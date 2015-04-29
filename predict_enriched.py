#Author : Lewis Mervin lhm30@cam.ac.uk
#Supervisor : Dr. A. Bender
#All rights reserved 2014
#Protein Target Prediction Tool trained on SARs from PubChem (Mined 08/04/14) and ChEMBL18
#Molecular Descriptors : 2048bit Morgan Binary Fingerprints (Rdkit) - ECFP4
#Dependencies : rdkit, sklearn, numpy

#libraries
import pymysql
import random
import time
import getpass
random.seed(2)
from functools import wraps
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.naive_bayes import BernoulliNB
import cPickle
import glob
from lxml import etree as ET
from collections import Counter
import os
import sys
import urllib2
import socket
import numpy as np
from multiprocessing import Pool
import multiprocessing
multiprocessing.freeze_support()
N_cores = 25
socket.setdefaulttimeout(400)
opener = urllib2.build_opener()
opener.addheaders = [('User-agent', 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_8_4) AppleWebKit/536.30.1 (KHTML, like Gecko) Version/6.0.5 Safari/536.30.1')]
opener.addheaders = [('Connection', 'keep-alive')]
opener.addheaders = [('Accept', 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8')]
opener.addheaders = [('Accept-Language', 'en-us')]

def introMessage():
	print '=============================================================================================='
	print ' Author: Lewis Mervin\n Email:  lhm30@cam.ac.uk\n Supervisor: Dr. A. Bender.  Number of cores: ' + str(N_cores)
	print ' Address: Centre For Molecular Informatics, Dept. Chemistry, Lensfield Road, Cambridge CB2 1EW'
	print '==============================================================================================\n'
	return

def retry(ExceptionToCheck, tries=20, delay=3, backoff=4, logger=None):
	def deco_retry(f):
		@wraps(f)
		def f_retry(*args, **kwargs):
			mtries, mdelay = tries, delay
			while mtries > 1:
				try:
					return f(*args, **kwargs)
				except ExceptionToCheck, e:
					msg = time.ctime() + " %s, Retrying in %d seconds...\n" % (str(e), mdelay)
					if logger:
						logger.warning(msg)
					else:
						print msg
					mtries -= 1
					mdelay *= backoff
			return f(*args, **kwargs)
		return f_retry
	return deco_retry

#attach retry decorator for GET links
@retry(Exception)
def urlopen_with_retry(url): 
	return opener.open(url, timeout=15).read()

def login():
	user = raw_input(" Enter Username for PIDGIN DB [%s]: " % getpass.getuser())
	if not user:
		user = getpass.getuser()

	pprompt = lambda: (getpass.getpass(' Enter Password for DB: '), getpass.getpass(' Retype password: '))

	p1, p2 = pprompt()
	while p1 != p2:
		print(' Passwords do not match. Try again')
		p1, p2 = pprompt()
	samples = raw_input(" Enter Number of Samples: ")

	return user, p1, int(samples)

def printprog(size,count,message):
	count = count+1
	percent = (float(count)/float(size))*100
	sys.stdout.write(message + ' : %3d%%\r' % percent)
	sys.stdout.flush()

#import user query
def importQuery(name):
	query = open(name).read().splitlines()
	matrix = []
	for q in query:
		matrix.append(calcFingerprints(q))
	matrix = np.array(matrix, dtype=np.uint8)
	return matrix
	
#calculate 2048bit morgan fingerprints, radius 2
def calcFingerprints(smiles):
	m1 = Chem.MolFromSmiles(smiles)
	fp = AllChem.GetMorganFingerprintAsBitVect(m1,2, nBits=2048)
	binary = fp.ToBitString()
	return list(binary)
	
def arrayFP(input):
	outfp = []
	for i in input:
		outfp.append(calcFingerprints(i[0]))
	return np.array(outfp, dtype=np.uint8)

def getRandomFP(needed):
	global usr, pw
	conn = pymysql.connect(db='pidgin', user=usr, passwd=pw, host='localhost', port=3306)
	cur = conn.cursor()
	cur.execute("SELECT stdsmiles FROM compounds WHERE RAND(12)<(SELECT (("+str(needed)+"/COUNT(*))*10) FROM compounds) ORDER BY RAND(12) LIMIT "+str(needed)+";")
	smiles = cur.fetchall()
	print ' Number of BG mols : ' + str(len(smiles)) + ' ' * 35
	return chunks(smiles, len(querymatrix))
	
#get names of uniprots
def getUpName():
	global u_name
	t_file = open('classes_in_model.txt').read().splitlines()
	t_file.pop(0)
	for t in t_file:
		t = t.split('\t')
		u_name[t[1]] = t[0]
	return

#import thresholds as specified by user	
def importThresholds():
	global thresholds
	global metric
	if metric == 'p':
		m = 1	
	if metric == 'f':
		m = 2
	if metric == 'r':
		m = 3	
	if metric == 'a':
		m = 4
	if metric == '0.5':
		m = 5
	t_file = open('thresholds.txt').read().splitlines()
	for t in t_file:
		t = t.split('\t')
		thresholds[t[0]] = float(t[m])
	return
		
#parallel train models	
def trainModels():
	models = dict()
	pool = Pool(processes=N_cores)  # set up resources
	train_tasks = [modelFile for modelFile in glob.glob('models/*.pkl')] #create queue
	jobs = pool.imap_unordered(trainer, train_tasks)
	t_job = len(train_tasks)
	for i, result in enumerate(jobs):
		models[result[0]] = result[1]
	pool.close()
	pool.join()
	return models
	
#trainer worker	
def trainer(x):
	with open(x, 'rb') as fid:
			loaded = cPickle.load(fid)
	return [x[7:-4], loaded]

#predict worker
def predict(x):
	global models
	mod, input = x
	results = dict()	
	hits = 0
	probs = models[mod].predict_proba(input)
	for prob in probs:
		#if the probability of activity is above threshold then active
		if prob[1] >= thresholds[mod]:
			hits+=1
	return [mod, hits]

#chunk random fp into jobs
def chunks(l, n):
	for i in xrange(0, len(l), n):
		yield l[i:i+n]

#calculate enriched target metrics and calculate background pw array
def calculateEnrichmentT(poolcomp):
	global positives
	global bsid_a
	lwin = dict((el,0) for el in positives.keys())
	avr = dict((el,0) for el in positives.keys())
	bgpw = []
	for i, slice in enumerate(poolcomp):
		printprog(samples,i,' Calculating Enriched Targets vs BG ')
		pw = []
		selected = arrayFP(slice)
		pool = Pool(processes=N_cores)  # set up resources
		background_prediction_tasks = [[mod, selected] for mod in models.keys()] #create queue
		jobs = pool.imap_unordered(predict, background_prediction_tasks)
		for i, result in enumerate(jobs):
			unip, hits = result
			if hits > 0:
				pw = pw + (bsid_a[unip] * hits)
			#update count of hits for target (for average-ratio)
			avr[unip] = avr[unip] + hits
			#update times that query was larger than background (for e-ratio)
			if positives[unip] > hits:
				lwin[unip] +=1
		pool.close()
		pool.join()
		bgpw.append(pw)
	return lwin, avr, bgpw
	
def calculateEnrichmentPW():
	global positivespw
	global bgpw
	lwin = dict()
	avr = dict()
	pool = Pool(processes=N_cores)  # set up resources
	tasks = [[bsid, count] for bsid, count in positivespw.items()] #create queue
	jobs = pool.imap_unordered(process, tasks)
	for i, result in enumerate(jobs):
		printprog(len(tasks),i,' Calculating Enriched PW vs BG')
		lwin[result[0]]= result[1]
		avr[result[0]]= result[2]
	aratiopw = calcAR(avr,positivespw)
	return lwin,avr,aratiopw
	
def process(input):
	global bgpw
	lwin = 0
	avr = 0
	bsid, count = input
	for split in bgpw:
		split = Counter(split)
		try:
			if count > split[bsid]:
				lwin +=1
			avr = avr + split[bsid]
		except KeyError:
			lwin +=1
	return [bsid,lwin,avr]
	
#calculate enrichment
def calcAR(avr,positiv):
	global samples
	aratio = dict()
	for annotation, bhits in avr.items():
		#average background hit ratio
		normhit = float(bhits)/float(samples)
		#number of predictions
		numpreds = float(len(querymatrix))
		try:
			#normed positive hit ratio / normed background hits
			aratio[annotation] = (float(normhit)/float(numpreds))/(float(positiv[annotation])/float(numpreds))
		except:
			if float(bhits) == 0.0:
				aratio[annotation] = 0.0
			if positiv[annotation] == 0.0:
				aratio[annotation] = 999.0
	return aratio	

def getBSID(term):
	bsidlist = []
	#send bsid term to esearch
	url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=biosystems&term=' + term + '&retmode=xml&RetMax=500'
	response = urlopen_with_retry(url)
	root = ET.fromstring(response)
	#search the xml for the gids
	for information in root.findall('IdList'):
		for bsids in information.findall('Id'):
			bsid=(bsids.text)
			bsidlist.append(bsid)
	#remove duplicates
	bsidlist = list(set(bsidlist))
	return bsidlist
	
def getBSIDname(bsids):
	BSID_n = dict()
	#send bsid term to esearch
	chunks = [bsids[x:x+120] for x in xrange(0, len(bsids), 120)]
	for chunk in chunks:
		url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=biosystems&id=' + ','.join(map(str,chunk))
		response = urlopen_with_retry(url)
		root = ET.fromstring(response)
		ids = root.xpath('//Id')
		info = [root.xpath('//Item[@Name="biosystemname"]'), root.xpath('//Item[@Name="externalid"]')]
		for i, id in enumerate(ids):
			BSID_n[(id.text)] = [info[0][i].text, info[1][i].text]
	return BSID_n
	
#main
introMessage()
usr, pw, samples = login()
metric = sys.argv[1]
print ' Using Class Specific Cut-off Thresholds of : ' + metric
thresholds = dict()
importThresholds()
file_name = sys.argv[2]
output_name, output_name2 = ['out_targets_enriched.txt', 'out_pathways_enriched.txt']
models = trainModels()
u_name = dict()
getUpName()
t_count = len(models.keys())
print ' Total Number of Classes : ' + str(t_count)
querymatrix = importQuery(file_name)
print ' Total Number of Library Molecules : ' + str(len(querymatrix))

positives = dict()
positivespw = []
bsid_a = dict()
pool = Pool(processes=N_cores)  # set up resources
test_prediction_tasks = [[mod, querymatrix] for mod in models.keys()] #create queue
jobs = pool.imap_unordered(predict, test_prediction_tasks)
for i, result in enumerate(jobs):
	printprog(len(test_prediction_tasks),i,' Calculating Targets and Pathways for ' + file_name)
	bsid_a[result[0]] = getBSID(result[0])
	positives[result[0]] = result[1]
	#update list of hit pw
	if result[1] > 0:
		positivespw = positivespw + (bsid_a[result[0]] * result[1])
pool.close()
pool.join()
positivespw = Counter(positivespw)

#import background db
poolcomp = getRandomFP(samples*len(querymatrix))
#predict for random background, calculating number of times enriched in lib
lwin, avr, bgpw = calculateEnrichmentT(poolcomp)
#calculate average ratio
aratio = calcAR(avr,positives)
numpreds = float(len(querymatrix))

#write to target file
file = open(output_name, 'w')
file.write('uniprot\tname\tquery_hits\tenriched_count\te_score\tbg_rate\taverage_ratio\n')
for uniprot, hit in positives.items():
	file.write(uniprot + '\t' + u_name[uniprot]  + '\t' +  str(hit)  + '\t' +  str(lwin[uniprot])  + '\t' + str(1.0-(float(lwin[uniprot])/float(samples))) + '\t' + str(float(avr[uniprot])) + '\t' + str(aratio[uniprot]) + '\n')
print '\n Wrote Target Results to : ' + output_name
file.close()
	
#write to pw file
file = open(output_name2, 'w')
file.write('bsid\tname\torigin\tquery_hits\tenriched_count\te_score\tbg_rate\taverage_ratio\n')
lwin, avr, aratiopw = calculateEnrichmentPW()
BSID_n = getBSIDname(positivespw.keys())
count = 0
for bsid, count in positivespw.items():
	count +=1
	printprog(count,i,' Writing Pathways to file ' + file_name)
	file.write(str(bsid) + '\t' + '\t'.join(map(str,BSID_n[bsid]))  + '\t' +  str(count)  + '\t' +  str(lwin[bsid])  + '\t' + str(1.0-(float(lwin[bsid])/float(samples))) + '\t' + str(float(avr[bsid])) + '\t' + str(aratiopw[bsid]) + '\n')
print ' Wrote Pathway Results to : ' + output_name
file.close()

