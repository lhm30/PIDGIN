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
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.naive_bayes import BernoulliNB
import cPickle
import glob
import gc
from collections import Counter
import os
import sys
import numpy as np
from multiprocessing import Pool
import multiprocessing
multiprocessing.freeze_support()
N_cores = 10

def introMessage():
	print '=============================================================================================='
	print ' Author: Lewis Mervin\n Email:  lhm30@cam.ac.uk\n Supervisor: Dr. A. Bender.  Number of cores: ' + str(N_cores)
	print ' Address: Centre For Molecular Informatics, Dept. Chemistry, Lensfield Road, Cambridge CB2 1EW'
	print '==============================================================================================\n'
	return

def login():
	user = raw_input(" Enter Username for PIDGIN & BIOSYSTEMS DB [%s]: " % getpass.getuser())
	if not user:
		user = getpass.getuser()

	pprompt = lambda: (getpass.getpass(' Enter Password for DB: '), getpass.getpass(' Retype password: '))

	p1, p2 = pprompt()
	while p1 != p2:
		print(' Passwords do not match. Try again')
		p1, p2 = pprompt()
	samples = raw_input(" Enter Number of Samples: ")

	return user, p1, int(samples)
	
def ispwneeded():
	msg = " Calculate Pathway Enrichment from BioSystems? [y/n]: "
	pwneeded = raw_input(msg)
	while pwneeded not in ['y','n']:
		print(' Please type y for yes, or n for no. Try again')
		pwneeded = raw_input(msg)
	return pwneeded

def printprog(size,count,message):
	count = count+1
	percent = (float(count)/float(size))*100
	sys.stdout.write(message + ' : %3d%%\r' % percent)
	sys.stdout.flush()

#import user query
def importQuery(name):
	outproblem = open('problematic_smiles.smi','w')
	query = open(name).read().splitlines()
	matrix = []
	problem = 0
	for q in query:
		try:
			fp = calcFingerprints(q)
			gc.disable()
			matrix.append(fp)
			gc.enable()
		except:
			problem +=1
			outproblem.write(q + '\n')
	matrix = np.array(matrix, dtype=np.uint8)
	if problem > 0:
		print 'WARNING: ' + str(problem) + ' SMILES HAVE ERRORS'
		outproblem.close()
	else:
		outproblem.close()
		os.remove('problematic_smiles.smi')
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
		gc.disable()
		outfp.append(calcFingerprints(i[0]))
		gc.enable()
	return np.array(outfp, dtype=np.uint8)

def getRandomCompoundPredictions(metric):
	global usr, pw
	conn = pymysql.connect(db='pidgin', user=usr, passwd=pw, host='localhost', port=3306)
	cur = conn.cursor()
	cur.execute("SELECT "+metric+" FROM preds limit 100000;")
	preds = np.array(cur.fetchall())[:,0]
	return preds
	
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
	m = None
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
	if m is None:
		print ' ERROR: Please enter threshold!'
		quit()
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

def getPW():
	global models
	bsid_a = dict()
	conn = pymysql.connect(db='biosystems', user=usr, passwd=pw, host='localhost', port=3306)
	cur = conn.cursor()
	for m in models.keys():
		cur.execute("SELECT bsid FROM target_bsid WHERE target ='"+str(m)+"';")
		bsids = np.array(cur.fetchall(),dtype=int)
		try:
			bsid_a[m] = bsids[::,0]
		except IndexError:
			bsid_a[m] = []
	return bsid_a

#predict worker
def predict(x):
	global models
	global thresholds
	mod, input = x
	hits = 0
	probs = models[mod].predict_proba(input)[::,1]
	hits = probs > [thresholds[mod]]*len(probs)
	return [mod, hits.sum()]

#calculate enriched target metrics and calculate background pw array
def calculateEnrichmentT(bgpred):
	global bsid_a
	global positives
	print
	lwin = dict((el,0) for el in positives.keys())
	avr = dict((el,0) for el in positives.keys())
	bgpw = []
	#for each comparison
	for _ in range(samples):
		try:
			chunk = random.sample(bgpred,len(querymatrix))
		except ValueError:
			chunk = [random.choice(bgpred) for r in range(len(querymatrix))]
		printprog(samples,_,' Calculating Enriched Targets vs BG ')
		chunk = np.matrix(map(list,chunk),dtype=np.uint8)
		pw = dict()
		for i,mod in enumerate(sorted(models.keys())):
			hits = np.sum(chunk[:,i])
			if hits >= 1:
				#update count of hits for target (for average-ratio)
				avr[mod] = avr[mod] + hits
				for b in bsid_a[mod]:
					try:
						pw[b] += hits
					except KeyError:
						pw[b] = hits
			#update times that query was larger than background (for e-ratio)
			if positives[mod] > hits:
				lwin[mod] +=1
		bgpw.append(pw)
	return lwin, avr, bgpw
	
def calculateEnrichmentPW():
	global positivespw
	lwin = dict()
	avr = dict()
	pool = Pool(processes=N_cores)  # set up resources
	tasks = [[bsid, count] for bsid, count in positivespw.items()] #create queue
	jobs = pool.imap_unordered(processPW, tasks)
	for i, result in enumerate(jobs):
		lwin[result[0]]= result[1]
		avr[result[0]]= result[2]
	aratiopw = calcAR(avr,positivespw)
	return lwin,avr,aratiopw
	
def processPW(input):
	global bgpw
	lwin = 0
	avr = 0
	bsid, count = input
	for split in bgpw:
		try:
			split = split[bsid]
		except:
			split = 0
		if count > split:
			lwin +=1
		avr = avr + split
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
	
#main
introMessage()
usr, pw, samples = login()
metric = sys.argv[1]
print ' Using Class Specific Cut-off Thresholds of : ' + metric
thresholds = dict()
importThresholds()
file_name = sys.argv[2]
output_name, output_name2 = [file_name + 'out_targets_enriched.txt', file_name + 'out_pathways_enriched.txt']
models = trainModels()
u_name = dict()
getUpName()
bsid_a = getPW()
t_count = len(models.keys())
print ' Total Number of Classes : ' + str(t_count)
querymatrix = importQuery(file_name)
print ' Total Number of Library Molecules : ' + str(len(querymatrix))


positives = dict()
positivespw = dict()
pool = Pool(processes=N_cores)  # set up resources
test_prediction_tasks = [[mod, querymatrix] for mod in models.keys()] #create queue
jobs = pool.imap_unordered(predict, test_prediction_tasks)
for i, result in enumerate(jobs):
	mod, hit = result
	printprog(len(test_prediction_tasks),i,' Calculating Targets and Pathways for ' + file_name)
	positives[mod] = hit
	#update list of hit pw
	if hit >= 1:
		for b in bsid_a[mod]:
			try:
				positivespw[b] += hit
			except KeyError:
				positivespw[b] = hit
pool.close()
pool.join()

#import background db
bgpred = getRandomCompoundPredictions(metric)
#predict for random background, calculating number of times enriched in lib
lwin, avr, bgpw = calculateEnrichmentT(bgpred)
bgpred = None
#calculate average ratio
aratio = calcAR(avr,positives)
numpreds = float(len(querymatrix))

#write to target file
file = open(output_name, 'w')
file.write('uniprot\tname\tquery_hits\te_score\taverage_ratio\n')
for uniprot, hit in positives.items():
	if hit >=1:
		file.write(uniprot + '\t' + u_name[uniprot]  + '\t' +  str(hit)  + '\t' + str(1.0-(float(lwin[uniprot])/float(samples))) + '\t' + str(aratio[uniprot]) + '\n')
print '\n Wrote Target Results to : ' + output_name
file.close()

#run pathway analysis?
if ispwneeded() == 'n': quit()
	
#write to pw file
file = open(output_name2, 'w')
file.write('bsid\tname\tdatabase\texternal_id\tclass\tquery_hits\te_score\taverage_ratio\n')
lwin, avr, aratiopw = calculateEnrichmentPW()
conn = pymysql.connect(db='biosystems', user=usr, passwd=pw, host='localhost', port=3306)
cur = conn.cursor()
for bsid, count in positivespw.items():
	cur.execute("SELECT * FROM bsid_info WHERE bsid ='"+str(bsid)+"';")
	BSID_n = cur.fetchall()[0]
	file.write('\t'.join(map(str,BSID_n))  + '\t' +  str(count)  + '\t' + str(1.0-(float(lwin[bsid])/float(samples))) + '\t' + str(aratiopw[bsid]) + '\n')
print ' Wrote Pathway Results to : ' + output_name2
file.close()

