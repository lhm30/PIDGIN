#Author : Lewis Mervin lhm30@cam.ac.uk
#Supervisor : Dr. A. Bender
#All rights reserved 2014
#Protein Target Prediction Tool trained on SARs from PubChem (Mined 08/04/14) and ChEMBL18
#Molecular Descriptors : 2048bit Morgan Binary Fingerprints (Rdkit) - ECFP4
#Dependencies : rdkit, sklearn, numpy

#libraries
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.naive_bayes import BernoulliNB
import cPickle
import glob
import os
import sys
import pymysql
import numpy as np
from multiprocessing import Pool
import multiprocessing
import getpass
multiprocessing.freeze_support()	# Allows entire job to die
N_cores = 25

def login():
    user = raw_input(" Enter Username for PIDGIN DB [%s]: " % getpass.getuser())
    if not user:
        user = getpass.getuser()

    pprompt = lambda: (getpass.getpass(), getpass.getpass(' Retype password: '))

    p1, p2 = pprompt()
    while p1 != p2:
        print(' Passwords do not match. Try again')
        p1, p2 = pprompt()
    samples = raw_input(" Enter Number of Samples: ")

    return user, p1, int(samples)

def introMessage():
	print '=============================================================================================='
	print ' Author: Lewis Mervin\n Email:  lhm30@cam.ac.uk\n Supervisor: Dr. A. Bender.  Number of cores: ' + str(N_cores)
	print ' Address: Centre For Molecular Informatics, Dept. Chemistry, Lensfield Road, Cambridge CB2 1EW'
	print '==============================================================================================\n'
	return

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
	
def getRandomFP(needed):
	global usr, pw
    conn = pymysql.connect(db='pidgin', user=usr, passwd=pw, host='localhost', port=3306)
    cur = conn.cursor()
	#cur.execute("SELECT stdsmiles FROM compounds limit "+str(needed)+";")
	cur.execute("SELECT stdsmiles FROM compounds WHERE RAND(12)<(SELECT (("+str(needed)+"/COUNT(*))*10) FROM compounds) ORDER BY RAND(12) LIMIT "+str(needed)+";")
        #cur.execute("SELECT stdsmiles FROM compounds ORDER BY RAND(123) LIMIT "+str(needed)+";")
	smiles = cur.fetchall()
	print ' Number of BG mols : ' + str(len(smiles))
	return chunks(smiles, len(querymatrix))

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

def arrayFP(input):
	outfp = []
	for i in input:
		outfp.append(calcFingerprints(i[0]))
	return np.array(outfp, dtype=np.uint8)

def calculateEnrichment(poolcomp):
	global positives
	lwin = dict((el,0) for el in positives.keys())
	avr = dict((el,0) for el in positives.keys())
	background = dict()
	for i, slice in enumerate(poolcomp):
		selected = arrayFP(slice)
		percent = (float(i)/float(samples))*100 	
		sys.stdout.write('\r Calculating Enrichment : %3d%%' % percent)
		sys.stdout.flush()
		pool = Pool(processes=N_cores)  # set up resources
		background_prediction_tasks = [[mod, selected] for mod in models.keys()] #create queue
		jobs = pool.imap_unordered(predict, background_prediction_tasks)
		for i, result in enumerate(jobs):
			unip, hits = result
			background[unip] = hits
			#update count of hits for target (for average-ratio)
			avr[unip] = avr[unip] + hits
			#update times that query was larger than background (for e-ratio)
			if positives[unip] > hits:
				lwin[unip] +=1
		pool.close()
		pool.join()
	return lwin, avr

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
	
#calculate enrichment
def calcAR(avr):
	aratio = dict()
	for uniprot, bhits in avr.items():
		#average background hit ratio
		normhit = float(bhits)/float(samples)
		#number of predictions
		numpreds = float(len(querymatrix))
		try:
			#normed positive hit ratio / normed background hits
			aratio[uniprot] = (float(normhit)/float(numpreds))/(float(positives[uniprot])/float(numpreds))
		except:
			if float(bhits) == 0.0:
				aratio[uniprot] = 0.0
			if positives[uniprot] == 0.0:
				aratio[uniprot] = 999.0
	return aratio

	
#main
introMessage()
usr, pw, samples = login()
conn = pymysql.connect(db='pidgin', user=usr, passwd=pw, host='localhost', port=3306)
metric = sys.argv[1]
file_name = sys.argv[2]
print ' Using Class Specific Cut-off Thresholds of : ' + metric
models = trainModels()
t_count = len(models.keys())
print ' Total Number of Classes : ' + str(t_count)
output_name = 'out_results_enriched.txt'
thresholds = dict()
importThresholds()

#import and predict for library
querymatrix = importQuery(file_name)
print ' Total Number of Library Molecules : ' + str(len(querymatrix))
positives = dict()
pool = Pool(processes=N_cores)  # set up resources
test_prediction_tasks = [[mod, querymatrix] for mod in models.keys()] #create queue
jobs = pool.imap_unordered(predict, test_prediction_tasks)
t_job = len(test_prediction_tasks)
for i, result in enumerate(jobs):
	positives[result[0]] = result[1]
	percent = 1+((float(i)/float(t_job))*100.0)
	sys.stdout.write('\r Performing Classification on Library Molecules: %3d%%' % percent)
	sys.stdout.flush()
pool.close()
pool.join()

print '\n Number of Samples : ' + str(samples)
#import background db
poolcomp = getRandomFP(samples*len(querymatrix))
#predict for random background, calculating number of times enriched in lib
lwin, avr = calculateEnrichment(poolcomp)
#calculate average ratio
aratio = calcAR(avr)
numpreds = float(len(querymatrix))

#get uniprot labels for file
u_name = dict()
getUpName()

#write to file
file = open(output_name, 'w')
file.write('\nuniprot\tname\tquery_hits\tenriched_count\te_score\tbg_rate\taverage_ratio\n')
for uniprot, hit in positives.items():
	file.write(uniprot + '\t' + u_name[uniprot]  + '\t' +  str(hit)  + '\t' +  str(lwin[uniprot])  + '\t' + str(1.0-(float(lwin[uniprot])/float(samples))) + '\t' + str(float(avr[uniprot])) + '\t' + str(aratio[uniprot]) + '\n')
print '\n\n Wrote Results to : ' + output_name
file.close()
