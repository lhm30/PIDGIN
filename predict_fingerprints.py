#Author : Lewis Mervin lhm30@cam.ac.uk
#Supervisor : Dr. A. Bender
#All rights reserved 2014
#Protein Target Prediction Tool trained on SARs from PubChem (Mined 08/04/14) and ChEMBL18
#Molecular Descriptors : 2048bit Morgan Binary Fingerprints (Rdkit) - ECFP4
#Dependencies : rdkit, sklearn, numpy

#libraries
import pymysql
import random
random.seed(2)
import time
import getpass
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

	return user, p1
	
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
	return matrix, query
	
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
def predict(q):
	global models
	global thresholds
	bioact_profile = []
	pwfp = []
	for name, m in sorted(models.iteritems()):
		prob = m.predict_proba(q)[:,1]
		hit = prob > [thresholds[name]]
		bioact_profile.append(int(hit))
		if hit == True:
			try:
				for pw in bsid_a[name]:
					pwfp.append(pw)
			except KeyError: pass
	return bioact_profile, pwfp
	
#main
introMessage()
usr, pw = login()
metric = sys.argv[1]
file_name = sys.argv[2]
print ' Using Class Specific Cut-off Thresholds of : ' + metric
thresholds = dict()
importThresholds()
output_name, output_name2 = [file_name + 'out_targets_fingerprints.txt', file_name + 'out_pathways_fingerprints.txt']
models = trainModels()
u_name = dict()
getUpName()
bsid_a = getPW()
t_count = len(models.keys())
print ' Total Number of Classes : ' + str(t_count)
querymatrix, smiles = importQuery(file_name)
print ' Total Number of Library Molecules : ' + str(len(querymatrix))

allpw = []
pwfp = dict()
pool = Pool(processes=N_cores)  # set up resources
prediction_tasks = [q for q in querymatrix] #create queue
jobs = pool.imap(predict, prediction_tasks)
outf=open(output_name,'w')
outf.write('SMILES\t' + '\t'.join(map(str,sorted(models.keys()))) + '\n')
for i, result in enumerate(jobs):
	printprog(len(prediction_tasks),i,' Calculating Targets and Pathways for ' + file_name)
	bioact, pws = result
	outf.write(smiles[i] + '\t' + '\t'.join(map(str,bioact)) + '\n')
	pwfp[i] = pws
	allpw += pws
pool.close()
pool.join()
print ' Wrote Target Results to : ' + output_name
outf.close()

allpw = list(set(allpw))
allpwnames = []
conn = pymysql.connect(db='biosystems', user=usr, passwd=pw, host='localhost', port=3306)
cur = conn.cursor()
for pw in sorted(allpw):
	cur.execute("SELECT * FROM bsid_info WHERE bsid ='"+str(pw)+"';")
	allpwnames.append(cur.fetchall()[0])
outf2 = open(output_name2, 'w')
outf2.write('SMILES\t' + '\t'.join(map(str,sorted(allpw))) + '\n')
outf2.write('SMILES\t' + '\t'.join(map(str,sorted(allpwnames))) + '\n')
for smilescount,bsids in sorted(pwfp.iteritems()):
	bsidcount = Counter(bsids)
	hits = []
	for pw in sorted(allpw):
		try:
			hits.append(bsidcount[pw])
		except KeyError:
			hits.append(0)
	outf2.write(smiles[smilescount] + '\t' + '\t'.join(map(str,hits)) + '\n')
print ' Wrote Pathway Results to : ' + output_name2
outf2.close()


# conn = pymysql.connect(db='biosystems', user=usr, passwd=pw, host='localhost', port=3306)
# cur = conn.cursor()
# cur.execute("SELECT * FROM bsid_info WHERE bsid ='"+str(bsid)+"';")
# BSID_n = cur.fetchall()[0]

