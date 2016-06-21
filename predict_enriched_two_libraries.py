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
import operator
import numpy as np

def introMessage():
    print '=============================================================================================='
    print ' Author: Lewis Mervin\n Email:  lhm30@cam.ac.uk\n Supervisor: Dr. A. Bender'
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

def predict(input, name):
	results = dict()
	count=0
	#for each model
	for filename in glob.glob('models/*.pkl'):
		hits = 0
		count +=1
		#unpickle model
		with open(filename, 'rb') as fid:
			bnb = cPickle.load(fid)
			probs = bnb.predict_proba(input)
			for prob in probs:
				#if the probability of activity is above threshold then active
				if prob[1] >= thresholds[filename[7:-4]]:
					hits+=1
		results[filename[7:-4]] = hits
		#update precent finished
		percent = (float(count)/float(t_count))*100
		sys.stdout.write(' Performing Classification on '+name+' Molecules: %3d%%\r' % percent)
		sys.stdout.flush()
	print
	return results

def calculateEnrichment(positives,background):
	out = dict()
	for uniprot, hits in positives.items():
		if hits == 0:
			out[uniprot] = 999.0
			continue
		try:
			out[uniprot] = (float(hits)/float(len(querymatrix)))/(float(background[uniprot])/float(len(querymatrix2)))
		except ZeroDivisionError:
			out[uniprot] = 0.0
	return out

#main
introMessage()
file_name = sys.argv[1]
file_name2 = sys.argv[2]
metric = sys.argv[3]
print ' Using Class Specific Cut-off Thresholds of : ' + metric
t_count = len(glob.glob('models/*.pkl'))
print ' Total Number of Classes : ' + str(t_count)
outf = open(file_name + '_vs_' + file_name2 + '_out_results_enriched.txt','w')
thresholds = dict()
importThresholds()
u_name = dict()
getUpName()
querymatrix = importQuery(file_name)
querymatrix2 = importQuery(file_name2)
print ' Total Number of Library Molecules : ' + str(len(querymatrix))
print ' Total Number of Background Molecules : ' + str(len(querymatrix2))
positives = predict(querymatrix, file_name)
background = predict(querymatrix2, file_name2)
enrichedTargets = calculateEnrichment(positives,background)
#write to file
outf.write('Uniprot\tName\tHits\tBG_Hits\tOdds_Ratio\n')
for uniprot, rate in sorted(enrichedTargets.items(), key=operator.itemgetter(1)):
	if positives[uniprot] == 0: continue
	outf.write(uniprot + '\t' + u_name[uniprot] + '\t' + str(round(positives[uniprot]/len(querymatrix),2)) + '\t' + str(round(background[uniprot]/len(querymatrix),2)) + '\t' + str(rate) + '\n')
outf.close()
