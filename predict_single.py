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
import numpy as np

def introMessage():
    print '=============================================================================================='
    print ' Author: Lewis Mervin\n Email:  lhm30@cam.ac.uk\n Supervisor: Dr. A. Bender'
    print ' Address: Centre For Molecular Informatics, Dept. Chemistry, Lensfield Road, Cambridge CB2 1EW'
    print '==============================================================================================\n'
    return

#import user query
def importQuery():
    query = open(file_name).read().splitlines()
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
def getName():
    global u_name
    t_file = open('classes_in_model.txt').read().splitlines()
    t_file.pop(0)
    for t in t_file:
        t = t.split('\t')
        u_name[t[1]] = t[0]
    return

#import thresholds as specified by user    
def importThresholds(uniprot):
    t_file = open('thresholds.txt').read().splitlines()
    for t in t_file:
        t = t.split('\t')
        if t[0] == uniprot:
            thresholds = (float(t) for t in t[1:])
    return thresholds


#main
introMessage()
file_name = sys.argv[1]
t_count = len(glob.glob('models/*.pkl'))
print ' Total Number of Classes : ' + str(t_count)
output_name = 'out_result_single.txt'
file = open(output_name, 'w')
querymatrix = importQuery()
u_name = dict()
getName()
print ' Query Molecule : ' + file_name
file.write('NAME\tUNIPROT\tRAW_SCORE\tPRECISION\tF_SCORE\tRECALL\tACCURACY\t0.5\n')

count=0
#for each model
for filename in glob.glob('models/*.pkl'):
    row = []
    count +=1
    #unpickle model
    with open(filename, 'rb') as fid:
        row = [u_name[filename[7:-4]],filename[7:-4]]
        bnb = cPickle.load(fid)
        prob = bnb.predict_proba(querymatrix)[0][1]
        row.append(prob)
        #if the probability of activity is above threshold then active
        thresholds = importThresholds(filename[7:-4])
        for thresh in thresholds: 
            if prob >= thresh:
                row.append('1')
            else:
                row.append('0')
        file.write('\t'.join(map(str,row)) + '\n')
    #update precent finished
    percent = (float(count)/float(t_count))*100
    sys.stdout.write(' Performing Classification on Query Molecules: %3d%%\r' % percent)
    sys.stdout.flush()
print '\n Wrote Results to: ' + output_name
file.close()