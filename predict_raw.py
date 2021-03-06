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

#main
introMessage()
file_name = sys.argv[1]
t_count = len(glob.glob('models/*.pkl'))
print 'Total Number of Classes : ' + str(t_count)
output_name = 'out_results.txt'
file = open(output_name, 'w')
querymatrix = importQuery()
u_name = dict()
getName()
print 'Total Number of Query Molecules : ' + str(len(querymatrix))

count=0
#for each model
for filename in glob.glob('models/*.pkl'):
    count +=1
    #unpickle model
    with open(filename, 'rb') as fid:
        bnb = cPickle.load(fid)
        probs = bnb.predict_proba(querymatrix)
        row = [u_name[filename[7:-4]],filename[7:-4]]
        for prob in probs:
            row.append(prob[1])
        file.write('\t'.join(map(str,row)) + '\n')
    #update precent finished
    percent = (float(count)/float(t_count))*100
    sys.stdout.write('Performing Classification on Query Molecules: %3d%%\r' % percent)
    sys.stdout.flush()
print '\nWrote Results to: ' + output_name
file.close()