from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from sklearn.naive_bayes import BernoulliNB
import cPickle
import glob
import os
import sys
from pylab import save
import numpy as np
import pymysql
import getpass

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
	
def calcNormalFingerprints(smiles):
	m1 = Chem.MolFromSmiles(smiles)
	fp = AllChem.GetMorganFingerprintAsBitVect(m1,2, nBits=2048)
	return fp 

def importQuery(name):
	query = open(name).read().splitlines()
	matrix = []
	for q in query:
		matrix.append(calcNormalFingerprints(q))
	return matrix
	
#get names of uniprots
def getUpName():
	global u_name
	t_file = open('classes_in_model.txt').read().splitlines()
	t_file.pop(0)
	for t in t_file:
		t = t.split('\t')
		u_name[t[1]] = t[0]
	return
	
def introMessage():
	print '=============================================================================================='
	print ' Author: Lewis Mervin\n Email:  lhm30@cam.ac.uk\n Supervisor: Dr. A. Bender'
	print ' Address: Centre For Molecular Informatics, Dept. Chemistry, Lensfield Road, Cambridge CB2 1EW'
	print '==============================================================================================\n'
	return 
	  
#main
introMessage()
output_name = 'out_ad.txt'
file = open(output_name, 'w')
t_count = len(glob.glob('models/*.pkl'))
print 'Total Number of Classes : ' + str(t_count)
u_name = dict()
getUpName()
fp = importQuery(sys.argv[1])
print 'Total Number of Query Molecules : ' + str(len(fp))
usr, pw, samples = login()
conn = pymysql.connect(db='pidgin', user=usr, passwd=pw, host='localhost', port=3306)
count = 0
firstcols = []
matrix = []
for filename in glob.glob('models/*.pkl'):
	querynn = []
	pool = []
	cur = conn.cursor()
	cur.execute("SELECT stdsmiles FROM actives WHERE UNIPROT = '"+filename[7:-4]+"';")
	for row in cur:
		pool.append(calcNormalFingerprints(row[0]))
	with open(filename, 'rb') as fid:
		firstcols.append([u_name[filename[7:-4]],filename[7:-4]])
		bnb = cPickle.load(fid)
	#for each wombat compound get nn
	ncount = int(round(bnb.class_count_[1]*0.1))
	sim_array = []
	for f in fp:
		sim_array.append(DataStructs.BulkTanimotoSimilarity(f,pool))
	sims = []
	for sim in sim_array:
		sims.append(np.average(np.sort(sim)[-ncount:]))
	matrix.append(sims)
	count +=1
	#update precent finished
	percent = (float(count)/float(t_count))*100
	sys.stdout.write(' Performing NN search on Query Molecule: %3d%%\r' % percent)
	sys.stdout.flush()
matrix = np.concatenate((firstcols, matrix), axis=1)
headings = ['Name','Uniprot']
for i in range(len(fp)):
	headings.append("C"+str(i+1))
file.write('\t'.join(headings))
for row in matrix:
	file.write('\t'.join(map(str,row)) + '\n')
sys.stdout.write('Wrote Results to: ' + output_name + 20*' ')
file.close()