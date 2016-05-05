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
import math
import getpass
from multiprocessing import Pool
import multiprocessing
import getpass
multiprocessing.freeze_support()
N_cores = 25

def introMessage():
	print '=============================================================================================='
	print ' Author: Lewis Mervin\n Email:  lhm30@cam.ac.uk\n Supervisor: Dr. A. Bender'
	print ' Address: Centre For Molecular Informatics, Dept. Chemistry, Lensfield Road, Cambridge CB2 1EW'
	print '==============================================================================================\n'
	return 

def login():
    user = raw_input("Enter Username for PIDGIN DB [%s]: " % getpass.getuser())
    if not user:
        user = getpass.getuser()
    pprompt = lambda: (getpass.getpass(), getpass.getpass('Retype password: '))
    p1, p2 = pprompt()
    while p1 != p2:
        print('Passwords do not match. Try again')
        p1, p2 = pprompt()
    return user, p1
	
def calcNormalFingerprints(smiles):
	m1 = Chem.MolFromSmiles(smiles)
	fp = AllChem.GetMorganFingerprintAsBitVect(m1,2, nBits=2048)
	return fp 

def fp_array(x):
	mod, sm = x
	ret = []
	for s in sm:
		ret.append(calcNormalFingerprints(s))
	return [mod,ret]
		

def importQuery(name):
	smis = []
	query = open(name).read().splitlines()
	matrix = []
	for q in query:
		try:
			matrix.append(calcNormalFingerprints(q))
			smis.append(q)
		except:
			print 'err ' + q
			pass
	return matrix, smis
	
def getUpName():
	global u_name
	t_file = open('classes_in_model.txt').read().splitlines()
	t_file.pop(0)
	for t in t_file:
		t = t.split('\t')
		u_name[t[1]] = t[0]
	return
	  
def processtarget(x):
	filename,fps = x
	ret = []
	for i,fp in enumerate(fps):
		line = [u_name[filename[7:-4]],filename[7:-4],smis[i]]
		res = sorted(zip(DataStructs.BulkTanimotoSimilarity(fp,modfps[filename]),s_dict[filename]),reverse=True)
		for i in range(req):
			smi = res[i][1]
			line.append(smi)
			sim = round(res[i][0],3)
			line.append(sim)
		ret.append(line)
	return ret

#main
introMessage()
mods = glob.glob('models/*.pkl')
print 'Total Number of Classes : ' + str(len(mods))
u_name = dict()
getUpName()
fps,smis = importQuery(sys.argv[1])
req = int(sys.argv[2])
of = open(sys.argv[1] + '_out_ad_detailed_' + str(req) + '_nn.txt', 'w')
of.write('Name\tTarget\t' + '\t'.join(map(str,list(np.ravel([['Comp_' + str(i+1),'Sim_' + str(i+1)] for i in range(req)])))) + '\n')
print 'Total Number of Query Molecules : ' + str(len(fps))
usr, pw = login()
conn = pymysql.connect(db='pidgin', user=usr, passwd=pw, host='localhost', port=3306)
s_dict = dict()
print 'Gathering active compounds for all targets'
for j, mod in enumerate(mods):
	cur = conn.cursor()
	cur.execute("SELECT stdsmiles FROM actives WHERE UNIPROT = '"+mod[7:-4]+"';")
	s_dict[mod] = np.array(cur.fetchall())[:,0]

print 'Calculating fingerprints for all actives'
modfps = dict()
pool = Pool(processes=N_cores)  # set up resources
jobs = pool.imap_unordered(fp_array, [[mod,smiles] for mod, smiles in s_dict.iteritems()])
for i, result in enumerate(jobs):
	modfps[result[0]] = result[1]
pool.close()
pool.join()

print 'Calculating near-neighours for all input compounds'
ad_tasks = [[mod,fps] for mod in sorted(mods)]
pool = Pool(processes=N_cores)  # set up resources
jobs = pool.imap(processtarget, ad_tasks)
for i, result in enumerate(jobs):
	for r in result:
		of.write('\t'.join(map(str,r)) + '\n')
of.close()