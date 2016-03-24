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
from multiprocessing import Pool
import multiprocessing
import getpass
multiprocessing.freeze_support()
N_cores = 3

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
	  
def processtarget(x):
	filename,fps,inpsmi = x
	ret = []
	pool = modpool[filename]
	smiles = smilespool[filename]
	with open(filename, 'rb') as fid:
		bnb = cPickle.load(fid)
	#for each wombat compound get nn
	for i,fp in enumerate(fps):
		sims = sorted(zip(DataStructs.BulkTanimotoSimilarity(fp,pool),smiles,pool),reverse=True)[:req]
		for sim in sims:
			pred = round(bnb.predict_proba([map(int,list(sim[2].ToBitString()))])[:,1],3)
			ret.append([u_name[filename[7:-4]],filename[7:-4],inpsmi[i],sim[1],str(round(sim[0],3)),str(pred)])
	return ret

#main
introMessage()
output_name = sys.argv[1] + '_out_ad.txt'
req = int(sys.argv[2])
of = open(output_name, 'w')
mods = open('mods.txt').read().splitlines()
print 'Total Number of Classes : ' + str(len(mods))
u_name = dict()
getUpName()
fps,smis = importQuery(sys.argv[1])
print 'Total Number of Query Molecules : ' + str(len(fps))
usr, pw, samples = 'lhm30', 'change', 0
conn = pymysql.connect(db='pidgin', user=usr, passwd=pw, host='localhost', port=3306)
matrix = []
modpool = dict()
smilespool = dict()
for mod in mods:	
	pool=[]
	smiles= []
	cur = conn.cursor()
	cur.execute("SELECT stdsmiles FROM actives WHERE UNIPROT = '"+mod[7:-4]+"';")
	for row in cur:
		pool.append(calcNormalFingerprints(row[0]))
		smiles.append(row[0])
	modpool[mod] = pool
	smilespool[mod] = smiles
#mods = [mod for mod in glob.glob('models/*.pkl')]
ad_tasks = [[mod,fps,smis] for mod in mods]
pool = Pool(processes=N_cores)  # set up resources
jobs = pool.imap_unordered(processtarget, ad_tasks)
for i, result in enumerate(jobs):
	for res in result:
		of.write('\t'.join(map(str,res)) + '\n')
pool.close()
pool.join()
of.close()