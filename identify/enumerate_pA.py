import snappy
import csv
import sqlite3
import subprocess
import itertools

import sys
sys.path.insert(0,"/home/jonathan/Dropbox/repo/Veering/scripts")
sys.path.insert(0,"/home/jonathan/Dropbox/repo/Veering")
sys.setrecursionlimit(100000)

import veering
import veering.taut
import sage
import sage.all
import ast


sys.path.insert(0,"/home/jonathan/Dropbox/jonathan/transversefol")
import prepare

def basic_hash(manifold, digits=6):
	return to_str_at_prec(manifold.volume(), digits) + " " + repr(manifold.homology())

def to_str_at_prec(x, d):
	return ('%.' + repr(d) + 'f') % x

def generate_census():
	subprocess.run(["rm", "veering_census.sq3"]) 
	conn = sqlite3.connect("veering_census.sq3")

	cur = conn.cursor()
	reader = csv.reader(open("veering_census_with_data.txt"), delimiter=" ")
	data = []
	for i, row in enumerate(reader):
		isosig = row[0].split("_")[0]
		M=snappy.Manifold(isosig)
		if i%500==0:
			print(i)

		cusps = M.cusp_info('is_complete').count(True)
		H = M.homology()
		betti = H.betti_number()
		torsion = [c for c in H.elementary_divisors() if c != 0]
		h=basic_hash(M)
		data.append((row[0], isosig, int(i), float(M.volume()), int(cusps), int(betti), str(torsion), h))

	cur.execute("CREATE TABLE census (name TEXT, triangulation TEXT, id INTEGER, volume REAL, cusps INT, betti INT, torsion TEXT, hash INT)")
	cur.executemany("INSERT INTO census VALUES(?,?,?,?,?,?,?,?)", data)
	conn.commit()

#generate_census()


VeeringDB=snappy.database.ManifoldTable(table='census', mfld_hash=basic_hash, db_path="veering_census.sq3")

def find(isosig):
	reader = csv.reader(open("veering_census_with_data.txt"), delimiter=" ")
	for row in reader:
		if isosig.split("_")[0] == row[0].split("_")[0]:
			yield row




"""
def iterate_pA_flows(M, maxdepth=2, max_segments=6):
	x = VeeringDB.identify(M)
	if x:
		yield x
	if maxdepth > 0:
		for i in M.length_spectrum(cutoff=5.0,grouped=False, include_words=True):
		#for i in range(len(M.dual_curves(max_segments=max_segments))):
			try:
				N=M.drill_word(i.word)
				for i in iterate_pA_flows(N, maxdepth=maxdepth-1, max_segments = max_segments):
					yield i
			except Exception as e:
				pass

"""

"""
yields N (drilled manifold), x (veering isosig), n (number of drillings)
"""
def iterate_pA_flows2(M, count=10, ndrill=2, maxlength=3):
	short_words = [(i.word,i.length.real()) for i in M.length_spectrum_alt(count=count, bits_prec=100)]

	#print("short words:", short_words)

	for i in range(ndrill+1):
		for words in itertools.combinations(short_words,i):
			#print(words)
			if sum(i[1] for i in words) < maxlength:
				N=M.drill_words([i[0] for i in words], bits_prec=100)
				x = VeeringDB.identify(N)
				if x:
					yield N, x.name(), len(words)

L=[
"s137(5,4)",
"s460(6,1)", 
"s593(6,1)", 
"s614(5,1)", 
"s753(6,-1)",
"s956(4,1)",
"v1333(5,-1)",
"v3045(4,-1)", 
"t06114(5,1)", 
"t08155(5,-1)", 
"o9_12518(6,-1)",
"o9_12544(4,-3)",
"o9_13679(6,-1)",
"o9_14675(1,-5)",
"o9_15066(5,-1)",
"o9_22743(7,1)", 
"o9_30634(6,1)", 
"o9_36699(7,1)"]


def pA_flows(MM, count=10, ndrill=2, maxlength=3):
	D=dict()
	for N,isosig,n in iterate_pA_flows2(MM, count=count, ndrill=ndrill, maxlength=maxlength):
		tri, angle = veering.taut.isosig_to_tri_angle(isosig)
		assert tri.isOriented() #important to do it this way, because veering fixes the orientation, and then we can pass to snappy without any problems
		M=snappy.Manifold(tri)
		isoms = N.is_isometric_to(M, return_isometries=True)
		assert len(isoms) >= 1
		for isom in isoms:
			l = [(0,0) for i in range(n)]
			for i in range(n):
				l[isom.cusp_images()[i]] = isom.cusp_maps()[i]*sage.all.vector([1,0])
			for i in range(n):
				if l[i][0] < 0:
					l[i] = (-l[i][0], -l[i][1])

			""" Some checks
			Mtmp = snappy.Manifold(tri)
			assert n == Mtmp.num_cusps()
			Mtmp.dehn_fill(l)
			assert Mtmp.is_isometric_to(MM)
			"""


			yield isosig + "_" + str(l).replace(" ", "")

#Takes a veering isosig with filling curves and checks whether prong counts are >= 2
def honest_pA(isosig):
	d=prepare.prepare_by_isosig(isosig)
	return all(x>=2 for x in d["prong_counts"])


#M=snappy.OrientableClosedCensus[0]
#M=snappy.Manifold("K9_48(0,1)")
#M=snappy.Manifold(L[4])

name = "K12n242(6,1)"
M=snappy.Manifold(name)


with open("batch/" + str(name) + "_pAflows.txt","w") as f:
	isosigs = list(dict.fromkeys(filter(honest_pA, pA_flows(M, count=6, ndrill=3, maxlength=3))))
	for isosig in isosigs:
		print(isosig, file=f)
#print(pA_flows(snappy.Manifold("L6a4"), maxdepth=1))

"""
import pandas as pd
import csv
df = pd.read_csv('/home/jonathan/Downloads/conjecture_data/floer/final_data/QHSpheres.csv')

for i in range(10):
	print(snappy.Manifold(df.iloc[i]["name"]),df.iloc[i]["L_space"])
	M=snappy.Manifold(df.iloc[i]["name"])
	isosig = df.iloc[i]["name"]
	with open("batch/" + isosig + "_pAflows.txt","w") as f:
		isosigs = list(dict.fromkeys(filter(honest_pA, pA_flows(M, count=6, ndrill=3, maxlength=4))))
		for isosig in isosigs:
			print(isosig, file=f)
			"""
