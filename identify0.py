import snappy
import csv
import sqlite3
import subprocess

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




def iterate_pA_flows(M, maxdepth=2, max_segments=6):
	x = VeeringDB.identify(M)
	if x:
		yield x
	if maxdepth > 0:
		for i in range(len(M.dual_curves(max_segments=max_segments))):
			N=M.drill(i)
			for i in iterate_pA_flows(N, maxdepth=maxdepth-1, max_segments = max_segments):
				yield i

#s137(5, 4) s460(6, 1) s593(6, 1) s614(5, 1) s753(−6, 1) s956(4, 1)
#v1333(−5, 1) v3045(−4, 1) t06114(5, 1) t08155(−5, 1) o912518(−6, 1) o912544(−4, 3)
#o913679(−6, 1) o914675(−1, 5) o915066(−5, 1) o922743(7, 1) o930634(6, 1) o936699(7, 1)


def pA_flows(M, maxdepth=1, max_segments=6):
	D=dict()
	for s in iterate_pA_flows(M, maxdepth=maxdepth, max_segments=max_segments):
		print(s)
		D[s.name()]=s

	ret = []
	for isosig in sorted(list(D.values()),key=lambda x: (x.num_cusps(),x.volume())):
		for row in find(isosig.name()):
			ret.append(row)
	return ret
