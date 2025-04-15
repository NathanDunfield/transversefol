from flipper import *
from flipper.kernel import *
from collections import defaultdict
from unionfind import unionfind
from functools import reduce
import operator
import snappy
import snappy.snap.peripheral.peripheral as periph
import veering
import veering.veering_tri
import regina
import sys
import numpy as np
import random
import math



#sys.path.append("/home/jonathan/Dropbox/repo/Veering/scripts")
sys.set_int_max_str_digits(0)
#import boundary_triangulation
#import prepare

def mod(i,n):
	return i%n

def make_triangulation(n,k): #n component chain link. Want to compose a bunch of vertical dehn twists and one horizontal dehn twist.
#We will build a big strip like this:
#		 i	   i+1
# -----------------------------
#		 |	  / |		 
#		 |   /  |		  
#		 |  /   |		  
#		 | /	|		  
# -----------------------------
#		 i	  i+1
#
#There are two negative Dehn twists separated by distance k.

	edge_count=[0]
	def new_edge():
		edge_count[0]+=1
		return edge_count[0]-1
	vertices = [flipper.kernel.Vertex(i) for i in range(n)]
	vertedges = [Edge(vertices[i], vertices[i], new_edge()) for i in range(n)]
	horedges = [Edge(vertices[i], vertices[mod(i+1,n)], new_edge()) for i in range(n)]
	diagedges = [Edge(vertices[i], vertices[mod(i+1,n)], new_edge()) for i in range(n)]

	triangles = []
	for i in range(n):
		triangles.append((vertedges[i], diagedges[i], horedges[i].reversed_edge)) #counterclockwise
		triangles.append((horedges[i], vertedges[mod(i+1,n)].reversed_edge, diagedges[i].reversed_edge))

	triangulation = Triangulation([Triangle(t) for t in triangles])

	def make_twist(edges, k=1):#raised to the power k
		indices=set([e.index for e in edges])
		weights=[1 if i in indices else 0 for i in range(edge_count[0])]
		assert triangulation.lamination(weights).is_curve()
		return triangulation.lamination(weights).encode_twist(k)
	
	m=math.floor(n/2-k/2)
	vertical_dts = [make_twist([diagedges[i], horedges[i]], -1 if (i == m or i == mod(m+k,n)) else 1) for i in range(n)] #first two twists left handed
	horizontal_dt = make_twist(vertedges + diagedges, -1)

	encoded_twists = vertical_dts + [horizontal_dt]
	
	phi = reduce(operator.mul, encoded_twists, triangulation.id_encoding())
	
	return triangulation, phi, encoded_twists

def find_s3_slope(M,unfilled_index, n=20):
	current_fillings=[x.filling for x in M.cusp_info()]
	if M.num_cusps()==1:
		slope_candidates = M.short_slopes()[0]
	else:
		slope_candidates = [(p,q) for p in range(-n,n) for q in range(0,n)]

	for ss in slope_candidates:
		M.dehn_fill(ss,unfilled_index)
		G=M.fundamental_group()
		if len(G.generators())==0:
			return ss
	assert False
	return None

def find_s3_slope2(M):
	ret = []
	for i in range(M.num_cusps()):
		M.dehn_fill((0,0),i)
	for ss in itertools.product(*M.short_slopes(length=6)):
		M.dehn_fill(list(ss))
		G=M.fundamental_group()
		print(ss)
		if len(G.generators())==0:
			ret.append(ss)
	return ret


def analyze(n, k, name="test"):
	print()
	tau,phi,twists = make_triangulation(n,k)
	print(str(len(tau.triangles)) + " triangles")
	#print(phi.nielsen_thurston_type())
	L=phi.invariant_lamination()
	
	#for i in L.geometric:
	#	print(i,i.minpoly())
	
	scale_factor = 5/np.max(np.array(L.geometric, dtype=np.float64))
	draw_triangulation(n,k, L.geometric, name="traintracks/traintrack({},{}).svg".format(n,k), scale_factor = scale_factor)

	ell=L
	draw_triangulation(n,k, ell.geometric, name="traintracks/folding_sequence({},{})_{}.svg".format(n,k,0), scale_factor = scale_factor)
	for (i,t) in enumerate(reversed(twists),1):
		ell = t(ell)
		draw_triangulation(n,k, ell.geometric, name="traintracks/folding_sequence({},{})_{}.svg".format(n,k,i), scale_factor = scale_factor)

	print(phi.stratum())
	bun=phi.bundle()
	s=bun.snappy_string()
	M=snappy.Manifold(bun)
	#print(M)
	print("volume: ", M.volume())
	print(M.identify())
	degen = bun.degeneracy_slopes()
	fibre = bun.fibre_slopes()

	print("degeneracy slopes: ", degen)
	print("fiber slopes: ", fibre)


	"""
	current_fillings = [x.filling for x in M.cusp_info()]
	unfilled_index=0
	for i in range(len(current_fillings)):
		if current_fillings[i] == (0.0,0.0):
			unfilled_index=i
			break


	try:
		find_s3_slope(M,unfilled_index)
		meridian = [x.filling for x in M.cusp_info()]


		for i in range(M.num_cusps()):
			print("cusp "+str(i))
			if i==unfilled_index:
				A=np.linalg.inv(np.transpose(np.array([fibre[i],meridian[i]])))
			else:
				A=np.linalg.inv(np.transpose(np.array([fibre[i],degen[i]])))

			print("degen_slope " + str(np.matmul(A,np.array(degen[i]))))
			print("fiber_slope " + str(np.matmul(A,np.array(fibre[i]))))
			print("s3_slope " + str(np.matmul(A,np.array(meridian[i]))))
	except:
		print("failed to find s3 slope")

	M=snappy.Manifold(bun)
	tri=regina.SnapPeaTriangulation(s)
	angle=[1 for i in range(M.num_tetrahedra())] #flipper arranges that all the tetrahedra are flattened in the same way
	longitude =[(x[0].label,x[1](3)) for x in bun.immersion.values()]

	v=veering.veering_tri.veering_triangulation(tri,angle)
	prepare.prepare_example(v, isosig=name, longitude=longitude)

	info_file = open("batch/" + name + ".info.txt",'w')
	print(phi.stratum(), file=info_file)
	print(bun.snappy_string(), file=info_file)

	info_file.close()
	"""
	return bun
import drawsvg as dw
import numpy as np

def distance(A,B):
	return np.sqrt(np.sum((A-B)**2))

def draw_triangle(d, A,B,C, wA, wB, wC, text=False):
	#A,B,C are the positions of the vertices
	#wA,wB,and wC are weights of the opposite edges

	#mA, mB, and mC are masses of the vertices

	semiperimeter = (distance(A,B) + distance(B,C) + distance(A,C))/2
	mA = 1/(semiperimeter - distance(B,C))
	mB = 1/(semiperimeter - distance(A,C))
	mC = 1/(semiperimeter - distance(A,B))

	fC = (mA * A + mB * B)/(mA+mB)
	fB = (mA * A + mC * C)/(mA+mC)
	fA = (mB * B + mC * C)/(mB+mC)


	def stroke(start, end, radius, width):
		print(start,end,radius,width)
		if width > 0.001:
			p=dw.Path(stroke='green',fill='none', stroke_width=width)
			p=p.M(*start).A(radius,radius, 0, False, True, *end)
			if text:
				d.append(dw.Text('%.2f' % width, 18, path=p))
			else:
				d.append(p)

	
	W=wA+wB+wC
	stroke(fB,fC, 1/mA, W-2*wA)
	stroke(fC,fA, 1/mB, W-2*wB)
	stroke(fA,fB, 1/mC, W-2*wC)

def draw_triangulation(n,k, weights, name, scale_factor = None):
	# Create an SVG drawing object
	drawing = dw.Drawing(200*(n+1),200)

	weights = np.array(weights, dtype=np.float64)

	if scale_factor == None:
		scale_factor = 20/np.max(weights)

	weights = scale_factor*np.array(weights)

	vert= weights[0:n]
	hor= weights[n:2*n]
	diag= weights[2*n:3*n]

#		 i	   i+1
# -----------------------------
#		 |	  / |		 
#		 |   /  |		  
#		 |  /   |		  
#		 | /	|		  
# -----------------------------
#		 i	  i+1

	upper_vertices = [200*np.array([i,0]) for i in range(n+1)]
	lower_vertices = [200*np.array([i+0.5,np.sqrt(3)/2]) for i in range(n+1)]


	for i in range(n):
		draw_triangle(drawing,upper_vertices[i+1], upper_vertices[i], lower_vertices[i], vert[i], diag[i], hor[i], text=False)
		draw_triangle(drawing,upper_vertices[i+1], lower_vertices[i], lower_vertices[i+1], hor[i], vert[mod(i+1,n)], diag[i], text=False)	
	for i in range(n):
		draw_triangle(drawing,upper_vertices[i+1], upper_vertices[i], lower_vertices[i], vert[i], diag[i], hor[i], text=True)
		draw_triangle(drawing,upper_vertices[i+1], lower_vertices[i], lower_vertices[i+1], hor[i], vert[mod(i+1,n)], diag[i], text=True)	

	drawing.save_svg(name)

# Call the function with an example value of N
for n in range(3,10):
	for k in range(1,math.floor(n/2)):
		analyze(n,k)
