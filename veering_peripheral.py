import veering
import regina
import veering.veering_tri
import veering.taut as taut
import veering.transverse_taut

import snappy
import snappy.snap.peripheral.peripheral as periph
import sys
sys.path.append("/home/jonathan/Dropbox/repo/Veering/scripts")
sys.set_int_max_str_digits(0)
import boundary_triangulation

from snappy.snap import t3mlite as t3m
from snappy.snap.peripheral import link, dual_cellulation

def peripheral_curve_from_snappy(dual_cell, snappy_data):
	D = dual_cell
	T = D.dual_triangulation
	M = T.parent_triangulation
	data = snappy_data
	weights = len(D.edges)*[0]
	print("ntetrahedra", len(M.Tetrahedra))
	for tet_index, tet in enumerate(M.Tetrahedra): #for each tetrahedron
		for x in t3m.EdgeFacePairs:
			print(x)
		for vert_index, V in enumerate(t3m.ZeroSubsimplices): #for each vertex in the tetrahedron
			triangle = tet.CuspCorners[V]
			sides = triangle.oriented_sides()
			for tri_edge_index, tet_edge in enumerate(link.TruncatedSimplexCorners[V]): #for each edge in the corresponding truncation triangle
				tet_face_index = t3m.ZeroSubsimplices.index(tet_edge ^ V)
				side = sides[tri_edge_index]
				global_edge = side.edge()
				if global_edge.orientation_with_respect_to(side) > 0:
					dual_edge = D.from_original[global_edge]
					weight = data[tet_index][4*vert_index + tet_face_index]
					weights[dual_edge.index] = -weight

	# Sanity check
	total_raw_weights = sum([sum(abs(x) for x in row) for row in data])
	assert 2*sum(abs(w) for w in weights) == total_raw_weights
	return weights


isosig="iLMzMPcbcdefghhhhhhhxxqdl_12211002"
x = taut.isosig_to_tri_angle(isosig)
v= veering.veering_tri.veering_triangulation(*x)
M=snappy.Manifold(isosig.split("_")[0])
N = t3m.Mcomplex(M)
C = link.LinkSurface(N)
D = dual_cellulation.DualCellulation(C)
cusp_indices, data = M._get_cusp_indices_and_peripheral_curve_data()
meridian = peripheral_curve_from_snappy(D, [data[i] for i in range(0, len(data), 4)])
longitude = peripheral_curve_from_snappy(D, [data[i] for i in range(2, len(data), 4)])

print(len(D.edges))





print(len(meridian))
print(len(longitude))


