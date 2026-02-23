import snappy

M=snappy.Manifold("s580")

curves = M.dual_curves()
for curve in M.dual_curves(max_segments=15):
	print(curve)
	N=M.drill(curve)
	print(N.identify())


