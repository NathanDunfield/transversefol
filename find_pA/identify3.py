import snappy
import spherogram as sp


M=snappy.Manifold("s593")
M=snappy.Manifold("s776")
M=snappy.Manifold("o9_30634")
M.browse()

curves = M.dual_curves()
for i in range(8):
	print(M.drill(i).identify())

