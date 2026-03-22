import snappy
import spherogram as sp
import identify0

T = sp.RationalTangle(-1, 2) + sp.RationalTangle(1, 5) + sp.RationalTangle(1, 5)
L = T.numerator_closure()
print((L.exterior().identify()))

M=L.exterior()

for i in identify0.pA_flows(M,maxdepth=1):
	print(i)
