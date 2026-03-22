import snappy
import spherogram as sp
import identify0

L=sp.Link(braid_closure=9*[1,2,3,4]+4*[4,3,2])
print((L.exterior().identify()))
M=L.exterior()
print(M.volume())

for i in identify0.pA_flows(M,maxdepth=2):
	print(i)
