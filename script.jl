import Pkg
Pkg.activate(@__DIR__; io=devnull)

using Revise
using TransverseFol
using TransverseFol.Envelopes
using TransverseFol.VeeringCensus

using DataFrames
using StatProfilerHTML
using Dates
using Profile, PProf
using JSON

#isosig = "siddhi2"
#isosig = "eLMkbcddddedde_2100"
#isosig = "gvLQQcdeffeffffaafa_201102" #L6a5
L6a5 = "gvLQQcdeffeffffaafa_201102"
#isosig = "gLLAQcdecfffhsermws_122201"
#isosig = "fLLQcbecdeepuwsua_20102"
#isosig = "fLLQcbeddeehhbghh_01110"
#isosig = "challenge2"
#isosig = "gLLPQbefefefhhhhhha_011102"
#isosig = "gLLPQcdfefefuoaaauo_022110"
#isosig = "gLLPQcdfefefuoaaauo_022110"
#setup(isosig)

#=
isosig = "eLMkbcddddedde_2100"
setup(isosig)

Profile.init(n=10^7, delay=0.01)
@profile setup(isosig)
using ProfileView
if isinteractive()
	ProfileView.view()
end
=#



#add_trace!(p, _plotjs(L6a2E, fill=true))
#add_trace!(p, _plotjs(accurate_E,fill=true))
#scatter_envelope!(p, accurate_E)
#scatter!(p, [-2,-1,-1/2,0],[1/2,1/3,1/6,0])
#scatter!(p,[-2,-1,-1/2],[1/3,1/6,0])


#=
if false
	c=subdivide(subdivide(longitude_to_candidate(bt,longitude)))
	println(all_slopes(c, time=50000))
	#=
	Profile.init(n=10^7, delay=0.01)
	vals1, candidate = @profile annealing(c->objective(all_slopes(c)), c, x->jiggle(x,0.005), 900, 1800, 2000000, verbose=true)
	println(all_slopes(candidate, time=50000))
	println(objective(all_slopes(candidate, time=10000)))
	=#
end
=#


#=
Profile.init(n=10^7, delay=0.01)
vals2, candidate = @profile annealing(c->objective2(all_slopes(c)), c, x->jiggle(x,0.03), 3000, 3000, 100000)
println(all_slopes(candidate, time=10000))
=#

#=
Profile.init(n=10^7, delay=0.01)
@profile random_trials(bt)
=#

#=
xs=0:0.01:1
p=plot(xs, [[f(x) for x in xs] for (J,f) in candidate.d if !J.inv])
display(p)
=#


#siddhi's example
#longitude  = 1/36
#meridian = 0
#46 triangles => 69 edges => H_1 has rank 24 => genus 12
#So the most we should expect is (36,1) - (2*12-1,0) = (13,1)
#1/13
#Looks like we're getting (1/3,1/3)
#Question is whether we can get (1/3+eps, 1/13-eps)
#Records: (1/13, 4/11) (1/3,1/3)



#bojun's example
#longitude = 1/48,   -1/4
#genus = 14
#prediction that we can get to 1/21, -1/4
#But we're getting 1/18



#L6a2
#Records
#(-2,1/2)   (-1, 1/3)    (-1/2, 1/6)    (-1/3, 1/9)    (-1/6, 1/18)
#
#
#pts = [(-2,1/2),   (-1, 1/3) ,   (-1/2, 1/6),    (-1/3, 1/9),    (-1/6, 1/18)]

#(0,0) is triangulated with 4 ideal triangles, has at least 2 punctures
# V - 3/2 *4 + 4 = V - 6 + 4 = V - 2. So V=2 => torus, V=4 => sphere. It's either a twice punctured torus or a 4 times punctured sphere
#
#
#



# m125
# Records
# [[0.5, -2.0], [0.014492753623188406, -1.3714285714285714], [-1.375, 0.06666666666666667], [0.25, -1.5], [-1.3703703703703705, 0.014492753623188406], [0.0136986301369863, -1.368421052631579], [-2.0, 0.5], [-1.4, 0.1111111111111111], [0.1111111111111111, -1.4], [0.07142857142857142, -1.375], [0.0, 0.0], [-1.5, 0.25]]




#=
dummy_candidate=random_candidate(bt,0)
L6a2E=Envelope()
push!(L6a2E, (T[0,0], dummy_candidate))
for pt in [(-4,1/2+0.000001), (-2,1/2), (-1, 1/3), (-1/2, 1/6), (-1/3, 1/9), (-1/4, 1/12), (-1/5, 1/15), (-1/6, 1/18)]
	push!(L6a2E, (T[pt...], dummy_candidate))
end
=#


#problematic examples
#11 - endpoint
#38 - disaster
#35 - endpoint
#49 - endpoint
#71 - 1-prong surgery
#83 - not fibered
#91 - endpoint
#95 - disaster
#26 - not fibered
#160
#194 - 1-prong surgery
#214
#244 - 1-prong surgery

#bad_examples = [38, 95, 160, 214, 278, 338, 356, 370, 406, 448, 453, 470, 473, 485]


#For L6a5:
#Know from magic.py that in snappy coordinates
#
#degeneracy = (0,-1), (0,1), (-1,1)
#fiber = (4,1), (1,-1), (-1,0)
#meridian = (-1,0), (0,1), (-1,1)
#
#So in my coordinates, we know that the meridian is (0, infty, infty)
#Now we want to change coordinates, so that 0->infty, and
