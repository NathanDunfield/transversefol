using Serialization
using DataFrames

include("search.jl")
include("find_surface.jl")

function mathematica_print(l::Union{Array,Tuple})
	print("{")
	for i in l[1:end-1]
		mathematica_print(i)
		print(",")
	end
	mathematica_print(l[end])
	print("}")
end

function mathematica_print(l::Real)
	print(Float64(l))
end

function find_s2_longitudes(isosig; nlongs=300)
	D=Set()
	include("batch/$(isosig).txt")

	global bt=BoundaryTriangulation(fans, face_coorientations, firstrungs, alledges, rungs)	
	global ncusps = length(bt.firstrungs)
	global Elong=PEnvelope()
	global longitudeDF = DataFrame()
	global long_dict = DefaultDict(()->[])
	global longitudes = []

	for l in find_longitudes_iterative(fans)
		if any([sum(abs.(x))==0 for x in slopes(Longitude(bt,l))])
			continue
		end
		ss=[y//x for (x,y) in slopes(Longitude(bt,l))]
		if normalizedchi(Longitude(bt,l))==2 && all(abs.(ss) .< 100)
			push!(D,ss)
		end
		if length(D)>=nlongs
			break
		end
	end
	mathematica_print(collect(D))
end

function compute_longitudes!(isosig; nlongs=10)
	global Elong=PEnvelope()
	global longitudeDF = DataFrame()
	global long_dict = DefaultDict(()->[])
	global longitudes = []

	for l in find_longitudes_random(fans)
		if any([sum(abs.(x))==0 for x in slopes(Longitude(bt,l))])
			continue
		end
		ss = [y//x for (x,y) in slopes(Longitude(bt,l))]

		if !haskey(long_dict, ss)
			@show length(long_dict)
			flush(stdout)
		end
		push!(long_dict[ss], l)
		if length(long_dict) >= nlongs
			break
		end
	end

	for (ss, ls) in long_dict
		_, i = findmin(x-> (count(y->y==0, x), sum(x.^2)), ls)
		c=longitude_to_candidate(bt,ls[i])
		push!(Elong, (ss, c))
		push!(longitudes, ls[i])
		push!(longitudeDF, (ss=slopes(Longitude(bt,ls[i])), l=ls[i], chi = -sum(ls[i])//2, normalizedchi = normalizedchi(Longitude(bt,ls[i]))))
	end
end
function setup(isosig)
	println("setting up $(isosig), search for $(nlongs) points on fibered face")
	flush(stdout)
	if !isfile("batch/$(isosig).txt")
		println("batch/$(isosig).txt not found, preparing it now")
		flush(stdout)
		run(`python3 prepare.py $(isosig)`)
	end
	include("batch/$(isosig).txt")

	global bt=BoundaryTriangulation(fans, face_coorientations, firstrungs, alledges, rungs)	
	global ncusps = length(bt.firstrungs)

	Eupper = Envelope{Upper}(copy(Elong.A))
	Elower = Envelope{Lower}(copy(Elong.A))

	serialize("batch/$(isosig).jls", (bt=bt, Eupper=Eupper, Elower=Elower, Elong=Elong, longitudes=longitudes))

	println("done setup")
	flush(stdout)
end
function regimen(E::Envelope{S}; ncusps = length(E.A[1][1])) where {S}
	@show ncusps
	if S==Upper
		target=[CLIP for i in 1:ncusps]
	elseif S==Lower
		target=[-CLIP for i in 1:ncusps]
	else
		@assert false
	end
	E = try_improve(E; nsubdivide=1, iters=30000, time=1000, target=target, radius=0.001)
	flush(stdout)
	if isinteractive()
		add_trace!(p, _plotjs(E))
	end
	E = try_improve(E; nsubdivide=0, iters=100000, time=2000, target=target, radius=0.001, beta=800)
	flush(stdout)
	if isinteractive()
		add_trace!(p, _plotjs(E))
	end
	E = try_improve(E; nsubdivide=1, iters=2000000, time=2000, target=target, radius=0.001, beta=1300)
	flush(stdout)
	if isinteractive()
		add_trace!(p, _plotjs(E))
	end
	return E
end

function runjob(i::Int; rt=false, ex=false, nlongs=10)
	include("batch/2cusp_manifest.txt")
	isosig=isosigs[i]
	runjob(isosig; rt=rt, ex=ex, nlongs=nlongs)
end

function runjob3d(isosig)
	setup(isosig, nlongs=100)
	global p=PlotlyJS.plot()
	if isinteractive()
		display(p)
	end

	Elower2, Eupper2 = extreme_candidates(bt)

	Eupper = Envelope{Upper}(copy(Elong.A))
	Elower = Envelope{Lower}(copy(Elong.A))
	Eupper = regimen(Eupper)
	Elower = regimen(Elower)


	#=
	randE = random_trials(bt, ntrials=1000000)
	randE2 = PEnvelope()
	for (x,c) in randE.A
		push!(Eupper, (approximant_all_slopes(c::Candidate; time=10000), c))
	end
	=#

	for x in Eupper2.A
		push!(Eupper, x)
	end
	for x in Elower2.A
		push!(Elower, x)
	end

	longitudeDF = DataFrame()
	for (ss, ls) in long_dict
		_, i = findmin(x-> (count(y->y==0, x), sum(x.^2)), ls)
		c=longitude_to_candidate(bt,ls[i])
		normchi = normalizedchi(Longitude(bt,ls[i]))
		push!(longitudeDF, (ss=slopes(Longitude(bt,ls[i])), x=ss[1], y=ss[2], l=ls[i], normchi = normchi, text=string((normchi=normchi,ss=ss))))
	end

	add_trace!(p, _plotjs(Eupper))

	PlotlyJS.savefig(p, "batch/$(isosig).html")
	serialize("batch/$(isosig).jls", (bt=bt, Eupper=Eupper, Elower=Elower, Elong=Elong, longitudes=longitudes))
	flush(stdout)
	return p
end

function runjob(isosig::String; rt=false, ex=false, reg=true, nlongs=10)
	#setup(isosig)

	setup(isosig, nlongs=nlongs)
	global p=PlotlyJS.plot()
	if isinteractive()
		display(p)
	end

	

	Eupper = Envelope{Upper}(copy(Elong.A))
	Elower = Envelope{Lower}(copy(Elong.A))

	if reg
		Eupper = regimen(Eupper)
		Elower = regimen(Elower)
	end

	if ex
		Elower2, Eupper2 = extreme_candidates(bt)

		for x in Elower2.A
			push!(Elower, x)
		end
		for x in Eupper2.A
			push!(Eupper, x)
		end
	end
	if rt
		randE = random_trials(bt,ntrials=200000)
		for (x,c) in randE.A
			push!(Eupper, c)
			push!(Elower, c)
		end
	end

	longitudeDF = DataFrame()
	for (ss, ls) in long_dict
		_, i = findmin(x-> (count(y->y==0, x), sum(x.^2)), ls)
		c=longitude_to_candidate(bt,ls[i])
		normchi = normalizedchi(Longitude(bt,ls[i]))
		push!(longitudeDF, (ss=slopes(Longitude(bt,ls[i])), x=ss[1], y=ss[2], l=ls[i], normchi = normchi, text=string((normchi=normchi,ss=ss))))
	end

	addtraces!(p, _plotjs(Elower, Eupper)...)
	#addtraces!(p, _plotjs(Elower2, Eupper2)...)
	#add_trace!(p, _plotjs(Elong; color=LONGITUDE_COLOUR))
	if ncusps == 2
		add_trace!(p, PlotlyJS.scatter(longitudeDF,x=:x, y=:y, marker=attr(line=attr(width=0), size=25 ./ log.(4 .- longitudeDF[!,:normchi]), color=LONGITUDE_COLOUR), text=:text, mode="markers"))

		long_ranges = [[minimum(filter(r -> abs(r) < CLIP, collect(x[i] for x in keys(long_dict)))),
							   maximum(filter(r -> abs(r) < CLIP, collect(x[i] for x in keys(long_dict))))] for i in 1:ncusps]

		xe = 0.5 * (long_ranges[1][2]-long_ranges[1][1])
		ye = 0.5 * (long_ranges[2][2]-long_ranges[2][1])
		update_xaxes!(p,range=[long_ranges[1][1]-xe, long_ranges[1][2]+xe], autorange=false, title="cusp 1 slope")
		update_yaxes!(p,range=[long_ranges[2][1]-ye, long_ranges[2][2]+ye], autorange=false, title="cusp 2 slope")
	end

	PlotlyJS.savefig(p, "batch/$(isosig).html")
	serialize("batch/$(isosig).jls", (bt=bt, Eupper=Eupper, Elower=Elower, Elong=Elong, longitudes=longitudes))
	flush(stdout)
	return p
end

function dump_points(isosig)
	tup = deserialize("/home/jonathan/Dropbox/jonathan/transversefol/batch/$(isosig).jls")

	print("upperpts=")
	mathematica_print([x for (x,y) in tup.Eupper.A])
	print("lowerpts=")
	mathematica_print([x for (x,y) in tup.Elower.A])

end

function quickview(i::Int)
	try
		include("batch/2cusp_manifest.txt")
		isosig=isosigs[i]
		quickview(isosig; index=i)
	catch e
		@show e
	end
end

function load_isosig(isosig)
	locations=[]
	push!(locations, "/home/jonathan/Dropbox/jonathan/transversefol/batch/$(isosig).jls")
	push!(locations, "/home/jonathan/engaging_sshfs/transversefol/batch/$(isosig).jls")

	locations = sort(filter(isfile, locations),by=mtime)
	if length(locations)==0
		println("not found")
	else
		println("loading from $(locations[end])")
	end

	return deserialize(locations[end])
end

function improve(isosig::String; index=0)
	tup = load_isosig(isosig)
	Eupper = tup.Eupper
	Elower = tup.Elower

	bt = tup.bt
	ncusps = length(bt.firstrungs)
end

function quickview(isosig::String; index=0)
	tup = load_isosig(isosig)

	Eupper = tup.Eupper
	Elower = tup.Elower


	dummy_candidate=random_candidate(tup.bt,0)
	if isosig == "eLMkbcddddedde_2100"
		for pt in [(-2,1/2), (-1, 1/3), (-1/2, 1/6), (-1/3, 1/9), (-1/4, 1/12), (-1/5, 1/15), (-1/6, 1/18)]
			push!(Eupper, (pt, dummy_candidate))
			push!(Elower, (map(x->-x,pt), dummy_candidate))
		end
	end


	bt = tup.bt
	ncusps = length(bt.firstrungs)

	Econstr=[PEnvelope() for i in 1:4]
	Econstr_all = PEnvelope()
	longitudeDF = DataFrame()
	long_slopes = []
	for l in tup.longitudes
		c=longitude_to_candidate(bt,l)
		L=Longitude(bt,l)
		ss=slopes(L)
		sss = [y//x for (x,y) in ss]
		push!(long_slopes, sss)
		push!(longitudeDF, (ss=slopes(L), l=l,  x=sss[1], y=sss[2],text=string((normchi=normalizedchi(L),ss=ss)), nchi = normalizedchi(L)))
		for (k,s) in enumerate(constraints(L))
			if all(!isnan(x) for x in s) && all(!isinf(x) for x in s)
				push!(Econstr[k], (s,c))
				push!(Econstr_all, (s,c))
			end
		end
	end

	Econstr_upper = Envelope{Upper}()
	Econstr_lower = Envelope{Lower}()

	push!(Econstr_lower, Econstr[1])
	push!(Econstr_lower, Econstr[3])

	push!(Econstr_upper, Econstr[2])
	push!(Econstr_upper, Econstr[4])


	long_ranges = [[minimum(filter(r -> abs(r) < CLIP, collect(x[i] for x in long_slopes))),
						   maximum(filter(r -> abs(r) < CLIP, collect(x[i] for x in long_slopes)))] for i in 1:ncusps]

	#=
	config = PlotConfig(modeBarButtonsToAdd=[
    "drawline",
    "drawopenpath",
    "drawclosedpath",
    "drawcircle",
    "drawrect",
    "eraseshape"
	])
	=#
	config = PlotConfig(modeBarButtonsToAdd=[
	])

	#info = run(pipeline(`cat veering_census_with_data.txt`, `grep $(isosig)`))
	layout = Layout(title="$(isosig)")

	p=PlotlyJS.plot(layout,config=config)

	dummy_candidate=random_candidate(tup.bt,0)

	addtraces!(p, _plotjs(Elower, Eupper)...)
	if length(Econstr_lower.A) > 0
		addtraces!(p, _plotjs(Econstr_lower, Envelope{Upper}([([CLIP,CLIP],dummy_candidate)]), color=OBSTRUCTION_COLOUR)...)
	end
	if length(Econstr_upper.A) > 0
		addtraces!(p, _plotjs(Envelope{Lower}([([-CLIP,-CLIP],dummy_candidate)]),Econstr_upper, color=OBSTRUCTION_COLOUR)...)
	end

	#add_trace!(p, _plotjs(tup.Elong, color=LONGITUDE_COLOUR))
	add_trace!(p, PlotlyJS.scatter(longitudeDF,x=:x, y=:y, marker=attr(line=attr(width=0), size=25 ./ log.(4 .- longitudeDF[!,:nchi]), color=LONGITUDE_COLOUR), text=:text, mode="markers"))


	#=
	for k in 1:4
		add_trace!(p, _plotjs(Econstr[k],color=OBSTRUCTION_COLOUR))
	end
	=#

	add_trace!(p, _plotjs(Econstr_all,color=OBSTRUCTION_COLOUR))


	#add_trace!(p, _plotjs(Econstr_upper, color=OBSTRUCTION_COLOUR))
	#add_trace!(p, _plotjs(Econstr_lower, color=OBSTRUCTION_COLOUR))

	#=
	randE = random_trials(bt, ntrials=10000)
	randE2 = PEnvelope()
	for (x,c) in randE.A
		push!(randE2, (approximant_all_slopes(c::Candidate; time=10000), c))
	end

	add_trace!(p, _plotjs(randE, color=LONGITUDE_COLOUR))
	add_trace!(p, _plotjs(randE2, color=LONGITUDE_COLOUR))
	=#

	xe = 0.5 * (long_ranges[1][2]-long_ranges[1][1])
	ye = 0.5 * (long_ranges[2][2]-long_ranges[2][1])


	update_xaxes!(p,range=[clip(long_ranges[1][1]-xe,4), clip(long_ranges[1][2]+xe,4)],autorange=false, title="cusp 1 slope")
	update_yaxes!(p,range=[clip(long_ranges[2][1]-ye,4), clip(long_ranges[2][2]+ye,4)], autorange=false, title="cusp 2 slope")



	PlotlyJS.savefig(p, "batch/$(index).html")
	#serialize("batch/$(isosig).jls", (bt=bt, Eupper=Eupper, Elower=Elower, Elong=Elong))
	flush(stdout)
	p

	#interesting example isosigs[63]
end

function review()
	include("batch/2cusp_manifest.txt")
	missing_isosigs = []
	interesting_isosigs = []
	for i in 1:10
		try
			tup = deserialize("/home/jonathan/engaging_sshfs/transversefol/batch/$(isosigs[i]).jls")
			L=length(tup.Elower) + length(tup.Eupper)
			println("$i total weight $(L)")
			if L > 2
				push!(interesting_isosigs,i)
			end
		catch error
			@show error
			println("$i missing")
			push!(missing_isosigs,i)
		end
	end
	@show missing_isosigs
	@show interesting_isosigs

end

function verify(isosig::String)
	tup = load_isosig(isosig)
	@threads for (x,c) in tup.Elower.A
		@show x, approximant_all_slopes(c)
	end

	@threads for (x,c) in tup.Eupper.A
		@show x, approximant_all_slopes(c)
	end
end

function aggregate_bounds()
	include("batch/2cusp_manifest.txt")
	D=DefaultDict(()->[])
	
	bad_D=DefaultDict(()->[])
	sharp_D=DefaultDict(()->[])

	for i in 1:100
		isosig = isosigs[i]
		f="/home/jonathan/engaging_sshfs/transversefol/batch/$(isosig).jls"
		if isfile(f)
			tup = deserialize(f)
			for l in sort(tup.longitudes, by=sum)
				L=Longitude(tup.bt,l)
				ss=slopes(L)

				q,p = ss[1]
				s,r = ss[2]

				g=gcd(p,q,r,s)

				npunctures = [gcd(a,b) for (a,b) in ss]
				chi = -sum(L.weights)//2
				closedchi = (chi + npunctures[1])//npunctures[2]
				#=
				if denominator(closedchi) !=1
					@show isosig
					@show L.weights
					@show ss, chi
					@show chi + npunctures[1] + npunctures[2]
					@show chi + npunctures[1]
					@show closedchi
					#@assert denominator(closedchi)==1
				end
				=#


				b1=rationalize(bound(tup.Elower, p//q))
				b2=rationalize(bound(tup.Eupper, p//q))
				if s!= 0 && q!=0

				for constr in constraints(L)[1:2]
					info = (i,closedchi=closedchi,b1=b1,b2=b2,ss=[y//x for (x,y) in ss],constr=constr)
					push!(D[(mod(r//s,1),closedchi)], info)
					if abs(constr[2]-b1) < 0.00001 || abs(constr[2]-b2) < 0.00001
						push!(sharp_D[mod(r//s,1)],info)
					end
					if  b1 < constr[2] < b2
						@show (b1, ss[2], b2)
						@show ss
						@show chi, closedchi, ss
						@show (mod(r//s,1),closedchi)
						println(i)
						push!(bad_D[(mod(r//s,1),closedchi)],info)
					end
				end
			end


			end
		end
	end

	for (k,v) in D
		sort!(v)
	end
	for (k,v) in bad_D
		sort!(v)
	end
	for (k,v) in sharp_D
		sort!(v)
	end
	return D, bad_D, sharp_D
end

#isosig = "siddhi2"
#isosig = "eLMkbcddddedde_2100"
#isosig = "gvLQQcdeffeffffaafa_201102"
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
#



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
