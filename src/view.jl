
using Blink
#=
@eval AtomShell begin
    function init(; debug = false)
        electron() # Check path exists
        p, dp = port(), port()
        debug && inspector(dp)
        dbg = debug ? "--debug=$dp" : []
        proc = (debug ? run_rdr : run)(
            `$(electron()) --no-sandbox $dbg $mainjs port $p`; wait=false)
        conn = try_connect(ip"127.0.0.1", p)
        shell = Electron(proc, conn)
        initcbs(shell)
        return shell
    end
end
=#

function quickview(tup::NamedTuple; longitudes=true, obstructions=true, contact_structures=false, font_size=30, save_html=true, save_png=false,png_width::Int=1920, png_height::Int=1080, png_scale::Real=1)
    isosig = tup.isosig
    index = VeeringCensus.index(isosig)

    #include("batch/$(isosig).txt")
	Eupper = tup.Eupper
	Elower = tup.Elower


    #=
	if isosig == "eLMkbcddddedde_2100"
        dummy_candidate=random_cand(tup.bt, 1, DOWN)
		for pt in [(-2,1/2), (-1, 1/3), (-1/2, 1/6), (-1/3, 1/9), (-1/4, 1/12), (-1/5, 1/15), (-1/6, 1/18)]
			push!(Eupper, (pt, dummy_candidate))
			push!(Elower, (map(x->-x,pt), dummy_candidate))
		end
	end
    =#


    @show length(tup.longitudes)


	bt = tup.bt
    ncusps=bt.ncusps

    function namedtuple(slopes)
        @assert length(slopes)==ncusps
        if ncusps==1
            return (x=slopes[1], y=0)
        elseif ncusps==2
            return (x=slopes[1], y=slopes[2])
        elseif ncusps==3
            return (x=slopes[1], y=slopes[2], z=slopes[3])
        end
        @assert false
    end

    function plotting_directives()
        if ncusps ==1 || ncusps == 2
            return (x=:x, y=:y, type="scatter")
        elseif ncusps==3
            return (x=:x, y=:y, z=:z, type="scatter3d")
        end
    end


    Econstr = PEnvelope()
    Econstr_upper = Envelope{Upper,Float64,Cand{DiscreteHomeo{Tuple{Int,Int}}}}()
    Econstr_lower = Envelope{Lower,Float64,Cand{DiscreteHomeo{Tuple{Int,Int}}}}()
    longitudeDF = DataFrame()
    constrDF = DataFrame()
    long_slopes = []
    for l in tup.longitudes
        c=longitude_to_candidate(bt,l)
        L=Longitude(bt,l)
        ss=slopes(L)
        sss = map(slope_to_rat, ss)


        #b1=bound(tup.Eupper, sss[1])
        #b2=bound(tup.Elower, sss[1])
        push!(long_slopes, sss)
        push!(longitudeDF, (ss=slopes(L), l=l,  namedtuple(sss)..., text=string((normchi=normalizedchi(L),ss=ss,weights=l)), nchi = normalizedchi(L)))

        if is_fiber(l,tup.prep.top_bot_pairs)
            #@assert connected_components(l,fans)==1
            for (s,info) in constraints(L)
                if all(!isnan(x) for x in s) && all(!isinf(x) for x in s) && info.npunc==1# && info.interior_prong >= 2
                    push!(Econstr, (s,c))
                    push!(constrDF, (namedtuple(s)...,text=string(info)))
                    if info.dir[2] == -1
                        push!(Econstr_upper, (s,c))
                    else
                        @assert info.dir[2] == 1
                        push!(Econstr_lower, (s,c))
                    end

                end
            end
        else
            #println("nonfiber: $(sss)")
        end
    end



    long_ranges = [[minimum(filter(r -> abs(r) < CLIP/5, collect(x[i] for x in long_slopes)), init=CLIP),
                           maximum(filter(r -> abs(r) < CLIP/5, collect(x[i] for x in long_slopes)), init=-CLIP)] for i in 1:ncusps]

    paddings = [0.5 * (y-x) for (x,y) in long_ranges]
    trimmed_ranges = [[range[1]-pad, range[2]+pad] for (range,pad) in zip(long_ranges, paddings)]

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
    config = PlotConfig(modeBarButtonsToAdd=[])

    #info = run(pipeline(`cat veering_census_with_data.txt`, `grep $(isosig)`))
    
    axes = if ncusps == 1
        (
        xaxis=attr(
            showticklabels=true,
            range=trimmed_ranges[1]
        ),
        yaxis=attr(
            showticklabels=false,
            range=[-1,1]
           ))
    elseif ncusps == 2
        (
        xaxis=attr(
            showticklabels=false,
            range=trimmed_ranges[1],
            title="cusp 1 slope"
        ),
        yaxis=attr(
            showticklabels=false,
            range=trimmed_ranges[2],
            title="cusp 2 slope"
           ))
    elseif ncusps == 3
        attr(scene=(
        xaxis=attr(
            showticklabels=false,
            range=trimmed_ranges[1]
        ),
        yaxis=attr(
            showticklabels=false,
            range=trimmed_ranges[2]
           ),

        zaxis=attr(
            showticklabels=false,
            range=trimmed_ranges[3]
           )))

    end

    data = VeeringCensus.lookup_row(index)

    layout = Layout(title=attr(text="#$(index)   $(data[:isosig])   $(data[:depth])    $(data[:names])", font=attr(size=font_size));
                    axes...,
                    font=attr(size=font_size))

    p=PlotlyJS.plot(layout)

    #contact structures
    if contact_structures
        dummy_candidate = random_cand(bt, 1, DOWN)
        addtraces!(p, _plotjs(Elower, Envelope{Upper,Float64,Cand{DiscreteHomeo}}([([CLIP for i in 1:ncusps], dummy_candidate)]), color=NEG_CONTACT_COLOUR, name="negative contact structures")...)

        addtraces!(p, _plotjs(Envelope{Lower,Float64,Cand{DiscreteHomeo}}([([-CLIP for i in 1:ncusps],dummy_candidate)]), Eupper, color=POS_CONTACT_COLOUR, name="positive contact structures")...)
    end

    if length(Elower.A) > 0 && length(Eupper.A) > 0
        addtraces!(p, _plotjs(Elower, Eupper, name="Foliation multislopes")...)
    end

    if haskey(tup, :Elowerbound)
        addtraces!(p, _plotjs(tup.Elowerbound, tup.Eupperbound, name="bound")...)
    end

    if obstructions
        if length(Econstr_lower.A) > 0
            addtraces!(p, _plotjs(Econstr_lower, Envelope{Upper}([([CLIP for i in 1:ncusps], nothing)]), color=OBSTRUCTION_COLOUR, name="obstructions")...)
        end
        if length(Econstr_upper.A) > 0
            addtraces!(p, _plotjs(Envelope{Lower}([([-CLIP for i in 1:ncusps],nothing)]), Econstr_upper, color=OBSTRUCTION_COLOUR, name="obstructions")...)
        end
    end


    function clip_df(df)
        return df
        return subset(df, :x => x->abs.(x).<=CLIP, :y => y->abs.(y).<=CLIP)
    end

    if length(tup.Elong.A) > 0
        #add_trace!(p, _plotjs(tup.Elong, color=LONGITUDE_COLOUR))
        if longitudes
            add_trace!(p, PlotlyJS.scatter(clip_df(longitudeDF); plotting_directives()..., marker=attr(line=attr(width=0), size=(ncusps<=2 ? 25 : 10) ./ log.(4 .- longitudeDF[!,:nchi]), color=LONGITUDE_COLOUR), text=:text, mode="markers", name="Fibration multislopes"))
        end
        if obstructions
            add_trace!(p, PlotlyJS.scatter(clip_df(constrDF); plotting_directives()..., marker=attr(color=OBSTRUCTION_COLOUR, size=(ncusps<=2 ? 5 : 3)), text=:text, mode="markers", name="obstructions"))
        end
    else
        println("no longitudes")
    end



    #=
    crevices = PEnvelope()
    for y in [(Vector{T}(x).+0.01, dummy_candidate) for x in crevices_general(Econstr_lower)]
        push!(crevices, y)
    end
    add_trace!(p, _plotjs(crevices, color=OBSTRUCTION_COLOUR))
    =#

    #add_trace!(p, _plotjs(Econstr_all, color=OBSTRUCTION_COLOUR))


    #add_trace!(p, _plotjs(Econstr_upper, color=OBSTRUCTION_COLOUR))
    #add_trace!(p, _plotjs(Econstr_lower, color=OBSTRUCTION_COLOUR))

    if false
        randE = random_trials(bt, nsubdivide=3, ntrials=1000000)
        randE2 = PEnvelope()
        for (x,c) in randE.A
            push!(randE2, (approximant_all_slopes(c::Candidate; time=10000), c))
        end

        #add_trace!(p, _plotjs(randE, color=TAUT_COLOUR))
        add_trace!(p, _plotjs(randE2, color=TAUT_COLOUR, name="random sample"))
    end

    if save_html
        PlotlyJS.savefig(p, joinpath(BATCH_DIR, "$(index).html"))
    end
    if save_png
        PlotlyJS.savefig(p, joinpath(BATCH_DIR, "$(index).png"), width=png_width, height=png_height, scale=png_scale)
    end
	flush(stdout)

    on(p["click"]) do data
        @show data
        for point in data["points"]
            @show point
            if ncusps == 1
                coords = [point["x"]]
            elseif ncusps == 2
                coords = [point["x"],point["y"]]
            elseif ncusps == 3
                coords = [point["x"],point["y"],point["z"]]
            end

            println(rationalize.(coords))

            for l in sort(tup.longitudes, by=sum)
                L=Longitude(tup.bt,l)

                ss=slopes(L)
                if rationalize.(coords)==map(slope_to_rat, ss)
                    println(constraints(L))
                end
            end


            for (s,cand) in Iterators.flatten([tup.Elong.A, tup.Eupper.A, tup.Elower.A])
                if s == coords
                    global lastcand = cand
                    for i in 1:ncusps
                        draw(cand,i) |> display
                    end
                    break
                end
            end
        end
    end

    return p

	#interesting example isosigs[63]
end



function viewladderpole(i::Int)
    run(`evince $(joinpath(BATCH_DIR, "$(VeeringCensus.lookup(i)).pdf"))`)
end

function viewladderpole(i::Int, ncusps::Int)
    run(`evince $(joinpath(BATCH_DIR, "$(VeeringCensus.lookup(i,ncusps)).pdf"))`)
end

function viewladderpole(isosig::String)
    run(`evince $(joinpath(BATCH_DIR, "$(isosig).pdf"))`)
end

function quickview(i::Int, ncusps::Int; kwargs...)
    quickview(load(VeeringCensus.lookup(i, ncusps)); kwargs...)
end

function quickview(i::Int; kwargs...)
    quickview(load(VeeringCensus.lookup(i)); kwargs...)
end

function quickview(isosig::String; kwargs...)
    quickview(load(isosig); kwargs...)
end
