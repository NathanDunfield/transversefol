using Base.Threads
using VeeringCensus
using Random

N=11031
r = collect(0:100:N)
push!(r,N)
@threads for (x,y) in shuffle([(r[i],r[i+1]) for i in 1:(length(r)-1)])
    try
        run(`/home/jonathan/miniconda3/envs/sage/bin/python3 enumerate_pA.py $(x) $(y)`)
    catch e
        print(e)
    end 
end
