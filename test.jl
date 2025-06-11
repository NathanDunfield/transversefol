using Revise
using JSON
includet("script.jl")

@show JSON.parse("batch/$(VeeringCensus.lookup(1,2)).json")
