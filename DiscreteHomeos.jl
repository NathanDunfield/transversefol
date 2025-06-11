struct DiscreteHomeo{T} <: Homeo #todo: precompute the output heights
    ordering::Vector{2,T}
    dir::Int #1 or 2
    roundmode::RoundMode
end

function (f::DiscreteHomeo{T})(r::T) where {T}
    r = searchsorted(f.ordering[f.dir,:])
    @assert length(r) > 0
    if f.roundmode == UP
        return f.ordering[2-f.dir,r[end]]
    else
        return f.ordering[2-f.dir,r[1]]
    end
end

function inv(f::DiscreteHomeo)
    return DiscreteHomeo(f.ordering, 2-f.dir, f.roundmode)
end

function inject(s::State)

end
