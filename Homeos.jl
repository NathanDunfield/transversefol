abstract type Homeo end
#represents a piecewise linear function from [0,1] to [0,1]
#
#=
#really, the identity function
struct Linear <: Homeo
end
=#

@enum Dir LEFT=1 RIGHT=2 BOTH=3 #idea: make BOTH take up two slots. Then we never need to change the length of the array.
@enum RoundMode DOWN=1 UP=2
struct DiscreteHomeo <: Homeo #todo: precompute the output heights
    left_heights::Vector{Rational{Int}} #sorted list of heights on the left
    right_heights::Vector{Rational{Int}} #sorted list of heights on the right
    ordering::Vector{Dir} #ordering of the left_heights and right_heights. We should have ordering[1] == ordering[end]==BOTH
    dir::Dir
    roundmode::RoundMode
end

function inv(d::DiscreteHomeo) #uses the same underlying array
    return DiscreteHomeo(d.right_heights, d.left_heights, d.ordering, inv(d.dir), d.roundmode)
end

function inv(d::Dir)
    if d == LEFT
        return RIGHT
    elseif d==RIGHT
        return LEFT
    else
        @assert d==BOTH
        return BOTH
    end
end


function assign_height_left(f::Homeo, r::Real)
    return f(r)
end

function assign_height_right(f::Homeo, r::Real)
    return r
end

function assign_height_left(f::DiscreteHomeo, r::Rational{Int})
    index_range = searchsorted(f.left_heights, r)
    @assert length(index_range)==1
    index = index_range.stop

    count=0
    for k in 1:length(f.ordering)
        if f.ordering[k] == f.dir || f.ordering[k] == BOTH
            count+=1
        end
        if count==index
            return k//length(f.ordering)
        end
    end
    @assert false
end

function assign_height_right(f::DiscreteHomeo, r::Rational{Int})
    index_range = searchsorted(f.right_heights, r)
    @assert length(index_range)==1
    index = index_range.stop

    count=0
    for k in 1:length(f.ordering)
        if f.ordering[k] == inv(f.dir) || f.ordering[k] == BOTH
            count+=1
        end
        if count==index
            return k//length(f.ordering)
        end
    end
    @assert false

end

function (f::DiscreteHomeo)(r::Rational{Int})
    verbose =  false #rand() < 0.01

    index_range = searchsorted(f.left_heights, r)
    if length(index_range) != 1
        @show f
        @show r
    end
    @assert length(index_range)==1

    index = index_range.stop

    if f.roundmode == DOWN
        curr_leftindex = 0
        curr_rightindex = 0

        for k in 1:length(f.ordering)
            if f.ordering[k] == f.dir
                curr_leftindex += 1
            elseif f.ordering[k] == inv(f.dir)
                curr_rightindex += 1
            else 
                @assert f.ordering[k] == BOTH
                curr_leftindex += 1
                curr_rightindex += 1
            end
            if curr_leftindex == index
                return f.right_heights[curr_rightindex]
            end
        end
    elseif f.roundmode == UP
        curr_leftindex = length(f.left_heights)+1
        curr_rightindex = length(f.right_heights)+1

        for k in length(f.ordering):-1:1
            if f.ordering[k] == f.dir
                curr_leftindex -= 1
            elseif f.ordering[k] == inv(f.dir)
                curr_rightindex -= 1
            else
                @assert f.ordering[k] == BOTH
                curr_leftindex-=1
                curr_rightindex-=1
            end
            if curr_leftindex == index
                return f.right_heights[curr_rightindex]
            end
        end
    end 
    @show f.roundmode
    @show f.ordering
    @show curr_leftindex
    @show curr_rightindex
    @assert false
end


function test_discrete_homeo()
    thickness=5
    left_heights = Rational{Int}[i//(thickness) for i in 0:thickness]
    right_heights = Rational{Int}[i//(thickness) for i in 0:thickness]
    ordering = shuffle(vcat(Dir[LEFT for i in 1:length(left_heights)-2], Dir[RIGHT for i in 1:length(right_heights)-2]))
    @show left_heights
    @show right_heights
    @show ordering
    f = DiscreteHomeo(left_heights, right_heights, ordering, LEFT, UP)

    @show [f(r) for r in left_heights]
end
