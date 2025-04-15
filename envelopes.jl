#defining crevices
#order all points by x,y, and z coordinate.
#a crevice has the property that it is not in the interior of envelope, but if you move slightly down, left, or right, you enter the envelope. So it is a triple of points in the envelope, p(_x, p_y, p_z, such that 
#

#We'll start with the points [Inf, -Inf, -Inf], [-Inf, -Inf, Inf], [-Inf, Inf, -Inf].
#We'll maintain a tree-like structure, where each node is a triangle
#Whenever we insert a point, it restricts both the 
#

mutable struct Crevice{N} #can also be thought of as an octant in space
    faces::MVector{N,SVector{N,T}} #should be the points giving rise to this crevice
    pivot::Union{SVector{N,T},Nothing}
    children::Vector{Crevice{N}}
end
function Crevice(faces::AbstractVector{SVector{N,T}}) where {N}
    return Crevice{N}(MVector{N,SVector{N,T}}(faces), nothing, Crevice{N}[])
end

function contains(c::Crevice{N}, p::SVector{N,T}) where {N}
    all(c.faces[i][i] < p[i] for i in 1:N)
end

function is_valid_crevice(c::Crevice{N}) where {N}
    return all( (i==j || c.faces[i][i] < c.faces[j][i]+0.000001)
               for i in 1:N, j in 1:N
              )
end

function push!(c::Crevice{N}, pt::SVector{N,T}) where {N}
    if contains(c,pt)
        if c.pivot == nothing
            c.pivot=pt
            for i in 1:N
                childfaces = copy(c.faces)
                childfaces[i] = pt
                child = Crevice(childfaces)
                if is_valid_crevice(child)
                    push!(c.children, child)
                end
            end 
        else
            for child in c.children
                push!(child, pt)
            end
        end
    end 
end

function leaves(c::Crevice{N}) where {N}
    if c.pivot == nothing
        return [SVector{N}([c.faces[i][i] for i in 1:N])]
    else
        return Iterators.flatten(map(leaves, c.children))
    end
end

function all_leaves(c::Crevice{N}) where {N}
    return unique(collect(leaves(c)))
end

function crevices_general(e::Envelope{Upper})
    N=length(e.A[1][1])
    c=Crevice([SVector{N,T}([(i==j ? -CLIP : CLIP) for j in 1:N]) for i in 1:N])
    for (x,_) in e.A
        push!(c, SVector{N,T}(x))
    end
    return [Vector{T}(x) for x in all_leaves(c)]
end

function crevices_general(e::Envelope{Lower})
    N=length(e.A[1][1])
    c=Crevice([SVector{N,T}([(i==j ? -CLIP : CLIP) for j in 1:N]) for i in 1:N])
    for (x,_) in e.A
        push!(c, -SVector{N,T}(x))
    end
    return [Vector{T}(-x) for x in all_leaves(c)]
end
