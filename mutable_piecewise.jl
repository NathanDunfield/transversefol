mutable struct MPiecewise <: Homeo
	left::Union{MPiecewise, Linear}
	right::Union{MPiecewise, Linear}
	x::T
	fx::T #x maps to fx
end #todo: Piecewise should just have a flag letting you lazily invert it.

struct InvMPiecewise
    h::MPiecewise
end


inv(h::MPiecewise) = InvMPiecewise(h)
inv(h::Linear) = Linear()
inv(h::InvMPiecewise) = h.h

struct MCandidate
	bt::BoundaryTriangulation
	d::ArrayDict{Junction, Union{InvMPiecewise, MPiecewise,Linear},2}
end

function (f::MPiecewise)(y::T)
	@assert 0 <= y <=1
	if y < f.x
		ret = f.left(y/f.x) * f.fx
	else
		ret = f.right((y-f.x)/(1-f.x)) * (1-f.fx) + f.fx
	end
	@assert 0 <= ret <= 1
	return ret
end

function inv_apply(f::MPiecewise, y::T)
	@assert 0 <= y <=1
	if y < f.fx
		ret = inv_apply(f.left, y/f.fx) * f.x
	else
		ret = inv_apply(f.right, (y-f.fx)/(1-f.fx)) * (1-f.x) + f.x
	end
	@assert 0 <= ret <= 1
	return ret
end

function inv_apply(f::Linear, y::T)
    @assert 0 <= y <= 1
    return y
end

function to_mutable(h::Piecewise)
    return MPiecewise(to_mutable(h.left), to_mutable(h.right), h.x, h.fx)
end

function to_static(h::MPiecewise)
    return Piecewise(to_static(h.left), to_static(h.right), h.x, h.fx)
end

function to_static(h::InvMPiecewise)
    return inv(to_static(h.h))
end

to_static(l::Linear) = Linear()
to_mutable(l::Linear) = Linear()

function (f::InvMPiecewise)(y::T)
    inv_apply(f.h, y)
end

function to_mutable(c::Candidate)
    D=Dict{Junction, Union{InvMPiecewise, MPiecewise, Linear}}(x=>to_mutable(y) for (x,y) in c.d)
    for J in c.bt.junctions
        @assert !J.inv
        D[inv(J)] = inv(D[J])
    end
    MCandidate(c.bt,
               ArrayDict(D)
              )
end

function to_static(c::MCandidate)
    D=Dict{Junction, Union{Piecewise, Linear}}()
    for J in c.bt.junctions
        D[J]=to_static(c.d[J])
    end
    return Candidate(c.bt, D)
end


function jiggle!(p::MPiecewise, r::T)
    jiggle!(p.left,r)
    jiggle!(p.right,r)
    p.x = jiggle(p.x, r)
    p.fx = jiggle(p.fx, r)
    return p
end

function jiggle!(p::InvMPiecewise, r::T)
    jiggle!(p.h, r)
end

function jiggle!(c::MCandidate, r::T)
    for J in c.bt.junctions
        jiggle!(c.d[J],r)
    end 
end

function jiggle!(p::Linear, r::T)
    return p
end

function copy_to!(c1::MCandidate, c2::MCandidate) #copy all the weights of c1 into c2
    for J in c1.bt.junctions
        copy_to!(c1.d[J], c2.d[J])
    end
end

function copy_to!(c1::Linear, c2::Linear)
    return
end

function copy_to!(c1::MPiecewise, c2::MPiecewise)
    c2.x = c1.x
    c2.fx = c1.fx
    copy_to!(c1.left, c2.left)
    copy_to!(c1.right, c2.right)
end
