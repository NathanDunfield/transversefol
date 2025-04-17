
function prunings(p::Piecewise)
	chnl = Channel{Homeo}(3)
	put!(chnl, Linear())
	@async begin
		for i1 in prunings(p.left)
			push!(chnl, Piecewise(i1,p.right,p.x,p.fx))
		end
		for i2 in prunings(p.right)
			push!(chnl, Piecewise(p.left,i2,p.x,p.fx))
		end
		close(chnl)
	end
	return chnl
end

function complexity(p::Piecewise)
	return complexity(p.left) + complexity(p.right)
end
function complexity(p::Linear)
	return 1
end

function prunings(l::Linear)
	return (Linear(),)
end

function inv(p::Piecewise)
	Piecewise(inv(p.left),inv(p.right),p.fx,p.x)
end

function inv(x::Linear)
	Linear()
end

function (f::Piecewise)(y::T)
	@assert 0 <= y <=1
	if y < f.x
		ret = f.left(y/f.x) * f.fx
	else
		ret = f.right((y-f.x)/(1-f.x)) * (1-f.fx) + f.fx
	end
	@assert 0 <= ret <= 1
	return ret
end

function (f::Linear)(y::T)
	@assert 0 <= y <= 1
	return y
end

function prunings(c::Candidate)
	ch = Channel{Candidate}(3)
	D=Dict{Junction, Union{Piecewise,Linear}}()
	for J in c.bt.junctions
		D[J] = c[J]
	end
	@async begin
		for J in c.bt.junctions
			for p in prunings(D[J])
				if complexity(p) < complexity(D[J])
					cnew=Candidate(c.bt,Dict{Junction,Union{Piecewise,Linear}}())
					for _J in c.bt.junctions
						cnew[_J]=c[_J]
					end
					cnew[J]=p
					put!(ch, cnew)
				end
			end
		end
		close(ch)
	end
	return ch
end

function complexity(c::Candidate)
	return sum(complexity(c[J]) for J in c.bt.junctions)
end

function prune(c::Candidate)
	val = approximant_all_slopes(c)
	old_complexity = complexity(c)

	@label here
	for cnew in prunings(c)
		if approximant_all_slopes(cnew,time=10000) == val
			c=cnew
			@goto here
		end
	end

	@show (old_complexity, complexity(c))
	return c
end

function prune(E::Envelope{S}) where {S}
    E2 = Envelope{S}(copy(E.A))
    prune!(E2)
    return E2
end

function prune!(E::Envelope)
	@threads for i in eachindex(E.A)
		val,c = E.A[i]
		E.A[i] = (val, prune(c))
	end
    return E
end


function piecewise(inputs::Vector, outputs::Vector)
	@assert length(inputs)==length(outputs)
	for x in inputs
		@assert 0 <= x <= 1
	end
	for x in outputs 
		@assert 0 <= x <= 1
	end
	@assert inputs == sort(inputs)
	@assert outputs == sort(outputs)

	if length(inputs)==0
		return Linear()
	else
		i = round(Int,ceil(length(inputs)/2))
		ret = Piecewise(piecewise(inputs[1:i-1] ./ inputs[i], outputs[1:i-1] ./ outputs[i]),
						 piecewise((inputs[i+1:end] .- inputs[i]) ./ (1-inputs[i]), (outputs[i+1:end] .- outputs[i]) ./ (1-outputs[i])),
						 inputs[i],
						 outputs[i])
		for (i,j) in zip(inputs, outputs)
			@assert abs(j-ret(i))<0.0001
		end
		return ret
	end
end


function random_candidate(bt,depth)
	c=Candidate(bt, Dict{Junction, Union{Linear,Piecewise}}())
	for j in bt.junctions
		#@show depth + round(Int, log2(j.left_len + j.right_len))
		c.d[j]=random_piecewise(depth + round(Int, log2(j.left_len + j.right_len)))
	end
	return c
end
function subdivide(p::Piecewise)
	return Piecewise(subdivide(p.left), subdivide(p.right), p.x, p.fx)
end
function subdivide(p::Linear)
	return Piecewise(Linear(), Linear(), 0.5, 0.5)
end
function subdivide(c::Candidate)
	return Candidate(c.bt, Dict{Junction, Union{Piecewise,Linear}}(x=>subdivide(y) for (x,y) in c.d))
end
function random_piecewise(depth)
	if depth<=0
		Linear()
	else
		Piecewise(random_piecewise(depth-1), random_piecewise(depth-1), uniform(0.01,0.99), uniform(0.01,0.99))
	end
end
function uniform(x,y)
	@assert y>=x
	return rand()*(y-x) + x
end

struct Piecewise <: Homeo
	left::Union{Piecewise, Linear}
	right::Union{Piecewise, Linear}
	x::T
	fx::T #x maps to fx
end #todo: Piecewise should just have a flag letting you lazily invert it.
