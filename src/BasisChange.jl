# BasisChange: slice, permute, and apply per-coordinate SL(2,Z) transforms.
# Included into the Envelopes module.
#
# Slice occurs first (filled cusps removed), then basis_change applied in
# survival order, then perm reorders to output coordinates.
# Convention: perm[k] = output position for the k-th surviving coordinate.
struct BasisChange
    slice        :: Vector{Tuple{Int,Int}}   # (0,0) = survives; nonzero = filling slope
    basis_change :: Vector{Matrix{Int}}       # one 2x2 matrix per surviving coordinate
    perm         :: Vector{Int}               # perm[k] = output position of surviving coord k
end

function _int_inv(A::Matrix{Int})
    d = A[1,1]*A[2,2] - A[1,2]*A[2,1]
    @assert abs(d) == 1
    return Int[A[2,2] -A[1,2]; -A[2,1] A[1,1]] ./ d
end

function (B::BasisChange)(s::Nothing)
    return nothing
end

function (B::BasisChange)(s::Vector{Rational{Int}})
    tups = Tuple{Int,Int}[(denominator(x), numerator(x)) for x in s]
    result = B(tups)
    result === nothing && return nothing
    return Rational{Int}[r[2]//r[1] for r in result]
end

function (B::BasisChange)(s::Vector{Tuple{Int,Int}})
    @assert length(s) == length(B.slice)
    sliced = Tuple{Int,Int}[]
    for (x, y) in zip(B.slice, s)
        if x == (0,0)
            push!(sliced, y)
        elseif x[1]*y[2] - x[2]*y[1] != 0
            return nothing
        else
            @debug "nontrivial slice"
        end
    end
    # output[perm[k]] = basis_change[k] * sliced[k]
    # i.e. output[l] = basis_change[invperm[l]] * sliced[invperm[l]]
    ip = invperm(B.perm)
    return Tuple{Int,Int}[tuple(B.basis_change[ip[l]] * Int[sliced[ip[l]][1], sliced[ip[l]][2]]...) for l in 1:length(B.perm)]
end

function *(B::BasisChange, s::Vector{Tuple{Int,Int}})
    return B(s)
end

function *(B::BasisChange, s::Vector{Rational{Int}})
    return B(s)
end

function inv(B::BasisChange)
    @assert all(s == (0,0) for s in B.slice) "inv requires no sliced dimensions"
    inv_perm = invperm(B.perm)
    inv_basis = [_int_inv(B.basis_change[inv_perm[j]]) for j in 1:length(inv_perm)]
    return BasisChange(B.slice, inv_basis, inv_perm)
end

# Apply B2 first, then B1.
function *(B1::BasisChange, B2::BasisChange)
    n2    = length(B2.slice)
    n_mid = length(B2.perm)
    @assert length(B1.slice) == n_mid

    unfilled2 = [j for j in 1:n2 if B2.slice[j] == (0,0)]
    @assert length(unfilled2) == n_mid

    unfilled1 = [i for i in 1:n_mid if B1.slice[i] == (0,0)]

    comp_slice = Vector{Tuple{Int,Int}}(undef, n2)
    for j in 1:n2
        if B2.slice[j] != (0,0)
            comp_slice[j] = B2.slice[j]
        else
            k = findfirst(==(j), unfilled2)
            i = B2.perm[k]               # new convention: surviving k -> output i
            if B1.slice[i] == (0,0)
                comp_slice[j] = (0,0)
            else
                v = round.(Int, _int_inv(B2.basis_change[k]) * collect(B1.slice[i]))
                comp_slice[j] = (v[1], v[2])
            end
        end
    end

    comp_unfilled = [j for j in 1:n2 if comp_slice[j] == (0,0)]

    n_comp_unfilled = length(comp_unfilled)
    comp_perm         = Vector{Int}(undef, n_comp_unfilled)
    comp_basis_change = Vector{Matrix{Int}}(undef, n_comp_unfilled)

    for (p, j) in enumerate(comp_unfilled)
        k = findfirst(==(j), unfilled2)
        i = B2.perm[k]                   # intermediate position
        m = findfirst(==(i), unfilled1)
        comp_perm[p] = B1.perm[m]
        comp_basis_change[p] = B1.basis_change[m] * B2.basis_change[k]
    end

    return BasisChange(comp_slice, comp_basis_change, comp_perm)
end

# Apply a BasisChange to an (Elower, Eupper) pair.
# Slices at the filling slopes, permutes coordinates, then calls basis_change.
function *(B::BasisChange, p::Tuple{Envelope{Lower,T,D}, Envelope{Upper,T,D}}) where {T<:Rational, D}
    Elower, Eupper = p

    Elower_s = slice(Elower, B.slice)
    Eupper_s = slice(Eupper, B.slice)

    n_out = length(B.perm)
    ip = invperm(B.perm)

    # Reorder surviving coords to output positions: output[l] = surviving[invperm[l]]
    Elower_p = typeof(Elower_s)([(v[ip], c) for (v, c) in Elower_s.A], SpinLock())
    Eupper_p = typeof(Eupper_s)([(v[ip], c) for (v, c) in Eupper_s.A], SpinLock())
    transforms = Matrix{Int}[B.basis_change[ip[i]] for i in 1:n_out]

    return basis_change(Elower_p, Eupper_p, transforms)
end
