import Base: *
using LinearAlgebra: I
using Random

# --- Random helpers for testing BasisChange ---

function rand_sl2()
    A = Matrix{Int}(I, 2, 2)
    for _ in 1:4
        n = rand(-3:3)
        if rand(Bool)
            A = [1 n; 0 1] * A
        else
            A = [1 0; n 1] * A
        end
    end
    if rand(Bool)
        A = [-1 0; 0 -1] * A
    end
    return A
end

function rand_basis_change(n_in::Int, n_out::Int)
    @assert n_out <= n_in
    surviving = sort(shuffle(1:n_in)[1:n_out])
    slice = Vector{Tuple{Int,Int}}(undef, n_in)
    for j in 1:n_in
        if j in surviving
            slice[j] = (0, 0)
        else
            p, q = rand(1:4), rand(-3:3)
            slice[j] = (p, q)
        end
    end
    bc = [rand_sl2() for _ in 1:n_out]
    perm = shuffle(1:n_out)
    return BasisChange(slice, bc, perm)
end

function rand_valid_input(B::BasisChange)
    s = Vector{Tuple{Int,Int}}(undef, length(B.slice))
    for (j, x) in enumerate(B.slice)
        if x == (0,0)
            p, q = rand(1:4), rand(-3:3)
            s[j] = (p, q)
        else
            s[j] = x
        end
    end
    return s
end

function test_basis_change(; n_trials=100000)
    all_passed = true
    for _ in 1:n_trials
        n2    = rand(2:5)
        n_mid = rand(1:n2)
        n_out = rand(1:n_mid)

        B2  = rand_basis_change(n2, n_mid)
        B1  = rand_basis_change(n_mid, n_out)
        B12 = B1 * B2

        s   = rand_valid_input(B2)
        lhs = B12(s)
        rhs = B1(B2(s))

        if lhs != rhs
            @warn "Associativity failed" B1 B2 s lhs rhs
            all_passed = false
        end
    end

    for _ in 1:n_trials
        n    = rand(2:5)
        B    = rand_basis_change(n, n)
        Binv = inv(B)

        s         = rand_valid_input(B)
        recovered = Binv(B(s))

        if recovered != s
            @warn "inv(B)*B*s != s" B s recovered
            all_passed = false
        end
    end

    if all_passed
        println("All tests passed.")
    end
    return all_passed
end
