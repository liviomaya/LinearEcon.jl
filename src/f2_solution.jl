"""
    LREMSolution(var, eig, flagrank)

Stores the solution to an `LREM` model.

### Arguments

- `var::VAR`: VAR object with the solution. The system incorportes all variables of the LREM, not just the state variables.

- `eig::Vector{<:Number}`: vector with eigenvalues of the `LREM` model.

- `flagrank::Bool`: indicate if rank condition is verified
"""
struct LREMSolution{T<:Number}
    var::VAR
    eig::Vector{T}
    flagrank::Bool
end

noranksolution(eig) = LREMSolution(setvar(0), eig, false)
noranksolution() = noranksolution([NaN])

function solve(lrem::LREM)

    # unload parameters
    A = lrem.A
    B = lrem.B
    C = lrem.C
    Σ = lrem.Σ
    n = lrem.n
    m = lrem.m

    # schur decomposition
    F = schur(B, A)
    eig = F.α ./ F.β # eigenvalues
    Istable = abs.(F.α) .< abs.(F.β)
    # n̄ = count(Istable)
    m̄ = count(.!Istable)
    (m̄ < m) && return noranksolution(eig) # infinite solutions
    (m̄ > m) && return noranksolution(eig) # no stationary solution

    noforward = (m == 0)
    nostate = (n == 0)

    T, S, Q, Z, ~, ~ = ordschur(F, Istable)
    Ct = Q' * C

    if noforward
        (det(Z) == 0) && return noranksolution(eig)
        P = Z * inv(S) * T * inv(Z)
        Q = Z * inv(S) * Ct
    elseif nostate
        M = -inv(T) * Ct
        P = zeros(m, m)
        Q = Z * M
    else
        Z11 = Z[1:n, 1:n]
        (det(Z11) == 0) && return noranksolution(eig)

        Z12 = Z[1:n, n+1:end]
        Z21 = Z[n+1:end, 1:n]
        Z22 = Z[n+1:end, n+1:end]

        S11 = S[1:n, 1:n]
        # S12 = S[1:n, n+1:end]
        # S22 = S[n+1:end, n+1:end]

        T11 = T[1:n, 1:n]
        T12 = T[1:n, n+1:end]
        T22 = T[n+1:end, n+1:end]

        Ct1 = Ct[1:n, :]
        Ct2 = Ct[n+1:end, :]

        M = -inv(T22) * Ct2
        Py = Z21 * inv(Z11)
        Qy = (Z22 - Z21 * inv(Z11) * Z12) * M

        Px = Z11 * inv(S11) * T11 * inv(Z11)
        Qx = Z11 * inv(S11) * ((T12 - T11 * inv(Z11) * Z12) * M + Ct1)

        P = [Px zeros(n, m); Py zeros(m, m)]
        Q = [Qx; Qy]
    end

    flagrank = true

    P, Q = map(x -> real.(x), [P, Q])
    P, Q = map(x -> round.(x, digits=12), [P, Q])
    v = setvar(zeros(n + m), P, Q, Σ)

    sol = LREMSolution(v, eig, flagrank)

    return sol
end