function ConvertArray2(x)
        if typeof(x) in (Int64,Float64)
            return reshape([Float64(x)],1,1)
        elseif typeof(x) in (Array{Int64,1},Array{Float64,1})
            n = length(x)
            return reshape(Float64.(x),n,1)
        elseif typeof(x) in (Array{Int64,2},Array{Float64,2})
            return Float64.(x)
        else
            error("Type error.")
        end
end

function ConvertReal(x, flag)
        flag = flag || maximum(abs.(imag.(x))) .> 1e-5
        x = real.(x)
        return x, flag
end
ConvertReal(x) = ConvertReal(x, true)[1]


function SolutionFailedRank(λ)
        S = EmptyVARModel()
        flag_rank = false
        flag_complex = false
        λ, flag_complex = ConvertReal(λ, flag_complex)
        return λ, S, flag_rank, flag_complex
end

function SolveModel(A, B, C, Σ, n, m)
        
        F = schur(B, A)
        λ = F.α ./ F.β
        Istable = abs.(λ) .< 1
        n̄ = count(Istable)
        m̄ = count(.!Istable)
        (m̄ < m) && return SolutionFailedRank(λ) # infinite solutions
        (m̄ > m) && return SolutionFailedRank(λ) # no stationary solution

        noforward = (m == 0)
        nostate = (n == 0)

        T, S, Q, Z, ~, ~ = ordschur(F, Istable)
        Ct = Q' * C

        if noforward
                (det(Z) == 0) && return SolutionFailedRank(λ)
                P = Z * inv(S) * T * inv(Z)
                Q = Z * inv(S) * Ct
        elseif nostate
                M = -inv(T) * Ct
                P = zeros(0,0)
                Q = Z * M
        else
                Z11 = Z[1:n,1:n]
                (det(Z11) == 0) && return SolutionFailedRank(λ)

                Z12 = Z[1:n,n+1:end]
                Z21 = Z[n+1:end,1:n]
                Z22 = Z[n+1:end,n+1:end]
        
                S11 = S[1:n,1:n]
                S12 = S[1:n,n+1:end]
                S22 = S[n+1:end,n+1:end]

                T11 = T[1:n,1:n]
                T12 = T[1:n,n+1:end]
                T22 = T[n+1:end,n+1:end]

                Ct1 = Ct[1:n,:]
                Ct2 = Ct[n+1:end,:]

                M = -inv(T22)*Ct2
                Py = Z21 * inv(Z11)
                Qy = (Z22 - Z21*inv(Z11)*Z12) * M

                Px = Z11*inv(S11)*T11*inv(Z11)
                Qx = Z11*inv(S11)*((T12 - T11*inv(Z11)*Z12)*M + Ct1)

                P = [Px zeros(n,m); Py zeros(m, m)]
                Q = [Qx; Qy]
        end

        flag_rank = true

        P, Q = map(x -> round.(x, digits=10), [P, Q])
        S = VARModel(P, Q, Σ)

        flag_complex = false
        P, flag_complex = ConvertReal(P, flag_complex)
        Q, flag_complex = ConvertReal(Q, flag_complex)

        return λ, S, flag_rank, flag_complex
end

function model(A, B, C, Σ, n)
        m = size(A,1) - n
        q = size(Σ,1)
        A, B, C, Σ = map(ConvertArray2, (A, B, C, Σ))
        λ, S, flag_rank, flag_complex = SolveModel(A, B, C, Σ, n, m)
        m = Model(A, B, C, Σ, n, m, q, λ, S, flag_rank, flag_complex)
end
