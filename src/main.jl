function solvemodel(m0::model)
        A, B, C, Σ = m0.A, m0.B, m0.C, m0.Σ
        n, m, q = m0.n, m0.m, m0.q
        
        F = schur(B, A)
        Istable = abs.(F.α) .< abs.(F.β)
        n̄ = count(Istable)
        m̄ = count(.!Istable)

        # Check eigenvalues
        println("$n̄ stable eigenvalues to $n state variables.")
        if m̄ < m
                error("More forward looking variables than
                unstable eigenvalues: infinite solutions.")
        elseif m̄ > m
                error("More unstable eigenvalues than
                forward looking variables: no solution.")
        end

        T, S, Q, Z, t, s = ordschur(F, Istable)

        Z11 = Z[1:n,1:n]
        if det(Z11) == 0
                error("Rank condition not verified. Z₁₁ is singular.")
        end
        println("Rank condition verified.")

        Z12 = Z[1:n,n+1:end]
        Z21 = Z[n+1:end,1:n]
        Z22 = Z[n+1:end,n+1:end]
        
        S11 = S[1:n,1:n]
        S12 = S[1:n,n+1:end]
        S22 = S[n+1:end,n+1:end]

        T11 = T[1:n,1:n]
        T12 = T[1:n,n+1:end]
        T22 = T[n+1:end,n+1:end]

        Ct = Q'*C
        Ct1 = Ct[1:n,:]
        Ct2 = Ct[n+1:end,:]

        M = -inv(T22)*Ct2

        Py = Z21 * inv(Z11)
        Qy = (Z22 - Z21*inv(Z11)*Z12) * M

        Px = Z11*inv(S11)*T11*inv(Z11)
        Qx = Z11*inv(S11)*((T12 - T11*inv(Z11)*Z12)*M + Ct1)

        P = [Px; Py]
        Q = [Qx; Qy]

        flag_complex = false
        P, flag_complex = convert_real(P, flag_complex)
        Q, flag_complex = convert_real(Q, flag_complex)
        flag_complex && println("Imaginary component found in the solution.")

        sol = solution(P,Q)
        return sol
end

function solution(m0::model)
        m1, Iaux = addauxiliaries(m0)
        sol1 = solvemodel(m1)
        # Remove auxiliaries
        if m0.n == 0
                P = [0.0]
        else
                P = sol1.P[.!Iaux, .!Iaux[1:m1.n]]
        end
        Q = sol1.Q[.!Iaux, :]

        P,Q = map(x -> round.(x, digits=10), [P,Q])
        P,Q = map(x -> fix_type(x), [P,Q])
        sol = solution(P,Q)
        return sol
end
