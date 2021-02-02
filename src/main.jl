function solvemodel(m0::model)

        #dm = dynmodel(m0)
        #cx,cy,Px,Py,Qx,Qy,Λ2,Lt2,Mx0,My0,αx,αy = dynsolution(dm)
        #A11,A12,A21,A22,A31,A32, B1,B2,B3, C11,C12,C13,C21,C22,C23,C31,C32,C33, D1,D2,D3, F1,F2,F3 = submatrices(m0)

        A, B, C, Φ, Ω, Σ = m0.A, m0.B, m0.C, m0.Φ, m0.Ω, m0.Σ
        n, m, p, q = m0.n, m0.m, m0.p, m0.q
        
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

        vecM = inv( ( kron(Φ', S22) - kron(collect(I(p)), T22) )  ) * Ct2[:]
        M = reshape(vecM, m, p)

        Py = Z21 * inv(Z11)
        Qy = (Z22 - Z21*inv(Z11)*Z12) * M

        Px = Z11*inv(S11)*T11*inv(Z11)
        Qx = (Z12 - Z11*inv(S11)*S12)*M*Φ + (Z11*inv(S11)*T12 - Z11*inv(S11)*T11*inv(Z11)*Z12)*M + Z11*inv(S11)*Ct1

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
