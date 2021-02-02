function dynsolution(H::Array{Float64,2},
        J::Array{Float64,2},
        K::Array{Float64,2},
        L::Array{Float64,2},
        n::Int64)

        H1, H2 = H[1:n, :], H[n+1:end, :]
        K1, K2 = K[1:n, :], K[n+1:end, :]
        L1, L2 = L[1:n, :], L[n+1:end, :]

        m = size(H,1) - n
        Vi = eigen(J).vectors
        Λ_vec = eigen(J).values
        Id = sortperm(abs.(Λ_vec))
        Λ_vec = Λ_vec[Id]
        Vi = Vi[:,Id]
        V = inv(Vi)
        Λ = diagm(0 => Λ_vec)

        n̄ = sum(abs.(Λ_vec).<=1)
        m̄ = sum(abs.(Λ_vec).>1)
        V11 = V[1:n̄,1:n]
        V12 = V[1:n̄,n+1:n+m]
        V21 = V[n̄+1:n̄+m̄,1:n]
        V22 = V[n̄+1:n̄+m̄,n+1:n+m]
        R11 = Vi[1:n,1:n̄]
        R12 = Vi[1:n,n̄+1:n̄+m̄]
        R21 = Vi[n+1:n+m,1:n̄]
        R22 = Vi[n+1:n+m,n̄+1:n̄+m̄]
        Λ1 = Λ[1:n̄,1:n̄]
        Λ2 = Λ[n̄+1:n̄+m̄,n̄+1:n̄+m̄]
        Ht1 = V11 * H1 + V12 * H2
        Ht2 = V21 * H1 + V22 * H2
        Kt1 = V11 * K1 + V12 * K2
        Kt2 = V21 * K1 + V22 * K2
        Lt1 = V11 * L1 + V12 * L2
        Lt2 = V21 * L1 + V22 * L2

        # Check rank conditions
        println("$m̄ unstable equations to $m forward looking variables.")
        if m̄ < m
                error("More forward looking variables than
                unstable eigenvalues: infinite solutions.")
        elseif m̄ > m
                error("More unstable eigenvalues than
                forward looking variables: no solution.")
        end
        if det(V22) == 0
                error("Rank condition not verified. V₂₂ not invertible.")
        end
        println("Rank condition verified.")

        In = diagm(0 => ones(n))
        Im = diagm(0 => ones(m))
        cx = R11*Ht1 - (In - R11*Λ1*inv(R11)) * R12 * inv(Im-inv(Λ2)) * inv(Λ2) * Ht2
        cy = - inv(V22) * inv(Im - inv(Λ2)) * inv(Λ2) * Ht2
        Px = R11 * Λ1 * inv(R11)
        Py = R21 * inv(R11)
        Qx = R11*Kt1 + R11*Λ1*inv(R11)*R12*inv(Λ2)*Kt2
        Qy = -inv(V22)*inv(Λ2)*Kt2
        Mx0 = R11 * Lt1 + R11 * Λ1 * inv(R11) * R12 * inv(Λ2) * Lt2
        My0 = -inv(V22)*inv(Λ2)*Lt2
        αx = R11 * Λ1 * inv(R11) * R12 * inv(Λ2) - R12
        αy = -inv(V22)*inv(Λ2)

        return cx,cy,Px,Py,Qx,Qy,Λ2,Lt2,Mx0,My0,αx,αy
end

dynsolution(dm::dynmodel) = dynsolution(dm.H,
        dm.J, dm.K, dm.L, dm.n)

function solvemodel(m0::model)

        dm = dynmodel(m0)
        cx,cy,Px,Py,Qx,Qy,Λ2,Lt2,Mx0,My0,αx,αy = dynsolution(dm)
        A11,A12,A21,A22,A31,A32, B1,B2,B3, C11,C12,C13,C21,C22,C23,C31,C32,C33, D1,D2,D3, F1,F2,F3 = submatrices(m0)

        auxA = A31 + A32 * Py
        cz = inv(C33) * (auxA*cx + A32 * cy - B3 - C32 * cy)
        Pz = inv(C33) * (auxA*Px - C31 - C32 * Py)
        Qz = inv(C33) * (auxA*Qx - C32 * Qy - D3 * m0.Σ)
        Mz0 = inv(C33) * (auxA*Mx0 - C32*My0 - F3)
        αz =  inv(C33) * (auxA*αx + A32*αy*Λ2 - C32*αy)

        c = [cx; cy; cz]
        P = [Px; Py; Pz]
        Q = [Qx; Qy; Qz]

        flag_complex = false
        c, flag_complex = convert_real(c, flag_complex)
        P, flag_complex = convert_real(P, flag_complex)
        Q, flag_complex = convert_real(Q, flag_complex)
        Lt2, flag_complex = convert_real(Lt2, flag_complex)
        Mx0, flag_complex = convert_real(Mx0, flag_complex)
        My0, flag_complex = convert_real(My0, flag_complex)
        Mz0, flag_complex = convert_real(Mz0, flag_complex)

        flag_complex && println("Imaginary component found in the solution.")

        # the following lines are problematic:
        Λ2 = convert_real(Λ2)
        αx = convert_real(αx)
        αy = convert_real(αy)
        αz = convert_real(αz)

        dp = dpsolution(Λ2,Lt2,Mx0,My0,Mz0,αx,αy,αz)
        sol = solution(c,P,Q,dp)
        return sol
end

function solution(m0::model)
        m1, Iaux = addauxiliaries(m0)
        sol1 = solvemodel(m1)
        # Remove auxiliaries
        c = sol1.c[.!Iaux, 1]
        if m0.n == 0
                P = [0.0, 0.0]
        else
                P = sol1.P[.!Iaux, .!Iaux[1:m1.n]]
        end
        Q = sol1.Q[.!Iaux, :]

        if m0.m == 0
                Mx0 = sol1.dp.Mx0
                My0 = zeros(1,m0.r)
                αx = zeros(m0.n, 1)
                αy = zeros(1, 1) # irrelevant
                Λ2 = 1.5 # irrelevant
                Lt2 = zeros(1, m0.r)
                if m0.p > 0
                        Mz0 = sol1.dp.Mz0
                        αz = zeros(m0.p, 1)
                elseif m0.p == 0
                        Mz0 = zeros(1,m0.r)
                        αz = 0.0
                end
        elseif m0.n == 0
                Mx0 = zeros(1,m0.r)
                My0 = sol1.dp.My0
                αx = zeros(1, m0.m)
                αy = sol1.dp.αy
                Λ2 = sol1.dp.Λ2
                Lt2 = sol1.dp.Lt2
                if m0.p > 0
                        Mz0 = sol1.dp.Mz0
                        αz = sol1.dp.αz
                elseif m0.p == 0
                        Mz0 = zeros(1, m0.r)
                        αz = zeros(1, m0.m)
                end
        else
                Mx0 = sol1.dp.Mx0
                My0 = sol1.dp.My0
                αx = sol1.dp.αx
                αy = sol1.dp.αy
                Λ2 = sol1.dp.Λ2
                Lt2 = sol1.dp.Lt2
                if m0.p > 0
                        Mz0 = sol1.dp.Mz0
                        αz = sol1.dp.αz
                elseif m0.p == 0
                        Mz0 = zeros(1, m0.r)
                        αz = zeros(1, m0.m)
                end
        end

        Λ2,Lt2,Mx0,My0,Mz0,αx,αy,αz,c,P,Q = map(x -> round.(x, digits=8), [Λ2,Lt2,Mx0,My0,Mz0,αx,αy,αz,c,P,Q])
        Λ2,Lt2,Mx0,My0,Mz0,αx,αy,αz,c,P,Q = map(x -> fix_type(x), [Λ2,Lt2,Mx0,My0,Mz0,αx,αy,αz,c,P,Q])

        dp1 = dpsolution(Λ2,Lt2,Mx0,My0,Mz0,αx,αy,αz)
        sol = solution(c,P,Q,dp1)

        return sol
end
