# convert array into Array{Float64,2}
function fix_type(x)
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

# build model structure
function model(A,B,C,D,F,Σ,n)
    m = size(A,2) - n
    p = size(A,1) - n - m
    q = size(Σ,1)
    r = size(F,2)
    A,B,C,D,F,Σ = map(fix_type, (A,B,C,D,F,Σ))
    m = model(A,B,C,D,F,Σ,n,m,p,q,r)
    return m
end
# m = model(A,B,C,D,F,Σ,n)

model(A,B,C,D,Σ,n) = model(A,B,C,D,zeros(size(A,1),1),Σ,n)

function addauxiliaries(m0::model)
        n,m,p = m0.n, m0.m, m0.p
        nmp = n + m + p
        nm = n + m
        A, B, C, D, F = m0.A, m0.B, m0.C, m0.D, m0.F
        Iaux = [false for i in 1:nmp]

        # If no static variable, add independent one
        if m0.p == 0
                A = [A; 1 zeros(1,nm-1)]
                B = [B; 0]
                C = [C zeros(nmp); zeros(1,nmp) 1]
                D = [D; zeros(1,m0.q)]
                F = [F; 0]
                Iaux = [Iaux; true]
                p += 1
                nmp += 1
        end

        # If no state variable, include an independent one
        if m0.n == 0
                A = [1 zeros(1,nm); zeros(nmp) A]
                B = [0.0; B]
                C = [0.50 zeros(1,nmp); zeros(nmp) C]
                D = [zeros(1,m0.q); D]
                F = [0; F]
                Iaux = [true; Iaux]
                n += 1
                nm += 1
                nmp += 1
        end

        # If no forward looking variable, include an independent one
        if m0.m == 0
                A = [A[1:nm,:] zeros(nm); zeros(1,nm) 1; A[nm+1:end,:] zeros(p)]
                B = [B[1:nm,1]; 0.0; B[nm+1:end,1]]
                C = [C[1:nm,1:nm] zeros(nm) C[1:nm,nm+1:end]; zeros(1,nm) 1.5 zeros(1,p);
                C[nm+1:end,1:nm] zeros(p) C[nm+1:end,nm+1:end]]
                D = [D[1:nm,:]; zeros(1,m0.q); D[nm+1:end,:]]
                F = [F[1:nm,:]; zeros(1,m0.r); F[nm+1:end,:]]
                Iaux = [Iaux[1:nm,1]; true; Iaux[nm+1:end,1]]
                m += 1
                nm += 1
                nmp += 1
        end
        m1 = model(A,B,C,D,F,m0.Σ,n)
        return m1, Iaux
end

function submatrices(m0::model)
        n,m,p,q,r = m0.n, m0.m, m0.p, m0.q, m0.r
        A11 = m0.A[1:n, 1:n]
        A12 = m0.A[1:n, n+1:end]
        A21 = m0.A[n+1:n+m, 1:n]
        A22 = m0.A[n+1:n+m, n+1:end]
        A31 = m0.A[n+m+1:end, 1:n]
        A32 = m0.A[n+m+1:end, n+1:end]
        B1 = m0.B[1:n,1]
        B2 = m0.B[n+1:n+m,1]
        B3 = m0.B[n+m+1:end,1]
        C11 = m0.C[1:n, 1:n]
        C12 = m0.C[1:n, n+1:n+m]
        C13 = m0.C[1:n, n+m+1:end]
        C21 = m0.C[n+1:n+m, 1:n]
        C22 = m0.C[n+1:n+m, n+1:n+m]
        C23 = m0.C[n+1:n+m, n+m+1:end]
        C31 = m0.C[n+m+1:end, 1:n]
        C32 = m0.C[n+m+1:end, n+1:n+m]
        C33 = m0.C[n+m+1:end, n+m+1:end]
        D1 = m0.D[1:n,:]
        D2 = m0.D[n+1:n+m,:]
        D3 = m0.D[n+m+1:end,:]
        F1 = m0.F[1:n,:]
        F2 = m0.F[n+1:n+m,:]
        F3 = m0.F[n+m+1:end,:]

        return A11,A12,A21,A22,A31,A32, B1,B2,B3, C11,C12,C13,C21,C22,C23,C31,C32,C33, D1,D2,D3, F1,F2,F3
end

function dynmodel(m0::model)

        A11,A12,A21,A22,A31,A32, B1,B2,B3, C11,C12,C13,C21,C22,C23,C31,C32,C33, D1,D2,D3, F1,F2,F3 = submatrices(m0)

        Ah1 = [A11 A12; A21 A22]
        Ah2 = [A31 A32]
        Bh1 = [B1; B2]
        Bh2 = B3
        Ch11 = [C11 C12; C21 C22]
        Ch12 = [C13; C23]
        Ch21 = [C31 C32]
        Ch22 = C33
        Dh1 = [D1; D2]
        Dh2 = D3
        Fh1 = [F1; F2]
        Fh2 = F3

        det(Ch22)==0 && error("Can't express static variables as functions of dynamic variables. Please re-order equations. (Ĉ₂₂ not invertible).")

        Gh = Ah1 - Ch12 * inv(Ch22) * Ah2
        Hh = Bh1 - Ch12 * inv(Ch22) * Bh2
        Jh = Ch11 - Ch12 * inv(Ch22) * Ch21
        Kh = Dh1 - Ch12 * inv(Ch22) * Dh2
        Lh = Fh1 - Ch12 * inv(Ch22) * Fh2

        H = inv(Gh) * Hh
        J = inv(Gh) * Jh
        K = inv(Gh) * Kh * m0.Σ
        L = inv(Gh) * Lh
        H,J,K,L = map(fix_type, (H,J,K,L))
        dm = dynmodel(H,J,K,L,m0.n,m0.m,m0.q,m0.r)
        return dm
end
# dm = dynmodel(m)
