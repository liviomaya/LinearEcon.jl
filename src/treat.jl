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
function model(A,B,C,Σ,n)
    m = size(A,1) - n
    q = size(Σ,1)
    A,B,C,Σ = map(fix_type, (A,B,C,Σ))
    m = model(A,B,C,Σ,n,m,q)
    return m
end
# m = model(A,B,C,Σ,n)

function addauxiliaries(m0::model)
        n, m, q = m0.n, m0.m, m0.q
        neq = n + m
        A, B, C, Σ = m0.A, m0.B, m0.C, m0.Σ
        Iaux = [false for i in 1:nmp]

        # If no state variable, include an independent one
        if n == 0
            A = [1 zeros(1,neq); zeros(neq) A]
            B = [0.20 zeros(1,neq); zeros(neq) B]
            C = [zeros(1, q); C]
            Iaux = [true; Iaux]
            n += 1
            neq += 1
        end

        # If no forward looking / static variable, include an independent one
        if m == 0
            A = [A zeros(neq); zeros(1,neq) 1]
            B = [B zeros(neq); zeros(1,neq) 1.5]
            C = [C; zeros(1, q)]
            Iaux = [Iaux; true]
            m += 1
            neq += 1
        end
        m1 = model(A,B,C,Σ,n)
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
