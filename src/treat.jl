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
    m = Model(A,B,C,Σ,n,m,q)
    return m
end
# m = model(A,B,C,Σ,n)

function addauxiliaries(m0::Model)
        n, m, q = m0.n, m0.m, m0.q
        neq = n + m
        A, B, C, Σ = m0.A, m0.B, m0.C, m0.Σ
        Iaux = [false for i in 1:neq]

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
        m1 = Model(A,B,C,Σ,n)
        return m1, Iaux
end
