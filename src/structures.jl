struct model
    A::Array{Float64,2}
    B::Array{Float64,2}
    C::Array{Float64,2}
    Σ::Array{Float64,2}
    n::Int64
    m::Int64
    q::Int64
end
# m = model(A,B,C,Σ,n,m,q)
struct solution
    P::Array{Float64,2}
    Q::Array{Float64,2}
    flag_rank::Bool
    flag_complex::Bool
end
# sol = solution(P, Q, flag_rank, flag_complex)
