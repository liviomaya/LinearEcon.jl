struct model
    A::Array{Float64,2}
    B::Array{Float64,2}
    C::Array{Float64,2}
    D::Array{Float64,2}
    F::Array{Float64,2}
    Σ::Array{Float64,2}
    n::Int64
    m::Int64
    p::Int64
    q::Int64
    r::Int64
end
# m = model(A,B,C,D,F,Σ,n,m,p,q,r)

struct dynmodel
    H::Array{Float64,2}
    J::Array{Float64,2}
    K::Array{Float64,2}
    L::Array{Float64,2}
    n::Int64
    m::Int64
    q::Int64
    r::Int64
end
# dm = dynmodel(H,J,K,L,n,m,q,r)

struct dpsolution
    Λ2::Array{Float64,2}
    Lt2::Array{Float64,2}
    Mx0::Array{Float64,2}
    My0::Array{Float64,2}
    Mz0::Array{Float64,2}
    αx::Array{Float64,2}
    αy::Array{Float64,2}
    αz::Array{Float64,2}
end

struct solution
    c::Array{Float64,2}
    P::Array{Float64,2}
    Q::Array{Float64,2}
    dp::dpsolution
end
# sol = solution(c,P,Q,dp)
