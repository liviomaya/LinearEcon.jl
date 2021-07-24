
# VAR Model
# xₜ = P xₜ₋₁ + Q ϵₜ, 
# cov(ϵₜ) = Σ
struct VARModel
    P::Array{Float64,2}
    Q::Array{Float64,2}
    Σ::Array{Float64,2}
    n::Int64
    q::Int64
end

VARModel(P,Q,Σ) = VARModel(P, Q, Σ, size(P, 1), size(Σ, 1))
empty_VARModel() = VARModel(zeros(0, 0), zeros(0, 0), zeros(0, 0))

# Linear Rational Expectations Model
# A (xₜ ; Eₜ yₜ₊₁) = B (xₜ₋₁ ; yₜ) + C ϵₜ
# cov(ϵₜ) = Σ
struct Model
    A::Array{Float64,2} 
    B::Array{Float64,2}
    C::Array{Float64,2}
    Σ::Array{Float64,2}
    n::Int64 # number of states
    m::Int64 # number of forward-looking or static variables
    q::Int64 # number of shocks
    S::VARModel # solution in VAR form
    flag_rank::Bool
    flag_complex::Bool
end
# m = Model(A,B,C,Σ,n,m,q,S,flag_rank,flag_complex)

# Model(A,B,C,Σ,n)
