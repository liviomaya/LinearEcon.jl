# promote functions
mat(x::T) where {T<:Real} = fill(float(x), 1, 1)
mat(x::Vector{T}) where {T<:Real} = reshape(float(x), length(x), 1)
mat(x::Matrix{T}) where {T<:Real} = float(x)
mat(x::Diagonal) = Matrix(float(x))
makehermitian(X::Matrix{Float64}) = (X .+ X') / 2 |> Hermitian |> Matrix


"""
    LREM(A,B,C,Σ,n)
    LREM(A,B,C,Σ,n,m,q)

`LREM` stores the parameters that define the linear rational expectations model:

        A[xₜ, Eₜyₜ₊₁] = B [xₜ₋₁, yₜ] + C ϵₜ             ϵₜ ∼ 𝑁(0,Σ)

State variables `x` should be indexed on top.

### Sizes

- `n::Int64`: number of state variables (necessary argument)

- `m::Int64`: number of forward-looking or static variables

- `q::Int64`: number of shocks
"""
struct LREM
    A::Matrix{Float64}
    B::Matrix{Float64}
    C::Matrix{Float64}
    Σ::Matrix{Float64}
    n::Int64 # number of states
    m::Int64 # number of forward-looking or static variables
    q::Int64 # number of shocks
    LREM(A, B, C, Σ, n, m, q) = new(mat(A), mat(B), mat(C), mat(Σ), n, m, q)
    LREM(A, B, C, Σ, n) = new(mat(A), mat(B), mat(C), mat(Σ), n,
        size(A, 2) - n, size(Σ, 1))
end
