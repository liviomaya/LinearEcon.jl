
"""
    M = moving_average(S, T)

Compute the moving average representation of VAR model `xₜ = P xₜ₋₁ + Q ϵₜ` stored in `S`, going `T` periods back. The MA representation is `xₜ = ∑ᵢ₌₀ Cᵢ ϵₜ₋ᵢ` with `Cᵢ = Pⁱ Q`. Output `M[:,:,i]` stores `Cᵢ₋₁`.
"""
function moving_average(S::VARModel, T)
    MA = zeros(S.n, S.q, T)
    for t in 1:T
        MA[:, :, t] = (S.P^(t - 1)) * S.Q
    end
    return MA
end

"""
    V = var_decomp(S, T)

Compute the variance decomposition of `T` period ahead forecast errors for VAR model `S`. Output `V[i,j]` stores `i`-th variable variance due to shock `j`.
"""
    function var_decomp(S::VARModel, T)
    @assert isdiag(S.Σ)
    MA = moving_average(S, T)
    V = zeros(S.n, S.q)
    for j in 1:S.q
        Ij = zeros(S.q, S.q)
        Ij[j,j] = 1
        V[:,j] = diag(sum(MA[:,:,t] * Ij * S.Σ * MA[:,:,t]' for t in 1:T))
    end
    return V
end


function path_plot_template(X, label)
    T = size(X, 1)

    Fig = plot(xlim=(1, Inf), gridalpha=0.05, tickfontsize=12)

    plot!(Fig, 1:T, X, markershape=:circle, lw=2, ms=6, label=label)

    plot!(Fig, 1:T, zeros(T), color=:black, alpha=0.15, label=:none)

    return Fig
end

"""
    Fig = irf(S, e)
    Fig = irf(S, e; T, Vars, label)

Returns a pre-formatted figure `Fig` with the impulse response function (the MA presentation) of VAR model `S` to the shock indexed by `e`.

# Keyword Arguments
`T`: number of periods. Default = 25

`id`: array with indices of variables to be displayed. Default = All variables

`label`: figure labels. Default = :none
"""
    function irf(S::VARModel, e::Int64; T=25, id=1:S.n, label=:none)
    MA = moving_average(S, T)
    X = MA[id,e,:]'
    Fig = path_plot_template(X, label)
    return Fig
end


"""
    Fig = simulate(S)
    Fig = simulate(S; T, id, label, X0)

Returns a pre-formatted figure `Fig` with a simulated path of the VAR model `S`.

# Keyword Arguments
`T`: number of periods. Default = 25

`id`: array with indices of variables to be displayed. Default = All variables

`label`: figure labels. Default = `nothing`

`X0`: initial state of the system. Default = vector of zeros
"""
function simulate(S::VARModel; T=25, id=1:S.n, label=nothing, X0=zeros(S.n))
    X = zeros(S.n, T)
    ShockDist = MvNormal(zeros(S.q), S.Σ)
    for t in 1:T
        State = (t == 1) ? X0 : X[:,t - 1]
        X[:,t] = S.P * State + S.Q * rand(ShockDist)
    end
    Y = X[id, :]'
    Fig = path_plot_template(Y, label)
    return Fig
end


"""
    V, FlagConverged = cov(S)
    V, FlagConverged = cov(S; IterMax, Tol)

Compute the covariance matrix `V` of the VAR model `S`. Boolean `FlagConverged` indicates if convergence of Lyapunov operator was achieved.

# Keyword Arguments
`IterMax`: maximum number of iterations. Default = 10000.

`Tol`: tolerated approximation error. Default = 1e-8.
"""
    function cov(S::VARModel; IterMax=10000, Tol=1e-8)
    damp = 1.0
    FlagConverged = false
        V = zeros(S.n, S.n)
    for i in 1:IterMax
        TV = S.P * V * S.P' + S.Q * S.Σ * S.Q'
        Dist = maximum(abs.(TV .- V))
        FlagConverged = (Dist < Tol)
        FlagConverged && break
        V = damp * TV + (1 - damp) * V
    end
    V = round.(V, digits=10)
    V = (V .+ V') / 2
    return V, FlagConverged
end

"""
    R, FlagConverged = corr(S)
    R, FlagConverged = corr(S; IterMax, Tol)

Compute the correlation matrix `R` of the VAR model `S`. Boolean `FlagConverged` indicates if convergence of Lyapunov operator was achieved.

# Keyword Arguments
`IterMax`: maximum number of iterations. Default = 10000.

`Tol`: tolerated approximation error. Default = 1e-8.
"""
function corr(S::VARModel; IterMax=10000, Tol=1e-8)
    V, FlagConverged = cov(S, IterMax=IterMax, Tol=Tol)
    SD = sqrt.(diag(V))
    SDM = diagm(0 => SD)
    iSDM = inv(SDM)
    iSDM[SDM .== 0] .= 0.0 # set correlation = 0 to constant variables
R = iSDM * V * iSDM
return R, FlagConverged
end


"""
    M, R = cholesky(S::VarModel)

Orthogonalize innovations of system `S` and return new VAR model `M`. Each column of array `R` contains the linear combination of the original system's IRFs that lead to system `M` IRFs. 

`S` model: xₜ = P xₜ₋₁ + Q ϵₜ       cov(ϵ) = Σ

`M` model: xₜ = P xₜ₋₁ + Q*R ηₜ     cov(η) = I
"""
function cholesky(S::VARModel)
    # ηₜ = R⁻¹ ϵₜ       cov(η) = I
    G = Matrix(cholesky(S.Σ).L) # Σ = G G'
    # I = cov(η) = R⁻¹ Σ R⁻¹' = R⁻¹ G G' R⁻¹' => R = G
    R = G 
    # MA: xₜ = ∑ Cᵢ ϵₜ₋ᵢ = ∑ Cᵢ R R⁻¹ ϵₜ₋ᵢ = ∑ Cᵢ R ηₜ₋ᵢ
    # State Space : xₜ = P xₜ₋₁ + Q ϵₜ = P xₜ₋₁ + (Q R) ηₜ
    NewQ = S.Q * R
    Id = Matrix(I(S.q))
    M = VARModel(S.P, NewQ, Id)
    return M, R
end
