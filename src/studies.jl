"""
    x = ss(m, sol)

Compute the steady state for the model `m` with solution `sol`.
"""
function ss(m::Model, sol::Solution)
    if m.n == 0
        SS = sol.c
    else
        Px = sol.P[1:m.n,1:m.n]
        eig = eigen(Px).values
        if maximum(abs.(eig)) < 1
            Id = diagm(0 => ones(m.n))
            xbar = inv(Id - sol.P[1:m.n,1:m.n]) * sol.c[1:m.n,1]
            #SS = sol.c + sol.P * xbar
        else
            Vi = eigen(Px).vectors
            V = inv(Vi)
            Id1 = abs.(eig) .≈ 1
            Vc = V * sol.c[1:m.n]
            if !all(Vc[Id1] .≈ 0)
                error("Unit root found. No stationary steady state.")
            end
            x̂s = Vc[.!Id1] ./ (1 .- eig[.!Id1])

            x̂ = zeros(m.n)
            x̂[.!Id1] = x̂s
            xbar = Vi * x̂
        end
        SS = sol.c + sol.P * xbar
    end
    SS = round.(SS, digits=7)
    return SS
end

"""
    V, R = covariance(m, sol)

Compute the covariance and correlation matrices `V` and `R`, respectivelly, of the model `m` with solution `sol`.

# Keyword Arguments
`itermax`: maximum number of iterations. Default = 10000.

`guess`: `n+m` × `n+m` matrix to be initial guess for the fixed point.
"""
function covariance(m::Model, sol::Solution; itermax = 10000, guess = zeros(m.n + m.m, m.n + m.m))
    !sol.flag_rank && error("Rank condition not verified")

    Σ = m.Σ
    P, Q = sol.P, sol.Q
    if m.n == 0
        cov = Q * Σ * Q'
    else
        damp = 1.0
        converged = false
        cov = guess
        if m.m > 0
            P̄ = [P zeros(m.m+m.n,m.m)]
        else
            P̄ = P
        end
        for i in 1:itermax
            covUpd = P̄ * cov * P̄' + Q * Σ * Q'
            dist = maximum(abs.(covUpd - cov))
            (log10(dist) < -10) && (converged = true)
            converged && break
            cov = damp * covUpd + (1-damp) * cov
        end
        !converged && println("Covariance matrix: failed convergence.")
    end
    cov = round.(cov, digits=10)
    stdev = sqrt.(diag(cov))
    stdev_mat = diagm(0 => sqrt.(diag(cov)))
    inv_stdev_mat = inv(stdev_mat)
    inv_stdev_mat[stdev_mat .== 0] .= 0.0
    cor = inv_stdev_mat * cov * inv_stdev_mat

    return cov, cor
end

"""
    X = path(m, sol)
    X = path(m, sol; T, ϵ, x0, displayFigure, varIndex, labels, title)

Compute a path for the model `m` with solution `sol`. Returns an `n`+`m` × `T` array `X`.

# Keyword Arguments
`T`: number of periods. Default = 25.

`ϵ`: `q` × `T` vector of structural shocks. Default drawn from i.i.d. Gaussian.

`x0`: initial condition for state variables. Zeros if not specified.

`displayFigure`: `TRUE` to display figure. Default = `TRUE`.

`varIndex`: array with indices of variables to be displayed in the plot. All variables assigned if not specified.

`labels`: one-dimensional string array with labels for the variables ploted.

`title`: title of the plot.

`saveName`: if passed, save figure under the name `saveName`.
"""
function path(m0::Model, sol::Solution; T::Int64=25, ϵ=randn(m0.q,T), x0=standardinitialx0(m0), displayFigure = true, varIndex = 1:(m0.n+m0.m), labels=nothing, title=nothing, saveName=nothing)
    !sol.flag_rank && error("Rank condition not verified")

    n, m, q = m0.n, m0.m, m0.q
    P, Q = sol.P, sol.Q
    Iaux = addauxiliaries(m0)[2]
    if size(ϵ,1) != q
        error("Matrix ϵ must have q = $(q) rows.")
    end
    if (n > 0) && (size(x0,1) != n)
        error("Initial state x0 must have n = $(n) rows.")
    end
    if (size(ϵ,2) < T)
        error("Shock sequence ϵ must have at least T = $(T) columns.")
    end

    Path = zeros(n + m, T)
    no_state = (n==0)
    for t in 1:T
        if no_state
            Path[:,t] = Q * ϵ[:,t]
        else
            state = (t==1) ? x0 : Path[1:n, t-1]
            Path[:, t] = P * state + Q * ϵ[:,t]
        end
    end

    if displayFigure || !isnothing(saveName)
        label = isnothing(labels) ? :none : permutedims(labels)
        fig = plot(xlim=(1,Inf), gridalpha=0.05, legend=:topright)

        plot!(fig, 1:T, Path[varIndex,:]', markershape =:circle, lw = 2, ms = 6, label = label)

        !isnothing(title) && plot!(fig, title=title)

        plot!(fig, 1:T, zeros(T), color =:black, alpha=0.15, label = :none)
    end

    if displayFigure
        plot!(fig, tickfontsize = 12)
        display(fig)
    end

    if !isnothing(saveName)
        plot!(fig, tickfontsize = 14, legendfontsize = 14)
        savefig(fig, "$(saveName)")
    end

    return Path
end

function movavg(m::Model, sol::Solution; T::Int64 = 25)

    nn = m.n
    nm = m.m
    nq = m.q
    Ix = 1:nn
    Iy = nn+1:nn+nm
    
    Px = sol.P[Ix,:]
    Py = sol.P[Iy,:]

    Qx = sol.Q[Ix,:]
    Qy = sol.Q[Iy,:]

    MA = zeros(nn+nm, nq, T)
    MA[Ix,:,1] = Qx
    MA[Iy,:,1] = Qy
    for t in 2:T
        MA[Ix,:,t] = (Px^(t-1))*Qx
        MA[Iy,:,t] = Py*(Px^(t-2))*Qx
    end
    return MA
end


"""
    X = irf(m, sol, i)
    X = irf(m, sol, i; T, displayFigure, varIndex, labels, saveName)

Compute the impulse response functions for model `m` with solution `sol` to a one sd. dev. shock in exogenous variable i. Returns an `n`+`m` × `T` array `X`. All variables represented as deviations from their steady state values.

# Keyword Arguments
`T`: number of periods. 25 if not specified.

`displayFigure`: `TRUE` to display figure. `TRUE` if not specified.

`varIndex`: array with indices of variables to be displayed in the plot. All variables assigned if not specified.

`labels`: one-dimensional string array with labels for the variables ploted.

`title`: title of the plot.

`saveName`: if passed, save figure under the name `saveName`.

"""
function irf(m0::Model, sol::Solution, i::Int64; T::Int64=25, displayFigure = true, varIndex = 1:(m0.n+m0.m), labels=nothing, title = nothing, saveName = nothing)
    n, m, q = m0.n, m0.m, m0.q
    Irf = zeros(n+m, T)
    ϵ = zeros(q,T)
    ϵ[i,1] = 1.
    Irf = path(m0, sol; T=T, ϵ=ϵ, displayFigure = displayFigure, varIndex = varIndex, labels=labels, title=title, saveName = saveName)
    return Irf
end


"""
    X = irf(m, sol)
    X = irf(m, sol; T, displayFigure, varIndex, labels)

Compute the impulse response functions for model `m` with solution `sol` to every innovation in the model. Returns an `n`+`m` × `T` × `q` array `X`. All variables represented as deviations from their steady state values. No `saveName` option.

"""
function irf(m0::Model, sol::Solution; T::Int64=25, displayFigure = true, varIndex = 1:(m0.n+m0.m), labels=nothing)
    n, m, q = m0.n, m0.m, m0.q
    Irf = zeros(n+m, T, q)
    for iq in 1:q
        ϵ = zeros(q,T)
        ϵ[iq,1] = 1.0
        Irf[:,:,iq] = irf(m0, sol, iq; T=T, displayFigure = displayFigure, varIndex = varIndex, labels=labels, title="Shock #$(iq)")
    end
    return Irf
end


"""
    X = vardecomp(m, sol)
    X = vardecomp(m, sol; T)

Compute the variance decomposition of `T` period ahead forecast errors for model `m` with solution `sol`. Returns an `n`+`m` × `q` array `V`. The `i`,`j` element contains the `i`-th variable variance due to shock `j`. If `Σ` is a diagonal matrix and `T` is large, the elements of each row combine to the total variance of the corresponding variable.

# Keyword Argument
`T`: forecast error horizon. (default = 1)

"""
function vardecomp(m0::Model, sol::Solution; T = 1)
    MA = movavg(m0, sol, T = T)
    nn, nm, nq = m0.n, m0.m, m0.q 
    neq = nn + nm

    V = zeros(neq, nq)
    for j in 1:nq
        v = zeros(neq)
        Iq = zeros(nq, nq)
        Iq[j,j] = m0.Σ[j,j]
        for t in 1:T
            C = MA[:,:,t]
            Term = C*Iq*C' 
            v .+= diag(Term)
        end
        V[:,j] = v
    end
    return V
end