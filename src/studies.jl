"""
    x = ss(m, sol)

Compute the steady state for the model `m` with solution `sol`.
"""
function ss(m::model, sol::solution)
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
"""
function covariance(m::model, sol::solution)

    if m.n == 0
        cov = sol.Q * sol.Q'
    else
        Px = sol.P[1:m.n, :]
        Qx = sol.Q[1:m.n, :]
        xcov = Qx * Qx'
        itermax = 1000
        iter = 1
        while iter <= itermax
            xcovUpd = Px * xcov * Px' + Qx * Qx'
            dist = maximum(abs.(xcovUpd - xcov))
            (log10(dist) < -7) && break
            iter += 1
            xcov = 0.50 * xcovUpd + 0.50 * xcov
        end
        if m.m + m.p > 0
            Pnx = sol.P[m.n+1:end, :]
            Qnx = sol.Q[m.n+1:end, :]
            ycov = Pnx * xcov * Pnx' + Qnx * Qnx'
            yxcov = Pnx * xcov * Px' + Qnx * Qx'
            cov = [xcov yxcov'; yxcov ycov]
        else
            cov = xcov
        end
    end
    stdev = sqrt.(diag(cov))
    stdev_mat = diagm(0 => sqrt.(diag(cov)))
    cor = inv(stdev_mat) * cov * inv(stdev_mat)

    return cov, cor
end

"""
    X = path(m, sol)
    X = path(m, sol; T, ϵ, dp, devSS, x0, displayFigure, varIndex, labels, title)

Compute a path for the model `m` with solution `sol`. Returns an `n`+`m`+`p` × `T` array `X`.

# Keyword Arguments
`T`: number of periods. 25 if not specified.

`ϵ`: `q` × `T` vector of structural shocks. Drawn from i.i.d. Gaussian if not specified.

`dp`: `r` × `T` deterministic process. Zeros if not specified. `path` assumes that the final value of the deterministic process is the last entry in `dp`.

`devSS`: `TRUE` to represent variables as deviation from their steady-state values. `TRUE` if not specified.

`x0`: initial condition for state variables. Steady state values if not specified.

`displayFigure`: `TRUE` to display figure. `TRUE` if not specified.

`varIndex`: array with indices of variables to be displayed in the plot. All variables assigned if not specified.

`labels`: one-dimensional string array with labels for the variables ploted.

`title`: title of the plot.

`saveName`: if passed, save figure under the name `saveName`.
"""
function path(m0::model, sol::solution; T::Int64=25, ϵ=randn(m0.q,T), dp=zeros(m0.r,T), devSS = true, x0=standardinitialx0(m0,sol), displayFigure = true, varIndex = 1:(m0.n+m0.m+m0.p), labels=nothing, title=nothing, saveName=nothing)

    n, m, p, q, r = m0.n, m0.m, m0.p, m0.q, m0.r
    Iaux = addauxiliaries(m0)[2]
    if size(ϵ,1) != q
        error("Matrix ϵ must have q = $(q) rows.")
    end
    if size(dp,1) != r
        error("Matrix dp must have r = $(r) rows.")
    end
    if (n > 0) && (size(x0,1) != n)
        error("Initial state x0 must have n = $(n) rows.")
    end
    if (size(ϵ,2) < T)
        error("Shock sequence ϵ must have at least T = $(T) columns.")
    end
    Path = zeros(n + m + p, T)
    onlypredet = (m + p == 0)
    nopredet = (n == 0)
    alltypes = !onlypredet && !nopredet
    !nopredet && (Path[1:n, 1] = x0)
    for t in 1:T
        dpfuturepath = (t >= size(dp,2)) ? dp[:,end] : dp[:,t:end]
        dpcomponent = dpdiscountedsum(sol.dp, dpfuturepath)[.!Iaux]
        if nopredet
            Path[:,t] = sol.c + sol.Q * ϵ[:,t] + dpcomponent
        else
            state = (t==1) ? x0 : Path[1:n, t-1]
            Path[1:n, t] = sol.c[1:n,1] + sol.P[1:n,:] * state + sol.Q[1:n, :] * ϵ[:,t] + dpcomponent[1:n,1]
            if !onlypredet
                Path[n+1:end, t] = sol.c[n+1:end,1] + sol.P[n+1:end,:] * state + sol.Q[n+1:end, :] * ϵ[:,t] + dpcomponent[n+1:n+m+p,1]
            end
        end
    end

    if devSS
        SS = ss(m0,sol)
        Path .-= SS*ones(1,T)
    end

    if displayFigure || !isnothing(saveName)
        label = isnothing(labels) ? :none : permutedims(labels)
        fig = plot(xlim=(1,Inf), gridalpha=0.05, legend=:topright)
        plot!(fig, 1:T, Path[varIndex,:]', markershape =:circle, lw = 2, ms = 6, title = title, label = label)
        devSS && plot!(fig, 1:T, zeros(T), color =:black, alpha=0.15, label = :none)
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




"""
    X = irf(m, sol, i)
    X = irf(m, sol, i; T, displayFigure, varIndex, labels, saveName)

Compute the impulse response functions for model `m` with solution `sol` to a one sd. dev. shock in exogenous variable i. Returns an `n`+`m`+`p` × `T` array `X`. All variables represented as deviations from their steady state values.

# Keyword Arguments
`T`: number of periods. 25 if not specified.

`displayFigure`: `TRUE` to display figure. `TRUE` if not specified.

`varIndex`: array with indices of variables to be displayed in the plot. All variables assigned if not specified.

`labels`: one-dimensional string array with labels for the variables ploted.

`title`: title of the plot.

`saveName`: if passed, save figure under the name `saveName`.

"""
function irf(m0::model, sol::solution, i::Int64; T::Int64=25, displayFigure = true, varIndex = 1:(m0.n+m0.m+m0.p), labels=nothing, title = nothing, saveName = nothing)
    n, m, p, q, r = m0.n, m0.m, m0.p, m0.q, m0.r
    Irf = zeros(n+m+p, T)
    ϵ = zeros(q,T)
    ϵ[i,1] = 1.
    Irf = path(m0, sol; T=T, ϵ=ϵ, displayFigure = displayFigure, varIndex = varIndex, labels=labels, title=title)
    return Irf
end


"""
    X = irf(m, sol)
    X = irf(m, sol; T, displayFigure, varIndex, labels)

Compute the impulse response functions for model `m` with solution `sol` to every innovation in the model. Returns an `n`+`m`+`p` × `T` × `q` array `X`. All variables represented as deviations from their steady state values. No `saveName` option.

"""
function irf(m0::model, sol::solution; T::Int64=25, displayFigure = true, varIndex = 1:(m0.n+m0.m+m0.p), labels=nothing)
    n, m, p, q, r = m0.n, m0.m, m0.p, m0.q, m0.r
    Irf = zeros(n+m+p, T, q)
    for iq in 1:q
        ϵ = zeros(q,T)
        ϵ[iq,1] = 1.0
        Irf[:,:,iq] = irf(m0, sol, iq; T=T, displayFigure = displayFigure, varIndex = varIndex, labels=labels, title="Shock #$(iq)")
    end
    return Irf
end
