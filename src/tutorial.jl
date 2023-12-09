using VectorAutoRegressions, Plots, LinearAlgebra, DataFrames
using LinearEcon

function stockprice()

    #=  PRICE OF A STOCK
        (1) dₜ = ρ dₜ₋₁ + ϵₜ
        (2) sₜ = dₜ + β Eₜsₜ₊₁

        st.dev(ϵ) = σ =#

    neq = 2 # number of equations/variables
    nn = 1 # number of state variables
    nm = 1 # number of non-state variables
    nq = 1 # number of exogenous variables

    # parameters
    ρ = 0.6
    β = 0.98
    σ = 1.0

    # define variable indices
    d, s = 1:2
    e = 1

    # define matrices
    A = zeros(neq, neq)
    B = zeros(neq, neq)
    C = zeros(neq, nq)
    Σ = zeros(nq, nq)

    # equation (1): dₜ = ρ dₜ₋₁ + ϵₜ 
    A[1, d] = 1
    B[1, d] = ρ
    C[1, e] = 1

    # equation (2): sₜ = dₜ + β Eₜsₜ₊₁
    A[2, s] = -β
    A[2, d] = -1
    B[2, s] = -1

    # covariance matrix
    Σ = [σ]

    m = LREM(A, B, C, Σ, nn) # define LREM
    sol = solve(m)

    R = cor(sol.var)

    # plot options
    T = 12 # periods in the irf
    label = ["Dividend" "Stock Price"]

    # response to surplus shock
    irf(sol.var, 1, T=T, options=Dict(:label => label))

    # response to surplus shock
    data = rand(sol.var, 100)
    display(plot(data))
    return
end
stockprice()

function oldkeynesian()

    #=  OLD KEYNESIAN MODEL
        (1) xₜ = xₜ₋₁ - γ (iₜ - πₜ₋₁) 
        (2) πₜ = β πₜ₋₁ + κ xₜ
        (3) iₜ = θ₁ xₜ + θ₂ πₜ + uₜ
        (4) uₜ = ρ uₜ₋₁ + ϵ

        st.dev.(ϵ) = σ
        2 states (x and π) =#

    neq = 4 # number of equations/variables
    nn = 3 # number of state variables
    nm = 1 # number of non-state variables
    nq = 1 # number of exogenous variables

    # parameters
    γ = 1.0
    β = 0.98
    κ = 0.5
    θ1 = 0.5
    θ2 = 1.5
    σ = 1.0
    ρ = 0.7

    # define variable indices
    x, pi, u = 1:nn
    ir = neq
    e = 1

    if true
        # define matrices
        A = zeros(neq, neq)
        B = zeros(neq, neq)
        C = zeros(neq, nq)
        Σ = zeros(nq, nq)

        # equation (1): xₜ = xₜ₋₁ - γ (iₜ - πₜ₋₁) 
        i = 1
        A[i, x] = 1
        B[i, x] = 1
        B[i, ir] = -γ
        B[i, pi] = γ

        # equation (2): πₜ = β πₜ₋₁ + κ xₜ
        i += 1
        A[i, pi] = 1
        A[i, x] = -κ
        B[i, pi] = β

        # equation (3): iₜ = θ₁ xₜ + θ₂ πₜ + uₜ
        i += 1
        A[i, x] = -θ1
        A[i, pi] = -θ2
        B[i, ir] = -1
        A[i, u] = -1

        # equation (4): uₜ = ρ uₜ₋₁ + ϵₜ
        i += 1
        A[i, u] = 1
        B[i, u] = ρ
        C[i, e] = 1

        # covariance matrix
        Σ = [σ]
    end

    lrem = LREM(A, B, C, Σ, nn) # define model object
    sol = solve(lrem)
    R = cor(sol.var)
    display(R)

    # plot options
    T = 12 # periods in the irf
    label = ["Output Gap" "Inflation" "Interest Rate"]

    # response to surplus shock
    irf(sol.var, 1, T=T, id=1:3, options=Dict(
        :label => label,
        :title => "Monetary Policy Shock"))
    return nothing
end
oldkeynesian()

function monetaryfrictionless()
    #= FRICTIONLESS MODEL UNDER MONETARY DOMINANCE
        (1) vₜ = ρ vₜ₋₁ + ϵ₁ₜ
        (2) rₜ = ν rₜ₋₁ + ϵ₂ₜ
        (3) iₜ = rₜ + Eₜπₜ₊₁
        (4) iₜ = θ πₜ + vₜ

        st.dev(ϵ₁) = σ₁
        st.dev(ϵ₂) = σ₂ =#

    neq = 4 # number of equations/variables
    nn = 2 # number of state variables
    nm = 2 # number of non-state variables
    nq = 2 # number of exogenous variables

    # parameters
    ρ = 0.75
    ν = 0.75
    θ = 1.5
    σ₁ = 0.01
    σ₂ = 0.02
    ψ = 0.95

    # define variable indices
    v, r, pi, ir = 1:4
    e1, e2 = 1:2

    if true
        # define matrices
        A = zeros(neq, neq)
        B = zeros(neq, neq)
        C = zeros(neq, nq)
        Σ = zeros(nq, nq)

        # equation (1): vₜ = ρ vₜ₋₁ + ϵ₁ₜ
        i = 1
        A[i, v] = 1
        B[i, v] = ρ
        C[i, e1] = 1

        # equation (2): rₜ = ν rₜ₋₁ + ϵ₂ₜ
        i += 1
        A[i, r] = 1
        B[i, r] = ν
        C[i, e2] = 1

        # equation (3): iₜ = rₜ + Eₜπₜ₊₁
        i += 1
        A[i, pi] = -1
        A[i, r] = -1
        B[i, ir] = -1

        # equation (4): iₜ = θ πₜ + vₜ
        i += 1
        A[i, v] = -1
        B[i, pi] = θ
        B[i, ir] = -1

        # covariance matrix
        Σ[e1, e1] = σ₁^2
        Σ[e2, e2] = σ₂^2
        Σ[e1, e2] = ψ * σ₁ * σ₂
        Σ[e2, e1] = ψ * σ₁ * σ₂
    end

    lrem = LREM(A, B, C, Σ, nn)
    sol = solve(lrem)

    # variance decomposition
    V = vardecomp(sol.var, 1) # variance decomposition
    display(pie(["Mon. Pol Shock", "Real Interest Shock"],
        V[pi, :], title="Variance Decomposition"))

    # cholesky
    vchol, comb = cholesky(sol.var) # cholesky factorization
    display(comb)
    Vchol = vardecomp(vchol, 1) # variance decomposition
    display(pie(["Shock 1", "Shock 2"],
        Vchol[pi, :], title="Variance Decomposition (Orthogonalized)"))

    # plot irf
    id = 2:4 # display real interest, inflation and nominal rate
    label = ["Real Interest" "Inflation" "Nominal Interest"]
    irf(sol.var, 1, id=id, options=Dict(:label => label))

    # simulation
    T = 100
    data = rand(sol.var, T)
    fig = plot(1:T, data[:, id], label=label, title="Simulation")
    display(fig)

    return nothing
end
monetaryfrictionless()

function nk3()

    #=  THREE EQUATION NK MODEL
        (1) xₜ = Eₜxₜ₊₁ - γ (iₜ - Eₜπₜ₊₁) 
        (2) πₜ = β Eₜπₜ₊₁ + κ xₜ
        (3) iₜ = θ₁ xₜ + θ₂ πₜ + uₜ
        (4) uₜ = ρ uₜ₋₁ + ϵₜ

        st.dev.(ϵ) = σ
        No states =#

    neq = 4 # number of equations/variables
    nn = 1 # number of state variables
    nm = neq - nn # number of non-state variables
    nq = 1 # number of exogenous variables

    # parameters
    γ = 1.0
    β = 0.98
    κ = 0.5
    θ1 = 0.0
    θ2 = 1.5
    σ = 1.0
    ρ = 0.0

    # define variable indices
    u = 1
    x, pi, ir = nn+1:neq
    e = 1

    if true
        # define matrices
        A = zeros(neq, neq)
        B = zeros(neq, neq)
        C = zeros(neq, nq)
        Σ = zeros(nq, nq)

        # equation (1): xₜ = Eₜxₜ₊₁ - γ (iₜ - Eₜπₜ₊₁) 
        i = 1
        A[i, x] = -1
        A[i, pi] = -γ
        B[i, x] = -1
        B[i, ir] = -γ

        # equation (2): πₜ = β Eₜπₜ₊₁ + κ xₜ
        i += 1
        A[i, pi] = -β
        B[i, x] = κ
        B[i, pi] = -1

        # equation (3): iₜ = θ₁ xₜ + θ₂ πₜ + uₜ
        i += 1
        B[i, x] = θ1
        B[i, pi] = θ2
        B[i, ir] = -1
        A[i, u] = -1

        # (4) uₜ = ρ uₜ₋₁ + ϵₜ
        i += 1
        A[i, u] = 1
        B[i, u] = ρ
        C[i, e] = 1

        # covariance matrix
        Σ = [σ]
    end

    lrem = LREM(A, B, C, Σ, nn) # define model object
    sol = solve(lrem)

    # plot options
    T = 3 # periods in the irf
    label = ["Output Gap" "Inflation" "Interest Rate"]
    id = [x, pi, ir]

    # response to surplus shock
    irf(sol.var, 1, T=T, id=id, options=Dict(
        :label => label,
        :title => "Monetary Policy Shocks"))
    return nothing
end
nk3()

function nkfiscal()

    #=  NK MODEL + FISCAL DOMINANCE
        (1) xₜ = Eₜxₜ₊₁ - γ (iₜ - Eₜπₜ₊₁) 
        (2) πₜ = β Eₜπₜ₊₁ + κ xₜ + e1ₜ
        (3) iₜ = ϕ πₜ + ρ iₜ₋₁ + e2ₜ
        (4) β vₜ = vₜ₋₁ + rₜ - πₜ
        (5) rₜ = ωβqₜ - qₜ₋₁
        (6) qₜ = -iₜ + ωβ Eₜqₜ₊₁
        (7) qlₜ = qₜ

        st.dev(e1) = 1
        st.dev(e2) = σ
        No states =#

    neq = 7 # number of equations/variables
    nn = 3 # number of state variables
    nq = 2 # number of exogenous variables

    # parameters
    γ = 1.0
    β = 0.98
    κ = 0.5
    ρ = 0.0
    ω = 0.7
    ϕ = 0.4
    σ = 1.0

    # define variable indices
    ir, v, ql = 1:nn
    x, πi, r, q = nn+1:neq
    e1, e2 = 1:nq

    if true
        # define matrices
        A = zeros(neq, neq)
        B = zeros(neq, neq)
        C = zeros(neq, nq)
        Σ = zeros(nq, nq)

        # equation (1): xₜ = Eₜxₜ₊₁ - γ (iₜ - Eₜπₜ₊₁) 
        i = 1
        A[i, x] = -1
        B[i, x] = -1
        A[i, ir] = γ
        A[i, πi] = -γ

        # equation (2): πₜ = β Eₜπₜ₊₁ + κ xₜ + e1ₜ
        i += 1
        B[i, πi] = -1
        A[i, πi] = -β
        B[i, x] = κ
        C[i, e1] = 1

        # equation (3): iₜ = ϕ πₜ + ρ iₜ₋₁ + e2ₜ
        i += 1
        A[i, ir] = 1
        B[i, πi] = ϕ
        B[i, ir] = ρ
        C[i, e2] = 1

        # equation (4): β vₜ = vₜ₋₁ + rₜ - πₜ
        i += 1
        A[i, v] = β
        B[i, v] = 1
        B[i, r] = 1
        B[i, πi] = -1

        # (5) rₜ = ωβqₜ - qₜ₋₁
        i += 1
        B[i, r] = -1
        B[i, q] = ω * β
        B[i, ql] = -1

        # (6) qₜ = -iₜ + ωβ Eₜqₜ₊₁
        i += 1
        B[i, q] = -1
        A[i, ir] = 1
        A[i, q] = -ω * β

        # (7) qlₜ = qₜ
        i += 1
        A[i, ql] = 1
        B[i, q] = 1

        # covariance matrix
        Σ = [1 0; 0 σ]
    end

    # define and solve
    lrem = LREM(A, B, C, Σ, nn)
    sol = solve(lrem)

    # plot options
    T = 12 # periods in the irf
    label = ["Interest Rate" "Public Debt" "Debt Return" "Inflation"]
    id = [ir, v, r, πi]

    # cholesky decomposition
    vchol, R = cholesky(sol.var)

    # response to cost-push shock
    irf(sol.var, e1, T=T, id=id, options=Dict(
        :title => "Cost-Push Shock",
        :label => label,
        :ylim => (-0.75, 0.75)
    ))

    # response to monetary policy shock
    irf(sol.var, e2, T=T, id=id, options=Dict(
        :title => "Monetary Policy Shock",
        :label => label,
        :ylim => (-1.0, 1.0)
    ))

    V = cov(sol.var)
    100 * sqrt(V[πi, πi]) |> display

    return nothing
end
nkfiscal()

function openecon()

    #=  SIMPLIFIED OPEN ECONOMY MODEL
        (1) πₜ = πHₜ + (λ/(1-λ)) qₜ
        (2) isₜ = ϕ bFₜ + uₜ
        (3) qₜ = Eₜqₜ₊₁ + isₜ + Eₜπₜ₊₁
        (4) nxₜ = λ (qₜ - cₜ)
        (5) β bFₜ = bFₜ₋₁ - nxₜ + Γ (isₜ + Δeₜ - πHₜ)
        (6) Δeₜ = Δqₜ + πₜ
        (7) cₜ = ρ1 cₜ₋₁ + ϵ1ₜ
        (8) πHₜ = ρ2 πHₜ₋₁ + ϵ2ₜ 
        (9) uₜ = ρ₃ uₜ₋₁ + ϵ3ₜ
        (10) qlₜ = qₜ =#

    neq = 10 # number of equations/variables
    nn = 5 # number of state variables
    nm = neq - nn # number of non-state variables
    nq = 3 # number of exogenous variables

    # parameters
    λ = 0.20
    ϕ = 0.03
    β = 0.98
    Γ = 1.0
    ρ1 = 0.9
    ρ2 = 0.9
    ρ3 = 0.0

    # define variable indices
    bF, c, πH, u, ql = 1:nn
    πi, q, nx, Δe, is = nn+1:neq
    e1, e2, e3 = 1:nq

    if true
        # define matrices
        A = zeros(neq, neq)
        B = zeros(neq, neq)
        C = zeros(neq, nq)
        Σ = zeros(nq, nq)

        # equation (1): πₜ = πHₜ + (λ/(1-λ)) Δqₜ
        i = 1
        B[i, πi] = -1
        A[i, πH] = -1
        B[i, q] = λ / (1 - λ)
        B[i, ql] = -λ / (1 - λ)

        # equation (2): isₜ = ϕ bFₜ + uₜ
        i += 1
        B[i, is] = -1
        A[i, bF] = -ϕ
        A[i, u] = -1

        # equation (3): qₜ = Eₜqₜ₊₁ + isₜ + Eₜπₜ₊₁
        i += 1
        B[i, q] = -1
        A[i, q] = -1
        B[i, is] = 1
        A[i, πi] = -1

        # equation (4): nxₜ = λ (qₜ - cₜ)
        i += 1
        B[i, nx] = -1
        B[i, q] = λ
        A[i, c] = λ

        # equation (5): β bFₜ = bFₜ₋₁ - nxₜ + Γ (isₜ + Δeₜ - πHₜ)
        i += 1
        A[i, bF] = β
        B[i, bF] = 1
        B[i, nx] = -1
        B[i, is] = Γ
        B[i, Δe] = Γ
        A[i, πH] = Γ

        # equation (6):  Δeₜ = Δqₜ + πₜ = qₜ - qlₜ₋₁ + πₜ
        i += 1
        B[i, Δe] = -1
        B[i, q] = 1
        B[i, ql] = -1
        B[i, πi] = 1

        # equation (7):  cₜ = ρ1 cₜ₋₁ + ϵ1ₜ
        i += 1
        A[i, c] = 1
        B[i, c] = ρ1
        C[i, e1] = 1
        # B[i,is] = -0.3

        # equation (8): πHₜ = ρ2 πHₜ₋₁ + ϵ2ₜ 
        i += 1
        A[i, πH] = 1
        B[i, πH] = ρ2
        C[i, e2] = 1

        # equation (9): uₜ = ρ₃ uₜ₋₁ + ϵ3ₜ
        i += 1
        A[i, u] = 1
        B[i, u] = ρ3
        C[i, e3] = 1

        # equation (10): qlₜ = qₜ
        i += 1
        A[i, ql] = 1
        B[i, q] = 1

        # covariance matrix
        Σ = diagm(ones(nq))
    end

    # define and solve
    lrem = LREM(A, B, C, Σ, nn) # define model object
    sol = solve(lrem)

    # plot options
    T = 12 # periods in the irf
    id = [c, nx, bF, Δe]
    label = ["c" "nx" "bF," "Δe"]

    # response to shocks
    for e in 1:nq
        irf(sol.var, e, T=T, id=id, options=Dict(
            :label => label,
            :title => "Response to shock $e"
        ))
    end

    return nothing
end
openecon()

function fiscalshocks()

    #= MULTIPLE FISCAL SHOCKS
    System linearized around SS with Y=1 and V=Y
    (1) cₜ = Eₜcₜ₊₁ - σ (iₜ - Eₜπₜ₊₁)
    (2) πₜ = β Eₜπₜ₊₁ + κ yₜ
    (3) yₜ = (1-γ) cₜ + γ gₜ
    (4) β vₜ = vₜ₋₁ + iₜ₋₁ - πₜ - sₜ
    (5) β vsₜ = vsₜ₋₁ + iₜ₋₁ - πsₜ - sₜ
    (6) πsₜ = Eₜ₋₁πₜ + ηₜ
    (7) sₜ = -(r+γ) trₜ - γ gₜ
    (8) trₜ = ρ trₜ₋₁ - α1 vsₜ + ϵ1ₜ
    (9) gₜ = ρ gₜ₋₁ - α2 vsₜ + ϵ2
    (10) iₜ = 0
    (11) Eₜ₋₁πₜ = Eₜ₋₁πₜ
    (12) iₜ = iₜ
    =#

    neq = 12# number of equations/variables
    nn = 6 # number of state variables
    nm = neq - nn # number of non-state variables
    nq = 3 # number of exogenous variables

    # parameters 
    σ = 0.5
    β = 0.98
    κ = 0.5
    γ = 0.2
    r = 1 / β - 1
    ρ = 0.90
    α1 = 0.8
    α2 = 0.0
    corvaltr = 0
    corvalg = 0
    cortrg = 0
    stdval = 1
    stdtr = 1
    stdg = 1

    # build Σ 
    stdarray = [stdval; stdtr; stdg]
    cormat = [
        1 corvaltr corvalg
        corvaltr 1 cortrg
        corvalg cortrg 1]
    Σ = diagm(0 => stdarray) * cormat * diagm(0 => stdarray)

    # define variables indices
    v, vs, tr, g, Epi, irlag = 1:nn
    c, pii, y, pis, s, ir = nn+1:neq
    eta, e1, e2 = 1:nq

    # model equations
    if true

        # define matrices
        A = zeros(neq, neq)
        B = zeros(neq, neq)
        C = zeros(neq, nq)
        Σ = zeros(nq, nq)

        # (1) cₜ = Eₜcₜ₊₁ - σ (iₜ - Eₜπₜ₊₁)
        i = 1
        A[i, [c, pii]] .= [-1, -σ]
        B[i, [c, ir]] .= [-1, -σ]

        # (2) πₜ = β Eₜπₜ₊₁ + κ yₜ
        i += 1
        A[i, [pii]] .= [-β]
        B[i, [pii, y]] .= [-1, κ]

        # (3) yₜ = (1-γ) cₜ + γ gₜ
        i += 1
        A[i, [g]] .= [-γ]
        B[i, [y, c]] .= [-1, 1 - γ]

        # (4) β vₜ = vₜ₋₁ + iₜ₋₁ - πₜ - sₜ
        i += 1
        A[i, [v]] .= [β]
        B[i, [v, irlag, pii, s]] .= [1, 1, -1, -1]

        # (5) β vsₜ = vsₜ₋₁ + iₜ₋₁ - πsₜ - sₜ
        i += 1
        A[i, [vs]] .= [β]
        B[i, [vs, irlag, pis, s]] .= [1, 1, -1, -1]

        # (6) πsₜ = Eₜ₋₁πₜ + ηₜ
        i += 1
        B[i, [pis, Epi]] .= [-1, 1]
        C[i, [eta]] .= [1]

        # (7) sₜ = -(r+γ) trₜ - γ gₜ
        i += 1
        A[i, [tr, g]] .= [γ + r, γ]
        B[i, [s]] .= [-1]

        # (8) trₜ = ρ trₜ₋₁ - α1 vsₜ + ϵ1ₜ
        i += 1
        A[i, [tr]] .= [1]
        B[i, [tr, vs]] .= [ρ, -α1]
        C[i, [e1]] .= [1]

        # (9) gₜ = ρ gₜ₋₁ - α2 vsₜ + ϵ2
        i += 1
        A[i, [g]] .= [1]
        B[i, [g, vs]] .= [ρ, -α2]
        C[i, [e2]] .= [1]

        # (10) iₜ = 0
        i += 1
        B[i, [ir]] .= [-1]

        # (11) Eₜ₋₁πₜ = Eₜ₋₁πₜ
        i += 1
        A[i, [pii, Epi]] .= [1, -1]

        # (12) iₜ = iₜ
        i += 1
        A[i, [irlag]] .= [1]
        B[i, [ir]] .= [1]
    end

    # define and solve 
    lrem = LREM(A, B, C, Σ, nn) # define model object
    sol = solve(lrem)
    display(sol.flagrank)
    !sol.flagrank && return

    # plot options
    T = 7 # periods in the irf
    id = [y, pii, v, s, g]
    label = ["y" "π" "v" "s" "g"]

    # response to shocks
    for e in 1:nq
        irf(sol.var, e, T=T, id=id, options=Dict(
            :label => label,
            :title => "Response to shock $e",
            :xlim => (1, T + 2),
            :ylim => (-2, 2),
            :xticks => 1:T))
    end

    shock = [1, inv(r + γ), 0]
    shock = [1, 0, inv(γ)]
    irf(sol.var, shock, T=T, id=id, options=Dict(
        :label => label,
        :title => "Response to shock",
        :xlim => (1, T + 2),
        :ylim => (-2, 2),
        :xticks => 1:T))

    # present discounted value
    eye = diagm(0 => ones(neq))
    Is = eye[:, s]
    Ipi = eye[:, pii]
    Ii = eye[:, ir]
    Ψ = sol.var.Ψ[1]
    Γ = sol.var.Γ
    tb = DataFrame()
    eye3 = diagm(0 => ones(nq))
    for e in 1:nq
        sPDV = Is' * inv(I(neq) .- β * Ψ) * Γ * eye3[:, e]
        piPDV = Ipi' * inv(I(neq) .- β * Ψ) * Γ * eye3[:, e]
        iPDV = -Ii' * inv(I(neq) .- β * Ψ) * Γ * eye3[:, e]
        df = DataFrame(
            "Shock" => [e],
            "s" => [round(sPDV, digits=5)],
            "-ir" => [round(iPDV, digits=5)],
            "pi" => [round(piPDV, digits=5)])
        tb = vcat(tb, df)
    end
    display(tb)

    # all entries should be zero:
    # (Is .- Ir .+ Ipi)' * inv(I(neq) .- β * Ψ) * Γ |> display

    return nothing
end
fiscalshocks()