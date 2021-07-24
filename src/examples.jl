using Plots:display
using Parameters, LinearAlgebra, Plots, LinearEcon

function price_stock()

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
    A[1,d] = 1
    B[1,d] = ρ
    C[1,e] = 1

    # equation (2): sₜ = dₜ + β Eₜsₜ₊₁
    A[2,s] = -β
    A[2,d] = -1
    B[2,s] = -1

    # covariance matrix
    Σ = [σ]

    m = Model(A, B, C, Σ, nn) # define model object
    R, FlagConverged = Correlation(m.S)
    display(FlagConverged)
    display(R)

    # plot options
    T = 12 # periods in the irf
    Labels = ["Dividend", "Stock Price"]

    # response to surplus shock
    Fig = IRF(m.S, 1, T=T, Labels=Labels)
    plot!(Fig, title="Dividend Shock")
    display(Fig)
    nothing
end
price_stock()

function old_keynesian()

    #=  OLD KEYNESIAN MODEL
        (1) xₜ = xₜ₋₁ - γ (iₜ - πₜ₋₁) 
        (2) πₜ = β πₜ₋₁ + κ xₜ
        (3) iₜ = θ₁ xₜ + θ₂ πₜ + ϵ

        st.dev.(ϵ) = σ
        2 states (x and π) =#

    neq = 3 # number of equations/variables
    nn = 2 # number of state variables
    nm = 1 # number of non-state variables
    nq = 1 # number of exogenous variables

    # parameters
    γ = 1.0
    β = 0.98
    κ = 0.5
    θ1 = 0.5
    θ2 = 1.5
    σ = 1.0

    # define variable indices
    x, pi, ir = 1:3
    e = 1

    # define matrices
    A = zeros(neq, neq)
    B = zeros(neq, neq)
    C = zeros(neq, nq)
    Σ = zeros(nq, nq)

    # equation (1): xₜ = xₜ₋₁ - γ (iₜ - πₜ₋₁) 
    A[1,x] = 1
    B[1,x] = 1
    B[1,ir] = -γ
    B[1,pi] = γ

    # equation (2): πₜ = β πₜ₋₁ + κ xₜ
    A[2,pi] = 1
    A[2,x] = -κ
    B[2,pi] = β

    # equation (3): iₜ = θ₁ xₜ + θ₂ πₜ + ϵ
    A[3,x] = -θ1
    A[3,pi] = -θ2
    B[3,ir] = -1
    C[3,e] = 1

    # covariance matrix
    Σ = [σ]

    m = Model(A, B, C, Σ, nn) # define model object
    R, FlagConverged = Correlation(m.S)
    display(FlagConverged)
    display(R)

    # plot options
    T = 12 # periods in the irf
    Labels = ["Output Gap", "Inflation", "Interest Rate"]

    # response to surplus shock
    Fig = IRF(m.S, 1, T=T, Labels=Labels)
    plot!(Fig, title="Monetary Policy Shock")
    display(Fig)
    nothing
end
old_keynesian()

function frictionless_mon_dominance()
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

    # define variable indices
    v, r, pi, ir = 1:4
    e1, e2 = 1:2

    # define matrices
    A = zeros(neq, neq)
    B = zeros(neq, neq)
    C = zeros(neq, nq)
    Σ = zeros(nq, nq)

    # equation (1): vₜ = ρ vₜ₋₁ + ϵ₁ₜ
    A[1,v] = 1
    B[1,v] = ρ
    C[1,e1] = 1

    # equation (2): rₜ = ν rₜ₋₁ + ϵ₂ₜ
    A[2,r] = 1
    B[2,r] = ν
    C[2,e2] = 1

    # equation (3): iₜ = rₜ + Eₜπₜ₊₁
    A[3,pi] = -1
    A[3,r] = -1
    B[3,ir] = -1

    # equation (4): iₜ = θ πₜ + vₜ
    A[4,v] = -1
    B[4,pi] = θ
    B[4,ir] = -1

    # covariance matrix
    Σ[e1,e1] = σ₁
    Σ[e2,e2] = σ₂

    m = Model(A, B, C, Σ, nn) # define model
    V = VarDecomp(m.S, 500) # variance decomposition
    display(pie(["Mon. Pol Shock", "Real Interest Shock"], V[pi,:], title="Variance Decomposition"))

    Vars = 2:4 # display real interest, inflation and nominal rate
    Labels = ["Real Interest", "Inflation", "Nominal Interest"]

    Fig = IRF(m.S, 2, Vars=Vars, Labels=Labels)
    plot!(Fig, title="Natural Rate Shock")
    display(Fig)

    Fig = Simulation(m.S, Vars=Vars, Labels=Labels)
    plot!(Fig, title="Simulation")
    display(Fig)
    
    nothing
end
frictionless_mon_dominance()

function three_eq_nk()

    #=  THREE EQUATION NK MODEL
        (1) xₜ = Eₜxₜ₊₁ - γ (iₜ - Eₜπₜ₊₁) 
        (2) πₜ = β Eₜπₜ₊₁ + κ xₜ
        (3) iₜ = θ₁ xₜ + θ₂ πₜ + ϵ

        st.dev.(ϵ) = σ
        No states =#

    neq = 3 # number of equations/variables
    nn = 0 # number of state variables
    nm = 3 # number of non-state variables
    nq = 1 # number of exogenous variables

    # parameters
    γ = 1.0
    β = 0.98
    κ = 0.5
    θ1 = 0.5
    θ2 = 1.5
    σ = 1.0

    # define variable indices
    x, pi, ir = 1:3
    e = 1

    # define matrices
    A = zeros(neq, neq)
    B = zeros(neq, neq)
    C = zeros(neq, nq)
    Σ = zeros(nq, nq)

    # equation (1): xₜ = Eₜxₜ₊₁ - γ (iₜ - Eₜπₜ₊₁) 
    A[1,x] = -1
    A[1,pi] = -γ
    B[1,x] = -1
    B[1,ir] = -γ

    # equation (2): πₜ = β Eₜπₜ₊₁ + κ xₜ
    A[2,pi] = -β
    B[2,x] = κ
    B[2,pi] = -1

    # equation (3): iₜ = θ₁ xₜ + θ₂ πₜ + ϵ
    B[3,x] = θ1
    B[3,pi] = θ2
    B[3,ir] = -1
    C[3,e] = 1

    # covariance matrix
    Σ = [σ]

    m = Model(A, B, C, Σ, nn) # define model object

    # plot options
    T = 12 # periods in the irf
    Labels = ["Output Gap", "Inflation", "Interest Rate"]

    # response to surplus shock
    Fig = IRF(m.S, 1, T=T, Labels=Labels)
    plot!(Fig, title="Monetary policy shock")
    display(Fig)
    nothing
end
three_eq_nk()

function fiscal_dom_nk()

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
    x, πi, r, q = nn + 1:neq
    e1, e2 = 1:nq

    if true
        # define matrices
        A = zeros(neq, neq)
        B = zeros(neq, neq)
        C = zeros(neq, nq)
        Σ = zeros(nq, nq)

        # equation (1): xₜ = Eₜxₜ₊₁ - γ (iₜ - Eₜπₜ₊₁) 
        i = 1
        A[i,x] = -1
        B[i,x] = -1
        A[i,ir] = γ
        A[i,πi] = -γ

        # equation (2): πₜ = β Eₜπₜ₊₁ + κ xₜ + e1ₜ
        i += 1
        B[i,πi] = -1
        A[i,πi] = -β
        B[i,x] = κ
        C[i,e1] = 1

        # equation (3): iₜ = ϕ πₜ + ρ iₜ₋₁ + e2ₜ
        i += 1
        A[i,ir] = 1
        B[i,πi] = ϕ
        B[i,ir] = ρ
        C[i,e2] = 1

        # equation (4): β vₜ = vₜ₋₁ + rₜ - πₜ
        i += 1
        A[i,v] = β
        B[i,v] = 1
        B[i,r] = 1
        B[i,πi] = -1

        # (5) rₜ = ωβqₜ - qₜ₋₁
        i += 1
        B[i,r] = -1
        B[i,q] = ω * β
        B[i,ql] = -1

        # (6) qₜ = -iₜ + ωβ Eₜqₜ₊₁
        i += 1
        B[i,q] = -1
        A[i,ir] = 1
        A[i,q] = -ω * β

        # (7) qlₜ = qₜ
        i += 1
        A[i,ql] = 1
        B[i,q] = 1

        # covariance matrix
        Σ = [1 0; 0 σ]

        m = Model(A, B, C, Σ, nn) # define model object
    end

    # plot options
    T = 12 # periods in the irf
    Labels = ["Interest Rate", "Public Debt", "Debt Return", "Inflation"]
    Vars = [ir, v, r, πi]

    # response to surplus shock
    Fig = IRF(m.S, e1, T=T, Labels=Labels, Vars=Vars)
    plot!(Fig, title="Monetary policy shock", ylim=(-0.75, 0.75))
    display(Fig)

    V = Covariance(m.S)[1]
    100 * V[πi, πi] |> display

    nothing
end
fiscal_dom_nk()

function open_econ()

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
    πi, q, nx, Δe, is = nn + 1:neq
    e1, e2, e3 = 1:nq
    
    if true
        # define matrices
        A = zeros(neq, neq)
        B = zeros(neq, neq)
        C = zeros(neq, nq)
        Σ = zeros(nq, nq)
        
        # equation (1): πₜ = πHₜ + (λ/(1-λ)) Δqₜ
        i = 1
        B[i,πi] = -1
        A[i,πH] = -1
        B[i,q] = λ / (1 - λ) 
        B[i,ql] = -λ / (1 - λ) 
        
        # equation (2): isₜ = ϕ bFₜ + uₜ
        i += 1
        B[i,is] = -1
        A[i,bF] = -ϕ
        A[i,u] = -1
        
        # equation (3): qₜ = Eₜqₜ₊₁ + isₜ + Eₜπₜ₊₁
        i += 1
        B[i,q] = -1
        A[i,q] = -1
        B[i,is] = 1
        A[i,πi] = -1

        # equation (4): nxₜ = λ (qₜ - cₜ)
        i += 1
        B[i,nx] = -1
        B[i,q] = λ
        A[i,c] = λ

        # equation (5): β bFₜ = bFₜ₋₁ - nxₜ + Γ (isₜ + Δeₜ - πHₜ)
        i += 1
        A[i,bF] = β
        B[i,bF] = 1
        B[i,nx] = -1
        B[i,is] = Γ
        B[i,Δe] = Γ
        A[i,πH] = Γ

        # equation (6):  Δeₜ = Δqₜ + πₜ = qₜ - qlₜ₋₁ + πₜ
        i += 1
        B[i,Δe] = -1
        B[i,q] = 1
        B[i,ql] = -1
        B[i,πi] = 1

        # equation (7):  cₜ = ρ1 cₜ₋₁ + ϵ1ₜ
        i += 1
        A[i,c] = 1
        B[i,c] = ρ1
        C[i,e1] = 1
        # B[i,is] = -0.3

        # equation (8): πHₜ = ρ2 πHₜ₋₁ + ϵ2ₜ 
        i += 1
        A[i,πH] = 1
        B[i,πH] = ρ2
        C[i,e2] = 1

        # equation (9): uₜ = ρ₃ uₜ₋₁ + ϵ3ₜ
        i += 1
        A[i,u] = 1
        B[i,u] = ρ3
        C[i,e3] = 1
        
        # equation (10): qlₜ = qₜ
        i += 1
        A[i,ql] = 1
        B[i,q] = 1

        # covariance matrix
        Σ = Float64.(collect(I(nq)))
    end
    
    m = Model(A, B, C, Σ, nn) # define model object
    
    # plot options
    T = 12 # periods in the irf
    Vars = [c, nx, bF, Δe]
    Labels = ["c", "nx", "bF,", "Δe"]
    
    # display(m.flag_complex)
    # display(m.S.P)
    
    # response to shocks
    for e in 1:nq
        Fig = IRF(m.S, e, T=T, Vars=Vars, Labels=Labels)
        plot!(Fig, title="Response to shock $e")
        display(Fig)
    end

    Vars = [bF, is, Δe, nx]
    Labels = ["bF", "is,", "Δe", "nx"]
    # Fig = IRF(m.S, e3, T=T, Vars=Vars, Labels=Labels)
    # plot!(Fig, title="Response to shock $e3")
    # display(Fig)

    nothing
end
open_econ()