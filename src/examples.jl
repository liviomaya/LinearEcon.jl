include("src\\RESOLVER.jl");

# comment
function price_stock()

    #=  PRICE OF A STOCK
        (1) dₜ = ρ dₜ₋₁ + ϵₜ
        (2) sₜ = dₜ + β Eₜsₜ₊₁

        st.dev(ϵ) = σ
    =#

    neq = 2 # number of equations/variables
    nn = 1 # number of state variables
    nm = 1 # number of forward looking variables
    np = 0 # number of static variables
    nq = 1 # number of exogenous variables

    # parameters
    ρ = 0.6
    β = 0.98
    σ = 1.0

    # define variable indices
    d,s = 1:2
    e = 1

    # define matrices
    A = zeros(neq,nn+nm)
    B = zeros(neq)
    C = zeros(neq,neq)
    D = zeros(neq,nq)
    Σ = zeros(nq,nq)

    # equation (1): dₜ = ρ dₜ₋₁ + ϵₜ 
    A[1,d] = 1
    C[1,d] = ρ
    D[1,e] = 1

    # equation (2): sₜ = dₜ + β Eₜsₜ₊₁
    A[2,s] = -β
    A[2,d] = -1
    C[2,s] = -1

    # covariance matrix
    Σ = [σ]

    m = model(A,B,C,D,Σ,nn) # define model object
    sol = solution(m) # solve model
    SS = ss(m, sol) # calculate steady state
    Cov, Cor = covariance(m, sol) # covariance and correlation matrices

    # plot options
    T = 12 # periods in the irf
    labels = ["Dividend", "Stock Price"]

    # response to surplus shock
    irf(m, sol, 1, T=T, labels=labels, title="Dividend shock")
    nothing
end
price_stock()

function old_keynesian()

    #=  OLD KEYNESIAN MODEL
        (1) xₜ = xₜ₋₁ - γ (iₜ - πₜ₋₁) 
        (2) πₜ = β πₜ₋₁ + κ xₜ
        (3) iₜ = θ₁ xₜ + θ₂ πₜ + ϵ

        st.dev.(ϵ) = σ
        2 states (x and π)
    =#

    neq = 3 # number of equations/variables
    nn = 2 # number of state variables
    nm = 0 # number of forward looking variables
    np = 1 # number of static variables
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
    A = zeros(neq,nn+nm)
    B = zeros(neq)
    C = zeros(neq,neq)
    D = zeros(neq,nq)
    Σ = zeros(nq,nq)

    # equation (1): xₜ = xₜ₋₁ - γ (iₜ - πₜ₋₁) 
    A[1,x] = 1
    C[1,x] = 1
    C[1,ir] = -γ
    C[1,pi] = γ

    # equation (2): πₜ = β πₜ₋₁ + κ xₜ
    A[2,pi] = 1
    A[2,x] = -κ
    C[2,pi] = β

    # equation (3): iₜ = θ₁ xₜ + θ₂ πₜ + ϵ
    A[3,x] = -θ1
    A[3,pi] = -θ2
    C[3,ir] = -1
    D[3,e] = 1

    # covariance matrix
    Σ = [σ]

    m = model(A,B,C,D,Σ,nn) # define model object
    sol = solution(m) # solve model
    SS = ss(m, sol) # calculate steady state
    Cov, Cor = covariance(m, sol) # covariance and correlation matrices

    # plot options
    T = 12 # periods in the irf
    labels = ["Output Gap", "Inflation", "Interest Rate"]

    # response to surplus shock
    irf(m, sol, 1, T=T, labels=labels, title="Monetary policy shock")
    nothing
end
old_keynesian()

function frictionless_mon_dominance()
    #= FRICTIONLESS MODEL UNDER MONETARY DOMINANCE
        (1) vₜ = ρ vₜ₋₁ + ϵ₁ₜ
        (2) rₜ = (1-ν) r̄ + ν rₜ₋₁ + ϵ₂ₜ
        (3) iₜ = rₜ + Eₜπₜ₊₁
        (4) iₜ = r̄ + πbar + θ (πₜ - πbar) + vₜ

        st.dev(ϵ₁) = σ₁
        st.dev(ϵ₂) = σ₂
    =#

    neq = 4 # number of equations/variables
    n = 2 # number of state variables
    m = 1 # number of forward looking variables
    p = 1 # number of static variables
    q = 2 # number of exogenous variables

    # parameters
    ρ = 0.75
    ν = 0.75
    r̄ = 0.01
    θ = 1.5
    πbar = 0.02
    σ₁ = 0.01
    σ₂ = 0.03

    # define variable indices
    v,r,pi,ir = 1:4
    e1, e2 = 1:2

    # define matrices
    A = zeros(neq,n+m)
    B = zeros(neq)
    C = zeros(neq,neq)
    D = zeros(neq,q)
    Σ = zeros(q,q)

    # equation (1): vₜ = ρ vₜ₋₁ + ϵ₁ₜ
    A[1,v] = 1
    C[1,v] = ρ
    D[1,e1] = 1

    # equation (2): rₜ = (1-ν) r̄ + ν rₜ₋₁ + ϵ₂ₜ
    A[2,r] = 1
    B[2] = (1-ν)*r̄
    C[2,r] = ν
    D[2,e2] = 1

    # equation (3): iₜ = rₜ + Eₜπₜ₊₁
    A[3,pi] = -1
    A[3,r] = -1
    C[3,ir] = -1

    # equation (4): iₜ = r̄ + πbar + θ (πₜ - πbar) + vₜ
    A[4,v] = -1
    B[4] = r̄ + (1-θ)*πbar
    C[4,pi] = θ
    C[4,ir] = -1

    # covariance matrix
    Σ[e1,e1] = σ₁
    Σ[e2,e2] = σ₂

    m = model(A,B,C,D,Σ,n) # define model
    sol = solution(m) # solve model
    SS = ss(m, sol) # calculate steady state

    varIndex = 2:4 # display real interest, inflation and nominal rate
    labels = ["Real Interest", "Inflation", "Nominal Interest"]

    irf(m, sol, varIndex=varIndex, labels=labels)
    #path(m, sol, varIndex=varIndex, labels=labels, title="Simulation", devSS = false)
    nothing
end
frictionless_mon_dominance()

function three_eq_nk()

    #=  THREE EQUATION NY MODEL
        (1) xₜ = Eₜxₜ₊₁ - γ (iₜ - Eₜπₜ₊₁) 
        (2) πₜ = β Eₜπₜ₊₁ + κ xₜ
        (3) iₜ = θ₁ xₜ + θ₂ πₜ + ϵ

        st.dev.(ϵ) = σ
        No states
    =#

    neq = 3 # number of equations/variables
    nn = 0 # number of state variables
    nm = 2 # number of forward looking variables
    np = 1 # number of static variables
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
    A = zeros(neq,nn+nm)
    B = zeros(neq)
    C = zeros(neq,neq)
    D = zeros(neq,nq)
    Σ = zeros(nq,nq)

    # equation (1): xₜ = Eₜxₜ₊₁ - γ (iₜ - Eₜπₜ₊₁) 
    A[1,x] = -1
    A[1,pi] = -γ
    C[1,x] = -1
    C[1,ir] = -γ

    # equation (2): πₜ = β Eₜπₜ₊₁ + κ xₜ
    A[2,pi] = -β
    C[2,x] = κ
    C[2,pi] = -1

    # equation (3): iₜ = θ₁ xₜ + θ₂ πₜ + ϵ
    C[3,x] = θ1
    C[3,pi] = θ2
    C[3,ir] = -1
    D[3,e] = 1

    # covariance matrix
    Σ = [σ]

    m = model(A,B,C,D,Σ,nn) # define model object
    sol = solution(m) # solve model
    SS = ss(m, sol) # calculate steady state
    Cov, Cor = covariance(m, sol) # covariance and correlation matrices

    # plot options
    T = 12 # periods in the irf
    labels = ["Output Gap", "Inflation", "Interest Rate"]

    # response to surplus shock
    irf(m, sol, 1, T=T, labels=labels, title="Monetary policy shock")
    nothing
end
three_eq_nk()

function nk_fisc_dominance()

    #=  NEW KEYNESIAN MODEL UNDER FISCAL DOMINANCE
        (1) ρ vₜ = vₜ₋₁ + iₜ - πₜ - sₜ
        (2) xₜ = Eₜxₜ₊₁ - σ (iₜ - Eₜπₜ₊₁)
        (3) πₜ = β Eₜπₜ₊₁ + κ xₜ
        (4) sₜ = τ xₜ + ϵₜ

        iₜ is a determistic process: agents anticipate interest path
    =#

    neq = 4 # number of equations/variables
    nn = 1 # number of state variables
    nm = 2 # number of forward looking variables
    np = 1 # number of static variables
    nq = 1 # number of exogenous variables
    nr = 1 # number of deterministic processes

    # parameters
    ρ = 0.98
    σ = 1.0
    β = 0.98
    κ = 0.5
    τ = 0.2
    σϵ = 1

    # define variable indices
    v, x, pi, s = 1:4
    e = 1
    ir = 1

    # define matrices
    A = zeros(neq,nn+nm)
    B = zeros(neq)
    C = zeros(neq,neq)
    D = zeros(neq,nq)
    Σ = zeros(nq,nq)
    F = zeros(neq,nr)

    # equation (1): ρ vₜ = vₜ₋₁ + iₜ - πₜ - sₜ
    A[1,v] = ρ
    C[1,v] = 1
    C[1,pi] = -1
    C[1,s] = -1
    F[1,ir] = 1

    # equation (2): xₜ = Eₜxₜ₊₁ - σ (iₜ - Eₜπₜ₊₁)
    A[2,x] = -1
    A[2,pi] = -σ
    C[2,x] = -1
    F[2,ir] = -σ

    # equation (3): πₜ = β Eₜπₜ₊₁ + κ xₜ
    A[3,pi] = -β
    C[3,pi] = -1
    C[3,x] = κ

    # equation (4): sₜ = τ xₜ + ϵₜ
    C[4,x] = τ
    C[4,s] = -1
    D[4,e] = 1

    # covariance matrix
    Σ = [σϵ]

    m = model(A,B,C,D,F,Σ,nn) # define model object
    sol = solution(m) # solve model
    SS = ss(m, sol) # calculate steady state
    Cov, Cor = covariance(m, sol) # covariance and correlation matrices

    # plot options
    varIndex = 1:4
    labels = ["Debt", "Output Gap", "Inflation", "Surplus"]
    T = 12 # periods in the irf

    # response to surplus shock
    irf(m, sol, 1, T=T, varIndex=varIndex, labels=labels, title="Surplus shock")

    # permanent increase in interest rates in period t > 1
    t = 1
    dp = zeros(nr,t)
    ϵ = zeros(nq,T)
    dp[ir,t] = 1
    path(m, sol, ϵ=ϵ, dp=dp, T=T, varIndex=varIndex, labels=labels, title = "Permanent increase in interest rate")
    nothing
end
nk_fisc_dominance()
