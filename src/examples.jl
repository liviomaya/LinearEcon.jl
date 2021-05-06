
function price_stock()

    #=  PRICE OF A STOCK
        (1) dₜ = ρ dₜ₋₁ + ϵₜ
        (2) sₜ = dₜ + β Eₜsₜ₊₁

        st.dev(ϵ) = σ
    =#

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
    A = zeros(neq,neq)
    B = zeros(neq,neq)
    C = zeros(neq,nq)
    Σ = zeros(nq,nq)

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

    m = Model(A,B,C,Σ,nn) # define model object
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
        2 states (x and π)
    =#

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
    A = zeros(neq,neq)
    B = zeros(neq,neq)
    C = zeros(neq,nq)
    Σ = zeros(nq,nq)

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

    m = Model(A,B,C,Σ,nn) # define model object
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
        st.dev(ϵ₂) = σ₂
    =#

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
    v,r,pi,ir = 1:4
    e1, e2 = 1:2

    # define matrices
    A = zeros(neq,neq)
    B = zeros(neq,neq)
    C = zeros(neq,nq)
    Σ = zeros(nq,nq)

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

    m = Model(A,B,C,Σ,nn) # define model
    V = VarDecomp(m.S, 500) # variance decomposition
    display( pie(["Mon. Pol Shock", "Real Interest Shock"], V[pi,:], title="Variance Decomposition")  )

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
        No states
    =#

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
    A = zeros(neq,neq)
    B = zeros(neq,neq)
    C = zeros(neq,nq)
    Σ = zeros(nq,nq)

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

    m = Model(A,B,C,Σ,nn) # define model object

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
