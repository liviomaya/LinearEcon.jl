function checkreal(x)
    (sum(abs.(imag.(x)) .> 1e-5) > 0) && println("Imaginary component found in the solution.")
    xt = real.(x)
    return xt
end

function standardinitialx0(m,sol)
    if m.n == 0
        return 0.0
    else
        SS = ss(m,sol)
        return SS[1:m.n,1]
    end
end

function dpdiscountedsum(dp::dpsolution, x)
    r = size(dp.Mx0, 2)
    (size(x,1) != r) && error("Deterministic process time series must be r x T.")
    M0 = [dp.Mx0; dp.My0; dp.Mz0]
    α = [dp.αx; dp.αy; dp.αz]
    Λ2 = dp.Λ2
    Lt2 = dp.Lt2
    Id = diagm(0 => ones(size(Λ2,1)))
    k = size(x,2)

    S = zeros(size(M0,1),1)
    S += M0 * x[:,1]
    if k == 1
        S += α * inv(Λ2) * inv(Id - inv(Λ2)) * Lt2 * x[:,k]
    else
        S += α * Λ2^(-(k-1)) * inv(Id - inv(Λ2)) * Lt2 * x[:,k]
        if k > 2
            S += sum([α*Λ2^(-(i-1))*Lt2*x[:,i] for i in 2:k-1])
        end
    end
    return S
end
