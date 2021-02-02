function convert_real(x, flag)
    flag = flag || maximum(abs.(imag.(x))) .> 1e-5
    x = real.(x)
    return x, flag
end
convert_real(x) = convert_real(x, true)[1]

function standardinitialx0(m)
    if m.n == 0
        return 0.0
    else
        return zeros(m.n)
    end
end

function standardinitialu0(m)
    if m.p == 0
        return 0.0
    else
        return zeros(m.p)
    end
end
