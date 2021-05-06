function ConvertReal(x, flag)
    flag = flag || maximum(abs.(imag.(x))) .> 1e-5
    x = real.(x)
    return x, flag
end
ConvertReal(x) = ConvertReal(x, true)[1]

function standardinitialx0(m)
    if m.n == 0
        return 0.0
    else
        return zeros(m.n)
    end
end
