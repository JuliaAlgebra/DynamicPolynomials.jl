function show(io::IO, m::Monomial)
    if isconstant(m)
        show(io, 1)
    else
        needsep = false
        for i in 1:nvars(m)
            if m.z[i] > 0
                if needsep
                    print(io, '*')
                end
                show(io, m.vars[i])
                if m.z[i] > 1
                    print(io, '^')
                    print(io, m.z[i])
                else
                    needsep = true
                end
            end
        end
    end
end

function show(io::IO, x::MonomialVector)
    print(io, typeof(x))
    print(io, "[ ")
    for (i, m) in enumerate(x)
        print(io, m)
        if i != length(x)
            print(io, ", ")
        end
    end
    print(io, " ]")
end
