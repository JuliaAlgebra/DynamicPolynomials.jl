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
