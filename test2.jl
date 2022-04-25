function main()
    a = 2
    for _ in 1:4
        a = 5 + a
        if a > 5
            return a 
        end
    end
    a
end
a = main()
println(a)