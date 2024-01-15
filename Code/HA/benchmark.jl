
using Pkg
Pkg.activate("..") ## because we have environment files in the parent directory
Pkg.instantiate() ## to download all missing packages
using Revise

using BenchmarkTools, TimerOutputs, Distributions, QuantEcon, IterTools, Plots, Optim, Interpolations, LinearAlgebra, Inequality, Statistics, ColorSchemes,PrettyTables, Roots
using Profile, PProf
### LOAD MODULE 
includet("aiyagari_module.jl")
using .Aiyagari

# some simple examples first
b = 1.0
function stupid(a)
    global b
    tmp = a
    for i in 1:1_000_000
        tmp = tmp + a + b
    end
    return tmp
end

function stupid_2(a)
    tmp = a
    for i in 1:1_000_000
        tmp = tmp + a + b
    end
    return tmp
end

function smart(a,b)
    tmp = a
    for i in 1:1_000_000
        tmp = tmp + a + b
    end
    return tmp
end

@btime stupid(1.0)
@btime stupid_2(1.0)
@btime smart(1.0,b)
@btime smart(1.0,1.0)

@benchmark stupid(1.0)
@benchmark stupid_2(1.0)
@benchmark smart(1.0,b)
@btime smart(1.0,1.0)

@code_native stupid(1.0)
@code_native smart(1.0,b)

# what is going on here? - types: do they help?
    function sum_float_array(x::AbstractVector{<:Number})
        sum = 0.0
        for i in eachindex(x)
            sum += x[i]
        end
        return sum
    end

    x_range = range(0, 1, length = 100_000)
    x = collect(x_range)
    typeof(x)

    @btime sum_float_array($x)

    function sum_array(x)
        sum = 0.0
        for i in eachindex(x)
            sum += x[i]
        end
        return sum
    end
    @btime sum_array($x)

    @btime sum($x)
    @btime sum($x_range)


# try one of our functions 

hh = create_ha_block(ϕ = 1.0)
a_grid = create_grid(hh,a_max = 100)
prices = (r=0.01, w=1.0, τ_w = 0.0, τ_r = 0.0, d = 0.0)

@benchmark solve_hh_block($a_grid,$hh,$prices)

for n in [30, 40, 50]
    println("n: $n")
    a_grid = create_grid(hh,N_a = n, a_max = 100)
    @benchmark solve_hh_block($a_grid,$hh,$prices)
end


#### functions to test

hh = create_ha_block(ϕ = 1.0)
a_grid = create_grid(hh,a_max = 100)
prices = (r=0.01, w=1.0, τ_w = 0.0, τ_r = 0.0, d = 0.0)

@profile solve_hh_block(a_grid,hh,prices)
pprof(;webport=58699)

