
using Pkg
Pkg.activate("..") ## because we have environment files in the parent directory
Pkg.instantiate() ## to download all missing packages
using Revise

using Distributions, QuantEcon, IterTools, Plots, Optim, Interpolations, LinearAlgebra, Inequality, Statistics, ColorSchemes,PrettyTables, Roots



### LOAD MODULE 
includet("aiyagari_module.jl")
using .Aiyagari

# access to all exported functions
hh = create_ha_block(ϕ = 1.0)
a_grid = create_grid(hh,a_max = 100)

# try to see if it works as intended
prices = (r=0.01, w=1.0, τ_w = 0.0, τ_r = 0.0, d = 0.0)



V, policies, error, iter, λ, λ_vector, λ_a, λ_z, A′, C = solve_hh_block(a_grid,hh,prices)


println("error in VFI = $error")


# test plotting
    lines_scheme = [get(ColorSchemes.thermal,LinRange(0.0,1.0,hh.N_z));]
    value_plot = plot(xlabel = "a", ylabel = "V", title = "Value function")

    for j in 1:hh.N_z
        plot!(a_grid.a_vec[1:20], V[1:20,j], label = false, color = lines_scheme[j], lw=3)
    end

    plot!(a_grid.a_vec[1:20],V[1:20], label = false, linestyle = :dash, color = :black)


# now do the example of the full thing

    firm = create_firm()
    govt = create_government()
    B = govt.B
    L = get_aggregate_labor(hh)

    r_init = 0.0
    K_L, w = solve_firm(firm,r_init)
    K = K_L * L
    asset_supply = K + B
    
    prices = (r=r_init, w=w, τ_w = govt.τ_w, τ_r = govt.τ_r, d = govt.d)

    function asset_demand(prices)
        V, policies, error, iter, λ, λ_vector, λ_a, λ_z, A′, C = solve_hh_block(a_grid,hh,prices)
        return A′
    end

    asset_demand(prices)

    excess_demand = asset_demand(prices) - asset_supply

# put all pieces together -- this is not really safe... 
    function aiyagari_residual(r_guess)
        
        L = get_aggregate_labor(hh)
        K_L, w = solve_firm(firm,r_guess)
        K = K_L * L
        asset_supply = K + B
        prices = (r=r_guess, w=w, τ_w = govt.τ_w, τ_r = govt.τ_r, d = govt.d)
        residual = asset_demand(prices) - asset_supply

        return residual
    end

# test 

    r_guess = 0.01
    aiyagari_residual(r_guess)

# plot 
    r_vec = LinRange(-0.06,0.03,10)
    plot(r_vec, aiyagari_residual.(r_vec), label = false, lw = 3, color = :black, title = "Excess demand for assets", xlabel = "r", ylabel = "Excess demand")

# solve for an equilibrium return 
    # note - the upper bound will also depend on τ_r
    r_star = find_zero(aiyagari_residual, (-0.025, hh.β^(-1)-1-0.01 ) ,verbose = true)
    println("equilibrium real rate = $r_star")

 
### APPLICATION: TAXATION AND INEQUALITY
# rewrite the functions in a safe way

function asset_demand(prices,hh,a_grid)
    V, policies, error, iter, λ, λ_vector, λ_a, λ_z, A′, C = solve_hh_block(a_grid,hh,prices)
    return A′
end
function aiyagari_residual(r_guess,hh,firm,govt,a_grid)
        L = get_aggregate_labor(hh)
        K_L, w = solve_firm(firm,r_guess)
        K = K_L * L
        asset_supply = K + govt.B
        prices = (r=r_guess, w=w, τ_w = govt.τ_w, τ_r = govt.τ_r, d = govt.d)
        residual = asset_demand(prices,hh,a_grid) - asset_supply

        return residual
end

govt_low = create_government(τ_w = 0.0, τ_r = 0.0, d = 0.0)
govt_high = create_government(τ_w = 0.2, τ_r = 0.2, d = 0.05)
    
r_star_high = find_zero(x -> aiyagari_residual(x, hh,firm,govt_high,a_grid), (-0.025, hh.β^(-1)-1-0.001 ) ,verbose = true)
r_star_low  = find_zero(x -> aiyagari_residual(x, hh,firm,govt_low,a_grid), (-0.025, hh.β^(-1)-1-0.001 ) ,verbose = true)

resid_high = aiyagari_residual(r_star_high, hh,firm,govt_high,a_grid)
resid_low  = aiyagari_residual(r_star_low, hh,firm,govt_low,a_grid)

println("equilibrium real rate with high taxes = $r_star_high")
println("residual of the asset market clearing equation  = $resid_high")
println("equilibrium real rate with low taxes = $r_star_low")
println("residual of the asset market clearing equation  = $resid_low")
# compare two economies 
prices_high = (r=r_star_high, w=solve_firm(firm,r_star_high)[2], τ_w = govt_high.τ_w, τ_r = govt_high.τ_r, d = govt_high.d)
prices_low = (r=r_star_low, w=solve_firm(firm,r_star_low)[2], τ_w = govt_low.τ_w, τ_r = govt_low.τ_r, d = govt_low.d)

V_high, policies_high, error_high, iter_high, λ_high, λ_vector_high, λ_a_high, λ_z_high, A′_high, C_high  = solve_hh_block(a_grid,hh,prices_high)
V_low, policies_low, error_low, iter_low, λ_low, λ_vector_low, λ_a_low, λ_z_low, A′_low, C_low            = solve_hh_block(a_grid,hh,prices_low)

L = get_aggregate_labor(hh)

K_high = solve_firm(firm,r_star_high)[1]*L
K_low = solve_firm(firm,r_star_low)[1]*L

w_high = solve_firm(firm,r_star_high)[2]
w_low = solve_firm(firm,r_star_low)[2]

resid_goods_high = firm.F(K_high,get_aggregate_labor(hh)) - C_high - firm.δ * K_high - get_G(r_star_high,  w_high, L, A′_high, govt_high)
resid_goods_low = firm.F(K_low,get_aggregate_labor(hh)) - C_low - firm.δ * K_low - get_G(r_star_low,  w_low, L, A′_low, govt_low)

println("Stock of capital:") 
println("High taxes = $(solve_firm(firm,r_star_high)[1]*get_aggregate_labor(hh))")
println("Low taxes = $(solve_firm(firm,r_star_low)[1]*get_aggregate_labor(hh))")



# plot wealth distribution 
    plot(a_grid.a_vec,  λ_a_high, label = "High taxes",alpha = 0.5, color = :red, linestyle = :dash,lw=3)
    plot!(a_grid.a_vec, λ_a_low,  label = "Low taxes",alpha = 0.5,  color = :blue, lw=3, title = "Marginal distribution of assets")


