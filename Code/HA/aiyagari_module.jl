
### MODULE 
module Aiyagari
using Distributions, QuantEcon, IterTools, Plots, Optim, Interpolations, LinearAlgebra, Inequality, Statistics, ColorSchemes,PrettyTables, Roots, Revise

export create_firm, solve_firm, get_aggregate_labor
export create_government, get_G
export create_grid, create_ha_block, getW!, vfi, T_operator
export get_q, get_transition, stationary_distribution_hh, show_statistics_hugget, solve_hh_block


function create_grid(ha_block;N_a=70,a_max=40);
    a_min = -ha_block.ϕ

    a_vec= collect(range(a_min, a_max, length=N_a))
    return (; a_vec, N_a, a_min, a_max)
end

function create_firm(;
    α = 1/3, # capital share
    δ = 0.1, # depreciation rate
    Z = 1, # productivity
    F   =  (K,L) ->  Z * K^α * L^(1-α), # production function
    F_K =  (K,L) ->  α * Z * K^(α-1) * L^(1-α), # marginal product of capital
    F_L =  (K,L) -> (1-α) * Z * K^α * L^(-α)) # marginal product of labor
    return (; α, δ, F, F_K, F_L)
end

function solve_firm(firm,r)
    (; α, δ, F, F_K, F_L) = firm
    K_L = (α/(r+δ))^(1/(1-α)) # capital to output ratio
    w = F_L(K_L,1) # wage
    return K_L, w
end

function get_aggregate_labor(ha_block)
    (; ρ_z, ν_z, γ,  u, ϕ, β, N_z, mc_z, z_vec, P_z , λ_z) = ha_block
    L =  sum(z_vec .* λ_z)
    return L
end

function create_government(;
    τ_w = 0.1, # labor tax
    τ_r = 0.1, # capital tax
    d = 0.0, # lump-sum transfer
    B = 0.0) # debt
    return (; τ_w, τ_r, d, B)
end

function get_G(r, w, L, A, government)
    (; τ_w, τ_r, d, B) = government
    G = τ_w * w * L + τ_r * r * A - d - r * B
    return G
end


function create_ha_block(; 
    ρ_z=0.96, # log of productivity persistence
    ν_z=sqrt(0.025), # log of productivity volatility
    γ = 2, # curvature parameter of utility function
    u = γ == 1 ? x -> log(x) : x -> (x^(1 - γ) - 1) / (1 - γ), # utility function
    ϕ = 0.0, # borrowing constraint
    β= 0.98, # discount factor
    N_z= 9, # grid size for Tauchen
    mc_z = tauchen(N_z, ρ_z, ν_z, 0),
    λ_z = stationary_distributions(mc_z)[1],
    P_z = mc_z.p, # transition matrix
    z_vec = exp.(mc_z.state_values) / sum(exp.(mc_z.state_values) .* λ_z)) # normalize so that mean is 1
    return (; ρ_z, ν_z, γ,  u, ϕ, β, N_z, mc_z, z_vec, P_z , λ_z)
end

function getW!(W,V,a_vec,z_vec,P)
    for (i, a′) in enumerate(a_vec)
        for (j, z) in enumerate(z_vec)
            W[i,j] = sum(P[j,j′] .* V[i,j′] for j′ in eachindex(z_vec))
        end
    end
return W
end

function vfi(ha_block,prices,grid;tol = 1e-7, maxiter = 2000)
    (; ρ_z, ν_z, γ,  u, ϕ, β, N_z, mc_z, z_vec, P_z , λ_z) = ha_block
    (; a_vec, N_a, a_min, a_max) = grid
    (; r, w, τ_w, τ_r, d) = prices
    
    V = zeros(N_a,N_z) # initialize value function
    W = similar(V) # initialize W
    σ_a′ = similar(V) # initialize policy function
    σ_c = similar(V) # initialize policy function
    
    # get W
    getW!(W,V,a_vec,z_vec,P_z)
    error = 1 + tol
    iter = 1;
    while error > tol && iter < maxiter
        V_new, (σ_a′,σ_c) = T_operator(W, ha_block,prices, grid)
        error = maximum(abs.(V_new .- V))
        getW!(W,V_new,a_vec,z_vec,P_z)
        V = V_new
        
        iter += 1
    end
    policies = (a′ = σ_a′, c= σ_c)

    return V, policies, error, iter
end


function T_operator(W, ha_block, prices, grid; tol = 1e-10)
    
    (; ρ_z, ν_z, γ,  u, ϕ, β, N_z, mc_z, z_vec, P_z , λ_z) = ha_block
    (; a_vec, N_a, a_min, a_max) = grid
    (; r, w, τ_w, τ_r, d) = prices
    W_hat = linear_interpolation((a_vec,z_vec),W,extrapolation_bc=Linear()) # linear interpolation of W

    TW = similar(W)
    σ_a′ = similar(W)
    σ_c = similar(W)
    for (i, a) in enumerate(a_vec)
        for (j, z) in enumerate(z_vec)
            # solve maximization for each point in (a,z), using a itself as initial condition.
            a_max = (1+(1-τ_r) * r) * a + (1 - τ_w ) * w*z + d
            results = optimize(a′ -> -u( (1+(1-τ_r) * r) * a + (1 - τ_w ) * w*z + d - a′ ) - β * W_hat(a′,z), a_min, a_max ,GoldenSection(),g_tol=1e-8)
            TW[i,j]   = -Optim.minimum(results)
            σ_a′[i,j] = Optim.minimizer(results)
            σ_c[i,j]  = (1+(1-τ_r) * r) * a + (1 - τ_w ) * w*z + d- σ_a′[i,j]
        end
    end
    policies = (a′ = σ_a′, c= σ_c)
    return TW, policies 
end

# function from https://discourse.julialang.org/t/findnearest-function/4143/4
function closest_index(a::Vector,x::Real)

    if isempty(a) == true
      error("xGrid is empty in function closest_index.")
    end

    if isnan(x) == true
      error("val is NaN in function closest_index.")
    end

   idx = searchsortedfirst(a,x)
   if (idx==1); return idx; end
   if (idx>length(a)); return length(a); end
   if (a[idx]==x); return idx; end
   if (abs(a[idx]-x) < abs(a[idx-1]-x))
      return idx
   else
      return idx-1
   end
end

function closest_value_and_index(xGrid::Vector, val::Real)

   # get index
   ibest = closest_index(xGrid, val)

   # Return best value on grid, and the corresponding index
   return xGrid[ibest], ibest

end

function get_q(w_grid,w) # code from https://julienpascal.github.io/post/young_2010/


    w_min = minimum(w_grid)
    w_max = maximum(w_grid)
    nW = length(w_grid)
    
    q = zeros(length(w_grid),length(w_grid))
    
    
        for (wIndexTrue, w_temp) in enumerate(w_grid)
    
    
            # Project true value on the grid:
            (wValue_proj, wIndex_proj) = closest_value_and_index(w_grid, w[wIndexTrue])
    
            # To store the location of the value below and above the true value:
            wIndex_below = 0
            wIndex_above = 0
    
            # If the true value is above the projection
            if w[wIndexTrue] >= wValue_proj
                wIndex_below = wIndex_proj
                wIndex_above = wIndex_proj + 1
            # If the true value is below the projection
            elseif w[wIndexTrue] < wValue_proj
                wIndex_below = wIndex_proj -1
                wIndex_above = wIndex_proj
            end
    
            # Boundary cases
            if wIndex_proj == 1
                wIndex_below = 1
                wIndex_above = 2
            elseif wIndex_proj == nW
                wIndex_below = nW - 1
                wIndex_above = nW
            end
    
            # Special case 1: w < w_min
            if w[wIndexTrue] <= w_min
                p = 1
            elseif w[wIndexTrue] >= w_max
            # Special case 2: w > w_max
                p = 0
            else
                p = 1.0 - ((w[wIndexTrue] - w_grid[wIndex_below])/(w_grid[wIndex_above] - w_grid[wIndex_below]))
                p = min(1.0, max(0.0, p))
            end
    
        
    
        q[wIndexTrue, wIndex_below] = p
        q[wIndexTrue, wIndex_above] = 1.0 - p
    end
    
        return q
    end

function get_transition(ha_block, policies, grid)

    (; ρ_z, ν_z, γ,  u, ϕ, β, N_z, mc_z, z_vec, P_z , λ_z) = ha_block
    (; a_vec, N_a, a_min, a_max) = grid
    
    Q = zeros(N_a * N_z,N_a * N_z)
    
    
    for (j, z) in enumerate(z_vec)
            for (j′, z′) in enumerate(z_vec)
                    Q[(j-1)*N_a+1:j*N_a,(j′-1)*N_a+1:j′*N_a] = get_q(a_vec,policies.a′[:,j]) .* P_z[j,j′]
            end
    end
    
    
    return Q
    end

function stationary_distribution_hh(ha_block, policies, grid)

    Q = get_transition(ha_block, policies, grid)

    N_a = grid.N_a
    N_z = ha_block.N_z    
    z_vec = ha_block.z_vec

    λ_vector = (Q^10000)[1,:]
    λ = zeros(N_a, N_z)

    for (j, z) in enumerate(z_vec)
        for (j, z′) in enumerate(z_vec)
            λ[:,j] = λ_vector[(j-1)*N_a+1:j*N_a]
        end
    end

    λ_a = sum(λ,dims=2)
    λ_z = sum(λ,dims=1)'
    return λ, λ_vector, λ_a, λ_z
end 

function show_statistics(ha_block,grid,λ_a,λ_z)
# warning - this can be misleading if we allow for negative values!
lorenz_a_pop,lorenz_a_share=lorenz_curve(grid.a_vec,vec(λ_a))
lorenz_z_pop,lorenz_z_share=lorenz_curve(ha_block.z_vec,vec(λ_z))



lorenz_a = LinearInterpolation(lorenz_a_pop, lorenz_a_share);
lorenz_z = LinearInterpolation(lorenz_z_pop, lorenz_z_share);


header = (["", "Assets", "Income"])

data = [           
                     "Bottom 50% share"         lorenz_a(0.5)        lorenz_z(0.5)    ;
                     "Top 10% share"            1-lorenz_a(0.9)         1-lorenz_z(0.9)     ;
                     "Top 1% share"             1-lorenz_a(0.99)        1-lorenz_z(0.99)    ;  
                     "Gini Coefficient"      wgini(grid.a_vec,vec(max.(0,λ_a)))      wgini(ha_block.z_vec,vec(max.(0.0,λ_z)))    ;]

return pretty_table(data;header=header,formatters=ft_printf("%5.3f",2:3))
end

function solve_hh_block(grid,ha_block,prices)
    V, policies, error, iter = vfi(ha_block,prices,grid)
    λ, λ_vector, λ_a, λ_z = stationary_distribution_hh(ha_block, policies, grid)
    A′ = sum(λ .* policies.a′)
    C  = sum(λ .* policies.c)


    return V, policies, error, iter, λ, λ_vector, λ_a, λ_z, A′, C
end



end # module ends

