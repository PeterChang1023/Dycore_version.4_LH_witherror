function HS_Forcing!(atmo_data::Atmo_Data, Δt::Int64, sinθ::Array{Float64, 1}, grid_u::Array{Float64, 3}, grid_v::Array{Float64, 3}, grid_p_half::Array{Float64, 3}, grid_p_full::Array{Float64, 3}, grid_t::Array{Float64, 3}, 
  grid_δu::Array{Float64, 3}, grid_δv::Array{Float64, 3}, grid_t_eq::Array{Float64, 3}, grid_δt::Array{Float64, 3}, physcis_params::Dict{String, Float64}, grid_tracers_c::Array{Float64, 3})
  # reset grid_δu, grid_δv, grid_t_eq, grid_δt
  
  #todo delete
  #@info "HS_Forcing", sum(abs.(grid_u)), sum(abs.(grid_v)), sum(abs.(grid_t)) , sum(abs.(grid_p_half)), sum(abs.(grid_p_full))

  σ_b = physcis_params["σ_b"]  
  k_f = physcis_params["k_f"]  #day^{-1}
  k_a = physcis_params["k_a"]  #day^{-1}
  k_s = physcis_params["k_s"]  #day^{-1}
  ΔT_y = physcis_params["ΔT_y"] #K
  Δθ_z = physcis_params["Δθ_z"] #K

  # rayleigh damping of wind components near the surface

  Rayleigh_Damping!(atmo_data, grid_p_half, grid_p_full, grid_u, grid_v, grid_δu, grid_δv, 
                    σ_b, k_f)

  #todo 
  grid_δt .= 0.0
  do_conserve_energy = true
  if (do_conserve_energy) 
    cp_air = atmo_data.cp_air
    grid_δt .= -((grid_u + 0.5*grid_δu*Δt).*grid_δu + (grid_v + 0.5*grid_δv*Δt).*grid_δv)/cp_air
 end
  
  
  #thermal forcing for held & suarez (1994) benchmark calculation
  
  Newtonian_Damping!(atmo_data, sinθ, grid_p_half, grid_p_full, grid_t, grid_t_eq, grid_δt, σ_b, k_a, k_s, ΔT_y, Δθ_z, grid_tracers_c)

end 



function Newtonian_Damping!(atmo_data::Atmo_Data, sinθ::Array{Float64,1}, grid_p_half::Array{Float64,3}, grid_p_full::Array{Float64,3}, grid_t::Array{Float64,3}, 
  grid_t_eq::Array{Float64,3}, grid_δt::Array{Float64,3}, 
  σ_b::Float64, k_a::Float64, k_s::Float64, ΔT_y::Float64, Δθ_z::Float64, grid_tracers_c::Array{Float64,3})
  
  #routine to compute thermal forcing for held & suarez (1994)
  
  nλ, nθ, nd = size(grid_δt)
  
  day_to_sec = 86400
  t_zero, t_strat = 315.0, 200.0
  k_a, k_s = k_a/day_to_sec,  k_s/day_to_sec#s^-1
  
  p_ref = 1.0e5

  grid_ps = grid_p_half[:,:,nd+1]
  
  sinθ_2d = repeat(sinθ', nλ, 1)
  sinθ2_2d = sinθ_2d .* sinθ_2d
  cosθ2_2d = 1.0  .- sinθ2_2d
  cosθ4_2d = cosθ2_2d .* cosθ2_2d
  
  kappa = atmo_data.kappa
  grid_p_norm  = zeros(Float64, nλ, nθ)
  σ = zeros(Float64, nλ, nθ)
  k_t = zeros(Float64, nλ, nθ)
  # compute equilibrium temperature (grid_t_eq)
  
  for  k = 1:nd
    grid_p_norm .= grid_p_full[:,:,k]/p_ref

    grid_t_eq[:,:,k] .= (t_zero .- ΔT_y*sinθ2_2d  .- Δθ_z*cosθ2_2d.*log.(grid_p_norm)) .*  grid_p_norm.^kappa
    grid_t_eq[:,:,k] .= max.( t_strat, grid_t_eq[:,:,k])
  end

  # #######
  # #debug
  # θ_deg = asin.(sinθ)*180/pi
  # X = repeat(θ_deg, 1, nd)
  # Y = grid_ps[1,:,1] .\ grid_p_full[1,:,:]
  # PyPlot.contourf(X, Y, grid_t_eq[1,:,:], levels = Array(LinRange(200, 310, 12)))
  # PyPlot.colorbar()
  # PyPlot.gca().invert_yaxis()
  # ######
  ### By CJY4
  es  = zeros(size(grid_tracers_c))
  extra_e = zeros(size(grid_tracers_c))
  qvs = zeros(size(grid_tracers_c))
  cp  = 1004.
  Lv = 2.5*10^6.
  Rv = 461.
  one_array  = zeros(size(grid_t))
  one_array .= 1.0
    
  new_array  = zeros(size(grid_t))
  new_array .= 273.15
  initial_RH  = 1.0
  es[:,:,:] .= 6.11*exp.(Lv/Rv*(one_array[:,:,:] ./ new_array[:,:,:]-one_array[:,:,:] ./ grid_t[:,:,:])) * initial_RH
  for k = 1:nd
        for j=1:nθ
            for i=1:nλ
                if grid_tracers_c[i,j,k] >= es[i,j,k]
                    #@info grid_tracers_c[i,j,k], es[i,j,k], extra_e[i,j,k]
                    extra_e[i,j,k] = grid_tracers_c[i,j,k] - es[i,j,k]
                    grid_tracers_c[i,j,k] -= extra_e[i,j,k]
                    
                end
            end
        end
  end
  ###
  for k=1:nd
    σ .= grid_p_full[:,:,k]./grid_ps

    #todo
    @assert(maximum(σ .- σ_b)/(1.0 - σ_b) <= 1.0)
    
    k_t .= k_a .+ (k_s - k_a)/(1.0-σ_b) * max.(0.0, σ .- σ_b) .* cosθ4_2d
        
    ### By CJY4
    #rad_plus_latent_heat .= -k_t .* (grid_t[:,:,k]  - grid_t_eq[:,:,k]) + (pre[:,:,k] * Lv / cp)

    grid_δt[:,:,k] .-= k_t .* (grid_t[:,:,k]  - grid_t_eq[:,:,k]) 
    grid_δt[:,:,k] .+= (extra_e[:,:,k]*0.622 ./(grid_p_full[:,:,k] / 100.) .* Lv ./ cp) 
    
    
  end
  @info (maximum(grid_tracers_c[:,:,:])), maximum(grid_p_full[:,:,:]/100.)
  #@info (maximum((extra_e[:,:,:].*0.622 ./(grid_p_full[:,:,:]./100) .* Lv ./ cp) ))
  #@info maximum(grid_tracers_c), minimum(es), maximum(extra_e), maximum(grid_t)
  # @info "teq = ", sum(grid_t_eq)
  # @info sum(grid_t), sum(grid_t_eq), sum(k_t[1,:,:]), sum(grid_δt)


end 

!#######################################################################

function Rayleigh_Damping!(atmo_data::Atmo_Data, grid_p_half::Array{Float64,3}, grid_p_full::Array{Float64,3}, 
  grid_u::Array{Float64,3}, grid_v::Array{Float64,3}, 
  grid_δu::Array{Float64, 3}, grid_δv::Array{Float64, 3}, 
  σ_b::Float64, k_f::Float64)
  
  #rayleigh damping of wind components near surface
  nλ, nθ, nd = size(grid_u)
  day_to_sec = 86400

  k_f = k_f/day_to_sec #s^-1

  grid_ps = grid_p_half[:,:,nd+1]
  
  
  σ = zeros(Float64, nλ, nθ)
  k_v = zeros(Float64, nλ, nθ)
  for k = 1:nd
    
    σ .= grid_p_full[:,:,k]./grid_ps

    #todo
    @assert(maximum(σ .- σ_b)/(1.0 - σ_b) <= 1.0)
    
    k_v = k_f/(1.0 - σ_b) * max.(0, σ .- σ_b)
    grid_δu[:,:,k]  .= -k_v .* grid_u[:,:,k]
    grid_δv[:,:,k]  .= -k_v .* grid_v[:,:,k]
    
  end
  
end 



