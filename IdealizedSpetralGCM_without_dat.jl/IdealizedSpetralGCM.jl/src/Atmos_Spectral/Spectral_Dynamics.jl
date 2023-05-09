export Compute_Corrections_Init, Compute_Corrections!, Four_In_One!, Spectral_Dynamics!, Get_Topography!, Spectral_Initialize_Fields!, Spectral_Dynamics_Physics!, Atmosphere_Update!



function Compute_Corrections_Init(vert_coord::Vert_Coordinate, mesh::Spectral_Spherical_Mesh, atmo_data::Atmo_Data,
  grid_u_p::Array{Float64, 3}, grid_v_p::Array{Float64, 3}, grid_ps_p::Array{Float64, 3}, grid_t_p::Array{Float64, 3}, 
  grid_δu::Array{Float64, 3}, grid_δv::Array{Float64, 3}, grid_δt::Array{Float64, 3},  
  Δt::Int64, grid_energy_temp::Array{Float64, 3}, grid_tracers_p::Array{Float64, 3}, grid_tracers_c::Array{Float64, 3}, grid_δtracers::Array{Float64,3}, grid_tracers_all::Array{Float64,3})
  
  do_mass_correction, do_energy_correction, do_water_correction = atmo_data.do_mass_correction, atmo_data.do_energy_correction, atmo_data.do_water_correction
  
  if (do_mass_correction) 
    mean_ps_p = Area_Weighted_Global_Mean(mesh, grid_ps_p)
  end
  
  if (do_energy_correction) 
    # due to dissipation introduced by the forcing
    cp_air, grav = atmo_data.cp_air, atmo_data.grav
    grid_energy_temp  .=  0.5*((grid_u_p + Δt*grid_δu).^2 + (grid_v_p + Δt*grid_δv).^2) + cp_air*(grid_t_p + Δt*grid_δt)
    mean_energy_p = Mass_Weighted_Global_Integral(vert_coord, mesh, atmo_data, grid_energy_temp, grid_ps_p)
  end
  
  if (do_water_correction)
    # By CJY3
    grid_tracers_all .= grid_tracers_p .+ Δt*grid_δtracers
    sum_tracers_p    = grav*Mass_Weighted_Global_Integral(vert_coord, mesh, atmo_data, grid_tracers_all, grid_ps_p)
  end
  
  return mean_ps_p, mean_energy_p, sum_tracers_p
end 

function Compute_Corrections!(vert_coord::Vert_Coordinate, mesh::Spectral_Spherical_Mesh, atmo_data::Atmo_Data,
  mean_ps_p::Float64, mean_energy_p::Float64, 
  grid_u_n::Array{Float64, 3}, grid_v_n::Array{Float64, 3},
  grid_energy_temp::Array{Float64, 3},
  grid_ps_n::Array{Float64, 3}, spe_lnps_n::Array{ComplexF64, 3}, 
  grid_t_n::Array{Float64, 3}, spe_t_n::Array{ComplexF64, 3}, grid_tracers_p::Array{Float64, 3}, grid_tracers_c::Array{Float64, 3}, grid_tracers_n::Array{Float64, 3},grid_δtracers::Array{Float64,3}, grid_tracers_all::Array{Float64,3},sum_tracers_p::Float64, spe_tracers_n::Array{ComplexF64, 3})
  
  do_mass_correction, do_energy_correction, do_water_correction = atmo_data.do_mass_correction, atmo_data.do_energy_correction, atmo_data.do_water_correction
  
  
  if (do_mass_correction) 
    mean_ps_n = Area_Weighted_Global_Mean(mesh, grid_ps_n)
    mass_correction_factor = mean_ps_p/mean_ps_n
    grid_ps_n .*= mass_correction_factor
    #P00 = 1 
    spe_lnps_n[1,1,1] += log(mass_correction_factor)
  end
  
  if (do_energy_correction) 
    cp_air, grav = atmo_data.cp_air, atmo_data.grav
    
    grid_energy_temp .=  0.5*(grid_u_n.^2 + grid_v_n.^2) + cp_air*grid_t_n
    mean_energy_n = Mass_Weighted_Global_Integral(vert_coord, mesh, atmo_data, grid_energy_temp, grid_ps_n)
    
    temperature_correction = grav*(mean_energy_p - mean_energy_n)/(cp_air*mean_ps_p)
    #@info grav, mean_energy_p , mean_energy_n, cp_air, mean_ps_p
    grid_t_n .+= temperature_correction
    spe_t_n[1,1,:] .+= temperature_correction
    
  end

  #@info mean_ps_p, mean_energy_p, mass_correction_factor, temperature_correction
  # error(6868)

   nλ = mesh.nλ
   nθ = mesh.nθ
   nd = mesh.nd
   grav = atmo_data.grav
  if (do_water_correction) 
    #while tracers_correction != 0 
    
    std = zeros(size(grid_tracers_n))
    total  = zeros(size(grid_tracers_n))
    for k in collect(1:nd)
        for j in collect(1:nθ)
            for i in collect(1:nλ)
                if (grid_tracers_n[i,j,k] .< std[i,j,k]) 
                    total[i,j,k] = (grid_tracers_n[i,j,k])
                    #print((total[i,j,k]))
                    grid_tracers_n[i,j,k] = max.(std[i,j,k], grid_tracers_n[i,j,k])
                else
                    grid_tracers_n[i,j,k] = grid_tracers_n[i,j,k]
                end
            end
        end
    end
    grid_tracers_all .= grid_tracers_n 
    sum_tracers_n    = grav*Mass_Weighted_Global_Integral(vert_coord, mesh, atmo_data, grid_tracers_all, grid_ps_n)

    #tracers_correction = (sum_tracers_p - (sum_tracers_n-abs(grav*Mass_Weighted_Global_Integral(vert_coord, mesh, atmo_data, total, grid_ps_n))))
    tracers_correction = (sum_tracers_p - (sum_tracers_n))
    #tracers_correction = (sum_tracers_p / (sum_tracers_n-abs(Mass_Weighted_Global_Integral(vert_coord, mesh, atmo_data, total, grid_ps_n))))
    # sum(total) is negative
    # meanly distribute tracers_correction
    for k in collect(1:nd)
        for j in collect(1:nθ)
            for i in collect(1:nλ)
                grid_tracers_n[i,j,k] += (grid_tracers_n[i,j,k]/sum_tracers_n) * tracers_correction
                 #grid_tracers_n[i,j,k] *=  tracers_correction
            end
        end
    end
    
    grid_tracers_all .= grid_tracers_n 
    sum_tracers_n    = grav*Mass_Weighted_Global_Integral(vert_coord, mesh, atmo_data, grid_tracers_all, grid_ps_n)
    #@info   (sum_tracers_p - sum_tracers_n)   # mean_tracers_p mean_tracers_n
    #print(sum(total))
    #end
    #error("water correction has not implemented")
    #grid_tracers_n .+= tracers_correction
  end
  return (sum_tracers_p - sum_tracers_n)
end 





"""
compute vertical mass flux and velocity 
grid_M_half[:,:,k+1] = downward mass flux/per unit area across the K+1/2
grid_w_full[:,:,k]   = dp/dt vertical velocity 

update residuals
grid_δps[:,:,k]  += -∑_{r=1}^nd Dr = -∑_{r=1}^nd ∇(vrΔp_r)
grid_δt[:,:,k]  += κTw/p 
(grid_δu[:,:,k], grid_δv[:,:,k]) -= RT ∇p/p 

!  cell boundary. This is the "vertical velocity" in the hybrid coordinate system.
!  When vertical coordinate is pure sigma: grid_M_half = grid_ps*d(sigma)/dt
"""

function Four_In_One!(vert_coord::Vert_Coordinate, atmo_data::Atmo_Data, 
  grid_div::Array{Float64,3}, grid_u::Array{Float64,3}, grid_v::Array{Float64,3}, 
  grid_ps::Array{Float64,3},  grid_Δp::Array{Float64,3}, grid_lnp_half::Array{Float64,3}, grid_lnp_full::Array{Float64,3}, grid_p_full::Array{Float64,3},
  grid_dλ_ps::Array{Float64,3}, grid_dθ_ps::Array{Float64,3}, 
  grid_t::Array{Float64,3}, 
  grid_M_half::Array{Float64,3}, grid_w_full::Array{Float64,3}, 
  grid_δu::Array{Float64,3}, grid_δv::Array{Float64,3}, grid_δps::Array{Float64,3}, grid_δt::Array{Float64,3}, grid_δtracers::Array{Float64,3})
  
  rdgas, cp_air = atmo_data.rdgas, atmo_data.cp_air
  nd, bk = vert_coord.nd, vert_coord.bk
  Δak, Δbk = vert_coord.Δak, vert_coord.Δbk
  vert_difference_option = vert_coord.vert_difference_option
  
  kappa = rdgas / cp_air
  
  # dmean_tot = ∇ ∑_{k=1}^{nd} vk Δp_k = ∑_{k=1}^{nd} Dk
  nλ, nθ, _ = size(grid_ps)
  dmean_tot = zeros(Float64, nλ, nθ)
  Δlnp_p = zeros(Float64, nλ, nθ)
  Δlnp_m = zeros(Float64, nλ, nθ)
  Δlnp = zeros(Float64, nλ, nθ)
  x1 = zeros(Float64, nλ, nθ)
  dlnp_dλ = zeros(Float64, nλ, nθ)
  dlnp_dθ = zeros(Float64, nλ, nθ)
  dmean = zeros(Float64, nλ, nθ)
  x5 = zeros(Float64, nλ, nθ)
    
  if (vert_difference_option == "simmons_and_burridge") 
    for k = 1:nd
      Δp = grid_Δp[:,:,k]
      
      Δlnp_p .= grid_lnp_half[:,:,k + 1] - grid_lnp_full[:,:,k]
      Δlnp_m .= grid_lnp_full[:,:,k]   - grid_lnp_half[:,:,k]
      Δlnp   .= grid_lnp_half[:,:,k + 1] - grid_lnp_half[:,:,k]
      
      # angular momentum conservation 
      #    ∇p_k/p =  [(lnp_k - lnp_{k-1/2})∇p_{k-1/2} + (lnp_{k+1/2} - lnp_k)∇p_{k+1/2}]/Δpk
      #         =  [(lnp_k - lnp_{k-1/2})B_{k-1/2} + (lnp_{k+1/2} - lnp_k)B_{k+1/2}]/Δpk * ∇ps
      #         =  x1 * ∇ps
      x1 .= (bk[k] * Δlnp_m + bk[k + 1] * Δlnp_p ) ./ Δp
      
      dlnp_dλ .= x1 .* grid_dλ_ps[:,:,1]
      dlnp_dθ .= x1 .* grid_dθ_ps[:,:,1]
      
      
      
      # (grid_δu, grid_δv) -= RT ∇p/p 
      grid_δu[:,:,k] .-=  rdgas * grid_t[:,:,k] .* dlnp_dλ
      grid_δv[:,:,k] .-=  rdgas * grid_t[:,:,k] .* dlnp_dθ
      
      # dmean = ∇ (vk Δp_k) =  divk Δp_k + vk  Δbk[k] ∇ p_s
      dmean .= grid_div[:,:,k] .* Δp + Δbk[k] * (grid_u[:,:,k] .* grid_dλ_ps[:,:,1] + grid_v[:,:,k] .* grid_dθ_ps[:,:,1])
      
  
      # energy conservation for temperature
      # w/p = dlnp/dt = ∂lnp/∂t + dσ ∂lnp/∂σ + v∇lnp
      # dσ ∂ξ_k/∂σ = [M_{k+1/2}(ξ_k+1/2 - ξ_k) + M_{k-1/2}(ξ_k - ξ_k-1/2)]/Δp_k
      # weight the same way (TODO)
      # vertical advection operator (M is the downward speed)
      # dσ ∂lnp_k/∂σ = [M_{k+1/2}(lnp_k+1/2 - lnp_k) + M_{k-1/2}(lnp_k - lnp_k-1/2)]/Δp_k
      # ∂lnp/∂t = 1/p ∂p/∂t = [∂p/∂t_{k+1/2}(lnp_k+1/2 - lnp_k) + ∂p/∂t_{k-1/2}(lnp_k - lnp_k-1/2)]/Δp_k
      # As we know
      # ∂p/∂t_{k+1/2} = -∑_{r=1}^k Dr - M_{k+1/2}
      
      # ∂lnp/∂t + dσ ∂lnp/∂σ =  [(-∑_{r=1}^k Dr)(lnp_k+1/2 - lnp_k) + (-∑_{r=1}^{k-1} Dr)(lnp_k - lnp_k-1/2)]/Δp_k
      #                      = -[(∑_{r=1}^{k-1} Dr)(lnp_k+1/2 - lnp_k-1/2) + D_k(lnp_k+1/2 - lnp_k)]/Δp_k
      
      x5 .= -(dmean_tot .* Δlnp + dmean .* Δlnp_p) ./ Δp .+ grid_u[:,:,k] .* dlnp_dλ + grid_v[:,:,k] .* dlnp_dθ
      # grid_δt += κT w/p
      grid_δt[:,:,k] .+=  kappa * grid_t[:,:,k] .* x5
      # grid_w_full = w
      grid_w_full[:,:,k] .= x5 .* grid_p_full[:,:,k]
      # update dmean_tot to ∑_{r=1}^k ∇(vrΔp_r)
      dmean_tot .+= dmean
      # M_{k+1/2} = -∑_{r=1}^k ∇(vrΔp_r) - B_{k+1/2}∂ps/∂t
      grid_M_half[:,:,k + 1] .= -dmean_tot
    end
    
  else
    error("vert_difference_option ", vert_difference_option, " is not a valid value for option")
    
  end
  # ∂ps/∂t = -∑_{r=1}^nd ∇(vrΔp_r) = -dmean_tot
  grid_δps[:,:,1] .-= dmean_tot
  
  for k = 1:nd-1
    # M_{k+1/2} = -∑_{r=1}^k ∇(vrΔp_r) - B_{k+1/2}∂ps/∂t
    grid_M_half[:,:,k+1] .+= dmean_tot * bk[k+1]
  end
  
  grid_M_half[:,:,1] .= 0.0
  grid_M_half[:,:,nd + 1] .= 0.0


  
end 



"""
The governing equations are
∂div/∂t = ∇ × (A, B) - ∇^2E := f^d                    
∂lnps/∂t= (-∑_k div_k Δp_k + v_k ∇ Δp_k)/ps := f^p    
∂T/∂t = -(u,v)∇T - dσ∂T∂σ + κTw/p + J:= f^t           
Φ = f^Φ                                               

implicit part: -∇^2Φ - ∇(RT∇lnp) ≈ I^d = -∇^2(γT + H2 ps_ref lnps) - ∇^2 H1 ps_ref lnps, here RT∇lnp ≈  H1 ps_ref ∇lnps
implicit part:  f^p              ≈ I^p = -ν div / ps_ref
implicit part:  - dσ∂T∂σ + κTw/p ≈ I^t = -τ div  
implicit part:  f^Φ              ≈ I^Φ = γT + H2 ps_ref lnps 

We have 
δdiv = f^d - I^d + I^d
δlnps = f^p - I^p + I^p
δT = f^t - I^t + I^t

"""

function Spectral_Dynamics!(mesh::Spectral_Spherical_Mesh,  vert_coord::Vert_Coordinate, 
  atmo_data::Atmo_Data, dyn_data::Dyn_Data, 
  semi_implicit::Semi_Implicit_Solver)
  
  
  # spectral equation quantities
  spe_lnps_p, spe_lnps_c, spe_lnps_n, spe_δlnps = dyn_data.spe_lnps_p, dyn_data.spe_lnps_c, dyn_data.spe_lnps_n, dyn_data.spe_δlnps
  spe_vor_p, spe_vor_c, spe_vor_n, spe_δvor = dyn_data.spe_vor_p, dyn_data.spe_vor_c, dyn_data.spe_vor_n, dyn_data.spe_δvor
  spe_div_p, spe_div_c, spe_div_n, spe_δdiv = dyn_data.spe_div_p, dyn_data.spe_div_c, dyn_data.spe_div_n, dyn_data.spe_δdiv
  spe_t_p, spe_t_c, spe_t_n, spe_δt = dyn_data.spe_t_p, dyn_data.spe_t_c, dyn_data.spe_t_n, dyn_data.spe_δt
  
  # grid quantities
  grid_u_p, grid_u, grid_u_n = dyn_data.grid_u_p, dyn_data.grid_u_c, dyn_data.grid_u_n
  grid_v_p, grid_v, grid_v_n = dyn_data.grid_v_p, dyn_data.grid_v_c, dyn_data.grid_v_n
  grid_ps_p, grid_ps, grid_ps_n = dyn_data.grid_ps_p, dyn_data.grid_ps_c, dyn_data.grid_ps_n
  grid_t_p, grid_t, grid_t_n = dyn_data.grid_t_p, dyn_data.grid_t_c, dyn_data.grid_t_n


  # related quanties
  grid_p_half, grid_lnp_half, grid_p_full, grid_lnp_full = dyn_data.grid_p_half, dyn_data.grid_lnp_half, dyn_data.grid_p_full, dyn_data.grid_lnp_full
  grid_dλ_ps, grid_dθ_ps = dyn_data.grid_dλ_ps, dyn_data.grid_dθ_ps
  grid_lnps = dyn_data.grid_lnps
  
  grid_div, grid_absvor, grid_vor = dyn_data.grid_div, dyn_data.grid_absvor, dyn_data.grid_vor
  grid_w_full, grid_M_half = dyn_data.grid_w_full, dyn_data.grid_M_half
  grid_geopots, grid_geopot_full, grid_geopot_half = dyn_data.grid_geopots, dyn_data.grid_geopot_full, dyn_data.grid_geopot_half
  
  grid_energy_full, spe_energy = dyn_data.grid_energy_full, dyn_data.spe_energy
  
  # By CJY2
  spe_tracers_n = dyn_data.spe_tracers_n
  spe_tracers_c = dyn_data.spe_tracers_c
  spe_tracers_p = dyn_data.spe_tracers_p 
    
  grid_tracers_n = dyn_data.grid_tracers_n
  grid_tracers_c = dyn_data.grid_tracers_c
  grid_tracers_p = dyn_data.grid_tracers_p 
    
  spe_δtracers   = dyn_data.spe_δtracers
  grid_δtracers  = dyn_data.grid_δtracers
  ### By CJY3 
  grid_tracers_full = dyn_data.grid_tracers_full

  ###
  ###
  # todo !!!!!!!!
  #  grid_q = grid_t
  
  # original 
  # pressure difference
  grid_Δp = dyn_data.grid_Δp
  # temporary variables
  grid_δQ = dyn_data.grid_d_full1
  ### By CJY4
  grid_δQ2 = dyn_data.grid_d_full2
  ###
    
  # incremental quantities
  grid_δu, grid_δv, grid_δps, grid_δlnps, grid_δt = dyn_data.grid_δu, dyn_data.grid_δv, dyn_data.grid_δps, dyn_data.grid_δlnps, dyn_data.grid_δt
  

  integrator = semi_implicit.integrator
  Δt = Get_Δt(integrator)

  mean_ps_p, mean_energy_p, sum_tracers_p = Compute_Corrections_Init(vert_coord, mesh, atmo_data,
  grid_u_p, grid_v_p, grid_ps_p, grid_t_p, 
  grid_δu, grid_δv, grid_δt,  
  Δt, grid_energy_full, grid_tracers_p, grid_tracers_c, grid_δtracers, grid_tracers_full)
  
  # compute pressure based on grid_ps -> grid_p_half, grid_lnp_half, grid_p_full, grid_lnp_full 
  Pressure_Variables!(vert_coord, grid_ps, grid_p_half, grid_Δp, grid_lnp_half, grid_p_full, grid_lnp_full)
  
  # compute ∇ps = ∇lnps * ps
  Compute_Gradients!(mesh, spe_lnps_c,  grid_dλ_ps, grid_dθ_ps)
  grid_dλ_ps .*= grid_ps
  grid_dθ_ps .*= grid_ps


  
  # compute grid_M_half, grid_w_full, grid_δu, grid_δv, grid_δps, grid_δt, 
  # except the contributions from geopotential or vertical advection
  Four_In_One!(vert_coord, atmo_data, grid_div, grid_u, grid_v, grid_ps, 
  grid_Δp, grid_lnp_half, grid_lnp_full, grid_p_full,
  grid_dλ_ps, grid_dθ_ps, 
  grid_t, 
  grid_M_half, grid_w_full, grid_δu, grid_δv, grid_δps, grid_δt, grid_δtracers)

  Compute_Geopotential!(vert_coord, atmo_data, 
  grid_lnp_half, grid_lnp_full,  
  grid_t, 
  grid_geopots, grid_geopot_full, grid_geopot_half)
  

  grid_δlnps .= grid_δps ./ grid_ps
  Trans_Grid_To_Spherical!(mesh, grid_δlnps, spe_δlnps)
  

  
  # compute vertical advection, todo  finite volume method 
  Vert_Advection!(vert_coord, grid_u, grid_Δp, grid_M_half, Δt, vert_coord.vert_advect_scheme, grid_δQ)
  grid_δu  .+= grid_δQ
  Vert_Advection!(vert_coord, grid_v, grid_Δp, grid_M_half, Δt, vert_coord.vert_advect_scheme, grid_δQ)
  grid_δv  .+= grid_δQ
  Vert_Advection!(vert_coord, grid_t, grid_Δp, grid_M_half, Δt, vert_coord.vert_advect_scheme, grid_δQ)
  grid_δt  .+= grid_δQ
  # Vert_Advection!(vert_coord, grid_tracers_c, grid_Δp, grid_M_half, Δt, vert_coord.vert_advect_scheme, grid_δQ)
  # grid_δtracers  .+= grid_δQ

  Add_Horizontal_Advection!(mesh, spe_t_c, grid_u, grid_v, grid_δt)
  Trans_Grid_To_Spherical!(mesh, grid_δt, spe_δt)
  ### By CJY2
  # Add_Horizontal_Advection!(mesh, spe_tracers_c, grid_u, grid_v, grid_δtracers)
  # Trans_Grid_To_Spherical!(mesh, grid_δtracers, spe_δtracers)
  ###

  grid_absvor = dyn_data.grid_absvor
  Compute_Abs_Vor!(grid_vor, atmo_data.coriolis, grid_absvor)
  
  
  grid_δu .+=  grid_absvor .* grid_v
  grid_δv .-=  grid_absvor .* grid_u
  
  
  Vor_Div_From_Grid_UV!(mesh, grid_δu, grid_δv, spe_δvor, spe_δdiv)
  

  grid_energy_full .= grid_geopot_full .+ 0.5 * (grid_u.^2 + grid_v.^2)
  Trans_Grid_To_Spherical!(mesh, grid_energy_full, spe_energy)
  Apply_Laplacian!(mesh, spe_energy)
  spe_δdiv .-= spe_energy
  
  
  Implicit_Correction!(semi_implicit, vert_coord, atmo_data,
  spe_div_c, spe_div_p, spe_lnps_c, spe_lnps_p, spe_t_c, spe_t_p, 
  spe_δdiv, spe_δlnps, spe_δt)
  
  Compute_Spectral_Damping!(integrator, spe_vor_c, spe_vor_p, spe_δvor)
  Compute_Spectral_Damping!(integrator, spe_div_c, spe_div_p, spe_δdiv)
  Compute_Spectral_Damping!(integrator, spe_t_c, spe_t_p, spe_δt)
  ### By CJY2
  #Compute_Spectral_Damping!(integrator, spe_tracers_c, spe_tracers_p, spe_δtracers)
  ###
    
  Filtered_Leapfrog!(integrator, spe_δvor, spe_vor_p, spe_vor_c, spe_vor_n)
  Filtered_Leapfrog!(integrator, spe_δdiv, spe_div_p, spe_div_c, spe_div_n)
  Filtered_Leapfrog!(integrator, spe_δlnps, spe_lnps_p, spe_lnps_c, spe_lnps_n)
  Filtered_Leapfrog!(integrator, spe_δt, spe_t_p, spe_t_c, spe_t_n)
  ### By CJY2
  #Filtered_Leapfrog!(integrator, spe_δtracers, spe_tracers_p, spe_tracers_c, spe_tracers_n)
  ###

  
  Trans_Spherical_To_Grid!(mesh, spe_vor_n, grid_vor)
  Trans_Spherical_To_Grid!(mesh, spe_div_n, grid_div)
  UV_Grid_From_Vor_Div!(mesh, spe_vor_n, spe_div_n, grid_u_n, grid_v_n)
  Trans_Spherical_To_Grid!(mesh, spe_t_n, grid_t_n)
  ### By CJ2
  #Trans_Spherical_To_Grid!(mesh, spe_tracers_n, grid_tracers_n)
  ###
  Trans_Spherical_To_Grid!(mesh, spe_lnps_n, grid_lnps)
  grid_ps_n .= exp.(grid_lnps)

  ### By CJY4
  # tracer_spe update
    
  Compute_Corrections!(vert_coord, mesh, atmo_data, mean_ps_p, mean_energy_p, 
  grid_u_n, grid_v_n,
  grid_energy_full,
  grid_ps_n, spe_lnps_n, 
  grid_t_n, spe_t_n, grid_tracers_p, grid_tracers_c, grid_tracers_n, grid_δtracers, grid_tracers_full, sum_tracers_p, spe_tracers_n)
    
  

  """
  1. horizontal advection
  2. veritical advection
  3. hole filling
  4. leapfrog and time filter 
  5. back to grid
  """
  """
  Add_Horizontal_Advection!(mesh, spe_tracers_c, grid_u, grid_v, grid_δtracers)
  #Trans_Grid_To_Spherical!(mesh, grid_δtracers, spe_δtracers)
  Vert_Advection!(vert_coord, grid_tracers_c, grid_Δp, grid_M_half, Δt, vert_coord.vert_advect_scheme, grid_δQ)
  grid_δtracers  .+= grid_δQ
  Trans_Grid_To_Spherical!(mesh, grid_δtracers, spe_δtracers)
  #Compute_Spectral_Damping!(integrator, spe_tracers_c, spe_tracers_p, spe_δtracers)
  Filtered_Leapfrog!(integrator, spe_δtracers, spe_tracers_p, spe_tracers_c, spe_tracers_n)
  Trans_Spherical_To_Grid!(mesh, spe_tracers_n, grid_tracers_n)
  """
  # tracer_grid update
  """
  1. create a temporary field with all physical tendency
  2. perform horizontal finite volume advection from t − ∆t to t + ∆t on this field to create yet another updated field, 
     then performing vertical advection on this last field
  3. , once again from t − ∆t to t + ∆t, using one of the finite volume options 
    (we do not recommend mixing finite volume horizontal advection and centered vertical advection), 
  4. then applying the time filter
  *** Grid tracers are never transformed to the spectral domain and back. ***
  *** The advecting velocities are taken from the centered time t.        ***
  """
  grid_t_p_2  = deepcopy(grid_t_p)
  grid_ps_p_2 = deepcopy(grid_ps_p)
  grid_u_p_2 = deepcopy(grid_u_p)
  grid_v_p_2 = deepcopy(grid_v_p)
  grid_u_2 = deepcopy(grid_u)
  grid_v_2 = deepcopy(grid_v)
  grid_t_2 = deepcopy(grid_t)
  grid_u_n_2 = deepcopy(grid_u_n)
  grid_v_n_2 = deepcopy(grid_v_n)
  grid_t_n_2 = deepcopy(grid_t_n)
  spe_vor_n_2 = deepcopy(spe_vor_n)
  spe_vor_c_2 = deepcopy(spe_vor_c)
  spe_vor_p_2 = deepcopy(spe_vor_p)
  spe_div_n_2 = deepcopy(spe_div_n)
  spe_div_c_2 = deepcopy(spe_div_c)
  spe_div_p_2 = deepcopy(spe_div_p)
  grid_δu_2 = deepcopy(grid_δv)
  grid_δv_2 = deepcopy(grid_δv)
  spe_δvor_2 = deepcopy(spe_δvor)
  spe_δdiv_2 = deepcopy(spe_δdiv)
  grid_δu_2 = deepcopy(grid_δv)
  grid_δv_2 = deepcopy(grid_δv)
  grid_δt_2 = deepcopy(grid_δt)
  grid_div_2 = deepcopy(grid_div)
    
  #spe_tracers_c_2  = deepcopy(spe_tracers_c)
  spe_lnps_c_2     = deepcopy(spe_lnps_c)
  spe_lnps_p_2     = deepcopy(spe_lnps_p)
  spe_lnps_n_2     = deepcopy(spe_lnps_n)
  #grid_δtracers_2  = deepcopy(grid_δtracers)
   nλ = mesh.nλ
   nθ = mesh.nθ
   nd = mesh.nd
  grid_δtracers_2  = deepcopy(grid_δtracers)

  grid_M_half_2 = deepcopy(grid_M_half)
  grid_w_full_2 = deepcopy(grid_w_full)
  grid_δps_2    = deepcopy(grid_δps)
  grid_ps_2     = deepcopy(grid_ps)
  grid_Δp_2     = deepcopy(grid_Δp)
  grid_lnp_half_2 = deepcopy(grid_lnp_half)
  grid_lnp_full_2 = deepcopy(grid_lnp_full)
  grid_p_full_2   = deepcopy(grid_p_full)
  grid_dλ_ps_2    = deepcopy(grid_dλ_ps)
  grid_dθ_ps_2    = deepcopy(grid_dθ_ps)

  spe_t_c_2 = deepcopy(spe_t_c)
  spe_t_p_2 = deepcopy(spe_t_p)
  spe_t_n_2 = deepcopy(spe_t_n)
    
  grid_geopots_2 = deepcopy(grid_geopots)
  grid_geopot_full_2 = deepcopy(grid_geopot_full)
  grid_geopot_half_2 = deepcopy(grid_geopot_half)
    
  grid_absvor_2 = deepcopy(grid_absvor)
    
  grid_p_half_2 = deepcopy(grid_p_half)
  grid_δlnps_2  = deepcopy(grid_δlnps)
  spe_δlnps_2   = deepcopy(spe_δlnps)
  grid_vor_2    = deepcopy(grid_vor)
  grid_div_2    = deepcopy(grid_div)

  grid_energy_full_2 = deepcopy(grid_energy_full)
  spe_energy_2  = deepcopy(spe_energy)
  spe_t_p_2     = deepcopy(spe_t_p)
  spe_δt_2     = deepcopy(spe_δt)
  spe_t_n_2    = deepcopy(spe_t_n)
  grid_lnps_2  = deepcopy(grid_lnps)
  grid_ps_n_2  = deepcopy(grid_ps_n)
  grid_ps_2  = deepcopy(grid_ps)
  grid_ps_p_2  = deepcopy(grid_ps_p)
  # todo 1. and 2. 
  #############################################################################################
  ### 1. create a temporary field with all physical tendency,                               ###
  ###    and perform horizontal finite volume advection from t − ∆t to t + ∆t on this field ###
  ###    (temporary field variables named ***_2)                                            ###
  #############################################################################################
  
  integrator_2 = semi_implicit.integrator
  Δt_2 = Get_Δt(integrator_2)
  """
  mean_ps_p_2, mean_energy_p_2, sum_tracers_p_2 = Compute_Corrections_Init(vert_coord, mesh, atmo_data,
  grid_u_p_2, grid_v_p_2, grid_ps_p_2, grid_t_p_2, 
  grid_δu_2, grid_δv_2, grid_δt_2,  
  Δt_2, grid_energy_full_2, grid_tracers_p, grid_tracers_c, grid_δtracers, grid_tracers_full)
  """
  # compute pressure based on grid_ps -> grid_p_half, grid_lnp_half, grid_p_full, grid_lnp_full 
  #Pressure_Variables!(vert_coord, grid_ps_2, grid_p_half_2, grid_Δp_2, grid_lnp_half_2, grid_p_full_2, grid_lnp_full_2)
  """
  # compute ∇ps = ∇lnps * ps
  Compute_Gradients!(mesh, spe_lnps_c_2,  grid_dλ_ps_2, grid_dθ_ps_2)
  grid_dλ_ps_2 .*= grid_ps_2
  grid_dθ_ps_2 .*= grid_ps_2
    
  Four_In_One!(vert_coord, atmo_data, grid_div_2, grid_u_2, grid_v_2, grid_ps_2, 
  grid_Δp_2, grid_lnp_half_2, grid_lnp_full_2, grid_p_full_2,
  grid_dλ_ps_2, grid_dθ_ps_2, 
  grid_t_2, 
  grid_M_half_2, grid_w_full_2, grid_δu_2, grid_δv_2, grid_δps_2, grid_δt_2, grid_δtracers)

  Compute_Geopotential!(vert_coord, atmo_data, 
  grid_lnp_half_2, grid_lnp_full_2,  
  grid_t_2, 
  grid_geopots_2, grid_geopot_full_2, grid_geopot_half_2)

  grid_δlnps_2 .= grid_δps_2 ./ grid_ps_2
  Trans_Grid_To_Spherical!(mesh, grid_δlnps_2, spe_δlnps_2)  
    
  Add_Horizontal_Advection!(mesh, spe_t_c_2, grid_u_2, grid_v_2, grid_δt_2)


  Compute_Abs_Vor!(grid_vor_2, atmo_data.coriolis, grid_absvor_2)
  
  grid_δu_2 .+=  grid_absvor_2 .* grid_v_2
  grid_δv_2 .-=  grid_absvor_2 .* grid_u_2
  
  Vor_Div_From_Grid_UV!(mesh, grid_δu_2, grid_δv_2, spe_δvor_2, spe_δdiv_2)
  
  grid_energy_full_2 .= grid_geopot_full_2 .+ 0.5 * (grid_u_2.^2 + grid_v_2.^2)
  Trans_Grid_To_Spherical!(mesh, grid_energy_full_2, spe_energy_2)
  Apply_Laplacian!(mesh, spe_energy_2)
  spe_δdiv_2 .-= spe_energy_2
  
  Compute_Spectral_Damping!(integrator_2, spe_vor_c_2, spe_vor_p_2, spe_δvor_2)
  Compute_Spectral_Damping!(integrator_2, spe_div_c_2, spe_div_p_2, spe_δdiv_2)
  Compute_Spectral_Damping!(integrator_2, spe_t_c_2, spe_t_p_2, spe_δt_2)
    
  Filtered_Leapfrog!(integrator_2, spe_δvor_2, spe_vor_p_2, spe_vor_c_2, spe_vor_n_2)
  Filtered_Leapfrog!(integrator_2, spe_δdiv_2, spe_div_p_2, spe_div_c_2, spe_div_n_2)
  Filtered_Leapfrog!(integrator_2, spe_δlnps_2, spe_lnps_p_2, spe_lnps_c_2, spe_lnps_n_2)
  Filtered_Leapfrog!(integrator_2, spe_δt_2, spe_t_p_2, spe_t_c_2, spe_t_n_2)

  Trans_Spherical_To_Grid!(mesh, spe_vor_n_2, grid_vor_2)
  Trans_Spherical_To_Grid!(mesh, spe_div_n_2, grid_div_2)
  UV_Grid_From_Vor_Div!(mesh, spe_vor_n_2, spe_div_n_2, grid_u_n_2, grid_v_n_2)
  Trans_Spherical_To_Grid!(mesh, spe_t_n_2, grid_t_n_2)

  Trans_Spherical_To_Grid!(mesh, spe_lnps_n_2, grid_lnps_2)
  grid_ps_n_2 .= exp.(grid_lnps_2)
  """
  Add_Horizontal_Advection!(mesh, spe_tracers_c, grid_u_2, grid_v_2, grid_δtracers)
  #grid_tracers_n .= grid_tracers_c + Δt_2*grid_δtracers
  #grid_tracers_c .= grid_tracers_p + Δt_2*grid_δtracers
  """
  Compute_Corrections!(vert_coord, mesh, atmo_data, mean_ps_p_2, mean_energy_p_2, 
  grid_u_n_2, grid_v_n_2,
  grid_energy_full_2,
  grid_ps_n_2, spe_lnps_n_2, 
  grid_t_n_2, spe_t_n_2, grid_tracers_p, grid_tracers_c, grid_tracers_n, grid_δtracers, grid_tracers_full, sum_tracers_p, spe_tracers_n)
  """
    ####################################################################################################### 
  """
  # time advance
      # update spectral variables
  spe_vor_p_2 .= spe_vor_c_2
  spe_vor_c_2 .= spe_vor_n_2
    
  spe_div_p_2 .= spe_div_c_2
  spe_div_c_2 .= spe_div_n_2

  spe_lnps_p_2 .= spe_lnps_c_2
  spe_lnps_c_2 .= spe_lnps_n_2

  spe_t_p_2 .= spe_t_c_2
  spe_t_c_2 .= spe_t_n_2
    
  spe_tracers_p .= spe_tracers_c
  spe_tracers_c .= spe_tracers_n
  # update spectral variables
  grid_u_p_2 .= grid_u_2
  grid_u_2 .= grid_u_n_2

  grid_v_p_2 .= grid_v_2
  grid_v_2 .= grid_v_n_2

  grid_ps_p_2 .= grid_ps_2
  grid_ps_2 .= grid_ps_n_2

  grid_t_p_2 .= grid_t_2
  grid_t_2 .= grid_t_n_2

  grid_tracers_p .= grid_tracers_c
  grid_tracers_c .= grid_tracers_n
  """
  ####################################################################################################### 
  #Pressure_Variables!(vert_coord, grid_ps_2, grid_p_half_2, grid_Δp_2, grid_lnp_half_2, grid_p_full_2, grid_lnp_full_2)
  ##########################################################################################################
  ### 2. performing vertical advection on this last field                                                ###
  ###    once again from t − ∆t to t + ∆t, using one of the finite volume options                        ###
  ###    (we do not recommend mixing finite volume horizontal advection and centered vertical advection) ###                                   
  ##########################################################################################################
  """
  mean_ps_p_2, mean_energy_p_2, sum_tracers_p_2 = Compute_Corrections_Init(vert_coord, mesh, atmo_data,
  grid_u_p_2, grid_v_p_2, grid_ps_p_2, grid_t_p_2, 
  grid_δu_2, grid_δv_2, grid_δt_2,  
  Δt, grid_energy_full_2, grid_tracers_p, grid_tracers_c, grid_δtracers, grid_tracers_full)
  # compute pressure based on grid_ps -> grid_p_half, grid_lnp_half, grid_p_full, grid_lnp_full 
  Pressure_Variables!(vert_coord, grid_ps_2, grid_p_half_2, grid_Δp_2, grid_lnp_half_2, grid_p_full_2, grid_lnp_full_2)
  
  # compute ∇ps = ∇lnps * ps
  Compute_Gradients!(mesh, spe_lnps_c_2,  grid_dλ_ps_2, grid_dθ_ps_2)
  grid_dλ_ps_2 .*= grid_ps_2
  grid_dθ_ps_2 .*= grid_ps_2
  
  # compute grid_M_half, grid_w_full, grid_δu, grid_δv, grid_δps, grid_δt, 
  # except the contributions from geopotential or vertical advection
  Four_In_One!(vert_coord, atmo_data, grid_div_2, grid_u_2, grid_v_2, grid_ps_2, 
  grid_Δp_2, grid_lnp_half_2, grid_lnp_full_2, grid_p_full_2,
  grid_dλ_ps_2, grid_dθ_ps_2, 
  grid_t_2, 
  grid_M_half_2, grid_w_full_2, grid_δu_2, grid_δv_2, grid_δps_2, grid_δt_2, grid_δtracers)

  Compute_Geopotential!(vert_coord, atmo_data, 
  grid_lnp_half_2, grid_lnp_full_2,  
  grid_t_2, 
  grid_geopots_2, grid_geopot_full_2, grid_geopot_half_2)

  grid_δlnps_2 .= grid_δps_2 ./ grid_ps_2
  Trans_Grid_To_Spherical!(mesh, grid_δlnps_2, spe_δlnps_2)  
  """
  #Vert_Advection!(vert_coord, grid_u_2, grid_Δp_2, grid_M_half_2, Δt_2, vert_coord.vert_advect_scheme, grid_δQ)
  #grid_δu_2  .+= grid_δQ
  #Vert_Advection!(vert_coord, grid_v_2, grid_Δp_2, grid_M_half_2, Δt_2, vert_coord.vert_advect_scheme, grid_δQ)
  #grid_δv_2  .+= grid_δQ
  #Vert_Advection!(vert_coord, grid_t_2, grid_Δp_2, grid_M_half_2, Δt_2, vert_coord.vert_advect_scheme, grid_δQ)
  #grid_δt_2  .+= grid_δQ
  """
  Compute_Abs_Vor!(grid_vor_2, atmo_data.coriolis, grid_absvor_2)
  
  grid_δu_2 .+=  grid_absvor_2 .* grid_v_2
  grid_δv_2 .-=  grid_absvor_2 .* grid_u_2
  
  Vor_Div_From_Grid_UV!(mesh, grid_δu_2, grid_δv_2, spe_δvor_2, spe_δdiv_2)
  
  grid_energy_full_2 .= grid_geopot_full_2 .+ 0.5 * (grid_u_2.^2 + grid_v_2.^2)
  Trans_Grid_To_Spherical!(mesh, grid_energy_full_2, spe_energy_2)
  Apply_Laplacian!(mesh, spe_energy_2)
  spe_δdiv_2 .-= spe_energy_2
  
  Compute_Spectral_Damping!(integrator_2, spe_vor_c_2, spe_vor_p_2, spe_δvor_2)
  Compute_Spectral_Damping!(integrator_2, spe_div_c_2, spe_div_p_2, spe_δdiv_2)
  Compute_Spectral_Damping!(integrator_2, spe_t_c_2, spe_t_p_2, spe_δt_2)
    
  Filtered_Leapfrog!(integrator_2, spe_δvor_2, spe_vor_p_2, spe_vor_c_2, spe_vor_n_2)
  Filtered_Leapfrog!(integrator_2, spe_δdiv_2, spe_div_p_2, spe_div_c_2, spe_div_n_2)
  Filtered_Leapfrog!(integrator_2, spe_δlnps_2, spe_lnps_p_2, spe_lnps_c_2, spe_lnps_n_2)
  Filtered_Leapfrog!(integrator_2, spe_δt_2, spe_t_p_2, spe_t_c_2, spe_t_n_2)

  Trans_Spherical_To_Grid!(mesh, spe_vor_n_2, grid_vor_2)
  Trans_Spherical_To_Grid!(mesh, spe_div_n_2, grid_div_2)
  UV_Grid_From_Vor_Div!(mesh, spe_vor_n_2, spe_div_n_2, grid_u_n_2, grid_v_n_2)
  Trans_Spherical_To_Grid!(mesh, spe_t_n_2, grid_t_n_2)

  Trans_Spherical_To_Grid!(mesh, spe_lnps_n_2, grid_lnps_2)
  grid_ps_n_2 .= exp.(grid_lnps_2)
  """
  Vert_Advection!(vert_coord, grid_tracers_c, grid_Δp_2, grid_M_half_2, Δt_2, vert_coord.vert_advect_scheme, grid_δQ)
  grid_δtracers  .+= grid_δQ
  grid_tracers_n .= grid_tracers_c + Δt_2*grid_δtracers
  grid_tracers_c .= grid_tracers_p + Δt_2*grid_δtracers
  """
  Compute_Corrections!(vert_coord, mesh, atmo_data, mean_ps_p_2, mean_energy_p_2, 
  grid_u_n_2, grid_v_n_2,
  grid_energy_full_2,
  grid_ps_n_2, spe_lnps_n_2, 
  grid_t_n_2, spe_t_n_2, grid_tracers_p, grid_tracers_c, grid_tracers_n, grid_δtracers, grid_tracers_full, sum_tracers_p, spe_tracers_n)
  """
  ####################################################################################################### \
  """
  # time advance
      # update spectral variables
  spe_vor_p_2 .= spe_vor_c_2
  spe_vor_c_2 .= spe_vor_n_2
    
  spe_div_p_2 .= spe_div_c_2
  spe_div_c_2 .= spe_div_n_2

  spe_lnps_p_2 .= spe_lnps_c_2
  spe_lnps_c_2 .= spe_lnps_n_2

  spe_t_p_2 .= spe_t_c_2
  spe_t_c_2 .= spe_t_n_2
    
  spe_tracers_p .= spe_tracers_c
  spe_tracers_c .= spe_tracers_n
  # update spectral variables
  grid_u_p_2 .= grid_u_2
  grid_u_2 .= grid_u_n_2

  grid_v_p_2 .= grid_v_2
  grid_v_2 .= grid_v_n_2

  grid_ps_p_2 .= grid_ps_2
  grid_ps_2 .= grid_ps_n_2

  grid_t_p_2 .= grid_t_2
  grid_t_2 .= grid_t_n_2
  
  grid_tracers_p .= grid_tracers_c
  grid_tracers_c .= grid_tracers_n
  """
  #######################################################################################################  
  #Pressure_Variables!(vert_coord, grid_ps_2, grid_p_half_2, grid_Δp_2, grid_lnp_half_2, grid_p_full_2, grid_lnp_full_2)
  ##

  ######################################################################################################
  ###      4. then applying the time filter                                                          ###
  ### *** Grid tracers are never transformed to the spectral domain and back. ***                    ###
  ### *** The advecting velocities are taken from the centered time t.        ***                    ###
  ###     This part I'm gonna do time integration from time t to t + ∆t for passive tracers          ###
  ###     Since variables_2 has integrated to t + ∆t, I did another temporary field from t - ∆t to t ###
  ###     This temporary field variables named ***_v                                                 ###
  ######################################################################################################
  """
  grid_u_v = deepcopy(grid_u)
  grid_v_v = deepcopy(grid_v)
  grid_t_v = deepcopy(grid_t)
  grid_u_n_v = deepcopy(grid_u_n)
  grid_v_n_v = deepcopy(grid_v_n)
  grid_t_n_v = deepcopy(grid_t_n)
  spe_vor_n_v = deepcopy(spe_vor_n)
  spe_vor_c_v = deepcopy(spe_vor_c)
  spe_vor_p_v = deepcopy(spe_vor_p)
  spe_div_n_v = deepcopy(spe_div_n)
  spe_div_c_v = deepcopy(spe_div_c)
  spe_div_p_v = deepcopy(spe_div_p)
  grid_δu_v = deepcopy(grid_δv)
  grid_δv_v = deepcopy(grid_δv)
  grid_δt_v = deepcopy(grid_δt)
  spe_δvor_v = deepcopy(spe_δvor)
  spe_δdiv_v = deepcopy(spe_δdiv)
  grid_δu_v = deepcopy(grid_δv)
  grid_δv_v = deepcopy(grid_δv)
  grid_δt_v = deepcopy(grid_δt)
  grid_div_v = deepcopy(grid_div)
    
  spe_tracers_c_v  = deepcopy(spe_tracers_c)
  spe_tracers_n_v  = deepcopy(spe_tracers_n)
  grid_tracers_c_v  = deepcopy(grid_tracers_c)
  grid_tracers_n_v  = deepcopy(grid_tracers_n)
  grid_tracers_p_v  = deepcopy(grid_tracers_p)
  spe_lnps_c_v     = deepcopy(spe_lnps_c)
  spe_lnps_p_v     = deepcopy(spe_lnps_p)
  spe_lnps_n_v     = deepcopy(spe_lnps_n)
  #grid_δtracers_2  = deepcopy(grid_δtracers)
  grid_δtracers_v  = zeros(Float64, nλ,  nθ, nd)

  grid_M_half_v = deepcopy(grid_M_half)
  grid_w_full_v = deepcopy(grid_w_full)
  grid_δps_v    = deepcopy(grid_δps)
  grid_ps_v     = deepcopy(grid_ps)
  grid_Δp_v     = deepcopy(grid_Δp)
  grid_lnp_half_v = deepcopy(grid_lnp_half)
  grid_lnp_full_v = deepcopy(grid_lnp_full)
  grid_p_full_v   = deepcopy(grid_p_full)
  grid_dλ_ps_v    = deepcopy(grid_dλ_ps)
  grid_dθ_ps_v    = deepcopy(grid_dθ_ps)

  spe_t_c_v = deepcopy(spe_t_c)
    
  grid_geopots_v = deepcopy(grid_geopots)
  grid_geopot_full_v = deepcopy(grid_geopot_full)
  grid_geopot_half_v = deepcopy(grid_geopot_half)
    
  grid_absvor_v = deepcopy(grid_absvor)
    
  grid_p_half_v = deepcopy(grid_p_half)
  grid_δlnps_v  = deepcopy(grid_δlnps)
  spe_δlnps_v   = deepcopy(spe_δlnps)
  grid_vor_v    = deepcopy(grid_vor)
  grid_div_v    = deepcopy(grid_div)

  grid_energy_full_v = deepcopy(grid_energy_full)
  spe_energy_v  = deepcopy(spe_energy)
  spe_t_p_v     = deepcopy(spe_t_p)
  spe_δt_v     = deepcopy(spe_δt)
  spe_t_n_v    = deepcopy(spe_t_n)
  grid_lnps_v  = deepcopy(grid_lnps)
  grid_ps_n_v  = deepcopy(grid_ps_n)
  """
  ## vvv
  """
  Pressure_Variables!(vert_coord, grid_ps_v, grid_p_half_v, grid_Δp_v, grid_lnp_half_v, grid_p_full_v, grid_lnp_full_v)
  
  # compute ∇ps = ∇lnps * ps
  Compute_Gradients!(mesh, spe_lnps_c_v,  grid_dλ_ps_v, grid_dθ_ps_v)
  grid_dλ_ps_v .*= grid_ps_v
  grid_dθ_ps_v .*= grid_ps_v
    
  Four_In_One!(vert_coord, atmo_data, grid_div_v, grid_u_v, grid_v_v, grid_ps_v, 
  grid_Δp_v, grid_lnp_half_v, grid_lnp_full_v, grid_p_full_v,
  grid_dλ_ps_v, grid_dθ_ps_v, 
  grid_t_v, 
  grid_M_half_v, grid_w_full_v, grid_δu_v, grid_δv_v, grid_δps_v, grid_δt_v, grid_δtracers)

  Compute_Geopotential!(vert_coord, atmo_data, 
  grid_lnp_half_v, grid_lnp_full_v,  
  grid_t_v, 
  grid_geopots_v, grid_geopot_full_v, grid_geopot_half_v)

  grid_δlnps_v .= grid_δps_v ./ grid_ps_v
  Trans_Grid_To_Spherical!(mesh, grid_δlnps_v, spe_δlnps_v)  
    
  Add_Horizontal_Advection!(mesh, spe_t_c_v, grid_u_v, grid_v_v, grid_δt_v)

  Compute_Abs_Vor!(grid_vor_v, atmo_data.coriolis, grid_absvor_v)
  
  grid_δu_v .+=  grid_absvor_v .* grid_v_v
  grid_δv_v .-=  grid_absvor_v .* grid_u_v
  
  Vor_Div_From_Grid_UV!(mesh, grid_δu_v, grid_δv_v, spe_δvor_v, spe_δdiv_v)
  
  grid_energy_full_v .= grid_geopot_full_v .+ 0.5 * (grid_u_v.^2 + grid_v_v.^2)
  Trans_Grid_To_Spherical!(mesh, grid_energy_full_v, spe_energy_v)
  Apply_Laplacian!(mesh, spe_energy_v)
  spe_δdiv_v .-= spe_energy_v
  
  Compute_Spectral_Damping!(integrator, spe_vor_c_v, spe_vor_p_v, spe_δvor_v)
  Compute_Spectral_Damping!(integrator, spe_div_c_v, spe_div_p_v, spe_δdiv_v)
  Compute_Spectral_Damping!(integrator, spe_t_c_v, spe_t_p_v, spe_δt_v)
    
  Filtered_Leapfrog!(integrator, spe_δvor_v, spe_vor_p_v, spe_vor_c_v, spe_vor_n_v)
  Filtered_Leapfrog!(integrator, spe_δdiv_v, spe_div_p_v, spe_div_c_v, spe_div_n_v)
  Filtered_Leapfrog!(integrator, spe_δlnps_v, spe_lnps_p_v, spe_lnps_c_v, spe_lnps_n_v)
  Filtered_Leapfrog!(integrator, spe_δt_v, spe_t_p_v, spe_t_c_v, spe_t_n_v)

  Trans_Spherical_To_Grid!(mesh, spe_vor_n_v, grid_vor_v)
  Trans_Spherical_To_Grid!(mesh, spe_div_n_v, grid_div_v)
  UV_Grid_From_Vor_Div!(mesh, spe_vor_n_v, spe_div_n_v, grid_u_n_v, grid_v_n_v)
  Trans_Spherical_To_Grid!(mesh, spe_t_n_v, grid_t_n_v)

  Trans_Spherical_To_Grid!(mesh, spe_lnps_n_v, grid_lnps_v)
  grid_ps_n_v .= exp.(grid_lnps_v)
  """
  ### vvv
  """
  robert_coef = integrator.robert_coef
  Δt = integrator.Δt
  Add_Horizontal_Advection!(mesh, spe_tracers_c_v, grid_u_2, grid_v_2, grid_δtracers_v)
  Vert_Advection!(vert_coord, grid_tracers_c_v, grid_Δp_2, grid_M_half_2, Δt, vert_coord.vert_advect_scheme, grid_δQ)
  grid_δtracers_v  .+= grid_δQ  
  grid_tracers_n .= grid_tracers_c + Δt*grid_δtracers_v
  Pressure_Variables!(vert_coord, grid_ps_v, grid_p_half_v, grid_Δp_v, grid_lnp_half_v, grid_p_full_v, grid_lnp_full_v)
  #grid_tracers_n .= grid_tracers_n_v
  #grid_tracers_c .= grid_tracers_c_v
  #grid_tracers_p .= grid_tracers_p_v
  """
 
  Time_Advance!(dyn_data)
  ###  
  #@info "sec: ", integrator.time+1200, sum(abs.(grid_u_n)), sum(abs.(grid_v_n)), sum(abs.(grid_t_n)) , sum(abs.(grid_ps_n))
  #@info "max: ", maximum(abs.(grid_u_n)), maximum(abs.(grid_v_n)), maximum(abs.(grid_t_n)) , maximum(abs.(grid_ps_n))
  #@info "loc", grid_u_n[100,30,10],  grid_t_n[100,30,10], grid_u_n[1,32,1],  grid_t_n[1,32,1]
  
  #@assert(maximum(grid_u) <= 100.0 && maximum(grid_v) <= 100.0)

  Pressure_Variables!(vert_coord, grid_ps, grid_p_half, grid_Δp, grid_lnp_half, grid_p_full, grid_lnp_full)

    
  ###

  

  
  return 
  
  
end 

function Get_Topography!(grid_geopots::Array{Float64, 3})
  #grid_geopots .= 0.0
  
  read_file = load("0404_300_50_8_variables.dat")
  initial_day = 150
  grid_geopots .= read_file["grid_geopots_xyzt"][:,:,1,initial_day]
  
  return
end 

function Spectral_Initialize_Fields!(mesh::Spectral_Spherical_Mesh, atmo_data::Atmo_Data, vert_coord::Vert_Coordinate, sea_level_ps_ref::Float64, init_t::Float64,
  grid_geopots::Array{Float64,3}, dyn_data::Dyn_Data)
  
  spe_vor_c, spe_div_c, spe_lnps_c, spe_t_c = dyn_data.spe_vor_c, dyn_data.spe_div_c, dyn_data.spe_lnps_c, dyn_data.spe_t_c
  spe_vor_p, spe_div_p, spe_lnps_p, spe_t_p = dyn_data.spe_vor_p, dyn_data.spe_div_p, dyn_data.spe_lnps_p, dyn_data.spe_t_p
  grid_u, grid_v, grid_ps, grid_t = dyn_data.grid_u_c, dyn_data.grid_v_c, dyn_data.grid_ps_c, dyn_data.grid_t_c
  grid_u_p, grid_v_p, grid_ps_p, grid_t_p = dyn_data.grid_u_p, dyn_data.grid_v_p, dyn_data.grid_ps_p, dyn_data.grid_t_p
  
  grid_lnps,  grid_vor, grid_div =  dyn_data.grid_lnps, dyn_data.grid_vor, dyn_data.grid_div
  
  grid_p_half, grid_Δp, grid_lnp_half, grid_p_full, grid_lnp_full = dyn_data.grid_p_half, dyn_data.grid_Δp, dyn_data.grid_lnp_half, dyn_data.grid_p_full, dyn_data.grid_lnp_full
  nλ, nθ, nd = mesh.nλ, mesh.nθ, mesh.nd
  
  ### By CJY2
  spe_tracers_n = dyn_data.spe_tracers_n
  spe_tracers_c = dyn_data.spe_tracers_c
  spe_tracers_p = dyn_data.spe_tracers_p 
    
  grid_tracers_n = dyn_data.grid_tracers_n
  grid_tracers_c = dyn_data.grid_tracers_c
  grid_tracers_p = dyn_data.grid_tracers_p 
  ###

  
  rdgas = atmo_data.rdgas
  #grid_t    .= init_t
  read_file = load("0404_300_50_8_variables.dat")
  initial_day = 150
    
  grid_t    .= read_file["grid_t_c_xyzt"][:,:,:,initial_day] 
  # dΦ/dlnp = -RT    Δp = -ΔΦ/RT

  grid_geopots .= read_file["grid_geopots_xyzt"][:,:,:,initial_day]
  #grid_lnps .= log(sea_level_ps_ref) .- grid_geopots / (rdgas * init_t)
  grid_lnps .= read_file["grid_lnps_xyzt"][:,:,1,initial_day]
  #grid_ps   .= exp.(grid_lnps)
  grid_ps    .= read_file["grid_ps_xyzt"][:,:,1,initial_day]
  
  """
  # By CJY
  spe_div_c .= 0.0
  spe_vor_c .= 0.0
  """
  # # initial perturbation
  num_fourier, num_spherical = mesh.num_fourier, mesh.num_spherical
  """
  initial_perturbation = 1.0e-7/sqrt(2.0)
  # initial vorticity perturbation used in benchmark code
  # In gfdl spe[i,j] =  myspe[i, i+j-1]*√2
  
  for k = nd-2:nd
    spe_vor_c[2,5,k] = initial_perturbation
    spe_vor_c[6,9,k] = initial_perturbation
    spe_vor_c[2,4,k] = initial_perturbation  
    spe_vor_c[6,8,k] = initial_perturbation
  end
  """

  ###
  
  # By CJY2
  spe_vor_c[:,:,:] .= read_file["spe_vor_c_xyzt"][:,:,:,initial_day]
  spe_div_c[:,:,:] .= read_file["spe_div_c_xyzt"][:,:,:,initial_day]
  grid_u[:,:,:]    .= read_file["grid_u_c_xyzt"][:,:,:,initial_day]
  grid_v[:,:,:]    .= read_file["grid_v_c_xyzt"][:,:,:,initial_day]  
  ###
  
  UV_Grid_From_Vor_Div!(mesh, spe_vor_c, spe_div_c, grid_u, grid_v)
  
  # initial spectral fields (and spectrally-filtered) grid fields
  Trans_Grid_To_Spherical!(mesh, grid_t, spe_t_c)
  Trans_Spherical_To_Grid!(mesh, spe_t_c, grid_t)

  Trans_Grid_To_Spherical!(mesh, grid_lnps, spe_lnps_c)
  Trans_Spherical_To_Grid!(mesh, spe_lnps_c,  grid_lnps)

  # By CJY ### grid_ps .= exp.(grid_lnps)
  grid_ps .= read_file["grid_ps_xyzt"][:,:,1,initial_day]
  # grid_ps .= exp.(grid_lnps)
  
  Vor_Div_From_Grid_UV!(mesh, grid_u, grid_v, spe_vor_c, spe_div_c)

  UV_Grid_From_Vor_Div!(mesh, spe_vor_c, spe_div_c, grid_u, grid_v)
  
  Trans_Spherical_To_Grid!(mesh, spe_vor_c, grid_vor)
  Trans_Spherical_To_Grid!(mesh, spe_div_c, grid_div)
  
  ### By CJY2
  # try qv = 0.01
  
  Lv = 2.5*10^6.
  Rv = 461.
  # initial vorticity perturbation used in benchmark code
  # In gfdl spe[i,j] =  myspe[i, i+j-1]*√2
  one_array  = zeros(size(grid_t))
  one_array .= 1.0
    
  new_array  = zeros(size(grid_t))
  new_array .= 273.15
  initial_RH  = 0.8
  """
  initial_RH2 = 0.6
  for k = nd-2:nd
      grid_tracers_c[2:50,4:20,k] = 6.11*exp.(Lv/Rv*(one_array[2:50,4:20,k] ./ new_array[2:50,4:20,k]-one_array[2:50,4:20,k] ./ grid_t[2:50,4:20,k])) * initial_RH
      grid_tracers_c[2:50,5:20,k] = 6.11*exp.(Lv/Rv*(one_array[2:50,5:20,k] ./ new_array[2:50,5:20,k]-one_array[2:50,5:20,k] ./ grid_t[2:50,5:20,k])) * initial_RH
      grid_tracers_c[6:30,8:20,k] = 6.11*exp.(Lv/Rv*(one_array[6:30,8:20,k] ./ new_array[6:30,8:20,k]-one_array[6:30,8:20,k] ./ grid_t[6:30,8:20,k])) * initial_RH2
      grid_tracers_c[6:30,9:20,k] = 6.11*exp.(Lv/Rv*(one_array[6:30,9:20,k] ./ new_array[6:30,9:20,k]-one_array[6:30,9:20,k] ./ grid_t[6:30,9:20,k])) * initial_RH2
  end
  """
  grid_tracers_c[:,:,:] .= 6.11*exp.(Lv/Rv*(one_array[:,:,:] ./ new_array[:,:,:]-one_array[:,:,:] ./ grid_t[:,:,:])) * initial_RH
  # By CJY2
  Trans_Grid_To_Spherical!(mesh, grid_tracers_c, spe_tracers_c)
  Trans_Spherical_To_Grid!(mesh, spe_tracers_c, grid_tracers_c)
  
  ###
  #update pressure variables for hs forcing
  ### By CJY2
  grid_p_half .= read_file["grid_p_half_xyzt"][:,:,:,initial_day]
  grid_Δp .= read_file["grid_Δp_xyzt"][:,:,:,initial_day]
  grid_lnp_half .= read_file["grid_lnp_half_xyzt"][:,:,:,initial_day]
  grid_p_full .= read_file["grid_p_full_xyzt"][:,:,:,initial_day]
  grid_lnp_full .= read_file["grid_lnp_full_xyzt"][:,:,:,initial_day]
  ###
  
  Pressure_Variables!(vert_coord, grid_ps, grid_p_half, grid_Δp,
  grid_lnp_half, grid_p_full, grid_lnp_full)
  
  
  spe_vor_p .= spe_vor_c
  spe_div_p .= spe_div_c
  spe_lnps_p .= spe_lnps_c
  spe_t_p .= spe_t_c


  grid_u_p .= grid_u
  grid_v_p .= grid_v
  grid_ps_p .= grid_ps
  grid_t_p .= grid_t
  
  # By CJY2
  spe_tracers_p  .= spe_tracers_c
  grid_tracers_p .= grid_tracers_c

end 


function Spectral_Dynamics_Physics!(atmo_data::Atmo_Data, mesh::Spectral_Spherical_Mesh, dyn_data::Dyn_Data, Δt::Int64, physics_params::Dict{String, Float64})
  
  
  
  grid_δu, grid_δv, grid_δps, grid_δt = dyn_data.grid_δu, dyn_data.grid_δv, dyn_data.grid_δps, dyn_data.grid_δt
  grid_u_p, grid_v_p,  grid_t_p = dyn_data.grid_u_p, dyn_data.grid_v_p, dyn_data.grid_t_p
  grid_p_half, grid_p_full = dyn_data.grid_p_half, dyn_data.grid_p_full
  grid_t_eq = dyn_data.grid_t_eq
  
  ### By CJY2
  spe_tracers_n = dyn_data.spe_tracers_n
  spe_tracers_c = dyn_data.spe_tracers_c
  spe_tracers_p = dyn_data.spe_tracers_p 
    
  grid_tracers_n = dyn_data.grid_tracers_n
  grid_tracers_c = dyn_data.grid_tracers_c
  grid_tracers_p = dyn_data.grid_tracers_p 
    
  spe_δtracers   = dyn_data.spe_δtracers
  grid_δtracers  = dyn_data.grid_δtracers
  ###
  grid_δps .= 0.0

  HS_Forcing!(atmo_data, Δt, mesh.sinθ, grid_u_p, grid_v_p, grid_p_half, grid_p_full, grid_t_p, grid_δu, grid_δv,
  grid_t_eq, grid_δt, physics_params, grid_tracers_c)

end


function Atmosphere_Update!(mesh::Spectral_Spherical_Mesh, atmo_data::Atmo_Data, vert_coord::Vert_Coordinate, semi_implicit::Semi_Implicit_Solver, 
                            dyn_data::Dyn_Data, physcis_params::Dict{String, Float64})

  Δt = Get_Δt(semi_implicit.integrator)
  Spectral_Dynamics_Physics!(atmo_data, mesh,  dyn_data, Δt, physcis_params)
  Spectral_Dynamics!(mesh,  vert_coord , atmo_data, dyn_data, semi_implicit)

  grid_ps , grid_Δp, grid_p_half, grid_lnp_half, grid_p_full, grid_lnp_full = dyn_data.grid_ps_c,  dyn_data.grid_Δp, dyn_data.grid_p_half, 
                                                                              dyn_data.grid_lnp_half, dyn_data.grid_p_full, dyn_data.grid_lnp_full 
  grid_t = dyn_data.grid_t_c
  grid_geopots, grid_z_full, grid_z_half = dyn_data.grid_geopots, dyn_data.grid_z_full, dyn_data.grid_z_half

  ### By CJY2
  spe_tracers_n = dyn_data.spe_tracers_n
  spe_tracers_c = dyn_data.spe_tracers_c
  spe_tracers_p = dyn_data.spe_tracers_p 
    
  grid_tracers_n = dyn_data.grid_tracers_n
  grid_tracers_c = dyn_data.grid_tracers_c
  grid_tracers_p = dyn_data.grid_tracers_p 
    
  spe_δtracers   = dyn_data.spe_δtracers
  grid_δtracers  = dyn_data.grid_δtracers
  ###
    
  Compute_Pressures_And_Heights!(atmo_data, vert_coord,     
  grid_ps, grid_geopots, grid_t, 
  grid_p_half, grid_Δp, grid_lnp_half, grid_p_full, grid_lnp_full, grid_z_full, grid_z_half)

  return
end 


# function Spectral_Dynamics_Main()
#   # the decay of a sinusoidal disturbance to a zonally symmetric flow 
#   # that resembles that found in the upper troposphere in Northern winter.
#   name = "Spectral_Dynamics"
#   #num_fourier, nθ, nd = 63, 96, 20
#   num_fourier, nθ, nd = 42, 64, 20
#   #num_fourier, nθ, nd = 21, 32, 20
#   num_spherical = num_fourier + 1
#   nλ = 2nθ
  
#   radius = 6371000.0
#   omega = 7.292e-5
#   sea_level_ps_ref = 1.0e5
#   init_t = 264.0
  
#   # Initialize mesh
#   mesh = Spectral_Spherical_Mesh(num_fourier, num_spherical, nθ, nλ, nd, radius)
#   θc, λc = mesh.θc,  mesh.λc
#   cosθ, sinθ = mesh.cosθ, mesh.sinθ
  
#   vert_coord = Vert_Coordinate(nλ, nθ, nd, "even_sigma", "simmons_and_burridge", "second_centered_wts", sea_level_ps_ref)
#   # Initialize atmo_data
#   do_mass_correction = true
#   do_energy_correction = true
#   do_water_correction = false
  
#   use_virtual_temperature = false
#   atmo_data = Atmo_Data(name, nλ, nθ, nd, do_mass_correction, do_energy_correction, do_water_correction, use_virtual_temperature, sinθ, radius,  omega)
  
#   # Initialize integrator
#   damping_order = 4
#   damping_coef = 1.15741e-4
#   robert_coef  = 0.04 
  
#   implicit_coef = 0.5
#   day_to_sec = 86400
#   start_time = 0
#   end_time = 2*day_to_sec  #
#   Δt = 1200
#   init_step = true
  
#   integrator = Filtered_Leapfrog(robert_coef, 
#   damping_order, damping_coef, mesh.laplacian_eig,
#   implicit_coef, Δt, init_step, start_time, end_time)
  
#   ps_ref = sea_level_ps_ref
#   t_ref = fill(300.0, nd)
#   wave_numbers = mesh.wave_numbers
#   semi_implicit = Semi_Implicit_Solver(vert_coord, atmo_data,
#   integrator, ps_ref, t_ref, wave_numbers)
  
#   # Initialize data
#   dyn_data = Dyn_Data(name, num_fourier, num_spherical, nλ, nθ, nd)
  
  
#   NT = Int64(end_time / Δt)
  
#   Get_Topography!(dyn_data.grid_geopots)
  
#   Spectral_Initialize_Fields!(mesh, atmo_data, vert_coord, sea_level_ps_ref, init_t,
#   dyn_data.grid_geopots, dyn_data)
  

#   Atmosphere_Update!(mesh, atmo_data, vert_coord, semi_implicit, dyn_data)

#   Update_Init_Step!(semi_implicit)
#   integrator.time += Δt 
#   for i = 2:NT

#     Atmosphere_Update!(mesh, atmo_data, vert_coord, semi_implicit, dyn_data)

#     integrator.time += Δt
#     @info integrator.time

#   end
  
# end


# #Spectral_Dynamics_Main()

