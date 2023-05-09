export Output_Manager, Update_Output!, Finalize_Output!
export Lat_Lon_Pcolormesh, Zonal_Mean, Sigma_Zonal_Mean_Pcolormesh, Sigma_Zonal_Mean_Contourf
                             
mutable struct Output_Manager
    nλ::Int64
    nθ::Int64
    nd::Int64
    n_day::Int64
    
    day_to_sec::Int64
    start_time::Int64
    end_time::Int64
    current_time::Int64
    spinup_day::Int64

    λc::Array{Float64, 1}
    θc::Array{Float64, 1}
    σc::Array{Float64, 1}
  
    # nθ × nd × n_day
    # The average is (start, end], namely it does not include the first snapshot.
    t_daily_zonal_mean::Array{Float64, 3}
    t_eq_daily_zonal_mean::Array{Float64, 3}
    u_daily_zonal_mean::Array{Float64, 3}
    v_daily_zonal_mean::Array{Float64, 3}
    
    ps_daily_mean::Array{Float64, 3}

    n_daily_mean::Array{Float64, 1}
    

    # The average from spinup_day+1 to n_day
    t_zonal_mean::Array{Float64, 2}
    t_eq_zonal_mean::Array{Float64, 2}
    u_zonal_mean::Array{Float64, 2}
    v_zonal_mean::Array{Float64, 2}
    ps_mean::Array{Float64, 2}
    
    ### By CJY
    grid_u_c_xyzt::Array{Float64, 4}
    grid_v_c_xyzt::Array{Float64, 4}
    grid_t_c_xyzt::Array{Float64, 4}
    grid_t_eq_xyzt::Array{Float64, 4}
    
    grid_geopots_xyzt::Array{Float64,4}
    grid_ps_xyzt::Array{Float64,4}
    spe_vor_c_xyzt::Array{ComplexF64, 4}
    spe_div_c_xyzt::Array{ComplexF64, 4}
    grid_lnps_xyzt::Array{Float64,4}
    
    grid_p_half_xyzt::Array{Float64,4}
    grid_Δp_xyzt::Array{Float64,4} 
    grid_lnp_half_xyzt::Array{Float64,4} 
    grid_p_full_xyzt::Array{Float64,4} 
    grid_lnp_full_xyzt::Array{Float64,4}
    
    spe_lnps_c_xyzt::Array{ComplexF64,4}
    spe_lnps_p_xyzt::Array{ComplexF64,4}
    
    ### By CJY2
    grid_tracers_n_xyz1t::Array{Float64,4}
    grid_tracers_c_xyz1t::Array{Float64,4}
    grid_tracers_p_xyz1t::Array{Float64,4}

    spe_tracers_n_xyz1t::Array{ComplexF64,4}
    spe_tracers_c_xyz1t::Array{ComplexF64,4}
    spe_tracers_p_xyz1t::Array{ComplexF64,4}
    ###
    ### By CJY3
    grid_tracers_c_all::Array{Float64,4}
    tracers_correction_all::Array{Float64,1}
end

function Output_Manager(mesh::Spectral_Spherical_Mesh, vert_coord::Vert_Coordinate, start_time::Int64, end_time::Int64, spinup_day::Int64, num_grid_tracters::Int64=1, num_spe_tracters::Int64=1) ### By CJY2
    nλ = mesh.nλ
    nθ = mesh.nθ
    nd = mesh.nd
    
    day_to_sec = 86400
    current_time = start_time

    λc = mesh.λc
    θc = mesh.θc

    #todo definition of sigma coordinate
    bk = vert_coord.bk
    σc = (bk[2:nd+1] + bk[1:nd])/2.0
  
    # nθ × nd × n_day
    # The average is (start, end], namely it does not include the first snapshot.
    n_day = Int64((end_time - start_time)/day_to_sec)
    t_daily_zonal_mean = zeros(Float64, nθ, nd, n_day)
    t_eq_daily_zonal_mean = zeros(Float64, nθ, nd, n_day)
    u_daily_zonal_mean = zeros(Float64, nθ, nd, n_day)
    v_daily_zonal_mean = zeros(Float64, nθ, nd, n_day)
    ps_daily_mean = zeros(Float64, nλ, nθ, n_day)
    n_daily_mean = zeros(Float64, n_day)

    # The average from spinup_day+1 to n_day
    t_zonal_mean = zeros(Float64, nθ, nd)
    t_eq_zonal_mean = zeros(Float64, nθ, nd)
    u_zonal_mean = zeros(Float64, nθ, nd)
    v_zonal_mean = zeros(Float64, nθ, nd)
    ps_mean = zeros(Float64, nλ, nθ)

    ### By CJY
    grid_u_c_xyzt  = zeros(Float64, nλ, nθ, nd, n_day) ####
    grid_v_c_xyzt  = zeros(Float64, nλ, nθ, nd, n_day) 
    grid_t_c_xyzt  = zeros(Float64, nλ,  nθ, nd, n_day)
    grid_t_eq_xyzt = zeros(Float64, nλ,  nθ, nd, n_day)
    
    grid_geopots_xyzt = zeros(Float64, nλ, nθ, 1, n_day)
    grid_ps_xyzt = zeros(Float64, nλ, nθ, 1, n_day)
    num_fourier, nθ, nd = 42, 64, 20
    num_spherical = num_fourier + 1
    spe_vor_c_xyzt = zeros(ComplexF64, num_fourier+1, num_spherical+1, nd, n_day)
    spe_div_c_xyzt = zeros(ComplexF64, num_fourier+1, num_spherical+1, nd, n_day)
    grid_lnps_xyzt = zeros(Float64, nλ,  nθ, 1, n_day)
    
    grid_p_half_xyzt   = zeros(Float64, nλ,  nθ, nd+1, n_day)
    grid_Δp_xyzt       = zeros(Float64, nλ,  nθ, nd, n_day)
    grid_lnp_half_xyzt = zeros(Float64, nλ,  nθ, nd+1, n_day)
    grid_p_full_xyzt   = zeros(Float64, nλ,  nθ, nd, n_day)
    grid_lnp_full_xyzt = zeros(Float64, nλ,  nθ, nd, n_day)
    
    spe_lnps_c_xyzt    = zeros(ComplexF64, num_fourier+1, num_spherical+1, 1, n_day)
    spe_lnps_p_xyzt    = zeros(ComplexF64, num_fourier+1, num_spherical+1, 1, n_day)
    
    ### By CJY2
    grid_tracers_n_xyz1t = zeros(Float64, nλ,  nθ, nd, n_day)
    grid_tracers_c_xyz1t = zeros(Float64, nλ,  nθ, nd, n_day)
    grid_tracers_p_xyz1t = zeros(Float64, nλ,  nθ, nd, n_day)
    
    spe_tracers_n_xyz1t = zeros(ComplexF64, num_fourier+1, num_spherical+1, nd, n_day)
    spe_tracers_c_xyz1t = zeros(ComplexF64, num_fourier+1, num_spherical+1, nd, n_day)
    spe_tracers_p_xyz1t = zeros(ComplexF64, num_fourier+1, num_spherical+1, nd, n_day)
    ###
    
    ### By CJY3
    spe_tracers_c_all = zeros(Float64, nλ,  nθ, nd, 1)
    tracers_correction_all = zeros(Float64, n_day)
    ###
    Output_Manager(nλ, nθ, nd, n_day,
    day_to_sec, start_time, end_time, current_time, spinup_day,
    λc, θc, σc,
    t_daily_zonal_mean, t_eq_daily_zonal_mean, u_daily_zonal_mean, v_daily_zonal_mean, 
    ps_daily_mean, n_daily_mean, 
    t_zonal_mean,t_eq_zonal_mean, u_zonal_mean, v_zonal_mean, ps_mean, grid_u_c_xyzt, grid_v_c_xyzt, grid_t_c_xyzt, grid_t_eq_xyzt, grid_geopots_xyzt, grid_ps_xyzt, spe_vor_c_xyzt, spe_div_c_xyzt, grid_lnps_xyzt, grid_p_half_xyzt, grid_Δp_xyzt,  grid_lnp_half_xyzt,  grid_p_full_xyzt, grid_lnp_full_xyzt, spe_lnps_c_xyzt, spe_lnps_p_xyzt, grid_tracers_n_xyz1t, grid_tracers_c_xyz1t, grid_tracers_p_xyz1t, spe_tracers_n_xyz1t, spe_tracers_c_xyz1t, spe_tracers_p_xyz1t, spe_tracers_c_all, tracers_correction_all)
end


function Update_Output!(output_manager::Output_Manager, dyn_data::Dyn_Data, current_time::Int64)
    @assert(current_time > output_manager.current_time)
    output_manager.current_time = current_time
    day_to_sec, start_time, n_day = output_manager.day_to_sec, output_manager.start_time, output_manager.n_day

    t_daily_zonal_mean, t_eq_daily_zonal_mean, u_daily_zonal_mean, v_daily_zonal_mean, ps_daily_mean, n_daily_mean = 
    output_manager.t_daily_zonal_mean, output_manager.t_eq_daily_zonal_mean,
    output_manager.u_daily_zonal_mean, output_manager.v_daily_zonal_mean, 
    output_manager.ps_daily_mean, output_manager.n_daily_mean

    ### By CJY
    grid_u_c_xyzt  = output_manager.grid_u_c_xyzt ###
    grid_v_c_xyzt  = output_manager.grid_v_c_xyzt
    grid_t_c_xyzt  = output_manager.grid_t_c_xyzt
    grid_t_eq_xyzt = output_manager.grid_t_eq_xyzt

    grid_geopots_xyzt = output_manager.grid_geopots_xyzt
    grid_ps_xyzt   = output_manager.grid_ps_xyzt
    spe_vor_c_xyzt = output_manager.spe_vor_c_xyzt
    spe_div_c_xyzt = output_manager.spe_div_c_xyzt
    grid_lnps_xyzt = output_manager.grid_lnps_xyzt
    
    grid_p_half_xyzt   = output_manager.grid_p_half_xyzt 
    grid_Δp_xyzt       = output_manager.grid_Δp_xyzt
    grid_lnp_half_xyzt = output_manager.grid_lnp_half_xyzt 
    grid_p_full_xyzt   = output_manager.grid_p_full_xyzt 
    grid_lnp_full_xyzt = output_manager.grid_lnp_full_xyzt
    
    spe_lnps_c_xyzt    = output_manager.spe_lnps_c_xyzt
    spe_lnps_p_xyzt    = output_manager.spe_lnps_p_xyzt
    ###
    ### By CJY2
    grid_tracers_n_xyz1t = output_manager.grid_tracers_n_xyz1t
    grid_tracers_c_xyz1t = output_manager.grid_tracers_c_xyz1t
    grid_tracers_p_xyz1t = output_manager.grid_tracers_p_xyz1t
    
    spe_tracers_n_xyz1t = output_manager.spe_tracers_n_xyz1t
    spe_tracers_c_xyz1t = output_manager.spe_tracers_c_xyz1t
    spe_tracers_p_xyz1t = output_manager.spe_tracers_p_xyz1t
    ###
    ### By CJY3
    tracers_correction_all = output_manager.tracers_correction_all
    ###
    
    i_day = div(current_time - start_time - 1, day_to_sec) + 1

    if(i_day > n_day)
        @info "Warning: i_day > n_day in Output_Manager:Update!"
        return 
    end
    
    t_daily_zonal_mean[:,:,i_day] .+= Zonal_Mean(dyn_data.grid_t_c)
    t_eq_daily_zonal_mean[:,:,i_day] .+= Zonal_Mean(dyn_data.grid_t_eq)
    u_daily_zonal_mean[:,:,i_day] .+= Zonal_Mean(dyn_data.grid_u_c)
    v_daily_zonal_mean[:,:,i_day] .+= Zonal_Mean(dyn_data.grid_v_c)

    ps_daily_mean[:,:,i_day] .+= dyn_data.grid_ps_c[:,:,1]
    ### By CJY
    grid_u_c_xyzt[:,:,:,i_day]  .= dyn_data.grid_u_c[:,:,:]
    grid_v_c_xyzt[:,:,:,i_day]  .= dyn_data.grid_v_c[:,:,:]
    grid_t_c_xyzt[:,:,:,i_day]  .= dyn_data.grid_t_c[:,:,:]
    grid_t_eq_xyzt[:,:,:,i_day] .= dyn_data.grid_t_eq[:,:,:]
    
    grid_geopots_xyzt[:,:,1,i_day] .= dyn_data.grid_geopots[:,:,1]
    grid_ps_xyzt[:,:,1,i_day]      .= dyn_data.grid_ps_c[:,:,1]
    spe_vor_c_xyzt[:,:,:,i_day]    .= dyn_data.spe_vor_c[:,:,:]
    spe_div_c_xyzt[:,:,:,i_day]    .= dyn_data.spe_div_c[:,:,:]
    grid_lnps_xyzt[:,:,1,i_day]    .= dyn_data.grid_lnps[:,:,1]
    
    grid_p_half_xyzt[:,:,:,i_day]    .= dyn_data.grid_p_half[:,:,:]
    grid_Δp_xyzt[:,:,:,i_day]        .= dyn_data.grid_Δp[:,:,:]
    grid_lnp_half_xyzt[:,:,:,i_day]  .= dyn_data.grid_lnp_half[:,:,:]
    grid_p_full_xyzt[:,:,:,i_day]    .= dyn_data.grid_p_full[:,:,:]
    grid_lnp_full_xyzt[:,:,:,i_day]  .= dyn_data.grid_lnp_full[:,:,:]
    
    spe_lnps_c_xyzt[:,:,1,i_day]     .= dyn_data.spe_lnps_c[:,:,1]
    spe_lnps_p_xyzt[:,:,1,i_day]     .= dyn_data.spe_lnps_p[:,:,1]
    ###
    ### By CJY2
    grid_tracers_n_xyz1t[:,:,:,i_day] .= dyn_data.grid_tracers_n[:,:,:]
    grid_tracers_c_xyz1t[:,:,:,i_day] .= dyn_data.grid_tracers_c[:,:,:]
    grid_tracers_p_xyz1t[:,:,:,i_day] .= dyn_data.grid_tracers_p[:,:,:]

    spe_tracers_n_xyz1t[:,:,:,i_day] .= dyn_data.spe_tracers_n[:,:,:]
    spe_tracers_c_xyz1t[:,:,:,i_day] .= dyn_data.spe_tracers_c[:,:,:]
    spe_tracers_p_xyz1t[:,:,:,i_day] .= dyn_data.spe_tracers_p[:,:,:]
    #print(maximum(grid_tracers_c_xyz1t))
    ### By CJY3
    tracers_correction_all[i_day] = sum(dyn_data.grid_tracers_c[:,:,:])
    ###
    n_daily_mean[i_day] += 1
end

function Finalize_Output!(output_manager::Output_Manager, save_file_name::String = "None", mean_save_file_name::String = "None")

    n_day = output_manager.n_day

    t_daily_zonal_mean, t_eq_daily_zonal_mean, u_daily_zonal_mean, v_daily_zonal_mean, ps_daily_mean, n_daily_mean = 
    output_manager.t_daily_zonal_mean, output_manager.t_eq_daily_zonal_mean,
    output_manager.u_daily_zonal_mean, output_manager.v_daily_zonal_mean, 
    output_manager.ps_daily_mean, output_manager.n_daily_mean
    ###
    grid_u_c_xyzt  = output_manager.grid_u_c_xyzt ###
    grid_v_c_xyzt  = output_manager.grid_v_c_xyzt
    grid_t_c_xyzt  = output_manager.grid_t_c_xyzt
    grid_t_eq_xyzt = output_manager.grid_t_eq_xyzt
    
    grid_geopots_xyzt = output_manager.grid_geopots_xyzt
    grid_ps_xyzt      = output_manager.grid_ps_xyzt
    spe_vor_c_xyzt    = output_manager.spe_vor_c_xyzt
    spe_div_c_xyzt    = output_manager.spe_div_c_xyzt
    grid_lnps_xyzt    = output_manager.grid_lnps_xyzt
    
    grid_p_half_xyzt   = output_manager.grid_p_half_xyzt 
    grid_Δp_xyzt       = output_manager.grid_Δp_xyzt
    grid_lnp_half_xyzt = output_manager.grid_lnp_half_xyzt 
    grid_p_full_xyzt   = output_manager.grid_p_full_xyzt 
    grid_lnp_full_xyzt = output_manager.grid_lnp_full_xyzt
    
    spe_lnps_c_xyzt    = output_manager.spe_lnps_c_xyzt
    spe_lnps_p_xyzt    = output_manager.spe_lnps_p_xyzt
    ###
    
    ### By CJY2
    grid_tracers_n_xyz1t = output_manager.grid_tracers_n_xyz1t
    grid_tracers_c_xyz1t = output_manager.grid_tracers_c_xyz1t
    grid_tracers_p_xyz1t = output_manager.grid_tracers_p_xyz1t
    
    spe_tracers_n_xyz1t = output_manager.spe_tracers_n_xyz1t
    spe_tracers_c_xyz1t = output_manager.spe_tracers_c_xyz1t
    spe_tracers_p_xyz1t = output_manager.spe_tracers_p_xyz1t
    ###
    grid_tracers_c_all = output_manager.grid_tracers_c_all
    ### By CJY3
    tracers_correction_all = output_manager.tracers_correction_all
    ###
    
    for i_day = 1:n_day
        grid_tracers_c_all[:,:,:,1] .+= grid_tracers_c_xyz1t[:,:,:,i_day]
    end
    
    """
    for i_day = 1:n_day
        t_daily_zonal_mean[:,:,i_day] ./= n_daily_mean[i_day]
        t_eq_daily_zonal_mean[:,:,i_day] ./= n_daily_mean[i_day]
        u_daily_zonal_mean[:,:,i_day] ./= n_daily_mean[i_day]
        v_daily_zonal_mean[:,:,i_day] ./= n_daily_mean[i_day]
        ps_daily_mean[:,:,i_day] ./= n_daily_mean[i_day]
        n_daily_mean[i_day] = 1.0
    end

    spinup_day = output_manager.spinup_day
    t_zonal_mean, t_eq_zonal_mean, u_zonal_mean, v_zonal_mean, ps_mean = 
    output_manager.t_zonal_mean, output_manager.t_eq_zonal_mean, 
    output_manager.u_zonal_mean, output_manager.v_zonal_mean, output_manager.ps_mean

    t_zonal_mean .= dropdims(mean(t_daily_zonal_mean[:,:,spinup_day+1:n_day], dims=3), dims=3)
    t_eq_zonal_mean .= dropdims(mean(t_eq_daily_zonal_mean[:,:,spinup_day+1:n_day], dims=3), dims=3)
    u_zonal_mean .= dropdims(mean(u_daily_zonal_mean[:,:,spinup_day+1:n_day], dims=3), dims=3)
    v_zonal_mean .= dropdims(mean(v_daily_zonal_mean[:,:,spinup_day+1:n_day], dims=3), dims=3)
    ps_mean .= dropdims(mean(ps_daily_mean[:,:,spinup_day+1:n_day], dims=3), dims=3)
    """ 
    
    if save_file_name != "None"
        @save save_file_name output_manager
    end

    if mean_save_file_name != "None"
        @save mean_save_file_name grid_u_c_xyzt grid_v_c_xyzt grid_t_c_xyzt grid_t_eq_xyzt grid_geopots_xyzt grid_ps_xyzt spe_vor_c_xyzt spe_div_c_xyzt  grid_p_full_xyzt  spe_lnps_c_xyzt  grid_tracers_c_xyz1t grid_tracers_c_all tracers_correction_all
    end
end





function Sigma_Zonal_Mean_Contourf(output_manager::Output_Manager, save_file_pref::String)
    
    θc = output_manager.θc
    σc = output_manager.σc
    nθ = length(θc)
    θc_deg = θc*180/pi
    nd = output_manager.nd
    
    X,Y = repeat(θc_deg, 1, nd), repeat(σc, 1, nθ)'

    t_zonal_mean, t_eq_zonal_mean, u_zonal_mean, v_zonal_mean = output_manager.t_zonal_mean, output_manager.t_eq_zonal_mean, output_manager.u_zonal_mean, output_manager.v_zonal_mean
   
    PyPlot.contourf(X, Y, t_zonal_mean, levels = 10)
    PyPlot.gca().invert_yaxis()
    PyPlot.xlabel("Latitude")
    PyPlot.ylabel("σ")
    PyPlot.colorbar()
    PyPlot.savefig(save_file_pref * "_T.png")
    PyPlot.close("all")

    PyPlot.contourf(X, Y, t_eq_zonal_mean, levels = 10)
    PyPlot.gca().invert_yaxis()
    PyPlot.xlabel("Latitude")
    PyPlot.ylabel("σ")
    PyPlot.colorbar()
    PyPlot.savefig(save_file_pref * "_Teq.png")
    PyPlot.close("all")

    PyPlot.contourf(X, Y, u_zonal_mean, levels = 10)
    PyPlot.gca().invert_yaxis()
    PyPlot.xlabel("Latitude")
    PyPlot.ylabel("σ")
    PyPlot.colorbar()
    PyPlot.savefig(save_file_pref * "_U.png")
    PyPlot.close("all")

    PyPlot.contourf(X, Y, v_zonal_mean, levels = 10)
    PyPlot.gca().invert_yaxis()
    PyPlot.xlabel("Latitude")
    PyPlot.ylabel("σ")
    PyPlot.colorbar()
    PyPlot.savefig(save_file_pref * "_V.png")
    PyPlot.close("all")
    
    
end


function Sigma_Zonal_Mean_Pcolormesh(output_manager::Output_Manager, save_file_pref::String)
    
    θc = output_manager.θc
    σc = output_manager.σc
    nθ = length(θc)
    θc_deg = θc*180/pi
    nd = output_manager.nd
    
    X,Y = repeat(θc_deg, 1, nd), repeat(σc, 1, nθ)'

    t_zonal_mean, u_zonal_mean, v_zonal_mean = output_manager.t_zonal_mean, output_manager.u_zonal_mean, output_manager.v_zonal_mean
   
    PyPlot.pcolormesh(X, Y, t_zonal_mean, shading= "gouraud", cmap="viridis")
    PyPlot.gca().invert_yaxis()
    PyPlot.colorbar()
    PyPlot.savefig(save_file_pref * "_T.png")
    PyPlot.close("all")

    PyPlot.pcolormesh(X, Y, u_zonal_mean, shading= "gouraud", cmap="viridis")
    PyPlot.gca().invert_yaxis()
    PyPlot.colorbar()
    PyPlot.savefig(save_file_pref * "_U.png")
    PyPlot.close("all")

    PyPlot.pcolormesh(X, Y, v_zonal_mean, shading= "gouraud", cmap="viridis")
    PyPlot.gca().invert_yaxis()
    PyPlot.colorbar()
    PyPlot.savefig(save_file_pref * "_V.png")
    PyPlot.close("all")
    
    
end


function Lat_Lon_Pcolormesh(output_manager::Output_Manager, grid_dat::Array{Float64,3}, level::Int64, save_file_name::String = "None")
    
    λc, θc = output_manager.λc, output_manager.θc
    nλ, nθ = length(λc), length(θc)
    λc_deg, θc_deg = λc*180/pi, θc*180/pi
    
    
    X,Y = repeat(λc_deg, 1, nθ), repeat(θc_deg, 1, nλ)'
    
    
    PyPlot.pcolormesh(X, Y, grid_dat[:,:,level], shading= "gouraud", cmap="viridis")
    PyPlot.axis("equal")
    PyPlot.colorbar()
    
    if save_file_name != "None"
        PyPlot.savefig(save_file_name)
        PyPlot.close("all")
    end
    
end


function Lat_Lon_Pcolormesh(mesh::Spectral_Spherical_Mesh, grid_dat::Array{Float64,3}, level::Int64, save_file_name::String = "None")
    
    λc, θc = mesh.λc, mesh.θc
    nλ, nθ = length(λc), length(θc)
    λc_deg, θc_deg = λc*180/pi, θc*180/pi
    
    X,Y = repeat(λc_deg, 1, nθ), repeat(θc_deg, 1, nλ)'
    
    
    PyPlot.pcolormesh(X, Y, grid_dat[:,:,level], shading= "gouraud", cmap="viridis")
    PyPlot.axis("equal")
    PyPlot.colorbar()
    
    if save_file_name != "None"
        PyPlot.savefig(save_file_name)
        PyPlot.close("all")
    end
    
end


function Zonal_Mean(grid_dat::Array{Float64,3})
    
    return dropdims(mean(grid_dat, dims=1), dims=1)
    
end


function Sigma_Zonal_Mean_Pcolormesh(output_manager::Output_Manager,
    zonal_mean_data::Array{Float64,2}, save_file_name::String = "None")
    
    θc = output_manager.θc
    σc = output_manager.σc
    nθ = length(θc)
    θc_deg = θc*180/pi
    nd = output_manager.nd
    
    X,Y = repeat(θc_deg, 1, nd), repeat(σc, 1, nθ)'
    
    
    PyPlot.pcolormesh(X, Y, zonal_mean_data, shading= "gouraud", cmap="viridis")
    PyPlot.colorbar()
    PyPlot.gca().invert_yaxis()
    
    if save_file_name != "None"
        PyPlot.savefig(save_file_name)
        PyPlot.close("all")
    end
    
end


function Sigma_Zonal_Mean_Contourf(output_manager::Output_Manager, 
    zonal_mean_data::Array{Float64,2}, save_file_name::String = "None")
    
    θc = output_manager.θc
    σc = output_manager.σc
    nθ = length(θc)
    θc_deg = θc*180/pi
    nd = output_manager.nd
    
    X,Y = repeat(θc_deg, 1, nd), repeat(σc, 1, nθ)'
    
    PyPlot.contourf(X, Y, zonal_mean_data, levels = 10)
    PyPlot.gca().invert_yaxis()
    
    if save_file_name != "None"
        PyPlot.savefig(save_file_name)
        PyPlot.close("all")
    end
    
end