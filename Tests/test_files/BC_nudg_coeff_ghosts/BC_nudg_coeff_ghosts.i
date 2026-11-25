# ------------------  INPUTS TO MAIN PROGRAM  -------------------
#
# Radiation conditions with nudging on u and v, with momentum nudging coefficients read from a
# file in 1/day. The boundary fill averages the coefficient onto a face from the cells on
# either side, which at the end of each box's part of a side is the outermost ghost layer, so
# the whole ghost region has to be converted to 1/s. On the outer ring, where the file holds
# 1/day, a cell left unconverted makes coefficient * dt about 1000 instead of 0.01: the fill
# overshoots by that factor every step, and the overflow trap stops the run. Each side is cut
# into several boxes so that such ends lie inside the domain.
#
# The grid is IdealMiniGrid's. The initial state, boundary values and climatology are uniform:
# u = 0.05, v = -0.02, temp 10, salt 35. The nudging coefficient is 1/day on the outer ring
# of the grid file, then 0.5, 0.1 and 0 inward.

# PROBLEM SIZE & GEOMETRY
remora.prob_name = IdealMiniGrid
remora.max_step = 100
remora.prob_lo     =   0.       0.  -5000.
remora.prob_hi     =   320000.  320000.       0.
remora.n_cell      =  10 16 4
remora.is_periodic = 0 0 0
amr.max_grid_size  = 4

amrex.fpe_trap_overflow = 1

# BOUNDARY CONDITIONS
remora.boundary_per_variable = 1
#                         West               South              East               North
remora.bc.temp.type   =   slipwall           slipwall           slipwall           slipwall
remora.bc.salt.type   =   slipwall           slipwall           slipwall           slipwall
remora.bc.scalar.type =   outflow            outflow            outflow            outflow
remora.bc.u.type      =   orlanski_rad_nudg  orlanski_rad_nudg  orlanski_rad_nudg  orlanski_rad_nudg
remora.bc.v.type      =   orlanski_rad_nudg  orlanski_rad_nudg  orlanski_rad_nudg  orlanski_rad_nudg
remora.bc.w.type      =   outflow            outflow            outflow            outflow
remora.bc.ubar.type   =   slipwall           slipwall           slipwall           slipwall
remora.bc.vbar.type   =   slipwall           slipwall           slipwall           slipwall
remora.bc.zeta.type   =   slipwall           slipwall           slipwall           slipwall
remora.bc.tke.type    =   outflow            outflow            outflow            outflow
remora.obcfac = 1.0

# TIME STEP CONTROL
remora.fixed_dt = 2000.0
remora.ndtfast  = 20

# CHECKPOINT AND PLOT FILES
remora.check_int     = -1
remora.plot_file     = plt
remora.plot_int      = 100
remora.plot_vars_3d  = salt temp x_velocity y_velocity z_velocity
remora.plotfile_type = amrex

# SOLVER CHOICE
remora.tracer_horizontal_advection_scheme = "upstream3"
remora.use_coriolis  = true
remora.coriolis_type = netcdf

remora.theta_s = 5.0
remora.theta_b = 0.5
remora.tcline  = 100.

# Linear EOS parameters
remora.R0    = 1027.0
remora.S0    = 32.0
remora.T0    = 10.0
remora.Tcoef = 0.0
remora.Scoef = 0.0
remora.rho0  = 1025.0

# INITIAL, BOUNDARY AND CLIMATOLOGY DATA
remora.ic_bc_type         = "netcdf"
remora.time_ref           = 20200101
remora.nc_init_file_0     = "nudg_ghosts_ini.nc"
remora.nc_grid_file_0     = "nudg_ghosts_grd.nc"
remora.nc_bdry_file       = "nudg_ghosts_bry.nc"
remora.do_m3_clim_nudg    = true
remora.nc_clim_his_file   = "nudg_ghosts_clm.nc"
remora.nc_clim_coeff_file = "nudg_ghosts_nudge.nc"
