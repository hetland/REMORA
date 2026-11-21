
# Every side reads the boundary file but none uses Flather, so no side gets tides and the run
# must match the same run without tides. ROMS would add tides on these clamped sides.
remora.nc_bdry_file = ../TidalBoundary/subtidal.nc
remora.nc_tide_file = ../TidalBoundary/tides.nc
remora.bc.zeta.type = clamped clamped clamped clamped
remora.bc.ubar.type = clamped clamped clamped clamped
remora.bc.vbar.type = clamped clamped clamped clamped
