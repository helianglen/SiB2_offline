!
! 'common_inc' module
!
! This module replaces the 'common.inc' file
!
! Notes:
!
! Lines beginning with '!CCC' mean (commented) original F77 source code.
! Lines begginnng with '!!!' mean F77 to F90 replacement.
!
! Last revision by nelsonvn on 2015/05/06
!

module common_inc

   implicit none

!
!  Parameters for sub-basin9 simulation
!
   integer, parameter :: max_nc = 89		! max number of ncols in dem map
   integer, parameter :: max_nr = 50		! max number of nrows in dem map
   integer, parameter :: max_stn = 3		! max number of all gauged weather stations
   integer, parameter :: max_stno = 3	! max number of all output stations
   integer, parameter :: max_soilo = 3	! max number of all output soil moisture grids
   integer, parameter :: max_time = 731	! max time steps of model input data 

   integer, parameter :: max_sub = 81		! max number of sub-catchments 
   integer, parameter :: max_flow = 28    	! max number of flow interval in a sub-catchment
   integer, parameter :: max_grid = 25    	! max number of grids in a flow-interval
   integer, parameter :: max_irrg = 5	! max number of irrigated grids for one flow interval
   integer, parameter :: max_year = 2	! max number of years in modeling

!
!  Parameters for Mogi-Guacu simulation
!
!!!   integer, parameter :: max_nc = 270		! max number of ncols in dem map
!!!   integer, parameter :: max_nr = 230		! max number of nrows in dem map
!!!   integer, parameter :: max_stn = 118		! max number of all gauged weather stations
!!!   integer, parameter :: max_stno = 15	! max number of all output stations
!!!   integer, parameter :: max_soilo = 3	! max number of all output soil moisture grids
!!!   integer, parameter :: max_time = 5114	! max time steps of model input data 
!!!
!!!   integer, parameter :: max_sub = 633		! max number of sub-catchments 
!!!   integer, parameter :: max_flow = 28    	! max number of flow interval in a sub-catchment
!!!   integer, parameter :: max_grid = 25    	! max number of grids in a flow-interval
!!!   integer, parameter :: max_irrg = 10	! max number of irrigated grids for one flow interval
!!!   integer, parameter :: max_year = 2	! max number of years in modeling



   integer, parameter :: ndvi_nc = 400		! max number of ncols in NDVI map (After patched up)
   integer, parameter :: ndvi_nr = 300		! max number of nrows in NDVI map (After patched up)
   integer, parameter :: ndvi_nt = 31		! max number of Times in one NDVI file

   integer, parameter :: fpar_nc = 70		! max number of nrows in fpar map 
   integer, parameter :: fpar_nr = 70		! max number of nrows in one fpar file

   integer, parameter :: LAI_nc = 70		! max number of nrows in LAI map 
   integer, parameter :: LAI_nr = 70		! max number of nrows in one LAI file

   integer, parameter :: nc_s = 2630			! max number of nrows in SIB2 vegetaion map 
   integer, parameter :: nr_s = 1820			! max number of nrows in SIB2 vegetaion map

!!!   real, parameter :: pi = 3.141592653589793238462643
   integer, parameter :: iout = 35
   integer, parameter :: iout1 = 38
   integer, parameter :: iout2 = 42
   integer, parameter :: iout3 = 45
   integer, parameter :: iout4 = 48
   integer, parameter :: iout5 = 69

   integer, parameter :: icho1 = 21
   integer, parameter :: icho2 = 22
   integer, parameter :: icho3 = 23
   integer, parameter :: icho4 = 24
   integer, parameter :: icho5 = 25
   integer, parameter :: icho6 = 26

!   integer :: icho1, icho2, icho3, icho4, icho5, icho6
!   integer :: icho32, icho33, icho34   ! NOT USED
!   integer :: iout, iout1, iout2, iout3, iout4, iout5
!   integer :: inpt_lai, inpt_fpar, inpt_tm, inpt_tmax, inpt_min   ! NOT USED
!   integer :: inpt_um, inpt_et, inpt_rsum, inpt_sun, inpt_fsm   ! NOT USED
!   integer :: nc_cl, nc_ia, nc_vg   ! NOT USED

end module common_inc



