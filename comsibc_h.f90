!
! This module replaces 'COMSIBC.H' file - It was tested with IUNDEF and RUNDEF
!

module stepv
   implicit none
   real :: tc
   real :: tg
   real :: td
   real :: capac(2)
   real :: snoww(2)
   real :: www(3)
end module stepv

module multicamada ! Adicionado RHV 9/06/2016
  implicit none
  integer :: nlayers
  real :: dzmultilayer
  real :: depth(100)
  real :: dzwww(100)
  real :: extfrac(100)
  real :: extfrac_new(100)
  real :: satcoz_d(100)
  real :: temw_d(100)
  real :: temwp_d(100)
  real :: temwpp_d(100)
  integer :: hr_proc
end module multicamada

module const
   implicit none
   real, parameter :: pie = 3.14159265358979323846
   real, parameter :: g = 9.81
   real, parameter :: asnow = 13.2
   real, parameter :: clai = 4.2 * 1000.0 * 0.2
   real, parameter :: cpair = 1010.0
   real, parameter :: cw = 4.2 * 1000.0 * 1000.0
   real, parameter :: epsfac = 0.622
   real, parameter :: kappa = 0.286
   real, parameter :: rhoair = 1.225
   real, parameter :: snomel = 370518.5 * 1000.0
   real, parameter :: stefan = 5.669 * 10e-9
   real, parameter :: tf = 273.16
   real, parameter :: vkc = 0.41
   real, parameter :: rcp = rhoair * cpair
   real, parameter :: timcon = pie / 86400.0
   real :: psy
   real :: hlat
   real :: snofac
end module 



module  atchem
   implicit none
   real, parameter :: po2m = 20900.0
   real, parameter :: pco2m = 34.0
end module 



module gridij
   implicit none
   integer :: ivtype 
   integer :: istype 
end module 



module vstate
   implicit none
   real :: z2
   real :: z1
   real :: vcover
   real :: stem
   real :: chil
   real :: tran(2,2)
   real :: ref(2,2)
!!!cxx      rootd, ph1, ph2, &
   real :: rootd
   real :: phc
   real :: effcon
   real :: gradm
   real :: binter
   real :: respcp
   real :: atheta
   real :: btheta
!!!cxx      trda, trdm, trop, tpsa, tpsb, tpsc
   real :: trda 
   real :: trdm
   real :: trop
   real :: slti
   real :: hlti
   real :: shti
   real :: hhti
   real :: rootex
end module 



module vdyijt
   implicit none
   real :: zlt
   real :: green
   real :: fparc
end module 



module vderiv
   implicit none
   real :: z0d
   real :: dd
   real :: cc1
   real :: cc2
   real :: vmax0
   real :: gmudmu
end module 



module soilij
   implicit none
   real :: sodep
   real :: soref(2)
end module 



module soils
   implicit none
   real :: bee
   real :: phsat
   real :: poros
   real :: satco
   real :: slope
   real :: zdepth(3)   ! tatsch, 6 jun 2011
   real :: slpp        ! tatsch, 6 jun 2011
   real :: satcoz(3)   ! tatsch 15 set 2011
   real :: satcoz_a    ! tatsch 15 set 2011
   real :: zlay(3)     ! tatsch 15 set 2011
   real :: zlay_a      ! tatsch 15 set 2011
   real :: decay       ! tatsch 15 set 2011
   real :: anik        ! roilan 8/06/2016 
!!!c
!!!c     satcoz = satco * exp(-decay*zlay)
!!!c     satcoz_a = sum(satcoz(i, i =1,3))/3
!!!c     zlay ! Mid-level of soil layer (m)
!!!c     zlay_a! Mid-level of soil layer (m)
end module 



module atmos
   implicit none
   real, parameter :: bps = 1.0
   real, parameter :: psur = 1000.0
   real :: em
   real :: tm
   real :: um
   real :: zm
   real :: ppc
   real :: ppl
   real :: radn(3,2)
   real :: sunang
   real :: coszen
   real :: swdown
   real :: rnetm
   real :: cloud
end module 



module caerod
   implicit none
   real :: corb1
   real :: corb2
   real :: ha
   real :: g1
   real :: g2
   real :: g3
   real :: ztz0
   real :: zwind
   real :: zmet
   real :: ht
end module 



module site
   implicit none
   real :: zlong
   real :: zlat
   real :: salb(2,2)
   real :: rab(2,3,2)
end module 



module steps
   implicit none
   integer :: itrunk 
   integer :: ilw 
   integer :: niter 
   integer :: iter 
!   integer :: ispare   ! NOT USED
   real :: dtt
end module 



!!! Note that, by default, TIME, YEAR, DAY and HOUR are REAL and MONTH is INTEGER
!!! Variables MONTH and HOUR are NEVER USED
module govern
   implicit none
!!!   integer :: time, year, month, day, hour
   integer :: year
   integer :: month 
   integer :: day
   integer :: hour
!   real :: year
!   real :: day
!   real :: hour
   real :: time
   real :: realday
   real :: sols
   real :: season
   real :: dec
   real :: sindec
   real :: cosdec
   real :: tandec
   real :: coshr
end module 



module donor
   implicit none
   real :: etmass
   real :: hflux
   real :: roff
   real :: zlwup
   real :: drag
end module 



module rause
   implicit none
   real :: z0
   real :: d
   real :: rbc
   real :: rdc
end module 



module aerorx
   implicit none
   real :: ra
   real :: rb
   real :: rd
   real :: rbbest              ! adicionado por tatsch
end module 



module grads
   implicit none
   real :: tgs
   real :: ta
   real :: ea
   real :: etc
   real :: etgs
   real :: getc
   real :: getgs
   real :: u2
   real :: ustar
end module 



module radabs
   implicit none
   real :: albedo(2,2,2)
   real :: radfac(2,2,2)
   real :: radt(2)
   real :: thermk
   real :: exrain
   real :: tgeff
end module 



module surfrs
   implicit none
   real :: rst
   real :: rstfac(4)
   real :: rsoil
   real :: cog1
   real :: cog2
   real :: hr
   real :: fc
   real :: fg
end module 



module hydrol
   implicit none
   real :: satcap(2)
   real :: wc
   real :: wg
   real :: canex
   real :: areas
   real :: satcap_c
   real :: satcap_g
end module 



module stores
   implicit none
   real :: ccx
   real :: cg
   real :: csoil
end module 



module delts
   implicit none
   real :: dtc
   real :: dtg
   real :: dtd
   real :: dth
   real :: dqm
end module 



module carbio
   implicit none
   real :: assimn
   real :: respc
   real :: respg
   real :: pco2i
   real :: gsh2o
end module 



module flux
   implicit none
   real :: ec
   real :: eg
   real :: hc
   real :: hg
   real :: chf
   real :: shf
   real :: gflux   ! adicionado por tatsch
   real :: ect
   real :: eci
   real :: egi
   real :: egs
   real :: ecmass
   real :: egmass
   real :: heaten
end module 



module snow
   implicit none
   real :: tsnow
   real :: rsnow
end module 



module initialv
   implicit none
   real :: tc_ini
   real :: tg_ini
   real :: td_ini
   real :: www_ini(3)
end module 



!!! Are these variables integer or real ones?
module readin
   implicit none
!!!      common/ readin/  iqcaws, zlwd, tprec, iqchyd, mevap, msensh, mustar    
   integer :: iqcaws
   integer :: iqchyd
   integer :: mevap
   integer :: msensh
   integer :: mustar
   real :: zlwd
   real :: tprec
end module 



module checkbal
   implicit none
   real :: totwb
end module 



module output
   implicit none
   real :: roff1
   real :: roff2
   real :: roff3
   real :: roff4
   real :: gwsoil
   real :: otest1
   real :: otest2
   real :: otest3
   real :: roffGA
   real :: roffq3g
   real :: roffq3g2   ! tatsch 11 set 2011
end module 



!!! Are these variables integer or real ones?
module temp
   implicit none
!!!      common/ temp/  item01,item02, GWdep
   integer :: item01
   integer :: item02
   real :: GWdep
end module 



module suroff
   implicit none
   real :: surdep
   real :: finfil
end module 



module preccoff
   implicit none
   real :: app
   real :: bpp
   real :: cpp
end module 



module irrgrid
   implicit none
   integer :: idirr
end module 


module precision
!-------------------------------------------------------------------------------
! Purpose:
!	Define the precision to use for floating point and integer operations
!	throughout the model.
!-------------------------------------------------------------------------------
  integer, parameter :: r4 = selected_real_kind(5)
  integer, parameter :: r8 = selected_real_kind(12)
  integer, parameter :: i8 = selected_int_kind(13)
end module precision

MODULE PhysicalConstants

!=======================================================================
! physical constants 
!=======================================================================

  use precision
  IMPLICIT NONE

  public
  real(r8), parameter :: denice = 917.      ! density of ice [kg/m3]
  real(r8), parameter :: denh2o = 1000.     ! density of liquid water [kg/m3]
  real(r8), parameter :: cpliq  = 4188.     ! Specific heat of water [J/kg-K]
  real(r8), parameter :: cpice  = 2117.27   ! Specific heat of ice [J/kg-K]
  real(r8), parameter :: cpair  = 1004.64   ! specific heat of dry air [J/kg/K]
  real(r8), parameter :: hfus   = 0.3336e6  ! latent heat of fusion for ice [J/kg]
  real(r8), parameter :: hvap   = 2.5104e6  ! latent heat of evap for water [J/kg]
  real(r8), parameter :: hsub   = 2.8440e6  ! latent heat of sublimation [J/kg]
  real(r8), parameter :: tkair  = 0.023     ! thermal conductivity of air [W/m/k]
  real(r8), parameter :: tkice  = 2.290     ! thermal conductivity of ice [W/m/k]
  real(r8), parameter :: tkwat  = 0.6       ! thermal conductivity of water [W/m/k]
  real(r8), parameter :: tfrz   = 273.16    ! freezing temperature [K]
  real(r8), parameter :: rgas   = 287.04    ! gas constant for dry air [J/kg/K]
  real(r8), parameter :: roverg = 4.71047e4 ! rw/g = (8.3144/0.018)/(9.80616)*1000. mm/K
  real(r8), parameter :: rwat   = 461.296   ! gas constant for water vapor [J/(kg K)]
  real(r8), parameter :: grav   = 9.80616   ! gravity constant [m/s2]
  real(r8), parameter :: vonkar = 0.4       ! von Karman constant [-]
  real(r8), parameter :: stefnc = 5.67e-8   ! Stefan-Boltzmann constant  [W/m2/K4]

END MODULE PhysicalConstants

