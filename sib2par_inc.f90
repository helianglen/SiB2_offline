!
! 'sib2par_inc.f90' file
!
! This file replaces the 'SiB2PAR.INC' file
!

module sib2_parms
   implicit none

   integer, parameter :: nv = 9   ! Number of vegetation types
   integer, parameter :: ns = 7   ! Number of soil types
end module sib2_parms



! ================================================
! Vegetation : static, dynamic, derived parameters                       
!  
! Morphological properties
! ================================================

! Type 1: time and biomes invariant parameters
! This module replace /morphology1/ common block
module morphology1
   implicit none

   real :: z0s_cst, g4_cst, dsfc_cst, g1_cst
end module morphology1



! Type 2: time-invariant biomes dependent parameters
! This module replace /morphology2/ common block
module morphology2
   use sib2_parms, only : nv
   implicit none

   real :: z2_v(nv), z1_v(nv), zc_v(nv), vcover_v(nv), chil_v(nv), &
            leafw_v(nv), leafl_v(nv), sodep_v(nv), rootd_v(nv),rootex_v(nv), &
            stem_v(nv)
   real :: laimax_v(nv), lais_v(nv), ndvi98_v(nv), ndvi5_v(nv), &
            fcl_v(nv)
end module morphology2

! Type 3: time and biomes dependent variables (derived)
! This module replace /morphology3/ common block
module morphology3
   implicit none

   real :: ha_var, z0d_var, dd_var, &
      g2_var, g3_var, cc1_var, cc2_var, corb1_var, corb2_var, &
      zlt_var, green_var, fparc_var, gmudmu_var
end module morphology3



! Type 4: location dependent variables
! This module replace /morphology4/ common block
module morphology4
   implicit none

   real :: slope_cst 
end module morphology4



! ==================
! Optical properties
! ==================

! This module replace /optical/ common block
module optical
   use sib2_parms, only : nv
   implicit none
   real :: reflv_v(nv),  refdv_v(nv),  refln_v(nv),  refdn_v(nv), &
      tranlv_v(nv), trandv_v(nv), tranln_v(nv), trandn_v(nv), &
      sorefv_v(nv), sorefn_v(nv)
end module optical



! ========================
! Physiological properties
! ========================

! Type 1: time and biomes invariant parameters
! This module replaces /physiology1/ common block
module physiology1
   implicit none

   real :: shti_cst, slti_cst, trda_cst, trdm_cst, trop_cst, &
      btheta_cst
end module physiology1



! Type 2: time-invariant biomes dependent parameters
! This module replaces /physiology2/ common block
module physiology2
   use sib2_parms, only : nv
   implicit none

   real :: vmax0_v(nv), effcon_v(nv), gradm_v(nv), binter_v(nv), &
      atheta_v(nv),hhti_v(nv),hlti_v(nv),phc_v(nv), respcp_v(nv)
end module physiology2



! ========================
! Soil : static parameters                       
! ========================
! This module replaces /soil/ common block

module soil
   use sib2_parms, only : ns
   implicit none

   real :: bee_s(ns), phsat_s(ns), poros_s(ns), satco_s(ns)
end module soil







!!!c
!!!c
!!!c     ##################################################################
!!!c     ##################################################################
!!!c     ######                                                      ######
!!!c     ######                     SiB2PAR.INC                      ######
!!!c     ######                     Developed by                     ######
!!!c     ######     River and Environmental Engineering Laboratory   ######
!!!c     ######                University of Tokyo                   ######
!!!c     ######                                                      ######
!!!c     ##################################################################
!!!c     ##################################################################
!!!c
!!!c
!!!c
!!!c#######################################################################
!!!c
!!!c     PURPOSE:
!!!c
!!!c     Include file 'SiB2par.inc', used to calculate the derived
!!!c     parameter in SiB2
!!!c
!!!C#######################################################################
!!!c      
!!!c
!!!c     AUTHOR: K. Yang
!!!c                                                                               
!!!c     MODIFICATION HISTORY:
!!!c     01/21/01 (K. Yang)
!!!c     Reform full documentation 
!!!c
!!!c#######################################################################
!!!
!!!      integer     nv            ! Number of vegetation types     
!!!      integer     ns            ! Number of soil types     
!!!      PARAMETER (nv = 9, ns = 7)
!!!c
!!!c######################################################################
!!!c
!!!c        Vegetation : static, dynamic, derived parameters                       
!!!c
!!!c######################################################################
!!!c 
!!!c
!!!C######################################################################
!!!c
!!!c     Morphological properties
!!!c
!!!C######################################################################
!!!c
!!!c     Type 1: time and biomes invariant parameters
!!!
!!!      real    z0s_cst, g4_cst, dsfc_cst, g1_cst
!!!      COMMON /morphology1/ z0s_cst, g4_cst, dsfc_cst, g1_cst
!!!
!!!c     Type 2: time-invariant biomes dependent parameters
!!!
!!!      real    z2_v(nv), z1_v(nv), zc_v(nv), vcover_v(nv), chil_v(nv), 
!!!     &        leafw_v(nv), leafl_v(nv), sodep_v(nv), rootd_v(nv)
!!!      real    laimax_v(nv), lais_v(nv), ndvi98_v(nv), ndvi5_v(nv),
!!!     &        fcl_v (nv)
!!!
!!!      COMMON /morphology2/ 
!!!     &        z2_v, z1_v, zc_v, vcover_v, chil_v, 
!!!     &        leafw_v, leafl_v, sodep_v, rootd_v,
!!!     &        laimax_v, lais_v, ndvi98_v, ndvi5_v, fcl_v
!!!
!!!c     Type 3: time and biomes dependent variables (derived)
!!!
!!!      real    ha_var, z0d_var, dd_var,  
!!!     &        g2_var, g3_var, cc1_var, cc2_var, corb1_var, corb2_var,
!!!     &        zlt_var, green_var, fparc_var, gmudmu_var
!!!     
!!!      COMMON /morphology3/
!!!     &        ha_var,z0d_var,dd_var,
!!!     &        g2_var,g3_var,cc1_var, cc2_var,corb1_var, corb2_var,
!!!     &        zlt_var, green_var, fparc_var, gmudmu_var
!!!
!!!
!!!
!!!c     Type 4: location dependent variables
!!!
!!!      real    slope_cst 
!!!
!!!      COMMON /morphology4/ slope_cst
!!!c
!!!C######################################################################
!!!c
!!!c     Optical properties
!!!c
!!!C######################################################################
!!!c
!!!
!!!      real    reflv_v(nv),  refdv_v(nv),  refln_v(nv),  refdn_v(nv), 
!!!     &        tranlv_v(nv), trandv_v(nv), tranln_v(nv), trandn_v(nv), 
!!!     &        sorefv_v(nv), sorefn_v(nv)
!!!
!!!      COMMON /optical/
!!!     &        reflv_v,  refdv_v,  refln_v,  refdn_v, 
!!!     &        tranlv_v, trandv_v, tranln_v, trandn_v, 
!!!     &        sorefv_v, sorefn_v
!!!c
!!!C######################################################################
!!!c
!!!c     Physiological properties
!!!c
!!!C######################################################################
!!!c
!!!
!!!c     Type 1: time and biomes invariant parameters
!!!
!!!      real    shti_cst, slti_cst, trda_cst, trdm_cst, trop_cst, 
!!!     &        btheta_cst       
!!!
!!!      COMMON /physiology1/
!!!     &        shti_cst, slti_cst, trda_cst, trdm_cst, trop_cst, 
!!!     &        btheta_cst       
!!!
!!!c     Type 2: time-invariant biomes dependent parameters
!!!
!!!      real    vmax0_v(nv), effcon_v(nv), gradm_v(nv), binter_v(nv), 
!!!     &        atheta_v(nv),hhti_v(nv),hlti_v(nv),phc_v(nv), respcp_v(nv)
!!!
!!!      COMMON /physiology2/
!!!     &        vmax0_v, effcon_v, gradm_v, binter_v, 
!!!     &        atheta_v,hhti_v,hlti_v,phc_v, respcp_v
!!!
!!!c
!!!c######################################################################
!!!c
!!!c        Soil : static parameters                       
!!!c
!!!c######################################################################
!!!c 
!!!
!!!      real    bee_s(ns), phsat_s(ns), poros_s(ns), satco_s(ns)
!!!      COMMON /soil/ bee_s, phsat_s, poros_s, satco_s

