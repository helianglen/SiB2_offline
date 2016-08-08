!
!#######################################################################
!
!     PURPOSE:
!
!     Read the static and dynamic parameters used in SiB2. The data include 
!     those for vegetation and soil and will be used for the calculation 
!     of derived aerodynamic parameters of SiB2. 
!
!     Included in subroutine SiB2derive_trans, SiB2derive_lai
!
!#######################################################################

subroutine sib2data

   use morphology1
   use morphology2, only : &
      z2_v, z1_v, zc_v, vcover_v, chil_v, leafw_v, leafl_v, &
      sodep_v, rootd_v,rootex_v
   use morphology3
   use morphology4
   use optical
   use physiology1
   use physiology2
   use soil



   implicit none



   !
   !#######################################################################
   !
   !     Definition of vegetation types:
   !     SiB2 type           Name		      ARPS type
   !
   !     1           Broadleaf-evergreen trees           7, 8  
   !     2           Broadleaf-deciduous trees		  6
   !     3           Broadleaf and needleleaf trees 
   !     4           Needleleaf-evergreen trees
   !     5           Needleleaf-deciduous trees
   !     6           Short vegetation/C4 grass	       1, 2, 3, 4, 5,11
   !     7           Broadleaf shrubs with bare soil
   !     8           Dwarf trees and shrubs		5, 12
   !     9           Agriculture/C3 grassland	        10
   !
   !#######################################################################
   !

   ! ************************
   ! Morphological properties
   ! ************************

   ! Type 1: time and biomes invariant parameters
   namelist /morpho1_nml/ &
      z0s_cst, &     ! Ground roughness length (m) 
      g4_cst, &      ! Transition height factor for mom. transfer
      dsfc_cst, &    ! Depth of surface layer (m)
      g1_cst         ! km(actual) : km(log-linear) at z2

   ! Type 2: time-invariant biomes dependent parameters
   namelist /morpho2_nml/ &
      z1_v, &        ! Canopy-base height (m)
      zc_v, &        ! Inflection height for leaf-area density (m)
      z2_v, &        ! Canopy-top height (m)
      vcover_v, &    ! Canopy cover fraction
      chil_v, &      ! Leaf area distribution factor
      leafw_v, &     ! Leaf width (m)
      leafl_v, &     ! Leaf length (m)
      sodep_v, &     ! Total depth oh three moisture layers (m)
      rootd_v, &        ! Rooting depth (m)
      rootex_v
   ! Type 3: time and biomes dependent variables
   !
   ! These variables should be derived from mophological and NDVI before 
   ! calling SiB2 and only need to calculate once for a not long-term
   ! simulation since their time-scale is quite long.
   ! They should be defined as permanent 2D arrary (nx,ny) in ARPS.
   namelist /morpho3_nml/ &
      ha_var, &
      z0d_var, &
      dd_var, &
      g2_var, &
      g3_var, &
      cc1_var, &
      cc2_var, &
      corb1_var, &
      corb2_var, &
      zlt_var, &
      green_var, &
      fparc_var, &
      gmudmu_var

   ! Type 4: location dependent variables
   namelist /morpho4_nml/ &
      slope_cst



   ! ******************
   ! Optical properties
   ! ******************
   namelist /optical_nml/ &
      reflv_v, &        ! Live leaf reflectance to visible
      refdv_v, &        ! Dead leaf reflectance to visible
      refln_v, &        ! Live leaf reflectance to near IR
      refdn_v, &        ! Dead leaf reflectance to near IR
      tranlv_v, &       ! Live leaf transmittance to visible
      trandv_v, &       ! Dead leaf transmittance to visible
      tranln_v, &       ! Live leaf transmittance to near IR
      trandn_v, &       ! Dead leaf transmittance to near IR
      sorefv_v, &       ! Soil reflectance to visible
      sorefn_v          ! Soil reflectance to near IR



   ! ************************
   ! Physiological properties
   ! ************************

   ! Type 1: time and biomes invariant parameters
   namelist /physio1_nml/ &
      shti_cst, &       ! Slope of high temp. inhibition func.(K-1)
      slti_cst, &       ! Slope of low temp. inhibition func. (K-1)
      trda_cst, &       ! Slope of high temp. inhibition function 
                        ! for respriation (K-1)
      trdm_cst, &       ! One-halp point of high temp. inhibition 
                        ! function for respriation (K)
      trop_cst, &       ! Temperature ceof. in GS-A model (K) 
                        ! using calculate qt=(Tc-Trop)/10
      btheta_cst        ! Coupling parameter for wp and ws

   ! time-invariant biomes dependent parameters
   namelist /physio2_nml/ &
      vmax0_v, &        ! Maximum rubsico capacity of top leaf (mol / m2 s)
      effcon_v, &       ! Intrinsic quantum efficiency (mol / mol)
      gradm_v, &        ! Conductance_photosynthesis slope parameter (mol / m2 s)
      binter_v, &       ! Minimum stamotal conductance (i.e., Conductance_photosynthesis intercept)
      atheta_v, &       ! Photosynthesis coupling coeffiicent for wc and we
                        ! Note: Photosynthesis coupling coeffiicent for wp,wc and we is independent
                        ! of vegetation and thus is a constant
      hhti_v, &         ! One half point of high temperature inhibition function (K)
      hlti_v, &         ! One half point of low temperature inhibition function (K)
      phc_v, &          ! One half critical leaf-water potential limit (m)
      respcp_v          ! Leaf respiration fraction of vmax0



   !
   !#######################################################################
   !
   !     Definition of soil types
   !     SiB2 type         Name                      ARPS type
   !
   !     1           Sand -> Loamy sand              1, 2 
   !     2           Sandy loam                      3
   !     3           Loam, silt loam                 4, 5 
   !     4           Clay loam -> sandy clay loam    6,7 
   !     5           Clay, Clay loam                 8,9,10,11
   !     6           Ice                             1
   !     7           Organic                         5
   !
   !#######################################################################
   !

   ! ***************
   ! Soil parameters
   ! ***************

   namelist /soil_nml/ &
      bee_s, &          ! Soil wetness exponent
      phsat_s, &        ! Soil tension at saturation (m)
      satco_s, &        ! Hydraulic conductivity at saturation (m/s)
      poros_s           ! Soil porosity  (m3/m3)



   ! Read parameters from configuration file
   !open(unit=1049, file='SIB2/sib2data.par', status='old')
   open(unit=1049, file='sib2data.par', status='old')
   read(1049, nml=morpho1_nml)
   read(1049, nml=morpho2_nml)
!   read(1049, nml=morpho3_nml)
!   read(1049, nml=morpho4_nml)
   read(1049, nml=optical_nml)
   read(1049, nml=physio1_nml)
   read(1049, nml=physio2_nml)
   read(1049, nml=soil_nml)
   close(unit=1049, status='keep')

end subroutine sib2data
