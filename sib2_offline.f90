!
! Program SIB2 (offline)
!
program sib2_offline

   use constants, only : pi, twopi, decmax

   use common_inc, only : icho1, icho2, icho3, icho4, icho5, icho6, &
      iout, iout1, iout2, iout3, iout4
!   use input_parms, only : startyear   ! ???

   use const   ! pie, g, asnow, clai, cpair, cw, epsfac, kappa, rhoair, 
   ! snomel, stefan, tf, vkc, rcp, timcon, psy, hlat, snofac
   use gridij   ! ivtype, istype
!!!   use sib2river_inc, only : isnow, ipbl, varcal_para
   use vdyijt, only : zlt, green
   use caerod, only : ha, g2, g3, corb1, corb2, zwind, zmet
   use vderiv, only : z0d, dd, cc1, cc2, gmudmu
   use initialv   ! tc_ini, tg_ini, td_ini, www_ini
   use steps, only : itrunk, ilw, dtt, iter
   use soils, only : poros, phsat, satco, bee, slope, slpp, decay,anik!, speyield
   use preccoff   ! app, bpp, cpp
   use site, only : zlong, zlat
   use stepv   ! tc, tg, td, capac, snoww, www
   use suroff, only : surdep
   use atmos, only : swdown, rnetm, em, tm, um
   use readin, only : tprec, zlwd !, mustar
   use irrgrid, only : idirr
   use temp, only : gwdep
!   use govern, only : year, month, day, hour, time, realday, sols
   use govern, only : year, month, day, hour, &
      sols, realday, time, season, dec, sindec, cosdec, tandec, coshr
!   use donor, only : etmass   ! TEST
   use morphology2
   use morphology1
   use optical
   use vstate
   use physiology1
   use physiology2
   use multicamada, only : nlayers,extfrac,dzmultilayer,depth,hr_proc,satcoz_d
   use hydrol , only : satcap_c,satcap_g
   implicit none


   integer, parameter :: in2 = 10
   integer, parameter :: out1 = 20
   integer, parameter :: ichmet = 7
   integer, parameter :: iu = 8
   integer :: date1, date2, nymd, jul, dayloc, sta_id
   integer :: i, j, ierr
   integer :: isnow, ipbl   ! Taken off from 'sib2river_inc' module
   real :: varcal_para(9,70,9)   ! Taken off from 'sib2river_inc' module
   character (len=120) :: infile2,dir_out
   character (len= 4) :: id
   integer, external :: cal2jul
   real :: vcov, greeness, greenpar
   real :: reflv,refdv,refln,refdn,tranlv,trandv,tranln,trandn,sorefv,sorefn
   real :: z0s, g4, g1
   real :: z2_c, z1_c, zc_c, chil_c, leafw_c, leafl_c
   real :: effcon_c, gradm_c, binter_c, respcp_c, atheta_c, btheta_c, rootd_c
   real :: phc_c, trop_c, slti_c, hlti_c, hhti_c, vmax_c,sodep_c
   real :: beta,froot_norm,cumextfrac,depwww
    integer :: kzat_pre
   namelist /sib_run/ date1, date2, infile2, dtt, itrunk, ilw
   namelist /sib_invars/ ivtype, istype, isnow, ipbl, idirr, decay, &
      zlong, zlat, poros, phsat, satco, bee, slope, slpp, tc, tg, td, &
      capac, snoww, www, surdep, gwdep, app, bpp, cpp
   namelist /sib_initialv/ tc_ini, tg_ini, td_ini, www_ini
   namelist /sib_caerod/ zwind, zmet
   namelist /sib_vderiv/ gmudmu

   ! Use the following subroutine to read a namelist input file or the read 
   ! command to read all input data from console
   !call read_config_file
   ! dir_out='./'
   ! id='999'

   read(*,*) date1, date2, infile2, dtt, itrunk, ilw, ivtype, istype, isnow,   &
      ipbl, idirr, decay, zlong, zlat, poros, phsat, satco, bee, slope, slpp,  &
      tc, tg, td, capac(1), capac(2), snoww(1), snoww(2), www(1),www(2),www(3),&
      surdep, gwdep, app, bpp, cpp, tc_ini,                                   &
      tg_ini, td_ini, www_ini(1), www_ini(2), www_ini(3), zwind, zmet, gmudmu,&
      vcov,greenpar, id, dir_out,                                             &
      ! propriedades opticas
      reflv,refdv,refln,refdn,tranlv,trandv,tranln,trandn,sorefv,sorefn,      & 
      ! propriedades morfologicas
      z0s, g4, g1,                                                            &
      z2_c, z1_c, zc_c, chil_c, leafw_c, leafl_c,                             &
      effcon_c, gradm_c, binter_c, respcp_c, atheta_c, btheta_c, rootd_c,     &
      phc_c, trop_c, slti_c, hlti_c, hhti_c, vmax_c,                          &
      sodep_c,rootex,anik,satcap_c,satcap_g,                                  &
      nlayers,beta,dzmultilayer,hr_proc,kzat_pre
      
      
   ! Write input data on screen
   write(*,nml=sib_run)
   write(*,nml=sib_invars)
   write(*,nml=sib_initialv)
   write(*,nml=sib_caerod)
   write(*,nml=sib_vderiv)
    
   call sib2data
   call cntrol(icho1, 0)
   ! Atualizando valores das variáveis 
   ! no proceso de calibração
   if(ivtype.eq.2) then
   vcover_v(ivtype)= vcov
   reflv_v(ivtype)= reflv 
   refdv_v(ivtype)= refdv
   refln_v(ivtype)= refln
   refdn_v(ivtype)= refdn
   tranlv_v(ivtype)= tranlv
   trandv_v(ivtype)= trandv
   tranln_v(ivtype)= tranln
   trandn_v(ivtype)= trandn
   sorefv_v(ivtype)= sorefv
   sorefn_v(ivtype)= sorefn*sorefv
   z0s_cst = z0s
   g1_cst = g1
   z2_v(ivtype) = z2_c
   z1_v(ivtype) = z1_c * z2_c ! 0 < z1 < 0.8 multiplicador
   zc_v(ivtype) = z1_c*z2_c + zc_c * (z2_c - z1_c * z2_c) ! 0 <=zc <= 1
   g4_cst = g4 
   chil_v(ivtype) = chil_c
   leafw_v(ivtype) = leafw_c
   leafl_v(ivtype) = leafl_c
!  Parametros no modelo de fotosintesis e condutancia 
   effcon_v(ivtype) = effcon_c
   gradm_v(ivtype) = gradm_c
   binter_v(ivtype) = binter_c
   respcp_v(ivtype) = respcp_c
   atheta_v(ivtype) = atheta_c
   btheta_cst = btheta_c
   rootd_v(ivtype) = rootd_c
   phc_v(ivtype) = phc_c
   trop_cst = trop_c
   slti_cst = slti_c
   hlti_v(ivtype) = hlti_c
   hhti_v(ivtype) = hhti_c * hlti_c
   vmax0_v(ivtype) = vmax_c
   ! soil 
   sodep_v(ivtype) = (int( (rootd_c + sodep_c)*10 ) + 0.0) /10.0
!   rootex_v(ivtype) = rootex
   end if
   

   !stop
   do j = 1, 9
      do i = 1, 70
         zlt = 0.1 * real(i)
         call derive_trans(j, zlt, rhoair, &
            ha, z0d, dd, g2, g3, cc1, cc2, corb1, corb2)
         varcal_para(1,i,j) = ha
         varcal_para(2,i,j) = z0d
         varcal_para(3,i,j) = dd
         varcal_para(4,i,j) = g2
         varcal_para(5,i,j) = g3
         varcal_para(6,i,j) = cc1
         varcal_para(7,i,j) = cc2
         varcal_para(8,i,j) = corb1
         varcal_para(9,i,j) = corb2
!         if(id == '001' .and. j .eq. 2) write(1765, '(i4,9(f15.9))')i,ha,z0d,dd,g2,g3,cc1,cc2,corb1,corb2
      end do
   end do
    
   ! stop
   ! Computes the location on the year of the initial date
   year = date1 / 1000000
   month = mod(date1, 1000000) / 10000
   day = mod(date1, 10000) / 100
   jul = cal2jul(year, month, day)
   dayloc = jul - cal2jul(year, 1, 1) + 1
   sols = (4141.0 / 24.0) + 0.25 * real(mod(year + 3, 4))

   open(unit=in2, file=trim(infile2), status='old')
   read(in2,*)
!   open(unit=200, file='output.dat', form='unformatted', access='stream', status='unknown')
   open(unit=iout1, file=trim(dir_out)//"sib2diag"//trim(id)//".txt", &
   status='unknown')
   !
   ! Main loop
   print*,trim(infile2)
   print*,trim(dir_out)//"sib2diag"//trim(id)//".txt"
   ! stop

   
   if(dzmultilayer .ne. 0.0)then
    !determinando numero de camadas a parir da espessura definida
    nlayers = int( (sodep_v(ivtype)-dsfc_cst) / dzmultilayer )
    ! atualizando a profundidade do solo para fechamento de 
    ! numero de camadas e espessura
    sodep_v(ivtype) = (nlayers * dzmultilayer) + dsfc_cst
   else
     ! se não eh definida espessura então o solo éh dividido 
     ! em nlayers camadas
    dzmultilayer = (sodep_v(ivtype)-dsfc_cst)/(nlayers+0.0)
   endif   
   
      extfrac = 0.0
      cumextfrac = 0.0
      froot_norm = 1 - beta**((sodep_v(ivtype)-dzmultilayer)*100)
    Do i=1, nlayers-1
      if(i .eq. 1)then
	depth(i) = (dsfc_cst + (dzmultilayer))   ! cm
	extfrac(i) =  1 - beta**(depth(i)* 100)
      else 
	depth(i) = depth(i-1) + dzmultilayer
	extfrac(i) = (1 - beta**(depth(i)* 100)) - (1 - beta**(depth(i-1)* 100))
      end if
      extfrac(i) = extfrac(i)/froot_norm
      cumextfrac =  cumextfrac + extfrac(i) 
    END DO
    
    if(kzat_pre .eq. 1)then
	   OPEN(9890,file='best_ksatz.txt',status='old')
      DO i=1,nlayers
	read(9890,*) satcoz_d(i)
      ENDDO
    
    else
	    depwww=dsfc_cst+(dzmultilayer/2.0)
      DO i=1,nlayers
	satcoz_d(i) = satco*exp(-decay * depwww)
	depwww=depwww+dzmultilayer
      ENDDO
    endif
    
   iter = 0
   do
!      read(in2,*,iostat=ierr) sta_id, nymd, swdown, rnetm, em, tm, um, tprec, zlt, green, zlwd
      read(in2,*,iostat=ierr) sta_id, nymd, swdown, rnetm, em, tm, um, tprec, zlt , greeness, zlwd
      if (ierr /= 0 .or. nymd > date2) exit

      if (nymd >= date1) then
         iter = iter + 1

         ! Extract YEAR, MONTH, DAY and HOUR from NYMD variable
         year = nymd / 1000000
         month = mod(nymd, 1000000) / 10000
         day = mod(nymd, 10000) / 100

         !!! TIRAR o '-1' DA LINHA ABAIXO POIS O ARQUIVO VAI DE 1 A 24 E DEVERIA IR DE 0 A 23 
         hour = mod(nymd, 100) - 1

         hour = mod(hour, 24)   ! Useful only if file have daily data

         ! Yearly update
         if (month == 1 .and. day == 1 .and. hour == 0) then
            dayloc = 1
            sols = (4141.0 / 24.0) + 0.25 * real(mod(year + 3, 4))
         end if
	 ! if(dayloc >31 ) dayloc = 1 !Se tirar o comentario perde a correção ao sunang
         realday = real(dayloc) + real(hour) / 24.0
         season = (realday - sols) / 365.2
         dec = decmax * cos(twopi * season)
         sindec = sin(dec)
         cosdec = cos(dec)
         tandec = tan(dec)
         time = real(hour) - 0.5
         coshr = cos(-pi + time / 24.0 * twopi)   ! Melhorar isto

         ! The following three variables must be read from data file
!         zlt = 3.54
        
         if(greenpar .lt. 0 ) green = amin1(max(0.19*zlt+0.0627, 0.4), 0.99)
         if(greenpar .gt. 1 ) green = zlt/7.0
         if(greenpar .ge. 0 .and. greenpar .le. 1) green = greeness*(0.98-greenpar) + greenpar
!         if(ivtype.eq.2) green = amin1(max(0.19*zlt+0.0627, 0.4), 0.99)
!         zlwd = 4.903E-3 / 24.0 * (tm ** 4) * (0.66 + 0.039 * sqrt(em / 10.0)) / 3600.0
         i = nint(10.0 * zlt)
         i = min(70, max(i, 1))

         ! Actually the following subroutine do nothing
!         call const2

         ! For use with 'SIB2sub.f' file (FORTRAN 77 version)
!         call vegpar
!         call varcal(varcal_para)
         ! For use with 'sib2_sub.f90' file (Fortran 90 version)
         call vegpar(ivtype)
         call varcal(ivtype, varcal_para(1,i,ivtype))
         ! The first three arguments are NOT used in the subroutine
!         call driver(iu, icho2, ichmet, isnow, nymd)
!         call driver(isnow, nymd)
         call driver(isnow)
         call balan(1)
         call inter2(idirr)
         call rada2
         call begtem
         call endtem(ipbl)
         call updat2
         call balan(2)
!         write(777,*)'sib2off.f90:206 tg:',tg, 'depois de balan2',iter
!         call outer(iout, iout1, iout2, iout3, iout4, nymd)
         call outer(iout1, nymd)

!!!      TEST
!!!         WRITE(100,'(I12,4F15.7)') NYMD, REALDAY, TIME, SUNANG, ETMASS

         if (hour == 23) dayloc = dayloc + 1
      end if   ! (nymd>=date1)
   end do   ! main loop

   close(unit=in2, status='keep')
   close(unit=200, status='keep')

contains

   subroutine read_config_file
      implicit none

      open(unit=1, file='sib2_offline.par', status='old')
      read(1, nml=sib_run)
      read(1, nml=sib_invars)
      read(1, nml=sib_initialv)
      read(1, nml=sib_caerod)
      read(1, nml=sib_vderiv)
      close(unit=1, status='keep')
   end subroutine read_config_file



end program sib2_offline
