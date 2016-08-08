!c=======================================================================
!c
!c                  SUBROUTINES
!c
!c=======================================================================
!c
!!!      SUBROUTINE driver( iu, icho2, ichmet, isnow, nymd )
!!!      SUBROUTINE driver(icho2, isnow, nymd)
subroutine driver(isnow)
!c-----------------------------------------------------------------------
!c     subroutines called : radc2
!c-----------------------------------------------------------------------
   use stepv, only : snoww, tc, tg
   use atmos, only : radn, cloud, em, ppc, ppl, sunang, swdown, tm, um
   use steps, only : iter
   use readin, only : mustar, tprec, zlwd
   use govern, only : sindec, cosdec, coshr
   use site, only : zlat

   implicit none

   ! Dummy arguments
   integer, intent(in) :: isnow

   ! Local variables
   real :: x, e, cold, difrat, rhair, ptot, ustarm, vnrat
!    real :: cloud_old
!c-----------------------------------------------------------------------

   e(x) = exp( 21.18123 - 5418.0 / x ) / 0.622

!c-----------------------------------------------------------------------
!c
!c     isnow = 0 : conventional run with received met data.
!c
!c     isnow = 1 : shock test with initial snow dump and freezing
!c                 temperatures at beginning of run warming to
!c                 normal over  5-day period.
!c
!c-----------------------------------------------------------------------
    
   ptot = 0.0

   if (isnow /= 0) then
      if (iter <= 1) then
         tc = 270.0
         tg = 270.0
         snoww(2) = 0.1
      end if
      cold = max (0.0, (120.0 - (1.0 * iter)) / 120.0)
      rhair = em / e(tm)
      tm = tm * (1.0 - cold) + (tm - 30.0) * cold
      em = e(tm) * rhair
      if (em < 0.0) em = 0.1
   end if

!c      if(nymd <= 85082523)then

   um = max(um, 0.25)
   ustarm = mustar / 100.0
   swdown = max(swdown, 0.1)
!!!         ihour = mod(nymd,100)
   ptot = ptot + tprec
   ppl = tprec
   ppc = tprec - ppl

!!!         call radc2_2
   call radc2_2(zlat, sindec, cosdec, coshr, sunang)

!   cloud = (1160.0 * sunang - swdown) / (963.0 * sunang)
   cloud = 2.33 - 3.33*(  swdown / (1160.0 * sunang) )
   cloud = max(cloud, 0.0)
   cloud = min(cloud, 1.0)
   cloud = max(0.38, cloud)
   
!   cloud_old = (1160.0 * sunang - swdown) / (963.0 * sunang)
!   cloud_old = max(cloud_old, 0.0)
!   cloud_old = min(cloud_old, 1.0)
!   cloud_old = max(0.58, cloud_old)
!   write(8987,'(i8,2(1x,f6.4),2(1x,f13.6))')iter, cloud,cloud_old, swdown, 1160.0 * sunang
   
   difrat = 0.0604 / (sunang - 0.0223) + 0.0683
   if (difrat < 0.0) difrat = 0.0
   if (difrat > 1.0) difrat = 1.0

   difrat = difrat + (1.0 - difrat) * cloud
   vnrat = (580.0 - cloud * 464.0) / ((580.0 - cloud * 499.0) &
      + (580.0 - cloud * 464.0))

   radn(1,1) = (1.0 - difrat) * vnrat * swdown
   radn(1,2) = difrat * vnrat * swdown
   radn(2,1) = (1.0 - difrat) * (1.0 - vnrat) * swdown
   radn(2,2) = difrat * (1.0 - vnrat) * swdown
   radn(3,2) = zlwd
!!!         RETURN
!c      END IF
!!!1000  CONTINUE
!!!      WRITE(icho2, 900) iu, nymd, iout
!!!      WRITE(icho2, 900) nymd, iout
!!! 900   FORMAT(5x,'eof encountered for unit= ',i2,' eof date= ',i8,1x,i2)
!!!900   FORMAT(5x,'eof encountered. EOF date= ',i8,1x,i2)
!!!      STOP 'in driver'

   
end subroutine driver




!c=======================================================================

subroutine balan(iplace)

!c=======================================================================
!c
!c     energy and water balance check.
!c
!c-----------------------------------------------------------------------

   use atmos
   use checkbal
   use donor, only : etmass, hflux, roff
   use flux, only : chf, eci, ect, egi, egs, hc, heaten, hg, shf
   use gridij
   use output, only : &
      gwsoil, roff1, roff2, roff3, roff4, &
      roffga, roffq3g, roffq3g2
   use radabs, only : radt
   use readin
   use soils, only : poros, zdepth
   use steps, only : dtt, iter
   use stepv, only : snoww, www, capac
   use temp
   use multicamada
   implicit none

   ! Dummy arguments
   integer :: iplace,  i

   ! Local variables
   real :: cbal, emeter, endwb, errore, errorw, gbal, pmeter
   real :: zlhs, zrhs,soilwat

   if (iplace == 1) then

      etmass = 0.0
      roff   = 0.0
      roff1	 = 0.0
      roff2	 = 0.0
      roffga = 0.0 ! tatsch
      roff3	 = 0.0
      roffq3g = 0.0 ! tatsch
      roffq3g2 = 0.0 ! tatsch
      roff4	 = 0.0
      gwsoil = 0.0
      if(iter .eq. 1)then
      totwb = www(1) * poros * zdepth(1) &
         + www(2) * poros * zdepth(2) &
         + www(3) * poros * zdepth(3) &
         + capac(1) + capac(2) + snoww(1) + snoww(2)
      else
	      soilwat = 0.0
	      do i=1,nlayers
		soilwat = soilwat + dzwww(i) * poros * dzmultilayer
	      enddo
!          + www(2) * poros * zdepth(2) &
!          + www(3) * poros * zdepth(3) &
         totwb = www(1) * poros * zdepth(1) &
         + soilwat                       &
         + capac(1) + capac(2) + snoww(1) + snoww(2)
      endif
   else
	soilwat = 0.0
      do i=1,nlayers
	soilwat = soilwat + dzwww(i) * poros * dzmultilayer
      enddo
      
      endwb = www(1) * poros * zdepth(1) &
!          + www(2) * poros * zdepth(2) &
!          + www(3) * poros * zdepth(3) &
	   + soilwat                     &
         + capac(1) + capac(2) + snoww(1) + snoww(2) &
!c.......................................................................
!c! 2005/9/28 changed by tangqh@iis.u-tokyo.ac.jp
!c old    &        - (ppl+ppc)/1000. + etmass/1000. + roff
         - (ppl+ppc)/1000. + etmass/1000. + roff !-gwsoil
!c.......................................................................
!c        errorw= totwb - endwb
      errorw = totwb - endwb
      pmeter = (ppl + ppc) / 1000.0
      emeter = etmass / 1000.0

!       if (abs(errorw) > 0.0001) then
!       !print*,"water balance violation"
!          write(242,*) iter,item01,item02
!          write(242,*) swdown,zlwd,em,tm,um
!          write(242,*) tprec,istype,ivtype
!          write(242,900) nymd, totwb, endwb, errorw, &
!             www(1), (dzwww(i),i=1,nlayers), &
!             capac(1), capac(2), snoww(1), snoww(2), &
!             pmeter, emeter, gwsoil,                &
!             roff,roff1,roff2,roff3, roff4
! 900      format(/,10x,'** warning: water balance violation **  ',/, &
!             /,1x,'date ', i12, &
!             /,1x,'begin, end, diff ', 3(f10.7,1x), &
!             /,1x,'www,1+1-10       ', 11(f7.4,1x), &
!             /,1x,'capac,1-2        ', 2(f10.8,1x), &
!             /,1x,'snoww,1-2        ', 2(f10.8,1x), &
!            /,1x,'p, et, roff,roff2, gwsoil ', 3(f10.8,1x) &
!            / ,1x,' roff, roff1 , roff2 , roff3 , roff4 : ',  5(f10.8,1x))
!       else
!         !write(242,*) iter,item01,item02 ,"NO water balance violation"
!       endif

      cbal = radt(1) - chf - (ect + hc + eci) / dtt
      gbal = radt(2) - shf - (egs + hg + egi) / dtt - heaten / dtt
      zlhs = radt(1) + radt(2) - chf - shf
      zrhs = hflux + (ect + eci + egi + egs) / dtt + heaten / dtt

      errore= zlhs - zrhs

      if(abs(errore) > 1.) then
!         write(icho2,*) 'iter,inr,inc:',iter,item01,item02
!         write(icho2,*) 'swdown, zlwd, vapor pressure:',swdown,zlwd,em
!         write(icho2,*) 'temperature, wind speed:',tm,um
!         write(icho2,*) 'precipitation, soiltype, vtype',tprec,istype,ivtype
!         write(icho2,910) nymd, zlhs, zrhs, radt(1), radt(2), chf, shf, &
!            hflux, ect, eci, egi, egs, hc, hg, heaten, cbal, gbal
!910      format(/,10x,'** warning: energy balance violation **',/, &
!            /,1x,'date ', i8, &
!            /,1x,'rhs, lhs              ', 2g12.5, &
!            /,1x,'rn1, rn2, chf, shf, h ', 5g12.5, &
!            /,1x,'ect, eci, egi, egs    ', 4g12.5, &
!            /,1x,'hc        hg          ',  g12.5, 12x, g12.5, &
!            /,1x,'heaten, c-bal, g-bal  ', 3g12.5 )
      endif

   end if   ! (iplace==1)
  
end subroutine balan




!c=======================================================================

!!!      SUBROUTINE outer ( iout, iout1, iout2, iout3, iout4, nymd )
subroutine outer(iout1, nymd)

!c=======================================================================
!c
!c     output of results to files.
!c
!c-----------------------------------------------------------------------

   use aerorx, only : ra
   use atmos, only : em, swdown, tm, um, radn
   use carbio, only : assimn
   use const, only : hlat, stefan, tf
   use donor, only : etmass, hflux, roff, zlwup
   use flux, only : chf, eci, ect, egi, egs, gflux, shf
   use grads, only : tgs, ustar
   use gridij, only : istype, ivtype
   use hydrol, only : areas
   use output, only : &
      gwsoil, roff1, roff2, roff3, &
      roffga, roffq3g, roffq3g2
   use radabs, only : tgeff, thermk, radt
   use readin, only : tprec
   use site, only : salb
   use soilij, only : sodep, soref
   use soils, only : bee, phsat, poros, satco, slope, zdepth
   use steps, only : dtt, iter
   use stepv, only : www, tc, tg
   use temp, only : item01, item02
   use vderiv, only : dd, vmax0, z0d
   use vdyijt, only : fparc, green, zlt
   use vstate, only : gradm, vcover
   use multicamada 
   
   implicit none

   ! Dummy arguments
   integer, intent(in) :: iout1, nymd
   integer :: i , j

   ! Local variables
   real :: calbe, canil, elat, evapg, evapg1, evapg2
   real :: gcstor, radswa, radswd, radtot, tc4, tg4
   real :: trant, zlwf
   real :: radswVa,radswVd,calbe_par
   real :: swc01 , swc02 , fracd, swcsfc
!   real :: levels(6)
   real :: niveis(8)
   integer :: out_layers(8)
   
   integer :: wcount20 ,wcount01 ,wcount02
   real :: vwc20 , vwc01 , vwc02
   
   data niveis/0.1,0.2,0.5,0.8,1.0,1.5,2.0,2.5/ 
          
   trant  = ect / hlat
   canil  = eci  /hlat
!c! tatsch 10/23/2009 mudanca para output de componentes de evaporacao
!c! apenas adicionei sub variaveis evapg1 e evapg2
   evapg = (egs + egi) / hlat
   evapg1 = egs / hlat
   evapg2 = egi / hlat

   radswd = radn(1,1) + radn(1,2) + radn(2,1) + radn(2,2)
   radswa = (1.0 - salb(1,1)) * radn(1,1) + (1.0 - salb(1,2)) * radn(1,2) &
      + (1.0 - salb(2,1)) * radn(2,1) + (1.0 - salb(2,2)) * radn(2,2)
   radtot = radt(1) + radt(2)
!    write(234,'(i12,8(1x,f12.7))')nymd,radn(1,1),radn(1,2), radn(2,1), radn(2,2), &
!	       salb(1,1),salb(1,2),salb(2,1),salb(2,2)
!c!tatsch 07/02/2010
!c!  calbe: surface albedo (dimensionless)
   calbe = 1.0 - (radswa / radswd)
!c!   assimn in umol m-2 s-1
!c!      assimn * 1.e06

    radswVd = radn(1,1) + radn(1,2)
    radswVa = (1.0 - salb(1,1)) * radn(1,1) + (1.0 - salb(1,2)) * radn(1,2)
    calbe_par = 1.0 - (radswVa / radswVd)
    
   elat = etmass / dtt * hlat
   gcstor = chf + shf

   tgs = min(tf,tg) * areas + tg * (1. - areas)
   tc4 = tc  * tc  * tc  * tc
   tg4 = tgs * tgs * tgs * tgs
   zlwf =  tc4 * stefan * vcover * (1.0 - thermk) &
      + (1.0 - vcover * (1.0 - thermk)) * stefan * tg4
   tgeff = sqrt(sqrt(zlwf / stefan))
  
!!!!!!!****************************************************************
 !     Gerando colunas de conteúdo de água para as camadas:  
 !           0.0 - 1.00 m   
 !           0.0 - 2.00 m
	swcsfc = www(1) * poros * zdepth(1)
	swc01 = 0.0
	swc02 = 0.0
      do i = 1, nlayers
        if(depth(i) .gt. 0.0 .and. depth(i) .le. 0.2 )then
          swcsfc = swcsfc + dzwww(i) * poros * dzmultilayer 
        else if( depth(i).gt. 0.2 .and. depth(i) .lt. (0.2 + dzmultilayer))then
          fracd = (depth(i)-0.2)/dzmultilayer
          swcsfc = swcsfc + (1-fracd) * dzwww(i)*poros*dzmultilayer
          swc01 = swc01 + fracd * dzwww(i)*poros*dzmultilayer
	else if(depth(i) .gt. (0.2 + dzmultilayer) .and. depth(i).le. 1.0 )then
	  swc01 = swc01 + dzwww(i) * poros * dzmultilayer 
	else if(depth(i).gt. 1.0 .and. depth(i) .lt. (1.0 + dzmultilayer))then
	  fracd = (depth(i)-1.0)/dzmultilayer
	  swc01 = swc01 + (1-fracd) * dzwww(i)*poros*dzmultilayer
	  swc02 = swc02 + fracd * dzwww(i)*poros*dzmultilayer
	else if( depth(i) .gt. (1.0 + dzmultilayer) .and. depth(i) .lt. 2.0)then
	  swc02 = swc02 + dzwww(i)*poros*dzmultilayer
	else if( depth(i) .lt. (2.0 + dzmultilayer))then
	  fracd = 1.0- ((depth(i)-2.0)/dzmultilayer)
	  swc02 = fracd * dzwww(i)*poros*dzmultilayer
	endif
      enddo
          swc01 = swc01 + swcsfc
	  swc02 = swc02 + swc01
 !******************************************************************** 
      ! Saída de umidade para as profundidades fixas da OBS PDG ~RHV
          do i=1, 8
	      out_layers(i) = 1
             do j=1, nlayers
 	      if( niveis(i) .gt. (depth(j)-dzmultilayer) .and.  &
 	          niveis(i) .lt. depth(j) ) then
		    out_layers(i) = j
	      endif
 	    enddo
 	 enddo
 	 
 !********************************************************************
 ! Saída para as umidades volumétricas das camadas:
 !            0.0 - 0.2 m
 !            0.0 - 1.0 m
 !            1.0 - 2.0 m
 ! Casos psrticuares das medidas de PDG, segundo dados fornecidos 
 ! pelo Osvaldo Cabral e JDT
	wcount20 = 0
	wcount01 = 0
	wcount02 = 0
	vwc20 = 0.0 
	vwc01 = 0.0
	vwc02 = 0.0
      do i=1,nlayers
	if( depth(i)-dzmultilayer .lt. 0.2 )then
	   vwc20 = vwc20 + dzwww(i)*poros
	   wcount20 = wcount20 + 1
	endif
	if( depth(i)-dzmultilayer .lt. 1.0 )then
	   vwc01 = vwc01 + dzwww(i)*poros
	   wcount01 = wcount01 + 1
	endif
	if( depth(i) .gt. 1.0 .and. depth(i)-dzmultilayer .lt. 2.0 )then
	   vwc02 = vwc02 + dzwww(i)*poros
	   wcount02 = wcount02 + 1
	endif	
      enddo
      vwc20 = vwc20 / wcount20
      vwc01 = vwc01 / wcount01
      vwc02 = vwc02 / wcount02
      
 !********************************************************************   
    
 
 
 
 
 
! NEW WAY FOR SIB2DIAG.TXT WRITE  ~RHV
      if (iter == 1) then                                           ! write labels
         write(iout1,'(16(1x,a12))') &
	  'istype','ivtype', &                                      !2
	  'zdepth1','zdepth2','zdepth3', 'sodep', &                 !6
          'bee','phsat','satco','poros','slope','vcover', &         !12
          'soilRefVis', 'soilRefIR', &                              !14
          'vmax0', 'gradm'                                          !16
      endif     
      
      if (iter == 1) then                            
      write(iout1,'(2i13,14(1x,f12.7))') &   
         istype, ivtype, &                                          !2
         zdepth(1),zdepth(2),zdepth(3), sodep, &                    !6
         bee, phsat, satco, poros, slope,vcover, &                  !12  
         soref(1), soref(2), &                                      !14
         vmax0, gradm                                               !16
      end if
      
      if (iter == 1) then                            ! write labels
         write(iout1,'(1a12,57(1x,a12))') &
            'nymd', &                                                !0
            'lai', 'fpar', 'green', 'z0d', 'dd',  &                  !5
            'ki','rn','ldwn','lupw','radswd','radswa', &             !11
            'albedo','albedoPAR','em','tm','um','tprec', &           !16
            'LE','eci','ect','egi','egs', &                          !21
            'ec','eg',&                                              !23  
            'H','tc','tg','G','gcstor', &                            !28
            'www1','www2', 'www3','roff','roff1', &                  !33
            'roff2','roff3','gwsoil', 'roffGA', 'roffq3g', &         !38
            'roffq3g2', &                                            !39
            'assimn', 'ra','ustar',         &                        !42
            'swc20','swc01','swc02',                           &    !44
            'dzwww10','dzwww20','dzwww50','dzwww80','dzwww100', &    !49
            'dzwww150','dzwww200','dzwww250', &                      !54
            'vwc20','vwc01','vwc02'                                  !57  
      endif    
      
      write(iout1,'(1(i12), 57(1x,f12.7))') &
         nymd, &                                                      !0
         zlt, fparc, green, z0d, dd, &                                !5
         swdown, radtot, radn(3,2), zlwup, radswd, radswa, &          !11
         calbe,calbe_par, em, tm, um, tprec, &                                  !16
         elat,eci/dtt,ect/dtt,egi/dtt,egs/dtt, &                      !21
         (trant+canil)/ dtt * hlat,(evapg1+evapg2)/ dtt * hlat, &     !23
         hflux,tc,tg, gflux, gcstor, &                                !28
         www(1), www(2), www(3), roff*1000, roff1*1000, &             !33
         roff2*1000, roff3*1000, gwsoil,roffga*1000,roffq3g*1000, &   !38
         roffq3g2 * 1000, &                                           !39
         assimn * 1.e06, ra, ustar,        &                          !42
         swcsfc*1000,swc01*1000,swc02*1000,                    &      !44
         (dzwww(out_layers(i)), i = 1, 8),                     &      !54
         vwc20,vwc01,vwc02                                            !57  
!c.......................................................................
!c! adicionado by tatsch para avaliar inputs horarios nos ptos com H<0
!c Jacutinga vtype 9
!      if (iter == 1) then                            ! write labels
!         write(iout1,'(2(a8), 1(a12), 1(a8), 1(a8), &
!            72(1x,a12))') &
!            'np','iter', 'nymd', &
!            'istype','ivtype', 'zdepth1','zdepth2','zdepth3', 'sodep', &
!            'bee','phsat','satco','poros','slope', 'soilRefVis', 'soilRefIR', &
!            'lai', 'fpar', 'green', 'vcover','z0d', 'dd', 'vmax0', 'gradm', &
!            'ki','rn','ldwn','lupw','radswd','radswa', &
!            'albedo','em','tm','um','tprec', &
!            'LE','eci','ect','egi','egs','ec','eg','H','tc','tg','G','gcstor', &
!            'www1','www2', 'www3','roff','roff1', &
!            'roff2','roff3','gwsoil', 'roffGA', 'roffq3g', 'roffq3g2', &
!            'assimn', 'ra','ustar'
!      endif
!c 52 vars
!c!       print *, swdown, radtot, radn(3,2), zlwup, em, tm
!c!       print*, satco
!c!       pause
 
   
!      write(iout1,'(2(i8), 1(i12), 1(i8), 1(i8), &
!         6(1x,f12.7), 1(1x,es14.7),54(1x,f12.7))') 1, iter, nymd, &   !3
!         istype, ivtype, zdepth(1),zdepth(2),zdepth(3), sodep, &      !9
!         bee, phsat, satco, poros, slope, soref(1), soref(2), &       !16
!         zlt, fparc, green, vcover, z0d, dd, vmax0, gradm, &          !24
!         swdown, radtot, radn(3,2), zlwup, radswd, radswa, &          !30
!         calbe, em, tm, um, tprec, &                                  !35
!         elat,eci/dtt,ect/dtt,egi/dtt,egs/dtt, &                      !37
!         (trant+canil)/ dtt * hlat,(evapg1+evapg2)/ dtt * hlat, &     !38
!         hflux,tc,tg, gflux, gcstor, &                               !41
!         www(1), www(2), www(3), roff*1000, roff1*1000, &             !46
!         roff2*1000, roff3*1000, gwsoil,roffga*1000,roffq3g*1000, &   !51
!         roffq3g2 * 1000, &                                           !53
!         assimn * 1.e06, ra, ustar                                    !55



!   write(200) zdepth(1),zdepth(2),zdepth(3), sodep, bee, phsat, satco, &
!      poros, slope, soref(1), soref(2), zlt, fparc, green, vcover, &
!      z0d, dd, vmax0, gradm, swdown, radtot, radn(3,2), zlwup, radswd, &
!      radswa, calbe, em, tm, um, tprec, elat, (trant+canil)/dtt*hlat, &
!      (evapg1+evapg2)/dtt*hlat, hflux, gflux, gcstor, www(1), www(2), &
!      www(3), roff*1000.0, roff1*1000.0, roff2*1000.0, roff3*1000.0, &
!      gwsoil, roffGA*1000.0, roffq3g*1000.0, roffq3g2*1000.0, &
!      assimn*1.e06, ra, ustar



!!!   endif
!c       pause
!cc!---------
   if (item01 == 25 .and. item02 == 61) then  !item01=row, item02=col
!c Inconfidentes  vtype 6
      write(iout1,'(2(i8), 1(i10), 1(i8), 1(i8), &
         6(1x,f12.7), 1(1x,ES14.7),48(1x,f12.7))') 2, iter, nymd, &   !3
         istype, ivtype, zdepth(1),zdepth(2),zdepth(3), sodep, &      !9
         bee, phsat, satco, poros, slope, soref(1), soref(2), &       !16
         zlt, fparc, green, vcover, z0d, dd, vmax0, gradm, &          !24
         swdown, radtot, radn(3,2), zlwup, radswd, radswa, &          !30
         calbe, em, tm, um, tprec, &                                  !35
         elat, (trant+canil)/ dtt * hlat, &                           !37
         (evapg1+evapg2)/ dtt * hlat, &                               !38
         hflux, gflux, gcstor, &                                      !41
         www(1), www(2), www(3), roff*1000, roff1*1000, &             !46
         roff2*1000, roff3*1000, gwsoil,roffGA*1000,roffq3g*1000, &   !51
         roffq3g2 * 1000, &                                           !53
         assimn * 1.e06, ra, ustar                                    !55
   endif
!c!---------
   if (item01 == 27 .and. item02 == 5) then
!c pto prox Exutorio vtype 1
      write(iout1,'(2(i8), 1(i10), 1(i8), 1(i8), &
         6(1x,f12.7), 1(1x,ES14.7),48(1x,f12.7))') 3, iter, nymd, &   !3
         istype, ivtype, zdepth(1),zdepth(2),zdepth(3), sodep, &      !9
         bee, phsat, satco, poros, slope, soref(1), soref(2), &       !16
         zlt, fparc, green, vcover, z0d, dd, vmax0, gradm, &          !24
         swdown, radtot, radn(3,2), zlwup, radswd, radswa, &          !30
         calbe, em, tm, um, tprec, &                                  !35
         elat, (trant+canil)/ dtt * hlat, &                           !37
         (evapg1+evapg2)/ dtt * hlat, &                               !38
         hflux, gflux, gcstor, &                                      !41
         www(1), www(2), www(3), roff*1000, roff1*1000, &             !46
         roff2*1000, roff3*1000, gwsoil,roffGA*1000,roffq3g*1000, &   !51
         roffq3g2 * 1000, &                                           !53
         assimn * 1.e06, ra, ustar                                    !55
   endif
   
!c.......................................................................
!c
!c      WRITE(iout,900) nymd, www(1), www(2), www(3), ppl
!c900   FORMAT( 1x, i8, 1x, 3f7.4, 1x, f8.2)
!c
!c! Added by Tatsch
!c! output of all ETP components in w m-2
!c      WRITE(38,910) nymd, roff, trant/dtt*hlat, canil/dtt*hlat,
!c     & evapg1/dtt*hlat, evapg2/dtt*hlat, elat, ppl
!c910   FORMAT( 1x, i8, 1x, f12.7, 5f15.4, f8.2)
!c
!c      WRITE(iout1,910) nymd, roff, trant, canil, evapg, etmass, ppl
!c910   FORMAT( 1x, i8, 1x, f12.7, 4f10.4, f8.2)
!c
!c      WRITE(iout2,920) nymd, tc, tg, td, tm, tgeff, capac(1), capac(2)
!c     & ,snoww(1), snoww(2)
!c920   FORMAT( 1x, i8, 1x, 5f7.2, 4f7.5)
!c
!c      WRITE(iout3,930) nymd, radn(3,2), zlwup, radswd, radswa, radtot,
!c     & elat, hflux, gcstor
!c930   FORMAT(1x,i8,1x,8(f5.0,1x) )
!c
!       WRITE(iout4,940) nymd, rst,(rstfac(i),i=1,4),gsh2o,assimn
! 940   FORMAT(1x,i8,1x, f10.1, 1x, 4( f8.6,1x), 1x, f8.4,1x,e10.3 )
end subroutine outer



!c-----------------------------------------------------------------------
!!!	SUBROUTINE 
subroutine vegpar(ivtype)

!!!        use gridij, only : ivtype
   use soilij, only : soref, sodep
   use vstate, only : &
      tran, ref, atheta, btheta, binter, chil, effcon, &
      gradm, hhti, hlti, phc, respcp, rootd, shti, slti, &
      trda, trdm, trop, vcover, z1, z2, rootex
   use vderiv, only : vmax0

   use morphology2, only : z1_v, z2_v, vcover_v, chil_v, rootd_v, sodep_v, &
      rootex_v
   use optical, only : &
      tranlv_v, tranln_v, trandv_v, trandn_v, &
      reflv_v, refln_v, refdv_v, refdn_v, sorefv_v, sorefn_v
   use physiology1, only : btheta_cst, shti_cst, slti_cst, &
      trda_cst, trdm_cst, trop_cst
   use physiology2, only : &
      phc_v, effcon_v, gradm_v, binter_v, respcp_v, atheta_v, &
      hlti_v, hhti_v, vmax0_v

   implicit none

   ! Dummy arguments
   integer, intent(in) :: ivtype

   if (ivtype >= 0) then
      z2 = z2_v(ivtype)
      z1 = z1_v(ivtype)
      vcover = vcover_v(ivtype)
      chil = chil_v(ivtype)

      rootd = rootd_v(ivtype)
      phc = phc_v(ivtype)

      tran(1,1) = tranlv_v(ivtype)
      tran(2,1) = tranln_v(ivtype)
      tran(1,2) = trandv_v(ivtype)
      tran(2,2) = trandn_v(ivtype)

      ref(1,1) = reflv_v(ivtype)
      ref(2,1) = refln_v(ivtype)
      ref(1,2) = refdv_v(ivtype)
      ref(2,2) = refdn_v(ivtype)

      effcon = effcon_v(ivtype)
      gradm = gradm_v(ivtype)
      binter = binter_v(ivtype)
      respcp = respcp_v(ivtype)
      atheta = atheta_v(ivtype)
      btheta = btheta_cst  ! Constant

      trda = trda_cst  ! Constant
      trdm = trdm_cst  ! Constant
      trop = trop_cst  ! Constant
      slti = slti_cst  ! Constant
      hlti = hlti_v(ivtype)
      shti = shti_cst  ! Constant
      hhti = hhti_v(ivtype)

      vmax0 = vmax0_v(ivtype)
      rootex = rootex_v(ivtype)

   else !RUN the sample (ivtype<0)
      z2 = 35.0
      z1 = 1.0
      vcover = 0.98
      chil = 0.1

      rootd = 0.5
      phc = -200.0

      tran(1,1) = 0.0500
      tran(2,1) = 0.2500
      tran(1,2) = 0.0010
      tran(2,2) = 0.0010

      ref(1,1) = 0.1000
      ref(2,1) = 0.4500
      ref(1,2) = 0.1600
      ref(2,2) = 0.3900

      effcon = 0.08
      gradm = 9.00
      binter = 0.01
      respcp = 0.015
      atheta = 0.95
      btheta = 0.95  ! Constant

      trda = 1.30  ! Constant
      trdm = 328.16  ! Constant
      trop = 298.16  ! Constant
      slti = 0.2   ! Constant
      hlti = 288.16
      shti = 0.3   ! Constant
      hhti = 313.16

      vmax0 = 0.00006
   endif

   if (ivtype >= 0) then
      sodep = sodep_v(ivtype)
      soref(1) = sorefv_v(ivtype)
      soref(2) = sorefn_v(ivtype)
   else
      sodep = 2.0
      soref(1) = 0.10
      soref(2) = 0.20
   endif
end subroutine vegpar



!c	SUBROUTINE soipar(soilpara)
!c	INCLUDE 'COMSIBC.H'
!c	include 'SiB2par.inc'
!c	real soilpara(7000,6)
!c
!c	if (istype>0) then
!c		bee	=	soilpara(istype,4)
!c		phsat	=	soilpara(istype,2)
!c		satco	=	soilpara(istype,3)
!c		poros	=	soilpara(istype,1)
!c	else
!c		bee	=	7.797
!c		phsat	=	-0.2
!c		satco	=	3.5E-06
!c		poros	=	0.458
!c	endif
!c
!c	END
!c=======================================================================



!!!	SUBROUTINE varcal(varcal_para)
subroutine varcal(ivtype, aero_parms)

   use caerod, only : corb1, corb2, g1, g2, g3, ha, zmet, ztz0, zwind
!!!        use gridij, only : ivtype
   use soils, only : &
      zdepth, zlay, satcoz, decay, satco, satcoz_a, &
      zlay_a
   use soilij, only : sodep
   use vderiv, only : cc1, cc2, dd, gmudmu, z0d
   use vdyijt, only : green!, zlt
   use vstate, only : rootd, z2, ref, tran

   use morphology1, only : g1_cst, g4_cst

   implicit none

   ! Dummy arguments
   integer, intent(in) :: ivtype

   ! Local variables
   real :: park, scatp

   real, intent(in) :: aero_parms(9)

   if (ivtype >= 0) then
!c	CALL derive_trans
!c     $	(ivtype,zlt,rhoair,
!c     $	ha, z0d, dd,g2, g3, cc1, cc2, corb1, corb2)

      ha      = aero_parms(1)
      z0d     = aero_parms(2)
      dd      = aero_parms(3)
      g2      = aero_parms(4)
      g3      = aero_parms(5)
      cc1     = aero_parms(6)
      cc2     = aero_parms(7)
      corb1   = aero_parms(8)
      corb2   = aero_parms(9)
      g1 = g1_cst
      ztz0 = g4_cst
      
      if (z2 + 10.0 * dd + 0.1 < 2.0) then
         zwind = 2.0   ! (data from Chinese Weather Bureau)
         zmet  = 2.0   ! (usually it is 2.00 m)
!      else
!c! added by tatsch sep/2010, if vtype=1 this give unrealistic values
!c!   zwind = z2+10.*dd+0.1 ! (data from Chinese Weather Bureau)
!c!   zmet  = z2+10.*dd+0.1 ! (usually it is 2.00 m)
!         zwind = 45.00         ! (data from Chinese Weather Bureau)
!         zmet  = 45.00         ! (usually it is 2.00 m)
      endif

   else

      zwind = 45.00 ! (data from Chinese Weather Bureau)
      zmet  = 45.00 ! (usually it is 2.00 m)
      z0d = 2.02
      dd = 28.81
      cc1 = 5.59
      cc2 = 1177.14
      corb1 = 0.111
      corb2 = 19.112
      ha = 24.81
      g1 = 1.449
      g2 = 0.801
      g3 = 0.801
      ztz0 = 11.785

   endif   ! (ivtype>=0)

   rootd = min(rootd, sodep * 0.75)
   zdepth(1) = 0.02
   zdepth(2) = rootd - 0.02
   zdepth(3) = sodep - zdepth(1) - zdepth(2)
!c!----------------------------------------------------------------
!c! tatsch 15 set 2011, added variables to change satco with depth
   zlay(1) =                          0.5 * zdepth(1)  !mid-layer depth
   zlay(2) =             zdepth(1) + (0.5 * zdepth(2))
   zlay(3) = zdepth(1) + zdepth(2) + (0.5 * zdepth(3))

   zlay_a = zlay(1) + zlay(2) + zlay(3)
   zlay_a = zlay_a / 3.0                               !avg mid-layer depth from vadose zone

   satcoz(1) = satco * exp (-decay*zlay(1))            !satco for the 1st layer
   satcoz(2) = satco * exp (-decay*zlay(2))
   satcoz(3) = satco * exp (-decay*zlay(3))

!c!                satcoz_a  = satco * exp (-decay*zlay_a)             !avg satco from vadose zone

   satcoz_a  = satcoz(1)*(zdepth(1)/sodep)+ &
      satcoz(2)*(zdepth(2)/sodep)+ &
      satcoz(3)*(zdepth(3)/sodep)
!c               write(*,'(I8)') ivtype
!c               write(*,'(7f)') 0.02, rootd, sodep
!c               write(*,'(7f)') zdepth(1),zdepth(2), zdepth(3)
!c               write(*,'(7f)') zlay(1),zlay(2), zlay(3), zlay_a,decay
!c                write(*,'(7f)') satco
!c                write(*,'(7f)') satco*1000*3600
!c               write(*,'(4(1x,ES14.4))') satcoz(1),satcoz(2),satcoz(3),
!c     &                          satcoz_a
!c               write(*,'(7f)') satcoz(1)*1000*3600,
!c     &                          satcoz(2)*1000*3600,
!c     &                          satcoz(3)*1000*3600,
!c     &                          satcoz_a*1000*3600
!c                pause

!c!----------------------------------------------------------------
!c! If green and fparc are known, to get zlt(Leaf Area Index) (Original SiB2)
   scatp = green * (tran(1,1) + ref(1,1)) + (1.0 - green) * (tran(1,2) + ref(1,2))
   park = sqrt(1.0 - scatp) * gmudmu

!cc     fparc = 1. - exp ( -park*zlt )
!c  if ((1.-fparc)<=0.0) then
!c   print *,"SIB2sub(359):fprac larger than 1.0:",(1.-fparc)
!c  endif
!c  zlt = -1./park*log( 1.-fparc )

!c! If the zlt(Leaf Area Index) and fparc are known, to get green
!c! To do so, Please remove code before: zlt = -1./park*log(1.-fparc)
!c  scatp = 1. - (log(1.-fparc)/zlt/gmudmu)**2
!c  tranref1= tran(1,1) + ref(1,1)
!c  tranref2= tran(1,2) + ref(1,2)
!c  green = (scatp - tranref2) / ( tranref1 - tranref2 )
!c  green = min(0.95, green)
!c  green = max(0.05, green)
end subroutine varcal



!c=======================================================================

subroutine const2

!c=======================================================================
!c
!c     initialization of physical constants
!c
!c     subroutine const2 is called at every time step
!c     because in some applications constants may depend on
!c     environmental conditions.
!c
!c-----------------------------------------------------------------------


   use constants, only : pi

   use atchem, only : po2m, pco2m
   use atmos, only : bps, psur
   use const, only : &
      asnow, clai, cpair, cw, epsfac, g, kappa, pie, &
      rhoair, snomel, stefan, tf, vkc, rcp, timcon

   implicit none

!!!      INCLUDE 'COMSIBC.H'
!c
!!!      asnow    = 13.2
!!!      bps      = 1.
!!!      clai     = 4.2 * 1000. * 0.2
!!!      cpair    = 1010.
!!!      cw       = 4.2 * 1000. * 1000.
!!!      epsfac   = 0.622
!!!      g        = 9.81
!!!      kappa    = 0.286
!!!CCC      pie      = 3.14159265
!!!      pie      = pi
!!!      po2m     = 20900.
!!!      pco2m    = 34.
!!!      psur     = 1000.
!!!      rhoair   = 1.225
!!!      snomel   = 370518.5 * 1000.
!!!      stefan   = 5.669 * 10e-9
!!!      tf       = 273.16
!!!      vkc      = 0.41
!!!      rcp      = rhoair * cpair
!!!      timcon   = pie / 86400.
!c
!c-----------------------------------------------------------------------
!c     n.b. :  snomel is expressed in j m-1
!c-----------------------------------------------------------------------
end subroutine const2



!c=======================================================================
!c
!!!      SUBROUTINE inter2
subroutine inter2(idirr)

!c=======================================================================
!c
!c     calculation of  interception and drainage of rainfall and snow
!c     incorporating effects of patchy snow cover and temperature
!c     adjustments.
!c
!c----------------------------------------------------------------------
!c
!c     (1) non-uniform precipitation
!c         convective ppn. is described by area-intensity
!c         relationship :-
!c
!c                   f(x) = a*exp(-b*x)+c
!c
!c         throughfall, interception and infiltration
!c         excess are functional on this relationship
!c         and proportion of large-scale ppn.
!c         reference: sato et al.(1989b), appendix.
!c
!c     (2) reorganisation of snowmelt and rain-freeze procedures.
!c               subroutine adjust
!c
!c     (3) additional calculation for partial snow-cover case.
!c               subroutine patchs
!c
!c     (4) reorganisation of overland flow.
!c         reference: SA-89B, appendix.
!c
!c     (5) modified calaculation of soil heat capacity and
!c         conductivity.
!c
!c=======================================================================
!c
!c     subroutines in this block : snow1
!c     -------------------------   adjust
!c   patchs
!c   snow1
!c
!c++++++++++++++++++++++++++++++output+++++++++++++++++++++++++++++++++++
!c
!c       roff           runoff (mm)
!c       tc             canopy temperature (K)
!c       tg             ground surface temperature (K)
!c       www(1)         ground wetness of surface layer
!c       capac(2)       canopy/ground liquid interception store (m)
!c       snoww(2)       canopy/ground snow interception store (m)
!c
!c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   use atmos, only : ppc, ppl, tm
   use const, only : clai, cw, pie, tf
   use donor, only : roff
   use grads, only : tgs
   use hydrol, only : areas, satcap
!!!        use irrgrid, only : idirr
   use output, only : &
      gwsoil, roff1, roff2, roff3, roff4, &
      roffga, roffq3g, roffq3g2
   use preccoff, only : app, bpp, cpp
   use radabs, only : exrain
   use snow, only : tsnow
   use soils, only : poros, zdepth, satcoz
   use suroff, only : finfil
   use steps, only : dtt!,iter
   use stepv, only : capac, snoww, www, tc, tg
   use stores, only : csoil
   use vstate, only : chil, vcover
   use vdyijt, only : zlt
   use multicamada
   implicit none

   ! Dummy arguments
   integer, intent(in) :: idirr

   ! Local variables
   integer :: iveg!,i
   real :: aa, bb, ap, arg, bp, capacp, chiv, cp, equdep, fpi
   real :: p0, pinf, props, pxi, realc, realg, roffo
   real :: shcap, slamda, snowwp, spechc, tex, thru, totalp, ts
   real :: tti, xs, xsc, xss, xw1, zload

   real :: pcoefs(2,2)

!c      data pcoefs(1,1)/ 20. /, pcoefs(1,2)/ .206e-8 /,
!c     &     pcoefs(2,1)/ 0.0001 /, pcoefs(2,2)/ 0.9999 /, bp /20. /
!c!      data pcoefs(1,1)/ 25. /, pcoefs(1,2)/ 1.38879E-11 /,
!c!     &     pcoefs(2,1)/ 25. /, pcoefs(2,2)/ 1.38879E-11 /, bp /25. /
!c! tatsch
   data pcoefs(1,1)/ 0.0001 /, pcoefs(1,2)/ 0.9999 /, &
      pcoefs(2,1)/ 0.0001 /, pcoefs(2,2)/ 0.9999 /, bp /20. /

!c-----------------------------------------------------------------------

   call snow1

!c-----------------------------------------------------------------------
!c
!c     prec ( pi-x )   : equation (c.3), SA-89B
!c
!c-----------------------------------------------------------------------

	pcoefs(1,1) = app
	pcoefs(2,1) = app
	pcoefs(1,2) = cpp
	pcoefs(2,2) = cpp
	bp = bpp

   ap = pcoefs(2,1)
   cp = pcoefs(2,2)
   totalp = ppc + ppl
   if( snoww(1) > 0. .or. snoww(2) > 0. .or. tm < tf ) &
      ppc = 0.
   ppl = totalp - ppc
   if(totalp >= 1.e-8) then
      ap = ppc/totalp * pcoefs(1,1) + ppl/totalp * pcoefs(2,1)
      cp = ppc/totalp * pcoefs(1,2) + ppl/totalp * pcoefs(2,2)
   end if

   roff = 0.0
   roff1 = 0.0
   roff2 = 0.0
   roffga  = 0.0     ! tatsch
   roff3	= 0.0
   roffq3g = 0.0
   roffq3g2 = 0.0
   roff4 = 0.0
   gwsoil = 0.0
   thru = 0.0
   fpi = 0.0
   finfil = 0.0

!c----------------------------------------------------------------------
!c     heat capacity of the soil, as used in force-restore heat flux
!c     description. dependence of csoil on porosity and wetness is
!c     based on CS-81.
!c----------------------------------------------------------------------

   slamda = (1.5 * (1.0 - poros) + 1.3 * www(1) * poros) &
      / (0.75 + 0.65 * poros - 0.4 * www(1) * poros) * 0.4186
   shcap = (0.5 * (1.0 - poros) + www(1) * poros) * 4.186 * 1.e6
   csoil = sqrt(slamda * shcap * 86400.0 / pie) / 2.0

!c----------------------------------------------------------------------
!c     input precipitation is given in mm, converted to m to give p0.
!c----------------------------------------------------------------------

   p0 = totalp * 0.001

   do iveg = 1, 2

      realc = 2.0 - iveg
      realg = iveg - 1.0

      xsc = max(0.0, capac(iveg) - satcap(iveg))
      capac(iveg) = capac(iveg) - xsc
      xss = max(0.0, snoww(iveg) - satcap(iveg)) * realc
      snoww(iveg) = snoww(iveg) - xss
      if (iveg == 1) then  !Because freezing, snoww(1)>satcap, put to snoww(2)
         snoww(2) =snoww(2) + xss
         xss = 0.0
      endif
      p0 = p0 + xsc + xss			!tang@2005/12/21

      capacp = capac(iveg)
      snowwp = snoww(iveg)

      spechc = min(0.05, capac(iveg) + snoww(iveg)) * cw &
         + realc * zlt * clai + realg * csoil
      ts = tc * realc + tg * realg

!c----------------------------------------------------------------------
!c     proportional saturated area (xs) and leaf drainage(tex)
!c
!c     tex ( d-c )     : equation (c.8), SA-89B
!c     xs  ( x-s )     : equation (c.7), SA-89B
!c     tex ( d-c )     : equation (c.8), SA-89B
!c
!c-----------------------------------------------------------------------

      chiv = chil
      if (abs(chiv) <= 0.01) chiv = 0.01
      aa = 0.5 - 0.633 * chiv - 0.33 * chiv * chiv
      bb = 0.877 * (1.0 - 2.0 * aa)
      exrain = aa + bb

      zload = capac(iveg) + snoww(iveg)
      fpi = (1.0 - exp(-exrain * zlt / vcover)) * vcover * realc + realg
      tti = p0 * (1.0 - fpi)
      xs = 1.0
      if (p0 >= 1.e-9) then
         arg = (satcap(iveg) - zload) / (p0 * fpi * ap) - cp / ap
         if ( arg >= 1.e-9 ) then
            xs = -1.0 / bp * log(arg)
            xs = min(xs, 1.0)
            xs = max(xs, 0.0)
         end if
      end if
      pxi = max(0.0, satcap(iveg) - zload)
      tex = p0 * fpi * (ap / bp * (1.0 - exp(-bp * xs)) + cp * xs) &
           - pxi * xs
      tex = max(tex, 0.0)

!c----------------------------------------------------------------------
!c     total throughfall (thru) and store augmentation
!c----------------------------------------------------------------------

      if ( iveg == 1 ) then

         thru = tti + tex
         pinf = p0 - thru
         if (tm > tf) capac(iveg) = capac(iveg) + pinf
         if (tm <= tf) snoww(iveg) = snoww(iveg) + pinf

         call adjust(tc, spechc, capacp, snowwp, iveg)

         p0 = thru

      else if (tg > tf .and. snoww(2) > 0.0) then

         call patchs(p0)

      else

         thru = tti + tex
         if (tg <= tf .or. tm <= tf) thru = 0.0
         pinf = p0 - thru
         if (tm > tf) capac(iveg) = capac(iveg) + pinf
         if (tm <= tf) snoww(iveg) = snoww(iveg) + pinf
         if (tm > tf) then
!c            www(1) = www(1) + thru/ ( poros*zdepth(1) )
            finfil = finfil + thru
         end if

         call adjust ( tg, spechc, capacp, snowwp, iveg )

      end if
   end do

!c----------------------------------------------------------------------
!c
!c     instantaneous overland flow contribution ( roff )
!c     roff( r-i )     : equation (c.13), SA-89B
!c     finfil : through fall and snow melt etc which reach the soil top surface
!c-----------------------------------------------------------------------

   xw1 = max(0.0, www(1) - 1.0)
   www(1) = www(1) - xw1
   finfil = finfil + xw1 * (poros * zdepth(1))
!c----------------------------------------------------------------------
!c if froze, reduce infiltration
!c-----------------------------------------------------------------------
   tsnow = min(tf - 0.01, tg)
   tgs = tsnow * areas + tg * (1.0 - areas)
   props = (tgs - (tf - 1.0)) / 1.0
   props = max(0.0, min(1.0, props))

   if (idirr /= 1) then
      ap = app
      bp = bpp
      cp = cpp
!c! equdep = satco * dtt * props
      equdep = satcoz(1) * dtt * props
      xs = 1.0
      if( finfil >= 1.e-9 ) then
         arg = equdep / (finfil * ap) - cp / ap
         if (arg >= 1.e-9) then
            xs = -1.0 / bp * log(arg)
            xs = min(xs, 1.0)
            xs = max(xs, 0.0)
         end if
      end if
      roffo = finfil * (ap / bp * (1.0 - exp(-bp * xs)) + cp * xs) - equdep * xs
      roffo = max(roffo, 0.0)
      roff = roff + roffo
      roff1 = roff1 + roffo
      finfil = finfil - roffo
   endif
!c if (finfil>0.) then
!c  xw1 = min(finfil, (1.0-www(1) )* ( poros*zdepth(1) ) )
!c  www(1) = www(1) + xw1 / ( poros*zdepth(1) )
!c  finfil = finfil -xw1
!c endif
!     write(543,*)"antes",(www(i),i=1,3)
   call green_ampt(finfil)
!c thru_re= finfil
   roff = roff + finfil
   roff4 = roff4 + finfil
   finfil = 0.0
   
!   if(iter .lt. 10)  write(*,'(a8,i10, 50(f13.6))')'GREEN-A ',iter,www(1), &
!     (dzwww(i),i=1,nlayers)  , roff3           ! WRITING
end subroutine inter2



!c=======================================================================
subroutine green_ampt(dep)

   use donor, only : roff
   use output, only : roff3, roffga
   use soils, only : bee, phsat, poros, satcoz_a, zdepth
   use steps, only : dtt, iter
   use stepv, only : www
   use multicamada 
   use morphology1, only : dsfc_cst
   implicit none

   ! Dummy arguments
   real :: dep

   ! Local variables
   real :: dt, dwww, f_i, f_i2, ph_a, put, seli, totwww, xw, y, y0
   real :: se, a1, b1
     real :: xw0,depwww,fracdz3,drenagem,maxW  !locais RHV
     integer :: i                !locais RHV
     real :: pwww2,pwww3, count2, count3
!c  Calculate infiltration to soil and instantaneous overland flow (roff2)
    xw0 = 0.0
    y = 0.0
    maxW = 1.0
    
   if (dep > 0.0) then
    maxW = min(1.0, (dep / (satcoz_a*dtt))**(1.0 / (2.0 * bee + 3.0)) )
    maxW = max(maxW, 0.0)
   
      seli = 1.0
      se = (www(1) * zdepth(1) + www(2) * zdepth(2) + www(3) * zdepth(3)) &
         / (zdepth(1) + zdepth(2) +zdepth(3))
!c  se = www(1)
      se = max(0.03, se)
      se = min(se, 1.0)
      ph_a = -phsat * (se**(-bee) - 1.0) / (zdepth(1) + zdepth(2) + zdepth(3)) * 2.0
!c!  f_i = satco*(1.+ ph_a )
      f_i = satcoz_a*(1.0 + ph_a)
      if (dep > f_i) then   !case 1
	  maxW = 1.0
         a1 = -phsat * max(0.0, seli - se) * poros  !*se**(-bee)
!c!   b1 = dtt*satco
         b1 = dtt * satcoz_a
         y0 = a1 + b1
         call newton_i(a1, b1, y0, y)
      else
         y = dep
         a1 = -phsat * max(0.0, seli - se) * poros  !*se**(-bee)
!c!   f_i2= satco*( 1 + a1/y )
         f_i2 = satcoz_a * (1.0 + a1 / y)
         
         if (dep <= f_i2 * dtt) then !case 3
            y = dep
         else
!c!    dt= satco*a1/(f_i2-satco)/f_i2
            dt = satcoz_a * a1 / (f_i2 - satcoz_a) / f_i2
            a1 = -phsat * max(0.0, seli - se) * poros  !*se**(-bee)
!c!    b1 = (dtt-dt)*satco
            b1 = (dtt - dt) * satcoz_a
            call newton_i(a1, b1, y0, y)
         endif
      endif

!c  a1 = -phsat*max(0.,seli-se)*poros !*se**(-bee)
!c  b1 = dtt*satco
!c  y0 = a1 + b1
!c  call newton_i(a1,b1,y0,y)
!c  write(*,'(7f)') phsat,poros,se,satco,a1,b1,y

!!!     MODIFICATIONS - ROYLAN
!!!        WRITE(1049,'(6E20.8)') PHSAT, SATCOZ_A, SE, PH_A, DEP, Y

      y = max(0.0, y)
      put = min(dep, y)
      totwww = (www(1) * zdepth(1) + www(2) * zdepth(2) &
         + www(3) * zdepth(3)) * poros
      www(1) = www(1) + put / (poros * zdepth(1))
      xw = max(0.0, www(1) - maxW)
      www(1) = www(1) - xw
      xw0 = xw
   else !(dep > 0.0)
      put = 0.0
      totwww = ( www(1) * zdepth(1) + www(2) * zdepth(2) &
         + www(3) * zdepth(3) ) * poros         
   end if !(dep < 0.0)
   
   !   www(2) = www(2) + xw * zdepth(1) / zdepth(2)
   !   xw = max(0.0, www(2) - 1.0)
   !   www(2) = www(2) - xw
   !   www(3) = www(3) + xw * zdepth(2) / zdepth(3)
   !   xw = max(0.0, www(3) - 1.0)
   !   www(3) = www(3) - xw
   !   dwww = (www(1) * zdepth(1) + www(2) * zdepth(2) &
   !      + www(3) * zdepth(3)) * poros - totwww
   !   roff = roff + put - dwww
!c! variavel roffGA foi criada para guardar o roff gerado por GAmpt
   !   roffGA = roffGA + put - dwww
!c!                roff2  = roff2 + put - dwww  !(Experimento 2)
   !   roff3 = roff3 + put - dwww ! original
   !   dep = dep - put
   ! endif
   
   if( iter .eq. 1 )then
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Preparando Green_Ampt para modelo multicamadas. 
      ! Depois de preencher a 1ra camada o restante de infiltração 
      ! preenche sucessivamente as subcamadas definidas.
      depwww = dsfc_cst ! a profundidade de la primeira camada eh constante
       ! Limitante para o preenchimento com a infiltração***
      !igual a 0.02m e comeca a 2da camada.
      !dzmultilayer = (zdepth(2)+zdepth(3))/(nlayers+0.0) ! Sustituir aqui pelo fator de 
      ! discretizacao (nlayers)
      xw0 = xw0 * zdepth(1)/dzmultilayer
      do i=1,nlayers ! Substituir aqui pelo (nlayers - 1)
	depwww = depwww + dzmultilayer ! Profundidade das camadas...
	!... consideradas no esquema de infiltracao de GA, 
	! soma acumuluda do dz.
	!>    
	if(depwww .le. (zdepth(1)+zdepth(2)) ) then
	! Definendo para as subcamadas que originalmente faziam parte da 
	! 2da camada...
	  dzwww(i) = www(2) + xw0 ! * dzmultilayer/dzmultilayer = 1 
	  xw0 = max(0.0,dzwww(i)-maxW)
	  dzwww(i)=dzwww(i)-xw0
	else if( depwww .gt. (zdepth(1)+zdepth(2)) &
	         .and. &
	         depwww .lt. (zdepth(1)+zdepth(2)+dzmultilayer) ) then
	! Fracao de dzcamada que tem parte dela na 2da camada e parte na 3ra 
	  fracdz3 = (depwww - (zdepth(1)+zdepth(2)))/dzmultilayer
	  dzwww(i) = ( (fracdz3 * www(3)) + ((1-fracdz3)*www(2)) ) + xw0
	  xw0 = max(0.0,dzwww(i)-maxW)
	  dzwww(i) = dzwww(i) - xw0
	else
	! Definindo para as camadas que originalmente faziam parte da...
	! ... 3ra camada   
	  dzwww(i) = www(3) + xw0 ! * dzmultilayer/dzmultilayer = 1 
	  xw0 = max(0.0,dzwww(i)-maxW)
	  dzwww(i)=dzwww(i)-xw0
	endif
	! Necessário fazer todo o ciclo para que as camadas sejam atualizadas 
	! com o valor das umidades do solo no tempo de integracao, 
	! Se ocorre infiltracao a saturacao acontece ate ser distribuída totalmente
      enddo
      end if ! (iter == 1 )
      
      if( xw0 .gt. 0.0)then
! 	    depwww = dsfc_cst
	    xw0 = xw0 * zdepth(1)/dzmultilayer
	  do i=1,nlayers ! Substituir aqui pelo (nlayers - 1)
	      dzwww(i) = dzwww(i) + xw0 ! * dzmultilayer/dzmultilayer = 1 
	      xw0 = max(0.0,dzwww(i)-max(maxW,dzwww(i)-xw0 ))
	      dzwww(i)=dzwww(i)-xw0
	  enddo
      endif
      
      dwww = (www(1) * zdepth(1) + sum(dzwww(1:nlayers)) * dzmultilayer ) * poros - totwww
      if(dwww .lt. 1.0e-12) dwww = 0.0
      drenagem = put - dwww
      if(drenagem .lt. 1.0e-8) drenagem = 0.0
      roff = roff + drenagem
      roffGA = roffGA + drenagem
      roff3 = roff3 + drenagem ! original
      dep = dep - put
      
            ! Retornando a modelo de três camadas para conferir resultados
       DEPWWW = DSFC_CST ! RETORNANDO A PROFUNDIDADE DA PRIMEIRA CAMADA
       COUNT2 = 0.0
       COUNT3 = 0.0
       PWWW2  = 0.0
       PWWW3  = 0.0
 	  ! REAL :: PWWW2,PWWW3, COUNT2, COUNT3
       DO I=1,NLAYERS
 	DEPWWW = DEPWWW + DZMULTILAYER ! PROFUNDIDADE DA BASE DA SUBCAMADA ...!
 	  IF (DEPWWW .LE. (ZDEPTH(1) + ZDEPTH(2))) THEN
 	    COUNT2 = COUNT2 + 1.0
	    ! INCREMENTO NO GRAU DE SATURACAO DA 3RA CAMADA
 	    PWWW2 = PWWW2 + DZWWW(I)
 	  ELSE IF (DEPWWW .GT. (ZDEPTH(1) + ZDEPTH(2))  &
 	                  .AND.                         &
 	           DEPWWW .LT. (ZDEPTH(1) + ZDEPTH(2) + DZMULTILAYER)) THEN
 	      ! FRACAO DE DZCAMADA QUE TEM PARTE DELA NA 2DA CAMADA E PARTE NA 3RA 
 	      FRACDZ3 = (DEPWWW - (ZDEPTH(1)+ZDEPTH(2)))/DZMULTILAYER
 	      ! INCREMENTO AO NUMERO DE SUBCAMADAS DA SEGUNDA
	      COUNT2 = COUNT2 + (1 - FRACDZ3)
 	      PWWW2  = PWWW2  + (1 - FRACDZ3) * DZWWW(I)
 	      ! INCREMENTO AO NUMERO DE SUBCAMADAS NA 3RA CAMADA
 	      COUNT3 = COUNT3 + FRACDZ3
 	      ! INCREMENTO NO GRAU DE SATURACAO DA 3RA CAMADA
 	      PWWW3  = PWWW3  + FRACDZ3 * DZWWW(I)
 	  ELSE
 	    ! CONTADOR DAS CAMADAS PARA DETERMINAR O GRAU DE SATURACAO MÉDIO 
 	    COUNT3 = COUNT3 + 1.0
 	    ! INCREMENTO NO GRAU DE SATURACAO DA 3RA CAMADA
 	    PWWW3  = PWWW3  + DZWWW(I)
 	  END IF
       ENDDO
       ! DETERMINANDO GRAU DE SATURAÇÃO MÉDIO DE CADA CAMADA
        WWW(2) = PWWW2 / COUNT2
        WWW(3) = PWWW3 / COUNT3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
end subroutine green_ampt



!c=======================================================================
subroutine newton_i(a1, b1, y0, y)
   implicit none

   real :: fyn, fyn1

   real :: a1, b1, y0, y, cret, err
   integer :: nn

   cret = 1.0e-10
   nn = 0
   err = cret + 1.0
   y = y0
   do while (nn<=100.and.abs(err)>cret)
      fyn = a1 * log(1.0 + y / a1) - y + b1
      err = fyn
      fyn1 = 1.0 / (1.0 + y / a1) - 1.0
      y = y - fyn / fyn1
      nn = nn + 1
   end do
end subroutine newton_i



!c=======================================================================

subroutine adjust(ts, spechc, capacp, snowwp, iveg)

!c=======================================================================
!c
!c     temperature change due to addition of precipitation
!c
!c++++++++++++++++++++++++++++++output+++++++++++++++++++++++++++++++++++
!c
!c     roff           runoff (mm)
!c     tc             canopy temperature (K)
!c     tg             ground surface temperature (K)
!c     www(1)         ground wetness of surface layer
!c     capac(2)       canopy/ground liquid interception store (m)
!c     snoww(2)       canopy/ground snow interception store (m)
!c
!c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   use atmos, only : tm
   use const, only : cw, snomel, tf
   use hydrol, only : satcap
   use stepv, only : snoww, capac
   use suroff, only : finfil

   IMPLICIT NONE

   ! Dummy arguments
   integer, intent(in) :: iveg
   real, intent(in) :: spechc, capacp, snowwp
   real :: ts

   ! Local variables
   real :: cca, ccb, ccp, cct, diff, freeze, tsd, tta, ttb, xs, ccc

   freeze = 0.0
   diff = (capac(iveg) + snoww(iveg) - capacp - snowwp) * cw
   ccp = spechc
   cct = spechc + diff

   tsd = (ts * ccp + tm * diff) / cct

   if ((ts - tf) * (tm - tf) < 0.0) then

      tta = ts
      ttb = tm
      cca = ccp
      ccb = diff
      if (tsd <= tf) then

!c----------------------------------------------------------------------
!c    freezing of water on canopy or ground
!c----------------------------------------------------------------------

         if (ts >= tm) then
            ccc = capacp * snomel
         else
            ccc = diff * snomel / cw
         end if

         tsd = (tta * cca + ttb * ccb + ccc) / cct

         freeze = (tf * cct - (tta * cca + ttb * ccb))
         freeze = min(ccc, freeze) / snomel
         if (tsd > tf) tsd = tf - 0.01

      else

!c----------------------------------------------------------------------
!c    melting of snow on canopy or ground, water infiltrates.
!c----------------------------------------------------------------------

         if (ts <= tm) then
            ccc = - snoww(iveg) * snomel
         else
            ccc = -diff * snomel / cw
         end if

         tsd = (tta * cca + ttb * ccb + ccc) / cct

         freeze = tf * cct - (tta * cca + ttb * ccb)
         freeze = max(ccc, freeze) / snomel
         if (tsd <= tf) tsd = tf - 0.01

      end if
   end if

   snoww(iveg) = snoww(iveg) + freeze
   capac(iveg) = capac(iveg) - freeze

   xs = max(0.0, capac(iveg) - satcap(iveg))
   if (snoww(iveg) >= 0.0000001) xs = capac(iveg)
!c      www(1) = www(1) + xs / ( poros * zdepth(1) )
   finfil = finfil + xs
   capac(iveg) = capac(iveg) - xs
   ts = tsd
end subroutine adjust



!c=======================================================================

subroutine patchs(p0)

!c=======================================================================
!c
!c     marginal situation: snow exists in patches at temperature tf
!c     with remaining area at temperature tg > tf.
!c
!c----------------------------------------------------------------------
!c
!c     calculation of effect of intercepted snow and rainfall on ground.
!c     patchy snowcover situation involves complex treatment to keep
!c     energy conserved.
!c
!c----------------------------------------------------------------------

   use atmos, only : tm
   use const, only : asnow, cw, snomel, tf
   use hydrol, only : areas
   use stepv, only : snoww, capac, tg
   use stores, only : csoil
   use suroff, only : finfil

   implicit none

   ! Dummy arguments
   real, intent(in) :: p0

   ! Local variables
   real :: dareas, dcap, ex, pinf, rhs, snowhc, thru, tsd, zmelt

   pinf = p0
   thru = 0.0
   snowhc = min(0.05, snoww(2)) * cw
   areas = min(1.0, asnow * snoww(2))

   if (tm <= tf) then

      write(*,*) 'entrando em patchs (sib2sub.f), com tm < tf'
      zmelt = 0.0

!c----------------------------------------------------------------------
!c     snow falling onto area
!c----------------------------------------------------------------------

      rhs = tm * pinf * cw + tf * (snowhc + csoil * areas) &
         + tg * csoil * (1.0 - areas)
      dareas = min(asnow * pinf, 1.0 - areas)
      ex = rhs - tf * pinf * cw - tf * (snowhc + csoil * (areas + dareas)) &
         - tg * csoil * (1.0 - areas - dareas)
      if (areas + dareas >= 0.999) tg = tf - 0.01

      if (ex >= 0.0) then

!c----------------------------------------------------------------------
!c     excess energy is positive, some snow melts and infiltrates.
!c----------------------------------------------------------------------

         if (asnow * (snoww(2) + pinf - zmelt) > 1.0) then
            zmelt = ex / snomel
         else
            if (asnow * (snoww(2) + pinf) < 1.0) then
               zmelt = 0.0
            else
               zmelt = (asnow * (snoww(2) + pinf) - 1.0) / asnow
            end if
            zmelt = (ex - zmelt * snomel) &
               / (snomel + asnow * csoil * (tg - tf)) + zmelt
         end if
         snoww(2) = snoww(2) + pinf - zmelt
!c          www(1) = www(1) + zmelt/(poros*zdepth(1))
         finfil = finfil + zmelt

      else

!c----------------------------------------------------------------------
!c     excess energy is negative, bare ground cools to tf, then whole
!c     area cools together to lower temperature.
!c----------------------------------------------------------------------

         if (areas + dareas > 0.999) then
            tsd = 0.0
         else
            tsd = ex/(csoil*( 1.-areas-dareas)) + tg
         end if
         if( tsd <= tf ) then
            tsd = tf + ( ex - (tf-tg)*csoil*(1.-areas-dareas)) &
               / (snowhc+pinf*cw+csoil)
         end if
         tg = tsd
         snoww(2) = snoww(2) + pinf
      end if   ! (ex>=0.)

   else

!c----------------------------------------------------------------------
!c     rain falling onto area
!c----------------------------------------------------------------------
!c
!c
!c----------------------------------------------------------------------
!c     rain falls onto snow-free sector first.
!c----------------------------------------------------------------------

      if (areas >= 0.999) then
         tsd = tf - 0.01
      else
         tsd = (tm * pinf * cw + tg * csoil) / (pinf * cw + csoil)
      end if
      tg = tsd
!c        www(1)= www(1)+pinf*(1.-areas)/(poros*zdepth(1))
      finfil = finfil + pinf * (1.0 - areas)

!c----------------------------------------------------------------------
!c     rain falls onto snow-covered sector next.
!c----------------------------------------------------------------------

      ex = (tm - tf) * pinf * cw * areas
      dcap = -ex / (snomel + (tg - tf) * csoil * asnow)
      if (snoww(2) + dcap >= 0.0) then
!c          www(1) = www(1)+(pinf*areas-dcap)/(poros*zdepth(1))
			finfil = finfil + (pinf * areas - dcap)
         snoww(2) = snoww(2) + dcap
      else
         tg = (ex - snomel * snoww(2) - (tg - tf) * csoil * areas) & 
            / csoil + tg
!c          www(1)=www(1)+(snoww(2)+pinf*areas)/(poros*zdepth(1))
         finfil = finfil + (snoww(2) + pinf * areas)
         capac(2) = 0.0
         snoww(2) = 0.0
      end if
   end if
end subroutine patchs



!c=======================================================================

subroutine snow1

!c=======================================================================
!c
!c     calculation of effects of snow cover on surface morphology and
!c     maximum water storage values.
!c
!c++++++++++++++++++++++++++++++output+++++++++++++++++++++++++++++++++++
!c
!c       z0             roughness length (m)
!c       d              zero plane displacement (m)
!c       rbc            rb coefficient (c1) (s m-1)**1/2
!c       rdc            rd coefficient (c2)
!c       satcap(2)      interception capacities (m)
!c       canex          fraction of exposed canopy (snow-free)
!c       areas          ground snow cover fraction
!c
!c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!c
!c    alteration of aerodynamic transfer properties in case of snow
!c    accumulation. calculation of maximum water storage values.
!c
!c      canex       (fraction of canopy not covered by snow)
!c      d           (snow-modified value of dd, used in all calculations)
!c      z0          (snow-modified value of z0d, used in all calculations)
!c      rbc         (snow-modified value of cc1, used in all calculations)
!c      rdc         (snow-modified value of cc2, used in all calculations)
!c      areas       (fraction of ground covered by snow)
!c      satcap(1)   (s-c)   : equation (56) , SE-86, page 741 se-89
!c      satcap(2)   (s-g)   : 0.002, surface interception store
!c----------------------------------------------------------------------

   use const, only : asnow
   use hydrol, only : areas, canex, satcap,satcap_c,satcap_g
   use rause
   use stepv, only : snoww
   use vderiv, only : z0d, dd, cc1, cc2
   use vdyijt, only : zlt
   use vstate, only : z1, z2

   implicit none

   canex = 1.0 - (snoww(2) * 5.0 - z1) / (z2 - z1)
   canex = max(0.1, canex)
   canex = min(1.0, canex)
   d = z2 - (z2 - dd) * canex
   z0 = z0d / (z2 - dd) * (z2 - d)
   rbc = cc1 / canex
   rdc = cc2 * canex
   areas = min(1.0, asnow * snoww(2))
!   satcap(1) = zlt * 0.0001 * canex
    satcap(1) = zlt * satcap_c * canex
!c 3/96 changes
!c old      satcap(2) = 0.002
!   satcap(2) = 0.0002
   satcap(2) = satcap_g
end subroutine snow1




!c=======================================================================

subroutine rada2

!c=======================================================================
!c
!c     calculation of albedos via two stream approximation( direct
!c     and diffuse ) and partition of radiant energy
!c
!c-----------------------------------------------------------------------
!c
!c     subroutines  called  : snow1
!c     --------------------   longrn
!c
!c
!c++++++++++++++++++++++++++++++output++++++++++++++++++++++++++++++++
!c
!c       salb(2,2)      surface albedos
!c       tgeff          effective surface radiative temperature (k)
!c       radfac(2,2,2)  radiation absorption factors
!c       thermk         canopy gap fraction for tir radiation
!c
!c++++++++++++++++++++++++++diagnostics+++++++++++++++++++++++++++++++
!c
!c       albedo(2,2,2)  component reflectances
!c       closs          tir emission from the canopy (w m-2)
!c       gloss          tir emission from the ground (w m-2)
!c
!rhv-NOTA-  albedo(canopy(1)/ground(2), vis(1)/infraverm(2), direct(1)/diffuse(2))
!c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   use atmos, only : sunang, radn
   use const, only : stefan, tf
   use donor, only : zlwup
   use grads, only : tgs
   use hydrol, only : areas, canex, satcap
   use radabs, only : albedo, radfac, radt, thermk
   use site, only : salb
   use soilij, only : soref
   use stepv, only : tc, tg, snoww
   use vdyijt, only : green, zlt
   use vstate, only : chil, vcover, ref, tran

   implicit none

   integer :: irad, iveg, iwave
   real :: aa, acss, bb, be, betao, bot, ce, chiv, closs, de, den
   real :: ek, epsi, extkb, f, f1, fac1, fac2, facs, fe, fmelt
   real :: ge, gloss, hh1, hh10, hh2, hh3, hh4, hh5, hh6, hh7, hh8, hh9
   real :: power1, power2, proj, psi, reff1, reff2, scat, scov
   real :: tc4, tg4, tran1, tran2, upscat, zat, zkat, zmew, zmk, zp
   real :: tranc1(2), tranc2(2), tranc3(2)

   f = sunang

!c----------------------------------------------------------------------
!c
!c
!c     modification for effect of snow on upper story albedo
!c         snow reflectance   = 0.80, 0.40 . multiply by 0.6 if melting
!c         snow transmittance = 0.20, 0.54
!c
!c
!c-----------------------------------------------------------------------



   call snow1

   facs = (tg - tf) * 0.04
   facs = max(0.0, facs)
   facs = min(0.4, facs)
   fmelt = 1.0 - facs

!   do 1000 iwave = 1, 2
   do iwave = 1, 2
  
      scov = min(0.5, snoww(1) / satcap(1))
      reff1 = (1.0 - scov) * ref(iwave,1) + scov * (1.2 - iwave * 0.4) * fmelt
      reff2 = (1.0 - scov) * ref(iwave,2) + scov * (1.2 - iwave * 0.4) * fmelt
      tran1 = tran(iwave,1) * (1.0 - scov) + scov &
         * (1.0 - (1.2 - iwave * 0.4) * fmelt) * tran(iwave,1)
      tran2 = tran(iwave,2) * (1.0 - scov) + scov &
         * (1.0 - (1.2 - iwave * 0.4) * fmelt) * tran(iwave,2) * 0.9

!c-----------------------------------------------------------------------
!c
!c     calculate average scattering coefficient, leaf projection and
!c     other coefficients for two-stream model.
!c
!c      scat  (omega)        : equation (1,2) , SE-85
!c      proj  (g(mu))        : equation (13)  , SE-85
!c      extkb (k, g(mu)/mu)  : equation (1,2) , SE-85
!c      zmew  (int(mu/g(mu)) : equation (1,2) , SE-85
!c      acss  (a-s(mu))      : equation (5)   , SE-85
!c      extk  (k, various)   : equation (13)  , SE-85
!c      upscat(omega-beta)   : equation (3)   , SE-85
!c      betao (beta-0)       : equation (4)   , SE-85
!c      psi   (h)            : appendix       , SE-85
!c
!c-----------------------------------------------------------------------

      scat = green * (tran1 + reff1) + (1.0 - green) * (tran2 + reff2)
      chiv = chil

      if (abs(chiv) <= 0.01) chiv = 0.01
      aa = 0.5 - 0.633 * chiv - 0.33 * chiv * chiv
      bb = 0.877 * (1.0 - 2.0 * aa)

      proj = aa + bb * f
      extkb = (aa + bb * f) / f
      zmew = 1.0 / bb * (1.0 - aa / bb * log((aa + bb) / aa))
      acss = scat / 2.0 * proj / (proj + f * bb)
      acss = acss * (1.0 - f * aa / (proj + f * bb) * log((proj &
         + f * bb + f * aa) / (f * aa)))

      upscat = green * tran1 + (1.0 - green) * tran2
      upscat = 0.5 * (scat + (scat - 2.0 * upscat) &
         * ((1.0 - chiv) / 2.0) ** 2)
      betao = (1.0 + zmew * extkb) / (scat * zmew * extkb) * acss

!c----------------------------------------------------------------------
!c
!c     intermediate variables identified in appendix of SE-85.
!c
!c      be          (b)     : appendix      , SE-85
!c      ce          (c)     : appendix      , SE-85
!c      bot         (sigma) : appendix      , SE-85
!c      hh1         (h1)    : appendix      , SE-85
!c      hh2         (h2)    : appendix      , SE-85
!c      hh3         (h3)    : appendix      , SE-85
!c      hh4         (h4)    : appendix      , SE-85
!c      hh5         (h5)    : appendix      , SE-85
!c      hh6         (h6)    : appendix      , SE-85
!c      hh7         (h7)    : appendix      , SE-85
!c      hh8         (h8)    : appendix      , SE-85
!c      hh9         (h9)    : appendix      , SE-85
!c      hh10        (h10)   : appendix      , SE-85
!c      psi         (h)     : appendix      , SE-85
!c      zat         (l-t)   : appendix      , SE-85
!c      epsi        (s1)    : appendix      , SE-85
!c      ek          (s2)    : appendix      , SE-85
!c--------------------------------------------------------------------

      be = 1.0 - scat + upscat
      ce = upscat
      bot = (zmew * extkb) ** 2 + (ce**2 - be**2)
      if (abs(bot) <= 1.e-10) then
         scat = scat* 0.98
         be = 1.0 - scat + upscat
         bot = (zmew * extkb) ** 2 + (ce**2 - be**2)
      end if

      de = scat * zmew * extkb * betao
      fe = scat * zmew * extkb * (1.0 - betao)
      hh1 = - de * be + zmew * de * extkb - ce * fe
      hh4 = - be * fe - zmew * fe * extkb - ce * de

      if (be**2 - ce**2 < 0.0) then
         print*, be, ce, scat, upscat, bot, zmew, green, tran1, reff1, reff2
         print*, tran(iwave,1), scov, iwave, fmelt, tran(iwave,2)
         print*, snoww(1), satcap(1)
      endif
      psi = sqrt(be**2 - ce**2) / zmew

      zat = zlt / vcover * canex

      power1 = min(psi * zat, 50.0)
      power2 = min(extkb * zat, 50.0)
      epsi = exp(-power1)
      ek = exp (-power2)

      albedo(2,iwave,1) = soref(iwave) * (1.0 - areas) &
         + (1.2 - iwave * 0.4) * fmelt * areas
      albedo(2,iwave,2) = soref(iwave) * (1.0 - areas) &
         + (1.2 - iwave * 0.4) * fmelt * areas
      ge = albedo(2,iwave,1) / albedo(2,iwave,2)

!c-----------------------------------------------------------------------
!c     calculation of diffuse albedos
!c
!c      albedo(1,ir,2) ( i-up ) : appendix , SE-85
!c
!c-----------------------------------------------------------------------

      f1 = be - ce / albedo(2,iwave,2)
      zp = zmew * psi

      den = (be + zp) * (f1 - zp) / epsi - (be - zp) * (f1 + zp) * epsi
      hh7 = ce * (f1 - zp) / epsi / den
      hh8 = -ce * (f1 + zp) * epsi / den
      f1 = be - ce * albedo(2,iwave,2)
      den = (f1 + zp) / epsi - (f1 - zp) * epsi

      hh9   = (f1 + zp) / epsi / den
      hh10  = -(f1 - zp) * epsi / den
      tranc2(iwave) = hh9 * epsi + hh10 / epsi

      albedo(1,iwave,2) = hh7 + hh8

!c-----------------------------------------------------------------------
!c     calculation of direct albedos and canopy transmittances.
!c
!c      albedo(1,iw,1) ( i-up ) : equation(11)   , SE-85
!c      tranc (iw)   ( i-down ) : equation(10)   , SE-85
!c
!c-----------------------------------------------------------------------

      f1 = be - ce / albedo(2,iwave,2)
      zmk = zmew * extkb

      den = (be + zp) * (f1 - zp) / epsi - (be - zp) * (f1 + zp) * epsi
      hh2 = (de - hh1 / bot * (be + zmk)) * (f1 - zp) / epsi &
         - (be - zp) * (de - ce * ge - hh1 / bot * (f1 + zmk)) * ek
      hh2 = hh2 / den
      hh3 = (be + zp) * (de - ce * ge - hh1 / bot * (f1 + zmk)) * ek &
         - (de - hh1 / bot * (be + zmk)) * (f1 + zp) * epsi
      hh3 = hh3 / den
      f1 = be - ce * albedo(2,iwave,2)
      den = (f1 + zp) / epsi - (f1 - zp) * epsi
      hh5 = - hh4 / bot * (f1 + zp) / epsi &
         - (fe + ce * ge * albedo(2,iwave,2) + hh4 / bot * (zmk - f1)) * ek
      hh5 = hh5 / den
      hh6 = hh4 / bot * (f1 - zp) * epsi &
         + (fe + ce * ge * albedo(2,iwave,2) + hh4 / bot * (zmk - f1)) * ek
      hh6 = hh6 / den
      tranc1(iwave) = ek
      tranc3(iwave) = hh4 / bot * ek + hh5 * epsi + hh6 / epsi

      albedo(1,iwave,1) = hh1 / bot + hh2 + hh3

!c----------------------------------------------------------------------
!c     calculation of terms which multiply incoming short wave fluxes
!c     to give absorption of radiation by canopy and ground
!c
!c      radfac      (f(il,imu,iv)) : equation (19,20) , SE-86
!c
!c----------------------------------------------------------------------

      radfac(2,iwave,1) = (1.0 - vcover) * (1.0 - albedo(2,iwave,1)) &
         + vcover * (tranc1(iwave) * (1.0 - albedo(2,iwave,1)) &
         + tranc3(iwave) * (1.0 - albedo(2,iwave,2)))

      radfac(2,iwave,2) = (1.0 - vcover) * (1.0 - albedo(2,iwave,2)) &
         + vcover * tranc2(iwave) * (1.0 - albedo(2,iwave,2))

      radfac(1,iwave,1) = vcover * ((1.0 - albedo(1,iwave,1)) &
         - tranc1(iwave) * (1.0 - albedo(2,iwave,1)) &
         - tranc3(iwave) * (1.0 - albedo(2,iwave,2)))

      radfac(1,iwave,2) = vcover * ((1.0 - albedo(1,iwave,2)) &
         - tranc2(iwave) * (1.0 - albedo(2,iwave,2)))

!c----------------------------------------------------------------------
!c     calculation of total surface albedos ( salb ) with weighting
!c     for cover fractions.
!c----------------------------------------------------------------------
! albedo da sperficie eh o albedo de vegetacao e o solo as radiacoes diretas 
! e difusas no espectro visivel e infraverm.  
      do irad = 1, 2
         salb(iwave,irad) = (1.0 - vcover) * albedo(2,iwave,irad) &
            + vcover * albedo(1,iwave,irad)
      end do

! 1000  CONTINUE
   end do   ! (iwave=1,2)

!c----------------------------------------------------------------------
!c
!c     calculation of long-wave flux terms from canopy and ground
!c
!c      closs ( fc - rnc )     : equation (21),  SE-86
!c      gloss ( fg - rng )     : equation (22),  SE-86
!c
!c----------------------------------------------------------------------

   tgs = min(tf,tg) * areas + tg * (1.0 - areas)
   tc4 = tc * tc * tc * tc
   tg4 = tgs * tgs * tgs * tgs

   zkat = 1.0 / zmew * zlt / vcover
   zkat = min(50.0, zkat)
   zkat = max(1.e-5, zkat)
   thermk = exp(-zkat)

   fac1 = vcover * (1.0 - thermk)
   fac2 = 1.0
   closs = 2.0 * fac1 * stefan * tc4
   closs = closs - fac2 * fac1 * stefan * tg4
   gloss = fac2 * stefan * tg4
   gloss = gloss - fac1 * fac2 * stefan * tc4
   zlwup = fac1 * stefan * tc4 + (1.0 - fac1) * fac2 * stefan * tg4

   call longrn(tranc1, tranc2, tranc3)

!c-----------------------------------------------------------------------
!c
!c     calculation of absorption of radiation by surface
!c
!c-----------------------------------------------------------------------

   radt(1) = 0.0
   radt(2) = 0.0

   do iveg  = 1, 2
      do iwave = 1, 2
         do irad  = 1, 2
            !radt <- radiação total absorvida pelo dossel e o solo
            !radfac(iveg,iwave,irad) <- fração de radiação absorvida para cada elemento (dossel/solo)
            !                                                                           (visivel/infra)
            !                                                                           (directa/difussa)
            !radn <- radiação global observada particionada na (visível/infraverm) 
            !                                                  (direta/difussa)
            radt(iveg) = radt(iveg) + radfac(iveg,iwave,irad) * radn(iwave,irad)

         end do
      end do
   end do
    ! Adicionando a radiação térmica para cada um dos termos de radiação
   radt(1) = radt(1) + radn(3,2) * vcover * (1.0 - thermk) - closs
   radt(2) = radt(2) + radn(3,2) * (1.0 - vcover * (1.0 - thermk)) - gloss
end subroutine rada2



!c=======================================================================

subroutine longrn(tranc1, tranc2, tranc3)

!c=======================================================================
!c
!c     calculation of downward longwave. this is not required in gcm if
!c     downward longwave is provided by gcm-radiation code as radn(3,2).
!c
!c-----------------------------------------------------------------------

   use atmos, only : radn, cloud, em, rnetm, swdown, tm
   use const, only : stefan
   use donor, only : zlwup
   use radabs, only : albedo
   use site, only : rab
   use steps, only : ilw
   use vstate, only : vcover

   implicit none

   ! Dummy arguments
   real, intent(in) :: tranc1(2), tranc2(2), tranc3(2)

   ! Local variables
   integer :: iwave
   real :: esky, swab, swup

   if (ilw == 1) then

!c----------------------------------------------------------------------
!c     downward long-wave assumed to be provided as radn(3,2)
!c     no calculation is needed
!c----------------------------------------------------------------------

   else if (ilw == 2) then

!c----------------------------------------------------------------------
!c     downward long-wave from brunt's equation, Monteith(1973), p37.
!c----------------------------------------------------------------------

!      esky = 0.53 + 0.06 * sqrt(em)
!      radn(3,2) = esky * (1.0 + 0.2 * (cloud * cloud)) * stefan * tm ** 4

      esky = 0.5 + 0.06 * sqrt(em)
      radn(3,2) = esky * (1.0 + 0.22  * cloud) * stefan * tm ** 4

   else if (ilw == 3) then
 
!c----------------------------------------------------------------------
!c     downward long-wave flux calculated as residual from measured
!c     net radiation and outgoing longwave radiation.
!c
!c     calculation of absorbed fractions of radiation ( expendable )
!c----------------------------------------------------------------------

      do iwave = 1, 2
         rab(2,iwave,1) = (1.0 - vcover) * radn(iwave,1) &
            * (1.0 - albedo(2,iwave,1))
         rab(2,iwave,2) = (1.0 - vcover) * radn(iwave,2) &
            * (1.0 - albedo(2,iwave,2))

         rab(2,iwave,1) = rab(2,iwave,1) + vcover &
            * (radn(iwave,1) * (tranc1(iwave) * (1.0 - albedo(2,iwave,1)) &
            + tranc3(iwave) * (1.0 - albedo(2,iwave,2))))
         rab(2,iwave,2) = rab(2,iwave,2) + vcover &
            * radn(iwave,2) * tranc2(iwave) * (1.0 - albedo(2,iwave,2))

         rab(1,iwave,1) =  vcover &
            * radn(iwave,1) * ((1.0 - albedo(1,iwave,1)) &
            - tranc1(iwave) * (1.0 - albedo(2,iwave,1)) &
            - tranc3(iwave) * (1.0 - albedo(2,iwave,2)))
         rab(1,iwave,2) =  vcover &
            * radn(iwave,2) * ((1.0 - albedo(1,iwave,2)) &
            - tranc2(iwave) * (1.0 - albedo(2,iwave,2)))
      end do

      swab = rab(1,1,1) + rab(1,1,2) + rab(1,2,1) + rab(1,2,2) &
         + rab(2,1,1) + rab(2,1,2) + rab(2,2,1) + rab(2,2,2)
      swup = swdown - swab
      radn(3,2) = rnetm - swab + zlwup
!c----------------------------------------------------------------------

   else if (ilw == 4) then

!c----------------------------------------------------------------------
!c     downward long-wave from Idso and Jackson
!c----------------------------------------------------------------------
!c
      esky = 0.26 * exp(-0.00077 * ((273.0 - tm) ** 2))
      radn(3,2) = stefan * tm ** 4 * (1.0 - esky) * (1.0 + 0.2 * (cloud * cloud))

   else if (ilw == 5) then

!c----------------------------------------------------------------------
!c     downward long-wave from Swinbank
!c----------------------------------------------------------------------

      esky = 0.000009
      radn(3,2) = esky * stefan * tm ** 6 * (1.0 + 0.2 * (cloud * cloud))

   end if
end subroutine longrn



!c=======================================================================
!c
subroutine begtem
!c
!c-----------------------------------------------------------------------
!c     core routine: calculation of canopy and ground temperature
!c     increments over time step, fluxes derived.
!c
!c-----------------------------------------------------------------------
!c
!c
!c     subroutinescalled : snow1
!c     -----------------
!c
!c++++++++++++++++++++++++++++++output+++++++++++++++++++++++++++++++++++
!c
!c       rstfac(2)      soil moisture stress factor
!c       rsoil          soil surface resistance (s m-1)
!c       hr             soil surface relative humidity
!c       wc             canopy wetness fraction
!c       wg             ground wetness fraction
!c       ccx            canopy heat capacity (j m-2 k-1)
!c       cg             ground heat capacity (j m-2 k-1)
!c
!c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!c
!c----------------------------------------------------------------------

   use atmos, only : psur, tm
   use const, only : hlat, clai, cpair, cw, g, psy, snofac, snomel, tf
   use grads, only : etc, etgs, getc, getgs, tgs
   use hydrol, only : areas, wc, wg, satcap
   use snow, only : rsnow, tsnow
   use soils, only : bee, phsat, poros
   use stepv, only : tc, tg, capac, www, snoww
   use stores
!   use soilij, only : sodep
   use surfrs, only : rstfac, hr, rsoil
   use vdyijt, only : zlt
   use vstate, only : phc
   use multicamada
   implicit none

   ! Local variables
   real :: e, ge, x
   real :: argg, fac, phroot, psit
   real :: dzwww_perfil,rstfac2
   integer :: i
   real :: theta_fc, theta_wp,sum_extfrac_new
   integer :: baker

!c----------------------------------------------------------------------
!c     e(x) is vapour pressure in mbars as a function of temperature
!c     ge(x) is d e(x) / d ( temp )
!c----------------------------------------------------------------------

   e(x) = exp(21.18123 - 5418.0 / x) / 0.622
   ge(x) = exp(21.18123 - 5418.0 / x) * 5418.0 / (x * x) / 0.622

!c----------------------------------------------------------------------

   call snow1

!c----------------------------------------------------------------------

   hlat = (3150.19 - 2.378 * tm) * 1000.0
   psy = cpair / hlat * psur / 0.622
   snofac = hlat/ (hlat + snomel / 1.e3)

!c----------------------------------------------------------------------
!c    calculation of canopy and ground heat capacities.
!c    n.b. this specification does not necessarily conserve energy when
!c    dealing with very large snowpacks.
!c----------------------------------------------------------------------

   ccx = zlt * clai + (snoww(1) + capac(1)) * cw
   cg = csoil + min(0.05, snoww(2) + capac(2)) * cw

!c----------------------------------------------------------------------
!c
!c----------------------------------------------------------------------
!c      calculation of ground surface temperature and wetness fractions
!c
!c----------------------------------------------------------------------

   tsnow = min (tf - 0.01, tg)
   rsnow = snoww(2) / (snoww(2) + capac(2) + 1.e-10)

   tgs = tsnow * areas + tg * (1.0 - areas)

   etc = e(tc)
   etgs = e(tgs)
   getc = ge(tc)
   getgs = ge(tgs)

   wc = min(1.0, (capac(1) + snoww(1)) / satcap(1))
   wg = max(0.0, capac(2) / satcap(2)) * 0.25
   wc = max(0.0, wc)
   wg = min(1.0, wg)

!c-----------------------------------------------------------------------
!c     calculation of soil moisture stress factor.
!c     average soil moisture potential in root zone (layer-2) used as
!c     source for transpiration.
!c
!c      phroot      (psi-r) : equation (47) , SE-86
!c      rstfac(2)  f(psi-l) :    "     (12) , SE-89
!c-----------------------------------------------------------------------
    ! real :: theta_fc, theta_wp
    
    ! Calculo do termo de estresse hídrico no solo seguindo 
    ! "Seasonal drought stress in the Amazon: Reconciling models and observations"
    ! JOURNAL OF GEOPHYSICAL RESEARCH, VOL. 113, G00B01, doi:10.1029/2007JG000644, 2008
    ! I. T. Baker, L. Prihodko, A. S. Denning, M. Goulden, S. Miller, and H. R. da Rocha
    baker = 0
    if(baker .eq. 1)then
    theta_fc = poros*(abs(phsat)/100/340)**(1/bee)     ! Digman, L. Hydrological Physics
      theta_wp = poros*(abs(phsat)/100/15000)**(1/bee)   ! Digman, L. Hydrological Physics
!         theta_fc = poros* ( (-15.0/9.8) / phsat ) ** (-1.0 / bee)
!         theta_wp = poros*( (-1500.0/9.8) / phsat ) ** (-1.0 / bee)
 !     !Correção na fração de raízes pelo conteúdo de água no solo 
        rstfac2 = 0.0
        extfrac_new = 0.0
      DO i=1,nlayers-1
        extfrac_new(i) = extfrac(i) *                          &
  		       ( 1 - theta_wp/(dzwww(i)*poros)) /    &
  		       ( 1 - theta_wp/theta_fc )
        extfrac_new(i) = max(0.0, extfrac_new(i))
        rstfac2 = rstfac2 + extfrac_new(i)
      ENDDO
!        ! Determinando o novo vetor de distribuição de raízes considerando 
!        ! o grau de saturação de cada camada
  	sum_extfrac_new = sum(extfrac_new)
        DO i=1,nlayers-1
  	extfrac_new(i) = extfrac_new(i)/sum_extfrac_new
        ENDDO
     
    else ! baker
	dzwww_perfil = 0.0
	      !phroot = phsat * max(0.02, www(2)) ** (-bee)  ! old (pre RHV)
	    do i=1,nlayers-1                                        !rhv
	      dzwww_perfil = dzwww_perfil + extfrac(i)*dzwww(i)    
	    enddo
	  phroot = phsat * max(0.02, dzwww_perfil )** (-bee)        !rhv
	  extfrac_new = extfrac
   endif ! not baker 
     
     phroot = max (phroot, -2.e3) 
     rstfac(2) = 1.0 / (1.0 + exp(0.02 *(phc - phroot)))
     rstfac(2) = max(0.0001, rstfac(2))
     rstfac(2) = min(1.0, rstfac(2))
   
   ! if( dzwww_perfil .gt. 0.1) print*, www(2),dzwww_perfil , rstfac(2)
!c----------------------------------------------------------------------
!c
!c      rsoil function from fit to FIFE-87 data.  soil surface layer
!c      relative humidity.
!c
!c      rsoil      (rsoil) : SE-92B (personal communication)
!c      hr         (fh)    : equation (66) , SE-86
!c
!c----------------------------------------------------------------------

   fac = min(www(1), 1.0)
   fac = max(fac, 0.02)
!c 3/96 changes
!      rsoil =  max (0.1, 694. - fac*1500.) + 23.6  !c old SELLERS
    rsoil = exp(8.206 - 4.255 * fac)              !c TANG  
   ! parametrizacao da umidade do solo
   psit = phsat * fac ** (-bee)               ! psi top layer
   argg = max(-10.0, (psit * g / 461.5 / tgs))
   hr = exp(argg)               ! relative humidity of soil pore space
end subroutine begtem




!c======================================================================

subroutine endtem(ipbl)

!c----------------------------------------------------------------------
!c
!c      calculation of ea, ta, ra, rb, rd and soil moisture stress
!c      for the beginning of the time step
!c
!c----------------------------------------------------------------------
!c
!c                        modifications
!c
!c     (1)     : change in cog1,cog2,cogs1 to allow soil evaporation
!c               from beneath ground cover canopy.
!c
!c     (2)     : change in temperature tgs to account for patchy snow.
!c               a bulk (area-weighted) temperature is substituted for
!c               all energy budget calculations. energy used to warm
!c               exposed patches is subtracted from snow evaporation.
!c
!c     (3)     : inclusion of randall-sellers backward implicit scheme
!c               for calculating dtc, dtg, dth, dqm. option remains to
!c               use original dtc, dtg scheme only using parameter ipbl.
!c
!c     (4)     : inclusion of integrated canopy  photosynthesis -
!c               conductance model. note that soil moisture stress is
!c               modelled as a chronic condition, while the relative
!c               humidity term is solved within a feedback loop.
!c               reference : SE-92A
!c
!c=======================================================================
!c
!c     subroutines  called : rasite --> unstab,stab,rafcal
!c     -------------------   rbrd
!c  phosib --> cycalc-sortin
!c  delrn
!c  delhf
!c  delef
!c  sibslv --> gauss
!c  dtcdtg
!c  newton
!c-----------------------------------------------------------------------
!c
!c++++++++++++++++++++++++++++++output+++++++++++++++++++++++++++++++++++
!c       ect            canopy transpiration (j m-2)
!c       eci            canopy interception loss (j m-2)
!c       egs            ground evaporation (j m-2)
!c       egi            ground interception loss (j m-2)
!c       ec             ect + eci
!c       eg             egs + egi
!c       hc             canopy sensible heat flux (j m-2)
!c       hg             ground sensible heat flux (j m-2)
!c       chf            canopy heat storage flux (j m-2)
!c       shf            soil heat storage flux (j m-2)
!c       fc             canopy dew indicator
!c       fg             ground dew indicator
!c       heaten         heat loss to snow melt process (j m-2)
!c
!c++++++++++++++++++++++++++diagnostics++++++++++++++++++++++++++++++++++
!c
!c       tsnow          snow temperature (k)
!c
!c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   use atmos, only : em, tm
   use aerorx, only : ra, rb, rd
   use caerod, only : ht
   use const, only : asnow, hlat, psy, rcp, snofac, timcon
   use delts, only : dqm, dtc, dtg, dth
   use donor, only : zlwup
   use flux, only : &
      chf, ec, eci, ect, eg, egi, egs, gflux, &
      hc, heaten, hg, shf
   use grads, only : ta, ea, etc, etgs, getc, getgs, tgs
   use hydrol, only : areas, wc, wg
   use radabs, only : radt, thermk
   use snow, only : rsnow, tsnow
   use soils, only : poros, zdepth
   use steps, only : dtt, itrunk
   use stepv, only : www, tc, td, tg, capac, snoww
   use stores
   use surfrs, only : fc, fg, hr, rsoil, rst
   use vstate, only : vcover

   use flxdif, only : &
      deadqm, deadtc, deadtg, &
      hcdtc, hcdtg, hcdth, hgdtc, hgdtg, hgdth, &
      rncdtc, rncdtg, rngdtc, rngdtg

   implicit none

   ! Dummy arguments
   integer, intent(in) :: ipbl

   ! Local variables
   integer :: ifirst, icount, i, iwalk, lx, nonpos, nox
   real :: arean, aven, coct, cogs1, cogs2, d1, darea, dewc, dewg
   real :: ecf, ecidif, ecod, ecpot, egf, egidif, egit, egpot, egsmax
   real :: emin, finc, hend, hrr, t1, t2, taen, y

!!!     NOVAS VARIAVEIS PARA O NEWTON - TESTE
   integer :: iter_aux
   real :: zinc_aux, a2_aux, y1_aux

   ifirst = 1
   icount = 0

   fc = 1.0
   fg = 1.0
   ta = (tgs + tc) / 2.0
   ea = em
   ht = 0.0

   do while (icount <= 4)
      icount = icount + 1
      call rasite
      call rbrd
   end do

   call phosib

   ifirst = 0
   call delrn


!c----------------------------------------------------------------------
!c
!c     dew calculation : dew condition is set at beginning of time step.
!c     if surface changes state during time step, latent heat flux is
!c     set to zero.
!c
!c----------------------------------------------------------------------

   if (ea > etc) fc = 0.0
   if (ea > etgs) fg = 0.0

!c----------------------------------------------------------------------
!c     start of non-neutral resistance calculation loop
!c----------------------------------------------------------------------

   i = 0

!c----------------------------------------------------------------------
!c         initialize newton-raphson iterative routine
!c                    for ra calculation
!c----------------------------------------------------------------------


   nox = 0
   nonpos = 0
   iwalk = 0
   lx = 2
   finc = 50.0

   iter_aux = 0
   zinc_aux = 0.0
   a2_aux = 0.0
   y1_aux = 0.0

   do while (nox == 0)

      call rasite
      call delhf
      call delef

      if (ipbl == 0) then
         call dtcdtg
      else if (ipbl == 1) then
         call sibslv
      end if

!c-----------------------------------------------------------------------
!c     calculation of evaporative potentials (ecpot, egpot) and
!c        interception losses; fluxes in j m**-2.  ecidif and egidif
!c        hold the excess energy if all intercepted water is evaporated
!c        during the timestep.  these energies are later added to the
!c        sensible heat fluxes.
!c
!c      eci         (e-wc)  : equation (59) , SE-86
!c      egi         (e-wg)  : equation (60) , SE-86
!c-----------------------------------------------------------------------
!c
!c     check if interception loss term has exceeded canopy storage
!c----------------------------------------------------------------------

      ecpot = etc - ea + (getc - deadtc) * dtc - deadtg * dtg - deadqm * dqm
      eci = ecpot * wc / (2.0 * rb) * rcp / psy * dtt
      ecidif = max(0.0, eci - (snoww(1) + capac(1)) * 1.e3 * hlat)
      eci = min(eci, (snoww(1) + capac(1)) * 1.e3 * hlat)

      egpot = (etgs - ea) + (getgs - deadtg) * dtg - deadtc * dtc - deadqm * dqm
      egi = egpot / rd * rcp / psy * dtt * (wg * (1.0 - areas) + areas)
      egidif = max(0.0, egi - (snoww(2) + capac(2)) * 1.e3 * hlat) * (1.0 - rsnow)
      egit = min(egi, (snoww(2) + capac(2)) * 1.e3 * hlat) * (1.0 - rsnow)

!c----------------------------------------------------------------------
!c     calculation of interception loss from ground-snow. if snow patch
!c     shrinks, energy is taken from egi to warm exposed soil to tgs.
!c----------------------------------------------------------------------

      t1 = snoww(2) - 1.0 / asnow
      t2 = max(0.0, t1)
      aven = egi - t2 * hlat * 1.e3 / snofac
      if ((t1 - t2) * egi > 0.0) aven = egi
      darea = aven / ((tsnow - tg) * csoil - 1.0 / asnow * hlat * 1.e3 / snofac)
      arean = areas + darea
      egidif = egidif - min(0.0, arean) / asnow * hlat * 1.e3 / snofac * rsnow
      darea = max(darea, -areas)
      darea = min(1.0 - areas, darea)
      heaten = (tsnow - tg) * csoil * darea * rsnow
      egit = egit + (egi - heaten - egidif) * rsnow
      egi = egit

!c----------------------------------------------------------------------

      d1 = 1.0 / ra + 1.0 / rb + 1.0 / rd
      taen = ((tgs + dtg) / rd + (tc + dtc) / rb + tm / ra) / d1
      hend = (taen - tm) * rcp / ra + (ecidif + egidif) / dtt
      y = ht - hend
      i = i + 1
      if ( i > itrunk ) goto 44771
!!!        CALL newton(ht,y,finc,nox,nonpos,iwalk,lx)
!!!      call newton_new(ht, y, finc, nox, nonpos, iwalk, lx, &
      call newton_new(ht, y, finc, nox, nonpos, iwalk, &
               zinc_aux, a2_aux, y1_aux, iter_aux)

   end do   ! (nox==0)
44771 continue

!c----------------------------------------------------------------------
!c     exit from non-neutral calculation
!c     evapotranspiration fluxes calculated first ( j m-2 )
!c
!c----------------------------------------------------------------------
!c     calculation of transpiration and soil evaporation fluxes for the
!c        end of the timestep. see figure (2) of se-86.
!c
!c      ect         (e-c)   : equation (64) , SE-86
!c      egs         (e-s)   : equation (66) , SE-86
!c----------------------------------------------------------------------

   hrr = hr
   if (fg < 0.5) then
      hrr = 1.0
   else
      hrr = hr
   end if
!c.......................................
!c! Tatstch 10/22/2009 these piece is not in the original sib2
   emin = -1000.0
   ecod = eci
   eci = max(emin, eci)
   ecidif = ecidif + ecod - eci
   ecod = egi
   egi = max(emin, egi)
   egidif = egidif + ecod - egi
!c.......................................

   coct = (1.0 - wc) / (rst * fc + 2.0 * rb)
   cogs1 = (1.0 - areas) / (rd + rsoil * fg) * (1.0 - wg) * hrr
   cogs2 = cogs1 / hrr

   ect = ecpot * coct * rcp/psy * dtt
!c.......................................
!c! Tatstch 10/22/2009 these peace is not in the original sib2
   ecidif = ecidif + min(0.0, ect)                                ! add
   ect = max(ect, 0.0)                                           ! add
!c.......................................
   egs = (etgs + getgs * dtg) * cogs1 &
      - (ea + deadtg * dtg + deadtc * dtc + deadqm * dqm ) * cogs2
   egs = egs * rcp / psy * dtt
   egsmax = www(1) / 2.0 * zdepth(1) * poros * hlat * 1000.0
   egidif = egidif + max(0.0, egs - egsmax)
   egs = min(egs, egsmax)
!c.......................................
!c! Tatstch 10/22/2009 these peace is not in the original sib2
   egidif = egidif + min(0.0, egs)                                 ! add
   egs = max (egs, 0.0)                                            ! add

!c----------------------------------------------------------------------
!c     sensible heat flux calculated with latent heat flux correction
!c----------------------------------------------------------------------
!c
!c     calculation of sensible heat fluxes for the end of the timestep.
!c        see figure (2) of se-86.  note that interception loss excess
!c        energies (ecidif, egidif) are added.
!c
!c      hc          (hc)    : equation (63) , SE-86
!c      hg          (hgs)   : equation (65) , SE-86
!c----------------------------------------------------------------------

   hc = hc + (hcdtc * dtc + hcdtg * dtg + hcdth * dth) * dtt + ecidif
   hg = hg + (hgdtc * dtc + hgdtg * dtg + hgdth * dth) * dtt + egidif

!c----------------------------------------------------------------------
!c     test of dew condition. latent heat fluxes set to zero if sign
!c     of flux changes over time step.excess of energy donated to sensible
!c     heat flux.
!c     calculation of total latent heat fluxes,  see figure (2), se-86.
!c
!c      ec          (ec)    : equation (63) , SE-86
!c      eg          (eg)    : equation (65) , SE-86
!c----------------------------------------------------------------------

   ecf = sign(1.0, ecpot)
   egf = sign(1.0, egpot)
   dewc = fc * 2.0 - 1.0
   dewg = fg * 2.0 - 1.0

   if(dewc * ecf <= 0.0) then
      hc = hc + eci + ect
      eci = 0.0
      ect = 0.0
   end if
   if(dewg * egf <= 0.0) then
      hg = hg + egs + egi
      egs = 0.0
      egi = 0.0
   end if

!c! Tatsch 10/22/2009
!c! canopy evapotranpiration and soil evaporation
   ec = eci + ect
   eg = egi + egs

!c----------------------------------------------------------------------
!c     adjustment of : temperatures and vapor pressure
!c                     net radiation terms
!c                     storage heat fluxes
!c                     longwave loss and effective surface temperature
!c
!c----------------------------------------------------------------------

   ta  = taen
   ea = ea + deadtc*dtc + deadtg*dtg

   radt(1) = radt(1) + rncdtc * dtc + rncdtg * dtg
   radt(2) = radt(2) + rngdtc * dtc + rngdtg * dtg

!c----------------------------------------------------------------------
!c     calculation of storage heat fluxes
!c
!c----------------------------------------------------------------------

   chf = ccx / dtt * dtc
!c.......................................
!c! Tatstch 10/22/2009 these change there is not in the HRR's sib2 code
!c 3/96 changes
!c old    shf = cg / dtt * dtg + timcon*cg*2. * ( tgs+dtg - td )
   shf = cg / dtt * dtg + timcon * csoil * 2.0 * (tgs + dtg - td)
!c.......................................
!c! added by Tatsch based on HRR's sib2 code
!c! NOTA: sera que nao devia ser csoil ao inves de cg?

   gflux = timcon * cg * 2.0 * (tgs + dtg - td)

   zlwup = zlwup - rncdtc * dtc / 2.0 &
      - rngdtg * dtg * (1.0 - vcover * (1.0 - thermk))
end subroutine endtem



!c======================================================================
!c
!c        *********    auxiliary subroutine     **********
!c
!c=======================================================================
!c
subroutine rbrd
!c
!c=======================================================================
!c
!c      calculation of rb and rd as functions of u2 and temperatures
!c
!c
!c++++++++++++++++++++++++++++++output+++++++++++++++++++++++++++++++++++
!c
!c       rb (grb)       canopy to cas aerodynamic resistance (s m-1)
!c       rd (grd)       ground to cas aerodynamic resistance (s m-1)
!c       ta (gta)       cas temperature (k)
!c                      cas : canopy air space
!c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   use aerorx, only : ra, rb, rd
   use atmos, only : tm
   use caerod, only : ht
   use const, only : g, rcp
   use grads, only : ta, tgs, u2
   use hydrol, only : areas, wg
   use rause, only : rbc, rdc
   use stepv, only : tc
   use surfrs, only : cog1, cog2, fg, hr, rsoil
   use vdyijt, only : zlt
   use vstate, only : z2

   implicit none

   ! Local variables
   real :: cogs, cogr, d1, fac, fih, temdif

!c-----------------------------------------------------------------------
!c     rb : equation (a9), SE-86
!c-----------------------------------------------------------------------

   temdif = max(0.1, tc - tm)
   fac = zlt / 890.0 * sqrt(sqrt(temdif / 0.05))
   rb  = 1.0 / (sqrt(u2) / rbc + fac)

!c-----------------------------------------------------------------------
!c     rd : equation (a15), se-86
!c-----------------------------------------------------------------------

   temdif = max(0.1, tgs - tm)
   fih = sqrt(1.0 + 9.0 * g * temdif * z2 / tgs / (u2 * u2))
   rd  = rdc / u2 / fih

!c-----------------------------------------------------------------------
!c     calculation of ta, ht and effective soil aero-surface conductance,
!c        see equation (66) of SE-86 with snow area term added
!c
!c-----------------------------------------------------------------------

   d1 = 1.0 / ra + 1.0 / rb + 1.0 / rd
   ta = ( tgs/rd + tc/rb + tm/ra ) / d1
   ht = ( ta - tm ) * rcp / ra

   cogr = (1.0 - wg) / (rsoil * fg + rd)
   cogs =  wg / rd
   cog1 = (cogs + cogr * hr) * (1.0 - areas) + areas / rd
   cog2 = (cogs + cogr) * (1.0 - areas) + areas / rd
end subroutine rbrd



!c=======================================================================
subroutine phosib
!c=======================================================================
!c
!c     calculation of canopy photosynthetic rate using the integrated
!c     model relating assimilation and stomatal conductance.
!c     method uses numerical solution based on extrapolation from error
!c     versus co2i line.
!c     units are converted from mks to biological units in this routine.
!c     base reference :  SE-92A
!c
!c units
!c-------
!c
!c      pco2m, pco2a, pco2i, po2m                : pascals
!c      co2a, co2s, co2i, h2oa, h2os, h2oa       : mol mol-1
!c      vmax0, respn, assim, gs, gb, ga, pfd     : mol m-2 s-1
!c      effcon                                   : mol co2 mol quanta-1
!c      gcan, 1/rb, 1/ra, 1/rst                  : m s-1
!c      evapkg                                   : kg m-2 s-1
!c      q                                        : kg kg-1
!c
!c                       conversions
!c                      -------------
!c
!c      1 mol h2o           = 0.018 kg
!c      1 mol co2           = 0.044 kg
!c      h2o (mol mol-1)     = ea / psur ( mb mb-1 )
!c      h2o (mol mol-1)     = q*mm/(q*mm + 1)
!c      gs  (co2)           = gs (h2o) * 1./1.6
!c      gs  (mol m-2 s-1 )  = gs (m s-1) * 44.6*tf/t*p/po
!c      par (mol m-2 s-1 )  = par(w m-2) * 4.6*1.e-6
!c      mm  (molair/molh2o) = 1.611
!c
!c
!coutput
!c                      -------------
!c
!c      assimn              = canopy net assimilation rate
!c      ea                  = canopy air space vapor pressure
!c      1/rst               = canopy conductance
!c      pco2i               = internal co2 concentration
!c      respc               = canopy respiration
!c      respg               = ground respiration
!c      rstfac(4)      canopy resistance stress factors
!c
!c----------------------------------------------------------------------
!c
!c         rstfac(1) ( f(h-s) )               : equation (17,18), SE-92A
!c         rstfac(2) ( f(soil) )              : equation (12 mod), SE-89
!c         rstfac(3) ( f(temp) )              : equation (5b)   , CO-92
!c         rstfac(4) ( f(h-s)*f(soil)*f(temp))
!c
!c
!c----------------------------------------------------------------------

   use aerorx, only : ra, rb
   use atchem, only : pco2m, po2m
   use atmos, only : em, psur, tm, radn
   use carbio, only : assimn, gsh2o, pco2i, respc, respg
   use const, only : tf
   use grads, only : ea, etc, etgs, tgs
   use hydrol, only : wc
   use stepv, only : tc
   use surfrs, only : rstfac, cog1, cog2, rst
   use vderiv, only : gmudmu, vmax0
   use vdyijt, only : fparc, green, zlt
   use vstate, only : &
      atheta, binter, btheta, effcon, gradm, hhti, hlti, &
      respcp, shti, slti, trda, trdm, trop, ref, tran

   implicit none

   ! Local variables
   integer :: ic, ic2
   real :: c3, c4, fparkk, gah2o, gammas, gbh2o, gog1, gog2
   real :: bintc, h2oa, h2oi, h2om, h2os, h2osl, omss, par, park
   real :: pfd, qt, range, respn, rrkk, scatg, scatp, spfy
   real :: temph, templ, tprcor, vm, zkc, zko

   real :: pco2y(6), eyy(6)
   real :: aux

!c----------------------------------------------------------------------

   respg = 0.0

!c----------------------------------------------------------------------

   if (effcon <= 0.07) then
      c3 = 0.0
      c4 = 1.0
   else
      c3 = 1.0
      c4 = 0.0
   end if

!c----------------------------------------------------------------------
!c
!c
!c     calculation of canopy par use parameter.
!c
!c     fparkk      (pi)     : equation (31) , SE-92A
!c-----------------------------------------------------------------------

   scatp = green * (tran(1,1) + ref(1,1)) + (1.0 - green) * (tran(1,2) + ref(1,2))
   scatg = tran(1,1) + ref(1,1)
   park = sqrt(1.0 - scatp) * gmudmu
   fparc = 1.0 - exp(-park * zlt)
   fparkk = fparc / park * green

!c-----------------------------------------------------------------------
!c
!c     q-10 temperature effects :
!c      qt          (qt)    : table (2)     , SE-92A
!c      qt for vm changed to 2.1
!c
!c-----------------------------------------------------------------------

   qt = 0.1 * (tc - trop)
   respn = respcp * vmax0 * rstfac(2)
   respc = respn * 2.0 ** qt / (1.0 + exp(trda * (tc - trdm)))
   vm = vmax0 * 2.1 ** qt

   templ = 1.0 + exp(slti * (hlti - tc))
   temph = 1.0 + exp(shti * (tc - hhti))
   rstfac(3) = 1.0 / (templ * temph)
   vm = vm / temph * rstfac(2) * c3 + vm * rstfac(2) * rstfac(3) * c4

!c-----------------------------------------------------------------------
!c
!c     Michaelis-Menten constants for co2 and o2, co2/o2 specificity,
!c     compensation point
!c
!c      zkc          (kc)     : table (2)     , SE-92A
!c      zko          (ko)     : table (2)     , SE-92A
!c      spfy         (s)      : table (2)     , SE-92A
!c      gammas       (gamma-*): table (2)     , SE-92A
!c      omss         (omega-s): equation (13) , SE-92A
!c      bintc        (b*zlt)  : equation (35) , SE-92A
!c
!c-----------------------------------------------------------------------

   zkc = 30.0 * 2.1 ** qt
   zko = 30000.0 * 1.2 ** qt
   spfy = 2600.0 * 0.57 ** qt
   gammas = 0.5 * po2m / spfy * c3
   pfd = 4.6 * 1.e-6 * gmudmu * (radn(1,1) + radn(1,2))

   h2oi = etc / psur
   h2oa = ea / psur
   h2om = em / psur
   h2osl = etgs / psur

   tprcor = tf * psur * 100.0 / 1.013e5

!   gbh2o = 0.5 / rb * 44.6 * tprcor / tc
!   gah2o = 1.0 / ra * 44.6 * tprcor / tm
!   gog1 = cog1 * 44.6 * tprcor / tgs
!   gog2 = cog2 * 44.6 * tprcor / tgs
   aux = 44.6 * tprcor
   gbh2o = 0.5 / rb * aux / tc
   gah2o = 1.0 / ra * aux / tm
   gog1 = cog1 * aux / tgs
   gog2 = cog2 * aux / tgs

   rrkk = zkc * (1.0 + po2m / zko) * c3 + vmax0 / 5.0 * (1.8 ** qt) * c4
   par = pfd * effcon * (1.0 - scatg)
   bintc = binter * zlt * green * max(0.1, rstfac(2))

   omss = (vmax0 / 2.0) * (1.8 ** qt) / templ * rstfac(2) * c3 &
      + rrkk * rstfac(2) * c4

!c-----------------------------------------------------------------------
!c
!c     first guess is midway between compensation point and maximum
!c     assimilation rate.
!c
!c-----------------------------------------------------------------------

   range = pco2m * (1.0 - 1.6 / gradm) - gammas

   do ic = 1, 6
      pco2y(ic) = 0.0
      eyy(ic) = 0.0
   end do

   do ic = 1, 6
      ic2 = ic

      call sortin(eyy, pco2y, range, gammas, ic)

      call cycalc(fparkk, vm, gradm, bintc, atheta, btheta, &
         gah2o, gbh2o, gog1, gog2, wc, &
         h2oi, h2om, h2osl, par, pco2m, psur, &
         gammas, respc, respg, rrkk, omss, c3, c4, &
         pco2y(ic), eyy(ic), gsh2o, assimn, h2os, h2oa)

      if (abs(eyy(ic)) < 0.1) goto 44772
   end do
44772 continue

   pco2i = pco2y(ic2)

   rstfac(1) = h2os / h2oi
   rstfac(4) = rstfac(1) * rstfac(2) * rstfac(3)
!   rst = min(1.e6, 1.0 / (gsh2o * tc / (44.6 * tprcor) ))
   rst = min(1.e6, 1.0 / (gsh2o * tc / aux))
   ea = h2oa * psur
end subroutine phosib



!c=======================================================================

subroutine cycalc(fparkk, vm, gradm, bintc, atheta, btheta, &
      gah2o, gbh2o, gog1, gog2, wc, &
      h2oi, h2om, h2osl, par, pco2m, psur, &
      gammas, respc, respg, rrkk, omss, c3, c4, &
      pco2i, eyy, gsh2o, assimn, h2os, h2oa)

!c=======================================================================
!c
!c     calculation equivalent to steps in figure 4 of SE-92A
!c     c4 calculation based on CO-92.
!c
!c-----------------------------------------------------------------------
!c
!c++++++++++++++++++++++++++++++output+++++++++++++++++++++++++++++++++++
!c
!c       pco2i          canopy internal co2 concentration (mol mol-1)
!c       gsh2o          canopy conductance (mol m-2 s-1)
!c       h2os           canopy surface h2o concentration (mol mol-1)
!c
!c++++++++++++++++++++++++++diagnostics++++++++++++++++++++++++++++++++++
!c
!c       omc            rubisco limited assimilation (mol m-2 s-1)
!c                        (omega-c): equation (11) , SE-92A
!c       ome            light limited assimilation (mol m-2 s-1)
!c                        (omega-e): equation (12) , SE-92A
!c       oms            sink limited assimilation (mol m-2 s-1)
!c       co2s           canopy surface co2 concentration (mol mol-1)
!c                        equation (18c) , SE-92
!c       assimn         (a-n)    :  equation (14,15), SE-92A
!c
!c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   implicit none

   ! Dummy arguments
   real, intent(in) :: fparkk, vm, gradm, bintc, atheta, btheta
   real, intent(in) :: gah2o, gbh2o, gog1, gog2, wc
   real, intent(in) :: h2oi, h2om, h2osl, par, pco2m, psur
   real, intent(in) :: gammas, respc, respg, rrkk, omss, c3, c4, pco2i
   real :: assimn, gsh2o, h2os, h2oa, eyy

   ! Local variables
   real alpha, aquad, assim, assmt, beta, bquad, co2s, co2st
   real cquad, div2, div3, hcdma, omc, ome, omp, oms, pco2a
   real pco2in, sqrtin

   omc = vm  * (pco2i - gammas) / (pco2i + rrkk) * c3 + vm * c4
   ome = par * (pco2i - gammas) / (pco2i + 2.0 * gammas) * c3 + par * c4
   sqrtin= max(0.0, (ome + omc) ** 2 - 4.0 * atheta * ome * omc)
   omp = (ome + omc - sqrt(sqrtin)) / (2.0 * atheta)
   oms = omss * c3 + omss * pco2i * c4
   sqrtin = max(0.0, (omp + oms) ** 2 - 4.0 * btheta * omp * oms)
   assim = (oms + omp - sqrt(sqrtin)) / (2.0 * btheta)
   assimn = (assim - respc) * fparkk

!c-----------------------------------------------------------------------
!c     gah2o bottom stopped to prevent negative values of pco2a
!c-----------------------------------------------------------------------

   pco2a = pco2m - (1.4 / max(0.446, gah2o) * (assimn - respg) * psur *100.0)
   co2s  = pco2a / (psur * 100.0) - 1.4 * assimn / gbh2o

   assmt = max(1.e-12, assimn)
   co2st = min(co2s, pco2a / (psur * 100.0))
   co2st = max(co2st, 1.0 / (psur * 100.0))

   div2  = gah2o + gbh2o + gog2
   hcdma = h2oi * co2st / (gradm * assmt)
   alpha = hcdma * gbh2o * (1.0 - wc) / div2
   beta  = (-hcdma * bintc * gbh2o * (1.0 - wc) + h2oi * gbh2o * wc &
      + h2osl * gog1 + h2om * gah2o) / div2
   aquad = hcdma
   bquad = gbh2o * (hcdma - alpha) - h2oi - bintc * hcdma
   cquad = -gbh2o * (hcdma * bintc + beta)

   sqrtin = max(0.0, (bquad ** 2 - 4.0 * aquad * cquad))
   gsh2o = (-bquad + sqrt (sqrtin)) / (2.0 * aquad)
   h2os = (gsh2o - bintc) * hcdma
   h2os = min(h2os, h2oi)
   h2os = max(h2os, 1.e-7)
   gsh2o = h2os / hcdma + bintc
   div3 = gbh2o * gsh2o / (gbh2o + gsh2o) * (1.0 - wc) + gbh2o * wc + gog2 + gah2o
   h2oa = ((h2oi - h2oi * gsh2o / (gbh2o + gsh2o)) * gsh2o * (1.0 - wc) &
      + h2oi * gbh2o * wc + h2osl * gog1 + h2om * gah2o ) / div3
   h2os = (h2oa * gbh2o + h2oi * gsh2o) / (gbh2o + gsh2o)

!c-----------------------------------------------------------------------
!c     implied value of co2i derived from assimilation rate and stomatal
!c     conductance.
!c-----------------------------------------------------------------------

   pco2in =(co2s - 1.6 * assimn / gsh2o ) * psur * 100.0
   eyy = pco2i - pco2in
end subroutine cycalc




!c=======================================================================

subroutine sortin(eyy, pco2y, range, gammas, ic)

!c=======================================================================
!c
!c-----------------------------------------------------------------------

!c     arranges successive pco2/error pairs in order of increasing pco2.
!c     estimates next guess for pco2 using combination of linear and
!c     quadratic fits.
!c
!c-----------------------------------------------------------------------

   implicit none

   ! Dummy arguments
   integer, intent(in) :: ic
   real, intent(in) :: range, gammas
   real :: eyy(6), pco2y(6)

   ! Local variables
   integer :: i, i1, i2, i3, is, isp, ix, j, n
   real :: a, ac1, ac2, aterm, b, bc1, bc2, bterm, cc1, cc2, cterm
   real :: emin, pco2b, pco2yl, pco2yq, pmin

   if (ic < 4) then

      pco2y(1) = gammas + 0.5 * range
      pco2y(2) = gammas + range * (0.5 - 0.3 * sign(1.0, eyy(1)))
      pco2y(3) = pco2y(1) &
         - (pco2y(1) - pco2y(2)) / (eyy(1) - eyy(2) + 1.e-10) * eyy(1)

      pmin = min(pco2y(1), pco2y(2))
      emin = min(eyy(1), eyy(2))
      if (emin > 0.0 .and. pco2y(3) > pmin) pco2y(3) = gammas

   else

      n = ic - 1
      do j = 2, n
         a = eyy(j)
         b = pco2y(j)
         do i = j-1, 1, -1
            if (eyy(i) <= a) goto 44773
            eyy(i+1) = eyy(i)
            pco2y(i+1) = pco2y(i)
         end do
44773 continue
!c          i = 0
         eyy(i+1) = a
         pco2y(i+1) = b
      end do

      pco2b = 0.0
      is = 1
      do ix = 1, n
         if (eyy(ix) < 0.0) then
            pco2b = pco2y(ix)
            is = ix
         end if
      end do
      i1 = is - 1
      i1 = max0(1, i1)
      i1 = min0(n-2, i1)
      i2 = i1 + 1
      i3 = i1 + 2
      isp = is + 1
      isp = min0(isp, n)
      is = isp - 1

      pco2yl = pco2y(is) &
         - (pco2y(is) - pco2y(isp)) / (eyy(is) - eyy(isp)) * eyy(is)

!c----------------------------------------------------------------------
!c   method using a quadratic fit
!c----------------------------------------------------------------------

      ac1 = eyy(i1) * eyy(i1) - eyy(i2) * eyy(i2)
      ac2 = eyy(i2) * eyy(i2) - eyy(i3) * eyy(i3)
      bc1 = eyy(i1) - eyy(i2)
      bc2 = eyy(i2) - eyy(i3)
      cc1 = pco2y(i1) - pco2y(i2)
      cc2 = pco2y(i2) - pco2y(i3)
      bterm = (cc1 * ac2 - cc2 * ac1) / (bc1 * ac2 - ac1 * bc2)
      aterm = (cc1 - bc1 * bterm) / ac1
      cterm = pco2y(i2) - aterm * eyy(i2) * eyy(i2) - bterm * eyy(i2)
      pco2yq= cterm
      pco2yq= max(pco2yq, pco2b)
      pco2y(ic) = (pco2yl + pco2yq) / 2.0
   end if   ! (ic<4)

   pco2y(ic) = max (pco2y(ic), 0.01)
end subroutine sortin




!c======================================================================

subroutine delrn

!c======================================================================
!c
!c     partial derivatives of radiative and sensible heat fluxes
!c
!c----------------------------------------------------------------------

   use const, only : stefan
   use grads, only : tgs
   use radabs, only : thermk
   use stepv, only : tc
   use vstate, only : vcover

   use flxdif, only : rncdtc, rncdtg, rngdtg, rngdtc

   implicit none

   ! Local variables
   real :: fac1, fac2, tc3, tg3

   tc3 = tc * tc * tc
   tg3 = tgs * tgs * tgs
   fac1 = (1.0 - thermk) * vcover
   fac2 = 1.0

!   rncdtc = -2.0 * 4.0 * fac1 * stefan * tc3
   rncdtc = -8.0 * fac1 * stefan * tc3
   rncdtg = 4.0 * fac1 * fac2 * stefan * tg3

   rngdtg = -4.0 * fac2 * stefan * tg3
   rngdtc = 4.0 * fac1 * fac2 * stefan * tc3
end subroutine delrn



!c======================================================================

subroutine delhf

!c======================================================================
!c
!c     calculation of partial derivatives of canopy and ground sensible
!c     heat fluxes with respect to tc, tgs, and theta-m.
!c     calculation of initial sensible heat fluxes.
!c
!c========================================================================
!c
!c
!c       hc             canopy sensible heat flux (j m-2)
!c       hg             ground sensible heat flux (j m-2)
!c       hcdtc          dhc/dtc
!c       hcdtg          dhc/dtgs
!c       hcdth          dhc/dth
!c       hgdtc          dhg/dtc
!c       hgdtg          dhg/dtgs
!c       hgdth          dhg/dth
!c       aac            dh/dtc
!c       aag            dh/dtgs
!c       aam            dh/dth
!c
!c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!c
!c     fluxes expressed in joules m-2
!c
!c      hc  ( h-c ) : equation(63), SE-86
!c      hg  ( h-g ) : equation(65), SE-86
!c
!c----------------------------------------------------------------------

   use atmos, only : bps, tm
   use aerorx, only : ra, rb, rd
   use const, only : rcp
   use flux, only : hc, hg
   use grads, only : ta, tgs
   use steps, only : dtt
   use stepv, only : tc

   use flxdif, only : aac, aag, aam, hcdtc, hcdtg, hcdth, hgdtc, hgdtg, hgdth

   implicit none

   ! Local variables
   real :: d1, hcdqm, hgdqm

   d1 = 1.0 / ra + 1.0 / rb + 1.0 / rd
   ta = (tgs / rd + tc / rb + tm / ra) / d1

   hc = rcp * (tc - ta) / rb * dtt
   hg = rcp * (tgs - ta) / rd * dtt
!c----------------------------------------------------------------------
!c
!c      n.b.      fluxes expressed in joules m-2
!c
!c      hcdtc     (dhc/dtc) : equation (14) , SA-89B
!c      hcdtg     (dhc/dtgs): equation (14) , SA-89B
!c      hcdth     (dhc/dth) : equation (14) , SA-89B
!c      hgdtc     (dhg/dtc) : equation (15) , SA-89B
!c      hgdtg     (dhg/dtgs): equation (15) , SA-89B
!c      hgdth     (dhg/dth) : equation (15) , SA-89B
!c      aac       (dh/dtc)  : equation (12) , SA-89B
!c      aag       (dh/dtgs) : equation (12) , SA-89B
!c      aam       (dh/dth)  : equation (12) , SA-89B
!c----------------------------------------------------------------------

   hcdtc = rcp / rb * (1.0 / ra + 1.0 / rd) / d1
   hcdtg = -rcp / (rb * rd) / d1

   hgdtg = rcp / rd * (1.0 / ra + 1.0 / rb) / d1
   hgdtc = -rcp / (rd * rb) / d1   !!! OBS: HCDTG is equal to HGDTC

   hcdth= -rcp / (rb * ra) / d1 * bps
   hcdqm = 0.0

   hgdth = -rcp / (rd * ra) / d1 * bps
   hgdqm = 0.0

   aag = 1.0 / (rd * d1)
   aac = 1.0 / (rb * d1)
   aam = 1.0 / (ra * d1) * bps
end subroutine delhf



!c======================================================================

subroutine delef

!c======================================================================
!c
!c     calculation of partial derivatives of canopy and ground latent
!c     heat fluxes with respect to tc, tgs, theta-m, and qm.
!c     calculation of initial latent heat fluxes.
!c
!c      ec  ( e-c ) : equation(64), SE-86
!c      eg  ( e-gs) : equation(66), SE-86
!c
!c++++++++++++++++++++++++++++++output+++++++++++++++++++++++++++++++++++
!c
!c       ec             ect + eci
!c       eg             egs + egi
!c       ecdtc          dec/dtc
!c       ecdtg          dec/dtgs
!c       ecdqm          dec/dqm
!c       egdtc          deg/dtc
!c       egdtg          deg/dtgs
!c       egdqm          deg/dqm
!c       bbc            de/dtc
!c       bbg            de/dtgs
!c       bbm            de/dqm
!c
!c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   use aerorx, only : ra, rb, rd
   use atmos, only : em, psur
   use const, only : epsfac, psy, rcp
   use flux, only : ec, eg
   use grads, only : ea, etc, etgs, getc, getgs
   use hydrol, only : areas, wc, wg
   use steps, only : dtt
   use surfrs, only : cog1, cog2, fc, fg, hr, rsoil, rst

   use flxdif, only : &
      bbc, bbg, bbm, deadqm, deadtc, deadtg, demdqm, &
      ecdqm, ecdtc, ecdtg, egdqm, egdtc, egdtg

   implicit none

   ! Local variables
   real :: coc, cogr, cogs, d2, deadem, hrr, rcc, sh, top

!c----------------------------------------------------------------------
!c     modification for soil dryness : hr = rel. humidity in top layer
!c----------------------------------------------------------------------

   if (fg < 0.5) then
      hrr = 1.0
   else
      hrr = hr
   end if

!c-----------------------------------------------------------------------
!c
!c     calculation of surface resistance components, see equations (64,66)
!c       of SE-86
!c
!c-----------------------------------------------------------------------

   rcc = rst * fc + 2.0 * rb
   coc = (1.0 - wc) / rcc + wc / (2.0 * rb)

   cogr = (1.0 - wg) / (rsoil * fg + rd)
   cogs =  wg / rd
   cog1 = (cogs + cogr * hrr) * (1.0 - areas) + areas / rd
   cog2 = (cogs + cogr) * (1.0 - areas) + areas / rd

   d2 = 1.0 / ra + coc + cog2
   top = coc * etc + cog1 * etgs + em / ra
   ea = top / d2

   ec = (etc - ea) * coc * rcp / psy * dtt
   eg = (etgs * cog1 - ea * cog2) * rcp / psy * dtt

   deadtc = getc * coc / d2
   deadtg = getgs * cog1 / d2

!c-----------------------------------------------------------------------
!c      ecdtc     (dec/dtc) : equation (14) , SA-89B
!c      ecdtg     (dec/dtgs): equation (14) , SA-89B
!c      ecdqm     (dec/dqm) : equation (14) , SA-89B
!c      egdtc     (deg/dtc) : equation (15) , SA-89B
!c      egdtg     (deg/dtgs): equation (15) , SA-89B
!c      egdqm     (deg/dqm) : equation (15) , SA-89B
!c      bbc       (de/dtc)  : equation (13) , SA-89B
!c      bbg       (de/dtgs) : equation (13) , SA-89B
!c      bbm       (de/dqm)  : equation (13) , SA-89B
!c-----------------------------------------------------------------------

   ecdtc = (getc - deadtc) * coc * rcp / psy
   ecdtg = -deadtg * coc * rcp / psy

   egdtg = (getgs * cog1 - deadtg * cog2 ) * rcp / psy
   egdtc = -deadtc * cog2 * rcp / psy

   sh = epsfac * em / (psur - em)
   deadem = 1.0 / (ra * d2)
   demdqm = epsfac * psur / (epsfac + sh) ** 2
   deadqm = deadem * demdqm

   ecdqm = -deadqm * coc * rcp / psy
   egdqm = -deadqm * cog2 * rcp / psy

   bbg = (cog1 / d2) * getgs * epsfac * psur / (psur - etgs) ** 2
   bbc = (coc /d2) * getc * epsfac * psur / (psur - etc) ** 2
   bbm = 1.0 / (ra * d2)
end subroutine delef



!c======================================================================

subroutine sibslv

!c======================================================================
!c
!c     solve for time changes of pbl and sib variables,
!c     using a semi-implicit scheme.
!c
!c      dtc, dtg, dth, dqm  : equations(12-15) , SA-89B + radiation terms
!c
!c
!c++++++++++++++++++++++++++++++output+++++++++++++++++++++++++++++++++++
!c
!c       dtc            canopy temperature increment (K)
!c       dtg            ground surface temperature increment (K)
!c       dth            mixed layer potential temperature increment (K)
!c       dqm            mixed layer mixing ratio increment (kg kg-1)
!c
!c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   use aerorx, only : ra
   use atmos, only : bps
   use const, only : cpair, g, hlat, rhoair, timcon
   use delts, only : dqm, dtc, dtg, dth
   use flux, only : ec, eg, hc, hg
   use grads, only : tgs
   use radabs, only : radt
   use steps, only : dtt
   use stepv, only : td
   use stores, only : ccx, cg, csoil

   use flxdif, only : &
      aac, aag, aam, bbc, bbg, bbm, ecdqm, ecdtc, ecdtg, &
      egdqm, egdtc, egdtg, hcdtc, hcdtg, hcdth, hgdtc, &
      hgdtg, hgdth, rncdtc, rncdtg, rngdtc, rngdtg

   implicit none

   ! Local variables
   integer :: i, j
   real :: dthb, dwb, etem, fths, fws, grav2, gv, psb
   real :: pblsib(4,4), cvec(4), solvec(4), chin(4,5), work(4,5)

   grav2 = 0.01 * g
   gv = grav2 * rhoair / ra
   psb = 100.0

   etem = 0.0
   dthb = 0.0
   dwb = 0.0

!c     cvec uses provisional values of the fluxes.

   fths = (hc + hg) / (dtt * cpair * bps)
   fws = (ec + eg) / (dtt * hlat)

!c     tg equation

!c 3/5/96 changes
!c  old   pblsib(1,1)= cg/dtt + hgdtg+egdtg-rngdtg + timcon*cg*2.0
   pblsib(1,1) = cg / dtt + hgdtg + egdtg - rngdtg + timcon * csoil * 2.0
   pblsib(1,2) =            hgdtc + egdtc - rngdtc
   pblsib(1,3) =            hgdth
   pblsib(1,4) =                    egdqm

!c     tc equation

   pblsib(2,1) =           + hcdtg + ecdtg - rncdtg
   pblsib(2,2) = ccx / dtt + hcdtc + ecdtc - rncdtc
   pblsib(2,3) =             hcdth
   pblsib(2,4) =                     ecdqm

!c     theta equation

   pblsib(3,1) = -gv * aag
   pblsib(3,2) = -gv * aac
   pblsib(3,3) = -gv * (aam - 1.0) + etem + psb / dtt
   pblsib(3,4) = 0.0

!c     sh equation

   pblsib(4,1) = -gv * bbg
   pblsib(4,2) = -gv * bbc
   pblsib(4,3) = 0.0
   pblsib(4,4) = -gv * (bbm - 1.0) + etem + psb / dtt

!c 3/96 changes
!c  old  cvec(1) = radt(2) - hg/dtt - eg/dtt - timcon*(tgs-td)*cg*2.
   cvec(1) = radt(2) - hg / dtt - eg / dtt - timcon * (tgs - td) * csoil * 2.0
   cvec(2) = radt(1) - hc / dtt - ec / dtt
   cvec(3) = grav2 * fths + etem * dthb
   cvec(4) = grav2 * fws  + etem * dwb

!c     solve 4 x 4 matrix equation

   do j = 1, 4
      do i = 1, 4
         chin(i,j) = pblsib(i,j)
      end do
   end do

   do i = 1, 4
      chin(i,5) = cvec(i)
   end do

   call gauss(chin, 4, 5, solvec, work)

   dtg = solvec(1)
   dtc = solvec(2)
   dth = solvec(3)
   dqm = solvec(4)
end subroutine sibslv




!c======================================================================

subroutine dtcdtg

!c----------------------------------------------------------------------
!c
!c     calculation of temperature tendencies assuming no interaction
!c     with the pbl : equations(69,70), SE-86
!c
!c----------------------------------------------------------------------

   use const, only : timcon
   use delts, only : dtc, dtg
   use flux, only : ec, eg, hc, hg
   use grads, only : tgs
   use radabs, only : radt
   use steps, only : dtt
   use stepv, only : td
   use stores, only : ccx, cg, csoil

   use flxdif, only : &
      rncdtc, hcdtc, ecdtc, rncdtg, hcdtg, ecdtg, egdtc, &
      egdtg, hgdtc, hgdtg, rngdtc, rngdtg

   implicit none

   ! Local variables
   real :: ccodtc, ccodtg, ccorhs, gcodtg, gcodtc, gcorhs, denom

   ccodtc = ccx / dtt - rncdtc + hcdtc + ecdtc
   ccodtg = -rncdtg + hcdtg + ecdtg
   ccorhs = radt(1) - (hc + ec) / dtt

!c 3/96 changes
!c  old  gcodtg = cg / dtt + timcon*cg*2. - rngdtg + hgdtg + egdtg
   gcodtg = cg / dtt + timcon * csoil * 2.0 - rngdtg + hgdtg + egdtg
   gcodtc = -rngdtc + hgdtc + egdtc
!c 3/96 changes
!c old   gcorhs = radt(2) - timcon*cg*2. * ( tgs -td ) - ( hg + eg ) / dtt
   gcorhs = radt(2) - timcon * csoil * 2.0 * (tgs - td) &
      - (hg + eg) / dtt

   denom = ccodtc * gcodtg - ccodtg * gcodtc

   dtc = (ccorhs * gcodtg - ccodtg * gcorhs) / denom
   dtg = (ccodtc * gcorhs - ccorhs * gcodtc) / denom
end subroutine dtcdtg



!c=======================================================================

!!!      SUBROUTINE updat2
subroutine updat2

!c=======================================================================
!c
!c     updating of all prognostic variables.
!c
!c-----------------------------------------------------------------------
!c
!c     subroutines called   : updat2
!c     ------------------     snow2
!c   run2
!c
!c++++++++++++++++++++++++++++++output from this block++++++++++++++++++++
!c
!c       dtc            canopy temperature increment (K)
!c       dtd            deep soil temperature increment (K)
!c       dtg            ground surface temperature increment (K)
!c       www(3)         ground wetness
!c       capac(2)       canopy/ground liquid interception store (m)
!c       snoww(2)       canopy/ground snow interception store (m)
!c       roff           runoff (mm)
!c       etmass (fws)   evapotranspiration (mm)
!c       hflux (fss)    sensible heat flux (w m-2)
!c
!c++++++++++++++++++++++++++diagnostics++++++++++++++++++++++++++++++++++
!c
!c       ecmass         canopy evapotranspiration (mm)
!c       egmass         ground evapotranspiration (mm)
!c
!c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   use const, only : snofac, hlat
   use delts, only : dtc, dtd, dtg, dqm, dth
   use donor, only : etmass, hflux
   use flux, only : ect, eci, ecmass, egs, egi, egmass, hc, hg
   use snow, only : rsnow
   use soils, only : poros, zdepth
   use steps, only : dtt
   use stepv, only : snoww, capac, www, tc, tg, td
   use multicamada
   use morphology1, only : dsfc_cst
   IMPLICIT NONE

   ! Dummy arguments
!   integer, intent(in) :: idirr

   ! Local variables
   integer :: iveg,i
   real :: facks, dtdc, ectdif, egsdif, extrak, facl, qm, th
!   real :: facl3,extrak3,rootex
   real :: depwww,ectnew
   real :: ectPot
   th = 0.0
   qm = 0.0

   call snow1

!c----------------------------------------------------------------------
!c    interception losses.
!c    evaporation losses are expressed in j m-2 : when divided by
!c    ( hlat*1000.) loss is in m m-2. mass terms are in kg m-2 dt-1
!c    interception and drainage treated in inter2.
!c
!c----------------------------------------------------------------------

   rsnow = snoww(1) / (snoww(1) + capac(1) + 1.e-10)
   facks = 1.0 + rsnow * (snofac - 1.0)
   if (ect + eci <= 0.0) then
      eci = ect + eci
      ect = 0.0
      facks = 1.0 / facks
   end if
   capac(1) = capac(1) - (1.0 - rsnow) * eci * facks / hlat / 1.e3
   snoww(1) = snoww(1) - rsnow * eci * facks / hlat / 1.e3
   ecmass = eci * facks / hlat

   rsnow = snoww(2) / (snoww(2) + capac(2) + 1.e-10)
   facks = 1.0 + rsnow * (snofac - 1.0)
   if (egs + egi <= 0.0) then
      egi = egs + egi
      egs = 0.
      facks = 1.0 / facks
   end if
   capac(2) = capac(2) - (1.0 - rsnow) * egi * facks / hlat / 1.e3
   snoww(2) = snoww(2) - rsnow * egi * facks / hlat / 1.e3
   egmass = egi * facks / hlat

!c----------------------------------------------------------------------
!c    dumping of small capac values onto soil surface store
!c----------------------------------------------------------------------

   do iveg = 1, 2
      if (snoww(iveg) + capac(iveg) <= 0.00001) then
         www(1) = www(1) + (snoww(iveg) + capac(iveg)) / (poros * zdepth(1))
         capac(iveg) = 0.0
         snoww(iveg) = 0.0
      end if
   end do

!c----------------------------------------------------------------------
!c    snowmelt / refreeze calculation
!c----------------------------------------------------------------------

   call snow2

!c----------------------------------------------------------------------
!c    evapotranspiration losses,
!c    extraction of transpiration loss from root zone, soil evaporation.
!c
!c      ect         (e-dc)  : equation (5,6), SE-86
!c      egs         (e-s)   : equation (5)  , SE-86
!c----------------------------------------------------------------------
! CODIGO ORIGINAL SiB2
!   facl   = 1.0 / hlat / 1.e3 / (poros * zdepth(2))
!   extrak = ect * facl             !trans pot.
!   extrak = min(extrak, www(2))    !limitação da transp pela www(2)
!   ectdif = ect - extrak / facl    !transp pot - transp real
!   ect    = extrak / facl          !transp agora é a real 
!   hc     = hc + ectdif            !se teve energia não consumida em transp vira dH
!   ecmass = ecmass + ect / hlat    !somando a LE a transp final
!   www(2) = www(2) - ect * facl    !tirando da 2da camada água consumida pela transp
   
! MODIFICACAO PARA INCLUIR TRANSPIRACAO NO INVERNO DESDE A 
!3RA CAMADA ....
!    real :: facl3,extrak3

!   facl   = 1.0 / hlat / 1.e3 / (poros * zdepth(2))
!   extrak = ect * facl
!   extrak = min(extrak, www(2)) 
!   ectdif = ect - extrak / facl 
!   
!   facl3  = 1.0 / hlat / 1.e3 / (poros * zdepth(3))
!   extrak3 = ectdif * facl3
!   extrak3 = min(extrak3,rootex*www(3)) 
!   ectdif = ectdif - extrak3 / facl3 
!   
!   ect    = extrak / facl + extrak3 / facl3
!   hc     = hc + ectdif   
!   ecmass = ecmass + ect / hlat 
!   www(2) = www(2) - (extrak / facl) * facl 
!   www(3) = www(3) - (extrak3 / facl3) * facl3
   
!....FIN DA MODIFICACAO ROILAN
!   if(iter .lt. 10) print*,'before LE',(dzwww(i),i=1,10)
!   if(iter .lt. 10) print*,'before LE',(extfrac_new(i),i=1,10)
  ! Determinacao da transpiracao desde cada uma das subcamadas
  ! >
  ! Constantes:  real :: depthreal,extfrac(nlayers)
   ectnew = 0.0 ! Armazena transpiracao final
   ectPot = ect
   depwww = dsfc_cst ! Profundidade da 1ra camada
    facl   = 1.0 / hlat / 1.e3 / (poros * dzmultilayer)
    ! O mais importante do modelo multicamadas é a distribuição das...
    ! ... raízes. portanto EXTFRAC determina quanto da transpiração ...
    ! ... vai se extraída de cada camada. Um vetor com a fraçao de extração...
    ! ... com o mesmo numero de elementos que NLAYERS 
  DO i=1,nlayers
    ! convertendo a ECT em unidade de umidade para a camada pela fraçao 
    ! de ECTtotal que pode transpirar, segundo distribuiçaõ das raízes
    extrak = (extfrac_new(i) * ect) * facl
    ! calculando quanto transpirou para a camada (i)
    extrak = min(extrak, dzwww(i))
    ! se a umidade da camada satisfaz a demanda EXTRAK vira
    ! transpiracao, caso contrário limitado pela umidade da 
    ! camada.
    !>
    ! Determina, portanto, o incremento a Calor Sensível (H) 
    ectdif = (extfrac_new(i) * ect) - (extrak / facl)
    ! Adiciondo dH
    hc     = hc + ectdif
    ! Susbtraindo da umidade do solo o que foi transpirado
    dzwww(i) = dzwww(i) - extrak
    ! Adicionando transpiracao da camada na total
    ectnew = ectnew + (extrak / facl)
    ecmass = ecmass + (extrak / facl) / hlat
  ENDDO
    ! Traspiração total é a soma do total transpirado em cada camada
    !	print*,ect,ectnew,max(ect-ectnew,1.e-8)
    ect = ectnew
!RHV Fim das mudanças na tranpiração pelo modelo multicamadas
!       if(ect .gt. 800) print*, ectPot,ectnew
! Ir a fluxo intercamadas para (SUBROUTINE run2)
    
   facl   = 1.0 / hlat / 1.e3 / (poros * zdepth(1))
   extrak = egs * facl
   extrak = min(extrak, www(1))
   egsdif = egs - extrak / facl
   egs    = extrak / facl
   hg     = hg + egsdif
   egmass = egmass + egs / hlat
   www(1) = www(1) - egs * facl

!c-----------------------------------------------------------------------
!c    calculation of total moisture and sensible heat fluxes from surface.
!c-----------------------------------------------------------------------

   etmass = ecmass + egmass
   hflux  = (hc + hg) / dtt

!c----------------------------------------------------------------------
!c    calculation of interflow, infiltration excess and loss to
!c    groundwater .  all losses are assigned to variable 'roff' .
!c----------------------------------------------------------------------
!  write(9856,'(a8, i10, 50(f13.6))') 'Transp  ',iter , egs *3600/2454000/dtt, &
! 	(extfrac(i)*ectPot*3600/2454000/dtt,i=1,nlayers),  &
! 	(ectnew-ectPot)*3600/2454000/dtt,ectPot*3600/2454000/dtt ! WRITING
! write(9856,'(a8, i10, 50(f13.6))') 'EXTRAK- ',iter,www(1), &
! 				  (dzwww(i),i=1,nlayers) ! WRITING
! write(9856,*)	 '--'
!  if(iter .lt. 10) print*,'after LE',(dzwww(i),i=1,10) 
  call run2(ectPot)

!c----------------------------------------------------------------------
!c
!c    update of temperatures and pbl variables. note that tc and tg
!c    are modified in interc as a result of precipitation inputs.
!c    update of interception stores.
!c
!c----------------------------------------------------------------------

   dtg = max(-50.0, dtg)
   dtg = min(50.0, dtg)
   dtc = max(-50.0, dtc)
   dtc = min(50.0, dtc)
   dtd = max(-50.0, dtd)
   dtdc = min(50.0, dtd)

!c	dth=max(-50., dth)
!c	dth=min(50., dth)
!c	dqm=max(-50., dqm)
!c	dqm=min(50., dqm)

   tc  = tc + dtc
   tg  = tg + dtg
   td  = td + dtd
   th  = th + dth
   qm  = qm + dqm
end subroutine updat2



!c=======================================================================

subroutine snow2

!c=======================================================================
!c
!c    snowmelt / refreeze calculation
!c----------------------------------------------------------------------
!c
!c     calculation of snowmelt and modification of temperatures
!c
!c     modification deals with snow patches:
!c          ts < tf, tsnow = ts
!c          ts > tf, tsnow = tf
!c
!c-----------------------------------------------------------------------

   use const, only : asnow, cw, pie, snomel, tf
   use delts, only : dtc, dtd, dtg
   use flux, only : chf, shf
   use hydrol, only : areas
   use snow, only : tsnow
   use soils, only : poros, zdepth
   use steps, only : dtt
   use stepv, only : snoww, www, capac, tc, tg
   use stores, only : ccx, cg, csoil

   implicit none

   ! Local variables
   integer :: iveg
   real :: avex, avheat, avmelt, cct, cctt, cool, dts, dtsg
   real :: dtsg2, dtsg3, exheat, exmelt, fluxx, fluxef, freeze
   real :: heat, realc, realg, safe, snowhc, tbulk
   real :: tn, ts, zmelt, zmelt2

   do iveg = 1, 2

      realc = real(2 - iveg)
      realg = real(iveg - 1)

      cctt = realc * ccx + realg * cg
      cct = realc * ccx + realg * csoil
      ts = realc * tc + realg * tg
      dts = realc * dtc + realg * dtg
      fluxx = realc * chf + realg * dtg / dtt * cg

      tsnow = min (tf - 0.01, ts)
      snowhc = min(0.05, snoww(iveg)) * cw * realg
      zmelt = 0.0

      if (snoww(iveg) <= 0.0) then

         if (ts + dts <= tf) then

!c-----------------------------------------------------------------------
!c
!c     no snow  present, simple thermal balance with possible freezing.
!c
!c-----------------------------------------------------------------------

            freeze = min(0.0, fluxx * dtt - (tf - 0.01 - ts) * cctt)
            snoww(iveg) = min(capac(iveg), -freeze / snomel)
            zmelt = capac(iveg) - snoww(iveg)
            capac(iveg) = 0.0
            dts = dts + snoww(iveg) * snomel / cctt
         end if

      else

!c-----------------------------------------------------------------------
!c
!c     snow present
!c
!c-----------------------------------------------------------------------

         if (ts >= tf .or. ts + dts >= tf) then
            if (ts <= tf) then

!c-----------------------------------------------------------------------
!c
!c     snow present : ts < tf,  ts+dts > tf
!c
!c-----------------------------------------------------------------------

               avex = fluxx - (tf - 0.01 - ts) * cctt / dtt
               avmelt = (avex / snomel * (areas * realg + realc)) * dtt
               zmelt = min(avmelt, snoww(iveg))
               snoww(iveg) = snoww(iveg) - zmelt
               avheat = avex * (1.0 - areas) * realg &
                  + (avmelt - zmelt) * snomel / dtt

               safe = max(1.0 - areas * realg, 1.e-8)
               dts = tf - 0.01 - ts + avheat / (cctt * safe) * dtt

            else

!c-----------------------------------------------------------------------
!c
!c     snow present and ts > tf : ground only.
!c
!c-----------------------------------------------------------------------

               tbulk = tsnow * areas + ts * (1.0 - areas)
               tn = tbulk + dts
               exheat = cct * (1.001 - max(0.1, areas)) * dts
               exmelt = fluxx * dtt - exheat
               heat = exheat
               dtsg = exheat / (cct * (1.001 - areas))
               if (ts + dtsg <= tf) then
                  heat = (tf - 0.01 - ts) * cct * (1.0 - areas)
                  dtsg = tf - 0.01 - ts
               end if
               exmelt = exmelt + exheat - heat

               if (exmelt >= 0.0) then
                  zmelt = exmelt/snomel
                  if (asnow * (snoww(iveg) - zmelt) < 1.0) then
                     zmelt = max(0.0, snoww(iveg) - 1.0 / asnow)
                  end if
                  snoww(iveg) = snoww(iveg) - zmelt
                  exmelt = exmelt - zmelt * snomel
                  zmelt2 = exmelt / (cct * (ts - tf) * asnow + snomel)
                  zmelt2 = min(zmelt2, snoww(iveg))
                  zmelt = zmelt + zmelt2
                  snoww(iveg) = snoww(iveg) - zmelt2
                  exmelt = exmelt - zmelt2 * (cct * (ts - tf) * asnow + snomel)
                  dts  = dtsg + exmelt / cct
               else
                  cool = min(0.0, tf - 0.01 - (ts + dtsg)) * cct * (1.0 - areas)
                  dtsg2 = max(cool, exmelt) / (cct * (1.001 - areas))
                  exmelt = exmelt - dtsg2 * cct * (1.0 - areas)
                  dtsg3 =exmelt / cctt
                  dts = dtsg + dtsg2 + dtsg3
               end if
            end if   ! (ts<=tf)
         end if   ! (ts>=tf OR ts+dts>=tf)
      end if   ! (snoww(iveg)<=0.0)

      www(1) = www(1) + zmelt / (poros * zdepth(1))

      dtc = dtc * realg + dts * realc
      dtg = dtg * realc + dts * realg

   end do   ! (iveg=1, 2)

   fluxef = shf - cg * dtg / dtt
!c 3/96 changes
!c  old  dtd = fluxef / ( cg * 2. * sqrt ( pie*365. ) ) * dtt
   dtd = fluxef / (csoil * 2.0 * sqrt(pie * 365.0)) * dtt
end subroutine snow2



!c=======================================================================
!c  Consider Soil Groundwater interaction tangqh@iis.u-tokyo.ac.jp
!!!      SUBROUTINE run2
subroutine run2(ectmass_2)

!c=======================================================================
!c    calculation of interflow, infiltration excess and loss to
!c    groundwater .  all losses are assigned to variable 'roff' .
!c----------------------------------------------------------------------

   use const, only : tf, hlat
   use donor, only : roff
!   use grads, only : tgs
!   use hydrol, only : areas
   use output, only : roff2, roff3, roff4!, roffq3g, roffq3g2,gwsoil
!   use snow, only : tsnow
   use soils
   use steps, only : dtt!,iter
   use stepv, only : www, td!, tg
!   use temp, only : gwdep
   use multicamada
   use morphology1, only : dsfc_cst
   
   implicit none

   ! Dummy arguments
!   integer, intent(in) :: idirr

   ! Local variables
   integer :: i
   real :: deficit!,avk, avkmin, avkmax, d3g,  denom, div
   real :: pows, excess!,dpdw, dpdwdz,  pmin, pmax,  props
!   real :: q3g, q3gmin, q3gmax, qmin, qmax, rdenom, rsame
   real :: ectmass_2!,ts, wmin, wmax

   real :: temw(4), temwp(4), temwpp(4)!, &
!      aaa(3), bbb(3), ccc(3), ddd(3), qqq(3), avkii(3)
!   real :: roffpart ! tatsch 6 jun 2011
   !real :: anik  ! RHV 8/06/2016
   real :: pwww2,pwww3, count2, count3 , depwww,fracdz3
   real :: theta_fc,q2h, q2hv(101),etr
   integer :: later
   !anik = 4.0 ! RHV 8/06/2016
   
    real :: z_soisno(nlayers+1),dz_soisno(nlayers+1),zi_soisno(nlayers+2), &
    t_soisno(nlayers+1),vol_liq(nlayers+1),vol_ice(nlayers+1),icefrac(nlayers+1),&
    eff_porosity(nlayers+1),porsl(nlayers+1),hksati(nlayers+1),bsw(nlayers+1), &
    psi0(nlayers+1),rootr(nlayers+1),dwat(nlayers+1),qcharge
!c-----------------------------------------------------------------------
!c   OLD 
!c    calculation of gravitationally driven drainage from w(1-3) : taken
!c    as an integral of time varying conductivity.addition of liston
!c    baseflow term to original q3g to insure flow in
!c    dry season. modified liston baseflow constant scaled
!c    by available water.
!c		pows = 2.0 * bee + 2.0  originalmente
!c-----------------------------------------------------------------------
!   

!   if(iter .lt. 10) print*,'before run2',(dzwww(i),i=1,10) 
    q2hv = 0.0
   pows = 2.0 * bee + 3.0  ! RHV

!c-----------------------------------------------------------------------
!  Calculando o escoamento subsuperficial para cada uma das 
!  camadas do solo, implementando o modelo multicamada para
!  as camadas 2da e 3ra originales. 
!  RHV 11/06/2016
!c-----------------------------------------------------------------------

    
    i = 1
    temw(i) = max(0.03, www(i))                             ! teta/tetas
    temwp(i) = temw(i) ** (-bee)                            ! psi
    temwpp(i) = min(1.0, temw(i)) ** pows      ! k 
  later = 0
  if (later .eq. 1) then                                                       ! RUNOFF Horizontal
 
      theta_fc = poros*(abs(phsat)/100/340)**(1/bee)
        !theta_fc = min(theta_fc,0.99*poros )
        !theta_fc = max(0.01, theta_fc)

      Do i = 1, nlayers-1
    
          if(i.eq.1)then
              if(www(i)*zdepth(i)*poros .gt. theta_fc * zdepth(i))then
                  q2h = satcoz(i)*temwpp(i)*slope*dtt*0.3
                  q2h = min(q2h, www(i) * zdepth(i) * poros  -  theta_fc * zdepth(i))
                  roff = roff + q2h
                  roff2 = roff2 + q2h
                  www(i) = www(i) - q2h / (poros * zdepth(i))
                  q2hv(i) = q2h
              endif
          endif
      
          if( dzwww(i) * dzmultilayer * poros .gt. theta_fc * dzmultilayer )then
              q2h = satcoz_d(i)*dzwww(i)*pows*slope*dtt*0.3
              q2h = min(q2h, dzwww(i) * dzmultilayer * poros  -  theta_fc * dzmultilayer)
              roff = roff + q2h
              roff2 = roff2 + q2h
              dzwww(i) = dzwww(i) - q2h / (poros * dzmultilayer)
              q2hv(i+1) = q2h
          endif
    
      ENDDO
 
  endif
!     write(9856,'(a8,i10, 50(f13.6))')'fiel_cap',iter,&
!     (theta_fc - www(1)),& !Acima da Capacidade de campo WWW(1)
!     ( (theta_fc - dzwww(i)),i=1,nlayers), & !Acima da Capacidade de campo DZWWW(1:nlayers)
!     theta_fc ! Capacidade de campo
!     write(9856,'(a8,i10, 50(f13.6))')'roff2mm ',iter,(q2hv(i)*1000,i=1,nlayers+1), sum(q2hv)*1000
!     write(9856,'(a8,i10, 50(f13.6))')'ROFF2-H ',iter,www(1), & 
!     (dzwww(i),i=1,nlayers)  ! WRITING DZWWW depois de ROFF2
! !    write(9856,*) '--'                                                     !
!       write(1221,*) "roff2 ",roff2*1000  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!c----------------------------------------------------------------------
!c
!c    calculation of inter-layer exchanges of water due to gravitation
!c    and hydraulic gradient. the values of w(x) + dw(x) are used to
!c    calculate the potential gradients between layers.
!c    modified calculation of mean conductivities follows ME-82 ),
!c    reduces recharge flux to top layer.
!c
!c      dpdw           : estimated derivative of soil moisture potential
!c                       with respect to soil wetness. assumption of
!c                       gravitational drainage used to estimate likely
!c                       minimum wetness over the time step.
!c
!c      qqq  (q     )  : equation (61) , SE-86
!c             i,i+1
!c            -
!c      avk  (k     )  : equation (4.14) , ME-82
!c             i,i+1
!c
!c----------------------------------------------------------------------

	!call run2n
	!call run2clm
	
	! call run2Tang
   DO i=1,nlayers+1
     if(i.eq.1)then
     ! Variáveis ou elementos só da 1ra camada
       z_soisno(i) = dsfc_cst
       dz_soisno(i) = zdepth(i)
      zi_soisno(i) = 0.0
       zi_soisno(i+1) = z_soisno(i)
       vol_liq(i) = www(i)*poros!*zdepth(1)
       hksati(i) = satcoz(i)*1000          ! mm H2O s-1
       rootr(i) = 0.0
     endif
    ! Variáveis ou elementos comuns para todas as camadas 
     vol_ice(i) = 0.0
     icefrac(i) = 0.0
     eff_porosity(i) = poros
     porsl(i) = poros
     bsw(i) = bee
    t_soisno(i) = 290.0
     zi_soisno(i+1) = z_soisno(i)
     psi0(i) = phsat * 1000                  ! mm
     
     if(i.ge.2)then
    ! Variáveis ou elementos para o resto do perfil
       hksati(i) = satcoz_d(i-1)*1000       ! mm H2O s-1
       rootr(i) = extfrac(i-1)
       z_soisno(i) = dsfc_cst + (i-1)*dzmultilayer
       dz_soisno(i) = dzmultilayer
      vol_liq(i) = dzwww(i-1)*poros!*dzmultilayer
       
     endif
    rootr(nlayers+1) = 0.0
  ENDDO
 dwat = 0.0
 qcharge = 0.0
! 	if(iter .lt. 10) print*,'WWW1-3',www(1),www(2),www(3)
! 	if(iter .lt. 10) print*,'FIRST',(vol_liq(i),i=1,11)
  	etr=ectmass_2
!         write(234,*)iter, etr, ectmass_2
!   	call soilwater(nlayers+1,dtt,0.05,-1.e8,&
!                          0.0,etr,z_soisno,dz_soisno,zi_soisno,&
!                          t_soisno,vol_liq,vol_ice,icefrac,eff_porosity,&
!                          porsl,hksati,bsw,psi0,rootr,&
!                          gwdep,dwat,qcharge)
  	call soilwater_sib3(nlayers+1,hlat,vol_liq,hksati/1000,poros,phsat, &
  			  eff_porosity,bee,td,etr,                            &
  			  0.0,rootr, dz_soisno*1000, z_soisno*1000,  &
 			  qcharge, dwat)
                        
!      write(9867,'(a8,i10, 50(1x,f9.4))')"CLM   ",iter,&
!      www(1) + dwat(1)/poros,(dzwww(i)+dwat(i+1)/poros,i=1,nlayers) ,qcharge*3600             
                       

!c-----------------------------------------------------------------------
!c     update wetness of each soil moisture layer due to layer interflow
!c        and base flow.
!c-----------------------------------------------------------------------
   do i = 1,nlayers
      if(i.eq.1)then
	www(i)  = www(i) + dwat(i)/poros
      endif
        dzwww(i) = dzwww(i) + dwat(i+1)/poros
   end do
      roff = roff + qcharge
      roff3 = roff3 + qcharge
      
                             
!      write(9856,'(a8,i10, 50(f13.6))')"INTERLA ",iter,www(1),&
!      (dzwww(i),i=1,nlayers) , roff3*1000! WRITING DZWWW depois de INTERLAYER
!    write(9856,*) '--'
      
      excess = max(0.0, www(1) - 1.0)
      www(1) = www(1) - excess
      roff   = roff + excess * poros * zdepth(1)
      roff4  = roff4 + excess * poros * zdepth(1)
      
   do i = 1,nlayers
      excess   = max(0.0, dzwww(i) - 1.0)
      dzwww(i) = dzwww(i) - excess
      roff     = roff + excess * poros * dzmultilayer
      roff4    = roff4 + excess * poros * dzmultilayer
   end do
      !write(1221,*) "roff4 ",roff4*1000 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!c-----------------------------------------------------------------------
!c     prevent negative values of www(i)
!c-----------------------------------------------------------------------
      deficit = max(0.0, 1.e-12 - www(1))
      www(1) = www(1) + deficit
      dzwww(1) = dzwww(1) - deficit
   do i = 1,nlayers-1
      deficit = max(0.0, 1.e-12 - dzwww(i))
      dzwww(i) = dzwww(i) + deficit
      dzwww(i+1) = dzwww(i+1) - deficit
   end do
      dzwww(nlayers) = max(dzwww(nlayers), 1.e-12)
      
         ! Retornando a modelo de três camadas para conferir resultados
      depwww = dsfc_cst ! Retornando a profundidade da primeira camada
      count2 = 0.0
      count3 = 0.0
      pwww2  = 0.0
      pwww3  = 0.0
	  ! real :: pwww2,pwww3, count2, count3
      DO i=1,nlayers
	depwww = depwww + dzmultilayer ! Profundidade da base da subcamada ...
	  if (depwww .le. (zdepth(1) + zdepth(2))) then
	    count2 = count2 + 1.0
	    ! Incremento no grau de saturacao da 3ra camada
	    pwww2 = pwww2 + dzwww(i)
	  else if (depwww .gt. (zdepth(1) + zdepth(2))  &
	                  .and.                         &
	           depwww .lt. (zdepth(1) + zdepth(2) + dzmultilayer)) then
	      ! Fracao de dzcamada que tem parte dela na 2da camada e parte na 3ra 
	      fracdz3 = (depwww - (zdepth(1)+zdepth(2)))/dzmultilayer
	      ! Incremento ao numero de subcamadas da segunda
	      count2 = count2 + (1 - fracdz3)
	      pwww2  = pwww2  + (1 - fracdz3) * dzwww(i)
	      ! Incremento ao numero de subcamadas na 3ra camada
	      count3 = count3 + fracdz3
	      ! Incremento no grau de saturacao da 3ra camada
	      pwww3  = pwww3  + fracdz3 * dzwww(i)
	  else
	    ! Contador das camadas para determinar o grau de saturacao médio 
	    count3 = count3 + 1.0
	    ! Incremento no grau de saturacao da 3ra camada
	    pwww3  = pwww3  + dzwww(i)
	  end if
      ENDDO
      ! Determinando grau de saturação médio de cada camada
       www(2) = pwww2 / count2
       www(3) = pwww3 / count3 
   
  
	  
end subroutine run2

!c====================================================================
!   Subroutine from CLM model for resolve interlayer flux
!   June/15/2016 01:24AM
!   Roilan Hernandez
!c--------------------------------------------------------------------
!	subroutine run2Tang
!
!!c----------------------------------------------------------------------
!!c
!!c    calculation of inter-layer exchanges of water due to gravitation
!!c    and hydraulic gradient. the values of w(x) + dw(x) are used to
!!c    calculate the potential gradients between layers.
!!c    modified calculation of mean conductivities follows ME-82 ),
!!c    reduces recharge flux to top layer.
!!c
!!c      dpdw           : estimated derivative of soil moisture potential
!!c                       with respect to soil wetness. assumption of
!!c                       gravitational drainage used to estimate likely
!!c                       minimum wetness over the time step.
!!c
!!c      qqq  (q     )  : equation (61) , SE-86
!!c             i,i+1
!!c            -
!!c      avk  (k     )  : equation (4.14) , ME-82
!!c             i,i+1
!!c
!!c----------------------------------------------------------------------
!      use steps, only : dtt
!      use stepv, only : www, td, tg
!      use temp, only : gwdep
!      use multicamada
!      use soils
!      use snow, only : tsnow
!      use grads, only : tgs
!      use const, only : tf
!      use donor, only : roff
!      use hydrol, only : areas
!      use output, only : roff2, roff3, roff4, roffq3g, roffq3g2
!      use morphology1
!      implicit NoNE
!
!      ! Local variables
!      integer :: i
!      real :: avk, avkmin, avkmax, d3g, deficit, denom, div
!      real :: dpdw, dpdwdz, excess, pmin, pmax, pows, props
!      real :: q3g, q3gmin, q3gmax, qmin, qmax, rdenom, rsame
!      real :: ts, wmin, wmax
!
!      real :: temw(4), temwp(4), temwpp(4), &
!	  aaa(3), bbb(3), ccc(3), ddd(3), qqq(3), avkii(3)
!      real :: pwww2,pwww3, count2, count3,fracdz3,depwww
!
!      pows = 2.*bee + 3.
!          ! Retornando a modelo de três camadas para conferir resultados
!      depwww = dsfc_cst ! Retornando a profundidade da primeira camada
!      count2 = 0.0
!      count3 = 0.0
!      pwww2  = 0.0
!      pwww3  = 0.0
!	  ! real :: pwww2,pwww3, count2, count3
!      DO i=1,nlayers
!	depwww = depwww + dzmultilayer ! Profundidade da base da subcamada ...
!	  if (depwww .le. (zdepth(1) + zdepth(2))) then
!	    count2 = count2 + 1.0
!	    ! Incremento no grau de saturacao da 3ra camada
!	    pwww2 = pwww2 + dzwww(i)
!	  else if (depwww .gt. (zdepth(1) + zdepth(2))  &
!	                  .and.                         &
!	           depwww .lt. (zdepth(1) + zdepth(2) + dzmultilayer)) then
!	      ! Fracao de dzcamada que tem parte dela na 2da camada e parte na 3ra
!	      fracdz3 = (depwww - (zdepth(1)+zdepth(2)))/dzmultilayer
!	      ! Incremento ao numero de subcamadas da segunda
!	      count2 = count2 + (1 - fracdz3)
!	      pwww2  = pwww2  + (1 - fracdz3) * dzwww(i)
!	      ! Incremento ao numero de subcamadas na 3ra camada
!	      count3 = count3 + fracdz3
!	      ! Incremento no grau de saturacao da 3ra camada
!	      pwww3  = pwww3  + fracdz3 * dzwww(i)
!	  else
!	    ! Contador das camadas para determinar o grau de saturacao médio
!	    count3 = count3 + 1.0
!	    ! Incremento no grau de saturacao da 3ra camada
!	    pwww3  = pwww3  + dzwww(i)
!	  end if
!      ENDDO
!      ! Determinando grau de saturação médio de cada camada
!       www(2) = pwww2 / count2
!       www(3) = pwww3 / count3
!
!    DO i=1,3
!      temw(i) = max(0.03, www(i))                             ! teta/tetas
!      temwp(i) = temw(i) ** (-bee)                            ! psi
!      temwpp(i) = min(1.0, temw(i)) ** pows                   ! k
!    END DO
!
!   wmax = max(www(1), www(2), www(3), 0.05)
!   wmax = min(wmax, 1.0)
!   pmax = wmax ** (-bee)
!   wmin =(pmax - (zdepth(1) + zdepth(2) + zdepth(3)) / phsat) ** (-1.0 / bee)
!
!   wmin = min(www(1), www(2), www(3), wmin)
!   wmin = max(wmin, 0.02)
!   pmin = wmin ** (-bee)
!   dpdw = phsat * (pmax - pmin) / (wmax - wmin)
!!c	print *,dpdw,wmax,wmin
!
!   do i = 1, 2
!      rsame = 0.0
!      avk = temwp(i) * temwpp(i) - temwp(i+1) * temwpp(i+1)
!      div = temwp(i+1) - temwp(i)
!      if (abs(div) < 1.e-6) rsame = 1.0
!!c!        avk = satco*avk / ( ( 1. + 3./bee ) * div + rsame )
!      avk = satcoz(i) * avk / ((1.0 + 3.0 / bee) * div + rsame)
!!c!        avkmin = satco * min( temwpp(i), temwpp(i+1) )
!      avkmin = satcoz(i) * min(temwpp(i), temwpp(i+1))
!!c!        avkmax = satco * max( temwpp(i), temwpp(i+1) )*1.01
!      avkmax = satcoz(i) * max(temwpp(i), temwpp(i+1)) * 1.01
!      avk = max(avk, avkmin)
!      avk = min(avk, avkmax)
!!c 3/96 changes
!      if (www(i) < www(i+1)) avk = 0.1 * avk
!
!!c-----------------------------------------------------------------------
!!c     conductivities and base flow reduced when temperature drops below
!!c     freezing.
!!c-----------------------------------------------------------------------
!
!      tsnow = min(tf - 0.01, tg)
!      tgs = tsnow * areas + tg * (1.0 - areas)
!      ts = tgs * real(2 - i) + td * real(i - 1)
!      props = (ts - (tf - 10.0)) / 10.0
!      props = max(0.05, min(1.0, props))
!      avk = avk * props
!
!!c-----------------------------------------------------------------------
!!c     backward implicit calculation of flows between soil layers.
!!c-----------------------------------------------------------------------
!
!      dpdwdz = dpdw * 2.0 / (zdepth(i) + zdepth(i+1))
!      aaa(i) = 1.0 + avk * dpdwdz * (1.0 / zdepth(i) + 1.0 / zdepth(i+1)) &
!         * dtt / poros
!      bbb(i) = -avk * dpdwdz * 1.0 / zdepth(2) * dtt / poros
!      ddd(i) = avk * (dpdwdz * (www(i) - www(i+1)) + 1.0)
!      if (i == 2) then
!         ccc(i) = -avk * dpdwdz * 1.0 / zdepth(3) * dtt / poros
!      endif
!      avkii(i) = avk
!   end do
!   ccc(1) = 0.0
!   ccc(3) = 0.0
!
!   d3g = gwdep - (zdepth(1) + zdepth(2) + zdepth(3) / 2.0)
!   d3g = max(d3g, zdepth(3) / 2.0, 0.5)
!
!   avk = satcoz(3) * temw(3) ** (2.0 * bee + 3.0)                  !tatsch (i)
!
!   avkii(3) = avk
!
!   wmax = 1.0
!   pmax = 1.0
!   wmin = (pmax - d3g / phsat) ** (-1.0 / bee)
!   wmin = min( www(3), wmin)
!   wmin = max(wmin, 0.02)
!   pmin = wmin ** (-bee)
!   dpdw = phsat * (pmax - pmin) / (wmax - wmin)
!
!   dpdwdz = dpdw / d3g
!   aaa(3) = 1.0 + avk * dpdwdz * (1.0 / zdepth(3)) * dtt / poros
!   bbb(3) = -avk * dpdwdz * 1.0 / zdepth(3) * dtt / poros
!   ddd(3) = avk * (dpdwdz * www(3) + 1.0 - dpdwdz)
!
!   denom = (aaa(1) * aaa(2) * aaa(3) - aaa(3) * bbb(1) * bbb(2) &
!      - aaa(1) * bbb(3) * ccc(2))
!   rdenom = 0.0
!   if (abs(denom) < 1.e-6) rdenom = 1.0
!   rdenom = (1.0 - rdenom) / (denom + rdenom)
!   qqq(1) = (aaa(2) * ddd(1) * aaa(3) - bbb(3) * ccc(2) * ddd(1) &
!      - aaa(3) * ddd(2) * bbb(1) + ccc(2) * ddd(3) * bbb(1)) * rdenom
!   qqq(2) = (aaa(1) * aaa(3) * ddd(2) - aaa(3) * ddd(1) * bbb(2) &
!      - aaa(1) * ccc(2) * ddd(3)) * rdenom
!   q3g = (aaa(1) * aaa(2) * ddd(3) - bbb(1) * bbb(2) * ddd(3) &
!      - aaa(1) * bbb(3) * ddd(2) + bbb(3) * ddd(1) * bbb(2) ) * rdenom
!
!   q3gmax = www(3) * (poros * zdepth(3) / dtt)
!   q3gmin = (www(3) - 1.0) * (poros * zdepth(3) / dtt)
!   q3gmin = 0.0
!   q3g = min(q3g, q3gmax)  !, avkii(3)
!   q3g = max(q3g, q3gmin)	!, -avkii(3)
!
!!c-----------------------------------------------------------------------
!!c     update wetness of each soil moisture layer due to layer interflow
!!c        and base flow.
!!c-----------------------------------------------------------------------
!
!   www(3) = www(3) - q3g * dtt / (poros * zdepth(3))
!!c	gwsoil	= gwsoil-q3gn*dtt
!
!   do i = 1, 2
!      qmax = www(i) * (poros * zdepth(i) / dtt)
!      qmin = -www(i+1) * (poros * zdepth(i+1) / dtt)
!      qqq(i) = min(qqq(i), qmax)         !, avkii(i)
!      qqq(i) = max(qqq(i), qmin)         !, -avkii(i)
!      www(i) = www(i) - qqq(i) / (poros * zdepth(i) / dtt)
!      www(i+1) = www(i+1) + qqq(i) / (poros * zdepth(i+1) / dtt)
!   end do
!
!   roff = roff + q3g * dtt
!   roffq3g2 = roffq3g2 + q3g * dtt
!   roff3 = roff3 + q3g* dtt
!
!	return
!	end subroutine run2Tang
!
!!c====================================================================
!!   Subroutine from CLM model for resolve interlayer flux
!!   June/15/2016 01:24AM
!!   Roilan Hernandez
!!c--------------------------------------------------------------------
!	subroutine run2clm
!
!	use multicamada
!	use soils
!	use morphology1
!	use steps
!	use stepv
!	use output
!	use donor
!	implicit NoNE
!
!	real :: pows
!	real :: temw(1), temwp(1), temwpp(1)
!! 	real :: temw_d(nlayers),temwp_d(nlayers),temwpp_d(nlayers)
!	integer :: i,j,k
!	real :: smp0, dsmpdw0,hk0,alpha0,dhkdw10,dhkdw20,dqodw10,dqodw20
!	real :: qin0,qout0
!	real :: smp(nlayers+1), dsmpdw(nlayers),alpha(nlayers), &
!		hk(nlayers),dhkdw1(nlayers),dhkdw2(nlayers)
!	real :: qin(nlayers+1),qout(nlayers+1)
!	real :: amx(nlayers+1),bmx(nlayers+1),cmx(nlayers+1),rmx(nlayers+1),dwat(nlayers+1)
!	real:: dqidw0(nlayers) ,dqidw1(nlayers) ,dqodw1(nlayers) ,dqodw2(nlayers)
! 	real :: q3g,errorw
!	real :: phsatmm,dzmultilayermm,dsfc_cstmm,satcoz_dmm(100),pows_d
!	real :: smpmin
!
!	pows= 2. * bee + 3.
!	pows_d = 2. * bee + 2.
!	smpmin = -1.e8
!
!    i = 1
!    temw(i) = min(1.0,max(0.01, www(i)) )                   ! teta/tetas
!    temwp(i) = temw(i) ** (-bee)                            ! psi
!    temwpp(i) = min(1.0, temw(i)) ** pows                   ! k
!
!    ! Compute matric potential and derivative based on liquid water content only
!    smp0 = max(phsat*temwp(i),smpmin)
!    dsmpdw0 = -bee * smp0 / (www(i)*poros)
!    do i=1, nlayers
!	temw_d(i) = max(0.01, dzwww(i))                     ! teta/tetas
!	temwp_d(i) = temw_d(i) ** (-bee)                    ! psi
!	temwpp_d(i) = min(1.0, temw_d(i)) ** pows           ! k
!	  smp(i) = max(phsat * temwp_d(i),smpmin)
!	  dsmpdw(i) = -bee * smp(i) / (dzwww(i)*poros)
!    enddo
!	smp(nlayers+1) = phsat! Aqui tô sopondo que a fronteira inferior...
!	! ... é o aquífero e portanto, está saturado assim o ...
!	! ... smp == phsat (www(Aquifero) == 1.0)
!	!print*,phsat,(smp(k),k=1,nlayers+1)
!
!! Condutividade hidráulica real, potencial matricial e as derivadas....
!    j = 1
!    qin0 = 0.0 ! Infiltração foi resolvida pelo GA
!    ! Para a primeira camada:
!
!    if(www(1).lt.1.e-3)then
!      hk0 = 0.
!      dhkdw20 = 0.0
!      dhkdw10 = 0.0
!    else
!      hk0 = satcoz(j) * temwpp(j)
!      alpha0 = (smp(j)-smp0)/dsfc_cst -1.0
!	if(alpha0 .le. 0.0)then
!	  dhkdw10 = satcoz(j) * pows * temw(j)**pows_d/ poros
!	  dhkdw20 = 0.0
!	else
!	  dhkdw10 = 0.0
!	  dhkdw20 = satcoz_d(j)* pows * temw_d(j)**pows_d / poros
!	endif
!    endif
!
!    do i=1,nlayers-1
!
!      alpha(i) = (smp(i+1)-smp(i))/dzmultilayer -1.0
!
!      if(dzwww(i).lt.1.e-3)then
!	hk(i) = 0.
!	dhkdw1(i) = 0.
!	dhkdw2(i) = 0.0
!      else
!
!	if(alpha(i) .le. 0.0)then
!	  hk(i) = satcoz_d(i) * temwpp_d(i)
!	  dhkdw1(i) = satcoz_d(i) * pows * temw_d(i)**pows_d / poros
!	  dhkdw2(i) = 0.0
!	else
!	  hk(i) = satcoz_d(i+1) * temwpp_d(i+1)
!	  dhkdw1(i) = 0.0
!	  dhkdw2(i) = satcoz_d(i+1) * pows * temw_d(i+1)**pows_d /poros
!	endif
!
!      endif
!
!    enddo
!
!	i = nlayers
!	if(dzwww(i).lt.1.e-3)then
!	  hk(i)     = 0.
!	  dhkdw1(i) = 0.
!	  dhkdw2(i) = 0.0
!	else
!	  hk(i)     = satcoz_d(i) * temwpp_d(i)
!	  dhkdw1(i) = satcoz_d(i) * pows * temw_d(i)**pows_d /poros
!	  dhkdw2(i) = 0.0
!	endif
!! .... Fim de determinação das K_theta e Psi_theta
!
!! Fluxos de entrada e saída
!!   ...para primeira camada:
!    qout0 =  hk0 * alpha0
!    dqodw10 = -(alpha0 * dhkdw10) - hk0*dsmpdw0/dsfc_cst
!    dqodw20 = -(alpha0 * dhkdw20) + hk0*dsmpdw(j)/dsfc_cst
!
!!   ...elementos da matriz para a 1ra camada
!    amx(j) = 0.0
!    bmx(j) = dsfc_cst/dtt + dqodw10
!    cmx(j) = dqodw20
!    rmx(j) = qin0 - qout0  !- etr * rootr(j) Aqui o CLM determina a transpiração
!
!
!  ! Fluxos verticais e elementos da matriz para as camadas discretizadas:
!    do i= 1,nlayers-1
!       j=i+1
!      if(i.eq.1)then
!       qin(i) = -hk0*alpha0
!       dqidw0(i) = -(alpha0*dhkdw10 - hk0*dsmpdw0/dsfc_cst)
!       dqidw1(i) = -(alpha0*dhkdw20 + hk0*dsmpdw(i)/dsfc_cst)
!      else
!       qin(i) = -hk(i-1)*alpha(i-1)
!       dqidw0(i) = -(alpha(i-1)*dhkdw1(i-1) - hk(i-1)*dsmpdw(i-1)/dzmultilayer)
!       dqidw1(i) = -(alpha(i-1)*dhkdw2(i-1) + hk(i-1)*dsmpdw(i)/dzmultilayer)
!      endif
!
!       qout(i) = -hk(i)*alpha(i)
!       dqodw1(i) = -(alpha(i)*dhkdw1(i) - hk(i)*dsmpdw(i)/dzmultilayer)
!       dqodw2(i) = -(alpha(i)*dhkdw2(i) + hk(i)*dsmpdw(i+1)/dzmultilayer)
!
!       amx(j) = -dqidw0(i)
!       bmx(j) =  dzmultilayer/dtt - dqidw1(i) + dqodw1(i)
!       cmx(j) =  dqodw2(i)
!       rmx(j) =  qin(i) - qout(i) ! - etr*rootr(j) Aqui o CLM determina a transpiração
!    enddo
!
!     ! Node j=nl_soil (bottom)
!
!      i = nlayers
!    qin(i) = -hk(i-1)*alpha(i-1)
!    dqidw0(i) = -(alpha(i-1)*dhkdw1(i-1) - hk(i-1)*dsmpdw(i-1)/dzmultilayer)
!    dqidw1(i) = -(alpha(i-1)*dhkdw2(i-1) + hk(i-1)*dsmpdw(i)/dzmultilayer)
!
!    qout(i) = hk(i)
!    dqodw1(i) = dhkdw1(i)
!    dqodw2(i) = 0.
!       j = i+1
!    amx(j) = -dqidw0(i)
!    bmx(j) =  dzmultilayer/dtt - dqidw1(i) + dqodw1(i)
!    cmx(j) =  dqodw2(i)
!    rmx(j) =  qin(i) - qout(i) ! - etr*rootr(j) Aqui o CLM determina a transpiração
!
!    ! Solve for dwat
!
!    call tridia (nlayers+1, amx, bmx, cmx, rmx, dwat )
!	write(1865, '(i12,50(2x,e13.6))')iter, (dwat(i),i=1,nlayers+1)
!
!	DO i = 1, nlayers+1
!	  j=i-1
!	  if(i.eq.1)then
!	    www(i)  = www(1) + dwat(i)/zdepth(i)/poros
!	  else
!	    dzwww(j)  = dzwww(j) + dwat(j)/dzmultilayer/poros
!	  endif
!	ENDDO
!
!
!	errorw = -dtt*(qin(1)-qout(nlayers)-dqodw1(nlayers)*dwat(nlayers))
!	      errorw = errorw - dwat(1)*dsfc_cst
!	    do j = 1, nlayers
!	      errorw = errorw+dwat(j)*dzmultilayer
!	    enddo
!
!    if(abs(errorw) > 1.e-3)then
!       write(*,*) iter,'mass balance error in time step =',errorw
!    endif
!	! Adicionando valores aos escoamentos
!	q3g = (qout(nlayers) + dqodw1(nlayers)*dwat(nlayers))
!	roff= roff +  q3g*dtt
!	roffq3g2 = roffq3g2 + q3g*dtt
!	roff3 = roff3 + q3g*dtt
!
!
!
!        end subroutine run2clm
!  
!!c====================================================================
!!c
!        subroutine run2n
!!c
!!c     modified multi-layer scheme:  Humberto da Rocha 03/04/1996
!!c	modified for compatibility with SiB2 do DBHM
!!c         Roilan Hernandez 12/06/2016
!!c====================================================================
!    use multicamada
!    use soils
!    use output
!    use donor
!    use steps
!    use stepv
!    implicit NoNE
!
!    real :: qdowrd(nlayers), qupwrd(nlayers), qqq(nlayers)
!    real :: a(nlayers), b(nlayers),c(nlayers), d(nlayers), u(nlayers)
!    real :: qmin, qmax, pows, qng
!    real :: avk, dpdw , avdif,dpsidz, avdiff(nlayers)
!    integer :: i,j,jesq
!    real :: temw(nlayers), temwp(nlayers), temwpp(nlayers)
!
!      qdowrd = 0.	! fluxos downward (drenagem rapida + difusiva), i=1 (infiltracao superficie)
!      qupwrd = 0.	! fluxos upward (capilar)
!
!      avdiff = 0.0
!      avdif = 0.0
!
!      DO i=1,nlayers-1
!
!	  if( i .eq. 1) qqq(i)=0.0
!
!	call retec  (dpdw, 0, i)
!	call hydcon (avk, i)
!
!	  avdif = avk * dpdw / poros
!
!	  avdiff(i) = avdif   !temporal
!
!! 	if(i .eq. 1)then
!! 	  temw(i) = max(0.03, www(i))                             ! teta/tetas
!! 	  temwp(i) = temw(1) ** (-bee)                            ! psi
!! 	  temwpp(i) = min(1.0, temw(i)) ** (2.0 * bee + 3.0)      ! k
!! 	  dpsidz = phsat * temwp(i)
!! 	  dpsidz = dpsidz - phsat * dzwww(i)
!! 	  dpsidz = 2.* dpsidz / (zdepth(i) + dzmultilayer)
!! 	  qqq(i) = -avk * (dpsidz + 1.)
!! 	else
!! 	  temw(i-1) = max(0.03, dzwww(i-1))                           ! teta/tetas
!! 	  temwp(i-1) = temw(i-1) ** (-bee)                            ! psi
!! 	  temwpp(i-1) = min(1.0, temw(i-1)) ** (2.0 * bee + 3.0)      ! k
!! 	  dpsidz = phsat * temwp(i-1)
!! 	  dpsidz = dpsidz - phsat * dzwww(i)
!! 	  dpsidz = 2.* dpsidz / (2 * dzmultilayer)
!! 	  qqq(i) = -avk * (dpsidz + 1.)
!! 	endif
!! 	jesq = 1
!
!! 	if(jesq.eq.1)then
!!	  a(i) = -2.*dtt*avdif / ( 3 * dzmultilayer )
!!	  b(i) =  2.*dtt*avdif / ( 2 * dzmultilayer ) + 1.
!!	  c(i) = -2.*dtt*avdif / ( 3 * dzmultilayer )
!!	  d(i) = qqq(i)
!!	  if (i.eq.1) 	      a(i) = 0.
!!	  if (i.eq.nlayers-1) c(i) = 0.
!! 	endif
!
!!  	if(jesq.eq.2)then
!!  	  a(i) = -dtt*avdif / (3 * dzmultilayer )
!!  	  b(i) =  1. + dtt*avdif / ( 2 * dzmultilayer )
!!  	  c(i) = -dtt*avdif / ( 3 * dzmultilayer )
!!  	  if (i.eq.1) 	   	a(i) = 0.
!!  	  if (i.eq.nlayers-1) 	c(i) = 0.
!!  	  da(i) = -a(i)
!!  	  db(i) = 1. - dtt*avdif / ( 2 *dzmultilayer )
!!            dc(i) = -c(i)
!!           if (i.ge.2.and.i.le.(nlayers-2)) d(i) = da(i)*qqq(i-1) + db(i)*qqq(i) + dc(i)*qqq(i+1)
!!
!!  	if (i.eq.1)
!!       &  d(i) =  db(i)*qqq(i) + dc(i)*qqq(i+1)
!!
!! 	if (i.eq.(nlayer-1))
!!       &  d(i) = da(i)*qqq(i-1) + db(i)*qqq(i)
!!  	endif
!
!
!      ENDDO
!
!	!call tridia2 (nlayers-1, d, c, a, b )
! 	 call  tridia (nlayers-1, a, b, c, d, u)
!
!
!	 ! i
!!       write(1635,*) iter, qqq(1), (u(i),i=1,nlayers) !,qqq(1)*1000, (u(j),j=1,nlayers-1)
!     !'(i12,50(1x,f13.10))'
!     !write(1221,*),"BEFORE",dzwww(nlayers) , d(nlayers)
!      	DO i=1,nlayers-1
!
!	   ! qqq(i)= d(i)  ! use tridia2
! 	     qqq(i)= u(i) ! use tridia
!
!	    if (qqq(i).lt.0.) qdowrd(i+1) = qdowrd(i+1) + abs(qqq(i))
!	    if (qqq(i).gt.0.) qupwrd(i+1) = qupwrd(i+1) + qqq(i)
!
!	    qmin=  -dzwww(i)   * (poros  * dzmultilayer / dtt)
!	    qmax=   dzwww(i+1) * (poros  * dzmultilayer / dtt)
!	    qqq(i)= amin1( qqq(i), qmax )
!	    qqq(i)= amax1( qqq(i), qmin )
!	    dzwww(i)  = dzwww(i)   + qqq(i) / (poros * dzmultilayer / dtt)
!	    dzwww(i+1)= dzwww(i+1) - qqq(i) / (poros * dzmultilayer / dtt)
!
!	ENDDO
!  	
!	pows= 2. * bee + 2.
!	qng= temwpp_d(nlayers) + satcoz_d(nlayers)/dzmultilayer/poros*slope*pows*dtt
!	qng= qng ** (1./pows)
!	qng= -( 1./qng-dzwww(nlayers) )* poros * dzmultilayer / dtt
!	qng= amax1( 0., qng)
!	qng= amin1( qng, dzwww(nlayers)*poros * dzmultilayer / dtt)
!	dzwww(nlayers)= dzwww(nlayers)- qng*dtt/(poros*dzmultilayer)
!
!!     write(9856,'(a8,i10, 50(f13.6))')'FLUX_IO ',iter,qqq(1)*3.6e6,(d(j)*3.6e6,j=1,nlayers-1),qng*3.6e6     ! WRITING DZWWW antes de INTERLAYER
!!        write(9856,*) '--'
!	!write(1221,*)"AFTER",dzwww(nlayers),"FLUXO OUT: ",qng*dtt*1000
!
!	! Adicionando valores aos escoamentos
!	roff= roff + qng*dtt
!	roffq3g2 = roffq3g2 + qng * dtt
!	roff3 = roff3 + qng* dtt
!! 	write(1221,*) "roff3 ",roff3*1000 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!    return
!    end
!
!!c=======================================================================
!	subroutine retec (dpdw, jhort, i)
!!c=======================================================================
!!c
!!c..	output: cdpdw (derivada curva retencao) dpsi / dw
!!c
!!c----------------------------------------------------------------------
!	use multicamada
!	use soils
!	use stepv, only : www
!!c...
!      integer :: j,i,jdpsi
!      real :: w0,w1,w2,psi1,psi2
!      real :: wmax,pmax,pmin, wmin,dpdw
!      j = i + 1
!
!      jdpsi = 4
!
!      if (i.gt.1) then
!	w0 = amin1(www(1),1.)
!	w0 = amax1(www(1),0.05)
!      endif
!
!	w1 = amin1(dzwww(i),1.)
!	w1 = amax1(dzwww(i),0.05)
!	w2 = amin1(dzwww(j),1.)
!	w2 = amax1(dzwww(j),0.05)
!
!      psi1 = phsat * (w1**(-bee))
!      psi2 = phsat * (w2**(-bee))
!      if(i.eq.5) write(9898,*), phsat,bee,w1,w2,psi1,psi2
!      if(jdpsi.eq.1)then
!!c... calculo método original SiB2 SE-96
!	wmax= amax1 (dzwww(i), dzwww(i+1))
!	wmax= amax1 ( amin1(wmax,1.0),  0.05)
!	pmax= wmax**(-bee)
!
!	if(i.eq.1) then
!      wmin= ( pmax-  2./(phsat *(zdepth(1) + 3.0*dzmultilayer )) )
!!ch    wmin= (pmax-2./(phsat(1)*(2.*zdepth(1)+zdepth(2))))
!      wmin= amin1(w1,w2,wmin)
!      wmin= amax1(wmin,1.0e-02)
!	endif
!	if(i.gt.1) then
!      wmin=(pmax-2./(phsat*(4.*dzmultilayer)))
!      wmin= amin1(w0,w1,w2,wmin)
!      wmin= amax1(wmin,1.0e-02)
!	endif
!	if (wmin.eq.wmax) then
!      wmin= (pmax-1./(phsat*(4.*dzmultilayer)))
!      wmin= amax1(wmin,1.0e-03)
!	endif
!	pmin= wmin**(-bee)
!	dpdw= phsat * ( pmax-pmin )/( wmax-wmin )
!
!	return
!      endif
!
!      dpdw1 = - bee * phsat / (w1 ** (bee + 1.))
!      dpdw2 = - bee * phsat / (w2 ** (bee + 1.))
!
!      if (jdpsi.eq.2) then	! media ponderada espessura
!	if(i.eq.1) then
!	    dpdw = ( dpdw1  * zdepth(i) +  dpdw2  * dzmultilayer ) / &
!	           ( dzmultilayer + zdepth(i) )
!        else
!            dpdw = ( dpdw1  * dzmultilayer +  dpdw2  * dzmultilayer ) / &
!                   ( dzmultilayer + dzmultilayer )
!        endif
!      endif
!!c...
!      if (jdpsi.eq.3) then	! media geometrica
!      dpdw = sqrt( dpdw1 * dpdw2 )
!      endif
!!c
!      if (jdpsi.eq.4) then	! media aritmetica
!      dpdw = ( dpdw1 + dpdw2 )/2.
!      endif
!
!	return
!	end
!
!
!!c====================================================================
!	subroutine hydcon (avk, i)
!!c====================================================================
!!c
!!c     avk (condutividade hidraulica segmento camadas i,i+1)
!!c
!!c--------------------------------------------------------------------
!    use multicamada
!    use soils
!    use snow
!    use hydrol
!    use const
!    use stepv
!    use grads
!    implicit NoNE
!
!    real :: w1,w2,psi1,psi2,avk1,avk2,avb, avk, div,avkmin,avkmax
!    integer :: j,i
!    real :: props,rsame,qng
!    real :: ts
!
!    j = i + 1
!      w1 = amin1(dzwww(i),1.)
!      w1 = amax1(dzwww(i),0.05)
!      w2 = amin1(dzwww(j),1.)
!      w2 = amax1(dzwww(j),0.05)
!
!      psi1 = phsat * (w1**(-bee))
!      psi2 = phsat * (w2**(-bee))
!      avk1 = satcoz_d(i) * (w1**(2.*bee +3.))
!      avk2 = satcoz_d(i) * (w2**(2.*bee +3.))
!
!      rsame= 0.
!      avb = bee
!      div= psi2 - psi1
!      if ( abs(div) .lt. 1.e-6 ) rsame=1.
!      avk = (psi1*avk1 - psi2*avk2) / &
!           ( ( 1. +3./avb ) * div + rsame )
!
!      avkmin = amin1 (avk1, avk2)
!      avk    = amax1 (avk, avkmin)
!      avkmax = amax1 (avk1, avk2 )
!      avk    = amin1 (avk, avkmax )
!
!
!!c----------------------------------------------------------------
!!c      conductivites and baseflow reduced when temperature drops
!!c      below freezing
!!c----------------------------------------------------------------
!
!	tsnow = amin1( tf-0.01,tg)
!	tgs = tsnow * areas + tg* (1. - areas)
!	ts = tgs * (2 - j)  + td*(j-1)
!	props = ( ts - (tf - 10.) ) / 10.
!	props = amax1(0.05,amin1(1.0,props))
!	!avk = avk * props
!
!	return
!	end
!
!
!!c=======================================================================
!	subroutine tridia2 (n, d, c, a, b,u)
!!c=======================================================================
!!c
!!c  	HR  algoritmo original 1996, corrigido 2009
!!c
!!c             solution of linear system of n x n dimension
!!c             for a tridiagonal matrix
!!c             d    : independent coefficients (input)
!!c                         vector solution      (output)
!!c             b  : main diagonal coefs
!!c             c   :  suprior diagonal coefs
!!c             a   : inferior diagonal coefs
!!c--------------------------------------------------------------
!
!	dimension c(n+1), a(n+1), b(n+1), d(n+1),u(n+1)
!
!	c(1) = c(1)/b(1)
!	d(1) = d(1)/b(1)
!
!	do 10 i = 2, (n-1), 1
!	ii = (i-1)
!
!	b(i) = b(i) - c(ii)*a(i)
!
!	if (i .eq. n) go to 10
!	c(i) = c(i) / b(i)		! calculado
!   10   d(i) = (d(i) - d(ii)*a(i)) / b(i)
!
!	do 20 k = 1, (n-1), 1
!	i = n -k		! back substitution
!   20	d(i) = d(i) - c(i) * d(i+1)
!        u = d
!	return
!	end
!
!c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 	Subroutine from COLM
!	Solução de matriz tridiagonal
!c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine tridia (n, a, b, c, r, u)
      
      implicit none
      !integer, parameter :: r8 = selected_real_kind(12)
      integer :: n        !length of diagonal element vector
      real :: a(1:n)  !subdiagonal elements
      real :: b(1:n)  !diagonal elements
      real :: c(1:n)  !superdiagonal elements
      real :: r(1:n)  !right hand side
      real :: u(1:n) !solution vector

      integer :: j
      real ::  gam(1:n),bet
! -----------------------------------------------------------------

      bet = b(1)
      u(1) = r(1) / bet
      do j = 2, n
            gam(j) = c(j-1) / bet
            bet = b(j) - a(j) * gam(j)
            u(j) = (r(j) - a(j)*u(j-1)) / bet
      end do
      do j = n-1, 1, -1
            u(j) = u(j) - gam(j+1) * u(j+1)
      end do

      end subroutine tridia
!c=======================================================================
!c
!!!       SUBROUTINE cntrol(icho2,ichmet,iopt)
        subroutine cntrol(icho2, iopt)

!c=======================================================================
!c
!c      initialisation and switches.
!c
!c-----------------------------------------------------------------------

   use initialv
!   use steps, only : dtt, itrunk, ilw

   implicit none

   ! Dummy arguments
   integer, intent(in) :: icho2, iopt

    ! Essas variáveis já foram inicializadas*
!!!	dtt    = 3600.
!!!	itrunk = 20
!!!	ilw    = 5 !tatsch
!!!
!!!	tc_ini = 298.0
!!!	tg_ini = 298.0
!!!	td_ini = 297.0
!!!
!!!	www_ini(1) = 0.70
!!!	www_ini(2) = 0.85
!!!	www_ini(3) = 0.98
        print*,icho2
      IF (iopt /= 1) THEN

!      WRITE(icho2,800)
!  800 FORMAT(10x,32('*')/10x,'*Sib2 off-line simulation in MP*'/10x, &
!         32('*')/5x)
!      IF(itrunk == 1) WRITE(icho2,801)
!      IF(itrunk >= 2) WRITE(icho2,802) itrunk
!  801 FORMAT(5x,'resistances calculated from initial fluxes')
!  802 FORMAT(5x,'resistances calculated by iteration, itrunk=',i4)
!      IF(ilw == 1) WRITE(icho2,816)
!      IF(ilw == 2) WRITE(icho2,817)
!      IF(ilw == 3) WRITE(icho2,818)
!      IF(ilw == 4) WRITE(icho2,819)
!  816 FORMAT(5x,'downward longwave radiation read in as data'// )
!  817 FORMAT(5x,'downward longwave radiation computed from brunts', &
!         ' equation', // )
!  818 FORMAT(5x,'downward longwave radiation computed as residual in ene &
!         rgy balance', /, 5x,'net radiation read in as data ',// )
!  819 FORMAT(5x,'downward longwave radiation computed from Idso&Jackson', &
!         ' equation', // )

      END IF
end subroutine cntrol



!!!c=======================================================================
!!!c
!!!      SUBROUTINE radc2_1
!!!c
!!!c=======================================================================
!!!c
!!!c     solar zenith angle computation; downcoming radiation at bottom.
!!!c
!!!c-----------------------------------------------------------------------
!!!c
!!!        use govern, only : time, year, day
!!!        use steps, only : dtt, iter
!!!
!!!        IMPLICIT NONE
!!!
!!!        REAL dayspy
!!!
!!!CCC      INCLUDE 'COMSIBC.H'
!!!c
!!!      IF ( amod( year, 4. ) == 0. ) THEN
!!!          dayspy = 366.
!!!      ELSE
!!!          dayspy = 365.
!!!      END IF
!!!c
!!!c-----------------------------------------------------------------------
!!!c    julian day and time update; skip on 1st time step (initialized)
!!!c-----------------------------------------------------------------------
!!!      IF(iter > 1)THEN
!!!         time = time + dtt / 3600.
!!!         IF ( time >= 23.99 ) time = 0.0
!!!         day = day +  dtt / 86400.
!!!      END IF
!!!c
!!!      IF ( day > dayspy ) THEN
!!!           year = year + 1.
!!!		 day = day - dayspy
!!!      END IF
!!!
!!!	END



!!!      SUBROUTINE radc2_2
subroutine radc2_2(zlat, sindec, cosdec, coshr, sunang)

!c=======================================================================
!c
!c     solar zenith angle computation; downcoming radiation at bottom.
!c
!c-----------------------------------------------------------------------

!!!        use constants, only : decmax
!!!        use atmos, only : sunang
!!!        use const, only : pie
   use constants, only : s2r
!!!        use govern, only : time, year, day
!!!        use govern, only : sols, time, realday
!!!        use site, only : zlat
!!!        use steps, only : dtt

   implicit none

!!!        REAL decmax, sols, season, dec, rfd, sind, cosd, hac
!!!        REAL season, dec, rfd, sind, cosd, hac
   ! Dummy arguments
   real, intent(in) :: zlat, sindec, cosdec, coshr
   real :: sunang
!!!        REAL dawn, dusk, sr, ss, coshr, h
   ! Local variables
   real :: lat

   lat = s2r * zlat

   sunang = sin(lat) * sindec + cos(lat) * cosdec * coshr
   sunang = max(0.01, sunang)
   
!c-----------------------------------------------------------------------
!c    solar declination calculation
!c-----------------------------------------------------------------------

!!!      decmax = pie * ( 23.5 / 180.)
!!!      sols   = ( 4141./24. ) + amod( year+3., 4. ) * 0.25

!!!      season = ( day - sols ) / 365.2
!!!      season = (realday - sols) / 365.2
!!!      dec    = decmax * cos ( 2. * pie * season )

!!!      rfd  = pie / 180.
!!!      sind = sin( dec )
!!!      cosd = cos( dec )
!!!      hac  = -tan( zlat * rfd )*tan( dec )
!!!      hac  = -tan(lat) * tandec
!!!      hac  = min(hac,1.0)
!!!      hac  = max(hac,-1.0)

!c-----------------------------------------------------------------------
!c     h is the half-day length (in radians)
!c-----------------------------------------------------------------------

!!!      h   = acos(hac)
!!!      dawn= -h
!!!      dusk= +h
!!!      sr  = 12.-(h/(15.*rfd))
!!!      ss  = 12.+(h/(15.*rfd))
!!!      FAC = H / (s2r * 15.0)
!!!      SR = 12.0 - FAC
!!!      SS = 12.0 + FAC
!!!      coshr = cos( - pie + (time + 0.5*dtt/3600.) / 24. * 2. * pie )
!!!      coshr = cos(-pi + time / 24.0 * twopi)
!!!      sunang = sin( zlat*rfd ) * sind + cos ( zlat*rfd ) * cosd * coshr
end subroutine radc2_2



!c======================================================================

subroutine rasite

!c======================================================================
!c
!c     calculation of ustar, u2, ra and drag using Paulson's method.
!c
!c----------------------------------------------------------------------
!c
!c     (1) site parameters derived from momopt program suite.
!c
!c     (2) routine is not suitable for gcm applications; designed for
!c         use with forcing variables measured at a field site.
!c
!c     (3) paulson psi-coefficients are constrained under unsatble
!c         conditions to prevent unrealistically low ra values.
!c
!c     (4) wind speed (um) must be greater than or equal to 0.1 m/s
!c
!c----------------------------------------------------------------------
!c
!c     variables that must enter through comsibc
!c
!c      tm     : air temperature at zmet
!c      um     : wind speed at zwind, um >= 0.1
!c      ht     : sensible heat flux from surface
!c
!c     parameters that must enter through comsibc
!c
!c      z2     : height of canopy top
!c      z0     : roughness length
!c      d      : zero plane displacement
!c      vkc    : von karmans constant = 0.41
!c      rhoair : air density
!c      cpair  : air specific heat
!c
!c     other parameters
!c     ----------------
!c
!c      g1, g2, g3, ztz0, corb1, corb2, ha, zwind, zmet
!c
!c      g1     : ratio of km(actual) to km(log-linear) at z = z2
!c      g2     : ratio of ra(actual) to ra(log-linear) for momentum
!c               between: z = z2 and z = zx, where zx = min(zl,zwind)
!c      g3     : ratio of ra(actual) to ra(log-linear) for heat
!c               between: z = z2 and z = zx, where zx = min(zl,zmet)
!c      ztz0   : parameter to determine depth of transition layer above
!c               canopy, zl. zl = z2 + ztz0 * z0
!c      corb1  : non-neutral correction for calculation of aerodynamic
!c               resistance between ha and z2. when multiplied by
!c               h*rbb/(tm*u2*u2*(z2-ha)) gives bulk estimate of local
!c               richardson number.
!c               rbb = ra for heat between ha and z2.
!c               Rib = corb1*h*rbb/(tm*u2*u2*(z2-ha))
!c               corb1 = 9*g/( rhoair*cpair* (du/dz)**2 ) *u2*u2
!c      corb2  : neutral value of rbb*u2 ( squared ), equivalent to
!c               rdc**2 for upper canopy
!c      ha     : canopy source height for heat
!c      zwind  : reference height for wind measurement
!c      zmet   : reference height for temperature, humidity measurement
!c
!c        the above are generated from sibx + momopt output
!c
!c-----------------------------------------------------------------------
!c
!c     variables returned from this routine via comsibc
!c
!c      ustar  : friction velocity
!c      u2     : wind speed at canopy top
!c      ra     : aerodynamic resistance for heat flux between ha and zmet
!c      drag   : shear stress at canopy top
!c
!c-----------------------------------------------------------------------
!c
!c     references
!c     ----------
!c
!c         Paulson C.A. (1970) ' Mathematical representation of wind
!c         and temperature profiles in the unstable atmospheric surface
!c         layer', J. Appl. Met., 9, 129-861.
!c
!c         SE-89
!c-----------------------------------------------------------------------

   use aerorx, only : ra, rbbest
   use atmos, only : tm, um
   use caerod, only : g2, ztz0, corb1, corb2, ha, zwind, ht
   use const, only : vkc, rhoair, cpair, g
   use donor, only : drag
   use grads, only : ustar, u2
   use rause, only : z0, d
   use steps, only : itrunk
   use vstate, only : z2

   implicit none

   ! Local variables
   integer :: i, iwalk, lx, nonpos, nox
   real :: arg1, bot, coef3, finc, gfac, hm1, hm2, hrb, hress, hss
   real :: ps1, ps2, raf, raf1, rafmax, ram, top, uest, us1, us2, uss
   real :: y, zl, zx1, zx2

!!!     NOVAS VARIAVEIS PARA O NEWTON - TESTE
   integer :: iter_aux
   real :: zinc_aux, a2_aux, y1_aux

   us1 = 0.0

   hress = ht
   zl = z2 + ztz0 * z0
   uest = vkc * um / log((zwind - d) / z0)

!c-----------------------------------------------------------------------
!c
!c     calculation of u2 assuming neutral conditions
!c
!c-----------------------------------------------------------------------

   if (zwind <= zl) then
      top = 0.0
      zx1 = zwind - d
      zx2 = z2 - d
   else
      zx1 = zwind - d
      zx2 = zl - d
      top = log(zx1 / zx2)
      zx1 = zl - d
      zx2 = z2 - d
   end if
   bot = log(zx1 / zx2)
   ram = 1.0 / (vkc * uest) * (top + g2 * bot)
   u2 = um - ram * uest ** 2

!c-----------------------------------------------------------------------
!c
!c     calculation of ra for heat follows : non-neutrality assumed
!c
!c-----------------------------------------------------------------------

   zx1 = zwind - d
   zx2 = 0.0
   arg1 = log(zx1 / z0)

!c-----------------------------------------------------------------------
!c         initialize newton-raphson iterative routine
!c-----------------------------------------------------------------------

   nox = 0
   nonpos = 1
   iwalk = 0
   lx = 1
   finc = 0.2

   iter_aux = 0
   zinc_aux = 0.0
   a2_aux = 0.0
   y1_aux = 0.0

   if (ht > 0.0) then

!c-----------------------------------------------------------------------
!c
!c     unstable case : calculation of ustar followed by ra
!c
!c-----------------------------------------------------------------------

      i = 0
      do while (nox == 0)
         call unstab(uest, zx1, zx2, arg1, ht, ps1, ps2)

         y = um - uest / vkc * (arg1 - ps1)

         if (abs(y - uest) < 0.0000001) then
            y = uest - 0.0000001
         end if
         i = i + 1
         if (i > itrunk) goto 44441
!!!          CALL newton ( uest, y, finc, nox, nonpos, iwalk, lx )
!!!         call newton_new( uest, y, finc, nox, nonpos, iwalk, lx, &
         call newton_new( uest, y, finc, nox, nonpos, iwalk, &
            zinc_aux, a2_aux, y1_aux, iter_aux)
      end do
44441		continue

      if (nox == 2) write(6,900)
900     FORMAT( /,' convergence failure in rasite - unstable case' )

      call rafcal ( zl, uest, ht, raf )

   else
!c-----------------------------------------------------------------------
!c
!c      stable case : calculation of ustar
!c
!c-----------------------------------------------------------------------
!c
!c
!c----------------------------------------------------------------------
!c
!c      interpolation zone is defined: this spans negative sensible
!c      heat fluxes from positive side ( hence factor 0.95 ) of the
!c      following conditions:
!c              y = 0,    dy/du* = 0.
!c      see notes for details
!c
!c----------------------------------------------------------------------

      gfac = log((zwind - d)/z0)
      hm1 = -0.95 * tm * rhoair * cpair / (2.0 * 4.7 * g * (zwind - d)) &
         * (2.0 * um / 3.0) ** 3 * (vkc / gfac) ** 2
      hm2 = 5.0 * hm1
      us2 = vkc * um / (gfac + 4.7)
      if (ht >= hm2) then
         ht = max(hm1, ht)

!c----------------------------------------------------------------------
!c
!c      ustar calculated for slightly stable conditions : ht >= hm1
!c
!c----------------------------------------------------------------------

         i = 0
         do while (nox == 0)
!!!             CALL stab ( uest, zx1, zx2, arg1, ht, ps1, ps2)
            call stab(uest, zx1, zx2, ht, ps1, ps2)

            y = um - uest / vkc * (arg1 - ps1)

            i = i + 1
            if (i > itrunk) goto 44442
!!!             CALL newton ( uest, y, finc, nox, nonpos, iwalk, lx )
!!!            call newton_new(uest, y, finc, nox, nonpos, iwalk, lx, &
            call newton_new(uest, y, finc, nox, nonpos, iwalk, &
               zinc_aux, a2_aux, y1_aux, iter_aux)
         end do
44442			continue

         if (nox == 2) write(6,910)
910        FORMAT( /,' convergence failure in rasite - stable case' )

         ht = hress

!c-----------------------------------------------------------------------
!c      ustar calculation in interpolation zone
!c-----------------------------------------------------------------------

         if (ht <= hm1) then
            us1 = uest
            uest = (ht - hm2) / (hm1 - hm2) * (us1 - us2) + us2
         end if

!c-----------------------------------------------------------------------
!c      ustar calculation for collapsed profiles
!c-----------------------------------------------------------------------

      else
         uest = us2
      end if   ! (ht>=hm2)

!c-----------------------------------------------------------------------
!c
!c      calculation of ra for heat transfer between z2 and zmet
!c
!c-----------------------------------------------------------------------

      raf = 1.e5

      call rafcal(zl, us2, hm2, rafmax)

      if (ht >= hm2) then
         hss = max(hm1, ht)
         uss = max(us1, uest)

         call rafcal(zl, uss, hss, raf)

         if (ht <= hm1) then
            raf1 = raf
            raf = (ht - hm2) / (hm1 - hm2) * (raf1 - rafmax) + rafmax
         end if
      end if

      raf = min(raf, rafmax)

!c-----------------------------------------------------------------------
!c     above canopy variables calculated.
!c-----------------------------------------------------------------------

   end if   ! (ht>0)

!      hrb = (ht + sqrt(ht ** 2)) / 2.0 + 0.1
   hrb = (ht + abs(ht)) / 2.0 + 0.1

!c-----------------------------------------------------------------------
!c     corb1 and corb2 are calculated for between ha and z2 only.
!c-----------------------------------------------------------------------

   rbbest = sqrt(corb2) / u2

!c-----------------------------------------------------------------------
!c           initialize newton-raphson iterative routine
!c-----------------------------------------------------------------------

   nox = 0
   nonpos = 1
   iwalk = 0
   lx = 1
   finc = 0.2

   iter_aux = 0
   zinc_aux = 0.0
   a2_aux = 0.0
   y1_aux = 0.0

   i = 0
   do while (nox /= 1)
      coef3 = corb1 * hrb / tm / ( z2-ha )

      y = coef3 * rbbest**3 + ( u2*rbbest )**2 - corb2

      i = i + 1
      if (i > itrunk) goto 44443
!!!        CALL newton( rbbest , y, finc , nox, nonpos, iwalk, lx)
!!!      call newton_new(rbbest , y, finc , nox, nonpos, iwalk, lx, &
      call newton_new(rbbest , y, finc , nox, nonpos, iwalk, &
         zinc_aux, a2_aux, y1_aux, iter_aux)
   end do
44443		continue

   ra = raf + rbbest
   ustar = uest
   drag = rhoair * uest * uest
end subroutine rasite



!c=======================================================================

subroutine unstab(uest, a, b, argz, heat, psione, psitwo)

!c=======================================================================
!c
!c      calculation of Paulson psi-function for unstable condition
!c
!c-----------------------------------------------------------------------

   use atmos, only : tm
   use const, only : rhoair, cpair, g, vkc

   implicit none

   ! Dummy arguments
   real, intent(in) :: uest, a, b, argz, heat
   real :: psione, psitwo

   ! Local variables
   integer :: i
   real :: zin, zml, fac
   real :: x(2)

   zin = a

   do i = 1, 2
      zml = -uest ** 3 * rhoair * cpair * tm
      zml = zml / (vkc * g * heat)
      fac = 16.0 * zin / zml
      x(i) = (1.0 - fac) ** 0.25
      zin = b
   end do

   psione = 2.0 * log((1.0 + x(1)) / (1.0 + x(2))) + log((1.0 + x(1) ** 2) &
      / (1.0 + x(2) ** 2)) - 2.0 * atan(x(1)) + 2.0 * atan(x(2))
   psione = min(argz * 0.75, psione)

   psitwo = 2.0 * log((1.0 + x(1) ** 2) / (1.0 + x(2) ** 2))
   psitwo = min(argz * 0.75, psitwo)
end subroutine unstab



!c=======================================================================

!!!      SUBROUTINE stab ( uest, a, b, argz, heat, psione, psitwo )
subroutine stab(uest, a, b, heat, psione, psitwo)

!c=======================================================================
!c
!c      calculation of Paulson psi-function for stable condition
!c
!c-----------------------------------------------------------------------

   use atmos, only : tm
   use const, only : rhoair, cpair, g, vkc

   IMPLICIT NONE

   ! Dummy arguments
   real, intent(in) :: uest, a, b, heat
   real :: psione, psitwo

   ! Local variables
   real zml

   psione = 0.0
   psitwo = 0.0
   if (abs(heat) > 1.e-4) then
      zml = -uest ** 3 * rhoair * cpair * tm
      zml = zml / (vkc * g * heat)

      psione = -4.7 * (a - b) / zml
      psione = max(-4.7, psione)

      psitwo = psione
   end if
end subroutine stab




!c=======================================================================

subroutine rafcal(zl, uest, heat, raf)

!c=======================================================================
!c
!c      calculation of ra for heat between z2 and zmet
!c
!c-----------------------------------------------------------------------

   use caerod, only : g3, zmet
   use const, only : vkc
   use rause, only : d
   use vstate, only : z2

   implicit none

   ! Dummy arguments
   real, intent(in) :: zl, uest, heat
   real :: raf

   ! Local variables
   real arg, ps1, ps2, zx1, zx2, bot, top

   if (zmet <= zl) then

      top = 0.0
      zx1 = zmet - d
      zx2 = z2 - d

   else

      zx1 = zmet - d
      zx2 = zl - d
      arg = log(zx1 / zx2)
      if (heat > 0.0) then
         call unstab(uest, zx1, zx2, arg, heat, ps1, ps2)
      else
         call stab(uest, zx1, zx2, heat, ps1, ps2)
      end if
      top = arg - ps2

      zx1 = zl - d
      zx2 = z2 - d

   end if

   arg = log (zx1 / zx2)
   if (heat > 0.0) then
      call unstab(uest, zx1, zx2, arg, heat, ps1, ps2)
   else
      call stab(uest, zx1, zx2, heat, ps1, ps2)
   end if
   bot = arg - ps2

   raf = 1.0 / (vkc * uest) * (top + g3 * bot)
end subroutine rafcal



!c=======================================================================

       subroutine newton(a1, y, finc, nox, nonpos, iwolk, l)

!c=======================================================================
!c
!c      the newton raphson iterative routine will be used to generate new
!c      values of a1 if dabsolute value of y is greater than ertol;
!c      a1 is estimate, y is resultant error
!c      nex is exit condition  (0=no exit) or (1 when dabs(y) lt ertol)
!c      ertol is the dabsolute value of y necessary to obtain an exit
!c      finc is initial increment size for second estimate of a1
!c      nonpos=0 if quantity to be minimized can be less than zero;
!c      nonpos=1 if quantity can only be positive
!c      l identifies which quantity is being calculated.
!c
!c      control values: finc,ertol,nox,nonpos,l:must be set by user
!c-----------------------------------------------------------------------

        IMPLICIT NONE

        INTEGER nox, nonpos, iwolk, l
        REAL a1, y, finc

       integer iter(3), iwalk(3), nex(3)
       real zinc(3), a2(3), y1(3)

        REAL cons
       data cons/1.0/

        REAL a, ertol, step

        NAMELIST /VARS/ A1, Y, Y1


!!!        TESTE
!!!        ITER = 0
!!!        Y1 = 0

!!!        cons = 1.0

        WRITE(*,NML=VARS)

       ertol = 0.05 * finc
       iwalk(l) = iwolk
       nex(l)=nox

       if ( iter(l) >= 490 ) go to 160
       if (ertol < 0.000001) ertol=0.000001
       if (abs(y) <= ertol) go to 150
       if((abs(y-y1(l)))<=0.01*ertol .and. iwalk(l)==0 ) go to 8

       if(abs(y1(l))>ertol) go to 1
       a2(l)=a1
!c**    a1=a1-y
       step = min( abs(y), abs(10.*finc) ) * sign(cons,y)
       a1=a1-step
       nex(l)=0
       y1(l)=y
       iter(l)=1
       if (iwalk(l) == 3) go to 101
       iwalk(l)=0
       go to 101
   1   iter(l)=iter(l)+1
       if(iter(l) == 20) iwalk(l)=1
       if(iwalk(l) /= 0) go to 2
       if(abs(y) > ertol) go to 3
       nex(l)=1
       go to 150
   3   a=a1-y*(a1-a2(l))/(y-y1(l))
       if(abs(a-a1)>(10.0*finc)) &
                 a=a1+10.0*finc*sign(cons,(a-a1))
       a2(l)=a1
       a1=a
       y1(l)=y
       go to 101
   2   if(iwalk(l)==2)go to 4
       if(iwalk(l)==3) go to 6
       if(sign(cons,y)==sign(cons,y1(l))) go to  3
       zinc(l)=(a1-a2(l))/4.0
       a1=a2(l)+zinc(l)
       iwalk(l)=2
       nex(l)=0
       go to 101
   4   if(sign(cons,y) ==sign(cons,y1(l))) go to 5
       zinc(l)=-zinc(l)/4.0
       a2(l)=a1
       a1=a1+zinc(l)
       nex(l)=0
       y1(l)=y
       go to 101
   5   a2(l)=a1
       a1=a1+zinc(l)
       y1(l)=y
       nex(l)=0
       go to 101
   6   if(sign(cons,y)==sign(cons,y1(l))) go to 7
       iwalk(l)=1
       go to 2
   7   a2(l) = a1
       a1 = a1+finc
       y1(l)=y
       nex(l) = 0
       go to 101
   8   a1 = a1 + finc*2.0
       nex(l)=0
       go to 101
160    continue
       write(6,900) y, l
 900   format ( 3x,' failure to converge after 490 iterations', &
         /, 3x,' y = ',g12.5, ' lx =',i2 )

 150   nex(l) = 1
       if( iter(l) >= 490 ) nex(l) = 2
       zinc(l)=0.0
       iter(l) = 0
       iwalk(l)=0
       y1(l)=0.0
       y=0.0
       a2(l)=0.0
 101   continue
       if(nonpos==1.and.a1<0.0) a1=a2(l)/2.0
       nox = nex(l)
       iwolk = iwalk(l)

        WRITE(*,NML=VARS)
        WRITE(6,*) 'SAINDO. L', L, 'ITER', ITER(L)
      RETURN
      end subroutine newton




!c=======================================================================

!!!       subroutine newton_new(a1, y, finc, nox, nonpos, iwolk, l, &
       subroutine newton_new(a1, y, finc, nox, nonpos, iwolk, &
           ZINC, A2, Y1, ITER)

!c=======================================================================
!c
!!!     THIS ROUTINE WAS MODIFIED BY NELSONVN ON 20150618
!c
!c      the newton raphson iterative routine will be used to generate new
!c      values of a1 if dabsolute value of y is greater than ertol;
!c      a1 is estimate, y is resultant error
!c      nex is exit condition  (0=no exit) or (1 when dabs(y) lt ertol)
!c      ertol is the dabsolute value of y necessary to obtain an exit
!c      finc is initial increment size for second estimate of a1
!c      nonpos=0 if quantity to be minimized can be less than zero;
!c      nonpos=1 if quantity can only be positive
!c      l identifies which quantity is being calculated.
!c
!c      control values: finc,ertol,nox,nonpos,l:must be set by user
!c-----------------------------------------------------------------------

        IMPLICIT NONE

!!!     Dummy arguments
!!!        real, parameter :: cons = 1.0
!!!        integer, intent(in) :: nonpos, l
!!!        integer nonpos, l
        integer nonpos
        integer  nox, iwolk
!!!        real, intent(in) :: finc
        real finc
        real a1, y

!!!     Local variables
        integer iter, iwalk, nex
        real zinc, a2, y1
        real a, ertol, step

        REAL cons
       data cons/1.0/

!!!        NAMELIST /VARS/ A1, Y, Y1, ITER

!!!        TESTE
!!!        ITER = 0
!!!        Y1 = A1 / 2.0
!!!        Y1 = Y + 10.0*ertol

!!!        cons = 1.0

!!!        WRITE(*,NML=VARS)

       ertol = 0.05 * finc
       iwalk = iwolk
       nex = nox

       if (iter >= 490) go to 160
       if (ertol < 0.000001) ertol = 0.000001
       if (abs(y) <= ertol) go to 150
       if (abs(y - y1) <= 0.01 * ertol .and. iwalk == 0) go to 8

       if (abs(y1) > ertol) go to 1
       a2 = a1
!c!!        a1=a1-y   !CCC linha descomentada no WRF
       step = min(abs(y), abs(10.0 * finc)) * sign(cons,y)   !CCC linha inexistente no WRF
       a1 = a1 - step   !CCC linha inexistente no WRF
       nex = 0
       y1 = y
       iter = 1
       if (iwalk == 3) go to 101
       iwalk = 0
       go to 101

   1   iter = iter + 1
       if (iter == 20) iwalk = 1
       if (iwalk /= 0) go to 2
       if (abs(y) > ertol) go to 3
       nex = 1
       go to 150

   3   a = a1 - y * (a1 - a2) / (y - y1)
       if (abs(a - a1) > 10.0 * finc) &
                 a = a1 + 10.0 * finc * sign(cons, a - a1)
       a2 = a1
       a1 = a
       y1 = y
       go to 101

   2   if (iwalk == 2) go to 4
       if (iwalk == 3) go to 6
       if (sign(cons, y) == sign(cons, y1)) go to  3
       zinc = (a1 - a2) / 4.0
       a1 = a2 + zinc
       iwalk = 2
       nex = 0
       go to 101

   4   if (sign(cons, y) == sign(cons, y1)) go to 5
       zinc = -zinc / 4.0
       a2 = a1
       a1 = a1 + zinc
       nex = 0
       y1 = y
       go to 101

   5   a2 = a1
       a1 = a1 + zinc
       y1 = y
       nex = 0
       go to 101

   6   if (sign(cons, y) == sign(cons, y1)) go to 7
       iwalk = 1
       go to 2

   7   a2 = a1
       a1 = a1 + finc
       y1 = y
       nex = 0
       go to 101

   8   a1 = a1 + finc * 2.0
       nex = 0
       go to 101

160    continue
!!!       write(6,900) y, l   !CCC linha inexistente no WRF
!!! 900   format ( 3x,' failure to converge after 490 iterations',
!!!     & /, 3x,' y = ',g12.5, ' lx =',i2 )
        WRITE(6,'(A30,I8,A12)') 'FAILURE TO CONVERGE AFTER', ITER, &
       ' ITERATIONS'

 150   nex = 1
       if (iter >= 490) nex = 2   !CCC linha inexistente no WRF
       zinc = 0.0
       iter = 0
       iwalk = 0
       y1 = 0.0
       y = 0.0
       a2 = 0.0

 101   continue
       if (nonpos == 1 .and. a1 < 0.0) a1 = a2 / 2.0
       nox = nex
       iwolk = iwalk

!!!        WRITE(*,NML=VARS)
!!!        WRITE(6,*) 'SAINDO. ITER', ITER
      RETURN
      END subroutine newton_new



!c=======================================================================

subroutine gauss(a, n, np1, x, work)

!c=======================================================================
!c
!c     solve a linear system by gaussian elimination.  developed by
!c     dr. chin-hoh moeng.  a is the matrix of coefficients, with the
!c     vector of constants appended as an extra column.  x is the vector
!c     containing the results.  the input matrix is not destroyed.
!c
!c-----------------------------------------------------------------------

   implicit none

   ! Dummy arguments
   integer, intent(in) :: n, np1
   real, intent(in) :: a(4,5)
   real :: work(4,5), x(4)

   ! Local variables
   integer :: i, j, k, l
   real :: r

   do i = 1, n
      do j = 1, np1
         work(i,j) = a(i,j)
      end do
   end do

   do i = 2, n
      do j = i, n
         r = work(j,i-1) / work(i-1,i-1)
         do k = 1, np1
            work(j,k) = work(j,k) - r * work(i-1,k)
         end do
      end do
   end do

   do i = 2, n
      k = n - i + 2
      r = work(k,np1) / work(k,k)
      do j = i, n
         l = n - j + 1
         work(l,np1) = work(l,np1) - r * work(l,k)
      end do
   end do

   do i = 1, n
      x(i) = work(i,np1) / work(i,i)
   end do
end subroutine gauss

!cRHV===================================================================
!
!         subroutine runTn 
! 
! 
! !!!c----------------------------------------------------------------------
! !!!c
! !!!c    calculation of inter-layer exchanges of water due to gravitation
! !!!c    and hydraulic gradient. the values of w(x) + dw(x) are used to
! !!!c    calculate the potential gradients between layers.
! !!!c    modified calculation of mean conductivities follows ME-82 ),
! !!!c    reduces recharge flux to top layer.
! !!!c
! !!!c      dpdw           : estimated derivative of soil moisture potential
! !!!c                       with respect to soil wetness. assumption of
! !!!c                       gravitational drainage used to estimate likely
! !!!c                       minimum wetness over the time step.
! !!!c
! !!!c      qqq  (q     )  : equation (61) , SE-86
! !!!c             i,i+1
! !!!c            -
! !!!c      avk  (k     )  : equation (4.14) , ME-82
! !!!c             i,i+1
! !!!c
! !!!c----------------------------------------------------------------------
! 
! !    wmax = max(www(1), www(2), www(3), 0.05)
!     wmax = max(www(1), dzwww, 0.05)
!     wmax = min(wmax, 1.0)
!     pmax = wmax ** (-bee)
!     wmin = (pmax - (zdepth(1) + zdepth(2) + zdepth(3)) / phsat) **(-1.0 / bee)
! 
! !    wmin = min(www(1), www(2), www(3), wmin)
!     wmin = min(www(1),dzwww, wmin)
!     wmin = max(wmin, 0.02)
!     pmin = wmin ** (-bee)
!     dpdw = phsat * (pmax - pmin) / (wmax - wmin)
! 
! !!!c-----------------------------------------------------------------------
! !!!c     conductivities and base flow reduced when temperature drops below
! !!!c     freezing.
! !!!c-----------------------------------------------------------------------
!       tsnow  = min(tf - 0.01, tg)
!       tgs    = tsnow * areas + tg * (1.0 - areas)
!       ts     = tgs * (2 - i) + td * (i - 1)
!       props  = (ts - (tf - 10.0)) / 10.0
!       props  = max(0.05, min(1.0, props))
! !!!c-----------------------------------------------------------------------   
! 	      
! !    do i = 1, 2
!      DO i=1, nlayers-1
!       rsame  = 0.0
!       
! 	  temw_d(i)   = max(0.03, dzwww(i))               ! teta/tetas
! 	  temwp_d(i)  = temw(i) ** (-bee)                     ! psi
! 	  temwpp_d(i) = min(1.0, temw(i)) ** (2.0 * bee + 3.0) 
! 	  
!       if(i.eq.1)then
! 	
! 	  temw(i)   = max(0.03, www(i))               ! teta/tetas
! 	  temwp(i)  = temw(i) ** (-bee)                     ! psi
! 	  temwpp(i) = min(1.0, temw(i)) ** (2.0 * bee + 3.0) 
! 	  
! 	avk = temwp(i) * temwpp(i) - temwp_d(1) * temwpp_d(i)
! 	div = temwp(i+1) - temwp(i)
! 	if (abs(div) < 1.e-6) rsame = 1.0
! 	avk    = satcoz(i) * avk / ((1.0 + 3.0 / bee) * div + rsame)
! 	avkmin = satcoz(i) * min(temwpp(i), temwpp_d(i))
! 	avkmax = satcoz(i) * max(temwpp(i), temwpp_d(i)) * 1.01
! 	avk    = max(avk, avkmin)
!         avk    = min(avk, avkmax)
!         if (www(i) < dzwww(i)) avk = 0.1 * avk
!         avk    = avk * props
!         
! 	  dpdwdz = dpdw * 2.0 / (zdepth(i) + dzmultilayer)
! 	  aaa(i) = 1.0 + avk * dpdwdz * (1.0 / zdepth(i) +1.0 / dzmultilayer) &
! 	    * dtt / poros
! 	  bbb(i) = -avk * dpdwdz * 1.0 / dzmultilayer * dtt / poros
! 	  ddd(i) = avk * (dpdwdz * (www(i) - dzwww(i)) + 1.0)
!       avkii(i) = avk
!       endif ! i==1
!       
!       
! 	! estoy por aqui
!       avk    = temwp(i) * temwpp(i) - temwp(i+1) * temwpp(i+1)
!       div    = temwp(i+1) - temwp(i)
!       if (abs(div) < 1.e-6) rsame = 1.0
! !!!c!        avk = satco*avk / ( ( 1. + 3./bee ) * div + rsame )
!       avk    = satcoz(i) * avk / ((1.0 + 3.0 / bee) * div + rsame)
! !!!c!        avkmin = satco * amin1( temwpp(i), temwpp(i+1) )
!       avkmin = satcoz(i) * min(temwpp(i), temwpp(i+1))
! !!!c!        avkmax = satco * amax1( temwpp(i), temwpp(i+1) )*1.01
!       avkmax = satcoz(i) * max(temwpp(i), temwpp(i+1)) * 1.01
!       avk    = max(avk, avkmin)
!       avk    = min(avk, avkmax)
! !!!c 3/96 changes
!       if (www(i) < www(i+1)) avk = 0.1 * avk
! 
!       avk    = avk * props
! 
! !!!c-----------------------------------------------------------------------
! !!!c     backward implicit calculation of flows between soil layers.
! !!!c-----------------------------------------------------------------------
!       dpdwdz = dpdw * 2.0 / (zdepth(i) + zdepth(i+1))
!       aaa(i) = 1.0 + avk * dpdwdz * (1.0 / zdepth(i) +1.0 / zdepth(i+1)) &
!          * dtt / poros
!       bbb(i) = -avk * dpdwdz * 1.0 / zdepth(2) * dtt / poros
!       ddd(i) = avk * (dpdwdz * (www(i) - www(i+1)) + 1.0)
!       if (i == 2) then
!          ccc(i) = -avk * dpdwdz * 1.0 / zdepth(3) * dtt / poros
!       endif
!       avkii(i) = avk
!    end do   ! (i=1,2)
! 
!    ccc(1) = 0.0
!    ccc(3) = 0.0
! 
!    d3g = GWdep - (zdepth(1) + zdepth(2) + zdepth(3) / 2.0)
! 	d3g = max(d3g, zdepth(3) / 2.0, 0.5)
! 
! !!!c!	avk = satco * temw(3)**( 2.*bee+3.)
! 	avk = satcoz(3) * temw(3) ** (2.0 * bee + 3.0)                  !tatsch (i)
! !!!c!	avk = (satco* exp(-decay*d3g)) * temw(3)**( 2.*bee+3.)   !tatsch (ii)
! !!!c
! !!!c!	avk = (satco* exp(-decay*d3g)) * temw(3)**( 2.*bee+3.)   !tatsch (ii)
! 
! !!!c        print*, '-----dentro da run2-------'
! !!!c        print *, item01, item02, ivtype,'       d3g:',d3g
! !!!c        print *,'z1-3:',zdepth(1)+zdepth(2)+zdepth(3)/2.,'GWdep',GWdep
! !!!c        print*, 'satco:',satco*1000*3600,'satcoz3:',satcoz(3)*1000*3600,
! !!!c     &  'decay:', decay
! !!!c        print*,'avk          ', 1000*3600*avk
! !!!c        print*,'avk decay d3g',
! !!!c     &   1000*3600*satco*exp(-decay*d3g)*temw(3)**(2.*bee+3.)
! !!!c        pause
! 
!    avkii(3) = avk
! 
!    wmax     = 1.0
!    pmax     = 1.0
!    wmin     = (pmax- d3g / phsat) ** (-1.0 / bee)
!    wmin     = min(www(3), wmin)
!    wmin     = max(wmin, 0.02)
!    pmin     = wmin ** (-bee)
!    dpdw     = phsat * (pmax - pmin) / (wmax - wmin)
! 
!    dpdwdz   = dpdw / d3g
!    aaa(3)   = 1.0 + avk * dpdwdz * (1.0 / zdepth(3)) * dtt / poros
!    bbb(3)   = -avk * dpdwdz * 1.0 / zdepth(3) * dtt / poros
!    ddd(3)   = avk * (dpdwdz * www(3) + 1.0 - dpdwdz)
! 
!    denom    = (aaa(1) * aaa(2) * aaa(3) - aaa(3) * bbb(1) * bbb(2) - aaa(1) &
!       * bbb(3) * ccc(2))
!    rdenom   = 0.0
!    if ( abs(denom) < 1.e-6 ) rdenom = 1.0
!    rdenom   = (1.0 - rdenom) / (denom + rdenom)
!    qqq(1)   = (aaa(2) * ddd(1) * aaa(3) - bbb(3) * ccc(2) * ddd(1) - aaa(3) &
!       * ddd(2) * bbb(1) + ccc(2) * ddd(3) * bbb(1) ) * rdenom
!    qqq(2)   = (aaa(1) * aaa(3) * ddd(2) - aaa(3) * ddd(1) * bbb(2) - aaa(1) &
!       * ccc(2) * ddd(3)) * rdenom
!    q3g      = (aaa(1) * aaa(2) * ddd(3) - bbb(1) * bbb(2) * ddd(3) - aaa(1) &
!       * bbb(3) * ddd(2) + bbb(3) * ddd(1) * bbb(2)) * rdenom
! 
!    q3gmax   =  www(3) * (poros * zdepth(3) / dtt)
!    q3gmin   = (www(3) - 1.0) * (poros * zdepth(3) / dtt)
!    q3gmin   = 0.0
!    q3g      = min(q3g , q3gmax)  !, avkii(3)
!    q3g      = max(q3g , q3gmin)	!, -avkii(3)
! 
! 
! 	end subroutine runTn



!  subroutine soilwater(nl_soil,deltim,wimp,smpmin,&
!                       qinfl,etr,z_soisno,dz_soisno,zi_soisno,&
!                       t_soisno,vol_liq,vol_ice,icefrac,eff_porosity,&
!                       porsl,hksati,bsw,psi0,rootr,&
!                       zwt,dwat,qcharge)
!
!!-----------------------------------------------------------------------
!! Original author : Yongjiu Dai, 09/1999, 04/2014, 07/2014
!!
!! some new parameterization are added, which are based on CLM4.5
!!
!! Soil moisture is predicted from a 10-layer model (as with soil
!! temperature), in which the vertical soil moisture transport is governed
!! by infiltration, runoff, gradient diffusion, gravity, and root
!! extraction through canopy transpiration. The net water applied to the
!! surface layer is the snowmelt plus precipitation plus the throughfall
!! of canopy dew minus surface runoff and evaporation.
!!
!! The vertical water flow in an unsaturated porous media is described by
!! Darcy's law, and the hydraulic conductivity and the soil negative
!! potential vary with soil water content and soil texture based on the work
!! of Clapp and Hornberger (1978) and Cosby et al. (1984). The equation is
!! integrated over the layer thickness, in which the time rate of change in
!! water mass must equal the net flow across the bounding interface, plus the
!! rate of internal source or sink. The terms of water flow across the layer
!! interfaces are linearly expanded by using first-order Taylor expansion.
!! The equations result in a tridiagonal system equation.
!!
!! Note: length units here are all millimeter
!! (in temperature subroutine uses same soil layer
!! structure required but lengths are m)
!!
!! Richards equation:
!!
!! d wat     d     d psi
!! ----- =  -- [ k(----- - 1) ] + S
!!   dt     dz       dz
!!
!! where: wat = volume of water per volume of soil (mm**3/mm**3)
!! psi = soil matrix potential (mm)
!! dt  = time step (s)
!! z   = depth (mm) (positive downward)
!! dz  = thickness (mm)
!! qin = inflow at top (mm h2o /s)
!! qout= outflow at bottom (mm h2o /s)
!! s   = source/sink flux (mm h2o /s)
!! k   = hydraulic conductivity (mm h2o /s)
!!
!!                       d qin                  d qin
!! qin[n+1] = qin[n] +  --------  d wat(j-1) + --------- d wat(j)
!!                       d wat(j-1)             d wat(j)
!!                ==================|=================
!!                                  < qin
!!
!!                 d wat(j)/dt * dz = qin[n+1] - qout[n+1] + S(j)
!!
!!                                  > qout
!!                ==================|=================
!!                        d qout               d qout
!! qout[n+1] = qout[n] + --------- d wat(j) + --------- d wat(j+1)
!!                        d wat(j)             d wat(j+1)
!!
!!
!! Solution: linearize k and psi about d wat and use tridiagonal
!! system of equations to solve for d wat,
!! where for layer j
!!
!!
!! r_j = a_j [d wat_j-1] + b_j [d wat_j] + c_j [d wat_j+1]
!!
!!-----------------------------------------------------------------------
!    use precision
!    use PhysicalConstants , only : grav,hfus,tfrz,denh2o,denice
!    use steps
!!    use const, only : hlat
!    IMPLICIT NONE
!
!    integer :: nl_soil  ! number of soil layers
!    real :: deltim  ! land model time step (sec)
!    real :: wimp    ! water impremeable if porosity less than wimp
!    real :: smpmin  ! restriction for min of soil potential (mm)
!
!    real :: qinfl   ! infiltration (mm H2O /s)
!    real :: etr     ! vegetation transpiration (mm H2O/s) (+ = to atm)
!
!    real :: z_soisno (1:nl_soil) ! layer depth (m)
!    real :: dz_soisno(1:nl_soil) ! layer thickness (m)
!    real :: zi_soisno(0:nl_soil) ! interface level below a "z" level (m)
!
!    real :: t_soisno (1:nl_soil) ! soil temperature (Kelvin)
!    real :: vol_liq  (1:nl_soil) ! liquid volumetric water content
!    real :: vol_ice  (1:nl_soil) ! ice volumetric water content
!    real :: icefrac  (1:nl_soil)
!    real :: eff_porosity(1:nl_soil) ! effective porosity = porosity - vol_ice
!
!    real :: porsl  (1:nl_soil) ! volumetric soil water at saturation (porosity)
!    real :: hksati (1:nl_soil) ! hydraulic conductivity at saturation (mm H2O /s)
!    real :: bsw    (1:nl_soil) ! Clapp and Hornberger "b"
!    real :: psi0   (1:nl_soil) ! minimum soil suction (mm) [-]
!    real :: rootr  (1:nl_soil) ! effective fraction of roots in each soil layer
!    real :: zwt                ! the depth from ground (soil) surface to water table [m]
!
!    real :: dwat(1:nl_soil)   ! change of soil water [m3/m3]
!    real :: qcharge           ! aquifer recharge rate (positive to aquifer) (mm/s)
!!
!! local arguments
!!
!    integer  :: j,i               ! do loop indices
!    real :: amx(1:nl_soil)    ! "a" left off diagonal of tridiagonal matrix
!    real :: bmx(1:nl_soil)    ! "b" diagonal column for tridiagonal matrix
!    real :: cmx(1:nl_soil)    ! "c" right off diagonal tridiagonal matrix
!    real :: rmx(1:nl_soil)    ! "r" forcing term of tridiagonal matrix
!    real :: zmm(1:nl_soil)    ! layer depth [mm]
!    real :: dzmm(1:nl_soil)   ! layer thickness [mm]
!    real :: zimm(0:nl_soil)   ! layer interface depth [mm]
!    real :: den(1:nl_soil)    ! used in calculating qin, qout
!    real :: alpha(1:nl_soil)  ! used in calculating qin, qout
!    real :: qin(1:nl_soil)    ! flux of water into soil layer [mm h2o/s]
!    real :: qout(1:nl_soil)   ! flux of water out of soil layer [mm h2o/s]
!    real :: dqidw0(1:nl_soil) ! d(qin)/d(vol_liq(j-1))
!    real :: dqidw1(1:nl_soil) ! d(qin)/d(vol_liq(j))
!    real :: dqodw1(1:nl_soil) ! d(qout)/d(vol_liq(j))
!    real :: dqodw2(1:nl_soil) ! d(qout)/d(vol_liq(j+1))
!    real :: dsmpdw(1:nl_soil) ! d(smp)/d(vol_liq)
!    real :: s_node            ! soil wetness
!    real :: s1                ! "s" at interface of layer
!    real :: s2                ! k*s**(2b+2)
!    real :: smp(1:nl_soil)    ! soil matrix potential [mm]
!    real :: hk(1:nl_soil)     ! hydraulic conductivity [mm h2o/s]
!    real :: dhkdw1(1:nl_soil) ! d(hk)/d(vol_liq(j))
!    real :: dhkdw2(1:nl_soil) ! d(hk)/d(vol_liq(j+1))
!    real :: imped(1:nl_soil)  !
!    real :: errorw            ! mass balance error for this time step
!
!    integer  :: jwt               ! index of the soil layer right above the water table (-)
!
!    real, parameter :: e_ice=6.0      !soil ice impedance factor
!
!!-----------------------------------------------------------------------
!    ! Novas variáveis para as mudanças da redistribuição da água
!    real :: dtran,hr_term(1:nl_soil),srcon(1:nl_soil), mpot50,hrb
!    real :: rootx, mpot(1:nl_soil)
!    real :: vol
!!-----------------------------------------------------------------------
!
!    !compute jwt index
!    ! The layer index of the first unsaturated layer,
!    ! i.e., the layer right above the water table
!
!    jwt = nl_soil
!    ! allow jwt to equal zero when zwt is in top layer
!    do j = 1, nl_soil
!       if(zwt <= zi_soisno(j)) then
!          jwt = j-1
!          exit
!       end if
!    enddo
!
!    ! Because the depths in this routine are in mm, use local
!    ! variable arrays instead of pointers
!    do j = 1, nl_soil
!       zmm(j) = z_soisno(j)*1000.
!       dzmm(j) = dz_soisno(j)*1000.
!       zimm(j) = zi_soisno(j)*1000.
!    end do
!
!    zimm(0) = 0.0
!
!    ! Compute matric potential and derivative based on liquid water content only
!    do j = 1, nl_soil
!       if(t_soisno(j)>tfrz) then
!          if(porsl(j)<1.e-6)then     ! bed rock
!             s_node = 0.001
!             smp(j) = psi0(j)
!             dsmpdw(j) = 0.
!          else
!             s_node = max(vol_liq(j)/porsl(j),0.01)
!             s_node = min(1.0,s_node)
!             smp(j) = psi0(j)*s_node**(-bsw(j))
!             smp(j) = max(smpmin,smp(j))
!             dsmpdw(j) = -bsw(j)*smp(j)/(s_node*porsl(j))
!          endif
!       else
!          ! when ice is present, the matric potential is only related to temperature
!          ! by (Fuchs et al., 1978: Soil Sci. Soc. Amer. J. 42(3):379-385)
!          ! Unit 1 Joule = 1 (kg m2/s2), J/kg /(m/s2) ==> m ==> 1e3 mm
!          smp(j) = 1.e3 * 0.3336e6/9.80616*(t_soisno(j)-tfrz)/t_soisno(j)
!          smp(j) = max(smpmin, smp(j))        ! Limit soil suction
!          dsmpdw(j) = 0.
!       endif
!    end do
!
!    ! Hydraulic conductivity and soil matric potential and their derivatives
!    do j = 1, nl_soil
!
!       if(j < nl_soil)then
!          den(j) = (zmm(j+1)-zmm(j))
!          alpha(j) = (smp(j+1)-smp(j))/den(j) - 1.
!       else
!          den(j) = 0.        ! not used
!          alpha(j) = 0.      ! not used
!       endif
!
!       if((eff_porosity(j) < wimp) .OR. (eff_porosity(min(nl_soil,j+1)) < wimp) &
!                                   .OR. (vol_liq(j) <= 0.001))then   !valor limite original é 1.e-3
!          imped(j) = 0.
!          hk(j) = 0.
!          dhkdw1(j) = 0.
!          dhkdw2(j) = 0.
!       else
!          ! The average conductivity between two heterogeneous medium layers (j and j + 1),
!          ! are computed using different methods
!          if(j < nl_soil)then
!! Method I: UPSTREAM MEAN
!             if(alpha(j) <= 0.)then
!                hk(j) = hksati(j) * (vol_liq(j)/porsl(j))**(2.*bsw(j)+3.)
!                dhkdw1(j) = hksati(j) * (2.*bsw(j)+3.)*(vol_liq(j)/porsl(j))**(2.*bsw(j)+2.)/porsl(j)
!                dhkdw2(j) = 0.
!             else
!                hk(j) = hksati(j+1) * (vol_liq(j+1)/porsl(j+1))**(2.*bsw(j+1)+3.)
!                dhkdw1(j) = 0.
!                dhkdw2(j) = hksati(j+1) * (2.*bsw(j+1)+3.)*(vol_liq(j+1)/porsl(j+1))**(2.*bsw(j+1)+2.)/porsl(j+1)
!             endif
!          else
!             hk(j) = hksati(j) * (vol_liq(j)/porsl(j))**(2.*bsw(j)+3.)
!             dhkdw1(j) = hksati(j) * (2.*bsw(j)+3.)*(vol_liq(j)/porsl(j))**(2.*bsw(j)+2.)/porsl(j)
!             dhkdw2(j) = 0.
!          endif
!
!          ! replace fracice with impedance factor
!          imped(j)=10.**(-e_ice*(0.5*(icefrac(j)+icefrac(min(nl_soil,j+1)))))
!          hk(j) = imped(j) * hk(j)
!          dhkdw1(j) = imped(j) * dhkdw1(j)
!          dhkdw2(j) = imped(j) * dhkdw2(j)
!       endif
!    end do
!
!
!! Trecho do HR code from Baker in Amazonia paper
!!...baker HR
!!... rhv nessa subroutina foi eliminado o termo de transpiração devdo a que o
!! .. SiB2 já calculou antes e extraiu a transpiração das camadas do solo
!! .. aqui apenas entra para determinar se é feito o HR
!
!	mpot50 = -1.0 ! baker sibtype.F90:418:
!	hrb = 3.22    ! baker sibtype.F90:419:   ! empirical constant
!
!	if(etr/ dtt  > 1.0 ) then
!	  dtran = 0.0
!	else
!	  dtran = 1.0
!	endif
!
!    hr_term =0.0
!    vol = 0.0
!    srcon = 0.0
!    mpot = 0.0
!!
!!     DO i=1,nl_soil
!!        hr_term(i) = 0.0
!!        mpot(i)    = psi0(i)* 0.0098 *(vol_liq(i) / porsl(i) ) ** (- bsw(i))
!!        srcon(i)   = 1.0 / ( 1.0 +  (mpot(i) / mpot50) ** hrb)
!!        vol = vol + vol_liq(i) !/ porsl(i))/nl_soil
!!         do j = 1,nl_soil
!!             rootx = 0.0
!!            if(j .ne. i ) then
!! 	    if(rootr(i) > rootr(j) ) then
!! 	      rootx = rootr(i)
!! 	    else
!! 	      rootx = rootr(j)
!! 	    endif
!! 	    if(rootx .gt.0.01) hr_term(i) = hr_term(i) +              &
!! 					    ( (mpot(j) - mpot(i))*    &
!! 		                           MAX(srcon(i),srcon(j))*    &
!! 		                          (rootr(i)*rootr(j))/(rootx)*dtran)
!!
!!
!! 	  endif
!! 	  enddo
!!         hr_term(i) = hr_term(i) * 0.002778        ! convert units to mm/sec
!!     enddo
!!...baker HR
!     hr_term =0.0
!    ! Node j=1 (top)
!
!    j = 1
!    qin(j) = qinfl
!    qout(j) = -hk(j)*alpha(j)
!    dqodw1(j) = -(alpha(j)*dhkdw1(j) - hk(j)*dsmpdw(j)/den(j))
!    dqodw2(j) = -(alpha(j)*dhkdw2(j) + hk(j)*dsmpdw(j+1)/den(j))
!
!!itb_HR...
!    if(hr_term(j) > 0.0 ) then
!      qin(j) = qin(j) + hr_term(j)
!    else
!      qout(j) = qout(j) - hr_term(j)
!    endif
!
!    amx(j) = 0.
!    bmx(j) = dzmm(j)/deltim + dqodw1(j)
!    cmx(j) = dqodw2(j)
!    rmx(j) = qin(j) - qout(j) !- etr * rootr(j)
!
!    ! Nodes j=2 to j=nl_soil-1
!
!    do j = 2, nl_soil - 1
!       qin(j) = -hk(j-1)*alpha(j-1)
!       dqidw0(j) = -(alpha(j-1)*dhkdw1(j-1) - hk(j-1)*dsmpdw(j-1)/den(j-1))
!       dqidw1(j) = -(alpha(j-1)*dhkdw2(j-1) + hk(j-1)*dsmpdw(j)/den(j-1))
!
!       qout(j) = -hk(j)*alpha(j)
!
!          if(hr_term(j) > 0.0 ) then
!	    qin(j) = qin(j) + hr_term(j)
!	  else
!	    qout(j) = qout(j) - hr_term(j)
!	  endif
!
!       dqodw1(j) = -(alpha(j)*dhkdw1(j) - hk(j)*dsmpdw(j)/den(j))
!       dqodw2(j) = -(alpha(j)*dhkdw2(j) + hk(j)*dsmpdw(j+1)/den(j))
!
!       amx(j) = -dqidw0(j)
!       bmx(j) =  dzmm(j)/deltim - dqidw1(j) + dqodw1(j)
!       cmx(j) =  dqodw2(j)
!       rmx(j) =  qin(j) - qout(j)! - etr*rootr(j)
!    end do
!
!    ! Node j=nl_soil (bottom)
!    j = nl_soil
!    qin(j) = -hk(j-1)*alpha(j-1)
!    dqidw0(j) = -(alpha(j-1)*dhkdw1(j-1) - hk(j-1)*dsmpdw(j-1)/den(j-1))
!    dqidw1(j) = -(alpha(j-1)*dhkdw2(j-1) + hk(j-1)*dsmpdw(j)/den(j-1))
!
!       qout(j) = hk(j)
!          if(hr_term(j) > 0.0 ) then
!	    qin(j) = qin(j) + hr_term(j)
!	  else
!	    qout(j) = qout(j) - hr_term(j)
!	  endif
!       dqodw1(j) = dhkdw1(j)
!       dqodw2(j) = 0.
!
!    amx(j) = -dqidw0(j)
!    bmx(j) =  dzmm(j)/deltim - dqidw1(j) + dqodw1(j)
!    cmx(j) =  dqodw2(j)
!    rmx(j) =  qin(j) - qout(j) !- etr*rootr(j)
!
!    ! Solve for dwat
!!      if(iter .lt. 10) print*,'OUT  ', (rmx(i),i=1,nl_soil)
!    call tridia (nl_soil, amx, bmx, cmx, rmx, dwat )
!!     write(654,*) iter, (rmx(i)*dtt,i=1,11)
!!     write(654,*) iter, (dwat(i) * dzmm(i),i=1,11)
!! The mass balance error (mm) for this time step is
!    errorw = -deltim*(qin(1)-qout(nl_soil)-dqodw1(nl_soil)*dwat(nl_soil))
!    do j = 1, nl_soil
!       errorw = errorw+dwat(j)*dzmm(j)/1000!+etr*rootr(j)*deltim
!    enddo
!
!     if(abs(errorw) > 1.e-3)then
!       write(*,*) iter,'mass balance error in time step =',errorw
!     endif
!
!
!    ! Recharge rate qcharge to groundwater (positive to aquifer)
!    qcharge = qout(nl_soil) + dqodw1(nl_soil)*dwat(nl_soil)
!   ! if(iter .lt. 10) print*,'OUT  ', (dwat(i),i=1,nl_soil)
!
!  end subroutine soilwater
!
  
  

subroutine soilwater_sib3(nsoil,hltm,vol_liq,satco,poros,phsat, &
			   eff_porosity,bee,td,ect,             &
			  infil,rootf, dzmm, zmm,               &
			  qcharge, dw_liq)

!----------------------------------------------------------------------
!
!     subroutine based on  CLM_SOILWATER
!
!     CLM web info: http://clm.gsfc.nasa.gov
!
!     Description: (taken directly from CLM_SOILWATER.F90)
!     Soil moisture is predicted from a 10-layer model (as with soil
!     temperature), in which the vertical soil moisture transport is 
!     governed by infiltration, runoff, gradient diffusion, gravity 
!     and root extraction through canopy transpiration.  The net water
!     applied to the surface layer is the snowmelt plus precipitation 
!     plus the throughfall of canopy dew minus surface runoff and 
!     evaporation. 
!
!     The vertical water flow in an unsaturated porous media is 
!     described by Darcy's Law, and the hydraulic conductivity and the
!     soil negative potential vary with soil water content and soil 
!     texture based on the work of Clapp and Hornberger (1978) and Cosby 
!     et al. (1984). The equation is integrated over the layer thickness,
!     in which the time rate of change in water mass must equal the net 
!     flow across the bounding interface, plus the rate of internal source
!     or sink.  The terms of water flow across the layer interfaces are 
!     linearly expanded by using first-order Taylor expansion.  The 
!     equations result in a tridiagonal system of equations.

!     Note: lentgh units here are all millimeter
!
!     Richards Equation
!
!     d wat     d      d wat   d psi
!     ----- = - -- [ k(-----   ----- - 1) ] + S
!      dt       dz       dz    d wat
!
!     where:
!     wat = volumetric water content (m^3/m^3)
!     psi = soil matric potential (mm)
!     dt  = time step (sec)
!     z   = depth (mm)
!     qin = inflow at top of layer (mm H2O/sec)
!     qout= outflow at bottom of layer (mm H2O/sec)
!     s   = source/sink flux (mm H2O/sec)
!     k   = hydraulic conductivity (mm H2O/sec)
!
!     Solution:
!     linearize k and psi about d wat and use tridiagonal system of 
!     equations to solve for d wat, where for layer j
!
!     r_j = a_j[d wat_j-1] + b_j [d wat_j] + c_j [d wat_j+1]
!
!     Revision History
!     15 September 1999: Yongjiu Dai; Initial code
!     15 December  1999: Paul Houser and Jon Radakovich; F90 Revision
!     25 January   2002: Ian Baker; SiB integration
!----------------------------------------------------------------------
use precision
! use sibtype
! 
! use physical_parameters, only:  &
!     grav,  &
!     tice,  &
!     hltm
! 
! use sib_const_module, only:  &
!     nsoil,  &
!     wimp,   &
!     phmin,  &
!     dti,    &
!     dtt
!     
!     use multicamada 
     use steps , only : dtt , iter
     use atmos, only : swdown
     use multicamada, only : hr_proc 
implicit none

!----------------------------------------------------------------------

! type(sib_t), intent(inout) :: sib

!----------------------------------------------------------------------  

!...INPUT VARIABLES
real  :: infil   ! infiltration into soil
                                             !  (kg m^-2 sec^-1) 
real  ::dzmm(1:nsoil)     ! espessura das camadas [mm] 
real  ::zmm(1:nsoil) ! profundidade das cama [mm]


!...OUTPUT VARIABLES
real :: hk(1:nsoil)                     ! hydraulic conductivity (mm H2O sec^-1)  
real :: satco(1:nsoil)                   ! hydraulic conductivity  (mm H2O sec^-1)  
					     ! satco = (satco*exp(-decay*depth))*1000
real :: dhkdw(1:nsoil)                   ! d(hk)/d(water)     
real :: dw_liq(1:nsoil)                  ! change in layer liquid  (m^3/m^3)(volumetric)


!...local variables...
integer ::   j     ! loop variables

real  :: hltmi,hltm     ! 1/hltm
real  :: s_node    ! volumetric wetness of node
real  :: s1        ! wetness at interface of layer
real  :: s2        ! conductivity*wetness**(2b+2)
real  :: vol_liq(1:nsoil)       ! layer liquid content (m^3/m^3)
real  :: smp(1:nsoil)           ! soil matric potential (mm)
real  :: dsmpdw(1:nsoil)        ! d(smp)/d(wetness)
real  :: qin       ! flux of water into layer 
                                    !  (mm H2O sec^-1)
real  :: qout      ! flux of water out of layer 
                                    !  (mm H2O sec^-1)  
real  :: den       ! used in calculating qin,qout   
real  :: num       ! used in calculating qin,qout
real  :: dqidw0    ! d(qin)/d(vol_liq(j-1)) 
real  :: dqidw1    ! d(qin)/d(vol_liq(j))
real  :: dqodw1    ! d(qout)/d(vol_liq(j))
real  :: dqodw2    ! d(qout)/d(vol_liq(j+1))
real  :: rmx(1:nsoil), tmx(1:nsoil)
                                    ! "r" forcing term of tridiag matrix
real  :: amx(1:nsoil)
                                    ! "a" left off diagonal of tridiag mat
real  :: bmx(1:nsoil)
                                    ! "b" diagonal column for tridiag mat
real  :: cmx(1:nsoil)
                                    ! "c" right off diagonal of tridiag mat

!itb_HR...
integer :: i,nsoil   ! loop variable

real :: w_in,w_out,bal,wtemp, dtran,grav,mpot50,phmin, rootx,tice,wimp
real :: bee,hrb,phsat,poros,td, ect,qcharge,dti
real :: eff_porosity(1:nsoil)  !porosidade efetiva 
real :: hr_term(1:nsoil)
real :: mpot(1:nsoil),srcon(1:nsoil),rootf(1:nsoil)
!----------------------------------------------------------------------

    if(iter .lt. 10)print*, iter, (vol_liq(i),i=1,nsoil)
    
    hltmi = 1.0/hltm
    dti = 1.0/dtt
    grav = 9.800
    phmin = -1.e8
    wimp = 0.05
    tice = 273.15
    mpot50 = -1.0
    hrb = 3.22
    
    print*,"ect não é usado nessa versão",ect
    
    ! Inicalização de vetores ::
    hk     = 0.0
    dhkdw  = 0.0
    smp    = 0.0
    dsmpdw = 0.0 
    !...set hydraulic conductivity to zero if effective porosity 5% in 
    !...any two adjoining layers, or if volumetric water (liquid) content
    !...less than 0.001
    do j=1,nsoil
!         if( ((sib%diag%eff_poros(j) < wimp)            .or.   &
	if( (eff_porosity(j) .lt. wimp ) .or.                 &
!             (sib%diag%eff_poros(min(nsoil,j+1)) < wimp)) .or. &
              ( eff_porosity(min(j+1,nsoil)) .lt. wimp ) .or.   &
!             (sib%prog%vol_liq(j) <= 1.E-3)) then
               (vol_liq(j) .le. 1.e-3) ) then
            hk(j)    = 0.0
            dhkdw(j) = 0.0
        else
!             s1 = 0.5*(sib%prog%vol_liq(j)+sib%prog%vol_liq(min(nsoil,j+1)))/ &
!                 sib%param%poros
	      s1 = 0.5*(vol_liq(j)+vol_liq(min(nsoil,j+1)))/ &
	           poros
!             s2 = sib%param%satco*1000.0*s1**(2.0*sib%param%bee+2.0)
	      s2 = satco(j)*1000.0*s1**(2.0*bee+2.0)
	      hk(j) = s1*s2
!             dhkdw(j) = (2.0*sib%param%bee+3.0)*s2*0.5/sib%param%poros
	      dhkdw(j) = (2.0*bee+3.0)*s2*0.5/poros
            !itb...I don't like the 0.5 factor in the CLM code-don't understand.
            !itb...it's in the CLM manual-CLM does not follow Bonan exactly.
            if(j==nsoil) dhkdw(j) = dhkdw(j)*2.0
        endif

        !...evaluate hydraulic conductivity, soil matric potential,
        !...d(smp)/d(vol_liq) and d(hk)/d(vol_liq)
!         if(sib%prog%td(j) > tice) then
          if(td > tice) then
!             s_node = max(sib%prog%vol_liq(j)/sib%param%poros,0.01_dbl_kind)
	    s_node = max(vol_liq(j)/poros,0.01)
            s_node = min(1.0,s_node)
!             smp(j) = sib%param%phsat*1000.0*s_node**(-sib%param%bee)
	    smp(j) = phsat*1000.0*s_node**(-bee)
            smp(j) = max(phmin,smp(j))

!             dsmpdw(j) = -sib%param%bee*smp(j)/(s_node*sib%param%poros)
	    dsmpdw(j) = -bee*smp(j)/(s_node*poros)
        else

        !...when ice is present, the matric potential is only related to
        !...temperature by (Fuchs et. al., 1978; Soil Sci. Soc. of Amer. J.
        !...42(3); 379-385) 
        !...Unit 1 joule = 1kg m2/sec2 j/kg/(m/s2) ==> m ==>1e3mm

!           smp(j) = 1.e3*0.3336e6_dbl_kind/grav*(sib%prog%td(j)-tice)/sib%prog%td(j)
	    smp(j) = 1.e3*0.3336e6/grav*(td-tice)/td
	    smp(j) = max(phmin,smp(j))
            dsmpdw(j) = 0.0

        endif
    enddo
  
!itb_HR...set up Hydraulic Redistribution factors...
!itb_HR    (only works when there is no transpiration)

!     if(sib%diag%ect > 0.0 ) then
      if(swdown > 1.0 ) then
!        sib%diag%dtran = 0.0_dbl_kind
       dtran = 0.0
      else
!      sib%diag%dtran = 1.0_dbl_kind
       dtran = 1.0
      endif


      hr_term = 0.0
      mpot    = 0.0
      srcon   = 0.0
!     hr_loop: do i=1,nsoil
      do i=1,nsoil
!        sib%diag%hr_term(i) = 0.0_dbl_kind
        hr_term(i) = 0.0
!        sib%diag%mpot(i) = sib%param%phsat * 0.0098_dbl_kind *      &
!            ( sib%prog%vol_liq(i) / sib%param%poros ) ** (- sib%param%bee)
        mpot(i) = phsat * 0.0098 * ( vol_liq(i) / poros ) ** (- bee)
!        sib%diag%srcon(i) = 1.0_dbl_kind / ( 1.0_dbl_kind +    &
!                (sib%diag%mpot(i) / sib%diag%mpot50) ** sib%diag%hrb)
        srcon(i) = 1.0 / ( 1.0 + (mpot(i) / mpot50) ** hrb)
!         hr_inside: do j = 1,nsoil
         do j = 1,nsoil
!            if(j == i) cycle hr_inside
	   if(j .ne. i)then 
	     rootx = 0.0
! 	    if(sib%param%rootf(i) > sib%param%rootf(j) ) then 
	    if(rootf(i) > rootf(j) ) then 
! 	      sib%diag%rootx = sib%param%rootf(i)
	      rootx = rootf(i)
	    else
! 	      sib%diag%rootx = sib%param%rootf(j)
	      rootx = rootf(j)
	    endif
! 	    sib%diag%hr_term(i) = sib%diag%hr_term(i) + ( (sib%diag%mpot(j) - sib%diag%mpot(i))     &
! 		    * MAX(sib%diag%srcon(i),sib%diag%srcon(j)) *      &
! 		  (sib%param%rootf(i) * sib%param%rootf(j))/sib%diag%rootx * sib%diag%dtran)
	    if(rootx .gt. 0.0) then
	      hr_term(i) = hr_term(i) + ( (mpot(j) - mpot(i)) &
		           * MAX(srcon(i),srcon(j)) *      &
		           (rootf(i) * rootf(j))/rootx *     & 
		           dtran)
	    endif	  
           endif
        enddo 
!         sib%diag%hr_term(i) = sib%diag%hr_term(i) * 0.002778_dbl_kind  ! convert units to mm/sec
        hr_term(i) = hr_term(i) * 0.002778  ! convert units to mm/sec
!       print'(a,i3,3g15.5)','HR:',i,sib%diag%mpot(i),sib%diag%srcon(i),sib%diag%hr_term(i)
    enddo! hr_loop
     if(hr_proc .eq. 0) hr_term = 0.0
     
    !...set up a, b, c and r vectors for tridiagonal solver
    !...node j=1

    j = 1
    qin    = infil
    den    = zmm(j+1) - zmm(j)
    num    = (smp(j+1)-smp(j)) - den
    qout   = -hk(j)*num/den

!itb_HR...
!     if(sib%diag%hr_term(j) > 0.0 ) then
     if(hr_term(j) > 0.0 ) then
!       qin = qin + sib%diag%hr_term(j)
        qin = qin + hr_term(j)
    else
!       qout = qout - sib%diag%hr_term(j)
      qout = qout - hr_term(j)
    endif

!    qout = MIN(qout,vol_liq(j)*dti*0.1)
!itb_HR...restrict water removal to 1/2 the total 
!itb_HR...water in a layer for a single timestep
!     qout = MIN(qout,sib%prog%www_liq(j)*dti*0.5_dbl_kind)
!      qout = MIN(qout,vol_liq(j)*dti*dzmm(j)*0.01)
!itb_HR...

    dqodw1 = -(-hk(j)*dsmpdw(j)   + num*dhkdw(j))/den
    dqodw2 = -( hk(j)*dsmpdw(j+1) + num*dhkdw(j))/den
     rmx(j) = qin - qout  !+ ( ect * rootf(j) * dti * hltmi)
!    rmx(j) = qin - qout 
    amx(j) = 0.0
    bmx(j) = dzmm(j) * (dti) + dqodw1
    cmx(j) = dqodw2

    !...nodes 2 through nsoil-1
    do j = 2,nsoil-1 
        den    = zmm(j) - zmm(j-1)
        num    = (smp(j)-smp(j-1)) - den
        qin    = -hk(j-1)*num/den
        dqidw0 = -(-hk(j-1)*dsmpdw(j-1) + num*dhkdw(j-1))/den
        dqidw1 = -( hk(j-1)*dsmpdw(j)   + num*dhkdw(j-1))/den  
        den    = zmm(j+1)-zmm(j)
        num    = smp(j+1)-smp(j) - den
        qout   = -hk(j)*num/den

!itb_HR...
!     if(sib%diag%hr_term(j) > 0.0 ) then
     if(hr_term(j) > 0.0 ) then
!       qin = qin + sib%diag%hr_term(j)
        qin = qin + hr_term(j)
    else
!       qout = qout - sib%diag%hr_term(j)
      qout = qout - hr_term(j)
    endif
          ! qout = MIN(qout,vol_liq(j)*dti*0.1)
!itb_HR...restrict water removal to 1/2 the total 
!itb_HR...water in a layer for a single timestep
!    qout = MIN(qout,sib%prog%www_liq(j)*dti*0.5_dbl_kind)
!      qout = MIN(qout,vol_liq(j)*dti*dzmm(j)*0.01)
!itb_HR...

        dqodw1 = -(-hk(j)*dsmpdw(j)  +  num*dhkdw(j))/den
        dqodw2 = -( hk(j)*dsmpdw(j+1) + num*dhkdw(j))/den
        rmx(j) = qin - qout  !+ (ect * dti * rootf(j) * hltmi)
        amx(j) = -dqidw0
        bmx(j) = dzmm(j)*dti - dqidw1 + dqodw1
        cmx(j) = dqodw2
    enddo


    !...node j=nsoil
    j = nsoil
    den    = zmm(j) - zmm(j-1)
    num    = smp(j) - smp(j-1) - den
    qin    = -hk(j-1) * num/den
    dqidw0 = -(-hk(j-1)*dsmpdw(j-1) + num*dhkdw(j-1))/den
    dqidw1 = -( hk(j-1)*dsmpdw(j)   + num*dhkdw(j-1))/den
    qout   = hk(j)

!itb_HR...
!     if(sib%diag%hr_term(j) > 0.0 ) then
     if(hr_term(j) > 0.0 ) then
!       qin = qin + sib%diag%hr_term(j)
        qin = qin + hr_term(j)
    else
!       qout = qout - sib%diag%hr_term(j)
      qout = qout - hr_term(j)
    endif

!itb_HR...restrict water removal to 1/2 the total 
!itb_HR...water in a layer for a single timestep
!     qout = MIN(qout,vol_liq(j)*dti*0.5)
!       qout = MIN(qout,vol_liq(j)*dti*dzmm(j)*0.01)
!itb_HR...

    dqodw1 = dhkdw(j)
    rmx(j) = qin - qout  !+ (ect * dti * rootf(j) * hltmi)
    amx(j) = -dqidw0
    bmx(j) = dzmm(j)*dti - dqidw1 + dqodw1
    cmx(j) = 0.0

    ! Add qout out of the bottom layer to runoff
!     sib%diag%roff = sib%diag%roff + qout*dtt
    qcharge =  qout*dtt/1000.0
!itb...check water balance
!     w_in  = 0.0_dbl_kind
      w_in  = 0.0
    do j=1,nsoil
!      w_in = w_in + sib%prog%www_liq(j) + sib%prog%www_ice(j)
     w_in = w_in + vol_liq(j) !+ sib%prog%www_ice(j)
    enddo


!itb...balance

      tmx = rmx
    call  tridia (nsoil, amx, bmx, cmx, rmx, dw_liq)

    w_out = 0.0
    wtemp = 0.0

    do j=1,nsoil
!      w_out = w_out + sib%prog%www_liq(j) + sib%prog%www_ice(j) + dw_liq(j)
     w_out = w_out + vol_liq(j) + dw_liq(j)
     wtemp = wtemp + dw_liq(j)
    enddo

    bal = w_in - w_out

    if(ABS(bal) .gt. 0.1) then
     write(*,'(a,4g15.5)')'balance=',bal,infil*dtt,qout*dtt,wtemp
    endif

end subroutine soilwater_sib3

  
