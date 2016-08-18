!
! DATE ROUTINE JD(YYYY, MM, DD) CONVERTS CALENDER DATE TO
! JULIAN DATE.  SEE CACM 1968 11(10):657, LETTER TO THE
! EDITOR BY HENRY F. FLIEGEL AND THOMAS C. VAN FLANDERN.
! EXAMPLE JD(1970, 1, 1) = 2440588
!
! Taken from http://jblevins.org/mirror/amiller/datesub.f90
! See http://www.stiltner.org/book/bookcalc.htm, too.
!
function cal2jul(yy, mm, dd) result(res)
   implicit none
   integer, intent(in) :: yy
   integer, intent(in) :: mm
   integer, intent(in) :: dd
   integer :: res

   res = dd - 32075 + &
      1461 * (yy + 4800 + (mm - 14) / 12) / 4 +  &
      367 * (mm - 2 - ((mm - 14) / 12) * 12) / 12 - &
      3 * ((yy + 4900 + (mm - 14) / 12) / 100) / 4
end function cal2jul



FUNCTION calendarday_date(YEAR,MM,DAY,HOUR)result(day_date)

  implicit none
  integer, intent(in) :: YEAR,MM,DAY,HOUR
  real   :: day_date
  integer :: DOY
  integer, external :: cal2jul

  DOY = cal2jul(YEAR,MM,DAY) - cal2jul(YEAR,1,1)  +1
  day_date = float(DOY) + (float(HOUR)*3600.-1.)/86400.

END FUNCTION calendarday_date


!program teste
!   implicit none
!   integer, external :: cal2jul
!   real, external :: calendarday_date
!    integer :: date(3),day_year
!
!    date(1) = 2004
!    date(2) = 11
!    date(3) = 3600*23
!
!   write(*,*) cal2jul(2012, 1, 1)
!   write(*,*) cal2jul(2013, 1, 1)
!   day_year = cal2jul(2013, 3, 3) - cal2jul(2013, 1, 1) + 1
!    write(*,*) day_year, "day_year"
!   write(*,*) calendarday_date(2013,3,3,2)
!
!end program teste

