module constants
  implicit none
  integer, parameter :: r8 = kind(1.0d0)
!!!  integer, parameter :: r8 = kind(1.0)
!!!  integer, parameter :: UNDEF = -100000   ! UNDEF must be always negative
 
  integer, parameter :: nrank_max = 2048
  integer, parameter :: nref_max = 10
  integer, parameter :: nseg_max = 5000
  integer, parameter :: nchar_max = 120

!!!  real (kind = r8), parameter :: pi = 3.14159265358979323846_r8
  real, parameter :: pi = 3.14159265358979323846
  real, parameter :: twopi = 6.28318530717958647692
  real, parameter :: pihalf = 1.57079632679489661923
  real, parameter :: s2r = pi / 180.
  real, parameter :: r2s = 180. / pi
  real, parameter :: decmax = s2r * 23.5

!!!  real (kind = r8), parameter :: a = 6.371229E+06_r8   ! Earth radius in meters
  real, parameter :: g = 9.80616
  real, parameter :: c2k = 273.16

  character (len = nchar_max) :: tmpdir = '/var/tmp/nelsonvn/'
end module constants
