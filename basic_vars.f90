module basic_vars

   use constants, only : &
      nref_max, &
      nchar_max



   implicit none



   integer :: nrank = 1
   integer :: myrank = 0
   integer :: myunit = 1000
   integer :: root = 0
!!!  integer :: nthread
!!!  integer :: mythread

   integer :: sglob1(0:nref_max)
   integer :: sglob2(0:nref_max)

   ! RUN parameters with their default values
   integer :: nbytes_real = 8
   logical :: times_enabled = .true.
   logical :: downscale = .false.
   character (len = 7) :: sfx = '.bin'
   character (len = nchar_max) :: outdir = '/var/tmp/nelsonvn/zzz_wave'



   namelist /run_parms/ nbytes_real, times_enabled, downscale, sfx, outdir



contains



   subroutine read_config_file
      implicit none
      ! Local variables
      !integer :: k

      ! Read file configuration
      open(unit = 1, file = 'model.par', status = 'old')
      read(1, nml = run_parms)
      !read(1, nml = mg_parms)
      !read(1, nml = grid_parms)
      close(1, status = 'keep')

      write(*, nml = run_parms)
   end subroutine read_config_file



   subroutine init_basic_vars
      implicit none
   end subroutine init_basic_vars



end module basic_vars
