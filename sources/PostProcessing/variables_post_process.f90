MODULE variables_post_process
!
USE type
!
! Define the different common variables for post-processing
!
INTEGER :: n1,n2
REAL(RP), ALLOCATABLE, DIMENSION(:)   :: x,y
REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: eta, phis
!
! Input file
INTEGER               :: i_ana, i_card, tecplot
REAL(RP)              :: T_start, T_stop, x_min, x_max, y_min, y_max, t_min, t_max
CHARACTER(LEN=100)    :: file_3d, file_mod
!
INTEGER, PARAMETER    :: n_hdr = 34 ! Headerlines in '3d.dat' including line variables
REAL(RP), PARAMETER   :: HfoHs = 2.0_rp ! Freak wave height on Hs threshold for detection


!! For output in VP_card
!REAL(RP) :: x_min, x_max, y_min, y_max, z_min, z_max, T_init
!REAL(RP), ALLOCATABLE, DIMENSION(:)   :: xvect, yvect, zvect
!REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: eta_card
!INTEGER :: i_xvect, i_yvect, i_zvect, imin, imax, jmin, jmax
!
CONTAINS
!
LOGICAL FUNCTION iseven(n)
!
IMPLICIT NONE
!
INTEGER :: n
!
iseven = (MOD(n,2) == 0)
!
END FUNCTION iseven
!
END MODULE variables_post_process