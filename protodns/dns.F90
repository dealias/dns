PROGRAM protodns_Fortran
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : REAL64
    IMPLICIT NONE
        INTEGER                                          :: Nx, Ny
        INTEGER                                          :: mx, my
        REAL(KIND = REAL64)                              :: nu, dt
        REAL(KIND = REAL64), ALLOCATABLE, DIMENSION(:)   :: vector
        REAL(KIND = REAL64), ALLOCATABLE, DIMENSION(:,:) :: w, f0, f1
END PROGRAM protodns_Fortran
