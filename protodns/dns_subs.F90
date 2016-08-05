MODULE DNS_subsMod
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : REAL64
    IMPLICIT NONE
    PRIVATE
CONTAINS
    SUBROUTINE init(mx, my, w)
        IMPLICIT NONE
            INTEGER,                                                        INTENT(in)    :: mx, my
            REAL(KIND = REAL64), DIMENSION(0:(mx - 1), (-my + 1):(my + 1)), INTENT(inout) :: w
            INTEGER                                                                       :: i,j,start

        DO j = (-my + 1), (my + 1)
            IF (j <= 0) THEN
                start = 1
            ELSE
                start = 0
            ENDIF
            DO i = start, (mx - 1)
                w(i,j) = REAL(1, KIND=REAL64)/REAL(i**2 + j**2, KIND=REAL64)
            END DO
        END DO
        RETURN
    END SUBROUTINE init

    SUBROUTINE multadvection2(F, m)
        IMPLICIT NONE

        DO i = lower, upper, stride
            F(i, :) = [F(i,2)**2 - F(i,1)**2, &
                       F(i,1)*F(i,2)]
        END DO
        RETURN
    END SUBROUTINE multadvection2

    SUBROUTINE Source(aruments)
        IMPLICIT NONE

        f0(0,0) = REAL(0, KIND=REAL64)
        f1(0,0) = REAL(0, KIND=REAL64)

        DO j = lower, upper, stride
            DO i = lower, upper, stride
                k2_inv  = REAL(1, KIND=REAL64)/REAL(i**2 + j**2, KIND=REAL64)
                ik2_inv = REAL(i, KIND=REAL64)*k2_inv
                jk2_inv = REAL(j, KIND=REAL64)*k2_inv
                f0(i, j) = COMPLEX()
                f1(i, j) = COMPLEX()
            END DO
        END DO


        DO j = lower, upper, stride
            DO i = lower, upper, stride
                f0(i,j) = i*j*f0(i,j) + (i**2 - j**2)*f1(i,j) - nu*(i**2 + j**j)*w(i,j)
            END DO
        END DO
        RETURN
    END SUBROUTINE Source

    SUBROUTINE Spectrum(aruments)
        IMPLICIT NONE

            
        RETURN
    END SUBROUTINE Spectrum

    SUBROUTINE Output(aruments)
        IMPLICIT NONE
        code
        RETURN
    END SUBROUTINE Output


END MODULE DNS_subsMod
