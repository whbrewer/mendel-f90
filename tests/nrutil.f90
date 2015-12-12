! This set of modules is useful for verification tests 

MODULE nrutil
! Numerical Recipes in Fortran 90 Utilities (obtained on-line)
        USE nrtype
        IMPLICIT NONE
        
        PRIVATE
        
        INTEGER(I4B), PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
        INTEGER(I4B), PARAMETER :: NPAR_GEOP=4,NPAR2_GEOP=2
        INTEGER(I4B), PARAMETER :: NPAR_CUMSUM=16
        INTEGER(I4B), PARAMETER :: NPAR_CUMPROD=8
        INTEGER(I4B), PARAMETER :: NPAR_POLY=8
        INTEGER(I4B), PARAMETER :: NPAR_POLYTERM=8
        
        INTERFACE swap
                MODULE PROCEDURE swap_i,swap_r,swap_rv,swap_c, &
                        swap_cv,swap_cm,swap_z,swap_zv,swap_zm, &
                        masked_swap_rs,masked_swap_rv,masked_swap_rm
        END INTERFACE
        
        INTERFACE assert
                MODULE PROCEDURE assert1,assert2,assert3,assert4,assert_v
        END INTERFACE
        
        PUBLIC :: swap,assert

CONTAINS

        SUBROUTINE swap_i(a,b)
        INTEGER(I4B), INTENT(INOUT) :: a,b
        INTEGER(I4B) :: dum
        dum=a
        a=b
        b=dum
        END SUBROUTINE swap_i

        SUBROUTINE swap_r(a,b)
        REAL(SP), INTENT(INOUT) :: a,b
        REAL(SP) :: dum
        dum=a
        a=b
        b=dum
        END SUBROUTINE swap_r

        SUBROUTINE swap_rv(a,b)
        REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
        REAL(SP), DIMENSION(SIZE(a)) :: dum
        dum=a
        a=b
        b=dum
        END SUBROUTINE swap_rv

        SUBROUTINE swap_c(a,b)
        COMPLEX(SPC), INTENT(INOUT) :: a,b
        COMPLEX(SPC) :: dum
        dum=a
        a=b
        b=dum
        END SUBROUTINE swap_c

        SUBROUTINE swap_cv(a,b)
        COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: a,b
        COMPLEX(SPC), DIMENSION(SIZE(a)) :: dum
        dum=a
        a=b
        b=dum
        END SUBROUTINE swap_cv

        SUBROUTINE swap_cm(a,b)
        COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
        COMPLEX(SPC), DIMENSION(size(a,1),size(a,2)) :: dum
        dum=a
        a=b
        b=dum
        END SUBROUTINE swap_cm

        SUBROUTINE swap_z(a,b)
        COMPLEX(DPC), INTENT(INOUT) :: a,b
        COMPLEX(DPC) :: dum
        dum=a
        a=b
        b=dum
        END SUBROUTINE swap_z

        SUBROUTINE swap_zv(a,b)
        COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: a,b
        COMPLEX(DPC), DIMENSION(SIZE(a)) :: dum
        dum=a
        a=b
        b=dum
        END SUBROUTINE swap_zv

        SUBROUTINE swap_zm(a,b)
        COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
        COMPLEX(DPC), DIMENSION(size(a,1),size(a,2)) :: dum
        dum=a
        a=b
        b=dum
        END SUBROUTINE swap_zm

        SUBROUTINE masked_swap_rs(a,b,mask)
        REAL(SP), INTENT(INOUT) :: a,b
        LOGICAL(LGT), INTENT(IN) :: mask
        REAL(SP) :: swp
        if (mask) then
                swp=a
                a=b
                b=swp
        end if
        END SUBROUTINE masked_swap_rs

        SUBROUTINE masked_swap_rv(a,b,mask)
        REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
        LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
        REAL(SP), DIMENSION(size(a)) :: swp
        where (mask)
                swp=a
                a=b
                b=swp
        end where
        END SUBROUTINE masked_swap_rv

        SUBROUTINE masked_swap_rm(a,b,mask)
        REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a,b
        LOGICAL(LGT), DIMENSION(:,:), INTENT(IN) :: mask
        REAL(SP), DIMENSION(size(a,1),size(a,2)) :: swp
        where (mask)
                swp=a
                a=b
                b=swp
        end where
        END SUBROUTINE masked_swap_rm

        SUBROUTINE assert1(n1,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, INTENT(IN) :: n1
        if (.not. n1) then
                write (*,*) 'nrerror: an assertion failed with this tag:', &
                        string
                STOP 'program terminated by assert1'
        end if
        END SUBROUTINE assert1

        SUBROUTINE assert2(n1,n2,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, INTENT(IN) :: n1,n2
        if (.not. (n1 .and. n2)) then
                write (*,*) 'nrerror: an assertion failed with this tag:', &
                        string
                STOP 'program terminated by assert2'
        end if
        END SUBROUTINE assert2

        SUBROUTINE assert3(n1,n2,n3,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, INTENT(IN) :: n1,n2,n3
        if (.not. (n1 .and. n2 .and. n3)) then
                write (*,*) 'nrerror: an assertion failed with this tag:', &
                        string
                STOP 'program terminated by assert3'
        end if
        END SUBROUTINE assert3

        SUBROUTINE assert4(n1,n2,n3,n4,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, INTENT(IN) :: n1,n2,n3,n4
        if (.not. (n1 .and. n2 .and. n3 .and. n4)) then
                write (*,*) 'nrerror: an assertion failed with this tag:', &
                        string
                STOP 'program terminated by assert4'
        end if
        END SUBROUTINE assert4

        SUBROUTINE assert_v(n,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, DIMENSION(:), INTENT(IN) :: n
        if (.not. all(n)) then
                write (*,*) 'nrerror: an assertion failed with this tag:', &
                        string
                STOP 'program terminated by assert_v'
        end if
        END SUBROUTINE assert_v

END MODULE nrutil

