MODULE random_pkg

IMPLICIT NONE
REAL, PRIVATE      :: zero = 0.0, half = 0.5, one = 1.0, two = 2.0,   &
                      vsmall = TINY(1.0), vlarge = HUGE(1.0)
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)

CONTAINS

FUNCTION randomnum(iseed) RESULT(fn_val)

!...  This function returns a uniform random deviate between 0.0 and
!...  1.0.  To initialize or reinitialize the sequence, set iseed to
!...  to any negative integer.  This function is adapted from p. 197
!...  of Numerical Recipes by Press et al, 1989.

  real fn_val, rm
  integer ia, ic, ir, iseed, ix, iy, j, m
  parameter (m=714025, ia=1366, ic=150889, rm=1./m)
  common /rndm/ ir(97), ix, iy
  
  if(iseed .lt. 0) then
     
     ix = mod(ic-iseed, m)
     
     do j=1,97
        ix    = mod(ia*ix+ic, m)
        ir(j) = ix
     end do
     
     ix = mod(ia*ix+ic, m)
     iy = ix
     
  end if
  
  j = 1 + (97*iy)/m
  
  if(j.gt.97 .or. j.lt.1) stop
  
  iy    = ir(j)
  ix    = mod(ia*ix+ic, m)
  ir(j) = ix
  
  fn_val = iy*rm
  
END FUNCTION randomnum

FUNCTION random_normal() RESULT(fn_val)

! Adapted from the following Fortran 77 code
!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

!  The function random_normal() returns a normally distributed pseudo-random
!  number with zero mean and unit variance.

!  The algorithm uses the ratio of uniforms method of A.J. Kinderman
!  and J.F. Monahan augmented with quadratic bounding curves.

REAL :: fn_val

!     Local variables
REAL     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,           &
            r1 = 0.27597, r2 = 0.27846, u, v, x, y, q

!     Generate P = (u,v) uniform in rectangle enclosing acceptance region

DO

DO
  u = randomnum(1)
  v = randomnum(1)
  v = 1.7156 * (v - half)

!     Evaluate the quadratic form
  x = u - s
  y = ABS(v) - t
  q = x**2 + y*(a*y - b*x)

!     Accept P if inside inner ellipse
  IF (q < r1) EXIT
!     Reject P if outside outer ellipse
  IF (q > r2) CYCLE
!     Reject P if outside acceptance region
  IF (v**2 < -4.0*LOG(u)*u**2) CYCLE
END DO

!     Return ratio of P's coordinates as the normal deviate
  fn_val = v/u

  u = fn_val/fn_val

  IF (u == u) EXIT

END DO

RETURN

END FUNCTION random_normal

function poisson(poisson_mean) result(deviate)

! This routine returns an integer, deviate, that is a random deviate
! drawn from a Poisson distribution of mean poisson_mean, using the
! function randomnum as a source of uniform random deviates.
! Taken from Press, et al., Numerical Recipes, 1986, pp. 207-208.
implicit none
real poisson_mean, old_mean, pi, em, t, g, sq, lg, y
common /pssn/ old_mean, sq, lg, g
integer deviate
parameter (pi = 3.14159265)
data old_mean /-1./             ! Flag for whether poisson_mean
                                ! has changed since last call.
if(poisson_mean < 12.) then     ! Use direct method.

   if(poisson_mean /= old_mean) then
      old_mean = poisson_mean
      g = exp(-poisson_mean)    ! If mean is new, compute g.
   end if

   em = -1.
   t  =  1.

   do while(t > g)
      em = em + 1.
      t  = t*randomnum(1)
   end do

else                            ! Use rejection method.

   if(poisson_mean /= old_mean) then
      old_mean = poisson_mean
      sq = sqrt(2.*poisson_mean)
      lg = log(poisson_mean)
      g  = poisson_mean*lg - gammln(poisson_mean + 1.)
   end if

 1       continue

   y  = tan(pi*randomnum(1))
   em = sq*y + poisson_mean
   if(em < 0.) go to 1

   em = int(em)
   t  = 0.9*(1. + y**2)*exp(em*lg - gammln(em + 1.) - g)
   if(randomnum(1) > t) go to 1

end if

deviate = em

end function poisson

function gammln(arg)

!   This function returns the value of log(gamma(arg)) for arg > 0., 
!   where gamma is the gamma function.  Full accuracy is obtained
!   for arg > 1.
!   Taken from Press, et al., Numerical Recipes, 1986, p. 157.

implicit none
real gammln, arg
real*8 cof, stp, x, tmp, ser
integer j
common /gmln/ cof(6), stp

data cof, stp /76.18009173d0, -86.50532033d0, 24.01409822d0, &
              -1.231739516d0, 0.120858003d-2,  -0.536382d-5, &
             2.50662827465d0/

x   = arg - 1.d0
tmp =  x + 5.5d0
tmp = (x + 0.5d0)*log(tmp) - tmp
ser = 1.d0

do j=1,6
   x   = x   + 1.d0
   ser = ser + cof(j)/x
end do

gammln = tmp + log(stp*ser)

end function gammln

END MODULE random_pkg

