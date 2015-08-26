program test_poisson
use random_pkg
integer*8 num_gen

poisson_method = 1 ! 0=Numerical Recipes, 1=RANLIB

poisson_mean = 1.e-5
num_gen = 1e6
if (poisson_method==1) then
   new_mutn = random_Poisson(poisson_mean,.true.)
end if

do i = 1, num_gen
   if (poisson_method == 1) then
      new_mutn = random_Poisson(poisson_mean,.false.)
   else
      new_mutn = poisson(poisson_mean)
   endif
   if (new_mutn /= 0) print *, i, new_mutn, poisson_mean
end do

end program test_poisson
