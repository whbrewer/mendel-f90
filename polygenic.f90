MODULE polygenic

implicit none
! this corresponds to MNP in diagnostics_polymorphism_plot()
integer, parameter :: max_polys = 50000
integer, parameter :: A=1, C=2, G=3, T=4
integer :: num_polys_cumulative, num_polys_this_gen, poly_str_len
integer :: poly_gen_first_instance, poly_gen_last_instance
integer :: poly_stop_gen, duration, npoly
integer :: poly_not_selected, poly_fixed_fmutn
logical :: polygenic_fixed
logical, parameter :: poly_debug = .false.
real :: percent_pop_poly

type poly_mutn_table
  integer :: id
  integer :: gen_enter
  integer :: gen_exit
end type
type(poly_mutn_table), allocatable, dimension(:) :: pmutn

contains

character(len=poly_str_len) function poly_to_string(buffer)
character*1 nucl
integer, intent(in) :: buffer(*)
integer :: i
poly_to_string = ''
do i = 1, poly_str_len
  nucl = mutn_to_nucl(buffer(i))
  poly_to_string = trim(poly_to_string)//trim(nucl)
end do
return
end function poly_to_string

! return ASCII code which corresponds to either A,C,G,T
character function mutn_to_nucl(mutn)
integer, intent(in) :: mutn
real :: x
integer :: i
character*1 nucl(4)
nucl = (/ 'A', 'C', 'G', 'T' /)
mutn_to_nucl = nucl(mod(mutn-1,4)+1)
return
end function mutn_to_nucl

! given A,C,G,T return either 1,2,3,4 
integer function nucl_to_int(nucl)
use random_pkg
character*1 nucl
select case (nucl)
   case ('A') 
      nucl_to_int = A
   case ('C') 
      nucl_to_int = C
   case ('G') 
      nucl_to_int = G
   case ('T') 
      nucl_to_int = T
   case default
      ! if not A,T,C,G return random nucleotide
      nucl_to_int = min(4, 1 + int(4*randomnum(1)))
end select 
end function nucl_to_int

character*1 function int_to_nucl(x)
integer, intent(in) :: x
select case (x)
   case (A)
      int_to_nucl = 'A'
   case (C)
      int_to_nucl = 'C'
   case (G)
      int_to_nucl = 'G'
   case (T)
      int_to_nucl = 'T'
end select
end function int_to_nucl

! return ASCII code which corresponds to either A,C,G,T
character function random_nucl(n)
integer, intent(in) :: n
real :: x
integer :: i
character*1 nucl(4)
nucl = (/ 'A', 'C', 'G', 'T' /)
call random_number(x)
i = int(4*x)+1
random_nucl = nucl(i)
return
end function random_nucl

logical function poly_match(buffer)
use inputs
integer, intent(in) :: buffer(*)
! note: if changing the string length below, need to also change in inputs.f90
character*40 str1, str2
str1 = poly_to_string(buffer)
str2 = to_upper(polygenic_target)
if(str1==str2) then
   poly_match = .true.
else
   poly_match = .false.
endif 
return
end function poly_match

function to_upper(strIn) result(strOut)
! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)

     implicit none

     character(len=*), intent(in) :: strIn
     character(len=len(strIn)) :: strOut
     integer :: i,j

     do i = 1, len(strIn)
          j = iachar(strIn(i:i))
          if (j>= iachar("a") .and. j<=iachar("z") ) then
               strOut(i:i) = achar(iachar(strIn(i:i))-32)
          else
               strOut(i:i) = strIn(i:i)
          end if
     end do

end function to_upper

end module polygenic

