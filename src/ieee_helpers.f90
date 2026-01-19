module ieee_helpers
use, intrinsic :: ieee_arithmetic
use, intrinsic :: ieee_exceptions

contains

subroutine ieee_init()
   call ieee_set_flag(ieee_invalid, .false.)
   call ieee_set_flag(ieee_divide_by_zero, .false.)
   call ieee_set_flag(ieee_overflow, .false.)
   call ieee_set_flag(ieee_underflow, .false.)
   call ieee_set_flag(ieee_inexact, .false.)
end

subroutine ieee_report(tag)
   character(*), intent(in) :: tag
   logical :: flag

   call ieee_get_flag(ieee_invalid, flag)
   if (flag) write(*,*) 'IEEE flag set at ', trim(tag), ': invalid'

   call ieee_get_flag(ieee_divide_by_zero, flag)
   if (flag) write(*,*) 'IEEE flag set at ', trim(tag), ': divide_by_zero'

   call ieee_get_flag(ieee_overflow, flag)
   if (flag) write(*,*) 'IEEE flag set at ', trim(tag), ': overflow'

   call ieee_get_flag(ieee_underflow, flag)
   if (flag) write(*,*) 'IEEE flag set at ', trim(tag), ': underflow'

   call ieee_get_flag(ieee_inexact, flag)
   if (flag) write(*,*) 'IEEE flag set at ', trim(tag), ': inexact'

   call ieee_set_flag(ieee_invalid, .false.)
   call ieee_set_flag(ieee_divide_by_zero, .false.)
   call ieee_set_flag(ieee_overflow, .false.)
   call ieee_set_flag(ieee_underflow, .false.)
   call ieee_set_flag(ieee_inexact, .false.)
end

end module ieee_helpers
