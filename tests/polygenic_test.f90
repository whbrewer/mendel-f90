program test
integer a(7)
real x
character*1 nucl
character*2 y
character*20 str

a = (/ 2, 3, 5, 7, 11, 13, 17 /)
!call poly_to_string(a) 
call random_seed()
str = ''
do i = 1, 10
   call random_number(x)
   j = int(100000*x)
   nucl = char(mutn_to_nucl(j))
   print *, j, str//nucl
   str = trim(str) // trim(nucl) 
end do
print '(a)', str
end program test

subroutine poly_to_string(buffer)
integer, intent(in) :: buffer(*)
do i = 1, 7
  print *, i, buffer(i)
end do
end subroutine poly_to_string

! return ASCII code which corresponds to either A,C,G,T
function mutn_to_nucl(mutn)
integer, intent(in) :: mutn
integer :: nucl(4)
integer mutn_to_nucl
nucl = (/ 65, 67, 71, 84 /)
mutn_to_nucl = nucl(mod(mutn,3)+1)
return
end function mutn_to_nucl

function poly_match(str1,str2)


end function poly_match
