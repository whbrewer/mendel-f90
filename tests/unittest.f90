module unittest

!integer, parameter :: max_del_mutn_per_indiv = 10000
integer, parameter :: imax = 1000, nls = 1000
integer :: lb_modulo

contains

subroutine assert( b )

logical b

if ( .not. b) then
  write(*,*) 'assert failed'
  stop 1
endif

end subroutine assert

subroutine test_init(dmutn,nmutn,fmutn,lb_mutn_count,linkage_block_fitness,max_size)
!include 'common.h'
implicit none
integer, intent(in)    :: max_size
integer, intent(inout) :: dmutn(imax,2,max_size)
integer, intent(inout) :: nmutn(imax,2,max_size)
integer, intent(inout) :: fmutn(imax,2,max_size)
integer, intent(inout) :: lb_mutn_count(nls,2,3,max_size)
real*8,  intent(inout) :: linkage_block_fitness(nls,2,max_size)

lb_modulo  = (2**30-2)/nls
dmutn = nls*lb_modulo + 1
nmutn = nls*lb_modulo + 1
fmutn = nls*lb_modulo + 1
dmutn(1,:,:)  = 0
nmutn(1,:,:)  = 0
fmutn(1,:,:)  = 0
lb_mutn_count = 0
linkage_block_fitness = 1.d0
end subroutine test_init

end module unittest

