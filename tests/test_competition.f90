! Mendel unit testing framework

! Future tests to implement:
! make sure input file is read correctly
! if I set mutation rate = 0, fitness should be zero.

program test_main
use mpi
use mpi_helpers
use random_pkg
use inputs
include 'common.h'
real*8,  allocatable, dimension(:,:,:)   :: linkage_block_fitness
integer, allocatable, dimension(:,:,:)   :: dmutn, nmutn, fmutn
integer, allocatable, dimension(:,:,:,:) :: lb_mutn_count
integer :: max_size, nmax, pop_size_allocation
real :: fraction_selected_away
integer :: i,j,k,pop_size_winner,pop_size_loser,num_migrate
integer, allocatable, dimension(:) :: pop_size_array
logical winner

open (5, file='mendel.in',status='old')
call read_parameters(5)
close(5)
call write_parameters(6)

call mpi_myinit(myid,ierr)

if(myid.eq.0) then
  pop_size = 100 
else
  pop_size = 0
endif

! compute global population size
if(is_parallel) then
  call mpi_isum(pop_size,global_pop_size,1)
  call mpi_mybcasti(global_pop_size,1)
end if

if(is_parallel) then
   allocate( pop_size_array(num_tribes) )
   pop_size_array = pop_size
   pop_size_allocation = global_pop_size
   if (tribal_competition) then
      write(*,*) 'Allocating tribe ', myid, ' with max pop_size of:',&
                   pop_size_allocation
   endif
else
   pop_size_allocation = pop_size
endif

if(tribal_competition) then
   fraction_selected_away = 1. - 1./reproductive_rate
   tribal_fitness_factor = 1.d0
   max_size = int(0.55*12.*pop_size_allocation &
                      *(1. - fraction_random_death))
   nmax = 12.*(1. - fraction_random_death) + 0.999
!  Limit the minimum value of heritability to be 10**-20.
   group_heritability = max(1.e-20, group_heritability)
else
   max_size = int(1.1*reproductive_rate*pop_size_allocation &
                     *(1. - fraction_random_death))
   nmax = 2.*reproductive_rate*(1. - fraction_random_death) + 0.999
end if

allocate(         dmutn(max_del_mutn_per_indiv/2,2,max_size),     &
                  nmutn(max_neu_mutn_per_indiv/2,2,max_size),     &
                  fmutn(max_fav_mutn_per_indiv/2,2,max_size),     &
                 lb_mutn_count(num_linkage_subunits,2,3,max_size),&
         linkage_block_fitness(num_linkage_subunits,2,max_size))

call test_init(dmutn,nmutn,fmutn,lb_mutn_count,linkage_block_fitness,max_size)

if (is_parallel .and. num_tribes > 1) then
   current_pop_size = pop_size
   call mpi_gather(current_pop_size,1,mpi_integer, &
        pop_size_array,1,mpi_integer,0,mycomm,ierr)

   if(myid==0) then
      write(6,*) myid,'competing pop sizes:', pop_size_array
      pop_size_winner = maxval(pop_size_array)
      print *,'** maxloc:',maxloc(pop_size_array)
      pop_size_loser = minval(pop_size_array)
      print *,'** minloc:',minloc(pop_size_array)
      call mpi_send_int(pop_size_winner,1,msg_num,ierr)
      msg_num  = msg_num + 1
      call mpi_send_int(pop_size_loser,1,msg_num,ierr)
      msg_num  = msg_num + 1
   else
      call mpi_recv_int(pop_size_winner,0,msg_num,ierr)
      msg_num  = msg_num + 1
      call mpi_recv_int(pop_size_loser,0,msg_num,ierr)
      msg_num  = msg_num + 1
   end if
   winner = .false.

   if(myid.eq.0) then ! dying tribe
      call mpi_recv_int(current_pop_size,1-myid,msg_num,ierr)
      msg_num  = msg_num + 1
      call mpi_recv_int(num_migrate,1-myid,msg_num,ierr)
      msg_num  = msg_num + 1
   else ! winning tribe... send half of its pop size to dying tribe
      winner = .true.
      write(6,*)'<font color=red>*** FISSION TRIBE ***</font>', myid
      num_migrate = (pop_size_winner - pop_size_loser)/2
      current_pop_size = (pop_size_winner + pop_size_loser)/2
      ! if odd, round down half the time and round up half the time
      if(mod(current_pop_size,2)==1 .and. randomnum(1).gt.0.5) then
         num_migrate = num_migrate+1
         current_pop_size = current_pop_size + 1
      end if
      call mpi_send_int(current_pop_size,1-myid,msg_num,ierr)
      msg_num  = msg_num + 1
      call mpi_send_int(num_migrate,1-myid,msg_num,ierr)
      msg_num  = msg_num + 1
      write(6,*) 'migrating half the tribe from',myid+1,' to ',2-myid
      print *, 'num_migrate:',num_migrate
  end if

  do k = 1, num_migrate
     i = pop_size_winner-num_migrate + k
     j = pop_size_loser + k
     if(myid.eq.0) write(*,*) 'migrating: ',i, 'to:',j
     call migrate_individual(i,j,dmutn,fmutn,nmutn,lb_mutn_count, &
                          linkage_block_fitness,winner)
  end do

  print *, 'myid:',myid,'current_pop_size:',current_pop_size

end if

call mpi_myfinalize(ierr)

end program test_main

subroutine test_init(dmutn,nmutn,fmutn,lb_mutn_count,linkage_block_fitness,max_size)
use inputs
include 'common.h'
integer, intent(in)    :: max_size
integer, intent(inout) :: dmutn(max_del_mutn_per_indiv/2,2,max_size)
integer, intent(inout) :: nmutn(max_neu_mutn_per_indiv/2,2,max_size)
integer, intent(inout) :: fmutn(max_fav_mutn_per_indiv/2,2,max_size)
integer, intent(inout) :: lb_mutn_count(num_linkage_subunits,2,3,max_size)
real*8,  intent(inout) :: linkage_block_fitness(num_linkage_subunits,2,max_size)

dmutn = num_linkage_subunits*lb_modulo + 1
nmutn = num_linkage_subunits*lb_modulo + 1
fmutn = num_linkage_subunits*lb_modulo + 1
dmutn(1,:,:)  = 0
nmutn(1,:,:)  = 0
fmutn(1,:,:)  = 0
lb_mutn_count = 0
linkage_block_fitness = 1.d0

end subroutine test_init
