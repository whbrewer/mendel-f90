! Mendel unit testing framework

! Future tests to implement:
! make sure input file is read correctly
! if I set mutation rate = 0, fitness should be zero.

program test_main
use random_pkg
use selection_module
use inputs
use unittest
!include 'common.h'
implicit none
!include 'mpif.h'
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

!======= REUSABLE CODE =============================================
allocate(         dmutn(max_del_mutn_per_indiv/2,2,max_size),     &
                  nmutn(max_neu_mutn_per_indiv/2,2,max_size),     &
                  fmutn(max_fav_mutn_per_indiv/2,2,max_size),     &
                 lb_mutn_count(num_linkage_subunits,2,3,max_size),&
         linkage_block_fitness(num_linkage_subunits,2,max_size))

call test_init(dmutn,nmutn,fmutn,lb_mutn_count,linkage_block_fitness,max_size)
!======= REUSABLE CODE =============================================

!****** START UNIT TESTING *************
call test_restart(dmutn,nmutn,fmutn,lb_mutn_count,linkage_block_fitness,max_size)
!call test_selection()
!call test_competition()
!call test_migration()
!call test_mating()
!call test_offspring()
!call test_diagnostics_history_plot()
!****** END UNIT TESTING ***************

end program test_main

subroutine test_restart(dmutn,nmutn,fmutn,lb_mutn_count,linkage_block_fitness, &
                        max_size)
use unittest
implicit none
integer, intent(in)    :: max_size
integer, intent(inout) :: dmutn(imax,2,max_size)
integer, intent(inout) :: nmutn(imax,2,max_size)
integer, intent(inout) :: fmutn(imax,2,max_size)
integer, intent(inout) :: lb_mutn_count(nls,2,3,max_size)
real*8,  intent(inout) :: linkage_block_fitness(nls,2,max_size)
integer :: shutdown_gen, gen_0
real,    allocatable, dimension(:) :: initial_allele_effects
character*3 myid_str
allocate (initial_allele_effects(nls))
write(myid_str,'(i3.3)') 0 !myid+1

!call write_output_dump(dmutn,nmutn,fmutn,lb_mutn_count, &
!                       linkage_block_fitness,initial_allele_effects, &
!                       shutdown_gen,myid_str)

!call read_restart_dump(dmutn,nmutn,fmutn,lb_mutn_count, &
!                       linkage_block_fitness,     &
!                       initial_allele_effects,    &
!                       gen_0,max_size,myid_str)

! Note: need to put much more tests here, but mainly 
! needs to test if the routine is going to crash or not.
! Need to test parallel restarts
call assert(shutdown_gen == gen_0)

end subroutine test_restart

