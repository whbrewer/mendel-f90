module genome
! This module represents genomic data by using an array of pointers 
! along with a user derived data type.

type genotype
  integer, pointer :: dm(:,:), nm(:,:), fm(:,:), lbmc(:,:,:)
  real*8, pointer :: lbf(:,:)
end type
type(genotype), allocatable, dimension(:) :: gp ! gene pool

contains

subroutine init_genome(n,dmutn,nmutn,fmutn,linkage_block_fitness,lb_mutn_count)
  use inputs
  integer, intent(in) :: n
  integer, target, intent(inout) :: dmutn(max_del_mutn_per_indiv/2,2,*)
  integer, target, intent(inout) :: nmutn(max_neu_mutn_per_indiv/2,2,*)
  integer, target, intent(inout) :: fmutn(max_fav_mutn_per_indiv/2,2,*)
  integer, target, intent(inout) :: lb_mutn_count(num_linkage_subunits,2,3,*)
  real*8,  target, intent(inout) :: linkage_block_fitness(num_linkage_subunits,2,*)

  integer :: i

  ! Setup pointers to represent gene pool
  do i = 1, n
    gp(i)%dm => dmutn(:,:,i)
    gp(i)%nm => nmutn(:,:,i)
    gp(i)%fm => fmutn(:,:,i)
    gp(i)%lbf  => linkage_block_fitness(:,:,i)
    gp(i)%lbmc => lb_mutn_count(:,:,:,i)
  end do
end subroutine init_genome

subroutine print_genotype(id,num)
   integer :: id ! the number of the indidual
   integer, optional :: num  ! the number of data to print

   if(present(num)) then
      n = num
   else
      n = 5
   end if

   write(6,*) '===== printing genome data ====='
   write(6,*) 'del:',gp(id)%dm(1:n,1)
   write(6,*) 'fav:',gp(id)%fm(1:n,1)
   write(6,*) 'neu:',gp(id)%nm(1:n,1)
   write(6,*) 'lbf:',gp(id)%lbf(1:n,1)
   write(6,*) 'lbmc:',gp(id)%lbmc(1:n,1,1)
end subroutine print_genotype

end module genome
