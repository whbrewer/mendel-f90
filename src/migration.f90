subroutine migration(dmutn,nmutn,fmutn,linkage_block_fitness, &
           lb_mutn_count,gen,ierr,msg_num)
! This subroutine passes individuals from one tribe to another
! but only works in such a way: when one tribe sends N individuals
! it must also receive N individuals.  This subroutine supports
! three different migration models:
!  1. ring pass - all tribes are arranged as a ring, and each
!                 tribe passes to the tribe on its right
!  2. stepping stone - tribes also arranged in a ring, but
!                 each tribe both sends and receives from its
!                 left and right neighbor
!  3. island model - every tribe sends and receives with every
!                 other tribe.
use random_pkg
use inputs
use mpi
use mpi_helpers
include 'common.h'

integer dmutn(max_del_mutn_per_indiv/2,2,*)
integer nmutn(max_neu_mutn_per_indiv/2,2,*)
integer fmutn(max_fav_mutn_per_indiv/2,2,*)
real*8  linkage_block_fitness(num_linkage_subunits,2,*)
integer lb_mutn_count(num_linkage_subunits,2,3,*)
integer gen

integer dbuffs((max_del_mutn_per_indiv+4)*num_indiv_exchanged)
integer dbuffr((max_del_mutn_per_indiv+4)*num_indiv_exchanged)
integer nbuffs((max_neu_mutn_per_indiv+4)*num_indiv_exchanged)
integer nbuffr((max_neu_mutn_per_indiv+4)*num_indiv_exchanged)
integer fbuffs((max_fav_mutn_per_indiv+4)*num_indiv_exchanged)
integer fbuffr((max_fav_mutn_per_indiv+4)*num_indiv_exchanged)

integer i, j, k, l, m, n, ii, z, num_receiving_tribes, nie, id
integer dest, src, status(MPI_Status_size,8), requests(8)
integer max_num_dmutn, max_num_dmutn_recvd
integer max_num_nmutn, max_num_nmutn_recvd
integer max_num_fmutn, max_num_fmutn_recvd
real*8  lbuffs(2*num_linkage_subunits*num_indiv_exchanged)
real*8  lbuffr(2*num_linkage_subunits*num_indiv_exchanged)
real*8  t0, t1, time
integer cbuff1s(2*num_linkage_subunits*num_indiv_exchanged)
integer cbuff1r(2*num_linkage_subunits*num_indiv_exchanged)
integer cbuff2s(2*num_linkage_subunits*num_indiv_exchanged)
integer cbuff2r(2*num_linkage_subunits*num_indiv_exchanged)
integer cbuff3s(2*num_linkage_subunits*num_indiv_exchanged)
integer cbuff3r(2*num_linkage_subunits*num_indiv_exchanged)
integer msg_size_dbuff,msg_size_nbuff,msg_size_fbuff
integer msg_size_lbuff
integer msg_size_rdbuff,msg_size_rfbuff,msg_size_rnbuff
integer swap_list(num_indiv_exchanged*num_tribes)
integer num_linkage_subunits_passed
logical debug, timing

debug = .false.
timing = .false.

!  swap individuals between processors

if (migration_model == 1) then ! ring pass
   num_receiving_tribes = 1
   nie = num_indiv_exchanged
else if (migration_model == 2) then ! stepping stone
   num_receiving_tribes = 2
   nie = num_indiv_exchanged*2
else if (migration_model == 3) then ! island model
   num_receiving_tribes = num_tribes - 1
   nie = num_indiv_exchanged*(num_tribes -1)
else
   write(*,*) 'migration_model ',migration_model, 'not supported'
   stop
end if

!  randomly select individuals to be swapped

call randomly_select_individuals(nie,swap_list)

if(tracking_threshold < 1.) then

!  Find max number of mutations in dmutn, nmutn, and fmutn

   max_num_dmutn = 0
   max_num_nmutn = 0
   max_num_fmutn = 0

   do i=1, nie
      max_num_dmutn = max(max_num_dmutn,dmutn(1,1,swap_list(i)), &
           dmutn(1,2,swap_list(i)))
      max_num_nmutn = max(max_num_nmutn,nmutn(1,1,swap_list(i)), &
           nmutn(1,2,swap_list(i)))
      max_num_fmutn = max(max_num_fmutn,fmutn(1,1,swap_list(i)), &
           fmutn(1,2,swap_list(i)))
   end do

   if(max_num_dmutn .gt. max_del_mutn_per_indiv/2) then
      write(0,*) 'ERROR in migration subroutine'
      write(0,*) '      max_num_dmutn exceeds limit'
      call mpi_myabort()
   end if

   if(max_num_nmutn .gt. max_neu_mutn_per_indiv/2) then
      write(0,*) 'ERROR in migration subroutine'
      write(0,*) '      max_num_nmutn exceeds limit'
      call mpi_myabort()
   end if

   if(max_num_fmutn .gt. max_fav_mutn_per_indiv/2) then
      write(0,*) 'ERROR in migration subroutine'
      write(0,*) '      max_num_fmutn exceeds limit'
      call mpi_myabort()
   end if

   max_num_dmutn = max_num_dmutn + 2
   max_num_nmutn = max_num_nmutn + 2
   max_num_fmutn = max_num_fmutn + 2

   if (debug) then
      if(mod(gen,10)==0) then
         write(6,'(2x,"migrating",i3, "indivuals with max_dmutn= ",i5, &
                   2x," and max_fav_mutn = ", i3)') &
                   num_indiv_exchanged, 2*(max_num_dmutn), 2*(max_num_fmutn)
      end if
   end if

   if(max_num_dmutn > max_del_mutn_per_indiv/2) then
      write(6,*) 'ERROR: need to increase max_del_mutn_per_indiv'
      call mpi_myabort()
   end if

   if(max_num_nmutn > max_neu_mutn_per_indiv/2) then
      write(6,*) 'ERROR: need to increase max_neu_mutn_per_indiv'
      call mpi_myabort()
   end if

   if(max_num_fmutn > max_fav_mutn_per_indiv/2) then
      write(6,*) 'ERROR: need to increase max_fav_mutn_per_indiv'
      call mpi_myabort()
   end if

end if

!  populate buffer arrays for message passing

do m = 1, num_receiving_tribes ! loop over number of receiving tribes
! Note: there is a bug in the Intel Fortran compiler that sets
!       dbuff(1)=0 so start with k = 2.
   k = 2
   j = 2
   n = 2
   l = 1

   do i = 1, num_indiv_exchanged ! loop over number_individuals_exchanged

      id = (m-1)*num_indiv_exchanged + i

      if(tracking_threshold < 1.) then

         do ii = 1, 2 ! one for each haplotype
            if (k+max_num_dmutn-1 > size(dbuffs)) then
               write(*,*) 'ERROR: dbuffs overflow', k, max_num_dmutn, size(dbuffs)
               stop
            end if
            dbuffs(k:k+max_num_dmutn-1) = dmutn(1:max_num_dmutn,ii,swap_list(id))
            k = k + max_num_dmutn

            if (j+max_num_nmutn-1 > size(nbuffs)) then
               write(*,*) 'ERROR: nbuffs overflow', j, max_num_nmutn, size(nbuffs)
               stop
            end if
            nbuffs(j:j+max_num_nmutn-1) = nmutn(1:max_num_nmutn,ii,swap_list(id))
            j = j + max_num_nmutn

            if (n+max_num_fmutn-1 > size(fbuffs)) then
               write(*,*) 'ERROR: fbuffs overflow', n, max_num_fmutn, size(fbuffs)
               stop
            end if
            fbuffs(n:n+max_num_fmutn-1) = fmutn(1:max_num_fmutn,ii,swap_list(id))
            n = n + max_num_fmutn
         end do

      end if

      do ii = 1, 2
         if (l+num_linkage_subunits-1 > size(lbuffs)) then
            write(*,*) 'ERROR: lbuffs overflow', l, num_linkage_subunits, size(lbuffs)
            stop
         end if
         lbuffs(l:l+num_linkage_subunits-1) = linkage_block_fitness(:,ii,swap_list(id))
         cbuff1s(l:l+num_linkage_subunits-1) = lb_mutn_count(:,ii,1,swap_list(id))
         cbuff2s(l:l+num_linkage_subunits-1) = lb_mutn_count(:,ii,2,swap_list(id))
         cbuff3s(l:l+num_linkage_subunits-1) = lb_mutn_count(:,ii,3,swap_list(id))
         l = l + num_linkage_subunits
      end do

   end do

   msg_size_dbuff = k
   msg_size_nbuff = j
   msg_size_fbuff = n
   msg_size_lbuff = l - 1

!  Parallel debugging for two processors: writes the data sent to
!  fort.1 and the data received to fort.2.
   if (debug.and.myid==0) then
      write(1,*) '------ GENERATION:',gen,'----------------------'
      write(1,*) 'dmutn sent:',gen
      write(1,*) ((dmutn(ii,1,swap_list(1))),ii=1,max_num_dmutn)
      write(1,*) 'dbuff sent:',gen
      write(1,*) (dbuffs(ii),ii=1,max_num_dmutn+1)
      write(1,*) 'fmutn sent:',gen
      write(1,*) ((fmutn(ii,1,swap_list(1))),ii=1,max_num_fmutn)
      write(1,*) 'fbuff sent:',gen
      write(1,*) (fbuffs(ii),ii=1,max_num_fmutn+1)
      write(1,*) 'linkage_block_fitness sent - hap1:'
      write(1,*)((linkage_block_fitness(ii,1,swap_list(1))),ii=1,9)
      write(1,*) 'linkage_block_fitness sent - hap2:'
      write(1,*)((linkage_block_fitness(ii,2,swap_list(1))),ii=1,9)
      write(1,*) 'lb_mutn_count sent - hap1:'
      write(1,*) ((lb_mutn_count(ii,1,1,swap_list(1))),ii=1,9)
      write(1,*) 'lb_mutn_count - hap2:'
      write(1,*) ((lb_mutn_count(ii,2,1,swap_list(1))),ii=1,9)
   end if

   t0 = MPI_Wtime()

!  compute destination tribe based on migration model

   if(migration_model == 1 .or. migration_model == 3) then
      dest = mod(myid + m, num_tribes)
      src = MPI_ANY_SOURCE
   else if (migration_model == 2) then ! stepping stone
      if (mod(m,2)==1) then
         dest = mod(myid + m, num_tribes)
      else
         dest = mod(myid + num_tribes - 1, num_tribes)
      end if
      src = MPI_ANY_SOURCE
   else
      write(*,*) 'ERROR: migration_model not supported'
      stop
   end if

   if(tracking_threshold < 1.) then

!     send the maximum number of deleterious mutations
      call MPI_Isend(max_num_dmutn,1,MPI_INTEGER,dest,msg_num,MYCOMM,requests(1),ierr)
      call MPI_IRecv(max_num_dmutn_recvd,1,MPI_INTEGER,src,msg_num,MYCOMM,requests(2),ierr)
      msg_num = msg_num + 1

!     send the maximum number of neutral mutations
      call MPI_Isend(max_num_nmutn,1,MPI_INTEGER,dest,msg_num,MYCOMM,requests(3),ierr)
      call MPI_IRecv(max_num_nmutn_recvd,1,MPI_INTEGER,src,msg_num,MYCOMM,requests(4),ierr)
      msg_num = msg_num + 1

!     send the maximum number of favorable mutations
      call MPI_Isend(max_num_fmutn,1,MPI_INTEGER,dest,msg_num,MYCOMM,requests(5),ierr)
      call MPI_IRecv(max_num_fmutn_recvd,1,MPI_INTEGER,src,msg_num,MYCOMM,requests(6),ierr)
      msg_num = msg_num + 1

!     complete the nonblocking sends+receives before proceeding
      call MPI_Waitall(6,requests,status,ierr)

!     communicate buffer of deleterious mutations
      msg_size_rdbuff = 2*(max_num_dmutn_recvd+1)*num_indiv_exchanged
      call MPI_Isend(dbuffs,msg_size_dbuff,MPI_INTEGER, &
           dest,msg_num,MYCOMM,requests(1),ierr)
      call MPI_IRecv(dbuffr,msg_size_rdbuff,MPI_INTEGER, &
           src,msg_num,MYCOMM,requests(2),ierr)
      msg_num = msg_num + 1

!     communicate buffer of neutral mutations
      msg_size_rnbuff = 2*(max_num_nmutn_recvd+1)*num_indiv_exchanged
      call MPI_Isend(nbuffs,msg_size_nbuff,MPI_INTEGER, &
           dest,msg_num,MYCOMM,requests(3),ierr)
      call MPI_IRecv(nbuffr,msg_size_rnbuff,MPI_INTEGER, &
           src,msg_num,MYCOMM,requests(4),ierr)
      msg_num = msg_num + 1

!     communicate buffer of favorable mutations
      msg_size_rfbuff = 2*(max_num_fmutn_recvd+1)*num_indiv_exchanged
      call MPI_Isend(fbuffs,msg_size_fbuff,MPI_INTEGER, &
           dest,msg_num,MYCOMM,requests(3),ierr)
      call MPI_IRecv(fbuffr,msg_size_rfbuff,MPI_INTEGER, &
           src,msg_num,MYCOMM,requests(4),ierr)
      msg_num = msg_num + 1

      call MPI_Waitall(4,requests,status,ierr)

   end if

!  communicate buffer of linkage_block_fitness
   call MPI_Isend(lbuffs,msg_size_lbuff,MPI_DOUBLE_PRECISION, &
        dest,msg_num,MYCOMM,requests(1),ierr)
   call MPI_IRecv(lbuffr,msg_size_lbuff,MPI_DOUBLE_PRECISION, &
        src,msg_num,MYCOMM,requests(2),ierr)
   msg_num = msg_num + 1

!  communicate buffer of lb_mutn_count
   call MPI_Isend(cbuff1s,msg_size_lbuff,MPI_INTEGER, &
        dest,msg_num,MYCOMM,requests(3),ierr)
   call MPI_IRecv(cbuff1r,msg_size_lbuff,MPI_INTEGER, &
        src,msg_num,MYCOMM,requests(4),ierr)
   msg_num = msg_num + 1

!  communicate buffer of lb_mutn_count
   call MPI_Isend(cbuff2s,msg_size_lbuff,MPI_INTEGER, &
        dest,msg_num,MYCOMM,requests(5),ierr)
   call MPI_IRecv(cbuff2r,msg_size_lbuff,MPI_INTEGER, &
        src,msg_num,MYCOMM,requests(6),ierr)
   msg_num = msg_num + 1

!  communicate buffer of lb_mutn_count
   call MPI_Isend(cbuff3s,msg_size_lbuff,MPI_INTEGER, &
        dest,msg_num,MYCOMM,requests(7),ierr)
   call MPI_IRecv(cbuff3r,msg_size_lbuff,MPI_INTEGER, &
        src,msg_num,MYCOMM,requests(8),ierr)
   msg_num = msg_num + 1

   call MPI_Waitall(8,requests,status,ierr)

   t1 = MPI_Wtime()

   time = t1 - t0

   if(mod(gen,10)==0.and.timing) write(*,*) "  communication time = ",time,"seconds"

!  Unpack transferred buffers

   k = 2 ! see note at first k = 2 statement
   j = 2
   n = 2
   l = 1

   if (debug) then
      if(max_num_dmutn == max_num_dmutn_recvd) then
         write(*,*)'  WARNING: max_num_dmutn same for 2 tribes', max_num_dmutn
      end if
      if(max_num_nmutn == max_num_nmutn_recvd) then
         write(*,*)'  WARNING: max_num_nmutn same for 2 tribes', max_num_nmutn
      end if
      if(max_num_fmutn == max_num_fmutn_recvd) then
         write(*,*)'  WARNING: max_num_fmutn same for 2 tribes', max_num_fmutn
      end if
   end if

   do i = 1, num_indiv_exchanged

      id = (m-1)*num_indiv_exchanged + i

      if(tracking_threshold < 1.) then

         dmutn(:,:,swap_list(id)) = num_linkage_subunits*lb_modulo+1
         nmutn(:,:,swap_list(id)) = num_linkage_subunits*lb_modulo+1
         fmutn(:,:,swap_list(id)) = num_linkage_subunits*lb_modulo+1

         do ii = 1, 2
            dmutn(1:max_num_dmutn_recvd,ii,swap_list(id)) = &
                  dbuffr(k:k+max_num_dmutn_recvd-1)
            k = k + max_num_dmutn_recvd

            nmutn(1:max_num_nmutn_recvd,ii,swap_list(id)) = &
                  nbuffr(j:j+max_num_nmutn_recvd-1)
            j = j + max_num_nmutn_recvd

            fmutn(1:max_num_fmutn_recvd,ii,swap_list(id)) = &
                  fbuffr(n:n+max_num_fmutn_recvd-1)
            n = n + max_num_fmutn_recvd
         end do

      end if

      do ii = 1, 2
         linkage_block_fitness(:,ii,swap_list(id)) = lbuffr(l:l+num_linkage_subunits-1)
         lb_mutn_count(:,ii,1,swap_list(id)) = cbuff1r(l:l+num_linkage_subunits-1)
         lb_mutn_count(:,ii,2,swap_list(id)) = cbuff2r(l:l+num_linkage_subunits-1)
         lb_mutn_count(:,ii,3,swap_list(id)) = cbuff3r(l:l+num_linkage_subunits-1)
         l = l + num_linkage_subunits
      end do

   end do

   if (debug.and.myid==1) then
      write(2,*) '------ GENERATION:',gen,'----------------------'
      write(2,*) 'swap_list for proc:',myid
      write(2,*) ((swap_list(z)),z=1,num_indiv_exchanged)
      write(2,*) msg_num,'received',max_num_dmutn_recvd,max_num_fmutn_recvd
      write(2,*) 'dmutn received:',gen
      write(2,*) ((dmutn(ii,1,swap_list(1))),ii=1,max_num_dmutn_recvd)
      write(2,*) 'dbuff received:',gen
      write(2,*) (dbuffr(ii),ii=1,max_num_dmutn_recvd+1)
      write(2,*) 'fmutn received:',gen
      write(2,*) ((fmutn(ii,1,swap_list(1))),ii=1,max_num_fmutn_recvd)
      write(2,*) 'fbuff received:',gen
      write(2,*) (fbuffr(ii),ii=1,max_num_fmutn_recvd+1)
      write(2,*) 'linkage_block_fitness received - hap1:'
      write(2,*)((linkage_block_fitness(ii,1,swap_list(1))),ii=1,9)
      write(2,*) 'linkage_block_fitness received - hap2:'
      write(2,*)((linkage_block_fitness(ii,2,swap_list(1))),ii=1,9)
      write(2,*) 'lb_mutn_count received - hap1:'
      write(2,*) ((lb_mutn_count(ii,1,1,swap_list(1))),ii=1,9)
      write(2,*) 'lb_mutn_count received - hap2:'
      write(2,*) ((lb_mutn_count(ii,2,1,swap_list(1))),ii=1,9)
   end if

end do
end

! randomly select individuals to be swapped
subroutine randomly_select_individuals(nie,swap_list)
use random_pkg
use inputs
include 'common.h'
integer i, j, nie
logical available(current_pop_size)
integer swap_list(num_indiv_exchanged*num_tribes)
available = .true.

do i = 1, nie
   j = min(current_pop_size,1 + int(current_pop_size*randomnum(1)))
   do while (.not. available(j))
      j = mod(j,current_pop_size) + 1
   end do
   available(j) = .false.
   swap_list(i) = j
end do
end

subroutine migrate_individual(other,sid,did,dmutn,fmutn,nmutn,lb_mutn_count, &
           linkage_block_fitness,sender)
! This subroutine just passes a single individual to another tribes
! in a unidirectional sense (i.e. one tribe sends, the other tribe
! receives).  sid is the number of the individual that is to be sent
! did is the destination number where to unpack that individual.
! NOTE: One current limitation of this routine is that the individual
! is not removed from the originating tribe.  The way tribal fissioning
! is setup, we deal with this by copying half the individuals, and then
! reduce the population size by half.  However, this should be generalized
! in the future, so that any individual can be migrated, and also deleted
! from the originating tribe.
! Credits: The skeleton of this code came from:
! http://stackoverflow.com/questions/13211990
use inputs
use mpi
include 'common.h'

integer, parameter :: n = 5000, m = 1000, NV = 5
integer, intent(in) :: sid, did  ! source id, destination id
integer, intent(in) :: other
integer, intent(inout) :: dmutn(max_del_mutn_per_indiv/2,2,*)
integer, intent(inout) :: fmutn(max_fav_mutn_per_indiv/2,2,*)
integer, intent(inout) :: nmutn(max_neu_mutn_per_indiv/2,2,*)
integer, intent(inout) :: lb_mutn_count(num_linkage_subunits,2,3,*)
real*8,  intent(inout) :: linkage_block_fitness(num_linkage_subunits,2,*)
logical, intent(in) :: sender

type individual
  integer :: del(n,2)
  integer :: fav(n,2)
  integer :: neu(n,2)
  integer :: lmc(m,2,3)
  real*8  :: lbf(m,2)
end type

type (individual) :: bob
integer :: datatype, oldtypes(NV), blockcounts(NV)
integer(kind=MPI_ADDRESS_KIND) :: offsets(NV)
integer :: i, myint
integer :: status(MPI_Status_size)

if (num_linkage_subunits > m) then
  write(6,*) "ERROR in migrate_individual: need to increase m array bounds"
  call exit(1)
endif

if (dmutn(1,1,sid) > n .or. fmutn(1,1,sid) > n .or. nmutn(1,1,sid) > n) then
  write(6,*) "ERROR in migrate_individual: need to increase n array bounds"
  call exit(1)
endif

if (sender) then
  bob%del(:n,:) = dmutn(:n,:,sid)
  bob%fav(:n,:) = fmutn(:n,:,sid)
  bob%neu(:n,:) = nmutn(:n,:,sid)
  bob%lmc(:m,:,:) = lb_mutn_count(:m,:,:,sid)
  bob%lbf(:m,:) = linkage_block_fitness(:m,:,sid)
else
  bob%del = 0
  bob%fav = 0
  bob%neu = 0
  bob%lmc = 0
  bob%lbf = 0
endif

call mpi_get_address(bob%del, offsets(1), ierr)
call mpi_get_address(bob%fav, offsets(2), ierr)
call mpi_get_address(bob%neu, offsets(3), ierr)
call mpi_get_address(bob%lmc, offsets(4), ierr)
call mpi_get_address(bob%lbf, offsets(5), ierr)

do i = 2, size(offsets)
  offsets(i) = offsets(i) - offsets(1)
end do
offsets(1) = 0

oldtypes = (/mpi_integer, mpi_integer, mpi_integer, mpi_integer, mpi_real8/)
blockcounts = (/2*n, 2*n, 2*n, 6*m, 2*m/)

! Note: mpi_type_struct is deprecated, should use mpi_type_create_struct
! instead. However, mpi_type_create_struct crashes for some reason
!call mpi_type_create_struct(NV, blockcounts, offsets, oldtypes, datatype, ierr)
call mpi_type_struct(NV, blockcounts, offsets, oldtypes, datatype, ierr)
call mpi_type_commit(datatype, ierr)

if (sender) then
  call mpi_send(bob, 1, datatype, other, 0, MYCOMM, ierr)
else ! receiver
  call mpi_recv(bob, 1, datatype, other, 0, MYCOMM, status, ierr)
end if

! receiver
if (.not. sender) then
  dmutn(:n,:,did) = bob%del(:n,:)
  fmutn(:n,:,did) = bob%fav(:n,:)
  nmutn(:n,:,did) = bob%neu(:n,:)
  lb_mutn_count(:m,:,:,did) = bob%lmc(:m,:,:)
  linkage_block_fitness(:m,:,did) = bob%lbf(:m,:)
endif

end
