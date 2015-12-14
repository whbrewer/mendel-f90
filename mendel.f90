program mendel

use init
use inputs
use genome
use profiler
use polygenic
use random_pkg
use selection_module
include 'common.h'
! START_MPI
include 'mpif.h'
! END_MPI

! Data structures to maintain genetic information
real*8,  allocatable, target, dimension(:,:,:)   :: linkage_block_fitness
integer, allocatable, target, dimension(:,:,:)   :: dmutn, nmutn, fmutn
integer, allocatable, target, dimension(:,:,:,:) :: lb_mutn_count

real,    allocatable, dimension(:) :: initial_allele_effects
real*8,  allocatable, dimension(:) :: fitness, pheno_fitness
real*8,  allocatable, dimension(:) :: work_fitness, sorted_score
logical, allocatable, dimension(:) :: available
integer, allocatable, dimension(:) :: global_run_status
integer, allocatable, dimension(:) :: pop_size_array

integer*4 :: chmod
integer :: i, j, k, lb, m, nmax, gen, run_status
integer :: npath, max_size, this_size, gen_0, shutdown_gen
integer :: other_run_status, bottleneck_modulo, ibuff
! START_MPI
integer :: status2(MPI_Status_size,2), requests(2)
! END_MPI
integer :: total_offspring, offspring_count, empty, red, green, blue
integer :: ica_count(3), cumulative_offspring, mutn_indx
integer :: pop_size_allocation, current_global_pop_size
integer :: global_pop, mean_pop, delta_pop, delta, nie, mm, id
integer :: pop_size_winner, pop_size_loser, num_migrate
integer :: id_winner, id_loser
integer :: num_dmutns, num_fmutns, encode_mutn, string(40)
integer :: OLDGROUP,NEWGROUP,ranks(1),num_tribes_at_start

real*8 accum(50), reproductive_advantage_factor
real selection_coefficient, aoki, migration_rate, x
real tribal_score, random_effects, genetic_effects, social_effects
real fraction_elimination, fraction_selected_away
real random, num_offspring, fav_mutn_per_gen, d
real real_pop_size
real tin_migration, tout_migration, tin_gen, tout_gen, tin_run
real tin_offspring, tout_offspring, tgen, par_tgen
real tin_diagnostics, tout_diagnostics, tin_selection, tout_selection
real total_time_offspring, time_offspring, time_selection
real par_time_offspring, par_time_selection, tsub
logical found, print_flag, am_parallel, file_exists, winner, create_fav_mutn
character*3 myid_str
character(len=128) :: arg, filename

call second(tin_run)

! Get command-line arguments
do i = 1, iargc()
   call getarg(i, arg)
   if (arg == "-h") then
      print *, "usage: mendel {-c | -f} [filename]"
      print *, "-c creates a new mendel.in file"
      print *, "-f uses the specified filename"
      print *, "   e.g. mendel -f /home/bob/mendel.in"
      stop
   elseif (arg == "-d") then
      call set_default_parameters
      print *, "Using default parameters"
      exit
   elseif (arg == "-f") then
      call getarg(i+1, filename)
      print *, "filename is: ", filename
      exit
   elseif (arg == "-c") then
      print *, "create input file"
      open (5, file='./mendel.in.new')
      call set_default_parameters
      call write_parameters(5)
      close(5)
      print *, "Parameter file written: mendel.in.new" 
      stop
   else ! web interface
      filename = './mendel.in'
      ! print *, "ERROR: argument not supported. Try mendel -h", arg
      ! stop
   end if
end do

if ( arg /= "-d" ) then
!  Default filename
   if ( iargc() == 0 )  then
      filename = './mendel.in'
   end if
   print *, "Using filename: ", filename

!  Open file containing run parameters
   inquire(file=filename, exist=file_exists)
   if ( file_exists ) then
      open (5, file=filename, status='old')
   else
      print *, 'ERROR: could not find ',filename,len(filename)
      stop
   end if

!  Read input parameters from input file mendel.in.
   call read_parameters(5)
   close(5)
end if

! Perform certain initializations including opening files.

call initialize(myid_str)
run_status = 0

if(is_parallel) then 
! Since we may turn off is_parallel in case a tribe dies,
! remember the original state.
  am_parallel = .true.
  num_tribes_at_start = num_tribes
! compute global population size
  !START_MPI
  call mpi_isum(pop_size,global_pop_size,1)
  call mpi_mybcasti(global_pop_size,1)
  !END_MPI
  allocate( global_run_status(num_tribes) )
  allocate( pop_size_array(num_tribes) )
else
  am_parallel = .false.
endif
 
if(is_parallel .and. tribal_competition) then
   global_run_status = 0
   run_status = 0
   other_run_status = 0
   pop_size_array = pop_size
   pop_size_allocation = global_pop_size
   write(*,*) 'Allocating tribe ', myid, ' with max pop_size of:',&
                pop_size_allocation
elseif (pop_growth_model > 0) then
   pop_size_allocation = carrying_capacity*100.
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

! Allocate memory for large arrays.
allocate(         dmutn(max_del_mutn_per_indiv/2,2,max_size),     &
                  nmutn(max_neu_mutn_per_indiv/2,2,max_size),     &
                  fmutn(max_fav_mutn_per_indiv/2,2,max_size),     &
                 lb_mutn_count(num_linkage_subunits,2,3,max_size),&
         linkage_block_fitness(num_linkage_subunits,2,max_size),  &
        initial_allele_effects(num_linkage_subunits),             &
         pheno_fitness(max_size),      fitness(max_size),         &
          work_fitness(max_size), sorted_score(max_size),         &
                        available(pop_size_allocation) )
allocate( gp(max_size) )
if(polygenic_beneficials) then
   allocate( pmutn(max_polys) )
   poly_gen_first_instance = -1
   num_polys_cumulative = 0
endif

!call init_genome(max_size,dmutn,nmutn,fmutn,linkage_block_fitness,lb_mutn_count)

x = 1E6 ! Megabyte
print *, '-------------------------------------'
print *, 'MEMORY USAGE INFORMATION (MBYTES):'
print *, 'genotype pointer:', sizeof(gp)/x
print *, 'dmutn size:'   , sizeof(dmutn)/x
print *, 'fmutn size:'   , sizeof(fmutn)/x
print *, 'nmutn size:'   , sizeof(nmutn)/x
print *, 'lb_mutn_count size:',sizeof(lb_mutn_count)/x
print *, 'linkage_block_fitness size:',sizeof(linkage_block_fitness)/x
print *, '-------------------------------------'
print *
print *, 'Initializing data arrays... please wait....'

! If this is a restart case, read the restart dump file and
! set the current dump number to the restart dump number.
! Otherwise, set it to zero.  The variable gen_0 is the initial
! generation number, retrieved from the restart dump in a restart
! case and zero otherwise.

if(restart_case) then
   call read_restart_dump(dmutn,nmutn,fmutn,lb_mutn_count, &
                          linkage_block_fitness,     &
                          initial_allele_effects,    &
                          gen_0,max_size,myid_str)
   dump_number = restart_dump_number
else
   gen_0 = 0
   dump_number = 0
end if

! If the bottleneck flag, bottleneck_yes, is false, set the value
! of bottleneck_generation beyond the generation range for this run.

if(.not.bottleneck_yes) &
   bottleneck_generation = 1 + gen_0 + num_generations

! Initialize the population size to be equal to the parameter
! pop_size unless the parameter bottleneck_generation has the
! value zero.  In the latter case, initialize the population size 
! to bottleneck_pop_size.

if(abs(bottleneck_generation) > 0) then
   current_pop_size = pop_size
else
   current_pop_size = bottleneck_pop_size
end if

! Setup cyclic bottlenecking

if(bottleneck_yes.and.bottleneck_generation < 0) then
   if(num_bottleneck_generations >= abs(bottleneck_generation))   &
   then
      write(*,*) 'ERROR: num_bottleneck_generations ',            &
                 '>= cyclic bottleneck_generations'
      stop
   end if  
   bottleneck_modulo = abs(bottleneck_generation)
   bottleneck_generation = abs(bottleneck_generation)
   cyclic_bottlenecking = .true.
end if

! If not a restart case, initialize entire population to have no 
! initial mutations.   

! Initialize the linkage block fitness such that all individuals
! in the population have identical haplotypes.  If initial 
! contrasting alleles are to be included, generate them here. 

if(.not.restart_case) then
   dmutn = num_linkage_subunits*lb_modulo + 1
   nmutn = num_linkage_subunits*lb_modulo + 1
   fmutn = num_linkage_subunits*lb_modulo + 1
   dmutn(1,:,:)  = 0
   fmutn(1,:,:)  = 0

   ! With polygenic beneficials each linkage block represents a single
   ! nucleotide.  Correspondingly, we do not accumulate mutations as 
   ! in the traditional sense, but rather we maintain a single mutation
   ! for each linkage block.  Therefore, number of mutations in each
   ! haplotype is constant, equal to the number of linkage subunits.
   ! We also use the array fmutn to store whether the string in nmutn
   ! matches the target or not and, if so, identify that instance with
   ! a unique positive identification number in fmutn(2,:,:).

   if (polygenic_beneficials) then
      nmutn(1,:,:) = num_linkage_subunits
      fmutn(1,:,:) = 1
      fmutn(2,:,:) = 0
      write(*,'(/A,$)') 'POLYGENIC STRING INITIALIZATION: '
      do lb=1,num_linkage_subunits
         nmutn(1+lb,:,:) = nucl_to_int(polygenic_init(lb:lb))
         write(*,'(a,$)') int_to_nucl(nucl_to_int(polygenic_init(lb:lb)))
      enddo
      write(*,*)
!     if(recombination_model==full_sexual) recombination_model = suppressed
   else
      nmutn(1,:,:) = 0
      lb_mutn_count = 0 
   endif
   linkage_block_fitness = 1.d0
   if(num_contrasting_alleles > 0)                                &
      call gen_initial_contrasting_alleles(dmutn, fmutn,          &
         linkage_block_fitness, initial_allele_effects, max_size)
end if

! Read in a file containing a specific set of mutations.

if(upload_mutations) then
   call read_mutn_file(dmutn,nmutn,fmutn,lb_mutn_count,           &
                       linkage_block_fitness,max_size)
end if

! Generate num_initial_fav_mutn random initial favorable mutations.

do k=1,num_initial_fav_mutn
   call favorable_mutn(fmutn,lb_mutn_count,linkage_block_fitness)
end do

post_sel_fitness  = 1.d0 
ica_count         = 0

call second(tout)
sec(1) = sec(1) + tout - tin

! Step population through num_generations generations.

do gen=gen_0+1,gen_0+num_generations

   !call print_genotype(1)
   !call print_genotype(1,10)

   msg_num = 1

   call second(tin_gen)

!  If the generation number lies within the bottleneck interval,
!  set the current population size to bottleneck_pop_size.

   if(cyclic_bottlenecking.and.                                   &
      (mod(gen,bottleneck_modulo)==0                              &
      .and.gen>gen_0+1+bottleneck_modulo)) then                   
      bottleneck_generation = bottleneck_generation +             &
                              bottleneck_modulo
   end if

   if(bottleneck_yes .and. gen >= bottleneck_generation .and.     &
      gen <  bottleneck_generation + num_bottleneck_generations)  &
   then
      current_pop_size = bottleneck_pop_size
      do i=6,9,3
         write(i,'(/"BOTTLENECK down to ",                        &
                   i6," individual(s) at generation = ",i6)')     &
                   bottleneck_pop_size, gen
      end do
   end if

   fertility_factor = 1. - fraction_random_death

!  For competing tribes compute tribal fertility factor.

   if(is_parallel .and. tribal_competition .and.                  &
      gen > gen_0+1) then

      reproductive_advantage_factor = tribal_fitness_factor-1.d0

      ! Compute the tribal fertility factor.  The number is used
      ! to control the size of the various tribes.
      fertility_factor = (1.d0 - fraction_random_death)           &
                         *(1.d0 + tc_scaling_factor               &
                                  *reproductive_advantage_factor)

      !selection_coefficient = max(1.d0 - fertility_factor,0.d0)
      selection_coefficient = 1.d0 - fertility_factor

      if(mod(gen,10)==0) then
         write(*,'("myid =",i2," fertility_factor =",f7.4," selection_coefficient =",f7.4)')      &
               myid, fertility_factor, selection_coefficient
      end if

   endif

!  Move individuals between tribes/processors.

!START_MPI
   if (is_parallel.and.mod(gen,migration_generations)==0.and.     &
       num_indiv_exchanged > 0) then
      call second(tin_migration)
      if(mod(gen,10)==0 .and. myid==0) then
         write(6,'(/"migrating ",i4," individual(s) every",       &
             i4," generation(s) between",i4," tribes")')          &
             num_indiv_exchanged, migration_generations,          &
             num_tribes
         if(tribal_competition) then
            write(*,*) 'competing pop sizes:', pop_size_array
         endif
             
      end if
      call migration(dmutn,nmutn,fmutn,linkage_block_fitness, &
           lb_mutn_count,gen,ierr,msg_num)
      call second(tout_migration)
      sec(4) = sec(4) + tout_migration - tin_migration
   end if
!END_MPI

   call second(tin_offspring)

   num_offspring    = 2.d0*reproductive_rate*fertility_factor

   if(mod(gen - gen_0, 10) == 1 .or. gen - gen_0 <= 10) then
      cumulative_offspring = 0
      new_mutn_count       = 0
   end if

!  Re-initialize random number generator using PID xor Time.
   if(reseed_rng) call init_random_seed() 

!  Randomly mate one half of the population with members
!  from the other half.

   time_offspring = 0

   call mating(dmutn,nmutn,fmutn,lb_mutn_count,linkage_block_fitness, &
        num_offspring,available,pop_size_allocation,nmax,offspring_count, &
        total_offspring,tsub,gen)

   ! create a beneficial mutation for each individual which received 
   ! the polygenic target during this generation

   if(polygenic_beneficials) then
      ! For waiting time experiments we map the mutations to nucleotides
      ! (e.g. the mutation number 732853 represents a mutation such as A->C)
      ! therefore, we need to maintain one mutation per linkage block
      ! in the neutral mutation array.

      ! Loop through all individuals in the population and count the number
      ! of matches with the target nucleotide string. 

      num_polys_this_gen = 0
      do id=1,current_pop_size
         do j=1,2 ! haplotype
            if(fmutn(2,j,id) > 0) then
               num_polys_this_gen = num_polys_this_gen + 1 
            endif 
         enddo
      enddo

      percent_pop_poly = real(num_polys_this_gen)/real(current_pop_size)/2.*100.
      if(recombination_model == clonal) percent_pop_poly = 2.*percent_pop_poly

      if (percent_pop_poly >= 99) then

          polygenic_fixed = .true.
          poly_stop_gen = gen
          print *
          print *, 'POLYGENICS: SHUTTING DOWN BECAUSE 99% OF POPULATION HAS ALLELE'

          call diagnostics_history_plot(dmutn, nmutn, fmutn, lb_mutn_count, &
               ica_count, gen, .true., current_global_pop_size)

          goto 20 ! shutdown
      endif

   endif

   time_offspring = time_offspring + tsub

!  Because the Poisson random number generator does not yield
!  the specified mean number of new mutations to sufficient
!  accuracy, to improve accuracy make an adjustment to the value
!  fed to the generator.

   cumulative_offspring = cumulative_offspring + offspring_count

   if((mod(gen - gen_0, 10) == 0 .or. gen - gen_0 < 10) .and. &
      mutn_rate >= 1.) then 
      d = real(new_mutn_count)/real(cumulative_offspring)
      poisson_mean = poisson_mean + 0.3*(mutn_rate - d)
   end if

   if(is_parallel .and. tribal_competition) then

!     Modify the tribal population size such that selection 
!     intensity depends only on the default fertility.

      real_pop_size = offspring_count/(reproductive_rate &
                      *(1. - fraction_random_death)) 
      if(fitness_dependent_fertility) real_pop_size =  &
         real_pop_size/sqrt(min(1.d0, post_sel_fitness))
      current_pop_size = int(real_pop_size)
      if(real_pop_size - current_pop_size > randomnum(1)) &
         current_pop_size = current_pop_size + 1
      migration_rate = num_indiv_exchanged/real(current_pop_size)
!     Following is the k value, what Aoki calls "group selection intensity"
!     ref: Aoki, Kenichi, "A condition for group selection to prevail over 
!     counteracting selection" by Kenichi Aoki, Evolution 36(4), 1982, 
!     pp. 832-842.
      aoki = 2*selection_coefficient*current_pop_size*migration_rate
      if(myid.eq.0) write(*,'(a,i5,x,a,f7.4,x,a,f7.4,x,a,i7)') & 
         'gen:', gen, 'group_selection_intensity: ', aoki,     &
         'migration_rate:', migration_rate,  &
         'deme_size:', current_pop_size

!     Modify the tribal population size to keep the global
!     population size nearly constant.

      !START_MPI
      call mpi_mybarrier()

      call mpi_isum(current_pop_size,global_pop,1)
      call mpi_mybcasti(global_pop,1)
      !END_MPI

      mean_pop  = global_pop_size/num_tribes
      delta_pop = global_pop_size-global_pop
      delta     = 2*delta_pop/num_tribes

      if(delta_pop > 0 .and. current_pop_size > mean_pop) then
         k = min(delta, int(delta*randomnum(1) + 0.5))
         current_pop_size = current_pop_size + k
      elseif(delta_pop<0 .and. current_pop_size<mean_pop) then
         k = max(delta, int(delta*randomnum(1) - 0.5))
         current_pop_size = current_pop_size + k
      end if

! START_MPI
      call mpi_gather(current_pop_size,1,mpi_integer, &
                        pop_size_array,1,mpi_integer,0,mycomm,ierr)
! END_MPI

   end if

!  If there is tribal competition and all tribes but one go
!  extinct, let the remaining tribe grow to the maximum size.

   if(tribal_competition .and. .not.is_parallel) then
      k = real(global_pop_size - current_pop_size)*0.1 + 0.9
      current_pop_size = current_pop_size + k
   end if

   current_pop_size = min(current_pop_size, offspring_count)
   !print *, '***',current_pop_size, offspring_count, pop_size

   call second(tout_offspring)
   sec(5) = sec(5) + tout_offspring - tin_offspring

!  Impose selection based on fitness to reduce the population 
!  size to a value not to exceed the parameter pop_size.
    
   call second(tin_selection)
   call selection(dmutn, nmutn, fmutn, lb_mutn_count,  &
        linkage_block_fitness, fitness, pheno_fitness, &
        work_fitness, sorted_score, initial_allele_effects, &
        max_size, total_offspring, gen, lb_modulo, current_pop_size)
   !call selection2(fitness, pheno_fitness, &
   !     work_fitness, sorted_score, initial_allele_effects, &
   !     max_size, total_offspring, gen, lb_modulo, current_pop_size)
   call second(tout_selection)
   time_selection = tout_selection - tin_selection
   sec(6) = sec(6) + time_selection

! START_MPI
   if(tribal_competition) then
      call compute_tribal_fitness(dmutn, fmutn, pop_size_array, &
           current_global_pop_size,gen)
   end if
! END_MPI

!  If the population size or the mean fitness has collapsed,
!  print message and shutdown all processors.

   if(post_sel_fitness < extinction_threshold) then

      ! If one tribe goes extinct, set trigger for shutdown
      if(is_parallel) then
         ! Must shutdown parallel before calling diagnostics_history_plot
         ! otherwise will hang waiting for communications
         if(num_tribes == 2) is_parallel = .false.
         run_status = -1
      else 
         do i=6,9,3
            write(i,'(/"** SHUTDOWN DUE TO EXTINCTION **"/)')
         end do
         goto 20
      end if

      call diagnostics_history_plot(dmutn, nmutn, fmutn,  &
           lb_mutn_count, ica_count, gen, .true.,  &
           current_global_pop_size)

      do i=6,9,3
         write(i,*)
         write(i,*)'*** SHUTDOWN TRIBE',myid+1, &
                   ' DUE TO EXTINCTION AT GEN:',gen
      end do

   end if

!  Since clonal cells such as bacteria can multiply from a single cell
!  we need to provide some means to allow a population size of one
!  to recover from a bottleneck. This is done by increasing it to
!  a value of 2 just before the mutational check, which will cause
!  it to rebound correctly. If it is left as 1, it will never rebound 
!  correctly, but stay at a fixed value of 1 for the remainder of the 
!  simulation.

   if(bottleneck_yes .and. bottleneck_pop_size == 1 .and. &
      recombination_model == clonal .and. current_pop_size == 1) then
      current_pop_size = 2
   end if

   !  Mutational meltdown scenario
   if((is_parallel .and. current_pop_size <  &
      extinction_threshold*pop_size) .or. &
      current_pop_size <= 1) then
      
      ! If one tribe melts down, set trigger for shutdown
      if(is_parallel) then
         do i=6,9,3
            write(i,*)
            write(i,*)'*** SHUTDOWN TRIBE',myid+1, &
                    'DUE TO MUTATIONAL MELTDOWN AT GEN:',gen
            write(i,*) myid+1,'Population size:', current_pop_size
         end do
         run_status = -(myid + 1)
      else
         do i=6,9,3
            write(i,'(/"** SHUTDOWN DUE TO MUTATIONAL MELTDOWN **")')
            write(i,*) 'Population size:', current_pop_size
         end do
         goto 20
      end if

   end if

   ! In case one tribe is set to run less generations than the other
   if(is_parallel .and. .not.homogenous_tribes .and. &
      gen==gen_0+num_generations) then
      do i=6,9,3 
         write(i,*) 'TRIBE',myid+1,'IS SHUTTING DOWN. GEN:',gen
      end do
      run_status = -(myid + 1)
   end if

   ! START_MPI
   ! For the limiting case of two tribes, we must turn off the parallel
   ! switch if one of the tribes goes extinct.  So, every generation 
   ! communicate the status of each tribe to the other.

   if(tribal_competition) then

      if(num_tribes == 2) then

         if(myid == 1) then
            call mpi_send_int(run_status,0,msg_num,ierr)
            msg_num = msg_num + 1
            call mpi_recv_int(other_run_status,0,msg_num,ierr)
         else 
            call mpi_recv_int(other_run_status,1,msg_num,ierr)
            msg_num = msg_num + 1
            call mpi_send_int(run_status,1,msg_num,ierr)
         end if
         msg_num = msg_num + 1
            
         if(tribal_fission) then

            ! Simple tribal fission - if one tribe goes extinct, split the
            ! surviving tribe into two and send half of its population
            ! to the tribe that went extinct  
            
            if (run_status < 0 .or. other_run_status < 0 ) then
               
               if(myid==0) then
                  write(6,*) myid,'competing pop sizes:', pop_size_array
                  pop_size_winner = maxval(pop_size_array)
                  pop_size_loser = minval(pop_size_array)
                  id_winner = maxloc(pop_size_array,1) - 1 
                  call mpi_send_int(pop_size_winner,1-myid,msg_num,ierr)
                  call mpi_send_int(pop_size_loser,1-myid,msg_num,ierr)
                  call mpi_send_int(id_winner,1-myid,msg_num,ierr)
               else
                  call mpi_recv_int(pop_size_winner,1-myid,msg_num,ierr)
                  call mpi_recv_int(pop_size_loser,1-myid,msg_num,ierr)
                  call mpi_recv_int(id_winner,1-myid,msg_num,ierr)
               end if
               winner = .false.
               
               if(run_status < 0) then ! dying tribe
                  run_status = 0 ! resurrect tribe - reset the status as ok
                  is_parallel = .true.
                  call mpi_recv_int(current_pop_size,1-myid,msg_num,ierr)
                  call mpi_recv_int(num_migrate,1-myid,msg_num,ierr)
               else ! winning tribe... send half of its pop size to dying tribe
                  winner = .true.
                  write(6,*)'<font color=red>*** FISSION TRIBE ***</font>', myid
                  other_run_status = 0
                  num_migrate = (pop_size_winner - pop_size_loser)/2
                  current_pop_size = (pop_size_winner + pop_size_loser)/2
                  ! if odd, round down half the time and round up half the time
                  if(mod(current_pop_size,2)==1 .and. randomnum(1).gt.0.5) then
                     num_migrate = num_migrate+1 
                     current_pop_size = current_pop_size + 1
                  end if
                  call mpi_send_int(current_pop_size,1-myid,msg_num,ierr)
                  call mpi_send_int(num_migrate,1-myid,msg_num,ierr)
                  write(6,*) 'migrating half the tribe from',myid+1,' to ',2-myid
               end if
               
               ! now migrate half the population
               do k = 1, num_migrate
                  i = pop_size_winner-num_migrate + k 
                  j = pop_size_loser + k
                  !if(myid.eq.0) write(*,*) 'migrating: ',i, 'to:',j
                  id_loser = 1 - id_winner
                  call migrate_individual(id_winner,id_loser,i,j,dmutn,fmutn,nmutn, & 
                       lb_mutn_count,linkage_block_fitness,winner)
               end do
               
            end if
            
         else ! .not.tribal_fission
            
            ! This is the alternative treatment to fission
            ! turn off parallel and let the one run to completion
            if(other_run_status < 0)  then
               num_tribes = num_tribes - 1
               is_parallel = .false.
            end if
            
            if(run_status < 0) goto 30
            
         end if ! tribal_fission

      else if (num_tribes > 2) then

         ! receive status from every process in group 
         ! to check if one process has died

         call mpi_allreduce(run_status,global_run_status,1, &
              mpi_integer,mpi_sum,mycomm,ierr)
  
         if(sum(global_run_status) < 0) then

            ranks(1) = -sum(global_run_status)

            call mpi_comm_group(mycomm,oldgroup,ierr)
            call mpi_group_excl(oldgroup,1,ranks,newgroup,ierr)
            call mpi_comm_create(mycomm,newgroup,mycomm,ierr)

            if(run_status < 0) goto 30

            call mpi_comm_size(mycomm,num_tribes,ierr)

         end if
      end if
   end if
! END_MPI

!  Write diagnostic information to output files.

   current_pop_size = max(1, current_pop_size)

   call second(tin_diagnostics)

   if(mod(gen, 10) == 0 .or. (.not.bottleneck_yes .and. & 
      current_pop_size <= pop_size/20) .or. gen <= 10) then
      print_flag = .true.
   else
      print_flag = .false.
   end if

   if(num_contrasting_alleles > 0) &
      call diagnostics_contrasting_alleles(dmutn, nmutn, fmutn, &
           work_fitness,initial_allele_effects, ica_count, max_size, .false.)

   call diagnostics_history_plot(dmutn,nmutn,fmutn,lb_mutn_count, &
                ica_count,gen,print_flag,current_global_pop_size)

   if(gen <= 3 .or. mod(gen, 20) == 0) then

      if (verbosity == 2) then
   !     Output fitness of each individual in population.
         rewind(16)
         write(16,'("# individual",2x,"fitness")')
         do i=1,current_pop_size
            write(16,*) i, fitness(i)
         end do
         call flush(16)
      endif

   end if

   if(fitness_distrib_type == 1 .and. &
      mod(gen, diagnostic_gens) == 0 .and. verbosity > 0) then

      if(tracking_threshold /= 1.0) &
         call diagnostics_mutn_bins_plot(dmutn, fmutn, accum, gen)

      call diagnostics_near_neutrals_plot(dmutn, fmutn, &
                linkage_block_fitness, lb_mutn_count, gen)

      call diagnostics_selection(sorted_score,pheno_fitness, &
                                 total_offspring,gen)

   end if

   call second(tout_diagnostics)
   sec(7) = sec(7) + tout_diagnostics - tin_diagnostics
   call second(tin_diagnostics)

   ! shutdown if fixation reached
   if(polygenic_beneficials.and.polygenic_fixed.and.poly_stop_gen==gen) then
      goto 20
   endif

   if (polygenic_beneficials) then
      if (percent_pop_poly >= 99) then
         print*,'POLYGENICS: SHUTTING DOWN BECAUSE 99% OF POPULATION HAS ALLELE'        
         goto 20 ! shutdown
      endif
   endif

   if(mod(gen, plot_allele_gens)==0 .and. gen /= num_generations) &
       call diagnostics_polymorphisms_plot(dmutn, nmutn, fmutn, &
                             work_fitness, max_size, gen)

   call second(tout_diagnostics)
   sec(8) = sec(8) + tout_diagnostics - tin_diagnostics

   call second(tout_gen)
   tgen = tout_gen - tin_gen

   if (verbosity == 2) then
      write(22,'(i12,f17.7,2f19.7)') gen, tgen, time_offspring, &
                                     time_selection
      call flush(22)
   endif
  
   !START_MPI
   if(is_parallel) then
      call mpi_ravg(tgen,par_tgen,1)
      call mpi_ravg(time_offspring,par_time_offspring,3)
      call mpi_ravg(time_selection,par_time_selection,1)
      if (myid==0) then
         write(23,'(i12,f17.7,2f19.7)') gen, par_tgen, &
           par_time_offspring, par_time_selection
         if (myid==0.and.(mod(gen,10)==0.or.gen<4)) then
            write(*,'("iteration time: ",i6,"  milliseconds")') &
                 int(1000.*tgen)
         end if
         call flush(23)
      end if
   end if
   !END_MPI

   ! Monitor state file for shutdown flag.

   if (run_status >= 0) then
      npath = index(data_file_path,' ') - 1
      open(10, file=data_file_path(1:npath)//case_id//'.st8', &
               status='unknown')
      read(10,*) run_status
      close(10)

      ! Premature shutdown
      shutdown_gen = gen
      if (run_status == 1) then
         write(6,*) 'STATE: WRITING RESTART FILE & EXITING RUN'
         write_dump = .true.
         restart_dump_number = 8
         dump_number = restart_dump_number
         goto 20
      end if
   end if

   !      Write PPM data.

   ! do i=1,pop_size
   !    if (fitness(i) > 1) then 
   !       red = 255
   !       green = 0
   !       blue = 0
   !    else
   !       red = int(fitness(i)*255)
   !       if (red < 0) red = 0 
   !          green = red
   !          blue = red
   !       end if
   !    write(15,'(i4,$)') red,green,blue
   ! end do
   ! write(15,*)

!  For dynamic population sizes compute new pop_size

   if(pop_growth_model > 0) then
      if(pop_growth_model == 1) then
         pop_size = ceiling(pop_growth_rate*pop_size)
      else if (pop_growth_model == 2) then
         pop_size = ceiling(pop_size*(1. + pop_growth_rate* &
                    (1. - pop_size/carrying_capacity)))
      else 
         write(0,*) 'ERROR: pop_growth_model ', pop_growth_model, &
                    'not supported'
         stop   
      end if   

      this_size = int((1.1*reproductive_rate*pop_size &
                  *(1. - fraction_random_death)))
      if(this_size > max_size) then
         write(0,*) "OUT OF MEMORY! SHUTTING DOWN!"
         goto 20
      end if 

   end if

   call flush(6)
   call flush(9)

! End generation loop
end do ! gen

! Shutdown procedures

20 continue

! Perform diagnositics on initial contrasting alleles and
! polymorphisms and create file for plotting.

call second(tin_diagnostics)

if(polygenic_beneficials) then
   open (12, file=data_file_path(1:npath)//case_id//'.'//myid_str &
         //'.poly',status='unknown')
   do i=6,12,3
      write(i,*)
      write(i,*) '-----------------------------------------------------'
      write(i,*) 'POLYGENIC BENEFICIALS SUMMARY:'
      write(i,*)
      write(i,*) '# instance     gen_enter    gen_exit    duration    string_id'
      write(i,*)
      poly_not_selected = 0
      do j=1,num_polys_cumulative
         if(pmutn(j)%gen_exit==0) then
            duration = 0
         else
            duration = pmutn(j)%gen_exit - pmutn(j)%gen_enter
         endif
         ! If gen_exit is 0, it means that there was a target match
         ! but that offspring was not selected.
         if(pmutn(j)%gen_exit > 0) then
            write(i,*) j, pmutn(j)%gen_enter, pmutn(j)%gen_exit, &
                       duration, pmutn(j)%id
         elseif(pmutn(j)%gen_exit == 0) then
            poly_not_selected = poly_not_selected + 1
         endif
         if(plot_allele_gens > 1) then
            print *, "NOTE: gen_exit and duration parameters can ONLY be"
            print *, "      computed correctly when plot_allele_gens = 1"
         endif 
      end do
      write(i,*)
      write(i,'(a,a,f7.2)') 'Percentage of instances lost to drift within ',  &
            'the generation in which it occurred:',  &
            real(poly_not_selected)/num_polys_cumulative*100.
      write(i,*)
      if(polygenic_fixed) then
         write(i,*) '# fixed polygenic id:', poly_fixed_fmutn
      endif
      write(i,*)
      write(i,*) 'First_inst_gen, Last_inst_gen, Fix_gen, Total_inst'
      write(i,'(2i15,i9,i12)') poly_gen_first_instance, &
                               poly_gen_last_instance,  & 
                               poly_stop_gen-plot_allele_gens, & 
                               num_polys_cumulative
      write(i,*) '-----------------------------------------------------'
   end do
   close(12)
end if

if(num_contrasting_alleles > 0) then
   call diagnostics_contrasting_alleles(dmutn, nmutn, fmutn, &
          work_fitness, initial_allele_effects, ica_count, max_size, .true.)
   call tabulate_initial_alleles(dmutn, fmutn, &
          linkage_block_fitness, initial_allele_effects, max_size)
end if

if(tracking_threshold /= 1.0) then
   call diagnostics_polymorphisms_plot(dmutn,nmutn,fmutn, &
                       work_fitness, max_size, gen-1)
! if(recombination_model /= clonal)
! &      call diagnostics_heterozygosity(dmutn, fmutn)
end if

call second(tout_diagnostics)
sec(8) = sec(8) + tout_diagnostics - tin_diagnostics

! Write an output dump that contains the current set of parameter
! values, the stored mutation array mutn, and the linkage block
! fitness array linkage_block_fitness.

dump_number = dump_number + 1

if(write_dump) call write_output_dump(dmutn,nmutn,fmutn,lb_mutn_count, &
                    linkage_block_fitness,initial_allele_effects, &
                    shutdown_gen,myid_str)

30 continue

call second(tout)
sec(2) = tout - tin_run

call profile(6)
call profile(9)

! Close files.

close(4)
close(6)
close(7)
close(8)
close(9)
close(11)
! close(12)
close(13)
close(14)
! close(15)
if (verbosity == 2) then 
   close(16)
   close(19)
   close(22)
   close(26)
endif
close(25)
if (polygenic_beneficials) close(20)
if (num_contrasting_alleles > 0) close(30)

if (is_parallel) then
   close(14)
   close(17)
   close(18)
   if(verbosity==2) close(23)
   close(24)
   close(34)
   close(35)
end if

! START_MPI
! Wait for all processes to reach this point before shutting down all processes
if(am_parallel.and.tribal_competition) then
   call mpi_mybarrier()
   if(num_tribes_at_start > 2) then
      call MPI_GROUP_FREE(OLDGROUP,ierr)
      call MPI_GROUP_FREE(NEWGROUP,ierr)
   end if
   call MPI_COMM_FREE(MYCOMM)
   call mpi_myfinalize(ierr)
elseif (is_parallel) then
   call MPI_COMM_FREE(MYCOMM,ierr)
   call mpi_myfinalize(ierr)
endif
! END_MPI

! Write an empty file called success which is used by the unit testing
! framework to know if the run has successfully completed or not 
filename = 'success'
j = index(data_file_path,' ') - 1
open (99, file=data_file_path(1:npath)//filename,status='unknown')
write(99,*)
close(99)
i = chmod(filename,'644')

!deallocate(dmutn, nmutn, fmutn, lb_mutn_count, linkage_block_fitness, &
!           initial_allele_effects,pheno_fitness,work_fitness,fitness, &
!           sorted_score,available,global_run_status,pop_size_array)
deallocate(dmutn,nmutn,fmutn,lb_mutn_count,linkage_block_fitness)
deallocate(initial_allele_effects,pheno_fitness,work_fitness,fitness)
deallocate(sorted_score,available)
!deallocate(global_run_status,pop_size_array) ! both of these cause the error:
!forrtl: severe (153): allocatable array or pointer is not allocated
deallocate(gp)
stop

end program mendel

! START_MPI
subroutine compute_tribal_fitness(dmutn, fmutn, pop_size_array, &
           current_global_pop_size, gen)
use inputs
use random_pkg
use selection_module
include 'common.h'
include 'mpif.h'
integer dmutn(max_del_mutn_per_indiv/2,2,*)
integer fmutn(max_fav_mutn_per_indiv/2,2,*)
integer pop_size_array(num_tribes), current_global_pop_size
integer i, j, k, m, gen
real*8 fitness, decode_fitness_del, decode_fitness_fav
real*8 par_pre_sel_fitness, par_post_sel_fitness, &
       mod_par_post_sel_fitness, &
       post_sel_fitness_array(num_tribes)
real*8 tribal_fitness_variance, par_tribal_fitness

! Compute total social bonus.
 
if(upload_mutations .and. altruistic) then
   social_bonus = 0.d0
   do i=1, current_pop_size
      do j = 1, 2
         do k = 2, dmutn(1,j,i)
            do m = 1, num_uploaded_mutn
               if (dmutn(k,j,i) == uploaded_mutn(m)) then
                  fitness = abs(decode_fitness_del(dmutn(k,j,i)))
                  social_bonus = social_bonus + fitness
               end if
            end do
         end do
      end do
   end do
   social_bonus = social_bonus / current_pop_size
else
   social_bonus = 0.d0
end if

! Compute tribal fitness

if (is_parallel) then

   call mpi_davg(post_sel_fitness,par_post_sel_fitness,1)
   call mpi_davg(pre_sel_fitness,par_pre_sel_fitness,1)
   call mpi_mybcastd(par_post_sel_fitness,1)

   if(tribal_competition) then

!     Gather the fitnesses from each tribe into a single 
!     array called post_sel_fitness_array in order to compute
!     tribal fitness variance below.

      call MPI_GATHER(post_sel_fitness,1,MPI_DOUBLE_PRECISION, &
                      post_sel_fitness_array,1, &
                      MPI_DOUBLE_PRECISION,0,MYCOMM,ierr)

      call mpi_isum(current_pop_size,current_global_pop_size,1)

      if(myid==0) then

!        Compute global weighted genetic fitness.

         global_genetic_fitness=0.d0
         do i = 1, num_tribes
            global_genetic_fitness = &
                    pop_size_array(i)*post_sel_fitness_array(i) + &
                                     global_genetic_fitness
         end do 
         global_genetic_fitness = global_genetic_fitness/ &
                                  dble(current_global_pop_size)

         if(mod(gen,10)==0) then
            write(*,'("global genetic fitness = ",f7.4)') &
                       global_genetic_fitness
         end if

         tribal_fitness_variance=0.d0
         do i = 1, num_tribes
            tribal_fitness_variance=tribal_fitness_variance+ &
                           (par_post_sel_fitness- &
                            post_sel_fitness_array(i))**2
         end do

!        Compute the tribal noise variance required to yield the 
!        specified group heritability.  Take the square root to 
!        obtain the standard deviation.

         tribal_noise = sqrt(tribal_fitness_variance* &
                             (1.d0 - group_heritability) &
                              /group_heritability)
      end if

!     Broadcast this value to all tribes.

      call mpi_mybcastd(global_genetic_fitness,1)
      call mpi_mybcastd(tribal_noise,1)

!     Add noise component to tribal_fitness (post_sel_fitness)
!     Add a tiny variable positive increment to eliminate identical
!     fitness values when the noise is zero.

      tribal_noise = tribal_noise*random_normal() + 1.d-15*myid

!     Compute the total tribal fitness.

      tribal_fitness = post_sel_fitness + tribal_noise 

!     The tribal_fitness_factor relates the fitness of the
!     current tribe to the average tribal fitness from all tribes.  

      call mpi_davg(tribal_fitness,par_tribal_fitness,1)
      call mpi_mybcastd(par_tribal_fitness,1)
      
      tribal_fitness_factor=tribal_fitness/ &
                            global_genetic_fitness

!     Add social_bonus into the total tribal fitness.

      tribal_fitness = tribal_fitness + &
                       social_bonus_factor*social_bonus

!     Recompute tribal_fitness_factor which now includes 
!     social_bonus.

      call mpi_davg(tribal_fitness,par_tribal_fitness,1)
      call mpi_mybcastd(par_tribal_fitness,1)
      tribal_fitness_factor=tribal_fitness/par_tribal_fitness

   else

!     Add social_bonus to fitness for the case when parallel is
!     turned on but tribal competition is turned off.

      tribal_fitness = tribal_fitness +  &
                       social_bonus_factor*social_bonus

   endif

else 

!  Add social bonus to fitness for non-parallel/single-tribe case.
   tribal_fitness = tribal_fitness +  &
                    social_bonus_factor*social_bonus

end if

end subroutine compute_tribal_fitness
! END_MPI
