program mendel

use init
use inputs
use genome
use profiler
use polygenic
use random_pkg
use selection_module
include 'common.h'

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
integer :: total_offspring, offspring_count, empty, red, green, blue
integer :: ica_count(3), cumulative_offspring, mutn_indx
integer :: max_pop_size, pop_size_allocation, current_global_pop_size
integer :: global_pop, mean_pop, delta_pop, delta, nie, mm, id
integer :: pop_size_winner, pop_size_loser, num_migrate
integer :: id_winner, id_loser
integer :: num_dmutns, num_fmutns, encode_mutn, string(40)
integer :: OLDGROUP,NEWGROUP,ranks(1),num_tribes_at_start

real*8 accum(50), reproductive_advantage_factor
real selection_coefficient, aoki, migration_rate
real carrying_capacity, x
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
!  Pass in max_pop_size through num_generations input parameter 
   max_pop_size = num_generations
   pop_size_allocation = max_pop_size       
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


!  Write diagnostic information to output files.

   current_pop_size = max(1, current_pop_size)

   call second(tin_diagnostics)

   if(mod(gen, 10) == 0 .or. (.not.bottleneck_yes .and. & 
      current_pop_size <= pop_size/20) .or. gen < 4) then
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
!        pass carrying capacity through num_generations
         carrying_capacity = num_generations
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
      write(i,*) 'First_inst_gen, Last_inst_gen, Fix_gen, Total_inst'
      write(i,'(2i15,i9,i12)') poly_gen_first_instance, &
                               poly_gen_last_instance,  & 
                               poly_stop_gen-plot_allele_gens, & 
                               num_polys_cumulative
      write(i,*) '# instance     gen_enter    gen_exit    duration     string_id'
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

subroutine diagnostics_history_plot(dmutn, nmutn, fmutn, &
           lb_mutn_count, ica_count, gen, print_flag, &
           current_global_pop_size)
use selection_module
use random_pkg
use polygenic
use inputs
use init
include 'common.h'
integer, parameter :: MNP=100000
integer dmutn(max_del_mutn_per_indiv/2,2,*)
integer nmutn(max_neu_mutn_per_indiv/2,2,*)
integer fmutn(max_fav_mutn_per_indiv/2,2,*)
integer lb_mutn_count(num_linkage_subunits,2,3,*)
integer gen, i, j, k, m, lb, num_recessive, dth, jmax
integer ica_count(3)
integer current_global_pop_size, global_num_back_mutn
real*8 total_fav_mutn, par_total_fav_mutn
real*8 total_neu_mutn, par_total_neu_mutn
real*8 par_tracked_fav_mutn
real*8 tracked_neu_mutn, par_tracked_neu_mutn
real*8 par_total_del_mutn, par_tracked_del_mutn
real*8 frac_recessive, frac_accum, total_mutn, st
real*8 par_pre_sel_fitness, par_post_sel_fitness, &
       par_pre_sel_geno_sd, par_pre_sel_pheno_sd, &
       par_pre_sel_corr, par_post_sel_geno_sd,    &
       par_post_sel_pheno_sd, par_post_sel_corr,  &
       mod_par_post_sel_fitness,                  &
       post_sel_fitness_array(num_tribes)
real   suma, sumc, sumg, sumt
logical print_flag

total_del_mutn   = 0
total_neu_mutn   = 0
total_fav_mutn   = 0

tracked_del_mutn = 0
tracked_neu_mutn = 0
tracked_fav_mutn = 0

do i=1,current_pop_size
   do lb=1,num_linkage_subunits
      total_del_mutn = total_del_mutn          &
                     + lb_mutn_count(lb,1,1,i) &
                     + lb_mutn_count(lb,2,1,i)
      total_fav_mutn = total_fav_mutn          & 
                     + lb_mutn_count(lb,1,2,i) &
                     + lb_mutn_count(lb,2,2,i)
      total_neu_mutn = total_neu_mutn          &
                     + lb_mutn_count(lb,1,3,i) &
                     + lb_mutn_count(lb,2,3,i)
   end do
   tracked_del_mutn = tracked_del_mutn + dmutn(1,1,i) &
                                       + dmutn(1,2,i)
   tracked_fav_mutn = tracked_fav_mutn + fmutn(1,1,i) &
                                       + fmutn(1,2,i)
   tracked_neu_mutn = tracked_neu_mutn + nmutn(1,1,i) &
                                       + nmutn(1,2,i)
end do

tracked_del_mutn = tracked_del_mutn - ica_count(1)
tracked_fav_mutn = tracked_fav_mutn - ica_count(2)
tracked_neu_mutn = tracked_neu_mutn - ica_count(3)

! Compute averages across processors.


num_recessive  = 0
frac_recessive = 0.

do i=1,current_pop_size
   do j=2,dmutn(1,1,i)+1
      if(dmutn(j,1,i) < 0) num_recessive = num_recessive + 1
   end do
   do j=2,dmutn(1,2,i)+1
      if(dmutn(j,2,i) < 0) num_recessive = num_recessive + 1
   end do
end do

if(tracked_del_mutn > 0) frac_recessive = &
                         real(num_recessive)/tracked_del_mutn

! Output to a file for plotting the generation number, the mean 
! fitness, the average number of mutations per individual, and 
! the number of fixed favorable mutations.

if(mod(gen,hst_gens)==0) then
   if(tribal_competition) then
   !  Unit 7 are the .001, .002, etc. invidual tribal files
      write(7,'(i12,2e16.4,1p3e14.4,3e14.4)') gen, post_sel_fitness, &
         post_sel_geno_sd, &
         total_del_mutn/current_pop_size, &
         total_fav_mutn/current_pop_size, &
         total_neu_mutn/current_pop_size, &
         current_pop_size/real(global_pop_size)*100., &
         global_genetic_fitness, fertility_factor
   else 
      write(7,'(i12,2e16.4,1p3e14.4,i12)') gen,  &
         post_sel_fitness, post_sel_geno_sd, &
         total_del_mutn/current_pop_size, &
         total_fav_mutn/current_pop_size, &
         total_neu_mutn/current_pop_size, current_pop_size
   endif
   call flush(7)
endif


if (.not.print_flag) return

if(mod(gen,output_gens)==0 .or. percent_pop_poly >= 99.) then
   call write_status(9, gen, current_pop_size,            &
        frac_recessive, total_del_mutn, tracked_del_mutn, &
        total_fav_mutn, total_neu_mutn, pre_sel_fitness,  &
        pre_sel_geno_sd, pre_sel_pheno_sd, pre_sel_corr,  &
        post_sel_fitness, post_sel_geno_sd, post_sel_pheno_sd, &
        post_sel_corr, num_polys_this_gen, num_polys_cumulative)
endif
 
if(allow_back_mutn) write(9,"('mean number of back mutations', &
   '/indiv =',f10.2)") real(num_back_mutn)/real(current_pop_size)

if(is_parallel) then

   if (myid == 0) then 
      call write_status(6, gen, current_global_pop_size, &
           frac_recessive, par_total_del_mutn,           &
           par_tracked_del_mutn, par_total_fav_mutn,     &
           par_tracked_neu_mutn,                         &
           par_pre_sel_fitness, par_pre_sel_geno_sd,     &
           par_pre_sel_pheno_sd, par_pre_sel_corr,       &
           par_post_sel_fitness, par_post_sel_geno_sd,   &
           par_post_sel_pheno_sd, par_post_sel_corr,     &
           num_polys_this_gen)

      if(allow_back_mutn) write(6,"('mean number of back ' &
         'mutations/indiv =',f10.2)") real(global_num_back_mutn) &
         /real(current_global_pop_size)
      
   end if

   if((gen < 4 .or. mod(gen,10)==0) .and. altruistic) then
      write(*,'(a,i3,a,f12.8,$)') 'tribe:', myid+1, &
                                  ' social_bonus:', social_bonus
      if(tribal_competition) then
         write(*,'(a,f12.8)') ' tribal_noise:', tribal_noise
      else
         write(*,*)
      end if
   end if

else

   if(mod(gen,output_gens)==0 .or. percent_pop_poly >= 99.) then
      call write_status(6, gen, current_pop_size,                 &
           frac_recessive, total_del_mutn, tracked_del_mutn,      &
           total_fav_mutn, total_neu_mutn, pre_sel_fitness,       &
           pre_sel_geno_sd, pre_sel_pheno_sd, pre_sel_corr,       &
           post_sel_fitness, post_sel_geno_sd, post_sel_pheno_sd, &
           post_sel_corr, num_polys_this_gen)

      if(polygenic_beneficials) then
         suma = 0.
         sumc = 0.
         sumg = 0.
         sumt = 0.
         do i=1,current_pop_size
            jmax = 2
            if(recombination_model==clonal) jmax = 1
            do j=1,jmax
               do k=2,num_linkage_subunits + 1
                  if(nmutn(k,j,i) == A) suma = suma + 1.
                  if(nmutn(k,j,i) == C) sumc = sumc + 1.
                  if(nmutn(k,j,i) == G) sumg = sumg + 1.
                  if(nmutn(k,j,i) == T) sumt = sumt + 1.
               end do
            end do
         end do
         suma = suma/(num_linkage_subunits*jmax*current_pop_size)
         sumc = sumc/(num_linkage_subunits*jmax*current_pop_size)
         sumg = sumg/(num_linkage_subunits*jmax*current_pop_size)
         sumt = sumt/(num_linkage_subunits*jmax*current_pop_size)
         write(6,'(a,4f6.3)') 'Average nucleotide frequencies (A C G T): ',  &
               suma, sumc, sumg, sumt
         if(poly_debug) write(6,*) 'polygenics:', pmutn(1:num_polys_this_gen)%id
         i =   current_pop_size/2
         j = 9*current_pop_size/10
         if(recombination_model == clonal) then
            write(6,'(a,1x,a,1x,a)')  &
               'DNA of 50th and 90th percentile individuals: ', &
               poly_to_string(nmutn(2:,1,i)), poly_to_string(nmutn(2:,1,j))
         else
            write(6,'(a,1x,a,1x,a,3x,a,1x,a)')  &
               'DNA of 50th and 90th percentile individuals: ', &
               poly_to_string(nmutn(2:,1,i)), poly_to_string(nmutn(2:,2,i)), &
               poly_to_string(nmutn(2:,1,j)), poly_to_string(nmutn(2:,2,j))
         end if
         write(6,'(a,1x,f7.1)') 'Percentage of population with target:', &
            percent_pop_poly
      endif

   endif

   if(allow_back_mutn) write(6,"('mean number of back mutations', &
   '/indiv =',f10.2)") real(num_back_mutn)/real(current_pop_size)

   if((gen < 4 .or. mod(gen,10)==0) .and. altruistic) then
      write(*,'(a,f12.8)') 'social_bonus: ',social_bonus
   end if

end if

if(tracking_threshold == 1.0 .and. gen >= 200 .and. &
   mod(gen,diagnostic_gens) == 0) then
   total_mutn = gen*current_pop_size*mutn_rate
   if(frac_fav_mutn /= 1. .or. polygenic_beneficials) then
      frac_accum = total_del_mutn/(total_mutn*(1.-frac_fav_mutn))
      st = 1.3*exp(-alpha_del*(1. - frac_accum)**gamma_del)
      write(9, '("deleterious selection threshold  =",1pe10.3)') st
      write(9, '("deleterious fraction accumulated =",f7.4)') frac_accum
      if(.not. is_parallel) then
         write(6, '("deleterious selection threshold  =",1pe10.3)') st
         write(6, '("deleterious fraction accumulated =",f7.4)') frac_accum
      end if
      write(25, '(i10,1p1e15.3,2a15)') gen, st, 'NaN', 'NaN' 
      if(is_parallel .and. myid==0) then
         write(35, '(i10,1p1e15.3,2a15)') gen, st, 'NaN', 'NaN' 
      end if
   end if
end if

if(num_contrasting_alleles > 0) &
   write(9,'(19x,"Statistics for initial contrasting alleles"/,  &
             " mean fitness contrib =",f7.4,"  fav mean freq =", &
             f7.4,"  fixed =",i4,"  lost =",i4)')                &
             ica_mean_effect, fav_mean_freq, fav_fixed, fav_lost

if(num_contrasting_alleles > 0 .and. .not. is_parallel)          &
   write(6,'(19x,"Statistics for initial contrasting alleles"/,  &
             " mean fitness contrib =",f7.4,"  fav mean freq =", &
             f7.4,"  fixed =",i4,"  lost =",i4)')                &
             ica_mean_effect, fav_mean_freq, fav_fixed, fav_lost

end subroutine diagnostics_history_plot

subroutine diagnostics_mutn_bins_plot(dmutn, fmutn, accum, gen)
use inputs
include 'common.h'
integer dmutn(max_del_mutn_per_indiv/2,2,*)
integer fmutn(max_fav_mutn_per_indiv/2,2,*)
integer gen, i, j, k, k0, accum_gen
integer fid, oneortwo
real*8 fitness_bins(200,2), par_fitness_bins(100,2), work(100,2)
real*8 bin_fitness_boxwidth(101), bin_fitness_midpoint(101)
real*8 refr_bins(100), bin_fitness(101), del_bin_width, b0, b1
real*8 d, x, x0, x1, y0, y1, s, mutn_sum, fav_bin_width
real*8 av1, av2, fm1, fm2, sum, accum(50), del, del_no_sel, ratio
real*8 del_dom_thres, del_rec_thres, fav_dom_thres, fav_rec_thres
real*8 par_del_dom_thres, par_fav_dom_thres

! Compute the total number of mutations and the bin widths.
   
mutn_sum  = current_pop_size*gen*mutn_rate
del_bin_width = -log(tracking_threshold)/50
if(max_fav_fitness_gain > 0.) then
   fav_bin_width = -log(tracking_threshold/max_fav_fitness_gain)/50
else
   fav_bin_width = del_bin_width
end if

x0 = 0.
y0 = 0.
do k=1,50
   x1 = (del_bin_width*k/alpha_del)**(1./gamma_del)
   refr_bins(k) = (1. - frac_fav_mutn)*mutn_sum*(x1 - x0)
   y1 = (fav_bin_width*k/alpha_fav)**(1./gamma_fav)
   refr_bins(50+k) = frac_fav_mutn*mutn_sum*(y1 - y0)
   x0 = x1
   y0 = y1
end do

! Compute statistics on favorable and recessive mutations.

fitness_bins = 0.

do i=1,current_pop_size

   do j=2,dmutn(1,1,i)+1
      x = mod(abs(dmutn(j,1,i)),lb_modulo)*del_scale
      d = alpha_del*x**gamma_del
      k = 1 + int(d/del_bin_width)
      if(dmutn(j,1,i) < 0) then
         if(k <= 50) then
            fitness_bins(k,1) = fitness_bins(k,1) + 1.
         end if
      else
         if(k <= 50) then
            fitness_bins(k,2) = fitness_bins(k,2) + 1.
         end if
      end if
   end do

   do j=2,dmutn(1,2,i)+1
      x = mod(abs(dmutn(j,2,i)),lb_modulo)*del_scale
      d = alpha_del*x**gamma_del
      k = 1 + int(d/del_bin_width)
      if(dmutn(j,2,i) < 0) then
         if(k <= 50) then
            fitness_bins(k,1) = fitness_bins(k,1) + 1.
         end if
      else
         if(k <= 50) then
            fitness_bins(k,2) = fitness_bins(k,2) + 1.
         end if
      end if
   end do

   do j=2,fmutn(1,1,i)+1
      d = alpha_fav*(mod(abs(fmutn(j,1,i)),lb_modulo)*fav_scale)**gamma_fav
      k = 51 + int(d/fav_bin_width)

      if(fmutn(j,1,i) < 0) then
         if(k <= 100) then
            fitness_bins(k,1) = fitness_bins(k,1) + 1.
         end if
      else
         if(k <= 100) then
            fitness_bins(k,2) = fitness_bins(k,2) + 1.
         end if
      end if
   end do

   do j=2,fmutn(1,2,i)+1
      d = alpha_fav*(mod(abs(fmutn(j,2,i)),lb_modulo)*fav_scale)**gamma_fav
      k = 51 + int(d/fav_bin_width)

      if(fmutn(j,2,i) < 0) then
         if(k <= 100) then
            fitness_bins(k,1) = fitness_bins(k,1) + 1.
         end if
      else
         if(k <= 100) then
            fitness_bins(k,2) = fitness_bins(k,2) + 1.
         end if
      end if
   end do

end do

! Compute fitness values for bin boundaries and bin centers.

do k=1,51
   bin_fitness(k) = exp(-del_bin_width*(k - 1))
   if (k > 1) then
      bin_fitness_boxwidth(k-1) =  abs(bin_fitness(k) - bin_fitness(k-1))
      bin_fitness_midpoint(k-1) =  (bin_fitness(k) + bin_fitness(k-1))/2.
   end if
end do

do k=51,101
   bin_fitness(k) = max_fav_fitness_gain*exp(-fav_bin_width*(k - 51))
   if (k > 51) then
      bin_fitness_boxwidth(k-1) =  abs(bin_fitness(k) - bin_fitness(k-1))
      bin_fitness_midpoint(k-1) =  (bin_fitness(k) + bin_fitness(k-1))/2.
   end if
end do

! Compute the frequency of output of the mutation accumulation
! statistics.

accum_gen = 1000000
if(mod(gen, 100) == 0) then
   if(gen <= 500) then
      accum_gen = 100
   elseif(gen <=  1000) then
      accum_gen = 500
   elseif(gen <=  5000) then
      accum_gen = 1000
   elseif(gen <= 10000) then
      accum_gen = 5000
   else
      accum_gen = 10000
   end if
   if(gen == 100) accum = 0.
end if

! Output the number of accumulated deleterious dominant mutations 
! in each of 50 bins during the previous accum_gen generations at
! appropriate intervals, along with other relevant information.

if(mod(gen, accum_gen) == 0 .and. verbosity == 2) then 

   if(gen == 100) then
   write(26,'("#"/"#",11x, "Generation    Accumulation Interval"/ &
      2i19/"#"/"#           Accumulation Over Previous Interval"/ &
      "#"/"# bin  fitness effect    actual      expected       "  &
      "ratio  expected fraction")') gen, accum_gen
   else
   write(26,'("#"/"#",11x, "Generation    Accumulation Interval"/ &
      2i19/"#"/"#           Accumulation Over Previous Interval"/ &
      "#"/"# bin  fitness effect    actual      expected       "  &
      "ratio      total accum")') gen, accum_gen
   end if

   do k=1,50
      del = fitness_bins(k,2) - accum(k)
      del_no_sel = refr_bins(k)*(1. - fraction_recessive)*accum_gen/gen
      ratio = del/del_no_sel
      if(gen == 100) then
         x0    = refr_bins(k)*(1. - fraction_recessive)/mutn_sum
         write(26,'(i4,1pe14.3,i12,0pf14.2,f14.7,1pe14.3)') k,       &
            bin_fitness_midpoint(k), int(del), del_no_sel, ratio, x0
      else
      write(26,'(i4,1pe14.3,i12,0pf14.2,f14.7,i14)') k,           &
         bin_fitness_midpoint(k), int(del), del_no_sel, ratio,    &
         int(fitness_bins(k,2))
      end if
   end do

   accum = fitness_bins(1:50,2)

   call flush(26)

end if

! Normalize the binned mutations by the reciprocal of the expected
! number of mutations per bin in the absence of selection.

x = 1. - fraction_neutral
if (x.eq.0) x = 1. ! don't scale data if fraction_neutral = 1
do k=1,100

   if(refr_bins(k) > 0. .and. fraction_recessive > 0.) then
      fitness_bins(k,1) = fitness_bins(k,1) &                  
                        /(fraction_recessive*refr_bins(k))/x
   else
      fitness_bins(k,1) = 0.
   end if

   if(refr_bins(k) > 0. .and. fraction_recessive < 1.) then
      fitness_bins(k,2) = fitness_bins(k,2) &
                       /((1. - fraction_recessive)*refr_bins(k))/x
   else
      fitness_bins(k,2) = 0.
   end if

end do

! Perform an iteration of smoothing on the fitness_bin values 
! using a three-point average.  Iterate three times.

do i=1,3
fm1 = fitness_bins(1,1)
fm2 = fitness_bins(1,2)
do k=2,49
  av1 = fitness_bins(k,1) + 0.5*(fm1 + fitness_bins(k+1,1))
  fm1 = fitness_bins(k,1)
  work(k,1) = 0.5*av1
  av2 = fitness_bins(k,2) + 0.5*(fm2 + fitness_bins(k+1,2))
  fm2 = fitness_bins(k,2)
  work(k,2) = 0.5*av2
end do
fitness_bins(50,:) = 0.5*(fitness_bins(49,:) + fitness_bins(50,:))
fitness_bins(2:49,:) = work(2:49,:)
end do

! For favorable distribution, limit maximum to a value of 100.
! To increase the smoothness, iterate the smoothing two times.

fitness_bins(51:100,:) = min(100., fitness_bins(51:100,:))

do i=1,2
fm1 = fitness_bins(51,1)
fm2 = fitness_bins(51,2)
do k=52,99
  av1 = fitness_bins(k,1) + 0.5*(fm1 + fitness_bins(k+1,1))
  fm1 = fitness_bins(k,1)
  work(k,1) = 0.5*av1
  av2 = fitness_bins(k,2) + 0.5*(fm2 + fitness_bins(k+1,2))
  fm2 = fitness_bins(k,2)
  work(k,2) = 0.5*av2
end do
fitness_bins(52:99,:) = work(52:99,:)
end do

! Write the fitness bin information for deleterious mutations
! to standard output.

if (.not. is_parallel) then
   write(6, '(14x,"Fraction of mutations retained versus fitness effect")')
   write(6, '(" effect:",10(x,1pe7.0))') (bin_fitness_midpoint(k),k=3,50,5)
   write(6, '(" domint:",10f7.4)') (fitness_bins(k,2),k=3,50,5) 
   if(fraction_recessive > 0.) write(6, '(" recess:",0p10f7.4)') &
             (fitness_bins(k,1),k=3,50,5)
end if

write(9, '(14x,"Fraction of mutations retained versus fitness effect")')
write(9, '(" effect:",10(x,1pe7.0))') (bin_fitness_midpoint(k),k=3,50,5)
write(9, '(" domint:",10f7.4)') (fitness_bins(k,2),k=3,50,5) 
if(fraction_recessive > 0.) write(6, '(" recess:",0p10f7.4)') &
          (fitness_bins(k,1),k=3,50,5)

rewind (8)

write(8,'("# generation = ",i8)') gen
write(8,'("# deleterious mutations")')
write(8,'("# bin_fitness",3x,"recessive",2x,"dominant",3x,"box_width")')
do k=1,50
   write(8, '(1pe13.5,0p2f11.5,1pe13.5)') bin_fitness_midpoint(k), &
      fitness_bins(k,1), fitness_bins(k,2), bin_fitness_boxwidth(k)
end do

write(8,'("# favorable mutations")')
write(8,'("# bin_fitness",3x,"recessive",2x,"dominant",3x,"box_width")')
do k=51,100
   write(8, '(1pe13.5,0p2f11.5,1pe13.5,e13.5)') bin_fitness_midpoint(k), &
      fitness_bins(k,1), fitness_bins(k,2), bin_fitness_boxwidth(k),     &
      refr_bins(k)
end do
call flush(8)


! Compute the current values of the selection thresholds for both
! dominant and recessive deleterious mutations.

! Compute estimate for the current deleterious dominant threshold.

k   = 1
k0  = 0
sum = (1. - fraction_recessive)*refr_bins(50)
del_dom_thres = 0.
 
do while(k <= 50 .and. sum > 5000) 
   if(fitness_bins(k,2) > 0.25 .and. k0 == 0) k0 = k
   if(fitness_bins(k,2) > 0.75 .and. k > k0 + 1) then
      x0 = 0.
      y0 = 0.
      do i=k0,k-1
         x0 = x0 + i - 0.5
         y0 = y0 + fitness_bins(i,2)
      end do
      x0 = x0/(k - k0)
      y0 = y0/(k - k0)
      s = 0.
      d = 0.
      do i=k0,k-1
         s = s + (i - 0.5 - x0)*(fitness_bins(i,2) - y0)
         d = d + (i - 0.5 - x0)**2
      end do
      b1 = s/d
      b0 =   y0 - b1*x0
      x1 = (0.5 - b0)/b1
      del_dom_thres = exp(-x1*del_bin_width)
      k = 50
   end if
   k = k + 1
end do

! Compute estimate for the current deleterious recessive threshold.

k   = 1
k0  = 0
sum = fraction_recessive*refr_bins(50)
del_rec_thres = 0.
 
do while(k <= 50 .and. sum > 5000) 
   if(fitness_bins(k,1) > 0.25 .and. k0 == 0) k0 = k
   if(fitness_bins(k,1) > 0.75 .and. k > k0 + 1) then
      x0 = 0.
      y0 = 0.
      do i=k0,k-1
         x0 = x0 + i - 0.5
         y0 = y0 + fitness_bins(i,1)
      end do
      x0 = x0/(k - k0)
      y0 = y0/(k - k0)
      s = 0.
      d = 0.
      do i=k0,k-1
         s = s + (i - 0.5 - x0)*(fitness_bins(i,1) - y0)
         d = d + (i - 0.5 - x0)**2
      end do
      b1 = s/d
      b0 =   y0 - b1*x0
      x1 = (0.5 - b0)/b1
      del_rec_thres = exp(-x1*del_bin_width)
      k = 50
   end if
   k = k + 1
end do

! Compute estimate for the current favorable dominant threshold.

! First find the bin with the maximum ratio of actual to expected
! mutations if there were no selection, but restricted to ratios
! less than 3.

if (frac_fav_mutn > 0 .or. polygenic_beneficials) then

   y0 = fitness_bins(96,2)
   k0 = 4
   k  = 5
   do while(k <= 50)
      if(fitness_bins(101-k,2) > y0) then
         y0 = fitness_bins(101-k,2)
         k0 = k
      end if
      if(fitness_bins(101-k,2) >= 3.) k = 50
      k = k + 1
   end do

!  Now find the first bin with k < k0 that is backeted by ratios 
!  below and above 2.0.

   j = k0 - 1
   do k=k0-1,k0-5,-1
      if(fitness_bins(100-k,2) >  2.0 .and. &
         fitness_bins(101-k,2) <= 2.0) then
         j = k
      end if
   end do

end if

sum = (1. - fraction_recessive)*refr_bins(100)
fav_dom_thres = 0.

if(sum > 2000 .and. y0 > 2.5) then

!  Use simple linear interpolation to find the fitness effect
!  value corresponding to the ratio of 2.0.

   s  = (fitness_bins(100-j,2) - fitness_bins(101-j,2))/ &
        (bin_fitness_midpoint(100-j) &
      -  bin_fitness_midpoint(101-j))
   y0 =  fitness_bins(101-j,2) - s*bin_fitness_midpoint(101-j)
   fav_dom_thres = max(0., (2.0 - y0)/s)

end if

if(.not. is_parallel) then
   if(del_dom_thres > 0.) then
      x0 = (-log(del_dom_thres)/alpha_del)**(1/gamma_del)
      write(6, '("deleterious selection threshold   =",1pe10.3)') del_dom_thres
      write(6, '("deleterious fraction unselectable =",f6.3)') 1. - x0
   end if
   if(fav_dom_thres > 0.) then
      x0 = (-log(fav_dom_thres/max_fav_fitness_gain)/alpha_fav)**(1/gamma_fav)
      write(6, '("  favorable selection threshold   =",1pe10.3)') fav_dom_thres
      write(6, '("  favorable fraction unselectable =",f6.3)') 1. - x0
   end if
else
   if(myid == 0 .and. par_fav_dom_thres > 0.) then
      x0 = (-log(par_fav_dom_thres/max_fav_fitness_gain) &
           /alpha_fav)**(1/gamma_fav)
      write(6, '("  favorable selection threshold   =",1pe10.3)')  &
            par_fav_dom_thres
      write(6, '("  favorable fraction unselectable =",f6.3)') 1. - x0
   end if
end if

if(del_dom_thres > 0.) then
   x0 = (-log(del_dom_thres)/alpha_del)**(1/gamma_del)
   write(9, '("deleterious selection threshold   =",1pe10.3)') del_dom_thres
   write(9, '("deleterious fraction unselectable =",f6.3)') 1. - x0
end if
if(fav_dom_thres > 0.) then
   x0 = (-log(fav_dom_thres/max_fav_fitness_gain)/alpha_fav)**(1/gamma_fav)
   write(9, '("  favorable selection threshold   =",1pe10.3)') fav_dom_thres
   write(9, '("  favorable fraction unselectable =",f6.3)') 1. - x0
end if

if(is_parallel.and.myid==0) then
   oneortwo=2
else
   oneortwo=1
endif

fid=25

do i=1,oneortwo
   if(oneortwo==2) then
      fid=35
      del_dom_thres = par_del_dom_thres
      fav_dom_thres = par_fav_dom_thres
   end if
   
!  Do not write out zero threshold values, instead use NaN's
   if (del_dom_thres==0 .and. del_rec_thres==0 .and. fav_dom_thres==0) then
      write(fid, '(i10,3a15)') gen, 'NaN', 'NaN', 'NaN' 
   else if (del_dom_thres == 0 .and. del_rec_thres == 0) then
      write(fid, '(i10,2a15,1p1e15.3)') gen, 'NaN', 'NaN', fav_dom_thres
   else if (del_dom_thres == 0 .and. fav_dom_thres == 0) then
      write(fid, '(i10,a15,1p1e15.3,a15)') gen, 'NaN', del_rec_thres, 'NaN' 
   else if (del_rec_thres == 0 .and. fav_dom_thres == 0) then
      write(fid, '(i10,1p1e15.3,2a15)') gen, del_dom_thres, 'NaN', 'NaN' 
   else if(del_dom_thres == 0) then
      write(fid, '(i10,a15,1p2e15.3)') gen, 'NaN', del_rec_thres, fav_dom_thres
   else if(del_rec_thres == 0) then
      write(fid, '(i10,1p1e15.3,a15,1p1e15.3)') &
            gen, del_dom_thres, 'NaN', fav_dom_thres
   else if (fav_dom_thres == 0) then
      write(fid, '(i10,1p2e15.3,a15)') &
            gen, del_dom_thres, del_rec_thres, 'NaN'
   else
      write(fid, '(i10,1p3e15.3)') &
            gen, del_dom_thres, del_rec_thres, fav_dom_thres
   end if

   call flush(fid)

end do

end subroutine diagnostics_mutn_bins_plot

subroutine diagnostics_near_neutrals_plot(dmutn, fmutn, &
           linkage_block_fitness, lb_mutn_count, gen)
use selection_module
use inputs
include 'common.h'
integer dmutn(max_del_mutn_per_indiv/2,2,*)
integer fmutn(max_fav_mutn_per_indiv/2,2,*)
integer lb_mutn_count(num_linkage_subunits,2,3,*)
real*8 linkage_block_fitness(num_linkage_subunits,2,*)
integer gen, i, j, k, lb
real*8 bins_mutns(200,2), expn_bins(200)
real*8 haplotype_bins(200), haplotype_bin_width, avg_lb_effect
real*8 lb_fitness_frac_positive, num_del_lb, num_fav_lb
real*8 par_expn_bins(200)
real*8 par_haplotype_bins(200), par_avg_lb_effect
real*8 par_lb_fitness_frac_pos, y0
real*8 d, x, x0, x0r, x1, z0, z1, del_bin_width, fav_bin_width

! Generate the theoretical distribution curves for plotting.
   
haplotype_bin_width = 1.e-04
del_bin_width = haplotype_bin_width
fav_bin_width = 0.01*min(0.01, max_fav_fitness_gain)

x0 = 1.
do k=100,1,-1
   x1 = (-log(del_bin_width*(101-k))/alpha_del)**(1./gamma_del)
   expn_bins(k) = (1. - frac_fav_mutn)*(x0 - x1)
   x0 = x1
end do

z0 = 1.
do k=100,1,-1

   x0 = fav_bin_width*(101-k)/max_fav_fitness_gain
   x0 = min(1., x0)
   z1 = (-log(x0)/alpha_fav)**(1./gamma_fav)
   expn_bins(201-k) = frac_fav_mutn*(z0 - z1)
   if(x0 > fav_scale*lb_modulo) expn_bins(201-k) = 0.
   z0 = z1
end do

if(expn_bins(100) > 0.) &
   expn_bins(  1:100) = expn_bins(  1:100)/expn_bins(100)
if(expn_bins(101) > 0.) &
   expn_bins(101:200) = expn_bins(101:200)/expn_bins(101)

! Compute statistics on favorable and recessive mutations.

haplotype_bins = 0
bins_mutns     = 0

do i=1,current_pop_size

   do lb=1,num_linkage_subunits

      y0 = (linkage_block_fitness(lb,1,i) - 1.d0)/haplotype_bin_width

      if(y0 < 1.d-20) then
         k = 100 + int(y0)
         if(k > 0)    haplotype_bins(k) = haplotype_bins(k) + 1.
      else
         k = 101 + int(y0)
         if(k <= 200) haplotype_bins(k) = haplotype_bins(k) + 1.
      end if

      y0 = (linkage_block_fitness(lb,2,i) - 1.d0)/haplotype_bin_width

      if(y0 < 1.d-20) then
         k = 100 + int(y0)
         if(k > 0)    haplotype_bins(k) = haplotype_bins(k) + 1.
      else
         k = 101 + int(y0)
         if(k <= 200) haplotype_bins(k) = haplotype_bins(k) + 1.
      end if

   end do

   do j=2,dmutn(1,1,i)+1
      x = mod(abs(dmutn(j,1,i)),lb_modulo)*del_scale
      d = dexp(-alpha_del*x**gamma_del)
      k = 100 - int(d/del_bin_width)
      if(dmutn(j,1,i) < 0) then
         if(k > 0 .and. k <= 100) then
            bins_mutns(k,1) = bins_mutns(k,1) + 1.
         end if
      else
         if(k > 0 .and. k <= 100) then
            bins_mutns(k,2) = bins_mutns(k,2) + 1.
         end if
      end if
   end do

   do j=2,dmutn(1,2,i)+1
      x = mod(abs(dmutn(j,2,i)),lb_modulo)*del_scale
      d = dexp(-alpha_del*x**gamma_del)
      k = 100 - int(d/del_bin_width)
      if(dmutn(j,2,i) < 0) then
         if(k > 0 .and. k <= 100) then
            bins_mutns(k,1) = bins_mutns(k,1) + 1.
         end if
      else
         if(k > 0 .and. k <= 100) then
            bins_mutns(k,2) = bins_mutns(k,2) + 1.
         end if
      end if
   end do

   do j=2,fmutn(1,1,i)+1
      d = dexp(-alpha_fav*(mod(abs(fmutn(j,1,i)),lb_modulo) &
                            *fav_scale)**gamma_fav)
      k = 101 + int(d*max_fav_fitness_gain/fav_bin_width)

      if(fmutn(j,1,i) < 0) then
         if(k > 100 .and. k <= 200) then
            bins_mutns(k,1) = bins_mutns(k,1) + 1.
         end if
      else
         if(k > 100 .and. k <= 200) then
            bins_mutns(k,2) = bins_mutns(k,2) + 1.
         end if
      end if
   end do

   do j=2,fmutn(1,2,i)+1
      d = dexp(-alpha_fav*(mod(abs(fmutn(j,2,i)),lb_modulo) &
                            *fav_scale)**gamma_fav)
      k = 101 + int(d*max_fav_fitness_gain/fav_bin_width)

      if(fmutn(j,2,i) < 0) then
         if(k > 100 .and. k <= 200) then
            bins_mutns(k,1) = bins_mutns(k,1) + 1.
         end if
      else
         if(k > 100 .and. k <= 200) then
            bins_mutns(k,2) = bins_mutns(k,2) + 1.
         end if
      end if
   end do

end do

x0 = 1.e-10
do k=1,200
   x0 = max(x0, haplotype_bins(k))
end do

haplotype_bins = haplotype_bins/x0

if(tracking_threshold /= 1.0) then
   bins_mutns(  1:100,1) = bins_mutns(  1:100,1) &
              *expn_bins( 99)/(bins_mutns( 99,1) + 1.e-10)
   bins_mutns(  1:100,2) = bins_mutns(  1:100,2) &
              *expn_bins( 99)/(bins_mutns( 99,2) + 1.e-10)
   bins_mutns(101:200,1) = bins_mutns(101:200,1) &
              *expn_bins(102)/(bins_mutns(102,1) + 1.e-10)
   bins_mutns(101:200,2) = bins_mutns(101:200,2) &
              *expn_bins(102)/(bins_mutns(102,2) + 1.e-10)
end if

rewind (4)
write(4,'("# generation = ",i8)') gen
write(4,'("# effect-bin",2x,"theory(red)",1x,"lb-fitns(g)", &
          2x,"dominants",2x,"recessives")')

write(4, '(5e12.4)') ((k - 100.5)*haplotype_bin_width, &
        expn_bins(k), haplotype_bins(k),  &
        bins_mutns(k,2), bins_mutns(k,1), k=1,200)

avg_lb_effect = (post_sel_fitness - 1.d0)/(2*num_linkage_subunits)

num_del_lb = 0
num_fav_lb = 0

do k=1,100
   num_del_lb = num_del_lb + haplotype_bins(k)
end do

do k = 101,200
   num_fav_lb = num_fav_lb + haplotype_bins(k)
end do
  
lb_fitness_frac_positive = num_fav_lb/(num_del_lb + num_fav_lb)

write(4,'("# favorable x-axis scaling = ",1pe12.4)')  &
             fav_bin_width/del_bin_width 
write(4,'("# avg_linkage_block_effect = ",1pe12.4)')  &
             avg_lb_effect
write(4,'("# lb_fitness_percent_positive = ",f12.4)') &
             lb_fitness_frac_positive*100.
call flush(4)

 
end subroutine diagnostics_near_neutrals_plot

subroutine diagnostics_contrasting_alleles(dmutn, nmutn, fmutn, &
   cum_effect, initial_allele_effects, ica_count, max_size, list)

! This routine analyzes the distribution of the initial contrasting
! alleles and their effects on overall fitness.  When logical
! variable list is .true., routine outputs a list of allele 
! frequencies.
use inputs
include 'common.h'

integer max_size, count(num_linkage_subunits,3)
! note: we were sending in lb_offsprng_mutn_count as the count buffer
! to be able to reuse already existing memory.  However, now
! this routines have moved to the mating subroutine.  In the future,
! we should pass in an unused array to reuse the space.
integer dmutn(max_del_mutn_per_indiv/2,2,max_size)
integer nmutn(max_neu_mutn_per_indiv/2,2,max_size)
integer fmutn(max_fav_mutn_per_indiv/2,2,max_size)
real initial_allele_effects(num_linkage_subunits)
real w, effect, freq
real*8  cum_effect(pop_size)
integer i, lb, m, indx, ica_count(3)
integer zygous(num_linkage_subunits)
character dom*3
logical list

w = multiplicative_weighting

count = 0
cum_effect = 1.d0

do i=1,current_pop_size

   zygous = 0

   do m=2,dmutn(1,1,i)+1
      if(mod(dmutn(m,1,i), lb_modulo) == lb_modulo - 1) then 
         lb = dmutn(m,1,i)/lb_modulo + 1
         zygous(lb) = zygous(lb) + 1
      end if
   end do

   do m=2,dmutn(1,2,i)+1
      if(mod(dmutn(m,2,i), lb_modulo) == lb_modulo - 1) then 
         lb = dmutn(m,2,i)/lb_modulo + 1
         zygous(lb) = zygous(lb) + 1
      end if
   end do

   do lb=1,num_linkage_subunits
      if(zygous(lb) == 1) then
         effect = initial_allele_effects(lb)*recessive_hetero_expression
         count(lb,1) = count(lb,1) + 1
         cum_effect(i) = (cum_effect(i) - (1.d0 - w)*effect) &
                                        * (1.d0 - w *effect)
      elseif(zygous(lb) == 2) then
         effect = initial_allele_effects(lb)
         count(lb,1) = count(lb,1) + 2
         cum_effect(i) = (cum_effect(i) - (1.d0 - w)*effect) &
                                        * (1.d0 - w *effect)
      end if
   end do

   zygous = 0

   do m=2,fmutn(1,1,i)+1
      if(mod(fmutn(m,1,i), lb_modulo) == lb_modulo - 1) then 
         lb = fmutn(m,1,i)/lb_modulo + 1
         zygous(lb) = zygous(lb) + 1
      end if
   end do

   do m=2,fmutn(1,2,i)+1
      if(mod(fmutn(m,2,i), lb_modulo) == lb_modulo - 1) then 
         lb = fmutn(m,2,i)/lb_modulo + 1
         zygous(lb) = zygous(lb) + 1
      end if
   end do

   do lb=1,num_linkage_subunits
      if(zygous(lb) == 1) then
         effect = initial_allele_effects(lb)*dominant_hetero_expression
         count(lb,2) = count(lb,2) + 1
         cum_effect(i) = (cum_effect(i) + (1.d0 - w)*effect) &
                                        * (1.d0 + w *effect)
      elseif(zygous(lb) == 2) then
         effect = initial_allele_effects(lb)
         count(lb,2) = count(lb,2) + 2
         cum_effect(i) = (cum_effect(i) + (1.d0 - w)*effect) &
                                        * (1.d0 + w *effect)
      end if
   end do

end do

ica_count = 0

do lb=1,num_linkage_subunits
   ica_count(1) = ica_count(1) + count(lb,1)
   ica_count(2) = ica_count(2) + count(lb,2)
   ica_count(3) = ica_count(3) + count(lb,3)
end do

if(.not.list) then

!  Compute average effect of initial contrasting alleles.

   ica_mean_effect = 0.

   do i=1,current_pop_size
      ica_mean_effect = ica_mean_effect + (cum_effect(i) - 1.d0)
   end do

   ica_mean_effect = ica_mean_effect/current_pop_size 

!  Compute mean frequency of positive initial contrasting alleles
!  and the number of positive initial contrasting alleles fixed and lost.

   fav_lost  = 0
   fav_fixed = 0
   fav_mean_freq = 0

   do lb=1,num_linkage_subunits
      fav_mean_freq = fav_mean_freq + count(lb,2)
      if(abs(initial_allele_effects(lb)) > 0.) then
         if(count(lb,2) == 2*current_pop_size) fav_fixed = fav_fixed + 1
         if(count(lb,2) == 0) fav_lost = fav_lost + 1
      end if
   end do

   fav_mean_freq = fav_mean_freq &
                   /(2*current_pop_size*num_contrasting_alleles)

else

   if(.not. is_parallel) then
      write(6,'(/3x,"List of initial contrasting allele freqencies"/ &
            "     and fitness effect values at end of run:"// &
            6x,"allele   linkage  favorable  homozygous"/ &
            15x,"subunit  frequency    effect")')
   end if

   write(9,'(/3x,"List of initial contrasting allele freqencies"/ &
         "     and fitness effect values at end of run:"// &
         6x,"allele   linkage  favorable  homozygous"/ &
         15x,"subunit  frequency    effect")')

   indx = 0
   do lb=1,num_linkage_subunits
      if(abs(initial_allele_effects(lb)) > 0.) then
         indx = indx + 1
         freq = 0.5*real(count(lb,2))/real(current_pop_size)
         if(.not. is_parallel) &
            write(6,'(i10,i10,f12.4,f11.4,6x,a3)') indx, lb, freq, &
                  abs(initial_allele_effects(lb))
         write(9,'(i10,i10,f12.4,f11.4,6x,a3)') indx, lb, freq, &
               abs(initial_allele_effects(lb))
      end if
   end do

   if(.not. is_parallel) write(6,*)
   write(9,*)

end if

end subroutine diagnostics_contrasting_alleles

subroutine diagnostics_heterozygosity(dmutn, fmutn)
use inputs
include 'common.h'
integer dmutn(max_del_mutn_per_indiv/2,2,*)
integer fmutn(max_fav_mutn_per_indiv/2,2,*)
integer i,j,k
integer count_homozygous, count_heterozygous, total_analyzed
real mutn_per_individual, fraction_heterozygous
real individual_heterozygosity 

! Analyze percent homozygosity and heterozygosity of deleterious mutations

count_homozygous = 0
count_heterozygous = 0
do i=1, current_pop_size

   j = 2
   do k = 2, dmutn(1,j,i)+1

      do while(abs(dmutn(k,1,i))>abs(dmutn(j,2,i)) .and. j <= dmutn(1,2,i))
         j = j + 1
      end do

      if(dmutn(k,1,i) == dmutn(j,2,i)) then
         count_homozygous = count_homozygous + 1
      else
         count_heterozygous = count_heterozygous + 1
      end if

   end do 
end do

total_analyzed = count_heterozygous+count_homozygous
if(total_analyzed < 10) goto 10
fraction_heterozygous = count_heterozygous/real(total_analyzed)
mutn_per_individual = tracked_del_mutn/real(current_pop_size)
write(*,'(/"HETEROZYGOSITY ANALYSIS OF DELETERIOUS MUTATIONS:")')
write(*,'("Out of    ",i10," total deleterious mutations")') &
  total_analyzed
write(*,'("there are ",i10," homozygous mutations")')   &
  count_homozygous
write(*,'("and       ",i10," heterozygous mutations")') &
  count_heterozygous
write(*,'("resulting in a percent heterozygosity of: ",f5.2,"%")') &
  fraction_heterozygous*100
write(*,'("There are: ",f7.1," tracked deleterious " &
  "mutations per individual")') mutn_per_individual
write(*,*)

10   continue

! Analyze percent homozygosity and heterozygosity of favorable mutations
    
if(frac_fav_mutn > 0. .or. polygenic_beneficials) then

count_homozygous = 0
count_heterozygous = 0
do i=1, current_pop_size

   j = 2
   do k=2, fmutn(1,1,i)+1

      do while(abs(fmutn(k,1,i))>abs(fmutn(j,2,i)) .and. j<=fmutn(1,2,i))
         j = j + 1
      end do

      if(fmutn(k,1,i) == fmutn(j,2,i)) then
         count_homozygous = count_homozygous + 1
      else
         count_heterozygous = count_heterozygous + 1
      end if
      
   end do 
end do

total_analyzed = count_heterozygous+count_homozygous
if(total_analyzed < 10) goto 20

fraction_heterozygous = count_heterozygous/real(total_analyzed)
mutn_per_individual = tracked_fav_mutn/real(current_pop_size)
write(*,'(/"HETEROZYGOSITY ANALYSIS OF FAVORABLE MUTATIONS:")')
write(*,'("Out of    ",i10," total favorable mutations")') total_analyzed
write(*,'("there are ",i10," homozygous mutations")')      count_homozygous
write(*,'("and       ",i10," heterozygous mutations")')    count_heterozygous
write(*,'("resulting in a percent heterozygosity of: ",f5.2,"%")') &
  fraction_heterozygous*100
write(*,'("There are: ",f7.1, " tracked favorable mutation per individual")') & 
  mutn_per_individual
write(*,*)
   
end if

20   continue

end subroutine diagnostics_heterozygosity

subroutine diagnostics_polymorphisms_plot(dmutn, nmutn, fmutn, &
                                          mfirst, max_size, gen)
use polygenic
use inputs
use init
include 'common.h'
! MNP = max number of polymorphisms, NP = number of bins
integer, parameter :: MNP=100000, NB=100 
integer dmutn(max_del_mutn_per_indiv/2,2,*)
integer nmutn(max_neu_mutn_per_indiv/2,2,*)
integer fmutn(max_fav_mutn_per_indiv/2,2,*)
integer mutn_count(MNP), mutn_list(MNP)
integer mfirst(2,max_size)
integer i, j, k, m, n, lb, it, it0, ie, max_size, gen, list_count
integer dcount, fcount, par_dcount, par_fcount
integer global_mutn_count(MNP,num_tribes)
integer global_mutn_list(MNP,num_tribes)
integer global_list_count(num_tribes)
integer num_falleles(3), num_dalleles(3), lb_limit, dwarn, fwarn
integer par_num_falleles(3), par_num_dalleles(3), mutn_limit
integer ncount, par_ncount, num_nalleles(3), par_nalleles(3)
integer nwarn, mutn_thres, mutn, jmax
real*8 dpbin(NB), dpbin_count(NB), dpbin_max, pbin_width
real*8 npbin(NB), npbin_count(NB), npbin_max, x
real*8 fpbin(NB), fpbin_count(NB), fpbin_max, scale_factor
real*8 par_dpbin(NB), par_fpbin(NB), dsum, fsum
real*8 par_npbin(NB), par_npbin_count(NB), nsum
real*8 par_dpbin_count(NB),par_fpbin_count(NB)
real*8 febin(10,NB), fe_bin_width
real   bin_fitness(11), bin_center(10)

! For diploid organisms, mutations can be homozygous, and so the
! limiting mutation count of each member of the population is
! homozygous in a given mutation is twice the population size.
! However, for clonal haploid organisms or pure self-fertilization,
! the limiting mutation count when each member carries a given
! mutation is the population size.  The polymorphism bin width
! pbin_width value reflects this difference.

if (is_parallel) then 
   pbin_width = 2.*current_pop_size/NB*num_tribes
else
   pbin_width = 2.*current_pop_size/NB
end if

if(recombination_model == clonal) pbin_width = pbin_width/2.

fe_bin_width = -log(tracking_threshold)/10.

if(int(pbin_width) == 0) then
   write(*,*) 'Polymorphism analysis skipped because population', &
              ' size is too small.'
   return
end if

! Compute statistics on deleterious polymorphisms.  

! To keep the required computer time reasonable, limit the maximum 
! number of polymorphisms to be considered to be MNP.  This may 
! result in the analysis being performed only over a portion of 
! the total total population.

dsum   = 0.
dpbin  = 0.
febin  = 0.
dwarn  = 0

if(frac_fav_mutn /= 1.) then

   mfirst = 2
   
   ! Loop over recessives first and dominants second.
   
   it0 = 1
   if(fraction_recessive == 0.) it0 = 2
   
   do it=it0,2
      
      do lb=1,num_linkage_subunits

         if(it == 1) then
            mutn_limit = -lb_modulo*(num_linkage_subunits - lb)
         else
            mutn_limit =  lb_modulo*lb
         end if
         
         call count_alleles(dmutn,max_del_mutn_per_indiv,MNP,mutn_limit, &
                            mutn_list,mutn_count,list_count,mfirst,gen)

         call bin_alleles(MNP,NB,list_count,mutn_count,pbin_width, &
                          dpbin,dsum,dwarn)

         ! The following do statement is for the .pmd file, 
         ! which bins mutations based on their fitness effect
         do k=1,list_count
            j = min(NB, 1 + int(mutn_count(k)/pbin_width))
            x = mod(abs(mutn_list(k)),lb_modulo)*del_scale
            if(x < 1.d0) then
               ie = 1 + int(alpha_del*x**gamma_del/fe_bin_width)
               febin(ie,j) = febin(ie,j) + 1
            end if
         end do
         
      end do

   end do

end if

! Compute statistics on neutral polymorphisms.  

nsum   = 0.
npbin  = 0.
nwarn  = 0

if(track_neutrals) then
   
   if(polygenic_beneficials) then

      call count_string_alleles(nmutn,max_neu_mutn_per_indiv,MNP,  &
                                mutn_list,mutn_count,list_count)

      call bin_alleles(MNP,NB,list_count,mutn_count,pbin_width, &
              npbin,nsum,nwarn)

   else

      mfirst = 2
   
      do lb=1,num_linkage_subunits
      
         mutn_limit = lb_modulo*lb

         call count_alleles(nmutn,max_neu_mutn_per_indiv,MNP,mutn_limit, &
                 mutn_list,mutn_count,list_count,mfirst,gen)

         call bin_alleles(MNP,NB,list_count,mutn_count,pbin_width, &
                       npbin,nsum,nwarn)

      end do

   end if

end if

! Compute statistics on favorable polymorphisms.  

fsum   = 0.
fpbin  = 0.
fwarn  = 0

if(frac_fav_mutn > 0. .or. polygenic_beneficials) then

   mfirst = 2

   ! Loop over recessives first and dominants second.

   it0 = 1
   if(fraction_recessive == 0.) it0 = 2
   
   do it=it0,2
      
      do lb=1,num_linkage_subunits

         if(it == 1) then
            mutn_limit = -lb_modulo*(num_linkage_subunits - lb)
         else
            mutn_limit =  lb_modulo*lb
         end if
         
         call count_alleles(fmutn,max_fav_mutn_per_indiv,MNP,mutn_limit, &
                            mutn_list,mutn_count,list_count,mfirst,gen)

         call bin_alleles(MNP,NB,list_count,mutn_count,pbin_width, &
                          fpbin,fsum,fwarn)
         
      end do
      
   end do
   
   if(polygenic_beneficials) then
      do i=1,current_pop_size
         jmax = 2
         if(recombination_model==clonal) jmax = 1
         do j=1,jmax
            do k=1,num_polys_cumulative
               if(fmutn(2,j,i) == pmutn(k)%id) then
                  pmutn(k)%gen_exit = gen + 1
                  poly_fixed_fmutn = pmutn(k)%id
               end if
            end do
         end do
      end do
   end if

end if

if(polygenic_beneficials) then

! Remove favorable allele count from the neutral allele count in the polygenic
! case, as they are redundant.

   npbin = npbin - fpbin
   nsum  = nsum  - fsum

end if

dpbin_max = 1
npbin_max = 1
fpbin_max = 1
num_dalleles(2) = 0
num_nalleles(2) = 0
num_falleles(2) = 0

do j=1,NB
   dpbin_max = max(dpbin_max, dpbin(j))
   npbin_max = max(npbin_max, npbin(j))
   fpbin_max = max(fpbin_max, fpbin(j))
   if(j > 1 .and. j < NB) then
      num_dalleles(2) = num_dalleles(2) + dpbin(j)
      num_nalleles(2) = num_nalleles(2) + npbin(j)
      num_falleles(2) = num_falleles(2) + fpbin(j)
   end if
end do

num_dalleles(1) = dpbin(1)
num_nalleles(1) = npbin(1)
num_falleles(1) = fpbin(1)
num_dalleles(3) = dpbin(NB)
num_nalleles(3) = npbin(NB)
num_falleles(3) = fpbin(NB)

dcount = dpbin(1) + num_dalleles(2) + dpbin(NB)
ncount = npbin(1) + num_nalleles(2) + npbin(NB)
fcount = fpbin(1) + num_falleles(2) + fpbin(NB)

dpbin_count = dpbin
dpbin = dpbin/dpbin_max

npbin_count = npbin
npbin = npbin/npbin_max

fpbin_count = fpbin
fpbin = fpbin/fpbin_max

! Compute values for fitness effect bin centers.

do k=1,11
   bin_fitness(k) = exp(-fe_bin_width*(k - 1))
   if (k > 1) then
      bin_center(k-1) = 0.5*(bin_fitness(k) + bin_fitness(k-1))
   end if
end do

if (verbosity==2) then
   write(19,'("# generation = ",i8)') gen
   write(19,"('#',9x,'Table of polymorphism frequency vs. fitness' &
      ' effect category'/'#'/ &
      '#freq',17x,'fitness effect category center value')") 
   write(19,'(3x,1p4e7.0,6e8.1)') bin_center(1:10)
   write(19,'(i3,4i7,6i8)') (k,int(febin(1:10,k)),k=1,NB)
   call flush(19)
endif

! Output global statistics to file caseid.000.plm
if(is_parallel) then

   if(myid.eq.0 .and. mod(gen,diagnostic_gens)==0 ) then
      rewind(21)
      write(21,'("# Number of tribes = ",i4)') num_tribes
      write(21,'("# generation = ",i8)') gen
      write(21,'("# frequency del_normalized fav_normalized", &
           "  del_count fav_count neu_normalized  neu_count")')
      write(21,'(i11,2f15.11,2f11.0,f15.11,f11.0)') (k, dpbin(k), &
          fpbin(k), dpbin_count(k), fpbin_count(k), npbin(k),     &
                    npbin_count(k), k=1,NB)
      write(21,'("# Allele summary statistics:")')
      write(21,'("#  Very rare",6x,"Polymorphic",9x,"Fixed", &
                 11x,"Total")')
      write(21,'("#",3x,"(0-1%)",9x,"(1-99%)",11x,"(100%)")')
      write(21,'("#",4i10," deleterious")') &
          num_dalleles(1), num_dalleles(2), num_dalleles(3), dcount
      write(21,'("#",4i10," favorable")')   &
          num_falleles(1), num_falleles(2), num_falleles(3), fcount
      write(21,'("#",4i10," neutral")')     &
          num_nalleles(1), num_nalleles(2), num_nalleles(3), ncount
      call flush(21)
   end if 

else 

   rewind(11)
   if (mod(gen,plot_allele_gens)==0.and.verbosity>0) then
       write(11,'("# generation = ",i8)') gen
       write(11,'("# frequency del_normalized fav_normalized ", &
                  "  del_count fav_count neu_normalized  neu_count")')
       write(11,'(i11,2f15.11,2f11.0,f15.11,f11.0)')  (k, dpbin(k),  &
                 fpbin(k), dpbin_count(k), fpbin_count(k), npbin(k), &
                           npbin_count(k), k=1,NB)
   endif

end if

! Keep a "snapshot" file of latest polymorphism data.

if (verbosity > 0) then
   rewind(13)
   write(13,'("# generation = ",i8)') gen
   write(13,'("# frequency del_normalized fav_normalized ", &
              "  del_count fav_count neu_normalized  neu_count")')
   write(13,'(i11,2f15.11,2f11.0,f15.11,f11.0)')  (k, dpbin(k),  &
             fpbin(k), dpbin_count(k), fpbin_count(k), npbin(k), &
                       npbin_count(k), k=1,NB)
   call flush(13)
end if

if(.not.is_parallel .and. mod(gen,plot_allele_gens)==0.and.verbosity>0) then
   write(11,"('#',11x,'Allele summary statistics (tracked mutations only):'/ &
      '#   (Statistics are based on ',i12,' tracked deleterious mutations'/  &
      '#                            ',i12,' tracked   favorable mutations'/  &
      '#                        and ',i12,' tracked     neutral mutations.)'/ &
      '#    Very rare   Polymorphic     Fixed      Total'/   &
      '#      (0-1%)      (1-99%)      (100%)')") int(dsum), &
      int(fsum), int(nsum)
   write(11,"('#',4i12,' deleterious')") num_dalleles(1), &
              num_dalleles(2), num_dalleles(3), dcount
   write(11,"('#',4i12,' favorable')")   num_falleles(1), &
              num_falleles(2), num_falleles(3), fcount
   write(11,"('#',4i12,' neutral')")   num_nalleles(1),   &
              num_nalleles(2), num_nalleles(3), ncount
   if(dwarn == 1) write(11,'("# Warning: Number of deleterious " &
      "polymorhisms exceeded the linkage block limit of ",i8)') MNP
   if(fwarn == 1) write(11,'("# Warning: Number of   favorable " &
      "polymorhisms exceeded the linkage block limit of ",i8)') MNP
   if(nwarn == 1) write(11,'("# Warning: Number of     neutral " &
      "polymorhisms exceeded the linkage block limit of ",i8)') MNP
   call flush(11)
end if

if(myid == 0 .and. mod(gen,diagnostic_gens)==0 ) then
   write(6,"(/12x,'Allele summary statistics (tracked mutations only):'/ &
    '    (Statistics are based on ',i12,' tracked deleterious mutations'/ &
    '                             ',i12,' tracked   favorable mutations'/ &
    '                         and ',i12,' tracked     neutral mutations.)'/ &
    '     Very rare   Polymorphic     Fixed      Total'/   &
    '       (0-1%)      (1-99%)      (100%)')")            &
          int(dsum), int(fsum), int(nsum)
   write(6,"(' ',4i12,' deleterious')") num_dalleles(1), &
        num_dalleles(2), num_dalleles(3), dcount
   write(6,"(' ',4i12,' favorable')")   num_falleles(1), &
        num_falleles(2), num_falleles(3), fcount
   write(6,"(' ',4i12,' neutral')")   num_nalleles(1),   &
        num_nalleles(2), num_nalleles(3), ncount
   if(dwarn == 1) write(6,'("  Warning: Number of deleterious " &
        "polymorhisms exceeded the linkage block limit of ",i8)') MNP
   if(fwarn == 1) write(6,'("  Warning: Number of   favorable " &
        "polymorhisms exceeded the linkage block limit of ",i8)') MNP
   if(nwarn == 1) write(6,'("  Warning: Number of     neutral " &
        "polymorhisms exceeded the linkage block limit of ",i8)') MNP
end if

if(mod(gen,diagnostic_gens)==0) then
   write(9,"(/12x,'Allele summary statistics (tracked mutations only):'/ &
    '    (Statistics are based on ',i12,' tracked deleterious mutations'/ &
    '                             ',i12,' tracked   favorable mutations'/ &
    '                         and ',i12,' tracked     neutral mutations.)'/ &
    '     Very rare   Polymorphic     Fixed      Total'/   &
    '       (0-1%)      (1-99%)      (100%)')")            &
          int(dsum), int(fsum), int(nsum)
   write(9,"(' ',4i12,' deleterious')") num_dalleles(1), &
        num_dalleles(2), num_dalleles(3), dcount
   write(9,"(' ',4i12,' favorable')")   num_falleles(1), &
        num_falleles(2), num_falleles(3), fcount
   write(9,"(' ',4i12,' neutral')")   num_nalleles(1),   &
        num_nalleles(2), num_nalleles(3), ncount
   if(dwarn == 1) write(9,'("  Warning: Number of deleterious " &
        "polymorhisms exceeded the linkage block limit of ",i8)') MNP
   if(fwarn == 1) write(9,'("  Warning: Number of   favorable " &
        "polymorhisms exceeded the linkage block limit of ",i8)') MNP
   if(nwarn == 1) write(9,'("  Warning: Number of     neutral " &
        "polymorhisms exceeded the linkage block limit of ",i8)') MNP
endif

end subroutine diagnostics_polymorphisms_plot

subroutine diagnostics_selection(fitness_pre_sel,fitness_post_sel, &
                                 total_offspring,gen)
use selection_module
use inputs
include 'common.h'
integer, parameter :: NB=200 ! Number of bins
real*8 fitness_pre_sel(*), fitness_post_sel(*)
integer total_offspring, gen, i, j, jmax, max_bin
integer sel_bins(NB,2), par_sel_bins(NB,2)
integer fid, oneortwo
real select_ratio, srm, srp
real*8 par_pre_sel_fitness, par_post_sel_fitness, &
       par_pre_sel_geno_sd, par_pre_sel_pheno_sd, &
       par_pre_sel_corr, par_post_sel_geno_sd,    &
       par_post_sel_pheno_sd, par_post_sel_corr

sel_bins = 0

do i=1,total_offspring
   j = min(NB, int(100.*fitness_pre_sel(i)))
   sel_bins(j,1) = sel_bins(j,1) + 1
end do

do i=1,current_pop_size
   j = min(NB, int(100.*fitness_post_sel(i)))
   sel_bins(j,2) = sel_bins(j,2) + 1
end do
 
jmax    = 1
max_bin = 0

do j=1,NB
   if(sel_bins(j,2) > max_bin) then
      max_bin = sel_bins(j,2)
      jmax = j
   end if
end do


! If parallel is turned on, write two files
! caseid.001.sel which contains results for one tribe
! and caseid.000.sel which contains data averaged from all tribes.
if (is_parallel .and. myid==0) then
    oneortwo = 2
else 
    oneortwo = 1
end if 

fid = 24

do i = 1, oneortwo
   if(oneortwo==2) then
      fid = 34
      sel_bins = par_sel_bins
      pre_sel_fitness  = par_pre_sel_fitness
      post_sel_fitness = par_post_sel_fitness
      pre_sel_geno_sd  = par_pre_sel_geno_sd
      post_sel_geno_sd = par_post_sel_geno_sd 
      pre_sel_pheno_sd = par_pre_sel_pheno_sd
   end if

   rewind (fid)

   write(fid,'("# generation = ",i8)') gen
   write(fid,'("# heritability            =",f9.5)') heritability 
   write(fid,'("# non scaling noise       =",f9.5)') non_scaling_noise
   write(fid,'("# pre  selection fitness  =",f9.5)') pre_sel_fitness
   write(fid,'("# post selection fitness  =",f9.5)') post_sel_fitness
   write(fid,'("# pre  selection geno  sd =",f9.5)') pre_sel_geno_sd
   write(fid,'("# post selection geno  sd =",f9.5)') post_sel_geno_sd
   write(fid,'("# pre  selection pheno sd =",f9.5)') pre_sel_pheno_sd
   write(fid,'("# post selection pheno sd =",f9.5)') post_sel_pheno_sd
   write(fid,'("# Effect of Selection on Phenotypic Fitness Distribution")')
   write(fid,'("# Note: Individuals with phenotypic fitness > ", &
              "1.5 are all put into last bin")')
   write(fid,'("# max bin fitness",5x,"before selection",5x, &
            "after selection",5x,"ratio")')
   do j=1,NB-1
      srm = sel_bins(j  ,2)/(real(sel_bins(j  ,1)) + 0.000001)
      if(sel_bins(j  ,1) == 0) srm = -2.
      srp = sel_bins(j+1,2)/(real(sel_bins(j+1,1)) + 0.000001)
      if(sel_bins(j+1,1) == 0) srp = -2.
      select_ratio = 0.5*(srm + srp)
      if(sel_bins(j,2) + sel_bins(j+1,2) <= max(4, max_bin/50)) &
         select_ratio = -1.
      if(select_ratio > 0.) then
         write(fid,'(f10.3,2i20,f18.3)') j*0.01, sel_bins(j,:), select_ratio
      else
         write(fid,'(f10.3,2i20,15x,"?")') j*0.01, sel_bins(j,:) 
      end if
   end do    

   call flush(fid)

end do

end subroutine diagnostics_selection

subroutine count_alleles(xmutn,max_mutn_per_indiv,MNP,mutn_limit, &
                         mutn_list,mutn_count,list_count,mfirst,gen)
use inputs
use polygenic
include 'common.h'
include 'mpif.h'
integer, intent(in)    :: MNP
integer, intent(in)    :: max_mutn_per_indiv
integer, intent(in)    :: xmutn(max_mutn_per_indiv/2,2,*)
integer, intent(in)    :: mutn_limit
integer, intent(in)    :: gen
integer, intent(out)   :: mutn_list(MNP)
integer, intent(out)   :: mutn_count(MNP)
integer, intent(out)   :: list_count
integer, intent(inout) :: mfirst(2,*)

integer :: i, j, k, m
integer :: global_mutn_count(MNP,num_tribes)
integer :: global_mutn_list(MNP,num_tribes)
integer :: global_list_count(num_tribes)
logical :: new_mutn
integer :: mutn, jmax

! mutn_list is the list of mutation indices being analyzed for their frequency.
! mutn_count is the number of occurrences of each mutation in mutn_list.
! list_count is the total number of alleles being analyzed.

mutn_list  = 0
mutn_count = 0
list_count = 0

do i=1,current_pop_size
   jmax = 2
   if(polygenic_beneficials .and. recombination_model==clonal) jmax = 1
   do j=1,jmax
      m = mfirst(j,i)
      do while(xmutn(m,j,i) < mutn_limit .and. &
           m <=  xmutn(1,j,i)+1 .and. list_count < MNP)
         mutn = mod(abs(xmutn(m,j,i)), lb_modulo)
         if(mutn /= lb_modulo-1) then
            if(upload_mutations) then
               do k=1,num_uploaded_mutn
                  if(xmutn(m,j,i) == uploaded_mutn(k)) then
                     list_count = min(MNP, list_count + 1)
                     mutn_count(k) = mutn_count(k) + 1
                     mutn_list (list_count) = xmutn(m,j,i)
                  end if
               end do
            else
               new_mutn = .true.
               do k=1,list_count
                  if(xmutn(m,j,i) == mutn_list(k)) then
                     mutn_count(k) = mutn_count(k) + 1
                     new_mutn = .false.
                  end if
               end do
               if(new_mutn .and. .not.(polygenic_beneficials  &
                           .and. xmutn(m,j,i) == 0)) then
                  list_count = min(MNP, list_count + 1)
                  mutn_list (list_count) = xmutn(m,j,i)
                  mutn_count(list_count) = 1
               end if
            end if
         end if
         m = m + 1
      end do
      mfirst(j,i) = m
   end do
end do


end subroutine count_alleles

subroutine count_string_alleles(xmutn,max_mutn_per_indiv,MNP, &
                                mutn_list,mutn_count,list_count)
use inputs
use polygenic
include 'common.h'
include 'mpif.h'
integer, intent(in)    :: MNP
integer, intent(in)    :: max_mutn_per_indiv
integer, intent(in)    :: xmutn(max_mutn_per_indiv/2,2,*)
integer, intent(out)   :: list_count
integer, intent(out)   :: mutn_list(MNP)
integer, intent(out)   :: mutn_count(MNP)

integer :: i, j, k
logical :: new_mutn
integer :: mutn, jmax

! mutn_list is the list of mutation indices being analyzed for their frequency.
! mutn_count is the number of occurrences of each mutation in mutn_list.
! list_count is the total number of alleles being analyzed.

mutn_list  = 0
mutn_count = 0
list_count = 0

do i=1,current_pop_size
   jmax = 2
   if(recombination_model==clonal) jmax = 1
   do j=1,jmax
      mutn = 0
      do k=1,num_linkage_subunits
         mutn = mutn + (xmutn(1+k,j,i) - 1)*4**(k - 1)
      end do
      new_mutn = .true.
      do k=1,list_count
         if(mutn == mutn_list(k)) then
            mutn_count(k) = mutn_count(k) + 1
            new_mutn = .false.
         end if
      end do
      if(new_mutn) then
         list_count = min(MNP, list_count + 1)
         mutn_list (list_count) = mutn
         mutn_count(list_count) = 1
      end if
   end do
end do

end subroutine count_string_alleles

subroutine bin_alleles(MNP,NB,list_count,mutn_count,pbin_width,pbin,xsum,warn)
  implicit none

  integer, intent(in)    :: MNP, NB, list_count, mutn_count(MNP)
  real*8,  intent(in)    :: pbin_width
  real*8,  intent(inout) :: xsum, pbin(NB)
  integer, intent(inout) :: warn

  integer :: j,k

  ! Load the allele hit counts into the proper bins and count
  ! the total number of hits accumulated.
  
  do k=1,list_count
     j = min(NB, 1 + int(mutn_count(k)/pbin_width))
     pbin(j) = pbin(j) + 1
     xsum = xsum + mutn_count(k)
  end do
  
  if(list_count == MNP) warn = 1
  
end subroutine bin_alleles

