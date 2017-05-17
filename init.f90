module init

integer :: output_gens, diagnostic_gens, hst_gens
!integer:: lb_modulo

contains

subroutine initialize(myid_str)

use random_pkg
use polygenic
use profiler
use inputs
include 'common.h'

integer     :: values(8), npath, i, k
real        :: d1, d2, sum, del_mean, fav_mean, alpha, gamma
real        :: fraction_del_tracked
character*3 :: myid_str, zero_str
logical     :: file_exists

! Initialize the profile timer
sec = 0.0
cyclic_bottlenecking = .false.

! Output version information.  RCS will automatically update
! the following $Id string on check-in

write(6,*) 'VERSION >>> v2.6.3-9-ga744efb-dirty <<< VERSION'

call date_and_time(VALUES=values)

if(is_parallel) then
   !START_MPI
   call mpi_myinit(myid,ierr)

   write(myid_str,'(i3.3)') myid+1

   ! Open files containing run parameters for heterogeneous tribes
   if (.not.homogenous_tribes) then
      inquire(file='./mendel.in.'//myid_str, exist=file_exists)
      if ( file_exists ) then
         open (5, file='mendel.in.'//myid_str,status='old')
         call read_parameters(5)
         close(5)
      else
         print *, 'ERROR: could not find mendel.in.'//myid_str
         stop
      end if
   end if

   if (myid==0) write(*,*) 'subpopulation size is ',pop_size

   if(num_indiv_exchanged > pop_size) then
      write(6,*) 'ERROR: num_indiv_exchanged >= tribal pop_size'
      write(6,*) 'ERROR: decrease num_indiv_exchanged.'
      call mpi_myabort()
   end if
   !END_MPI
else
   write(myid_str,'(i3.3)') 0
   myid = 0
end if

npath = index(data_file_path,' ') - 1

open (9, file=data_file_path(1:npath)//case_id//'.'//myid_str &
         //'.out',status='unknown')
open(10, file=data_file_path(1:npath)//case_id//'.st8', &
        status='unknown')
write(10,*) 0
close(10)

if (restart_case) then
   open (7, file=data_file_path(1:npath)//case_id//'.'//myid_str &
      //'.hst',status='unknown',position='append')
else
   open (7, file=data_file_path(1:npath)//case_id//'.'//myid_str &
      //'.hst',status='unknown')
   if (tribal_competition) then
      write(7,'("# generation",4x,"fitness",9x,"fitness_sd",     &
                6x,"num_dmutns",4x,"num_fmutns",4x,"num_nmutns", &
                4x,"pop_size",4x,"global_fit",4x,"fert_fac")')
   else
      write(7,'("# generation",4x,"fitness",9x,"fitness_sd",     &
                6x,"num_dmutns",4x,"num_fmutns",4x,"num_nmutns", &
                4x,"pop_size")')
   end if
end if

if (verbosity > 0) then

   open (4, file=data_file_path(1:npath)//case_id//'.'//myid_str &
         //'.hap',status='unknown')

   open (8, file=data_file_path(1:npath)//case_id//'.'//myid_str &
         //'.dst',status='unknown')

   open(11, file=data_file_path(1:npath)//case_id//'.'//myid_str &
         //'.plm',status='unknown')

   if (.not.is_parallel) then
      open(12, file=data_file_path(1:npath)//case_id//'.'//myid_str &
           //'.plmcor',status='unknown')
   end if

   ! Maintain a separate file for windows which just has a snapshot
   ! of the most recent polymorphism information.
   open(13, file=data_file_path(1:npath)//case_id//'.'//myid_str &
            //'.pls',status='unknown')

   open(24, file=data_file_path(1:npath)//case_id//'.'//myid_str &
         //'.sel',status='unknown')

   open(25, file=data_file_path(1:npath)//case_id//'.'//myid_str &
         //'.thr',status='unknown')
   write(25,'("#  generation",20x,"selection thresholds"/"#",13x,&
     "del_dom_thres  del_rec_thres  fav_dom_thres  fav_rec_thres")')

   if (num_contrasting_alleles > 0) then
      open(30, file=data_file_path(1:npath)//case_id//'.'//myid_str &
           //'.ica',status='unknown')
   endif

endif

! Output rarely used ancillary files only when verbosity switch is set to 3
if (verbosity == 2) then
   ! open(15, file=data_file_path(1:npath)//case_id//'.'//myid_str &
   !          //'.ppm',status='unknown')
   ! Write header for PPM file
   !     write(15,'("P3")')
   !     write(15,'(2i)') pop_size, num_generations
   !     write(15,'("255")')
   open(16, file=data_file_path(1:npath)//case_id//'.'//myid_str &
         //'.fit',status='unknown')
   open(19, file=data_file_path(1:npath)//case_id//'.'//myid_str &
        //'.pmd',status='unknown')
   open(22, file=data_file_path(1:npath)//case_id//'.'//myid_str &
        //'.tim',status='unknown')
   write(22,'("# generation",2x,"time_per_gen(s)",2x, &
           "time_offspring(s)",2x,"time_selection(s)")')
   open(26, file=data_file_path(1:npath)//case_id//'.'//myid_str &
        //'.acc',status='unknown')
endif

! If parallel, write additional average files with name-tag 000.

if (is_parallel) then

   if (verbosity > 0) then

      write(zero_str,'(i3.3)') 0

      open (14, file=data_file_path(1:npath)//case_id//'.'//zero_str &
            //'.hap',status='unknown')

      open (17, file=data_file_path(1:npath)//case_id//'.'//zero_str &
            //'.hst',status='unknown')
      if (tribal_competition) then
         write(17,'("# generation",4x,"fitness",9x,"fitness_sd",     &
                   6x,"num_dmutns",4x,"num_fmutns",4x,"num_nmutns",  &
                   4x,"pop_size",2x,"global_fit")')
      else
         write(17,'("# generation",4x,"fitness",9x,"fitness_sd",     &
                   6x,"num_dmutns",4x,"num_fmutns",4x,"num_nmutns",  &
                   4x,"pop_size")')
      end if

      open (18, file=data_file_path(1:npath)//case_id//'.'//zero_str &
            //'.dst',status='unknown')

      open (21, file=data_file_path(1:npath)//case_id//'.'//zero_str &
            //'.plm',status='unknown')

      open (34, file=data_file_path(1:npath)//case_id//'.'//zero_str &
            //'.sel',status='unknown')

      open (35, file=data_file_path(1:npath)//case_id//'.'//zero_str &
            //'.thr',status='unknown')
      write(35,'("#  generation",20x,"selection thresholds"/"#",13x, &
      "del_dom_thres  del_rec_thres  fav_dom_thres  fav_rec_thres")')

   endif

   if(verbosity == 2) then
      open (23, file=data_file_path(1:npath)//case_id//'.'//zero_str &
            //'.tim',status='unknown')
      write(23,'("# generation",2x,"time_per_gen(s)",2x, &
                 "time_offspring(s)",2x,"time_selection(s)")')
   endif

end if

if(myid==0) &
   write(6, '(/" Case started ",i2,"/",i2,"/",i4," at ",i2,":", &
      i2,":",i2/)') values(2:3), values(1), values(5:7)

   write(9, '(/" Case started ",i2,"/",i2,"/",i4," at ",i2,":", &
      i2,":",i2/)') values(2:3), values(1), values(5:7)

! Check to ensure parameter fraction_random_death does not reduce
! actual fertility below 1.0.

if(reproductive_rate*(1. - fraction_random_death) < 1.) then
   write(6,'("ERROR: Input value of fraction_random_death ", &
             "implies a fertility of less than 1.0.")')
   write(9,'("ERROR: Input value of fraction_random_death ", &
             "implies a fertility of less than 1.0.")')
   stop
end if

! Limit the minimum value of heritability to be 10**-20.

heritability = max(1.e-20, heritability)

! Do not track neutrals if there are none to track.

if(fraction_neutral == 0.) track_neutrals = .false.

! If neutrals are to be tracked, set the tracking threshold to
! zero.

if(track_neutrals) tracking_threshold = 0.

! If back mutations are allowed, set the tracking threshold to zero
! so that all mutations are tracked.

if(allow_back_mutn) tracking_threshold = 0.

! Initialize number of cumulative back mutations to zero.

num_back_mutn = 0

if(.not.dynamic_linkage) haploid_chromosome_number = 1

! Initialize some quantities from the input values.

if(polygenic_beneficials) then
   fitness_distrib_type = 0
   polygenic_fixed = .false.
   if(len(polygenic_target) > 40) then
      write(*,*) 'ERROR: polygenic_target > 40. Need to increase bound.'
      stop
   endif
end if

if(clonal_haploid) recombination_model = clonal

! For clonal reproduction, set the number of chromosomes to one
! and the number of linkage blocks to one.

if(recombination_model == clonal) then
   haploid_chromosome_number  = 1
   num_linkage_subunits       = 1
   fraction_recessive         = 0.
   if(clonal_haploid) dominant_hetero_expression = 1.
end if

if(polygenic_beneficials) then
   poly_str_len = len_trim(polygenic_target)
   num_linkage_subunits = poly_str_len
   haploid_chromosome_number = 1
end if

if(recombination_model == suppressed) then
   dynamic_linkage = .false.
end if

! For the case of dynamic linkage, ensure that the number of linkage
! subunits is an integer times the haploid chromosome number.

if(dynamic_linkage) num_linkage_subunits = (num_linkage_subunits &
           /haploid_chromosome_number)*haploid_chromosome_number

! Echo the parameters for this case to the output file.

if(is_parallel) then
   if(myid == 0) call write_parameters(6)
else
   call write_parameters(6)
end if
call write_parameters(9)

! If the genome_size is zero, most likely means equal effect fitness
! mutation distribution is selected, in which case genome_size is ignored
! concerning the mutational data, however it is still used for computing
! the tracking threshold.
if(genome_size == 0) genome_size = 1.e8

! When tracking_threshold is input as zero, this means that all
! mutations are to be tracked.  In this case, set tracking_threshold
! to be equal to the minimum fitness value to prevent various
! numerical overflow problems that would arise otherwise.

tracking_threshold = max(1./genome_size,tracking_threshold)

lb_modulo  = (2**30-2)/num_linkage_subunits

alpha_del  = log(genome_size)
if(max_fav_fitness_gain > 0.) then
   alpha_fav  = log(genome_size*max_fav_fitness_gain)
else
   alpha_fav = alpha_del
end if

gamma_del  = log(-log(high_impact_mutn_threshold)/alpha_del) &
            /log(high_impact_mutn_fraction)

gamma_fav  = log(-log(high_impact_mutn_threshold)/alpha_fav) &
            /log(high_impact_mutn_fraction)

if(tracking_threshold == 1./genome_size) then
   del_scale = 1./(lb_modulo - 2)
   fav_scale = 1./(lb_modulo - 2)
elseif(tracking_threshold == 1.0) then ! dont track any mutations
   del_scale = 0.
   fav_scale = 0.
else
   del_scale = exp(log(-log(tracking_threshold)/alpha_del) &
               /gamma_del)/(lb_modulo-2)
   if(max_fav_fitness_gain > 0.) then
      fav_scale = exp(log(-log(tracking_threshold   &
                  /max_fav_fitness_gain)/alpha_fav) &
                  /gamma_fav)/(lb_modulo-2)
   else
      fav_scale = 0.
   end if
end if

! Compute mean absolute fitness effect for deleterious mutations.

sum = 0.
d2  = 1.

do i=1,1000000
   d1 = d2
   d2 = exp(-alpha_del*(0.000001*i)**gamma_del)
   sum = sum + d1 + d2
end do

del_mean = 0.0000005*sum

! Compute mean absolute fitness effect for favorable mutations.

sum = 0.
d2  = 1.

do i=1,1000000
   d1 = d2
   d2 = exp(-alpha_fav*(0.000001*i)**gamma_fav)
   sum = sum + d1 + d2
end do

fav_mean = 0.0000005*sum*max_fav_fitness_gain

if(max_fav_fitness_gain > 0.) then
   alpha = alpha_fav
   gamma = gamma_fav
else
   alpha = 0.
   gamma = 0.
end if

! Print some properties of the fitness effect distribution.
if(fitness_distrib_type == 1) then
   do i=1,2
      k = 3 + 3*i
      if((myid==0 .and. i==1) .or. i==2) then

      write(k, &
      '("  Properties of the Weibull fitness effect distribution ", &
        "function:"//                                               &
      "              e(x) = exp(-alpha*x**gamma), 0 < x < 1"//      &
      "  genome_size       = ",1pe9.2/                              &
      "  e_high_impact     = ",1pe9.2,"   defining value of *high " &
      "impact* mutation"/ "  frac_high_impact  = ",1pe9.2,          &
      "   fraction *high impact* mutations of total"//              &
      "  mutation type:        deleterious  favorable "/            &
      "  alpha             = ",0p2f12.5,"   log(genome_size) for"/  &
                                            52x,"deleterious case"/ &
      "  gamma             = ",0p2f12.6// &
      "  mean   effect     = ",1p2e12.2/  &
      "  median effect     = ",1p2e12.2,"   (x = 0.5)"/     &
      "   0th   percentile = ",1p2e12.2,"   (x = 1.0)"/     &
      "  10th   percentile = ",1p2e12.2,"   (x = 0.9)"/     &
      "  20th   percentile = ",1p2e12.2,"   (x = 0.8)"/     &
      "  30th   percentile = ",1p2e12.2,"   (x = 0.7)"/     &
      "  40th   percentile = ",1p2e12.2,"   (x = 0.6)"/     &
      "  50th   percentile = ",1p2e12.2,"   (x = 0.5)"/     &
      "  60th   percentile = ",1p2e12.2,"   (x = 0.4)"/     &
      "  70th   percentile = ",1p2e12.2,"   (x = 0.3)"/     &
      "  80th   percentile = ",1p2e12.2,"   (x = 0.2)"/     &
      "  90th   percentile = ",1p2e12.2,"   (x = 0.1)"/     &
      "  99th   percentile = ",1p2e12.2,"   (x = 0.01)"/    &
      "  99.9   percentile = ",1p2e12.2,"   (x = 0.001)"/   &
      "  99.99  percentile = ",1p2e12.2,"   (x = 0.0001)"/  &
      "  99.999 percentile = ",1p2e12.2,"   (x = 0.00001)"//&
      "  Notes:"/ &
      "  (1) The e(x) values above are for a homozygous pair of ", &
         "mutations."/ &
      "  (2) For favorables, e(x) also includes the factor ", &
         "max_fav_fitness_gain."/)') &
      genome_size, high_impact_mutn_threshold,     &
      high_impact_mutn_fraction, alpha_del, alpha, &
      gamma_del, gamma, -del_mean, fav_mean,       &
      -exp(-alpha_del*0.5**gamma_del),             &
       exp(-alpha_fav*0.5**gamma_fav)*max_fav_fitness_gain, &
      -exp(-alpha_del*1.0**gamma_del), &
       exp(-alpha_fav*1.0**gamma_fav)*max_fav_fitness_gain,&
      -exp(-alpha_del*0.9**gamma_del), &
       exp(-alpha_fav*0.9**gamma_fav)*max_fav_fitness_gain,&
      -exp(-alpha_del*0.8**gamma_del), &
       exp(-alpha_fav*0.8**gamma_fav)*max_fav_fitness_gain,&
      -exp(-alpha_del*0.7**gamma_del), &
       exp(-alpha_fav*0.7**gamma_fav)*max_fav_fitness_gain,&
      -exp(-alpha_del*0.6**gamma_del), &
       exp(-alpha_fav*0.6**gamma_fav)*max_fav_fitness_gain,&
      -exp(-alpha_del*0.5**gamma_del), &
       exp(-alpha_fav*0.5**gamma_fav)*max_fav_fitness_gain,&
      -exp(-alpha_del*0.4**gamma_del), &
       exp(-alpha_fav*0.4**gamma_fav)*max_fav_fitness_gain,&
      -exp(-alpha_del*0.3**gamma_del), &
       exp(-alpha_fav*0.3**gamma_fav)*max_fav_fitness_gain,&
      -exp(-alpha_del*0.2**gamma_del), &
       exp(-alpha_fav*0.2**gamma_fav)*max_fav_fitness_gain,&
      -exp(-alpha_del*0.1**gamma_del), &
       exp(-alpha_fav*0.1**gamma_fav)*max_fav_fitness_gain,&
      -exp(-alpha_del*0.01**gamma_del), &
       exp(-alpha_fav*0.01**gamma_fav)*max_fav_fitness_gain,&
      -exp(-alpha_del*0.001**gamma_del),&
       exp(-alpha_fav*0.001**gamma_fav)*max_fav_fitness_gain,&
      -exp(-alpha_del*0.0001**gamma_del), &
       exp(-alpha_fav*0.0001**gamma_fav)*max_fav_fitness_gain,&
      -exp(-alpha_del*0.00001**gamma_del), &
       exp(-alpha_fav*0.00001**gamma_fav)*max_fav_fitness_gain
      end if
   end do
endif

fraction_del_tracked = del_scale*(lb_modulo-2)
if(tracking_threshold == 1./genome_size) fraction_del_tracked = 1.
if(myid == 0) then
   write(6,*) " Tracking threshold = ", tracking_threshold
   write(6,*) " Fraction deleterious mutations tracked = ", &
                fraction_del_tracked
   write(6,*) " Fraction favorable   mutations tracked = ", &
                fav_scale*(lb_modulo-2)
end if
write(9,*) " Tracking threshold = ", tracking_threshold
write(9,*) " Fraction deleterious mutations tracked = ", &
             fraction_del_tracked
write(9,*) " Fraction favorable   mutations tracked = ", &
             fav_scale*(lb_modulo-2)

! Impose a reasonable limit of the number of tracked mutations
! based on the number of new mutations per offspring, the number
! of generations, and the fraction of deleterious mutations tracked.
! If the run is a restart run, double the number again.  Limit the
! number by the input value for max_tracted_mutn_per_indiv.

k = 1.8*mutn_rate*num_generations*fraction_del_tracked
if(restart_case) k = 2*k

!     max_del_mutn_per_indiv = min(k, max_del_mutn_per_indiv)
!     max_fav_mutn_per_indiv = min(k, max_fav_mutn_per_indiv)

max_del_mutn_per_indiv = max_del_mutn_per_indiv &
                       + 2*num_contrasting_alleles
max_fav_mutn_per_indiv = max_fav_mutn_per_indiv &
                       + 2*num_contrasting_alleles

! Prevent array overflow for cases with large numbers of initial
! favorable mutations.

max_fav_mutn_per_indiv = max(max_fav_mutn_per_indiv, &
                             30*num_initial_fav_mutn/pop_size)

! When the parameter tracking_threshold is set to one, it signals
! that no mutations are to be tracked.  In this case, the size of
! the dmutn and fmutn arrays can be reduced to a minimum.

if(tracking_threshold == 1.0) then
   max_del_mutn_per_indiv = 4
   max_fav_mutn_per_indiv = 4
   fraction_recessive     = 0.
   dominant_hetero_expression  = 0.5
   recessive_hetero_expression = 0.5
end if

if(myid == 0) &
   write(6,'(" Maximum  deleterious mutations tracked = ",i8/   &
             " Maximum  beneficial  mutations tracked = ",i8/   &
             " Maximum  neutral     mutations tracked = ",i8)') &
               max_del_mutn_per_indiv, max_fav_mutn_per_indiv,  &
               max_neu_mutn_per_indiv
   write(9,'(" Maximum  deleterious mutations tracked = ",i8/   &
             " Maximum  beneficial  mutations tracked = ",i8/   &
             " Maximum  neutral     mutations tracked = ",i8)') &
               max_del_mutn_per_indiv, max_fav_mutn_per_indiv,  &
               max_neu_mutn_per_indiv

! Initialize random number generator.

d1 = randomnum(-abs(random_number_seed+myid))
poisson_mean = mutn_rate
k = random_poisson(poisson_mean,.true.)

! Compute initial_alleles_mean_effect from input parameters
initial_alleles_mean_effect = max_total_fitness_increase &
                            / num_contrasting_alleles

! Compute how often to output information
if (num_generations <= 1000) then
   hst_gens = 1
   output_gens = 10
   diagnostic_gens = 20
elseif (num_generations <= 10000) then
   hst_gens = 10
   output_gens = 100
   diagnostic_gens = 200
elseif (num_generations <= 100000) then
   hst_gens = 100
   output_gens = 1000
   diagnostic_gens = 2000
elseif (num_generations <= 1000000) then
   hst_gens = 100
   output_gens = 1000
   diagnostic_gens = 2000
else
   hst_gens = 10000
   output_gens = 100000
   diagnostic_gens = 200000
endif

end subroutine initialize

subroutine gen_initial_contrasting_alleles(dmutn, fmutn, &
   linkage_block_fitness, initial_allele_effects, max_size)

! This routine generates a small number (no larger than the number
! of linkage subunits) of paired alleles, with a random fitness
! effect on one haplotype set and an effect with the same magnitude
! but the opposite sign on the other haplotype set.  Variation of
! of fitness effect is according to a uniform random distribution
! with a user-specified mean value.

use random_pkg
use inputs
include 'common.h'

integer max_size
integer dmutn(max_del_mutn_per_indiv/2,2,max_size)
integer fmutn(max_fav_mutn_per_indiv/2,2,max_size)
real*8 linkage_block_fitness(num_linkage_subunits,2,max_size)
real initial_allele_effects(num_linkage_subunits)
real w, effect, expressed
integer h1_id, h2_id, lb, mutn, mutn_indx
integer m, n, nskp
integer np

w = multiplicative_weighting

if(num_contrasting_alleles > 0) then
   num_contrasting_alleles = min(num_linkage_subunits, &
                                 num_contrasting_alleles)
   nskp = num_linkage_subunits/num_contrasting_alleles
else
   return
end if

initial_allele_effects = 0.

np = int(initial_alleles_pop_frac*pop_size)

do n=1,num_contrasting_alleles

   lb = 1 + (n - 1)*nskp

!    Use the same mutation effect index for all of these paired
!    alleles. This index, lb_modulo-1, is reserved exclusively for
!    these alleles.  When treated as an ordinary mutation, the
!    fitness effect it would imply is the smallest effect possible.
!    The fitness effects associated with these alleles, however,
!    are handled via the linkage_block_fitness array.  We tag these
!    alleles with a mutation index to be able to track them over
!    successive generations and obtain statistics on them at the
!    end of a run.

   mutn = lb_modulo - 1

!  Add an offset to assign it to the appropriate linkage block.

   mutn_indx = mutn + (lb - 1)*lb_modulo

!  Generate random haplotype identities.

   h1_id = min(2, 1 + int(2.*randomnum(1)))
   h2_id = 3 - h1_id

   m = dmutn(1,h1_id,1) + 1
   dmutn(m+1,h1_id,1:np) = mutn_indx
   dmutn(  1,h1_id,1:np) = m

   m = fmutn(1,h2_id,1) + 1
   fmutn(m+1,h2_id,1:np) = mutn_indx
   fmutn(  1,h2_id,1:np) = m

!  Generate the uniformly distributed random fitness effect
!  associated with the allele pair.

   effect = 2.*initial_alleles_mean_effect*randomnum(1)
   if(num_contrasting_alleles < 11) &
      effect = initial_alleles_mean_effect

!  Store the value of the fitness effect for each allele pair in
!  array initial_allele_effects.

   initial_allele_effects(lb) = effect

!  We assume deleterious alleles behave in a recessive manner
!  and when heterozygous have an effect given by the allele
!  fitness effect multiplied by recessive_hetero_expression.
!  Similarly, we assume favorable alleles behave in a dominant
!  manner and when heterozygous have an effect given by the
!  allele fitness effect times dominant_hetero_expression.  The
!  full allele fitness effect is realized only when the same
!  version of the allele occurs on both instances of its linkage
!  block, that is, is homozygous.

!  Apply the appropriate fitness effects to the appropriate
!  linkage blocks.

   expressed = recessive_hetero_expression*effect
   linkage_block_fitness(lb,h1_id,1:np) = (1.d0 - (1.-w)*expressed) &
                                         * (1.d0 - w *expressed)

   expressed =  dominant_hetero_expression*effect
   linkage_block_fitness(lb,h2_id,1:np) = (1.d0 + (1.-w)*expressed) &
                                         * (1.d0 + w *expressed)

end do

end subroutine gen_initial_contrasting_alleles

subroutine tabulate_initial_alleles(dmutn, fmutn, &
   linkage_block_fitness, initial_allele_effects, max_size)
use inputs
! This routine generates a file to report the allele distribution
! resulting from a small number (no larger than the number
! of linkage subunits) of paired alleles, with a random fitness
! effect on one haplotype set and an effect with the same magnitude
! but the opposite sign on the other haplotype set.  Variation of
! of fitness effect is according to a uniform random distribution
! with a user-specified mean value.

include 'common.h'
integer max_size
integer dmutn(max_del_mutn_per_indiv/2,2,max_size)
integer fmutn(max_fav_mutn_per_indiv/2,2,max_size)
real*8 linkage_block_fitness(num_linkage_subunits,2,max_size)
real initial_allele_effects(num_linkage_subunits)
real effect, expressed
integer i, lb, m, mutn, mutn_indx, n, ndel, nfav, npath, nskp

if(num_contrasting_alleles > 0) then
   num_contrasting_alleles = min(num_linkage_subunits, &
                                 num_contrasting_alleles)
   nskp = num_linkage_subunits/num_contrasting_alleles
else
   return
end if

!write(30,*) "#      Distribution of Contrasting Alleles"
!write(30,*) "#   block      del       fav    fitness effect  "

do n=1,num_contrasting_alleles

   lb = 1 + (n - 1)*nskp

!  The same mutation effect index is used for all of these paired
!  alleles. This index, lb_modulo-1, is reserved exclusively for
!  these alleles.

   mutn = lb_modulo - 1

!  Add an offset to assign it to the appropriate linkage block.

   mutn_indx = mutn + (lb - 1)*lb_modulo

!  Loop over entire population, count the number of each of the
!  alleles, and output these statistics.

   ndel = 0
   nfav = 0

   do i=1,current_pop_size

      do m=2,dmutn(1,1,i) + 1
         if (dmutn(m,1,i) == mutn_indx) ndel = ndel + 1
      end do

      do m=2,dmutn(1,2,i) + 1
         if (dmutn(m,2,i) == mutn_indx) ndel = ndel + 1
      end do

      do m=2,fmutn(1,1,i) + 1
         if (fmutn(m,1,i) == mutn_indx) nfav = nfav + 1
      end do

      do m=2,fmutn(1,2,i) + 1
         if (fmutn(m,2,i) == mutn_indx) nfav = nfav + 1
      end do

   end do

   !write(30,'(3i10,f15.10)') lb, ndel, nfav, &
   !                          initial_allele_effects(lb)

end do

end subroutine tabulate_initial_alleles

end module init
