subroutine mating(dmutn,nmutn,fmutn,lb_mutn_count,linkage_block_fitness, &
           num_offspring,available,pop_size_allocation,nmax,offspring_count, &
           total_offspring,tsub,gen)
! Randomly mate one half of the population with members
! from the other half.
use selection_module
use random_pkg
use profiler
use inputs
use polygenic
include 'common.h'
integer, intent(inout), dimension(:,:,:) :: dmutn(max_del_mutn_per_indiv/2,2,*)
integer, intent(inout), dimension(:,:,:) :: nmutn(max_neu_mutn_per_indiv/2,2,*)
integer, intent(inout), dimension(:,:,:) :: fmutn(max_fav_mutn_per_indiv/2,2,*)
integer, intent(inout) :: lb_mutn_count(num_linkage_subunits,2,3,*)
real*8,  intent(inout) :: linkage_block_fitness(num_linkage_subunits,2,*)
integer, intent(in)    :: nmax, pop_size_allocation, gen
integer, intent(out)   :: offspring_count, total_offspring
logical, intent(inout), dimension(:) :: available(pop_size_allocation)
real, intent(in)       :: num_offspring
real, intent(out)      :: tsub

integer, allocatable, dimension(:,:,:)   :: dmutn_offsprng
integer, allocatable, dimension(:,:,:)   :: nmutn_offsprng
integer, allocatable, dimension(:,:,:)   :: fmutn_offsprng
integer, allocatable, dimension(:,:,:,:) :: offsprng_lb_mutn_count
logical, allocatable, dimension(:)       :: replaced_by_offspring
real*8,  allocatable, dimension(:,:,:)   :: offsprng_lb_fitness
integer :: i, j, k, hap, dad, mom, actual_offspring, replace
integer :: parent, child, slot, empty
real    :: fitness_adjusted_offspring
logical :: match

allocate( dmutn_offsprng(max_del_mutn_per_indiv/2,2,nmax),     &
          nmutn_offsprng(max_neu_mutn_per_indiv/2,2,nmax),     &
          fmutn_offsprng(max_fav_mutn_per_indiv/2,2,nmax),     &
        offsprng_lb_mutn_count(num_linkage_subunits,2,3,nmax), &
           offsprng_lb_fitness(num_linkage_subunits,2,nmax),   &
           replaced_by_offspring(pop_size_allocation))

call second(tin)

available = .true.
replaced_by_offspring = .false.
post_sel_fitness = 1.d0
total_offspring  = current_pop_size
offspring_count  = 0

do i=1,current_pop_size/2
   
   call random_mate(dad,available,current_pop_size,pop_size_allocation)
   call random_mate(mom,available,current_pop_size,pop_size_allocation)
   ! call nonrandom_mate(dad,available,current_pop_size,pop_size_allocation)
   ! call nonrandom_mate(mom,available,current_pop_size,pop_size_allocation)
   
   ! Generate an appropriate number offspring from the two parents. 
   
   if(fitness_dependent_fertility) then
      fitness_adjusted_offspring = & 
           num_offspring*sqrt(min(1.d0, post_sel_fitness))
      actual_offspring = int(fitness_adjusted_offspring)
      if(fitness_adjusted_offspring - actual_offspring > &
           randomnum(1)) actual_offspring = actual_offspring + 1
   else
      actual_offspring = int(num_offspring)
      if(num_offspring - int(num_offspring) > randomnum(1)) &
           actual_offspring = actual_offspring + 1
   end if
   
   if(i == 1) actual_offspring = max(1, actual_offspring)
   
   actual_offspring = min(int(num_offspring) + 1, &
        actual_offspring)
   
   ! If the parameter fraction_self_fertilization is non-zero
   ! (as can be the case for many types of plants), implement
   ! self-fertilization for the appropriate portion of the
   ! population by calling routine offspring using 'dad' for
   ! both parents for half the offspring and 'mom' for both
   ! parents for the other half.
   
   do child=1,actual_offspring
      
      if(fraction_self_fertilization <= randomnum(1) .and.  &
         recombination_model /= clonal) then
         
         ! This call generates an offspring from a sexual
         ! union of the two individuals dad and mom.
         
         call offspring(dmutn_offsprng(1,1,child), &
              nmutn_offsprng(1,1,child), &
              fmutn_offsprng(1,1,child), &
              offsprng_lb_mutn_count(1,1,1,child), &
              offsprng_lb_fitness(1,1,child), &
              dmutn,nmutn,fmutn,lb_mutn_count, &
              linkage_block_fitness,dad,mom,gen)
         
      elseif(actual_offspring >= 2 .and. child == 1) then
         
         ! This call generates an offspring exclusively
         ! from the genetic makeup of individual -dad-.
         ! If the recombination_model parameter is clonal,
         ! the offspring is a clone of -dad- except for
         ! possible new mutations.  If not, the genotype
         ! of the offspring is the product of gametic 
         ! shuffling of the chromosomes of -dad- via
         ! self-fertilization.
         
         call offspring(dmutn_offsprng(1,1,child), &
              nmutn_offsprng(1,1,child), &
              fmutn_offsprng(1,1,child), &
              offsprng_lb_mutn_count(1,1,1,child), &
              offsprng_lb_fitness(1,1,child), &
              dmutn,nmutn,fmutn,lb_mutn_count, &
              linkage_block_fitness,dad,dad,gen)
         
      elseif(actual_offspring >= 2 .and. child == 2) then
         
         ! This call generates an offspring exclusively
         ! from the genetic makeup of individual -mom-.
         ! If the recombination_model parameter is clonal,
         ! the offspring is a clone of -mom- except for
         ! possible new mutations.  If not, the genotype
         ! of the offspring is the product of gametic 
         ! shuffling of the chromosomes of -mom- via
         ! self-fertilization.
         
         call offspring(dmutn_offsprng(1,1,child), &
              nmutn_offsprng(1,1,child), &
              fmutn_offsprng(1,1,child), &
              offsprng_lb_mutn_count(1,1,1,child), &
              offsprng_lb_fitness(1,1,child), &
              dmutn,nmutn,fmutn,lb_mutn_count, &
              linkage_block_fitness,mom,mom,gen)
         
      else
         
         if(randomnum(1) < 0.5) then
            parent = dad
         else
            parent = mom
         end if
         
         ! This call generates an offspring exclusively
         ! from the genetic makeup of individual parent.
         ! If the recombination_model parameter is clonal,
         ! the offspring is a clone of parent except for
         ! possible new mutations.  If not, the genotype
         ! of the offspring is the product of gametic 
         ! shuffling of the parent's chromosomes via
         ! self-fertilization.
         
         call offspring(dmutn_offsprng(1,1,child), &
              nmutn_offsprng(1,1,child), &
              fmutn_offsprng(1,1,child), &
              offsprng_lb_mutn_count(1,1,1,child), &
              offsprng_lb_fitness(1,1,child), &
              dmutn,nmutn,fmutn,lb_mutn_count, &
              linkage_block_fitness,parent,parent,gen)
         
      end if
      
   end do
   
   offspring_count = offspring_count + actual_offspring
   
   ! Copy mutation list arrays for each of the first two 
   ! offspring into locations of the two parents.  Update
   ! the linkage block mutation count and the linkage block
   ! fitness for these two offspring.
   
   if(actual_offspring >= 1) then
      do hap=1,2
         k = dmutn_offsprng(1,hap,1) + 1
         dmutn(1:k,hap,dad) = dmutn_offsprng(1:k,hap,1)
         
         k = nmutn_offsprng(1,hap,1) + 1
         nmutn(1:k,hap,dad) = nmutn_offsprng(1:k,hap,1)
         
         k = fmutn_offsprng(1,hap,1) + 1
         fmutn(1:k,hap,dad) = fmutn_offsprng(1:k,hap,1)
      end do
      
      lb_mutn_count(:,:,:,dad) =  offsprng_lb_mutn_count(:,:,:,1)
      linkage_block_fitness(:,:,dad) =  offsprng_lb_fitness(:,:,1)
      replaced_by_offspring(dad) = .true.
   end if
   if(actual_offspring >= 2) then
      do hap=1,2
         k = dmutn_offsprng(1,hap,2) + 1
         dmutn(1:k,hap,mom) = dmutn_offsprng(1:k,hap,2)
         
         k = nmutn_offsprng(1,hap,2) + 1
         nmutn(1:k,hap,mom) = nmutn_offsprng(1:k,hap,2)
         
         k = fmutn_offsprng(1,hap,2) + 1
         fmutn(1:k,hap,mom) = fmutn_offsprng(1:k,hap,2)
      end do
      
      lb_mutn_count(:,:,:,mom) =  offsprng_lb_mutn_count(:,:,:,2)
      linkage_block_fitness(:,:,mom) =  offsprng_lb_fitness(:,:,2)
      replaced_by_offspring(mom) = .true.
   end if
   
   ! Copy the mutation list for any other offspring into arrays
   ! dmutn and fmutn with an index greater than current_pop_size.
   ! Update array linkage_block_fitness appropriately.
   
   do slot=3,actual_offspring
      total_offspring = total_offspring + 1
      j = total_offspring 
      
      do hap=1,2
         k = dmutn_offsprng(1,hap,slot) + 1
         dmutn(1:k,hap,j) = dmutn_offsprng(1:k,hap,slot)
         
         k = nmutn_offsprng(1,hap,slot) + 1
         nmutn(1:k,hap,j) = nmutn_offsprng(1:k,hap,slot)
         
         k = fmutn_offsprng(1,hap,slot) + 1
         fmutn(1:k,hap,j) = fmutn_offsprng(1:k,hap,slot)
      end do
      
      lb_mutn_count(:,:,:,j) =       offsprng_lb_mutn_count(:,:,:,slot)
      linkage_block_fitness(:,:,j) = offsprng_lb_fitness(:,:,slot)
   end do
   
end do

! For slots not overwritten by new offspring, move data from
! offspring with higher index numbers to populate these slots.

if(offspring_count < current_pop_size) then
   
   replace = current_pop_size
   empty   = 1
   
   do while(empty  <= offspring_count .and.  &
        replace > offspring_count)
      
      do while(replaced_by_offspring(empty) .and. &
           empty < offspring_count)
         empty = empty + 1
      end do
      
      do while(replace <= current_pop_size .and. &
           .not.replaced_by_offspring(replace))
         replace = replace - 1
      end do
      
      if(empty <=  offspring_count .and.  &
           replace > offspring_count) then
         
         do hap=1,2
            k = dmutn(1,hap,replace) + 1
            dmutn(1:k,hap,empty) = dmutn(1:k,hap,replace)
            
            k = nmutn(1,hap,replace) + 1
            nmutn(1:k,hap,empty) = nmutn(1:k,hap,replace)
            
            k = fmutn(1,hap,replace) + 1
            fmutn(1:k,hap,empty) = fmutn(1:k,hap,replace)
         end do
         
         lb_mutn_count(:,:,:,empty)       = lb_mutn_count(:,:,:,replace)
         linkage_block_fitness(:,:,empty) = linkage_block_fitness(:,:,replace)
         replaced_by_offspring(empty) = .true.
         replace = replace - 1
      end if
      
   end do
   
   total_offspring = offspring_count
   
end if

call second(tout)
tsub = tout - tin
sec(5) = sec(5) + tsub

end subroutine mating


subroutine offspring(dmutn_offsprng,nmutn_offsprng,fmutn_offsprng, &
                     offsprng_lb_mutn_count, &
                     offsprng_lb_fitness,dmutn,nmutn,fmutn, &
                     lb_mutn_count,linkage_block_fitness,dad,mom,gen)

 ! This routine creates an offspring from a pair of individuals,
 ! indexed 'dad' and 'mom', in the current population.  This 
 ! offspring inherits one set of mutations from each linkage block
 ! pair from each parent.  This mutation information is loaded into
 ! array 'offsprng'.  Also, a set of new deleterious mutations, 
 ! equal in number to mutn_rate*(1. - frac_fav_mutn) and
 ! chosen randomly, are generated for this offspring.  The fitness
 ! associated with each of these new mutations is applied to update
 ! the fitness of the linkage block in which it occurs. 

use polygenic
use random_pkg
use inputs
include 'common.h'

integer del, fav, neu
parameter (del = 1, fav = 2, neu = 3)

integer dmutn_offsprng(max_del_mutn_per_indiv/2,2), &
        nmutn_offsprng(max_neu_mutn_per_indiv/2,2), &
        fmutn_offsprng(max_fav_mutn_per_indiv/2,2)
integer          dmutn(max_del_mutn_per_indiv/2,2,*), &
                 nmutn(max_neu_mutn_per_indiv/2,2,*), &
                 fmutn(max_fav_mutn_per_indiv/2,2,*), &
      offsprng_lb_mutn_count(num_linkage_subunits,2,3), & 
               lb_mutn_count(num_linkage_subunits,2,3,*)
real*8   offsprng_lb_fitness(num_linkage_subunits,2)
real*8 linkage_block_fitness(num_linkage_subunits,2,*)
real*8 fitness_effect, se_effect, ran1, ran2, ran3, ran4, rtmr
real w, y, chance_back_mutn
integer dad, mom, chr_length, ch, ls0, ls1, ls2, iseg, iseg_max, gen
integer md1, md2, mf1, mf2, mn1, mn2
integer mdd_off, mfd_off, mnd_off, mdm_off, mfm_off, mnm_off
integer lb, hap_id, new_mutn, mut, mutn_indx, num_mutn, i
integer lb_syn, hap_id_syn, mutn_sum, mutn_type, string(40)
logical is_back_mutn
logical match, had_target

w = multiplicative_weighting

dmutn_offsprng(1,:) = 0
fmutn_offsprng(1,:) = 0
if(.not.polygenic_beneficials) nmutn_offsprng(1,:) = 0

if(recombination_model == full_sexual) then

   ! Set the number of segments.  Three sections of the chromosome that are 
   ! involved in the crossover.  Form the gametes chromosome by chromosome.
   if(dynamic_linkage) then
      iseg_max = 3
   else
      iseg_max = 1 ! can come from any parent
   end if
   
   ! Number of linkage blocks on one chromosome
   chr_length = num_linkage_subunits/haploid_chromosome_number
   
   ! Initialize mutation counters.  The first slot of mutation arrays 
   ! is used to store the total number of mutations, so start with 2.
   md1     = 2
   md2     = 2
   mf1     = 2
   mf2     = 2
   mn1     = 2
   mn2     = 2
   mdd_off = 2
   mfd_off = 2
   mnd_off = 2
   
   ! Loop over the total number of chromosomes
   do ch=1,haploid_chromosome_number
      
      ! Markers along the chromosome
      ! total array chopped up into linkage blocks
      ! ls0 points to the first linkage block for a given chromosome
      ls0 = (ch - 1)*chr_length + 1
      ! ls1 and ls2 mark the locations within each chromosome where the beginning
      ! and ending of the crossover portion occurs, from the other parent
      ls1 = min(chr_length-1, int(chr_length*randomnum(1))) + ls0
      ls2 = min(chr_length-1, int(chr_length*randomnum(1))) + ls0
      
      ! Reverse the numbers in the case that ls1 is greater than ls2
      if(ls1 > ls2) then
         lb  = ls1
         ls1 = ls2
         ls2 = lb
      end if
      
      if(.not.dynamic_linkage) ls1 = num_linkage_subunits
      
      ! If dynamic_linkage loop 3 times, else loop 1 time
      do iseg=1,iseg_max
         
         if(iseg == 2) then
            ls0 = ls1 + 1
            ls1 = ls2
         elseif(iseg == 3) then
            ls0 = ls2 + 1
            ls1 = ch*chr_length
         end if
         
         ! If hap_id is 1, for this segment the mutations are going to come
         ! from the father, and not the mother, depending on the value of iseg
         ! If hap_id is 2, from the second chromosome of the father.
         if(dynamic_linkage) hap_id = min(2, 1 + int(2.*randomnum(1)))
         
         do lb=ls0,ls1
            
            ! Copy the mutation list from the randomly selected haplotype
            ! from the father to form gamete one for the offspring.  Also
            ! copy the corresponding linkage block mutation count and
            ! fitness.
            
            if(.not.dynamic_linkage) hap_id = min(2, 1 + int(2.*randomnum(1)))
            
            if(tracking_threshold /= 1.) then
               
               ! deleterious mutations chromosome 1
               do while(abs(dmutn(md1,1,dad)) < lb*lb_modulo .and. &
                    md1 <= dmutn(1,1,dad) + 1)
                  if(hap_id == 1) then
                     dmutn_offsprng(mdd_off,1) = dmutn(md1,1,dad)
                     mdd_off = mdd_off + 1
                  end if
                  md1 = md1 + 1
               end do
               
               ! deleterious mutations chromosome 2
               do while(abs(dmutn(md2,2,dad)) < lb*lb_modulo .and. &
                    md2 <= dmutn(1,2,dad) + 1)
                  if(hap_id == 2) then
                     dmutn_offsprng(mdd_off,1) = dmutn(md2,2,dad)
                     mdd_off = mdd_off + 1
                  end if
                  md2 = md2 + 1
               end do
               
               ! favorable mutations chromosome 1
               do while(abs(fmutn(mf1,1,dad)) < lb*lb_modulo .and. &
                    mf1 <= fmutn(1,1,dad) + 1)
                  if(hap_id == 1) then
                     fmutn_offsprng(mfd_off,1) = fmutn(mf1,1,dad)
                     mfd_off = mfd_off + 1
                  end if
                  mf1 = mf1 + 1
               end do
               
               ! favorable mutations chromosome 2
               do while(abs(fmutn(mf2,2,dad)) < lb*lb_modulo .and. &
                    mf2 <= fmutn(1,2,dad) + 1)
                  if(hap_id == 2) then
                     fmutn_offsprng(mfd_off,1) = fmutn(mf2,2,dad)
                     mfd_off = mfd_off + 1
                  end if
                  mf2 = mf2 + 1
               end do
               
               ! neutral mutations chromosome 1
               do while(abs(nmutn(mn1,1,dad)) < lb*lb_modulo .and. &
                    mn1 <= nmutn(1,1,dad) + 1)
                  if(hap_id == 1) then
                     nmutn_offsprng(mnd_off,1) = nmutn(mn1,1,dad)
                     mnd_off = mnd_off + 1
                  end if
                  mn1 = mn1 + 1
               end do
               
               ! neutral mutations chromosome 2
               do while(abs(nmutn(mn2,2,dad)) < lb*lb_modulo .and. &
                    mn2 <= nmutn(1,2,dad) + 1)
                  if(hap_id == 2) then
                     nmutn_offsprng(mnd_off,1) = nmutn(mn2,2,dad)
                     mnd_off = mnd_off + 1
                  end if
                  mn2 = mn2 + 1
               end do
               
            end if
            
            offsprng_lb_mutn_count(lb,1,1) = lb_mutn_count(lb,hap_id,1,dad)
            offsprng_lb_mutn_count(lb,1,2) = lb_mutn_count(lb,hap_id,2,dad)
            offsprng_lb_mutn_count(lb,1,3) = lb_mutn_count(lb,hap_id,3,dad)
            offsprng_lb_fitness(lb,1) = linkage_block_fitness(lb,hap_id,dad)
         end do
         
      end do
      
   end do
   
   md1     = 2
   md2     = 2
   mf1     = 2
   mf2     = 2
   mn1     = 2
   mn2     = 2
   mdm_off = 2
   mfm_off = 2
   mnm_off = 2
   
   do ch=1,haploid_chromosome_number
      
      ls0 = (ch - 1)*chr_length + 1
      ls1 = min(chr_length-1, int(chr_length*randomnum(1))) + ls0
      ls2 = min(chr_length-1, int(chr_length*randomnum(1))) + ls0
      
      if(ls1 > ls2) then
         lb  = ls1
         ls1 = ls2
         ls2 = lb
      end if
      
      if(.not.dynamic_linkage) ls1 = chr_length
      
      do iseg=1,iseg_max
         
         if(iseg == 2) then
            ls0 = ls1 + 1
            ls1 = ls2
         elseif(iseg == 3) then
            ls0 = ls2 + 1
            ls1 = ch*chr_length
         end if
         
         if(dynamic_linkage) hap_id = min(2, 1 + int(2.*randomnum(1)))
         
         do lb=ls0,ls1
            
            ! Copy the mutation list from the randomly selected haplotype
            ! from the mother to form gamete two for the offspring.  Also
            ! copy the corresponding linkage block mutation count and
            ! fitness.
            
            if(.not.dynamic_linkage) hap_id = min(2, 1 + int(2.*randomnum(1)))
            
            if(tracking_threshold /= 1.) then
               
               do while(abs(dmutn(md1,1,mom)) < lb*lb_modulo .and. &
                    md1 <= dmutn(1,1,mom) + 1)
                  if(hap_id == 1) then
                     dmutn_offsprng(mdm_off,2) = dmutn(md1,1,mom)
                     mdm_off = mdm_off + 1
                  end if
                  md1 = md1 + 1
               end do
               
               do while(abs(dmutn(md2,2,mom)) < lb*lb_modulo .and. &
                    md2 <= dmutn(1,2,mom) + 1)
                  if(hap_id == 2) then
                     dmutn_offsprng(mdm_off,2) = dmutn(md2,2,mom)
                     mdm_off = mdm_off + 1
                  end if
                  md2 = md2 + 1
               end do
               
               do while(abs(fmutn(mf1,1,mom)) < lb*lb_modulo .and. &
                    mf1 <= fmutn(1,1,mom) + 1)
                  if(hap_id == 1) then
                     fmutn_offsprng(mfm_off,2) = fmutn(mf1,1,mom)
                     mfm_off = mfm_off + 1
                  end if
                  mf1 = mf1 + 1
               end do
               
               do while(abs(fmutn(mf2,2,mom)) < lb*lb_modulo .and. &
                    mf2 <= fmutn(1,2,mom) + 1)
                  if(hap_id == 2) then
                     fmutn_offsprng(mfm_off,2) = fmutn(mf2,2,mom)
                     mfm_off = mfm_off + 1
                  end if
                  mf2 = mf2 + 1
               end do
               
               do while(abs(nmutn(mn1,1,mom)) < lb*lb_modulo .and. &
                    mn1 <= nmutn(1,1,mom) + 1)
                  if(hap_id == 1) then
                     nmutn_offsprng(mnm_off,2) = nmutn(mn1,1,mom)
                     mnm_off = mnm_off + 1
                  end if
                  mn1 = mn1 + 1
               end do
               
               do while(abs(nmutn(mn2,2,mom)) < lb*lb_modulo .and. &
                    mn2 <= nmutn(1,2,mom) + 1)
                  if(hap_id == 2) then
                     nmutn_offsprng(mnm_off,2) = nmutn(mn2,2,mom)
                     mnm_off = mnm_off + 1
                  end if
                  mn2 = mn2 + 1
               end do
               
            end if
            
            offsprng_lb_mutn_count(lb,2,1) = lb_mutn_count(lb,hap_id,1,mom)
            offsprng_lb_mutn_count(lb,2,2) = lb_mutn_count(lb,hap_id,2,mom)
            offsprng_lb_mutn_count(lb,2,3) = lb_mutn_count(lb,hap_id,3,mom)
            offsprng_lb_fitness(lb,2) = linkage_block_fitness(lb,hap_id,mom)
         end do
         
      end do
      
   end do
   
   ! Load the mutation count into the first location in the first 
   ! index of each of the mutation arrays.
   
   dmutn_offsprng(1,1) = mdd_off - 2
   dmutn_offsprng(1,2) = mdm_off - 2
   fmutn_offsprng(1,1) = mfd_off - 2
   fmutn_offsprng(1,2) = mfm_off - 2
   nmutn_offsprng(1,1) = mnd_off - 2
   nmutn_offsprng(1,2) = mnm_off - 2

elseif(recombination_model == clonal) then 
   
   mdd_off = dmutn(1,1,dad)
   mdm_off = dmutn(1,2,dad)
   mfd_off = fmutn(1,1,dad)
   mfm_off = fmutn(1,2,dad)
   mnd_off = nmutn(1,1,dad)
   mnm_off = nmutn(1,2,dad)

   dmutn_offsprng(1:mdd_off+1,1) = dmutn(1:mdd_off+1,1,dad)
   dmutn_offsprng(1:mdm_off+1,2) = dmutn(1:mdm_off+1,2,dad)
   fmutn_offsprng(1:mfd_off+1,1) = fmutn(1:mfd_off+1,1,dad)
   fmutn_offsprng(1:mfm_off+1,2) = fmutn(1:mfm_off+1,2,dad)
   nmutn_offsprng(1:mnd_off+1,1) = nmutn(1:mnd_off+1,1,dad)
   nmutn_offsprng(1:mnm_off+1,2) = nmutn(1:mnm_off+1,2,dad)

   offsprng_lb_mutn_count(:,:,:) = lb_mutn_count(:,:,:,dad)
   offsprng_lb_fitness(:,:) = linkage_block_fitness(:,:,dad)

elseif (polygenic_beneficials .or. recombination_model == suppressed) then

   ! Choose a random haplotype from the father
   hap_id = min(2, 1 + int(2.*randomnum(1)))

   mdd_off = dmutn(1,hap_id,dad)
   mfd_off = fmutn(1,hap_id,dad)
   mnd_off = nmutn(1,hap_id,dad)

   dmutn_offsprng(1:mdd_off+1,1) = dmutn(1:mdd_off+1,hap_id,dad)
   fmutn_offsprng(1:mfd_off+1,1) = fmutn(1:mfd_off+1,hap_id,dad)
   nmutn_offsprng(1:mnd_off+1,1) = nmutn(1:mnd_off+1,hap_id,dad)

   offsprng_lb_mutn_count(:,1,:) = lb_mutn_count(:,hap_id,:,dad)
   offsprng_lb_fitness(:,1) = linkage_block_fitness(:,hap_id,dad)

   ! Choose a random haplotype from the mother
   hap_id = min(2, 1 + int(2.*randomnum(1)))

   mdm_off = dmutn(1,hap_id,mom)
   mfm_off = fmutn(1,hap_id,mom)
   mnm_off = nmutn(1,hap_id,mom)

   dmutn_offsprng(1:mdm_off+1,2) = dmutn(1:mdm_off+1,hap_id,mom)
   fmutn_offsprng(1:mfm_off+1,2) = fmutn(1:mfm_off+1,hap_id,mom)
   nmutn_offsprng(1:mnm_off+1,2) = nmutn(1:mnm_off+1,hap_id,mom)

   offsprng_lb_mutn_count(:,2,:) = lb_mutn_count(:,hap_id,:,mom)
   offsprng_lb_fitness(:,2) = linkage_block_fitness(:,hap_id,mom)

else

   write(6,*) 'ERROR: reproduction_model not supported'
   stop

end if

! Generate new_mutn new random mutations in randomly chosen
! linkage block locations, where new_mutn is a random deviate
! drawn from a Poisson distribution with a mean value given by
! the parameter poisson_mean, which under most circumstances is
! very close to the value of mutn_rate. 
! For mutation rates smaller that 1.e-4, use a simpler method
! based on the probability of a random number occurring within
! a tiny interval between 0. and 1.

if(mutn_rate <= 1.e-4) then
   rtmr = sqrt(mutn_rate)
   ran1 = randomnum(1)
   ran2 = randomnum(1)
   ran3 = min(randomnum(1), 1. - rtmr)
   ran4 = min(randomnum(1), 1. - rtmr)
   if((ran1 > ran3 .and. ran1 <= ran3+rtmr) .and.  &
      (ran2 > ran4 .and. ran2 <= ran4+rtmr)) then
      new_mutn = 1
   else
      new_mutn = 0
   end if
else
   if (poisson_method == 1) then
      new_mutn = random_poisson(poisson_mean, .false.)
   else
      new_mutn = poisson(poisson_mean)
   end if
end if

new_mutn_count = new_mutn_count + new_mutn

do mut=1,new_mutn

   is_back_mutn = .false.

   if(allow_back_mutn) then

      ! Possibly substitute a back mutation for a normal one.
      ! The chance of a back mutation = total number of mutations 
      ! in an individual divided by genome size divided by three.

      mutn_sum = fmutn_offsprng(1,1) + fmutn_offsprng(1,2) &
               + dmutn_offsprng(1,1) + dmutn_offsprng(1,2)
      chance_back_mutn = mutn_sum/(2.*genome_size*3.)

      if(randomnum(1) < chance_back_mutn) then

         is_back_mutn = .true.

         call back_mutn(dmutn_offsprng, fmutn_offsprng, &
              offsprng_lb_fitness, offsprng_lb_mutn_count)

         num_back_mutn = num_back_mutn + 1

      end if

   end if

   if(.not.is_back_mutn) then

      call mutation(mutn_indx,mutn_type,lb,hap_id,fitness_effect)

      ! Track this mutation if its fitness effect exceeds the value 
      ! of tracking_threshold or if track_neutrals is true. 
      
      if(fitness_effect>tracking_threshold .or. track_neutrals) then
         
         if(mutn_type == fav) then
            
            ! Test to see if the storage limit of array fmutn_offsprng 
            ! has been exceeded.  (Note that we are using the first  
            ! slot to hold the actual mutation count.)
            
            num_mutn = fmutn_offsprng(1,hap_id) + 1 
            
            if(num_mutn + 1 > max_fav_mutn_per_indiv/2) then
               write(6,*) 'Favorable mutation count exceeds limit'
               write(9,*) 'Favorable mutation count exceeds limit'
               stop
            end if
            
            fmutn_offsprng(1,hap_id) = num_mutn
            
            ! Insert new mutation such that mutations are maintained
            ! in ascending order of their absolute value.
            
            i = num_mutn
            
            do while(abs(fmutn_offsprng(i,hap_id)) > abs(mutn_indx) .and. i > 1)
               fmutn_offsprng(i+1,hap_id) = fmutn_offsprng(i,hap_id)
               i = i - 1
            end do
            
            fmutn_offsprng(i+1,hap_id) = mutn_indx
            
         elseif(mutn_type == del) then
            
            ! Test to see if the storage limit of array dmutn_offsprng 
            ! has been exceeded.  (Note that we are using the first  
            ! slot to hold the actual mutation count.)
            
            num_mutn = dmutn_offsprng(1,hap_id) + 1 
            
            if(num_mutn + 1 > max_del_mutn_per_indiv/2) then
               write(6,*) 'Deleterious mutation count exceeds limit'
               write(9,*) 'Deleterious mutation count exceeds limit'
               stop
            end if
            
            dmutn_offsprng(1,hap_id) = num_mutn
            
            ! Insert new mutation such that mutations are maintained
            ! in ascending order of their absolute value.
            
            i = num_mutn
            
            do while(abs(dmutn_offsprng(i,hap_id)) > abs(mutn_indx) .and. i > 1)
               dmutn_offsprng(i+1,hap_id) = dmutn_offsprng(i,hap_id)
               i = i - 1
            end do
            
            dmutn_offsprng(i+1,hap_id) = mutn_indx
            
         else  ! neutral

            if(polygenic_beneficials) then

               ! For the case of clonal haploid, always choose hap_id = 1.

               if(recombination_model == clonal) hap_id = 1

               ! In the treatment of polygenic beneficials, don't accumulate
               ! mutations but instead maintain a single nucleotide in each 
               ! linkage block and allow that nucleotide to mutate.
                
               ! Generate a new mutation if the mutation is the same nucleotide
               ! as what already exists (at that specific position).  This
               ! implies that we are allowing back mutation. However, for an 
               ! individual whose haplotype matches the target string, do 
               ! not allow such back mutation. (fmutn_offsprng(2,hapid) == 0
               ! implies that the target string has not yet been matched.)

               if(fmutn_offsprng(2,hap_id) == 0) then

                  do while(mutn_indx == nmutn_offsprng(lb+1,hap_id))
                     mutn_indx = min(4, 1 + int(4*randomnum(1)))
                  end do

                  nmutn_offsprng(lb+1,hap_id) = mutn_indx

                  ! Test to see if this mutation causes the haplotype string
                  ! to match the target string.

                  string = nmutn_offsprng(2:num_linkage_subunits+1,hap_id)

                  match  = poly_match(string)

                  ! If the string for this haplotype matches the target,
                  ! set the mutation index for linkage block one in the
                  ! offspring favorable mutation array to a positive random
                  ! number.  The fact that the number is positive and not zero
                  ! signals that the haplotype string is a match. The random
                  ! number gives this haplotype a unique identity and enables
                  ! it and future offspring carrying to be tracked.  Also,
                  ! set the favorable lb_mutn_count to 1 for this offspring.

                  if(match) then
                     fmutn_offsprng(2,hap_id) = int(1.d6*randomnum(1))
                     num_polys_cumulative = num_polys_cumulative + 1
                     pmutn(num_polys_cumulative)%id = fmutn_offsprng(2,hap_id)
                     pmutn(num_polys_cumulative)%gen_enter = gen
                     if(poly_gen_first_instance < 0)  &
                        poly_gen_first_instance = gen
                     poly_gen_last_instance  = gen
                     offsprng_lb_mutn_count(lb,hap_id,2) = 1
                  else
                     fmutn_offsprng(2,hap_id) = 0
                     offsprng_lb_mutn_count(lb,hap_id,2) = 0
                  end if

               end if

            else

               num_mutn = nmutn_offsprng(1,hap_id) + 1 
                
               if(num_mutn + 1 > max_neu_mutn_per_indiv/2) then
                  write(6,*) 'Neutral mutation count exceeds limit'
                  write(9,*) 'Neutral mutation count exceeds limit'
                  stop
               end if
               
               nmutn_offsprng(1,hap_id) = num_mutn
               
               ! Insert new mutation such that mutations are maintained
               ! in ascending order of their absolute value.
               
               i = num_mutn
               
               do while(abs(nmutn_offsprng(i,hap_id)) > abs(mutn_indx) .and. i > 1)
                  nmutn_offsprng(i+1,hap_id) = nmutn_offsprng(i,hap_id)
                  i = i - 1
               end do
               
               nmutn_offsprng(i+1,hap_id) = mutn_indx

            end if

         end if
         
      end if
      
      ! Recessive mutations (identified as such with a negative
      ! mutation index) here incur only recessive_hetero_expression
      ! times of their fitness effect, while dominant mutations incur
      ! only dominant_hetero_expression times their fitness effect,
      ! relative to the case of heterozygous expression.  The full 
      ! fitness effect is realized only when a mutation occurs in  
      ! both instances of its linkage block, that is, is homozygous.
      
      if(mutn_indx < 0) then
         fitness_effect = recessive_hetero_expression*fitness_effect
      else
         fitness_effect =  dominant_hetero_expression*fitness_effect
      end if
      
      ! Update linkage subunit fitness.
      
      if(mutn_type == fav) then
         
         offsprng_lb_fitness(lb,hap_id) = &
              (offsprng_lb_fitness(lb,hap_id) + (1. - w)*fitness_effect) &
              * (1.d0 + w *fitness_effect)
         
         ! Increment the mutation count for the linkage subunit
         ! in which the mutation occurs.
         
         offsprng_lb_mutn_count(lb,hap_id,2) =  &
              offsprng_lb_mutn_count(lb,hap_id,2) + 1
         
      elseif(mutn_type == del) then
         
         offsprng_lb_fitness(lb,hap_id) =  &
              (offsprng_lb_fitness(lb,hap_id) - (1. - w)*fitness_effect) &
              * (1.d0 - w *fitness_effect)
         
         ! Increment the mutation count for the linkage subunit
         ! in which the mutation occurs.
         
         offsprng_lb_mutn_count(lb,hap_id,1) =  &
              offsprng_lb_mutn_count(lb,hap_id,1) + 1
         
      else  ! neutral 
         
         ! polygenic beneficials maintain a single neutral mutation for each
         ! linkage, so don't increment in such cases
         if (.not. polygenic_beneficials) then
            offsprng_lb_mutn_count(lb,hap_id,3) = &
               offsprng_lb_mutn_count(lb,hap_id,3) + 1
         endif
         
      end if
      
   end if
   
end do

end subroutine offspring

subroutine random_mate(id,available,pop,psa)
use random_pkg
integer, intent(out) :: id
integer, intent(in)  :: pop, psa
logical, intent(inout), dimension(:) :: available(psa)
id = min(pop, 1 + int(pop*randomnum(1)))
! if id selected is not available, then search for the next available id
do while(.not.available(id))
 id = mod(id, pop) + 1
end do
available(id) = .false.
end subroutine random_mate

subroutine nonrandom_mate(id,available,pop,psa)
use random_pkg
use inputs
integer, save :: y = 1
integer, intent(out) :: id
integer, intent(in)  :: pop, psa
logical, intent(inout), dimension(:) :: available(psa)
id = y
! if id selected is not available, then search for the next available id
do while(.not.available(id))
 id = mod(id, pop) + 1
end do
available(id) = .false.
end subroutine nonrandom_mate

! Thoughts on non-random mating
! 1. Simplest case "shoreline model" - each individual mates with the neighbor to its right
! 2. Wright-Fisher reproduction model
! 3. Tournament model - population is arranged as nxm grid, where each individual
!    mates with the closest neighbor with the highest fitness

subroutine wright_fisher(pop_size,off1)
!ref:
!http://www.stat.berkeley.edu/users/terry/Classes/s260.1998/Week13a/week13a/node9.html
! This subroutine is not yet finished
  use sort_module
  implicit none

  integer, intent(in)  :: pop_size
  integer, intent(out) :: off1(pop_size)
  
  integer :: i, j, n, flag(pop_size), nmax
  integer :: factorial
  real*8  :: p(pop_size,pop_size), binomial_coefficient, sorted(pop_size)
  real*8  :: success_probability, failure_probability
  
  !set flags
  flag = 0
  n = pop_size
  
  ! Pairing according probability value
  
  ! Calculate probability value
  
  do i = 1, pop_size
     do j = i+1, pop_size
        ! binomial = number of combinations of 2n objects taken j at a time
        binomial_coefficient = 2*factorial(n) &
                               /real(factorial(j)*factorial(2*n-j))
        ! j = number of success
        success_probability = (i/real(2*n))**j
        ! 2n-j = number of failures
        failure_probability = (1-i/real(2*n))**(2*n-j)
        p(i,j) = binomial_coefficient*success_probability*failure_probability
        !print '(2i,3f12.8,f15.10)',i,j,binomial_coefficient, &
        !   success_probability, failure_probability,p(i,j)
        sorted(j) = p(i,j)
     enddo
     call heapsort(sorted, pop_size)
     print *, sorted
  enddo

  off1 = 0 ! just to get it to compile without error

end subroutine wright_fisher

recursive function factorial(n) result(nfact)
implicit none
integer, intent(in) :: n
integer :: nfact
if(n > 0) then
   nfact = n * factorial(n-1)
else
   nfact = 1
end if
end function factorial
