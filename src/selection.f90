module selection_module

real*8 :: pre_sel_fitness, post_sel_fitness,                        &
          pre_sel_geno_sd, pre_sel_pheno_sd, pre_sel_corr,          &
          post_sel_geno_sd, post_sel_pheno_sd, post_sel_corr

contains

subroutine selection(dmutn, nmutn, fmutn, lb_mutn_count,   &
           linkage_block_fitness, fitness, pheno_fitness, work_fitness, &
           sorted_score, initial_allele_effects, max_size, &
           total_offspring, gen, lb_modulo, current_pop_size)

! This routine eliminates the least fit individuals in a new
! generation to reduce the population size to a level not to exceed
! pop_size.  If the population is recovering from a bottlenecking
! event, let half the excess reproduction be used to increase 
! population size and the other half be used for selection.

use random_pkg
use sort_module
use inputs
implicit none

integer, intent(inout) :: dmutn(max_del_mutn_per_indiv/2,2,*)
integer, intent(inout) :: nmutn(max_neu_mutn_per_indiv/2,2,*) 
integer, intent(inout) :: fmutn(max_fav_mutn_per_indiv/2,2,*) 
integer, intent(inout) :: lb_mutn_count(num_linkage_subunits,2,3,*)
real*8,  intent(inout) :: linkage_block_fitness(num_linkage_subunits,2,*)
real*8  :: fitness(*), pheno_fitness(*), work_fitness(*)
real*8  :: decode_fitness_del, decode_fitness_fav
real*8  :: sorted_score(max_size)
real    :: initial_allele_effects(num_linkage_subunits)
integer :: max_size, total_offspring, gen, lb_modulo, current_pop_size

integer :: i, j, k, lb, mutn, m, n, zygous(num_linkage_subunits), remaining
real*8  :: fitness_norm, homozygous_fitness_loss, noise
real*8  :: homozygous_fitness_gain, fitness_loss, covariance
real*8  :: max_work_fitness, min_work_fitness, score_cutoff
real*8  :: geno_fitness_variance, pheno_fitness_variance, hetero_effect
real*8  :: mean_pheno_fitness, se_linked, se_nonlinked, x, e, sumesq
real    :: w, p, effect, factor, bonus

w = multiplicative_weighting

! If the population is recovering from a bottlenecking event,
! compute the new population size that accounts for selection
! as well as population growth.  As a place holder here, let
! half the excess reproduction be used to increase population 
! size and the other half be used for selection.

if(bottleneck_yes) then
   if(gen > bottleneck_generation + num_bottleneck_generations &
      .and. current_pop_size < pop_size) current_pop_size = &
      min(pop_size, int(1. + current_pop_size*(1. + 0.5* &
      (reproductive_rate*(1. - fraction_random_death) - 1.0))))
end if

! Compute the fitness of each member of the new generation.

fitness(1:total_offspring) = 1.d0

if(fitness_distrib_type == 0) then ! All mutations have equal effect
   homozygous_fitness_loss = uniform_fitness_effect_del
   homozygous_fitness_gain = uniform_fitness_effect_fav
else if (fitness_distrib_type == 2) then ! All mutations neutral
   homozygous_fitness_loss = 0
   homozygous_fitness_gain = 0
end if

if (polygenic_beneficials) then

   ! Each individual with the nucleotide sequence matching the target string
   ! receives a fitness effect bonus.  For haploid organisms, the bonus equals
   ! polygenic_effect.  Diploid organisms homozygous with the target string
   ! also receive a bonus equal to polygenic_effect. Diploid organisms that 
   ! are heterozygous in regard to the target string receive a bonus equal to
   ! recessive_hetero_expression*polygenic_effect if the target string is
   ! declared to be recessive (fraction_recessive = 1.) or a bonus equal to
   ! dominant_hetero_expression*polygenic_effect if the target string is
   ! taken to be dominant (fraction_recessive /= 1.)
   ! This bonus is bestowed on a generation by generation basis to each
   ! individual carrying the target string.

   if(fraction_recessive == 1.) then
      hetero_effect = recessive_hetero_expression*polygenic_effect
   else
      hetero_effect =  dominant_hetero_expression*polygenic_effect
   end if

   do i=1,total_offspring
      if(fmutn(2,1,i) > 0 .or. fmutn(2,2,i) > 0) then
         if((fmutn(2,1,i) > 0 .and. fmutn(2,2,i) > 0) .or.  &
            recombination_model == clonal) then
            fitness(i) = fitness(i) + polygenic_effect
         else
            fitness(i) = fitness(i) + hetero_effect
         end if 
      end if
   end do

else

   do i=1,total_offspring
      do lb=1,num_linkage_subunits
         fitness(i) = (fitness(i) - (1. - w)*(2.d0          &
                        - linkage_block_fitness(lb,1,i)     &
                        - linkage_block_fitness(lb,2,i)))   &
         *(1.d0 - (1.d0 - linkage_block_fitness(lb,1,i))*w) &
         *(1.d0 - (1.d0 - linkage_block_fitness(lb,2,i))*w)
      end do
   end do

endif

do i=1,total_offspring
   
   ! Apply the appropriate fitness degradation adjustment for
   ! homozygous deleterious mutations.  Skip this step for the
   ! cases of clonal reproduction and co-dominance. 
   
   if(recombination_model /= clonal .and. .not.polygenic_beneficials &
      .and. dominant_hetero_expression /= 0.5) then
      
      j = 2 
      do k=2,dmutn(1,1,i)+1
         
         do while(abs(dmutn(k,1,i)) >  abs(dmutn(j,2,i)) .and. &
              j <= dmutn(1,2,i))
            j = j + 1
         end do
         
         if(dmutn(k,1,i) == dmutn(j,2,i)) then
            
            if(dmutn(k,1,i) == num_linkage_subunits*lb_modulo + 1) &
                 write(6,*) 'ERROR: dmutn range invalid'
            
            if(fitness_distrib_type == 1) then ! Natural mutation dist
               homozygous_fitness_loss = abs(decode_fitness_del(dmutn(k,1,i)))
            end if
            
            ! Apply the proper fitness decrease associated with a
            ! homozygous mutation, giving it 100% of the nominal 
            ! mutation effect.
            
            fitness(i) = (fitness(i) - (1. - w)*homozygous_fitness_loss) &
                 *(1.d0 - w*homozygous_fitness_loss)
            if(dmutn(k,1,i) < 0)  homozygous_fitness_loss = &
                 recessive_hetero_expression*homozygous_fitness_loss
            if(dmutn(k,1,i) > 0)  homozygous_fitness_loss = &
                 dominant_hetero_expression*homozygous_fitness_loss
            ! Remove the fitness decreases that were applied elsewhere 
            ! when it was assumed the mutation was heterozygous.  
            ! Remove the heterozygous effect by adding it back twice.
            ! since it was carried out on both haplotypes.
            fitness(i) = fitness(i) / (1.d0 - w*homozygous_fitness_loss)**2 &
                 + (1. - w) *homozygous_fitness_loss*2.
         end if
      end do
      
      ! Apply the appropriate fitness enhancement adjustment for
      ! homozygous favorable mutations. 
      
      j = 2 
      do k=2,fmutn(1,1,i)+1
         
         do while(abs(fmutn(k,1,i))>abs(fmutn(j,2,i)) .and. j<=fmutn(1,2,i))
            j = j + 1
         end do
         
         if(fmutn(k,1,i) == fmutn(j,2,i)) then
            
            if(fmutn(k,1,i) == num_linkage_subunits*lb_modulo + 1) &
                 write(6,*) 'ERROR: fmutn range invalid'
            
            if(fitness_distrib_type == 1) & ! Natural mutation dist
                 homozygous_fitness_gain = decode_fitness_fav(fmutn(k,1,i))
            fitness(i) = (fitness(i)  &
                 + (1. - w)*homozygous_fitness_gain) &
                 *(1.d0 + w*homozygous_fitness_gain)
            if(fmutn(k,1,i) < 0) homozygous_fitness_gain = &
                 recessive_hetero_expression*homozygous_fitness_gain
            if(fmutn(k,1,i) > 0) homozygous_fitness_gain = &
                 dominant_hetero_expression*homozygous_fitness_gain
            fitness(i) = (fitness(i) &
                 - (1. - w) *homozygous_fitness_gain*2.) &
                 / (1.d0 + w*homozygous_fitness_gain)**2
         end if
      end do
      
   end if
end do

if(synergistic_epistasis .and. recombination_model /= clonal) then
   
   !  In our synergistic epistasis (SE) treatment, we break its
   !  effect into two parts, one involving interactions between
   !  mutations occurring on the same linkage block (linked
   !  interactions) and the other part involving interactions of
   !  mutations on different linkage blocks (nonlinked interactions).
   !  SE effects from linked interactions are inherited perfectly,
   !  while those from nonlinked interactions are progressively
   !  scrambled by recombination generation to generation.
   !
   !  Let us first consider the linked interactions. We apply the
   !  following considerations.  First, we require amplitude of
   !  the SE effect to be _proportional_ to the non-epistatic fitness
   !  effects of each of the two interacting mutations. This means
   !  that if a mutation's effect on the non-mutant genome is small,
   !  then the SE contributions from its interactions with other
   !  mutations is assumed likewise to be small. If we use f to
   !  denote linkage block fitness, then (1 - f) represents the sum
   !  of the non-epistatic fitness effects of all the mutations on
   !  the linkage block.  The sum of the products of the fitness
   !  effects of all the mutations is then given by 0.5*(1-f)**2,
   !  corrected for the self-interaction contributions.
   !  We assume co-dominance, however, that is, we assume that only
   !  50% of each mutations base value is used in computing the SE
   !  contribution. Further, we allow the user to scale the SE
   !  contribution through the parameter se_linked_scaling. These
   !  considerations then imply that the SE effect from linked
   !  mutations on a given linkage block is given by
   !
   !    0.125*se_linked_scaling*((1-f)**2 - self_int_contribution)
   !
   !  Lets now consider the nonlinked SE interactions.  If M is the
   !  total number of mutations in the genome and n is the number of
   !  linkage blocks, then the total number of pairwise interactions
   !  between mutations is M(M-1)/2, the mean number of mutations
   !  per linkage block is M/n, and the approximate number of linked
   !  interactions is n(M/n)[(M/n)-1]/2.  Since SE contributions
   !  become significant only when M becomes moderately large, we
   !  approximate M-1 by M and (M/n)-1 by M/n.  With these
   !  approximations, the number of linked interactions becomes
   !  M**2/(2n) and the number of nonlinked interactions becomes
   !  (1 - 1/n)*M**2/2.
   !
   !  Let F denote the overall genotypic fitness.  The total
   !  nonlinked SE fitness contribution is then proportional to the
   !  sum of the non-epistatic fitness effects of all the individual
   !  mutations, (1-F), but scaled to account for the portion of the
   !  mutations which are linked with the factor (1 - 1/n), times
   !  the mean non-epistatic fitness effect of these mutations,
   !  (1-F)/M, times the number of unique pair-wise interactions,
   !  (1 - 1/n)M/2, that each non-linked mutation has with the
   !  others.  We assume co-dominance, which implies each haploid
   !  occurrence of a mutation gives 50% expression of the mutations
   !  non-epistatic value which reduces the overall contribution by
   !  a factor of 0.25.  We also subtract the self-interaction
   !  contribution implicit in the 0.5*(1-F)**2 formula.
   !
   !  We scale this non-linked SE contribution by the user-specified
   !  input parameter se_nonlinked_scaling.  In general one expects
   !  that interaction between mutations within the same linkage
   !  block will on average have much greater SE effects than
   !  mutations which are more distant to one another within the
   !  genome.  Hence, a value for this parameter much less than (say,
   !  by a factor of 0.001 times) the parameter se_linked_scaling
   !  used for the linked mutations is usually appropriate.  The
   !  resulting expression for the non-linked SE contribution to
   !  individual fitness, to be subtracted from F, is
   !
   !                   0.125*se_nonlinked_scaling
   !             *((1-F)**2 - self_int_term)*(1 - 1/n)**2.
   
   do i=1,total_offspring
      
      !     Compute self-interaction term, neglecting the non-tracked
      !     mutations because of their small values.  Assume
      !     co-dominance.
      
      sumesq = 0.d0
      do k=1,2
         do j=2,dmutn(1,k,i)+1
            e = 0.5*abs(decode_fitness_del(dmutn(j,k,i)))
            sumesq = sumesq + e**2
         end do
      end do
      
      ! Sum the linked SE contributions from each of the linkage blocks.
      
      se_linked = 0.
      if(recombination_model /= clonal) then
         do lb=1,num_linkage_subunits
            se_linked = se_linked &
                 + (2.d0 - linkage_block_fitness(lb,1,i) &
                 - linkage_block_fitness(lb,2,i))**2
         end do
      else
         do lb=1,num_linkage_subunits
            se_linked = se_linked &
                 + (1.d0 - linkage_block_fitness(lb,1,i))**2 &
                 + (1.d0 - linkage_block_fitness(lb,2,i))**2
         end do
      end if
      
      ! Subtract the self-interaction sum from the se_linked sum
      ! and scale the remainder appropriately.
      
      se_linked = 0.125*se_linked_scaling &
           * max(0., se_linked - sumesq)
      
      ! Compute the non-linked SE contribution, subtract the
      ! self-interaction sum, and scale.
      
      se_nonlinked = 0.125*se_nonlinked_scaling &
           * max(0., ((1.d0 - fitness(i))**2 - sumesq))
      
      ! If linked SE is being included, add the appropriate factor
      ! to account for the number of linked interactions.  Note that
      ! for diploid organisms the total number of linkage blocks is
      ! 2*num_linkage_subunits.
      
      if(se_linked_scaling > 0.) then
         se_nonlinked = se_nonlinked*(1. - 0.5/num_linkage_subunits)**2
      end if
      if(recombination_model == clonal) se_nonlinked = 0.
      
      fitness(i) = fitness(i) - se_linked - se_nonlinked
      
   end do
   
end if

! Account for possible homozygosity in initial contrasting alleles.

if(num_contrasting_alleles > 0) then
   
   do i=1,total_offspring
      
      zygous = 0

      do m=2,dmutn(1,1,i)+1
         if(mod(dmutn(m,1,i), lb_modulo) == lb_modulo-1) then 
            lb = dmutn(m,1,i)/lb_modulo + 1
            zygous(lb) = zygous(lb) + 1
         end if
      end do
      
      do m=2,dmutn(1,2,i)+1
         if(mod(dmutn(m,2,i), lb_modulo) == lb_modulo-1) then 
            lb = dmutn(m,2,i)/lb_modulo + 1
            zygous(lb) = zygous(lb) + 1
         end if
      end do
      
      do lb=1,num_linkage_subunits
         if(zygous(lb) == 2) then
            effect = initial_allele_effects(lb)
            fitness(i) = (fitness(i) - (1. - w)*effect) &
                 *(1.d0 - w*effect)
            effect = recessive_hetero_expression*effect
            fitness(i) = fitness(i) / (1.d0 - w*effect)**2 &
                 + (1. - w) *effect*2.
         end if
      end do
      
      zygous = 0
      
      do m=2,fmutn(1,1,i)+1
         if(mod(fmutn(m,1,i), lb_modulo) == lb_modulo-1) then 
            lb = fmutn(m,1,i)/lb_modulo + 1
            zygous(lb) = zygous(lb) + 1
         end if
      end do
      
      do m=2,fmutn(1,2,i)+1
         if(mod(fmutn(m,2,i), lb_modulo) == lb_modulo-1) then 
            lb = fmutn(m,2,i)/lb_modulo + 1
            zygous(lb) = zygous(lb) + 1
         end if
      end do
      
      do lb=1,num_linkage_subunits
         if(zygous(lb) == 2) then
            effect = initial_allele_effects(lb)
            fitness(i) = (fitness(i) + (1. - w)*effect) &
                                     *(1.d0 + w*effect)
            effect =  dominant_hetero_expression*effect
            fitness(i) = (fitness(i) - (1. - w) *effect*2.) &
                                     / (1.d0 + w*effect)**2
         end if
      end do

   end do

end if

! Compute the mean genotypic fitness of the new generation.

pre_sel_fitness = 0.d0

do i=1,total_offspring
   pre_sel_fitness = pre_sel_fitness + fitness(i)
end do

pre_sel_fitness = pre_sel_fitness/total_offspring

! Compute the genotypic fitness variance of the new generation.

geno_fitness_variance = 0.d0

do i=1,total_offspring
   geno_fitness_variance = geno_fitness_variance &
                    + (fitness(i) - pre_sel_fitness)**2
end do

geno_fitness_variance = geno_fitness_variance/total_offspring

pre_sel_geno_sd = sqrt(geno_fitness_variance)

! If population has collapsed to a single individual, skip the
! selection process and return.

if(total_offspring == 1) then
   current_pop_size = 1
   return
end if

! Compute the noise variance required to yield the specified
! heritability.  Add to this fitness-dependent noise a noise 
! component that is fitness independent. Take the square root
! to obtain the standard deviation. 

noise = sqrt(geno_fitness_variance*(1. - heritability) &
             /heritability + non_scaling_noise**2) 

! Add noise to the fitness to create a phenotypic fitness score.
! Add a tiny variable positive increment to eliminate identical 
! fitness values when the noise is zero.

do i=1,total_offspring
   pheno_fitness(i) = fitness(i) + random_normal()*noise + 1.d-15*i
end do

! Compute the mean phenotypic fitness of offspring.

mean_pheno_fitness = 0.d0

do i=1,total_offspring
   mean_pheno_fitness = mean_pheno_fitness + pheno_fitness(i)
end do

mean_pheno_fitness = mean_pheno_fitness/total_offspring

! Compute the phenotypic fitness variance, the covariance of
! genotypic and phenotypic fitness, and the genotype-phenotype
! correlation.

pheno_fitness_variance = 0.d0
covariance = 0.d0

do i=1,total_offspring
   pheno_fitness_variance = pheno_fitness_variance &
        + (pheno_fitness(i) - mean_pheno_fitness)**2
   covariance = covariance + fitness(i)*pheno_fitness(i)
end do

pheno_fitness_variance = pheno_fitness_variance/total_offspring

pre_sel_pheno_sd = sqrt(pheno_fitness_variance)

covariance = covariance/total_offspring &
           - pre_sel_fitness*mean_pheno_fitness

pre_sel_corr = 0.
effect = sqrt(geno_fitness_variance*pheno_fitness_variance)
if(effect > 0.) pre_sel_corr = covariance/effect

! Move, in effect, those offspring whose phenotypic fitness is 
! negative to the end of the list of offspring, and then, in effect,
! truncate the list so that these individuals cannot reproduce and
! do not even participate in the subsequent selection process.

remaining = total_offspring
do i=1,total_offspring
   if(pheno_fitness(i) < 0.d0) then
      do while(pheno_fitness(remaining) < 0.d0 .and. remaining > 1)
         remaining = remaining - 1
      end do
      k = dmutn(1,1,remaining) + 1
      dmutn(1:k,1,i) = dmutn(1:k,1,remaining)
      k = dmutn(1,2,remaining) + 1
      dmutn(1:k,2,i) = dmutn(1:k,2,remaining)
      k = nmutn(1,1,remaining) + 1
      nmutn(1:k,1,i) = nmutn(1:k,1,remaining)
      k = nmutn(1,2,remaining) + 1
      nmutn(1:k,2,i) = nmutn(1:k,2,remaining)
      k = fmutn(1,1,remaining) + 1
      fmutn(1:k,1,i) = fmutn(1:k,1,remaining)
      k = fmutn(1,2,remaining) + 1
      fmutn(1:k,2,i) = fmutn(1:k,2,remaining)
      lb_mutn_count(:,:,:,i) = lb_mutn_count(:,:,:,remaining)
      linkage_block_fitness(:,:,i) = linkage_block_fitness(:,:,remaining)
      fitness(i) = fitness(remaining)
      pheno_fitness(i) = pheno_fitness(remaining)
      if(remaining > 1) remaining = remaining - 1
   end if
end do

total_offspring = remaining

! Adjust the population size for the next generation such that it
! does not exceed the number of offspring after removal of those
! with negative phenotypic fitness.

current_pop_size = min(current_pop_size, remaining)

! Allow the population size for the next generation potentially 
! to rebound from an earlier reduction in previous generations 
! because of individuals with negative phenotypic fitness.

if(.not.tribal_competition .and. .not.bottleneck_yes .and. &
   remaining > current_pop_size) &
   current_pop_size = min(remaining, pop_size)

! Copy the phenotypic fitnesses into array work_fitness.

work_fitness(1:total_offspring) = pheno_fitness(1:total_offspring)

if (selection_scheme == 2) then

!  For unrestricted probability selection, divide the phenotypic  
!  fitness by a uniformly distributed random number prior to 
!  ranking and truncation.  This procedure allows the probability 
!  of surviving and reproducing in the next generation to be
!  directly related to phenotypic fitness and also for the correct
!  number of individuals to be eliminated to maintain a constant
!  population size.

   do i=1,total_offspring
      work_fitness(i) = work_fitness(i)/(randomnum(1) + 1.d-15)
   end do
   
end if

if (selection_scheme == 3) then

!  For strict proportionality probability selection, rescale the
!  phenotypic fitness values such that the maximum value is one.
!  Then divide the scaled phenotypic fitness by a uniformly
!  distributed random number prior to ranking and truncation.  
!  Allow only those individuals to reproduce whose resulting 
!  ratio of scaled phenotypic fitness to the random number value
!  exceeds one.  This approach ensures that no individual 
!  automatically survives to reproduce regardless of the value
!  of the random number.  But it restricts the fraction of the 
!  offspring that can survive.  Therefore, when the reproduction
!  rate is low, the number of surviving offspring may not be
!  large enough to sustain a constant population size.

   max_work_fitness = 0.d0
   do i=1,total_offspring
      max_work_fitness = max(max_work_fitness, work_fitness(i))
   end do
  
   do i=1,total_offspring
      work_fitness(i) = work_fitness(i)/(max_work_fitness + 1.d-15)
      work_fitness(i) = work_fitness(i)/(randomnum(1) + 1.d-15)
   end do
   
end if

if (selection_scheme == 4) then

!  For partial truncation selection, divide the phenotypic  
!  fitness by the sum of theta and (1. - theta) times a random 
!  number distributed uniformly between 0.0 and 1.0 prior to 
!  ranking and truncation, where theta is the parameter 
!  partial_truncation_value.  This selection scheme is 
!  intermediate between truncation selection and unrestricted 
!  probability selection.  The procedure allows for the correct 
!  number of individuals to be eliminated to maintain a constant 
!  population size.

   do i=1,total_offspring
      work_fitness(i) = work_fitness(i)/(partial_truncation_value &
                        + (1. - partial_truncation_value)*randomnum(1))
   end do
   
end if

! Sort the resulting work fitnesses in ascending order.

sorted_score(1:total_offspring) = work_fitness(1:total_offspring)

if (total_offspring > 1) call heapsort(sorted_score, total_offspring)

if (selection_scheme <= 4) then

!  Apply truncation selection to reduce the population size to
!  current_pop_size.

!  Compute the score cutoff value.

   if(total_offspring > current_pop_size) then
      score_cutoff = sorted_score(total_offspring - current_pop_size)
   else
      score_cutoff = -1000.d0
   end if

   if(selection_scheme == 3) score_cutoff = max(1.d0, score_cutoff)

!  Copy pheno_fitness into array sorted_score for diagnostics
!  purposes.

   sorted_score(1:total_offspring) =  &
                              pheno_fitness(1:total_offspring)

!  Remove those individuals whose score lies below the cutoff
!  value to reduce the population size to its appropriate value.

   current_pop_size = min(current_pop_size, total_offspring)
   remaining = total_offspring

   do i=1,current_pop_size

!     If the work fitness if individual i is below the cutoff 
!     value, find another individual in the pool of excess 
!     offspring whose work fitness is equal to or above the
!     cutoff value and replace the first individual with the 
!     second in the list of reproducing individuals for that
!     generation.

      if(work_fitness(i) < score_cutoff .and. i < remaining) then
         do while(work_fitness(remaining) < score_cutoff .and. &
                  remaining > 1)
            remaining = remaining - 1
         end do
         if(i < remaining) then
            k = dmutn(1,1,remaining) + 1
            dmutn(1:k,1,i) = dmutn(1:k,1,remaining)
            k = dmutn(1,2,remaining) + 1
            dmutn(1:k,2,i) = dmutn(1:k,2,remaining)
            k = nmutn(1,1,remaining) + 1
            nmutn(1:k,1,i) = nmutn(1:k,1,remaining)
            k = nmutn(1,2,remaining) + 1
            nmutn(1:k,2,i) = nmutn(1:k,2,remaining)
            k = fmutn(1,1,remaining) + 1
            fmutn(1:k,1,i) = fmutn(1:k,1,remaining)
            k = fmutn(1,2,remaining) + 1
            fmutn(1:k,2,i) = fmutn(1:k,2,remaining)
            lb_mutn_count(:,:,:,i) = lb_mutn_count(:,:,:,remaining)
            linkage_block_fitness(:,:,i) = &
            linkage_block_fitness(:,:,remaining)
            fitness(i) = fitness(remaining)
            pheno_fitness(i) = pheno_fitness(remaining)
            if(remaining > 1) remaining = remaining - 1
         end if
      end if
   end do

   current_pop_size = min(current_pop_size, remaining)

else

   write(6,*) 'ERROR: invalid selection scheme', selection_scheme
   write(9,*) 'ERROR: invalid selection scheme', selection_scheme
   stop

end if

! Compute the mean genotypic and phenotypic fitnesses of the new
! generation after selection.

post_sel_fitness   = 0.d0
mean_pheno_fitness = 0.d0

do i=1,current_pop_size
   post_sel_fitness   = post_sel_fitness   + fitness(i)
   mean_pheno_fitness = mean_pheno_fitness + pheno_fitness(i)
end do

post_sel_fitness   = post_sel_fitness/current_pop_size
mean_pheno_fitness = mean_pheno_fitness/current_pop_size

! Compute the genotypic and phenotypic fitness variances, the 
! covariance of genotypic and phenotypic fitness, and the 
! genotype-phenotype correlation of the new generation.

geno_fitness_variance  = 0.d0
pheno_fitness_variance = 0.d0
covariance = 0.d0

do i=1,current_pop_size
   geno_fitness_variance  = geno_fitness_variance &
                    + (fitness(i) - post_sel_fitness)**2
   pheno_fitness_variance = pheno_fitness_variance &
        + (pheno_fitness(i) - mean_pheno_fitness)**2
   covariance = covariance + fitness(i)*pheno_fitness(i)
end do

geno_fitness_variance  = geno_fitness_variance/current_pop_size
pheno_fitness_variance = pheno_fitness_variance/current_pop_size

post_sel_geno_sd  = sqrt(geno_fitness_variance)
post_sel_pheno_sd = sqrt(pheno_fitness_variance)

covariance = covariance/current_pop_size &
           - post_sel_fitness*mean_pheno_fitness

post_sel_corr = 0.
effect = sqrt(geno_fitness_variance*pheno_fitness_variance)
if(effect > 0.) post_sel_corr = covariance/effect

post_sel_fitness = max(0., post_sel_fitness)

fitness(current_pop_size+1:max_size) = 0.d0

end subroutine selection

subroutine selection2(fitness, pheno_fitness, work_fitness, sorted_score, &
                      initial_allele_effects, max_size,  &
                      total_offspring, gen, lb_modulo, current_pop_size)

! Selection2 is the same as selection, but now passing the
! genome data (dmutn, fmutn, etc.) via the genome module as pointers

! This routine eliminates the least fit individuals in a new
! generation to reduce the population size to a level not to exceed
! pop_size.  If the population is recovering from a bottlenecking
! event, let half the excess reproduction be used to increase 
! population size and the other half be used for selection.

use genome
use random_pkg
use sort_module
use inputs
use init
implicit none

!integer, intent(inout) :: dmutn(max_del_mutn_per_indiv/2,2,*)
!integer, intent(inout) :: fmutn(max_fav_mutn_per_indiv/2,2,*) 
!integer, intent(inout) :: lb_mutn_count(num_linkage_subunits,2,3,*)
!real*8,  intent(inout) :: linkage_block_fitness(num_linkage_subunits,2,*)
real*8  :: fitness(*), pheno_fitness(*), work_fitness(*)
real*8  :: decode_fitness_del, decode_fitness_fav
real*8  :: sorted_score(max_size)
real    :: initial_allele_effects(num_linkage_subunits)
integer :: max_size, total_offspring, gen, lb_modulo, current_pop_size

integer :: i, j, k, lb, mutn, m, n, zygous(num_linkage_subunits), remaining
real*8  :: fitness_norm, homozygous_fitness_loss, noise
real*8  :: homozygous_fitness_gain, fitness_loss, covariance
real*8  :: max_work_fitness, min_work_fitness, score_cutoff
real*8  :: geno_fitness_variance, pheno_fitness_variance
real*8  :: mean_pheno_fitness, se_linked, se_nonlinked, x, e, sumesq
real    :: w, p, effect, factor

w = multiplicative_weighting

! If the population is recovering from a bottlenecking event,
! compute the new population size that accounts for selection
! as well as population growth.  As a place holder here, let
! half the excess reproduction be used to increase population 
! size and the other half be used for selection.

if(bottleneck_yes) then
   if(gen > bottleneck_generation + num_bottleneck_generations &
      .and. current_pop_size < pop_size) current_pop_size = &
      min(pop_size, int(1. + current_pop_size*(1. + 0.5* &
      (reproductive_rate*(1. - fraction_random_death) - 1.0))))
end if

! Compute the fitness of each member of the new generation.

fitness(1:total_offspring) = 1.d0

if(fitness_distrib_type == 0) then ! All mutations have equal effect
   homozygous_fitness_loss = uniform_fitness_effect_del
   homozygous_fitness_gain = uniform_fitness_effect_fav
else if (fitness_distrib_type == 2) then ! All mutations neutral
   homozygous_fitness_loss = 0
   homozygous_fitness_gain = 0
end if

do i=1,total_offspring

   do lb=1,num_linkage_subunits
      fitness(i) = (fitness(i) - (1. - w)*(2.d0          &
                     - gp(i)%lbf(lb,1)     &
                     - gp(i)%lbf(lb,2)))   &
      *(1.d0 - (1.d0 - gp(i)%lbf(lb,1))*w) &
      *(1.d0 - (1.d0 - gp(i)%lbf(lb,2))*w)
   end do

!  Apply the appropriate fitness degradation adjustment for
!  homozygous deleterious mutations.  Skip this step for the
!  cases of clonal reproduction and co-dominance. 

   if(recombination_model /= clonal .and. &
        dominant_hetero_expression /= 0.5) then
      
      j = 2 
      do k=2,gp(i)%dm(1,1)+1
         
         do while(abs(gp(i)%dm(k,1)) >  abs(gp(i)%dm(j,2)) .and. &
              j <= gp(i)%dm(1,2))
            j = j + 1
         end do
         
         if(gp(i)%dm(k,1) == gp(i)%dm(j,2)) then
            
            if(gp(i)%dm(k,1) == num_linkage_subunits*lb_modulo + 1) &
                 write(6,*) 'ERROR: gp(i)%dm range invalid'
            
            if(fitness_distrib_type == 1) then ! Natural mutation dist
               homozygous_fitness_loss = abs(decode_fitness_del(gp(i)%dm(k,1)))
               if(x >= 1.d0) homozygous_fitness_loss = 0.d0
            end if

            ! Apply the proper fitness decrease associated with a
            ! homozygous mutation, giving it 100% of the nominal 
            ! mutation effect.
            
            fitness(i) = (fitness(i) - (1. - w)*homozygous_fitness_loss) &
                        *(1.d0 - w*homozygous_fitness_loss)
            if(gp(i)%dm(k,1) < 0)  homozygous_fitness_loss = &
                 recessive_hetero_expression*homozygous_fitness_loss
            if(gp(i)%dm(k,1) > 0)  homozygous_fitness_loss = &
                 dominant_hetero_expression*homozygous_fitness_loss
            ! Remove the fitness decreases that were applied elsewhere 
            ! when it was assumed the mutation was heterozygous.  
            ! Remove the heterozygous effect by adding it back twice.
            ! since it was carried out on both haplotypes.
            fitness(i) = fitness(i) / (1.d0 - w*homozygous_fitness_loss)**2 &
                 + (1. - w) *homozygous_fitness_loss*2.
         end if
      end do
      
      ! Apply the appropriate fitness enhancement adjustment for
      ! homozygous favorable mutations. 
      
      j = 2 
      do k=2,gp(i)%fm(1,1)+1
         
         do while(abs(gp(i)%fm(k,1))>abs(gp(i)%fm(j,2)) .and. j<=gp(i)%fm(1,2))
            j = j + 1
         end do
         
         if(gp(i)%fm(k,1) == gp(i)%fm(j,2)) then
            
            if(gp(i)%fm(k,1) == num_linkage_subunits*lb_modulo + 1) &
                 write(6,*) 'ERROR: gp(i)%fm range invalid'
            
            if(fitness_distrib_type == 1) & ! Natural mutation dist
                 homozygous_fitness_gain = decode_fitness_fav(gp(i)%fm(k,1))
            fitness(i) = (fitness(i)  &
                 + (1. - w)*homozygous_fitness_gain) &
                 *(1.d0 + w*homozygous_fitness_gain)
            if(gp(i)%fm(k,1) < 0) homozygous_fitness_gain = &
                 recessive_hetero_expression &
                 *homozygous_fitness_gain
            if(gp(i)%fm(k,1) > 0) homozygous_fitness_gain = &
                 dominant_hetero_expression &
                 *homozygous_fitness_gain
            fitness(i) = (fitness(i) &
                 - (1. - w) *homozygous_fitness_gain*2.) &
                 / (1.d0 + w*homozygous_fitness_gain)**2
         end if
      end do
      
   end if
   
end do

if(synergistic_epistasis .and. recombination_model /= clonal) then

!  In our synergistic epistasis (SE) treatment, we break its
!  effect into two parts, one involving interactions between
!  mutations occurring on the same linkage block (linked
!  interactions) and the other part involving interactions of
!  mutations on different linkage blocks (nonlinked interactions).
!  SE effects from linked interactions are inherited perfectly,
!  while those from nonlinked interactions are progressively
!  scrambled by recombination generation to generation.
!
!  Let us first consider the linked interactions. We apply the
!  following considerations.  First, we require amplitude of
!  the SE effect to be _proportional_ to the non-epistatic fitness
!  effects of each of the two interacting mutations. This means
!  that if a mutation's effect on the non-mutant genome is small,
!  then the SE contributions from its interactions with other
!  mutations is assumed likewise to be small. If we use f to
!  denote linkage block fitness, then (1 - f) represents the sum
!  of the non-epistatic fitness effects of all the mutations on
!  the linkage block.  The sum of the products of the fitness
!  effects of all the mutations is then given by 0.5*(1-f)**2,
!  corrected for the self-interaction contributions.
!  We assume co-dominance, however, that is, we assume that only
!  50% of each mutations base value is used in computing the SE
!  contribution. Further, we allow the user to scale the SE
!  contribution through the parameter se_linked_scaling. These
!  considerations then imply that the SE effect from linked
!  mutations on a given linkage block is given by
!
!    0.125*se_linked_scaling*((1-f)**2 - self_int_contribution)
!
!  Lets now consider the nonlinked SE interactions.  If M is the
!  total number of mutations in the genome and n is the number of
!  linkage blocks, then the total number of pairwise interactions
!  between mutations is M(M-1)/2, the mean number of mutations
!  per linkage block is M/n, and the approximate number of linked
!  interactions is n(M/n)[(M/n)-1]/2.  Since SE contributions
!  become significant only when M becomes moderately large, we
!  approximate M-1 by M and (M/n)-1 by M/n.  With these
!  approximations, the number of linked interactions becomes
!  M**2/(2n) and the number of nonlinked interactions becomes
!  (1 - 1/n)*M**2/2.
!
!  Let F denote the overall genotypic fitness.  The total
!  nonlinked SE fitness contribution is then proportional to the
!  sum of the non-epistatic fitness effects of all the individual
!  mutations, (1-F), but scaled to account for the portion of the
!  mutations which are linked with the factor (1 - 1/n), times
!  the mean non-epistatic fitness effect of these mutations,
!  (1-F)/M, times the number of unique pair-wise interactions,
!  (1 - 1/n)M/2, that each non-linked mutation has with the
!  others.  We assume co-dominance, which implies each haploid
!  occurrence of a mutation gives 50% expression of the mutations
!  non-epistatic value which reduces the overall contribution by
!  a factor of 0.25.  We also subtract the self-interaction
!  contribution implicit in the 0.5*(1-F)**2 formula.
!
!  We scale this non-linked SE contribution by the user-specified
!  input parameter se_nonlinked_scaling.  In general one expects
!  that interaction between mutations within the same linkage
!  block will on average have much greater SE effects than
!  mutations which are more distant to one another within the
!  genome.  Hence, a value for this parameter much less than (say,
!  by a factor of 0.001 times) the parameter se_linked_scaling
!  used for the linked mutations is usually appropriate.  The
!  resulting expression for the non-linked SE contribution to
!  individual fitness, to be subtracted from F, is
!
!                   0.125*se_nonlinked_scaling
!             *((1-F)**2 - self_int_term)*(1 - 1/n)**2.

   do i=1,total_offspring

!     Compute self-interaction term, neglecting the non-tracked
!     mutations because of their small values.  Assume
!     co-dominance.

      sumesq = 0.d0
      do k=1,2
         do j=2,gp(i)%dm(1,k)+1
            e = 0.5*abs(decode_fitness_del(gp(i)%dm(j,k)))
            sumesq = sumesq + e**2
         end do
      end do
   
!     Sum the linked SE contributions from each of the linkage blocks.

      se_linked = 0.
      if(recombination_model == clonal) then
         do lb=1,num_linkage_subunits
            se_linked = se_linked &
                      + (2.d0 - gp(i)%lbf(lb,1) &
                              - gp(i)%lbf(lb,2))**2
         end do
      else
         do lb=1,num_linkage_subunits
            se_linked = se_linked &
                      + (1.d0 - gp(i)%lbf(lb,1))**2 &
                      + (1.d0 - gp(i)%lbf(lb,2))**2
         end do
      end if

!     Subtract the self-interaction sum from the se_linked sum
!     and scale the remainder appropriately.

      se_linked = 0.125*se_linked_scaling &
                * max(0., se_linked - sumesq)

!     Compute the non-linked SE contribution, subtract the
!     self-interaction sum, and scale.

      se_nonlinked = 0.125*se_nonlinked_scaling &
                   * max(0., ((1.d0 - fitness(i))**2 - sumesq))

!     If linked SE is being included, add the appropriate factor
!     to account for the number of linked interactions.  Note that
!     for diploid organisms the total number of linkage blocks is
!     2*num_linkage_subunits.

      if(se_linked_scaling > 0.) then
         se_nonlinked = se_nonlinked*(1. - 0.5/num_linkage_subunits)**2
      end if
      if(recombination_model == clonal) se_nonlinked = 0.

      fitness(i) = fitness(i) - se_linked - se_nonlinked

   end do

end if

! Account for possible homozygosity in initial contrasting alleles.

if(num_contrasting_alleles > 0) then

   do i=1,total_offspring

      zygous = 0

      do m=2,gp(i)%dm(1,1)+1
         if(mod(gp(i)%dm(m,1), lb_modulo) == lb_modulo-1) then 
            lb = gp(i)%dm(m,1)/lb_modulo + 1
            zygous(lb) = zygous(lb) + 1
         end if
      end do

      do m=2,gp(i)%dm(1,2)+1
         if(mod(gp(i)%dm(m,2), lb_modulo) == lb_modulo-1) then 
            lb = gp(i)%dm(m,2)/lb_modulo + 1
            zygous(lb) = zygous(lb) + 1
         end if
      end do

      do lb=1,num_linkage_subunits
         if(zygous(lb) == 2) then
            effect = initial_allele_effects(lb)
            fitness(i) = (fitness(i) - (1. - w)*effect) &
                                     *(1.d0 - w*effect)
            effect = recessive_hetero_expression*effect
            fitness(i) = fitness(i) / (1.d0 - w*effect)**2 &
                                    + (1. - w) *effect*2.
         end if
      end do

      zygous = 0

      do m=2,gp(i)%fm(1,1)+1
         if(mod(gp(i)%fm(m,1), lb_modulo) == lb_modulo-1) then 
            lb = gp(i)%fm(m,1)/lb_modulo + 1
            zygous(lb) = zygous(lb) + 1
         end if
      end do

      do m=2,gp(i)%fm(1,2)+1
         if(mod(gp(i)%fm(m,2), lb_modulo) == lb_modulo-1) then 
            lb = gp(i)%fm(m,2)/lb_modulo + 1
            zygous(lb) = zygous(lb) + 1
         end if
      end do

      do lb=1,num_linkage_subunits
         if(zygous(lb) == 2) then
            effect = initial_allele_effects(lb)
            fitness(i) = (fitness(i) + (1. - w)*effect) &
                                     *(1.d0 + w*effect)
            effect =  dominant_hetero_expression*effect
            fitness(i) = (fitness(i) - (1. - w) *effect*2.) &
                                     / (1.d0 + w*effect)**2
         end if
      end do

   end do

end if

! Compute the mean genotypic fitness of the new generation.

pre_sel_fitness = 0.d0

do i=1,total_offspring
   pre_sel_fitness = pre_sel_fitness + fitness(i)
end do

pre_sel_fitness = pre_sel_fitness/total_offspring

! Compute the genotypic fitness variance of the new generation.

geno_fitness_variance = 0.d0

do i=1,total_offspring
   geno_fitness_variance = geno_fitness_variance &
                    + (fitness(i) - pre_sel_fitness)**2
end do

geno_fitness_variance = geno_fitness_variance/total_offspring

pre_sel_geno_sd = sqrt(geno_fitness_variance)

! If population has collapsed to a single individual, skip the
! selection process and return.

if(total_offspring == 1) then
   current_pop_size = 1
   return
end if

! Compute the noise variance required to yield the specified
! heritability.  Add to this fitness-dependent noise a noise 
! component that is fitness independent. Take the square root
! to obtain the standard deviation. 

noise = sqrt(geno_fitness_variance*(1. - heritability) &
             /heritability + non_scaling_noise**2) 

! Add noise to the fitness to create a phenotypic fitness score.
! Add a tiny variable positive increment to eliminate identical 
! fitness values when the noise is zero.

do i=1,total_offspring
   pheno_fitness(i) = fitness(i) + random_normal()*noise + 1.d-15*i
end do

! Compute the mean phenotypic fitness of offspring.

mean_pheno_fitness = 0.d0

do i=1,total_offspring
   mean_pheno_fitness = mean_pheno_fitness + pheno_fitness(i)
end do

mean_pheno_fitness = mean_pheno_fitness/total_offspring

! Compute the phenotypic fitness variance, the covariance of
! genotypic and phenotypic fitness, and the genotype-phenotype
! correlation.

pheno_fitness_variance = 0.d0
covariance = 0.d0

do i=1,total_offspring
   pheno_fitness_variance = pheno_fitness_variance &
        + (pheno_fitness(i) - mean_pheno_fitness)**2
   covariance = covariance + fitness(i)*pheno_fitness(i)
end do

pheno_fitness_variance = pheno_fitness_variance/total_offspring

pre_sel_pheno_sd = sqrt(pheno_fitness_variance)

covariance = covariance/total_offspring &
           - pre_sel_fitness*mean_pheno_fitness

pre_sel_corr = 0.
effect = sqrt(geno_fitness_variance*pheno_fitness_variance)
if(effect > 0.) pre_sel_corr = covariance/effect

! Move, in effect, those offspring whose phenotypic fitness is 
! negative to the end of the list of offspring, and then, in effect,
! truncate the list so that these individuals cannot reproduce and
! do not even participate in the subsequent selection process.

remaining = total_offspring
do i=1,total_offspring
   if(pheno_fitness(i) < 0.d0) then
      do while(pheno_fitness(remaining) < 0.d0 .and. remaining > 1)
         remaining = remaining - 1
      end do
      k = gp(remaining)%dm(1,1) + 1
      gp(i)%dm(1:k,1) = gp(remaining)%dm(1:k,1)
      k = gp(remaining)%dm(1,2) + 1
      gp(i)%dm(1:k,2) = gp(remaining)%dm(1:k,2)
      k = gp(remaining)%fm(1,1) + 1
      gp(i)%fm(1:k,1) = gp(remaining)%fm(1:k,1)
      k = gp(remaining)%fm(1,2) + 1
      gp(i)%fm(1:k,2) = gp(remaining)%fm(1:k,2)
      gp(i)%lbmc(:,:,:) = gp(remaining)%lbmc(:,:,:)
      gp(i)%lbf(:,:) = gp(remaining)%lbf(:,:)
      fitness(i) = fitness(remaining)
      pheno_fitness(i) = pheno_fitness(remaining)
      if(remaining > 1) remaining = remaining - 1
   end if
end do

total_offspring = remaining

! Adjust the population size for the next generation such that it
! does not exceed the number of offspring after removal of those
! with negative phenotypic fitness.

current_pop_size = min(current_pop_size, remaining)

! Allow the population size for the next generation potentially 
! to rebound from an earlier reduction in previous generations 
! because of individuals with negative phenotypic fitness.

if(.not.tribal_competition .and. .not.bottleneck_yes .and. &
   remaining > current_pop_size) &
   current_pop_size = min(remaining, pop_size)

! Copy the phenotypic fitnesses into array work_fitness.

work_fitness(1:total_offspring) = pheno_fitness(1:total_offspring)

if (selection_scheme == 2) then

!  For unrestricted probability selection, divide the phenotypic  
!  fitness by a uniformly distributed random number prior to 
!  ranking and truncation.  This procedure allows the probability 
!  of surviving and reproducing in the next generation to be
!  directly related to phenotypic fitness and also for the correct
!  number of individuals to be eliminated to maintain a constant
!  population size.

   do i=1,total_offspring
      work_fitness(i) = work_fitness(i)/(randomnum(1) + 1.d-15)
   end do
   
end if

if (selection_scheme == 3) then

!  For strict proportionality probability selection, rescale the
!  phenotypic fitness values such that the maximum value is one.
!  Then divide the scaled phenotypic fitness by a uniformly
!  distributed random number prior to ranking and truncation.  
!  Allow only those individuals to reproduce whose resulting 
!  ratio of scaled phenotypic fitness to the random number value
!  exceeds one.  This approach ensures that no individual 
!  automatically survives to reproduce regardless of the value
!  of the random number.  But it restricts the fraction of the 
!  offspring that can survive.  Therefore, when the reproduction
!  rate is low, the number of surviving offspring may not be
!  large enough to sustain a constant population size.

   max_work_fitness = 0.d0
   do i=1,total_offspring
      max_work_fitness = max(max_work_fitness, work_fitness(i))
   end do
  
   do i=1,total_offspring
      work_fitness(i) = work_fitness(i)/(max_work_fitness + 1.d-15)
      work_fitness(i) = work_fitness(i)/(randomnum(1) + 1.d-15)
   end do
   
end if

if (selection_scheme == 4) then

!  For partial truncation selection, divide the phenotypic  
!  fitness by the sum of theta and (1. - theta) times a random 
!  number distributed uniformly between 0.0 and 1.0 prior to 
!  ranking and truncation, where theta is the parameter 
!  partial_truncation_value.  This selection scheme is 
!  intermediate between truncation selection and unrestricted 
!  probability selection.  The procedure allows for the correct 
!  number of individuals to be eliminated to maintain a constant 
!  population size.

   do i=1,total_offspring
      work_fitness(i) = work_fitness(i)/(partial_truncation_value &
                        + (1. - partial_truncation_value)*randomnum(1))
   end do
   
end if

! Sort the resulting work fitnesses in ascending order.

sorted_score(1:total_offspring) = work_fitness(1:total_offspring)

if (total_offspring > 1) call heapsort(sorted_score, total_offspring)

if (selection_scheme <= 4) then

!  Apply truncation selection to reduce the population size to
!  current_pop_size.

!  Compute the score cutoff value.

   if(total_offspring > current_pop_size) then
      score_cutoff = sorted_score(total_offspring - current_pop_size)
   else
      score_cutoff = -1000.d0
   end if

   if(selection_scheme == 3) score_cutoff = max(1.d0, score_cutoff)

!  Copy pheno_fitness into array sorted_score for diagnostics
!  purposes.

   sorted_score(1:total_offspring) =  &
                              pheno_fitness(1:total_offspring)

!  Remove those individuals whose score lies below the cutoff
!  value to reduce the population size to its appropriate value.

   current_pop_size = min(current_pop_size, total_offspring)
   remaining = total_offspring

   do i=1,current_pop_size

!     If the work fitness if individual i is below the cutoff 
!     value, find another individual in the pool of excess 
!     offspring whose work fitness is equal to or above the
!     cutoff value and replace the first individual with the 
!     second in the list of reproducing individuals for that
!     generation.

      if(work_fitness(i) < score_cutoff .and. i < remaining) then
         do while(work_fitness(remaining) < score_cutoff .and. &
                  remaining > 1)
            remaining = remaining - 1
         end do
         if(i < remaining) then
            k = gp(remaining)%dm(1,1) + 1
            gp(i)%dm(1:k,1) = gp(remaining)%dm(1:k,1)
            k = gp(remaining)%dm(1,2) + 1
            gp(i)%dm(1:k,2) = gp(remaining)%dm(1:k,2)
            k = gp(remaining)%fm(1,1) + 1
            gp(i)%fm(1:k,1) = gp(remaining)%fm(1:k,1)
            k = gp(remaining)%fm(1,2) + 1
            gp(i)%fm(1:k,2) = gp(remaining)%fm(1:k,2)
            gp(i)%lbmc(:,:,:) = gp(remaining)%lbmc(:,:,:)
            gp(i)%lbf(:,:) = &
            gp(remaining)%lbf(:,:)
            fitness(i) = fitness(remaining)
            pheno_fitness(i) = pheno_fitness(remaining)
            if(remaining > 1) remaining = remaining - 1
         end if
      end if
   end do

   current_pop_size = min(current_pop_size, remaining)

else

   write(6,*) 'ERROR: invalid selection scheme', selection_scheme
   write(9,*) 'ERROR: invalid selection scheme', selection_scheme
   stop

end if

! Compute the mean genotypic and phenotypic fitnesses of the new
! generation after selection.

post_sel_fitness   = 0.d0
mean_pheno_fitness = 0.d0

do i=1,current_pop_size
   post_sel_fitness   = post_sel_fitness   + fitness(i)
   mean_pheno_fitness = mean_pheno_fitness + pheno_fitness(i)
end do

post_sel_fitness   = post_sel_fitness/current_pop_size
mean_pheno_fitness = mean_pheno_fitness/current_pop_size

! Compute the genotypic and phenotypic fitness variances, the 
! covariance of genotypic and phenotypic fitness, and the 
! genotype-phenotype correlation of the new generation.

geno_fitness_variance  = 0.d0
pheno_fitness_variance = 0.d0
covariance = 0.d0

do i=1,current_pop_size
   geno_fitness_variance  = geno_fitness_variance &
                    + (fitness(i) - post_sel_fitness)**2
   pheno_fitness_variance = pheno_fitness_variance &
        + (pheno_fitness(i) - mean_pheno_fitness)**2
   covariance = covariance + fitness(i)*pheno_fitness(i)
end do

geno_fitness_variance  = geno_fitness_variance/current_pop_size
pheno_fitness_variance = pheno_fitness_variance/current_pop_size

post_sel_geno_sd  = sqrt(geno_fitness_variance)
post_sel_pheno_sd = sqrt(pheno_fitness_variance)

covariance = covariance/current_pop_size &
           - post_sel_fitness*mean_pheno_fitness

post_sel_corr = 0.
effect = sqrt(geno_fitness_variance*pheno_fitness_variance)
if(effect > 0.) post_sel_corr = covariance/effect

post_sel_fitness = max(0., post_sel_fitness)

fitness(current_pop_size+1:max_size) = 0.d0

end subroutine selection2

end module selection_module
