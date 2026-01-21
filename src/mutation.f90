subroutine mutation(mutn_indx,mutn_type,lb,hap_id,fitness_effect)
use inputs
use random_pkg
include 'common.h'

integer, parameter :: del = 1, fav = 2, neu = 3
integer, intent(out) :: mutn_indx, mutn_type, lb, hap_id
real*8,  intent(out) :: fitness_effect

integer :: mutn, i
real    :: x

! Select a random linkage block for the new mutation.

lb = min(num_linkage_subunits, 1+int(num_linkage_subunits*randomnum(1)))

hap_id = min(2, 1 + int(2.*randomnum(1)))

! Determine whether new mutation is deleterious, neutral, or favorable.

x = randomnum(1)

if(x < frac_fav_mutn*(1. - fraction_neutral) ) then 
   mutn_type = fav
elseif(x < 1. - fraction_neutral) then 
   mutn_type = del
else
   mutn_type = neu
end if

! Compute the mutation index that is used to specify
! its fitness. mutn is an integer index which represents
! the mutations fitness.  This value is later offset 
! by the linkage block number as follows:
! mutn_indx = (lb - 1)*lb_modulo + mutn

10 continue

x = randomnum(1)

if(tracking_threshold /= 1.0) then  ! if we're tracking the mutations
   if(mutn_type == fav) then
      mutn = min(lb_modulo-2, int(x/fav_scale))
   elseif (mutn_type == del) then
      mutn = min(lb_modulo-2, int(x/del_scale))
   else  ! neutral mutations have zero effect on fitness
      if(polygenic_beneficials) then
         mutn = min(4, 1 + int(4.*randomnum(1)))
      else
         mutn = int(x*(lb_modulo-2))
      end if
   end if
else  ! not tracking mutations
   mutn = 1 
end if

! Check to make sure mutation that is being generated does not
! correspond with one that has already been uploaded, so
! when we track the mutations, they will just include the effects
! of the uploaded mutations.
if(upload_mutations) then
   do i = 1, num_uploaded_mutn
      if (uploaded_mutn(i) == mutn) then
         write(*,*)'WARNING: trying to generate a mutation', &
              ' with same id as uploaded mutation.   ', &
              'Computing new mutation id.'
         goto 10
      end if
   end do
end if

if(.not.polygenic_beneficials) then
   mutn_indx = (lb - 1)*lb_modulo + mutn
else
   mutn_indx = mutn
end if

! For parallel cases add the tribe number to the mutation index
! The reason for this is to create unique mutation numbers 
! arising from each tribe.
if(is_parallel) mutn_indx = mutn_indx + myid

! Compute the fitness effect associated with this new mutation.
     
! When parameter fitness_distrib_type is 1, the fitness
! effect e is obtained from the mutation index mutn using a
! distribution function of the form

! e = exp(-alpha_del*x**gamma) ,
 
! where alpha_del is log(genome_size) and x is a random
! number uniformly distributed between zero and one.
  
! When parameter fitness_distrib_type is 0, the fitness
! effect is constant and given by the expression
   
! e = uniform_fitness_effect_del
    
! For favorable mutations, 
! e = uniform_fitness_effect_fav
     
if(fitness_distrib_type == 1) then  ! Natural mutation distribution
   if(mutn_type == fav) then
      fitness_effect = max_fav_fitness_gain*dexp(-alpha_fav*x**gamma_fav)
   elseif(mutn_type == del) then
      fitness_effect = dexp(-alpha_del*x**gamma_del)
   else  ! neutral
      fitness_effect = 0.
   end if
else if (fitness_distrib_type == 2) then  ! All mutations neutral
   fitness_effect = 0.
else  ! All mutations have equal effect
   if(mutn_type == fav) then
      fitness_effect = uniform_fitness_effect_fav
   elseif(mutn_type == del) then
      fitness_effect = uniform_fitness_effect_del
   else  ! neutral mutations have zero fitness effect
      fitness_effect = 0.
   end if
end if

! If neutrals are not being tracked and the fraction of neutrals
! is non-zero, assign the appropriate fraction of mutations zero
! fitness effect.
     
! NOTE: Since we force the user to track all mutations when 
! including neutrals, the following code will only rarely,
! if ever, be executed.  There are some problems in accurately
! correctly computing polymorphisms when all the mutations are 
! not tracked.

if(.not. track_neutrals .and. fraction_neutral > 0. .and. &
     mutn_type /= fav) then
   if(randomnum(1) < fraction_neutral) fitness_effect = 0.
end if

! Identify the appropriate fraction of new mutations as
! recessive.  To distinguish recessive mutations from the
! dominant ones, label the recessives by assigning them a
! negative mutation index.  Make all neutral mutations 
! dominant.

if(fraction_recessive > 0. .and. fitness_effect > 0.) then
   if(randomnum(1) < fraction_recessive) mutn_indx = -mutn_indx
end if

end subroutine mutation

subroutine favorable_mutn(fmutn,lb_mutn_count,linkage_block_fitness,uid, &
                          effect,mutn_indx)

! This routine generates a random mutation in a randomly chosen
! linkage block with a randomly chosen haploid identity in a
! randomly chosen individual in the population and modifies the
! linkage block fitness to reflect the resulting fitness change.

use inputs
use random_pkg
include 'common.h'

integer, intent(in),  optional :: uid
integer, intent(out), optional :: mutn_indx
real,    intent(in),  optional :: effect
integer fmutn(max_fav_mutn_per_indiv/2,2,*)
integer lb_mutn_count(num_linkage_subunits,2,3,*)
real*8 linkage_block_fitness(num_linkage_subunits,2,*)
real*8 fitness_gain
real w, x
integer id, lb, hap_id, mutn, num_mutn, j

! Specify the new random favorable mutation. 

! Generate the index of the random individual.

if(uid > 0) then
   id = uid ! mainly for polygenic beneficials
else
   id = min(current_pop_size, 1 + int(current_pop_size*randomnum(1)))
endif

! Generate the linkage block index.

lb = min(num_linkage_subunits, &
     1 + int(num_linkage_subunits*randomnum(1)))

! Generate the haploid identity.

hap_id = min(2, 1 + int(2.*randomnum(1)))

! Generate a random index mutn to specify the fitness effect
! associated with the mutation.

x = randomnum(1)

if(tracking_threshold /= 1.0) then
   mutn = min(lb_modulo-2, int(x/fav_scale))
else
   mutn = 1
end if

! Add an offset to assign it to the appropriate linkage block.

mutn_indx = mutn + (lb - 1)*lb_modulo 

! Specify whether the mutation is dominant or recessive.
! (Recessives have a negative mutation index.) 

if(fraction_recessive > 0.) then
  if(randomnum(1) < fraction_recessive) mutn_indx = -mutn_indx
end if

! Increment the favorable mutation count for the appropriate
! individual, linkage block, and haploid index.

if (.not. polygenic_beneficials) then
   lb_mutn_count(lb,hap_id,2,id) = lb_mutn_count(lb,hap_id,2,id) + 1
end if

! Compute the fitness factor associated with this new mutation.
! Incorporate this fitness contribution into the fitness of the
! the appropriate linkage block. 

! When parameter fitness_distrib_type is 1, the fitness
! factor f is obtained from the mutation index mutn using a
! distribution function of the form
!
!   f = (1. + max_fav_fitness_gain*exp(-alpha_fav*x**gamma_fav)
!
! where max_fav_fitness_gain is an input parameter, alpha_fav is 
! is log(genome_size*max_fav_fitness_gain) and x is a 
! random number uniformly distributed between zero and one.
!
! When parameter fitness_distrib_type is 0, the fitness
! factor is constant and given by the expression
!
!   f = 1. + uniform_fitness_effect_fav

if(effect > 0) then ! case of polygenic_beneficials with fixed effect
   fitness_gain = effect
else
   if(fitness_distrib_type == 1) then ! Natural distribution
      fitness_gain = max_fav_fitness_gain*dexp(-alpha_fav*x**gamma_fav)
   else if (fitness_distrib_type == 2) then ! All mutn neutral
      fitness_gain = 0
   else ! All mutations have equal effect
      fitness_gain = uniform_fitness_effect_fav
   end if
end if

! Track this mutation if its fitness gain exceeds the value of
! tracking_threshold. 

if(fitness_gain > tracking_threshold) then

!  Test to see if the storage limit of array fmutn has been
!  exceeded.  (Note that we are using the first slot to hold the
!  actual mutation count.)

   num_mutn = fmutn(1,hap_id,id) + 1 

   if(num_mutn + 1 > max_fav_mutn_per_indiv/2) then
      write(6,*) 'Favorable mutations exceed the storage limit'
      write(9,*) 'Favorable mutations exceed the storage limit'
      stop
   end if

   fmutn(1,hap_id,id) = num_mutn

!  Insert new mutation such that mutations are maintained
!  in ascending order of their absolute value.

   j = num_mutn

   do while(abs(fmutn(j,hap_id,id)) > abs(mutn_indx) &
               .and. j > 1)
      fmutn(j+1,hap_id,id) = fmutn(j,hap_id,id)
      j = j - 1
   end do

   fmutn(j+1,hap_id,id) = mutn_indx

end if

! Recessive mutations (identified as such with a negative
! mutation index) here incur only recessive_hetero_expression
! times of their fitness gain, while dominant mutations incur
! only dominant_hetero_expression times their fitness gain.
! The full fitness gain is realized only when a mutation occurs 
! in both instances of its linkage block, that is, is homozygous.

if(mutn_indx < 0) then
   fitness_gain = recessive_hetero_expression*fitness_gain
else
   fitness_gain =  dominant_hetero_expression*fitness_gain
end if

w = multiplicative_weighting

linkage_block_fitness(lb,hap_id,id) = &
   (linkage_block_fitness(lb,hap_id,id) + (1. - w)*fitness_gain) &
                                      * (1.d0 + w *fitness_gain)

end subroutine favorable_mutn

subroutine back_mutn(dmutn,fmutn,lb_fitness,lb_mutn_count)

use inputs
use profiler
use random_pkg
include 'common.h'
integer dmutn(max_del_mutn_per_indiv/2,2)
integer fmutn(max_fav_mutn_per_indiv/2,2)
integer lb_mutn_count(num_linkage_subunits,2,3)
integer lb, hap_id, mutn, i, j, tries, idorf, random_index
integer num_fmutns, num_dmutns, decode_lb
real*8  lb_fitness(num_linkage_subunits,2)
real*8  fitness, decode_fitness_del, decode_fitness_fav
real    w
logical fav
call second(tin)

tries = 0
 10   continue
tries = tries + 1
if(tries > 10) return

! Determine whether back mutation is deleterious or favorable.

if(randomnum(1) < frac_fav_mutn) then
   fav = .true.
else
   fav = .false.       
end if

! Select either 1 or 2 for haplotype.      

hap_id = min(2, 1 + int(2.*randomnum(1)))

if(fav) then
   num_fmutns = fmutn(1,hap_id)
   if(num_fmutns == 0) goto 10
!  Choose a random mutation from the individual
   random_index = min(num_fmutns, &
                  int(num_fmutns*randomnum(1))+1)
   mutn = fmutn(1+random_index,hap_id)
   idorf = 2
else
   num_dmutns = dmutn(1,hap_id)
   if(num_dmutns == 0) goto 10
   random_index = min(num_dmutns, &
                  int(num_dmutns*randomnum(1))+1)
   mutn = dmutn(1+random_index,hap_id)
   idorf = 1
endif
   
! Compute the linkage block of the mutation

lb = decode_lb(mutn)

! Decrement the linkage block mutation counter

lb_mutn_count(lb,hap_id,idorf) = lb_mutn_count(lb,hap_id,idorf)-1

w = multiplicative_weighting

! Rebuild the mutation list excluding the random_index, compute the
! fitness associated with this mutation, and decrement the total
! number of mutations.

if(fav) then
   do i=1+random_index,num_fmutns
      fmutn(i,hap_id) = fmutn(i+1,hap_id)
   end do
   fmutn(num_fmutns+1,hap_id) = num_linkage_subunits*lb_modulo + 1
   fitness = decode_fitness_fav(mutn)
   fmutn(1,hap_id) = fmutn(1,hap_id) - 1
else
   do i=1+random_index,num_dmutns
      dmutn(i,hap_id) = dmutn(i+1,hap_id)
   end do
   dmutn(num_dmutns+1,hap_id) = num_linkage_subunits*lb_modulo + 1
   fitness = decode_fitness_del(mutn)
   dmutn(1,hap_id) = dmutn(1,hap_id) - 1
end if

! Make appropriate adjustments to linkage block fitness.

lb_fitness(lb,hap_id) = (lb_fitness(lb,hap_id) - (1. - w)*fitness) &
                       *(1.d0 - w*fitness)

call second(tout)
sec(12) = sec(12) + tout - tin

end subroutine back_mutn

integer function encode_mutn(fitness,lb,dominance)
include 'common.h'
real fitness
integer lb, dominance 
encode_mutn = dominance*((lb-1)*lb_modulo+lb_modulo*abs(fitness))
return
end function encode_mutn

subroutine decode_mutn_del(mutn,lb,dominance,fitness)
include 'common.h'
real*8 fitness, x
integer mutn, lb, dominance 
dominance = sign(1,mutn)
lb  = abs(mutn)/lb_modulo
x   = mod(abs(mutn), lb_modulo)*del_scale
fitness = -dexp(-alpha_del*x**gamma_del)
if(x >= 1.d0) fitness = 0.d0
return
end subroutine decode_mutn_del
 
integer function decode_lb(mutn)
include 'common.h'
integer mutn
decode_lb = abs(mutn)/lb_modulo
return
end function decode_lb

real*8 function decode_fitness_fav(mutn)
use inputs
include 'common.h'
integer mutn, mtn
mtn = mod(abs(mutn), lb_modulo)
decode_fitness_fav = max_fav_fitness_gain*dexp(-alpha_fav &
                     *(real(mtn)*fav_scale)**gamma_fav)
return
end function decode_fitness_fav

real*8 function decode_fitness_del(mutn)
include 'common.h'
integer mutn
real*8 x
x = mod(abs(mutn), lb_modulo)*del_scale
if(x >= 1.d0) then
   decode_fitness_del = 0.d0
else
   decode_fitness_del = -dexp(-alpha_del*x**gamma_del)
end if
return
end function decode_fitness_del

integer function decode_mutn_id(mutn, mutn_type) 
include 'common.h'
integer mutn, mutn_type ! -1 = deleterious, 0 = neutral, 1 = favorable
real*8 x
if (mutn_type < 0) then
   decode_mutn_id = mod(abs(mutn), lb_modulo)*del_scale
else if (mutn_type == 0) then
   decode_mutn_id = mod(abs(mutn), lb_modulo)
else
   decode_mutn_id = mod(abs(mutn), lb_modulo)*fav_scale
end if
end function decode_mutn_id
