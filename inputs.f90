MODULE inputs

implicit none

integer, parameter :: clonal = 1, suppressed = 2, full_sexual = 3

integer pop_size, num_generations, num_linkage_subunits,          &
        bottleneck_generation, bottleneck_pop_size,               &
        num_bottleneck_generations, fitness_distrib_type,         &
        max_del_mutn_per_indiv, max_fav_mutn_per_indiv,           &
        max_neu_mutn_per_indiv,                                   &
        num_initial_fav_mutn, num_indiv_exchanged,                &
        random_number_seed, restart_dump_number,                  &
        haploid_chromosome_number, grow_fission_threshold,        &
        selection_scheme, migration_generations,                  &
        migration_model, num_contrasting_alleles,                 &
        pop_growth_model, plot_allele_gens, verbosity,            &
        poisson_method, recombination_model, carrying_capacity

real    reproductive_rate, mutn_rate,                             &
        genome_size, high_impact_mutn_fraction,                   &
        high_impact_mutn_threshold, fraction_recessive,           &
        dominant_hetero_expression, max_fav_fitness_gain,         &
        recessive_hetero_expression, frac_fav_mutn,               &
        heritability, uniform_fitness_effect_del,                 &
        uniform_fitness_effect_fav, multiplicative_weighting,     &
        fraction_random_death, fraction_self_fertilization,       &
        initial_alleles_mean_effect, non_scaling_noise,           &
        partial_truncation_value, se_nonlinked_scaling,           &
        se_linked_scaling, pop_growth_rate,                       &
        tc_scaling_factor, group_heritability, fraction_neutral,  &
        social_bonus_factor, max_total_fitness_increase,          &
        polygenic_effect, initial_alleles_pop_frac,               &
        initial_alleles_amp_factor

real*8 :: tracking_threshold, extinction_threshold

logical :: fitness_dependent_fertility, dynamic_linkage,             &
           synergistic_epistasis, is_parallel, bottleneck_yes,       &
           restart_case, write_dump, homogenous_tribes,              &
           clonal_haploid, write_vcf,                                &
           upload_mutations, altruistic, allow_back_mutn,            &
           cyclic_bottlenecking, track_neutrals, tribal_competition, &
           polygenic_beneficials, tribal_fission, reseed_rng,        &
           grow_fission

! note: if changing the string length of polygenic_target below,
! need to make corresponding change in polygenic.f90 function poly_match
character case_id*6, data_file_path*80, polygenic_target*40, polygenic_init*40

contains

subroutine read_parameters(nf)
integer nf

namelist /basic/ case_id, mutn_rate, frac_fav_mutn, &
     reproductive_rate, pop_size, num_generations

namelist /mutations/ fitness_distrib_type, &
     genome_size, high_impact_mutn_fraction, &
     high_impact_mutn_threshold, uniform_fitness_effect_del, &
     uniform_fitness_effect_fav, &
     max_fav_fitness_gain, num_initial_fav_mutn, &
     multiplicative_weighting, fraction_recessive, &
     recessive_hetero_expression, dominant_hetero_expression, &
     upload_mutations, allow_back_mutn, se_nonlinked_scaling, &
     se_linked_scaling, synergistic_epistasis

namelist /selection/ fraction_random_death, heritability, &
     non_scaling_noise, fitness_dependent_fertility, &
     selection_scheme, partial_truncation_value

namelist /population/ recombination_model, clonal_haploid, &
     dynamic_linkage, haploid_chromosome_number, &
     fraction_self_fertilization, num_linkage_subunits, &
     pop_growth_model, pop_growth_rate, bottleneck_yes, &
     bottleneck_generation, bottleneck_pop_size, &
     num_bottleneck_generations, carrying_capacity

namelist /substructure/ is_parallel, homogenous_tribes, &
     num_indiv_exchanged, migration_model, migration_generations, &
     tribal_competition, tc_scaling_factor, group_heritability, &
     altruistic, social_bonus_factor, tribal_fission, grow_fission, &
     grow_fission_threshold

namelist /computation/ tracking_threshold, extinction_threshold, &
     max_del_mutn_per_indiv, max_fav_mutn_per_indiv, &
     max_neu_mutn_per_indiv, random_number_seed, reseed_rng, &
     write_dump, write_vcf, restart_case, &
     restart_dump_number, data_file_path, plot_allele_gens, &
     verbosity, poisson_method

namelist /special_apps/ num_contrasting_alleles, &
     max_total_fitness_increase, initial_alleles_pop_frac, &
     initial_alleles_amp_factor, track_neutrals, fraction_neutral, &
     polygenic_effect, polygenic_beneficials, polygenic_target, &
     polygenic_init

read (unit=nf, nml=basic)
read (unit=nf, nml=mutations)
read (unit=nf, nml=selection)
read (unit=nf, nml=population)
read (unit=nf, nml=substructure)
read (unit=nf, nml=computation)
read (unit=nf, nml=special_apps)

end subroutine read_parameters


subroutine write_parameters(nf)

! This routine writes the current parameter values to logical
! unit nf.
integer nf

write(nf,'("&basic")')
write(nf,'(a32,6xa6)')  ' case_id = '              , case_id
write(nf,'(a32,e12.3)') ' mutn_rate = '            , mutn_rate
write(nf,'(a32,f12.7)') ' frac_fav_mutn = '        , frac_fav_mutn
write(nf,'(a32,f12.7)') ' reproductive_rate = '    , reproductive_rate
write(nf,'(a32,i12)')   ' pop_size = '             , pop_size
write(nf,'(a32,i12)')   ' num_generations =   '    , num_generations
write(nf,'("/")')

write(nf,'(/"&mutations")')
write(nf,'(a32,i12)')   ' fitness_distrib_type = ' , fitness_distrib_type
write(nf,'(a32,f12.7)') ' fraction_neutral = '     , fraction_neutral
write(nf,'(a32,e12.3)') ' genome_size = '          , genome_size
write(nf,'(a32,f12.7)') ' high_impact_mutn_fraction = ',   &
                          high_impact_mutn_fraction
write(nf,'(a32,f12.7)') ' high_impact_mutn_threshold = ',  &
                          high_impact_mutn_threshold
write(nf,'(a32,i12)')   ' num_initial_fav_mutn = ' , num_initial_fav_mutn
write(nf,'(a32,f12.7)') ' max_fav_fitness_gain = ' , max_fav_fitness_gain
write(nf,'(a32,f12.7)') ' uniform_fitness_effect_del = ',  &
                          uniform_fitness_effect_del
write(nf,'(a32,f12.7)') ' uniform_fitness_effect_fav = ',  &
                          uniform_fitness_effect_fav
write(nf,'(a32,f12.7)') ' fraction_recessive = '   , fraction_recessive
write(nf,'(a32,f12.7)') ' recessive_hetero_expression = ', &
                          recessive_hetero_expression
write(nf,'(a32,f12.7)') ' dominant_hetero_expression = ',  &
                          dominant_hetero_expression
write(nf,'(a32,f12.7)') ' multiplicative_weighting = ',    &
                          multiplicative_weighting
write(nf,'(a32,l)')     ' synergistic_epistasis = ',       &
                          synergistic_epistasis
write(nf,'(a32,e12.5)') ' se_nonlinked_scaling = ' , se_nonlinked_scaling
write(nf,'(a32,e12.5)') ' se_linked_scaling = '    , se_linked_scaling
write(nf,'(a32,l)')     ' upload_mutations = '     , upload_mutations
write(nf,'(a32,l)')     ' allow_back_mutn = '      , allow_back_mutn
write(nf,'(a32,l)')     ' polygenic_beneficials = ', polygenic_beneficials
write(nf,'(a32,a)')     ' polygenic_init = '       , polygenic_init
write(nf,'(a32,a)')     ' polygenic_target = '     , polygenic_target
write(nf,'(a32,f12.7)') ' polygenic_effect = '     , polygenic_effect
write(nf,'("/")')

write(nf,'(/"&selection")')
write(nf,'(a32,f12.7)') ' fraction_random_death = ', fraction_random_death
write(nf,'(a32,f12.7)') ' heritability = '         , heritability
write(nf,'(a32,f12.7)') ' non_scaling_noise = '    , non_scaling_noise
write(nf,'(a32,l)')     ' fitness_dependent_fertility = ', &
                          fitness_dependent_fertility
write(nf,'(a32,i12)')   ' selection_scheme = '     , selection_scheme
write(nf,'(a32,f12.7)') ' partial_truncation_value = ',    &
                          partial_truncation_value
write(nf,'("/")')

write(nf,'(/"&population")')
write(nf,'(a32,i12)')   ' recombination_model = '  , recombination_model
write(nf,'(a32,l)')     ' clonal_haploid = '       , clonal_haploid
write(nf,'(a32,f12.7)') ' fraction_self_fertilization = ', &
                          fraction_self_fertilization
write(nf,'(a32,i12)')   ' num_contrasting_alleles = ',     &
                          num_contrasting_alleles
write(nf,'(a32,f12.7)') ' initial_alleles_mean_effect = ', &
                          initial_alleles_mean_effect
write(nf,'(a32,f12.7)') ' initial_alleles_pop_frac = ', &
                          initial_alleles_pop_frac
write(nf,'(a32,f12.7)') ' initial_alleles_amp_factor = ', &
                          initial_alleles_amp_factor
write(nf,'(a32,l)')     ' dynamic_linkage = '      , dynamic_linkage
write(nf,'(a32,i12)')   ' haploid_chromosome_number = ',   &
                          haploid_chromosome_number
write(nf,'(a32,i12)')   ' num_linkage_subunits = ' , num_linkage_subunits
write(nf,'(a32,i12)')   ' pop_growth_model = '     , pop_growth_model
write(nf,'(a32,f12.7)') ' pop_growth_rate = '      , pop_growth_rate
write(nf,'(a32,i12)')   ' carrying_capacity = '    , carrying_capacity
write(nf,'(a32,l)')     ' bottleneck_yes = '       , bottleneck_yes
write(nf,'(a32,i12)')   ' bottleneck_generation = ', bottleneck_generation
write(nf,'(a32,i12)')   ' bottleneck_pop_size = '  , bottleneck_pop_size
write(nf,'(a32,i12)')   ' num_bottleneck_generations  = ', &
                          num_bottleneck_generations
write(nf,'("/")')

write(nf,'(/"&substructure")')
write(nf,'(a32,l)')     ' is_parallel = '           , is_parallel
write(nf,'(a32,l)')     ' homogenous_tribes = '     , homogenous_tribes
write(nf,'(a32,i12)')   ' num_indiv_exchanged = '   , num_indiv_exchanged
write(nf,'(a32,i12)')   ' migration_generations = ' , migration_generations
write(nf,'(a32,i12)')   ' migration_model = '       , migration_model
write(nf,'(a32,l)')     ' tribal_competition = '    , tribal_competition
write(nf,'(a32,f12.7)') ' tc_scaling_factor = '     , tc_scaling_factor
write(nf,'(a32,f12.7)') ' group_heritability = '    , group_heritability
write(nf,'(a32,l)')     ' altruistic = '            , altruistic
write(nf,'(a32,f12.7)') ' social_bonus_factor = '   , social_bonus_factor
write(nf,'("/")')

write(nf,'(/"&computation")')
write(nf,'(a32,1pe12.3)') ' tracking_threshold = '  , tracking_threshold
write(nf,'(a32,f12.7)') ' extinction_threshold = '  , extinction_threshold
write(nf,'(a32,i12)')   ' max_del_mutn_per_indiv = ', max_del_mutn_per_indiv
write(nf,'(a32,i12)')   ' max_neu_mutn_per_indiv = ', max_neu_mutn_per_indiv
write(nf,'(a32,i12)')   ' max_fav_mutn_per_indiv = ', max_fav_mutn_per_indiv
write(nf,'(a32,i12)')   ' random_number_seed = '    , random_number_seed
write(nf,'(a32,i12)')   ' poisson_method = '        , poisson_method
write(nf,'(a32,l)')     ' reseed_rng = '            , reseed_rng
write(nf,'(a32,l)')     ' track_neutrals = '        , track_neutrals
write(nf,'(a32,l)')     ' write_dump = '            , write_dump
write(nf,'(a32,l)')     ' write_vcf = '             , write_vcf
write(nf,'(a32,l)')     ' restart_case = '          , restart_case
write(nf,'(a32,i12)')   ' restart_dump_number = '   , restart_dump_number
write(nf,'(a32,i12)')   ' plot_allele_gens = '      , plot_allele_gens
write(nf,'(a32,i12)')   ' verbosity = '             , verbosity
write(nf,'(a20,a,a)')   "  data_file_path = '", trim(data_file_path),"'"
write(nf,'("/")')

write(nf,'(/"&special_apps")')

write(nf,'("/")')

end subroutine write_parameters

subroutine set_default_parameters()

! basic parameters
case_id = 'test00'
mutn_rate = 10.0
frac_fav_mutn = 0.0
reproductive_rate = 2.0
pop_size = 1000
num_generations = 500

! mutations
fitness_distrib_type = 1 ! exponential_mutation_effect
genome_size = 3.e+8
high_impact_mutn_fraction = 0.001
high_impact_mutn_threshold = 0.001
num_initial_fav_mutn = 0
max_fav_fitness_gain = 0.01
uniform_fitness_effect_del = 0.0
uniform_fitness_effect_fav = 0.0
fraction_recessive = 0.0
recessive_hetero_expression = 0.5
dominant_hetero_expression = 0.5
multiplicative_weighting = 0.0
synergistic_epistasis = .false.
   se_nonlinked_scaling = 0
   se_linked_scaling = 0
upload_mutations = .false.
allow_back_mutn = .false.

! selection
fraction_random_death = 0.0
heritability = 0.2
non_scaling_noise = 0.05
fitness_dependent_fertility = .false.
selection_scheme = 2
   partial_truncation_value = 0.5

! population
recombination_model = full_sexual
clonal_haploid = .false.
fraction_self_fertilization = 0.0
dynamic_linkage = .true.
haploid_chromosome_number = 23
num_linkage_subunits = 989
pop_growth_model = 0 ! fixed_population
pop_growth_rate = 0
bottleneck_yes = .false.
   bottleneck_generation = 0
   bottleneck_pop_size = 0
   num_bottleneck_generations  = 0

! substructure
is_parallel = .false.
homogenous_tribes = .true.
num_indiv_exchanged = 1
migration_generations = 10
migration_model = 1
tribal_competition = .false.
tc_scaling_factor = 0.1
group_heritability = 0.0
altruistic = .false.
social_bonus_factor = 1.0
tracking_threshold = 0
extinction_threshold = 0

! computation
max_del_mutn_per_indiv = 10000
max_neu_mutn_per_indiv = 10000
max_fav_mutn_per_indiv = 10000
random_number_seed = 42
reseed_rng = .false.
poisson_method = 0 ! Numerical Recipes
write_dump = .false.
write_vcf = .false.
restart_case = .false.
restart_dump_number = 0
plot_allele_gens = 100
verbosity = 1
data_file_path = './'

! special applications
num_contrasting_alleles = 0
initial_alleles_mean_effect = 0.0
initial_alleles_amp_factor = 1
track_neutrals = .false.
fraction_neutral = 0.0
polygenic_beneficials = .false.
polygenic_init = 'AAAAAA'
polygenic_target = 'TCGTCG'
polygenic_effect = 0.001

end subroutine set_default_parameters

end module inputs
