implicit none

integer :: lb_modulo
integer :: new_mutn_count, uploaded_mutn(1000)
integer :: dump_number, current_pop_size, live_pop_size
integer :: mutn_per_indiv, fav_fixed, fav_lost
integer :: global_pop_size, new_pop_size
integer :: num_uploaded_mutn, num_back_mutn
integer :: myid, num_tribes, ierr, msg_num, MYCOMM

real :: del_scale, fav_scale, gamma_del, gamma_fav
real :: poisson_mean, ica_mean_effect, fav_mean_freq

real*8 :: alpha_del, alpha_fav
real*8 :: tribal_fitness_factor,                                    &
          tribal_fitness_factor_scaled, tribal_fitness,             &
          tribal_noise, social_bonus
real*8 :: global_genetic_fitness, fertility_factor
real*8 :: tracked_fav_mutn, total_del_mutn, tracked_del_mutn

common /mndl1/ new_mutn_count, uploaded_mutn, dump_number, lb_modulo, &
               current_pop_size, mutn_per_indiv, fav_fixed, fav_lost, &
               global_pop_size, new_pop_size, num_uploaded_mutn,      &
               num_back_mutn, myid, num_tribes, MYCOMM

common /mndl2/ del_scale, fav_scale, gamma_del, gamma_fav, &
               poisson_mean, ica_mean_effect, fav_mean_freq

common /mndl3/ alpha_del, alpha_fav, &
  	            tribal_fitness_factor, tribal_fitness_factor_scaled,     &
	            tribal_fitness, tribal_noise, social_bonus,              &
	            global_genetic_fitness, fertility_factor,                &
                tracked_fav_mutn, total_del_mutn, tracked_del_mutn
