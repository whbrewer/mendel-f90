module profiler

real :: sec(12), tin, tout

contains

subroutine second(time)

real time, times(2)

time = etime(times)
time = times(1)

end subroutine second

subroutine profile(unit)

! Routine to dump timing information at end of run

use inputs
integer i, unit

!if (myid == 0) then
   write(unit,10) sec(2)
   sec(3) = sec(3)/num_generations ! currently not used
   sec(4) = sec(4)/num_generations*migration_generations
   sec(5) = sec(5)/num_generations ! *actual_offspring
   sec(6) = sec(6)/num_generations
   sec(7) = sec(7)/num_generations ! every 20 generations
   sec(8) = sec(8)/num_generations ! every 100 generations
   write(unit,20) (1000.*sec(i),i=4,6)
   write(unit,30) (1000.*sec(i),i=7,12)
!end if

 10   format(/10x,"CPU SECONDS USED PER PROCESSOR:"// &
       10x,"TOTAL     ",f10.3)

 20   format(/10x,"PRIMARY SUBROUTINES (ms/gen/proc):"// &
       10x,"migration        ",f10.3/ &
       10x,"offspring        ",f10.3/ &
       10x,"selection        ",f10.3)

 30   format(/10x,"AUXILIARY ROUTINES (ms/gen/proc):"// &
       10x,"mutn statistics  ",f10.3/ &
       10x,"allele statistics",f10.3/ &
       10x,"read_restart_dump",f10.3/ &
       10x,"write_output_dump",f10.3/ &
       10x,"write_sample     ",f10.3/ &
       10x,"back_mutn        ",f10.3/10x)

end subroutine profile

end module profiler
