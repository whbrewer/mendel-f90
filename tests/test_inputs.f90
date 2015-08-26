program test_inputs
use unittest
use inputs

open (5, file='mendel.in',status='old')
call read_parameters(5)
close(5)
call write_parameters(6)

end program test_inputs
