PROGRAM RIXS

use m_functions, only: deallocate_global_arrays

call readinput
call readbnd
call readsdt
call readrixsdata
call geometry
call kpmesh
call setaux
call printdata
call readrixsmme
call rixsloss
call absorption
call readrix
call writerix
call print_time_statistics
call deallocate_global_arrays( success = .true. )

END PROGRAM RIXS

