PROGRAM RIXS

call prolog
call readinput
call readbnd
call readsdt
call read_RIXS_data
call geometry
call kpmesh
call set_aux_arrays
call printdata
call read_RIXS_mme
call rixsloss
call absorption
call readrix
call write_results
call print_time_statistics
call epilog

END PROGRAM RIXS
