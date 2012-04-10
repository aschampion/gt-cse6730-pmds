program pmds
  character(len=255) :: input_file

  call getarg(1, input_file)
  call input_parser(input_file)
end program pmds

subroutine run_simulation(num_timesteps, pair_style)
  implicit none

  include 'globals.inc'

  integer :: num_timesteps, i
  character(*), intent(in) :: pair_style

  #open(99, file='out.dump')
  
  Nstep = 0
    
  do while(Nstep.lt.num_timesteps)
    if (mod(Nstep, 100).eq.0) write(*,'(A I8)') 'Running timestep: ', Nstep+1
    if (mod(Nstep, 10000).eq.0) then
      #write(99, '(I5)') Natom
      #write(99, '(F13.3 F13.3)') (Xx(i), Yy(i), i=1,Natom)
      #flush(99)
    endif
    select case (pair_style)
      case ('soft')
        A_soft = 19.0*(Nstep/num_timesteps) + 1.0
        call force_soft()
        if (mod(Nstep, 100).eq.0) write(*,'(A I4)') 'Calculating soft force'
      case ('lj')
        call force()
        if (mod(Nstep, 100).eq.0) write(*,'(A I4)') 'Calculating LJ force'
    end select
    
    call integrate()
    if (mod(Nstep, 100).eq.0) write(*,'(A I4)') 'Integrating'

    Nstep = Nstep + 1
  end do

  #close(99)
end subroutine
