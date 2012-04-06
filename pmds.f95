program pmds
  character(len=255) :: input_file

  call getarg(1, input_file)
  call input_parser(input_file)
end program pmds

subroutine run_simulation(num_timesteps, pair_style)
  implicit none

  integer :: num_timesteps
  character(*), intent(in) :: pair_style

  do while(Nstep.lt.num_timesteps)
    select case (pair_style)
      case ('soft')
        A_soft = 19.0*(Nstep/num_timesteps) + 1.0
        call force_soft()
      case ('lj/cut')
        call force()
    end select

    call integrate()

    Nstep = Nstep + 1
  end do
end subroutine
