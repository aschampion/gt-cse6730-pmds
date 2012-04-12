program pmds
  character(len=255) :: input_file

  call getarg(1, input_file)
  open(99, file='out.dump', status='REPLACE')
  call input_parser(input_file)

  close(99)
end program pmds

subroutine run_simulation(num_timesteps, pair_style)
  implicit none

  include 'globals.inc'

  integer :: num_timesteps, i
  character(*), intent(in) :: pair_style
  REAL(KIND=8) :: Xo(Natom),Yo(Natom)
<<<<<<< .mine

<<<<<<< .mine
  !open(99, file='out.dump')
=======
  !#open(99, file='out.dump')
=======
>>>>>>> .r40
>>>>>>> .r39
  
  Nstep = 0
  
  DO i=1,Natom
  	Xo(i) = Xx(i)
  END DO
    
  do while(Nstep.lt.num_timesteps)
    if (mod(Nstep, 100).eq.0) write(*,'(A I8)') 'Running timestep: ', Nstep+1
<<<<<<< .mine
    if (mod(Nstep, 10000).eq.0) then
<<<<<<< .mine
      !write(99, '(I5)') Natom
      !write(99, '(F13.3 F13.3)') (Xx(i), Yy(i), i=1,Natom)
      !flush(99)
=======
!      #write(99, '(I5)') Natom
 !     #write(99, '(F13.3 F13.3)') (Xx(i), Yy(i), i=1,Natom)
  !    #flush(99)
=======
    if (mod(Nstep, 10).eq.0) then
      write(99, '(F13.3 F13.3)') (Xx(i), Yy(i), i=1,Natom)
      flush(99)
>>>>>>> .r40
>>>>>>> .r39
    endif
    select case (pair_style)
      case ('soft')
        A_soft = 19.0*Nstep/num_timesteps + 1.0
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
  
  DO i=1,Natom
  	IF(Xx(i) .LT. 0.0D0) WRITE (*,*) Xx(i),Xo(i)
  END DO
<<<<<<< .mine

<<<<<<< .mine
  !close(99)
=======
 ! #close(99)
=======
>>>>>>> .r40
>>>>>>> .r39
end subroutine
