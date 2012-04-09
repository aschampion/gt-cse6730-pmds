subroutine input_parser(filename)
  character(*), intent(in) :: filename

  include 'globals.inc'

  character(len=1024) :: line, instr, args
  character(len=256) :: s1, s2, pair_style
  integer :: i1, i2
  double precision :: v1, v2, v3

  open(10, file='in.micelle', access='SEQUENTIAL')

  do while(.true.)
    read(10, '(A)', end=99) line
    
    instr_end = 2
    do while (line(instr_end:instr_end).ne.char(9))
      instr_end = instr_end + 1
    end do

    instr = line(1:(instr_end - 1))
    do while(line(instr_end:instr_end).eq.char(9))
      instr_end = instr_end + 1
    end do
    args = adjustl(trim(line(instr_end:len(line))))

    select case (instr)
      case ('bond_coeff')
        ! Bond Type ID, K_bond, Rcut_bond
        read(args, *) i1, v1, v2
        K_bond = v1
        Rcut_bond = v2
        WRITE (*,*) "Reading harmonic bond data"
      case ('pair_coeff')
        select case (pair_style)
          case ('lj/cut')
            ! Type ID 1, Type ID 2, P1, P2, Radius Cutoff
            read(args, *) i1, i2, v1, v2, v3
            Rcut(i1, i2) = v3
            sigma_matrix(i1, i2) = v2
            epsilon_matrix(i1, i2) = v1
            WRITE (*,*) "Reading LJ data"
          case ('soft')
            ! Type ID 1, Type ID 2, A, Radius Cutoff
            read(args, *) s1, s2, v1, v2
            A_soft = v1
            Rcut_soft = v2
            WRITE (*,*) "Reading soft force data"
        end select
      case ('pair_style')
        read(args, *) pair_style
        WRITE (*,*) "Pair style",pair_style
      case ('read_data')
      	WRITE (*,*) "Reading initial data from ",args
        call data_parser(args)
        WRITE (*,*) Natom,"atoms"
        WRITE (*,*) Nbond,"bonds"
      case ('run')
        ! Timesteps
        read(args, *) i1
        WRITE (*,*) "Running simualation for",i1,"steps"
        call run_simulation(i1, pair_style)
      case ('velocity')
        read(args, *) s1, s2, v1, i1
        Temp_target = v1 
        call initial       
        WRITE (*,*) "Setting velocity for temperature",Temp_target
      case ('#')  
      case ('')
      case default
        !write(*,'(A A)') 'Unrecognized instruction: ', trim(instr)
    end select
  end do

99  close(10)
end subroutine input_parser

subroutine data_parser(filename)
  character(*), intent(in) :: filename

  include 'globals.inc'

  character(len=1024) :: line
  character(len=80) :: state = ''
  integer :: i1, i2, i3
  double precision :: v1, v2, v3

  open(20, file=filename, access='SEQUENTIAL')

  do while(.true.)
    read(20, '(A)', end=199) line

    select case (line)
      case (' Atoms')
        state = 'Atoms'
        Natom = 0
        WRITE (*,*) "Reading atom data"
      case (' Bonds')
        state = 'Bonds'
        Nbond = 0
        WRITE (*,*) "Reading bond data"
      case ('')
      case default
        select case (state)
          case ('Atoms')
            ! Atom ID, Molecule ID, Atom Type, X, Y, Z
            read(line, *) i1, i2, i3, v1, v2, v3
            At(i1) = i3
            Xx(i1) = v1
            Yy(i1) = v2
            if (i1.gt.Natom) Natom = i1
          case ('Bonds')
            ! Bond ID, Bond Type, Atom ID 1, Atom ID 2
            read(line, *) i1, i2, i3, i4
            bondlist(1, i1) = i3
            bondlist(2, i1) = i4
            if (i1.gt.Nbond) Nbond = i1
          case ('')
        end select
    end select
  end do  

199 close(20)
end subroutine data_parser
