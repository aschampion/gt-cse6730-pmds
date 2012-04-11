	SUBROUTINE input_parser(filename)
  		
 		USE Globals
  		
  		CHARACTER(*), INTENT(IN) :: filename
  		CHARACTER(LEN=1024) :: line, instr, args
  		CHARACTER(LEN=256) :: s1, s2, pair_style
  		INTEGER :: i1, i2
  		REAL(KIND=8) :: v1, v2, v3

  		OPEN(UNIT=10, FILE='in.micelle')

  		DO WHILE(.TRUE.)
    		READ(10, '(A)',END=99) line		
    
    		instr_end = 2
    		DO WHILE (line(instr_end:instr_end) .NE. char(9))
      			instr_end = instr_end + 1
    		END DO

    		instr = line(1:(instr_end - 1))
    		DO WHILE(line(instr_end:instr_end) .EQ. char(9))
      			instr_end = instr_end + 1
    		END DO
    		args = ADJUSTL(TRIM(line(instr_end:LEN(line))))

    		SELECT CASE (instr)
      			CASE ('bond_coeff')
        			!Bond Type ID, K_bond, Rcut_bond
        			READ(args, *) i1, v1, v2
        			K_bond = v1
        			Rcut_bond = v2
        			WRITE (*,*) "Reading harmonic bond data"
      			CASE ('pair_coeff')
        			SELECT CASE (pair_style)
          				CASE ('lj/cut')
            				!Type ID 1, Type ID 2, P1, P2, Radius Cutoff
            				READ(args, *) i1, i2, v1, v2, v3
            				Rcut(i1, i2) = v3
            				sigma_matrix(i1, i2) = v2
            				epsilon_matrix(i1, i2) = v1
            				WRITE (*,*) "Reading LJ data"
          				CASE ('soft')
            				!Type ID 1, Type ID 2, A, Radius Cutoff
	            			READ(args, *) s1, s2, v1, v2
            				A_soft = v1
            				Rcut_soft = v2
            				WRITE (*,*) "Reading soft force data"
        			END SELECT
      			CASE ('pair_style')
        			sigma_matrix(1:MaxNumtypes,1:MaxNumtypes) = 0
        			epsilon_matrix(1:MaxNumtypes,1:MaxNumtypes) = 0
        			READ(args, *) pair_style
        			WRITE (*,*) "Pair style",pair_style
      			CASE ('read_data')
      				WRITE (*,*) "Reading initial data from ",args
        			CALL data_parser(args)
        			WRITE (*,*) Natom,"atoms"
        			WRITE (*,*) Nbond,"bonds"
        			WRITE(99, '(I8)') Natom
      			CASE ('run')
        			!Timesteps
        			READ(args, *) i1
        			WRITE (*,*) "Running simualation for",i1,"steps"
        			CALL run_simulation(i1, pair_style)
      			CASE ('velocity')
        			READ(args, *) s1, s2, v1, i1
        			Temp_target = v1 
        			CALL initial       
        			WRITE (*,*) "Setting velocity for temperature",Temp_target
      			CASE ('#')  
      			CASE ('')
      			CASE DEFAULT
        	!write(*,'(A A)') 'Unrecognized instruction: ', trim(instr)
    		END SELECT
  		END DO

99  	CLOSE(10)
	END SUBROUTINE input_parser

	SUBROUTINE data_parser(filename)

  		USE Globals
  		CHARACTER(*), INTENT(IN) :: filename
  		CHARACTER(LEN=1024) :: line
  		CHARACTER(LEN=80) :: state = ''
  		INTEGER :: i1, i2, i3
  		REAL(KIND=8) :: v1, v2, v3

  		OPEN(UNIT=20,FILE=filename)

  		DO WHILE(.TRUE.)
    		READ(20, '(A)', END=199) line

    		SELECT CASE (line)
      			CASE (' Atoms')
        			state = 'Atoms'
        			Natom = 0
        			WRITE (*,*) "Reading atom data"
      			CASE (' Bonds')
        			state = 'Bonds'
        			Nbond = 0
        			WRITE (*,*) "Reading bond data"
      			case ('')
      			CASE DEFAULT
        			SELECT CASE (state)
          				CASE ('Atoms')
            				!Atom ID, Molecule ID, Atom Type, X, Y, Z
            				READ(line, *) i1, i2, i3, v1, v2, v3
            				At(i1) = i3
            				Xx(i1) = v1
            				Yy(i1) = v2
            				IF(i1 .GT. Natom) Natom = i1
          				CASE ('Bonds')
            				!Bond ID, Bond Type, Atom ID 1, Atom ID 2
            				READ(line, *) i1, i2, i3, i4
            				bondlist(1, i1) = i3
            				bondlist(2, i1) = i4
!WRITE (*,*) "parser bondlist(1, i1) , bondlist(2, i1) ",bondlist(1, i1),bondlist(2, i1)
            				IF(i1 .GT. Nbond) Nbond = i1
          				case ('')
        		END SELECT
    		END SELECT
  		END DO

199 	CLOSE(20)
	END SUBROUTINE data_parser
