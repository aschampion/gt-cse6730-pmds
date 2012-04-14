	PROGRAM pmds
		USE mpi
  		
  		CHARACTER(LEN=255) :: input_file
		INTEGER :: ierr

  		CALL GETARG(1, input_file)
  		OPEN(UNIT=99, FILE='out.dump')

		CALL MPI_INIT(ierr)
  		
  		CALL input_parser(input_file)

  		CLOSE(99)

		CALL MPI_FINALIZE(ierr)
  		
	END PROGRAM pmds

	SUBROUTINE run_simulation(num_timesteps, pair_style)
		USE mpi
  		
  		USE Globals
		IMPLICIT NONE
  		
  		INTEGER :: num_timesteps, i, nprocs, rank, ierr, atomsp
  		CHARACTER(*), INTENT(IN) :: pair_style
  		REAL(KIND=8) :: Xo(Natom),Yo(Natom)
		DOUBLE PRECISION :: Xsend(Maxatom)

		CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
		CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
		atomsp = CEILING(REAL(Natom)/nprocs)
		NAstart = atomsp*rank + 1;
		NAend = MIN(Natom, atomsp*(rank + 1))
  
  		Nstep = 0
  
  		DO i=NAstart,NAend
  			Xo(i) = Xx(i)
  		END DO
    
  		DO WHILE(Nstep .LT. num_timesteps)
    		IF(MOD(Nstep,100) .EQ. 0) WRITE(*,'(A I8)') 'Running timestep: ', Nstep+1
    		IF(MOD(Nstep,10) .EQ. 0) THEN
      			WRITE(99, '(F13.3 F13.3)') (Xx(i), Yy(i), i=1,Natom)
      			FLUSH(99)
    		ENDIF
    		SELECT CASE(pair_style)
      			CASE('soft')
        			A_soft = 19.0*Nstep/num_timesteps + 1.0
        			CALL force_soft()
        			IF(MOD(Nstep,100) .EQ. 0) THEN
                                        WRITE(*,'(A I4)') 'Calculating soft force'
                                        WRITE (*,*) 'Potential Soft',Upot/Natom
                                        WRITE (*,*) 'Ukin Soft',     Ukin/Natom
                                        WRITE (*,*) 'Tot Soft',     Ukin/Natom + Upot/Natom
                                END IF

      			CASE('lj')
        			CALL force()
        			IF(MOD(Nstep,100) .EQ. 0) THEN
                                        WRITE(*,'(A I4)') 'Calculating LJ force'
                                        WRITE (*,*) 'Potential LJ',Upot/Natom
                                        WRITE (*,*) 'Ukin LJ',     Ukin/Natom
                                        WRITE (*,*) 'Tot LJ',     Ukin/Natom + Upot/Natom
                                END IF	
               END SELECT
    
    		CALL integrate
    		IF(MOD(Nstep,100) .EQ. 0) WRITE(*,'(A I4)') 'Integrating'
		
		Xsend = Xx(NAstart:NAend)
		CALL MPI_ALLGATHER(Xsend, NAend-NAstart, MPI_DOUBLE_PRECISION,&
				   Xx, NAend-NAstart, MPI_DOUBLE, MPI_COMM_WORLD, ierr)
		Xsend = Yy(NAstart:NAend)
		CALL MPI_ALLGATHER(Xsend, NAend-NAstart, MPI_DOUBLE_PRECISION,&
				   Yy, NAend-NAstart, MPI_DOUBLE, MPI_COMM_WORLD, ierr)

    		Nstep = Nstep + 1
	  	END DO
  
  		DO i=1,Natom
			WRITE (*,*) Xx(i),Xo(i)
  		END DO
  		
	END SUBROUTINE run_simulation

	SUBROUTINE broadcast_velocity
		USE mpi
		USE globals
		INTEGER :: ierr
		DOUBLE PRECISION :: Xsend(NAend-NAstart)
		Xsend = Vx(NAstart:NAend)
		CALL MPI_ALLGATHER(Xsend, NAend-NAstart, MPI_DOUBLE_PRECISION,&
				   Vx, NAend-NAstart, MPI_DOUBLE, MPI_COMM_WORLD, ierr)
		Xsend = Vy(NAstart:NAend)
		CALL MPI_ALLGATHER(Xsend, NAend-NAstart, MPI_DOUBLE_PRECISION,&
				   Vy, NAend-NAstart, MPI_DOUBLE, MPI_COMM_WORLD, ierr)
	END SUBROUTINE broadcast_velocity
