	PROGRAM pmds
		USE mpi
                USE globals
  		
  		CHARACTER(LEN=255) :: input_file, dump_file
		CHARACTER(LEN=1) :: pstring
		INTEGER :: ierr, rank, nprocs
                DOUBLE PRECISION :: start_time
                start_time = MPI_WTIME()

  		CALL GETARG(1, input_file)

		CALL MPI_INIT(ierr)
		CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
		CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
		WRITE(UNIT=PSTRING, FMT='(I1)') rank
		dump_file = 'out.dump.' // pstring
  		OPEN(UNIT=99, FILE=dump_file)

  		
  		CALL input_parser(input_file)

  		CLOSE(99)

		CALL MPI_FINALIZE(ierr)

                IF (rank .EQ. 0) THEN
                  OPEN(UNIT=42, FILE='times', ACCESS='APPEND')
                  WRITE(42,*) Natom, nprocs, MPI_WTIME() - start_time
                  CLOSE(42)
                END IF
  		
	END PROGRAM pmds

	SUBROUTINE run_simulation(num_timesteps, pair_style)
		USE mpi
  		
  		USE Globals
		IMPLICIT NONE
  		
  		INTEGER :: num_timesteps, i, nprocs, rank, ierr, atomsp
                INTEGER, ALLOCABLE :: atoms_procs(:)
  		CHARACTER(*), INTENT(IN) :: pair_style
		DOUBLE PRECISION :: Xsend(Maxatom)

		CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
		CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
		atomsp = CEILING(REAL(Natom)/nprocs)
		NAstart = atomsp*rank + 1;
		NAend = MIN(Natom, atomsp*(rank + 1))
                ALLOCATE(atoms_procs(nprocs))
                atoms_procs(1:(nprocs-1)) = atomsp
                atoms_procs(nprocs) = Natom - atomsp*(nprocs-1) + 1
  
  		Nstep = 0
  		
  		OPEN(FILE = "Stats.dat",UNIT=35)
  		
  		WRITE (35,*) "Timestep","Ukin","Pot","Total","Pressure"
  
  		DO WHILE(Nstep .LT. num_timesteps)
    		IF(MOD(Nstep,100) .EQ. 0) WRITE(*,'(A I8)') 'Running timestep: ', Nstep
    		IF(rank .EQ. 0 .AND. MOD(Nstep,10) .EQ. 0) THEN
      			WRITE(99, '(F13.3 F13.3)') (Xx(i), Yy(i), i=1,Natom)
      			FLUSH(99)
    		ENDIF
    		SELECT CASE(pair_style)
      			CASE('soft')
        			A_soft = 19.0*Nstep/num_timesteps + 1.0
        			CALL force_soft_neighbor
        			IF(MOD(Nstep,100) .EQ. 0) THEN
                                        WRITE(*,'(A I4)') 'Calculating soft force'
                                        WRITE (*,*) 'Potential Soft',Upot/Natom
                                        WRITE (*,*) 'Ukin Soft',     Ukin/Natom
                                        WRITE (*,*) 'Tot Soft',     Ukin/Natom + Upot/Natom
                                        WRITE (*,*) 'Pressure', Press
                                END IF
      			CASE('lj')
        			CALL force_neighbor
        			IF(MOD(Nstep,100) .EQ. 0) THEN
                                        WRITE(*,'(A I4)') 'Calculating lj force'
                                        WRITE (*,*) 'Potential Soft',Upot/Natom
                                        WRITE (*,*) 'Ukin Soft',     Ukin/Natom
                                        WRITE (*,*) 'Tot Soft',     Ukin/Natom + Upot/Natom
                                        WRITE (*,*) 'Pressure',Press
                                END IF
                    WRITE (35,*) i,Ukin/Natom,Upot/Natom,(Ukin+Upot)/Natom,Press
               END SELECT
    
    		CALL integrate
    		
    		IF(MOD(Nstep,100) .EQ. 0) WRITE(*,'(A I4)') 'Integrating'
		
		Xsend(1:(NAend - NAstart + 1)) = Xx(NAstart:NAend)
		CALL MPI_ALLGATHER(Xsend, NAend-NAstart+1, MPI_DOUBLE_PRECISION,&
				   Xx, atoms_procs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
		Xsend(1:(NAend - NAstart + 1)) = Yy(NAstart:NAend)
		CALL MPI_ALLGATHER(Xsend, NAend-NAstart+1, MPI_DOUBLE_PRECISION,&
				   Yy, atoms_procs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)

    		Nstep = Nstep + 1
	  	END DO
	  	CLOSE(35)
	  	
	  	OPEN(FILE = 'FinalRes.dat',UNIT=27)
	  	WRITE (27,*) "AtomID","Type","Vx","Vy","Fx","Fy"
	  	DO i=1,Natom
	  		WRITE (27,*) i,AT(i),Vx(i),Vy(i),SQRT(Vx(i)**2+Vy(i)**2),Fx(i),Fy(i)
	  	END DO
	  	CLOSE(27)
	  	
	  	DEALLOCATE(atoms_procs)
  
	END SUBROUTINE run_simulation

	SUBROUTINE broadcast_velocity
		USE mpi
		USE globals
		INTEGER :: ierr
		DOUBLE PRECISION :: Xsend(NAend-NAstart+1)
		Xsend = Vx(NAstart:NAend)
		CALL MPI_ALLGATHER(Xsend, NAend-NAstart+1, MPI_DOUBLE_PRECISION,&
				   Vx, NAend-NAstart+1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
		Xsend = Vy(NAstart:NAend)
		CALL MPI_ALLGATHER(Xsend, NAend-NAstart+1, MPI_DOUBLE_PRECISION,&
				   Vy, NAend-NAstart+1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
	END SUBROUTINE broadcast_velocity
