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
		DOUBLE PRECISION :: Fsend(Maxatom)
		DOUBLE PRECISION :: PressSend

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
        			CALL force_soft_neighbor
        			IF(MOD(Nstep,100) .EQ. 0) THEN
                                        WRITE(*,'(A I4)') 'Calculating soft force'
                                        WRITE (*,*) 'Potential Soft',Upot/Natom
                                        WRITE (*,*) 'Ukin Soft',     Ukin/Natom
                                        WRITE (*,*) 'Tot Soft',      Ukin/Natom + Upot/Natom
! 					WRITE (*,*) 'Press',         Press
                                END IF

      			CASE('lj')
        			CALL force_neighbor
        			IF(MOD(Nstep,1000) .EQ. 0) THEN
                                        WRITE(*,'(A I4)') 'Calculating LJ force'
                                        WRITE (*,*) 'Potential LJ', Upot/Natom
                                        WRITE (*,*) 'Ukin LJ',      Ukin/Natom
                                        WRITE (*,*) 'Tot LJ',       Ukin/Natom + Upot/Natom
! 					WRITE (*,*) 'Press',        Press
                                END IF	

               END SELECT

! 		IF(MOD(Nstep,100) .EQ. 0) THEN
! 		DO i = 1, Natom
! 		  WRITE (*,*) 'rank', rank, 'timestep', Nstep, '	before', Fx(i),Fy(i)
! 		ENDDO
! 		ENDIF

		Fsend = Fx(1:Natom)
		CALL MPI_ALLREDUCE(Fsend, Fx, Natom, MPI_DOUBLE_PRECISION,&
				   MPI_SUM, MPI_COMM_WORLD, ierr)
		Fsend = Fy(1:Natom)
		CALL MPI_ALLREDUCE(Fsend, Fy, Natom, MPI_DOUBLE_PRECISION,&
				   MPI_SUM, MPI_COMM_WORLD, ierr)

! 		IF(MOD(Nstep,100) .EQ. 0) THEN
! 		    WRITE (*,*) 'before reduce Press',        Press
! 		ENDIF

		PressSend = Press
		CALL MPI_ALLREDUCE(PressSend, Press, 1, MPI_DOUBLE_PRECISION,&
				   MPI_SUM, MPI_COMM_WORLD, ierr)

! 		IF(MOD(Nstep,100) .EQ. 0) THEN
! 		    WRITE (*,*) 'after reduce Press',        Press
! 		ENDIF
		
! 		IF(MOD(Nstep,100) .EQ. 0) THEN
! 		DO i = 1, Natom
! 		  WRITE (*,*)'rank', rank, '	timestep', Nstep, '	after', Fx(i),Fy(i)
! 		ENDDO
! 		ENDIF

    		CALL integrate
    		IF(MOD(Nstep,1000) .EQ. 0) WRITE(*,'(A I4)') 'Integrating'

		IF(MOD(Nstep,1000) .EQ. 0) THEN
		    WRITE (*,*) 'after integrate Press',        Press
		ENDIF
		
    		Nstep = Nstep + 1
	  	END DO
  
!   		DO i=1,Natom
! 			WRITE (*,*) Xx(i),Xo(i)
!   		END DO
  		
	END SUBROUTINE run_simulation

	SUBROUTINE broadcast_velocity
		USE mpi
		USE globals
		INTEGER :: ierr
		DOUBLE PRECISION :: Xsend(NAend-NAstart+1)
		Xsend = Vx(NAstart:NAend)
		CALL MPI_ALLGATHER(Xsend, NAend-NAstart+1, MPI_DOUBLE_PRECISION,&
				   Vx, NAend-NAstart+1, MPI_DOUBLE, MPI_COMM_WORLD, ierr)
		Xsend = Vy(NAstart:NAend)
		CALL MPI_ALLGATHER(Xsend, NAend-NAstart+1, MPI_DOUBLE_PRECISION,&
				   Vy, NAend-NAstart+1, MPI_DOUBLE, MPI_COMM_WORLD, ierr)
	END SUBROUTINE broadcast_velocity
