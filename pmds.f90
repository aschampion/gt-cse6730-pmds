	PROGRAM pmds
  		
  		CHARACTER(LEN=255) :: input_file

  		CALL GETARG(1, input_file)
  		OPEN(UNIT=99, FILE='out.dump')
  		
  		CALL input_parser(input_file)

  		CLOSE(99)
  		
	END PROGRAM pmds

	SUBROUTINE run_simulation(num_timesteps, pair_style)
  		
  		USE Globals
		IMPLICIT NONE
  		
  		INTEGER :: num_timesteps,i
  		CHARACTER(*), INTENT(IN) :: pair_style
  
  		Nstep = 0
  		
  		OPEN(FILE = "Stats.dat",UNIT=35)
  		
  		WRITE (35,*) "Timestep","Ukin","Pot","Total","Pressure"
  
  		DO WHILE(Nstep .LT. num_timesteps)
    		IF(MOD(Nstep,100) .EQ. 0) WRITE(*,'(A I8)') 'Running timestep: ', Nstep
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

    		Nstep = Nstep + 1
	  	END DO
	  	CLOSE(35)
	  	
	  	OPEN(FILE = 'FinalRes.dat',UNIT=27)
	  	WRITE (27,*) "AtomID","Type","Vx","Vy","Fx","Fy"
	  	DO i=1,Natom
	  		WRITE (27,*) i,AT(i),Vx(i),Vy(i),SQRT(Vx(i)**2+Vy(i)**2),Fx(i),Fy(i)
	  	END DO
	  	CLOSE(27)
	  	
	  	
  
	END SUBROUTINE run_simulation
