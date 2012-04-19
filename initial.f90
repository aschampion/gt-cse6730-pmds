	SUBROUTINE Initial
		
		USE Globals
		IMPLICIT NONE

		! paramters that will be used

  		INTEGER:: I,J, seed(3)
  		REAL(KIND=8):: Accx, Accy, Adjust

  		Accx = 0.0d0
  		Accy = 0.0d0
  		Ukin = 0.0d0
  		dT = 0.005D0
  		Estep = 100
  		Cstep = 5
  		! Box = 35.85686D0
  		Rskin = 0.3D0
  		MMov = Rskin

  		CALL ITIME(seed)

  		CALL SRAND(seed(3))
		
		!Generate velocity based on uniform distribution
  		DO I = 1, Natom
    		Vx(I) = RAND(0) - 0.5
    		Vy(I) = RAND(0) - 0.5
  		END DO
  		
  		CALL Normalize_Velocity()

		!Calculate Previous Position Using The Generated Velocity
  		DO I = 1,Natom
    		Xp(I) = Xx(I) - dT*Vx(I)
    		Yp(I) = Yy(I) - dT*Vy(I)           
    		!WRITE (*,*) i,Xx(i),Yy(i),dT*Vx(i),dT*Vy(i)    
		END DO
	
                ! Calculate Ecut for different atom types
                 DO I =1 , MaxNumtypes
                       DO J = 1, MaxNumtypes
                           Ecut(I,J) = 4.0*epsilon_matrix(I,J)*   & 
                           (((sigma_matrix(I,J)/Rcut(I,J))**12.0) - ((sigma_matrix(I,J)/Rcut(I,J))**6.0))
                       END DO
                 END DO
	END SUBROUTINE Initial
