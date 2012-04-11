	SUBROUTINE Initial
		
		IMPLICIT NONE
		INCLUDE 'globals.inc'

		! paramters that will be used

  		INTEGER:: I, seed(3)
  		REAL(KIND=8):: Accx, Accy, Adjust

  		Accx = 0.0d0
  		Accy = 0.0d0
  		Ukin = 0.0d0
  		dT = 0.005D0
  		Estep = 1
  		Box = 35.85686D0

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

  !Call Neighbour

	END SUBROUTINE Initial
