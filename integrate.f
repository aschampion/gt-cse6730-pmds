      SUBROUTINE Integrate
      
      	IMPLICIT NONE
      	INCLUDE 'globals.inc'
 
C     Integrate The Equations Of Motion And Calculate The Total Impulse
 
      	INTEGER	:: I
      	DOUBLE PRECISION :: Xxn(Natom),Yyn(Natom),Zzn(Natom),scale

Ccccccccccccccccccccccccccccccccccccccccccccccccccc 
C     Set Kinetic Energy To Zero                  C
C     Verlet Integrator                           C
C                                                 C
C     Xxx/Yyy/Zzz = New Position (T + Delta T)    C
C     Rxx/Ryy/Rzz = Position     (T          )    C
C     Rxf/Ryf/Rzf = Old Position (T - Delta T)    C
C                                                 C
C     Vxx/Vyy/Vzz = Velocity     (T          )    C
Ccccccccccccccccccccccccccccccccccccccccccccccccccc

      	Ukin = 0.0d0
      	
      	DO i=1,natom
      		Xxn(i) = 2.0D0*Xx(i)-Xp(i)+Fx(i)*dt*dt
      		Yyn(i) = 2.0D0*Yy(i)-Yp(i)+Fy(i)*dt*dt
      		
      		Vx(i) = (Xxn(i)-Xp(i))/(2.0D0*dt)
      		Vy(i) = (Yyn(i)-Yp(i))/(2.0D0*dt)
      	
      		Ukin = Ukin + 0.5D0*(Vx(i)**2 + Vy(i)**2)
      	END DO
 
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     For Nstep < Ninit; Use The Velocity Scaling To Get The Exact Temperature  C
C     Otherwise: Scale = 1. Beware That The Positions/Velocities Have To Be     C
C     Recalculated !!!!                                                         C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

		IF(MOD(Step,Estep) .EQ. 0) THEN
			scale = DSQRT(Temp*DBLE(2*Natom-2)/(2.0D0*Ukin))
		ELSE
         	scale = 1.0D0
        ENDIF
 
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Scale Velocities And Put Particles Back In The Box      C
C     Beware: The Old Positions Are Also Put Back In The Box  C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

		Ukin = 0.0D0

		DO i=1,natom
			Vx(i) = scale*Vx(i)
			Vy(i) = scale*Vy(i)
			
			Xxn(i) = Xp(i) + 2.0D0*Vx(i)*dt
			Yyn(i) = Yp(i) + 2.0D0*Vy(i)*dt
			
			Ukin = Ukin + Vx(i)**2+Vy(i)**2
			
			Xp(i) = Xx(i)
			Xx(i) = Xxn(i)
			Yp(i) = Yy(i)
			Yy(i) = Yyn(i)
			
C     Put Particles Back In The Box

			IF (Xx(i) .GT. Box) THEN
				Xx(i) = Xx(i) - Box
				Xp(i) = Xp(i) - Box
			ELSE IF (Xx(i) .LT. 0.0D0) THEN
				Xx(i) = Xx(i) + Box
				Yy(i) = Yy(i) + Box
			END IF
			
			IF (Yy(i) .GT. Box) THEN
				Yy(i) = Yy(i) - Box
				Yp(i) = Yp(i) - Box
			ELSE IF (Yy(i) .LT. 0.0D0) THEN
				Yy(i) = Yy(i) + Box
				Yp(i) = Yp(i) + Box
			END IF
			
		END DO
		
	END SUBROUTINE Integrate
