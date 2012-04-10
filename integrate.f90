      SUBROUTINE Integrate
      
      	IMPLICIT NONE
      	INCLUDE 'globals.inc'
 
!     Integrate The Equations Of Motion And calculate The Total Impulse
 
      	INTEGER	:: I
      	DOUBLE PRECISION :: Xxn(Natom),Yyn(Natom),Zzn(Natom),scale

!cccccccccccccccccccccccccccccccccccccccccccccccccc 
!     Set Kinetic Energy To Zero                  !
!     Verlet Integrator                           !
!                                                 !
!     Xxx/Yyy/Zzz = New Position (T + Delta T)    !
!     Rxx/Ryy/Rzz = Position     (T          )    !
!     Rxf/Ryf/Rzf = Old Position (T - Delta T)    !
!                                                 !
!     Vxx/Vyy/Vzz = Velocity     (T          )    !
!cccccccccccccccccccccccccccccccccccccccccccccccccc

      	Ukin = 0.0d0
      	
      	DO i=1,natom
      		Xxn(i) = 2.0D0*Xx(i)-Xp(i)+Fx(i)*dt*dt
      		Yyn(i) = 2.0D0*Yy(i)-Yp(i)+Fy(i)*dt*dt
      		
      		Vx(i) = (Xxn(i)-Xp(i))/(2.0D0*dt)
      		Vy(i) = (Yyn(i)-Yp(i))/(2.0D0*dt)
      	
      		Ukin = Ukin + 0.5D0*(Vx(i)**2 + Vy(i)**2)
      	END DO
 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     For Nstep < Ninit; Use The Velocity Scaling To Get The Exact Temperature  !
!     Otherwise: Scale = 1. Beware That The Positions/Velocities Have To Be     !
!     Recalculated !!!!                                                         !
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

		IF(MOD(Nstep,Estep) .EQ. 0) THEN
			scale = DSQRT(Temp_Target*DBLE(2*Natom-2)/(2.0D0*Ukin))
		ELSE
         	scale = 1.0D0
        ENDIF
        
        !WRITE (*,*) scale
 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Scale Velocities And Put Particles Back In The Box      !
!     Beware: The Old Positions Are Also Put Back In The Box  !
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

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
			
!     Put Particles Back In The Box

			IF (Xx(i) .GT. Box) THEN
				Xx(i) = Xx(i) - Box
				Xp(i) = Xp(i) - Box
			ELSE IF (Xx(i) .LT. 0.0D0) THEN
				Xx(i) = Xx(i) + Box				
				Xp(i) = Xp(i) + Box
			END IF
			
			IF (Yy(i) .GT. Box) THEN
				Yy(i) = Yy(i) - Box
				Yp(i) = Yp(i) - Box
			ELSE IF (Yy(i) .LT. 0.0D0) THEN
				Yy(i) = Yy(i) + Box
				Yp(i) = Yp(i) + Box
			END IF
			
		END DO
		
		!WRITE (*,*) dT
		!DO i=1,Natom
		!	WRITE (*,*) Xp(i),Xx(i),Yp(i),Yy(i)
		!END DO
	
	END SUBROUTINE Integrate
