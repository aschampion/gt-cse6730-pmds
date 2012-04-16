      SUBROUTINE Integrate
      
      	USE Globals
	USE mpi

      	IMPLICIT NONE
 
		!Integrate The Equations Of Motion And calculate The Total Impulse
 
      	INTEGER	:: I,ierr
      	REAL(KIND=8) :: Xxn(Natom),Yyn(Natom),Zzn(Natom)
	DOUBLE PRECISION :: MvelSend

	DO i = NAstart,NAend
      		Xxn(i) = 2.0D0*Xx(i)-Xp(i)+Fx(i)*dt*dt
      		Yyn(i) = 2.0D0*Yy(i)-Yp(i)+Fy(i)*dt*dt
      		
      		Vx(i) = (Xxn(i)-Xp(i))/(2.0D0*dt)
      		Vy(i) = (Yyn(i)-Yp(i))/(2.0D0*dt)
      	END DO
 
		!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		!     For Nstep < Ninit; Use The Velocity Scaling To Get The Exact Temperature  !
		!     Otherwise: Scale = 1. Beware That The Positions/Velocities Have To Be     !
		!     Recalculated !!!!                                                         !
		!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

		IF(MOD(Nstep,Estep) .EQ. 0) THEN
			CALL Broadcast_velocity()
			CALL Normalize_velocity()
		ENDIF
		
		Ukin = 0.0D0

		Mvel = 0.0D0
 
		!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		!     Scale Velocities And Put Particles Back In The Box      !
		!     Beware: The Old Positions Are Also Put Back In The Box  !
		!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		DO i = NAstart,NAend
		
			Ukin = Ukin + 0.5D0*(Vx(i)**2+Vy(i)**2)
		
			Xxn(i) = Xp(i) + 2.0D0*Vx(i)*dt
			Yyn(i) = Yp(i) + 2.0D0*Vy(i)*dt
	
			IF(SQRT(Vx(i)**2+Vy(i)**2) .GT. Mvel) Mvel = SQRT(Vx(i)**2+Vy(i)**2)

			Xp(i) = Xx(i)
			Xx(i) = Xxn(i)
			Yp(i) = Yy(i)
			Yy(i) = Yyn(i)
			
			!Put Particles Back In The Box

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
		
		Press = Press + 2.0d0*Ukin*Dble(Natom)/(Box*Box*Dble(2*Natom-2))

		MvelSend = Mvel
		CALL MPI_ALLREDUCE(MvelSend, Mvel, 1, MPI_DOUBLE_PRECISION,&
				   MPI_MAX, MPI_COMM_WORLD, ierr)

		MMov = MMov + Mvel

! 		IF(MOD(Nstep,100) .EQ. 0) THEN
! 		   WRITE (*,*) 'Press', Press
! 		ENDIF
		
		!WRITE (*,*) dT
		!DO i=1,Natom
		!	WRITE (*,*) Xp(i),Xx(i),Yp(i),Yy(i)
		!END DO
	
	END SUBROUTINE Integrate
