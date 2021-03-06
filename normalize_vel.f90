SUBROUTINE Normalize_Velocity
  
  USE Globals
  IMPLICIT NONE
      
  REAL(KIND=8) :: Accx,Accy,Adjust
  INTEGER :: i

  Accx = 0.0D0
  Accy = 0.0D0
  Ukin = 0.0D0

  DO i=1,Natom
    Accx = Accx + Vx(I)
    Accy = Accy + Vy(I)
  END DO

  Accx = Accx/REAL(Natom)
  Accy = Accy/REAL(Natom)

  !Shift velocity to get zero momentum
  !Calculate The Kinetic Energy Ukin

  DO I = 1, Natom
    Vx(I) = Vx(I) - Accx
    Vy(I) = Vy(I) - Accy    
    Ukin = Ukin + Vx(I)**2 + Vy(I)**2
  END DO

  !Scale All Velocities To The Correct Temperature

  Adjust = SQRT(Temp_target*REAL(2*Natom-2)/Ukin)
  IF(Natom .EQ. 1) THEN
    Adjust = 1.0D0
  ENDIF
  
  !WRITE (*,*) Adjust

  DO I = 1,Natom
    Vx(I) = Adjust*Vx(I)
    Vy(I) = Adjust*Vy(I)
  END DO
  
END SUBROUTINE Normalize_Velocity
