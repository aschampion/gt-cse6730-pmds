Subroutine Initial
Implicit None



Include 'globals.inc'

! paramters that will be used

  Integer:: I, seed(3)
  Double Precision:: Accx, Accy, Adjust

  Accx = 0.0d0
  Accy = 0.0d0
  Ukin = 0.0d0
  dT = 0.005D0

  Call itime(seed)

  Call srand(seed(3))

!    Generate Velocities From A Gaussian; Set Impulse To Zero

  Do I = 1, Natom
    Vx(I) = rand(0) - 0.5
    Vy(I) = rand(0) - 0.5
    Accx = Accx + Vx(I)
    Accy = Accy + Vy(I)
 Enddo

  Accx = Accx/Dble(Natom)
  Accy = Accy/Dble(Natom)

!	Shift velocity to get zero momentum
!     Calculate The Kinetic Energy Ukin

  Do I = 1, Natom
    Vx(I) = Vx(I) - Accx
    Vy(I) = Vy(I) - Accy

    Ukin = Ukin + Vx(I)*Vx(I) + Vy(I)*Vy(I)
  Enddo

!     Scale All Velocities To The Correct Temperature

  Adjust = Dsqrt(Temp_target*Dble(2*Natom-2)/(2.0d0*Ukin))
  
  Do I = 1,Natom
    Vx(I) = Adjust*Vx(I)
    Vy(I) = Adjust*Vy(I)
  Enddo

!     Calculate Previous Position Using The Generated Velocity

  Do I = 1,Natom
    Xp(I) = Xx(I) - dT*Vx(I)
    Yp(I) = Yy(I) - dT*Vy(I)           
    !WRITE (*,*) i,Xx(i),Yy(i),dT*Vx(i),dT*Vy(i)    
  Enddo

  !Call Neighbour

Return
End

