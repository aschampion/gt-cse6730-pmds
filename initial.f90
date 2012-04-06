Subroutine Initial
Implicit None



Include 'globals.inc'

! paramters that will be used

  Integer:: I, seed(3)
  Double Precision:: Accx, Accy, Adjust

  Accx = 0.0d0
  Accy = 0.0d0
  Ukin = 0.0d0

  Call itime(seed)

  Call srand(seed(3))

!    Generate Velocities From A Gaussian; Set Impulse To Zero

  Do I = 1, Natom
    Vx(I) = rand() - 0.5
    Vy(I) = rand() - 0.5
    Accx = Accx + Vx(I)
    Accy = Accy + Vy(I)
  Enddo

  Accx = Accx/Dble(Natom)
  Accy = Accy/Dble(Natom)

!     Calculate The Kinetic Energy Ukin

  Do I = 1, Natom
    Vx(I) = Vx(I) - Accx
    Vy(I) = Vy(I) - Accy

    Ukin = Ukin + Vx(I)*Vx(I) + Vy(I)*Vy(I)
  Enddo

!     Scale All Velocities To The Correct Temperature

  Adjust = Dsqrt(Temp*Dble(2*Natom-2)/(2.0d0*Ukin))
  
  Do I = 1,Natom
    Vx(I) = Adjust*Vx(I)
    Vy(I) = Adjust*Vy(I)
  Enddo

!     Calculate Previous Position Using The Generated Velocity

  Do I = 1,Natom
    Xp(I) = Vx(I) - dT*Vx(I)
    Yp(I) = Vy(I) - dT*Vy(I)
  Enddo

  Call Neighbour

Return
End

