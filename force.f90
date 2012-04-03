      Subroutine Force
      Implicit None
 
!     Calculate The Forces And Potential Energy
       
      Include 'globals.inc'

      Integer I,J
      Double Precision Dx,Dy,Ff,R_square,R_square_i,R_six_i, Rcutsq
 
!**********Initialize the Forces, Potential Energy and Pressure to 0************
 
      DO I = 1,Natom
         Fx(I) = 0.0
         Fy(I) = 0.0
      END DO
 
      Upot = 0.0
      Press = 0.0
!**********Initialize the Forces, Potential Energy and Pressure to 0************
 
! **** START LOOP THROUGH ALL ATOM INTERACTIONS ******************************** 
      DO I = 1,Natom - 1
         DO J = I + 1,Natom
 
! Calculate the distance between the two atoms
            Dx = Rxx(I) - Rxx(J)
            Dy = Ryy(I) - Ryy(J)
 
! Apply the periodic boundary conditions
            Dx = Dx - Box*nint(Dx/Box)
            Dy = Dy - Box*nint(Dy/Box)
 
            R_square = Dx*Dx + Dy*Dy  
!NEED TO ADD IN IF STATEMENT TO DETERMINE THE ATOM INTERACTION
! & APPROPRIATE PARAMETERS FOR FORCE


!     Check If The Distance Is Within The Cutoff Radius for Lennard-Jones Potential
 
            IF (R_square .Lt. Rcutsq) THEN
               R_square_i = 1.0/R_square
               R_six_i = R_square_i**3.0
 
               Upot  = Upot + 4.00*R_six_i*(R_six_i - 1.00) - Ecut
               Ff    = 48.0*R_six_i*(R_six_i - 0.5)
               Press = Press + Ff
               Ff    = Ff*R_square_i
!  Update the total force on atoms
               Fx(I) = Fx(I) + Ff*Dx
               Fy(I) = Fy(I) + Ff*Dy
 
               Fx(J) = Fx(J) - Ff*Dx
               Fy(J) = Fy(J) - Ff*Dy
 
            END IF
         END DO
      END DO
 
!    Scale The Pressure
 
      Press = Press/(3.0d0*Box*Box*Box)
 
      Return
      END SUBROUTINE
