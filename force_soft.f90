      Subroutine Force_soft()
      Implicit None
 
!     Calculate The Forces And Potential Energy
       
      Include 'globals.inc'

      Integer:: I,J,Type1, Type2, atom1, atom2
      Double Precision:: Dx,Dy,Ff,R_square,Rr,sqR_square
      DOUBLE PRECISION:: Part1, Part2,Press, Upot
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
            Dx = Xx(I) - Xx(J)
            Dy = Yy(I) - Yy(J)
 
! Apply the periodic boundary conditions
            Dx = Dx - Box*nint(Dx/Box)
            Dy = Dy - Box*nint(Dy/Box)
 
            R_square = Dx*Dx + Dy*Dy
            Rr = sqrt(R_square)  
! Determine the Lennard-Jones parameters dependent on the two atom types             
             Type1 = AT(I)   !Array list of atom type
             Type2 = AT(J)  
            
! Check If The Distance Is Within The Cutoff Radius for Lennard-Jones Potential
! If it is calculate Force and update total force on atom I & J 
            IF (Rr .Lt. Rcut_soft) THEN
               Part1 = (A_soft*pi)*(Rcut_soft*Rr)
               Part2 = (sin(pi*Rr/Rcut_soft))
               Ff    = (Part1*Part2) 
               Upot  = Upot + A_soft*(1.0+cos(pi*Rr/Rcut_soft))  ! Potential at Rcut_soft = 0
               Press = Press + Ff
!  Update the total force on atoms
               Fx(I) = Fx(I) + Ff*Dx
               Fy(I) = Fy(I) + Ff*Dy
 
               Fx(J) = Fx(J) - Ff*Dx
               Fy(J) = Fy(J) - Ff*Dy
 
            END IF


     END DO  !END J
   END DO    !END I

!*************** Harmonic Bonds ************************************************
! ROUGH DRAFT OF USING BONDS!
!Harmonic Bonds             
          !Run through the bond list and grab the interacting atoms
      DO k = 1, MaxBonds
              atom1 = BondList(1,k)
              atom2   BondList(2,k)
                sqR_square = sqrt(R_square)
                Ff = -K_bond*(2.0-(Rcut_bond/sqR_square))
                Upot = Upot + K_bond*((sqR_square - Rcut_bond)**2.0)  !The bond potential at Rcut = 0 
                Press = Press + Ff
                Fx(atom1) = Fx(atom1) + Ff*Dx
                Fy(atom1) = Fy(atom1) + Ff*Dy

                Fx(atom2) = Fx(atom2) - Ff*Dx
                Fy(atom2) = Fy(atom2) - Ff*Dy

      END DO
 
!    Scale The Pressure
 
      Press = Press/(3.0d0*Box*Box*Box)
 
      Return
      END SUBROUTINE
