      Subroutine Force(f)
      Implicit None
 
!     Calculate The Forces And Potential Energy
       
      Include 'globals.inc'

      Double Precision:: f(Maxatom)
      Integer:: I,J,Type1, Type2
      Double Precision:: Dx, Dy, Ff, R_square, R_square_i, R_six_i, Rcut, Rcutsq,&
			 sigma_square, K_bond, Rcut_bond,Part1, Part2, A_soft&
			 Rcutsq_bond, eps, sigma, sqR_square 
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


! ROUGH DRAFT OF USING BONDS!
                ! HEAD          TAIL1         TAIL2
            IF (Type1==2 .AND. Type2==3 .OR. Type2==4) THEN
             !A bond calculation required
              K_bond = Energy_bond
              Rcut_bond = bond_Rcut
              Rcutsq_bond = Rcut_bond**2.0
              IF (R_square .Lt. Rcutsq_bond) THEN
                sqR_square = sqrt(R_square)
                Ff = -K_bond*(2.0-(Rcut_bond/sqR_square))
                Upot = Upot + K_bond*((sqR_square - Rcut_bond)**2.0)   
                Press = Press + Ff
                Fx(I) = Fx(I) + Ff*Dx
                Fy(I) = Fy(I) + Ff*Dy
 
                Fx(J) = Fx(J) - Ff*Dx
                Fy(J) = Fy(J) - Ff*Dy
             END IF
            END IF

         END DO
      END DO
 
!    Scale The Pressure
 
      Press = Press/(3.0d0*Box*Box*Box)
 
      Return
      END SUBROUTINE
