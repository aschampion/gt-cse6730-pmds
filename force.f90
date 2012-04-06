      Subroutine Force(f) 
      Implicit None
 
!     Calculate The Forces And Potential Energy
       
      Include 'globals.inc'
      Double Precision:: f(Maxatom)
      Integer:: I,J,Type1, Type2
      Double Precision:: Dx, Dy, Ff, R_square, R_square_i, R_six_i, Rcut, Rcutsq,&
			 sigma_square, K_bond, Rcut_bond,&
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
! Determine the Lennard-Jones parameters dependent on the two atom types          
             Type1 = AT(I)   !Array list of atom type
             Type2 = AT(J)  
            
             ! Type 1 and 2 will be 1-4 integer of atom type
             Rcut = R_cut_matrix(Type1,Type2)
             Rcutsq = Rcut**2.0
             sigma = sigma_matrix(Type1,Type2)
             eps = epsilon_matrix(Type1,Type2)

! Check If The Distance Is Within The Cutoff Radius for Lennard-Jones Potential
! If it is calculate Force and update total force on atom I & J 
            IF (R_square .Lt. Rcutsq) THEN
               R_square_i = 1.0/R_square
               sigma_square = sigma**2.0
               R_six_i = (R_square_i*sigma_square)**3.0
          
               Upot  = Upot + 4.00*eps*R_six_i*(R_six_i - 1.00) - Ecut
               Ff    = 48.0*eps*R_square_i*R_six_i*(R_six_i - 0.5)

               Press = Press + Ff
               Ff    = Ff*R_square_i
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
!QUESTION: should Energy_bond be Global?? Or, We can specify it here
! DONG: The Bonds I still don't quite understand but Energy_bond, bond_Rcut 
! should be global and defined from the read in of the data file. We 
! probably don't need to reassign them to K-bond, etc as I did here. 
! This is a rough draft. 
              K_bond = Energy_bond
              Rcut_bond = bond_Rcut
              Rcutsq_bond = Rcut_bond**2.0
              IF (R_square .Lt. Rcutsq_bond) THEN
                sqR_square = sqrt(R_square)
                Ff = -K_bond*(2.0-(Rcut_bond/sqR_square))
                Upot = Upot + K_bond*((sqR_square - Rcut_bond)**2.0)  !The bond potential at Rcut = 0 
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
