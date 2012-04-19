      SUBROUTINE Force_Neighbor
      		
      		USE Globals
      		IMPLICIT NONE
       		
      		INTEGER	:: I,J,K,Type1, Type2, atom1, atom2
      		REAL(KIND=8):: Dx, Dy, Ff, R_square, R_square_i, R_six_i, Rcutsq,&
			 sigma_square, R_cut,Rcutsq_bond, eps, sigma, sqR_square    
			
			!Initialize the Forces, Potential Energy and Pressure to 0 
      		DO I = 1,Natom
       			Fx(I) = 0.0
         		Fy(I) = 0.0
      		END DO
       		Upot = 0.0
      		Press = 0.0
      		
      		IF(MMov .GT. (0.75D0*Rskin)) THEN
				CALL Neighbour
				MMov = 0.0D0
			END IF

			!START LOOP THROUGH ALL ATOM INTERACTIONS
      		DO I = NAstart,NAend

      			DO k = 1,nlist(i)
      	!	DO J = 1,Natom
      				J = list(I,k)
      			
      				!IF(J .GT. I) THEN
      	
         		!DO J = (I + 1), Natom
 
			!Calculate the distance between the two atoms
            			Dx = Xx(I) - Xx(J)
            			Dy = Yy(I) - Yy(J)
 
					!Apply the periodic boundary conditions
            			Dx = Dx - Box*NINT(Dx/Box)
            			Dy = Dy - Box*NINT(Dy/Box)
             			R_square = Dx*Dx + Dy*Dy  
					
					!Determine the Lennard-Jones parameters dependent on the two atom types          
             			Type1 = AT(I)   !Array list of atom type
             			Type2 = AT(J)  
            
             		!Type 1 and 2 will be 1-4 integer of atom type
             			R_cut = Rcut(Type1,Type2)
             			Rcutsq = R_cut**2.0
             			sigma = sigma_matrix(Type1,Type2)
             			eps = epsilon_matrix(Type1,Type2)

					!Check If The Distance Is Within The Cutoff Radius for Lennard-Jones Potential
					!If it is calculate Force and update total force on atom I & J 
            			IF (R_square .Lt. Rcutsq) THEN
               				R_square_i = 1.0/R_square
               				sigma_square = sigma**2.0
               				R_six_i = (R_square_i*sigma_square)**3.0
          
               			 	Upot  = Upot + 4.00*eps*R_six_i*(R_six_i - 1.00) - Ecut(Type1,Type2)
               			 	Ff     = 48.0*eps*R_six_i*(R_six_i - 0.5)

               				sqR_square = sqrt(R_square)
										
               				Ff = Ff*R_square_i
               				Press = Press + Ff*sqR_square
						
						!Update the total force on atoms
               				Fx(I) = Fx(I) + Ff*Dx
               				Fy(I) = Fy(I) + Ff*Dy
 
            			END IF
            		!END IF

     			END DO  !END J

   		END DO    !END I


		!Harmonic Bonds             
        !Run through the bond list and grab the interacting atoms
     	DO k = 1, NBond
        	atom1 = BondList(1,k)
            atom2 = BondList(2,k)
 
	    !If neither atom belongs to this process, ignore the bond
	    IF ((atom1 .lt. NAstart .OR. atom1 .gt. NAend) .AND. &
                (atom2 .lt. NAstart .OR. atom2 .gt. NAend)) CYCLE 
            !Calculate the distance between the two atoms
            Dx = Xx(atom1) - Xx(atom2)
            Dy = Yy(atom1) - Yy(atom2)
 
            !Apply the periodic boundary conditions
            Dx = Dx - Box*nint(Dx/Box)
            Dy = Dy - Box*nint(Dy/Box)
 
            R_square = Dx*Dx + Dy*Dy  

            sqR_square = sqrt(R_square)
            Upot = Upot + K_bond*((sqR_square - Rcut_bond)**2.0)  !The bond potential at Rcut = 0 
            IF (sqR_square .gt. 0) THEN
            Ff = (-K_bond*2.0*(sqR_square - Rcut_bond))/sqR_square 
            ELSE
            Ff = 0.0
            END IF
            Press = Press + Ff*sqR_square
            Fx(atom1) = Fx(atom1) + Ff*Dx
            Fy(atom1) = Fy(atom1) + Ff*Dy

            Fx(atom2) = Fx(atom2) - Ff*Dx
            Fy(atom2) = Fy(atom2) - Ff*Dy
      	END DO

      	
!      	WRITE (*,*) 'Potential LJ',Upot/Natom
 
		!Scale The Pressure
      	Press = (Natom*Temp_Target+Press/6.0D0)/(Box*Box)

	!WRITE (*,*) '2222222222 Press LJ',Press
       	
	END SUBROUTINE Force_Neighbor

