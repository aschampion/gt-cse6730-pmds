	SUBROUTINE Force_soft_neighbor

		USE Globals
		IMPLICIT NONE       

      	INTEGER:: I,J,k, Type1, Type2, atom1, atom2
      	REAL(KIND=8):: Dx, Dy, Ff, R_square, Rr, sqR_square,&
			 Part1, Part2
		
		!Initialize the Forces, Potential Energy and Pressure to 0
       	DO I = 1,Natom
         	Fx(I) = 0.0
         	Fy(I) = 0.0
      	END DO
 
	IF(MMov .GT. (0.75D0*Rskin)) THEN
		CALL Neighbour
		MMov = 0.0D0
	END IF

      	Upot = 0.0
      	Press = 0.0
      
		!START LOOP THROUGH ALL ATOM INTERACTIONS
		!DO i=1,natom
		!	WRITE (*,*) i,Xx(i),Yy(i)
		!END DO
		!WRITE (*,*) "i am here 11111111111111111111111"
      	DO I = NAstart,NAend
         	DO k = 1,nlist(i)

      		J = list(I,k)
		IF (I .eq. J) CYCLE 
				
		IF (((I.lt.J) .and. ((MOD(I+J, 2).eq.1))) .or. ((I.gt.J) .and.& 
		    ((MOD(I+J, 2).eq.0)))) THEN
		    Dx = Xx(I) - Xx(J)
		    Dy = Yy(I) - Yy(J)
		
    !            WRITE (*,*) I,J,Xx(I),Xx(J),Dx,Yy(I),Yy(J),Dy
			    
				    !Apply the periodic boundary conditions
		    Dx = Dx - Box*NINT(Dx/Box)
		    Dy = Dy - Box*NINT(Dy/Box)
    
		    R_square = Dx*Dx + Dy*Dy
		    Rr = SQRT(R_square)  
	      
				    !Check If The Distance Is Within The Cutoff Radius for Lennard-Jones Potential
				    !If it is calculate Force and update total force on atom I & J 
		    IF (Rr .LT. Rcut_soft) THEN
			    Part1 = (A_soft*pi)/(Rcut_soft*Rr)
			    Part2 = (sin(pi*Rr/Rcut_soft))
			    Ff    = (Part1*Part2) 
			    Upot  = Upot + A_soft*(1.0+cos(pi*Rr/Rcut_soft))  ! Potential at Rcut_soft = 0
			    Press = Press + Ff
				    
					    !Update the total force on atoms
			    Fx(I) = Fx(I) + Ff*Dx
			    Fy(I) = Fy(I) + Ff*Dy

			    Fx(J) = Fx(J) - Ff*Dx
			    Fy(J) = Fy(J) - Ff*Dy
		  
		    END IF
!WRITE (*,*) "i am here 222222222222222"
		  ENDIF
     		END DO  !END J
	!WRITE(*,*) i,Fx(i),Fy(i)
   		END DO    !END I
   
!   WRITE (*,*) Natom,Nbond

		!Harmonic Bonds             
        !Run through the bond list and grab the interacting atoms
		!WRITE (*,*) " NBond: ",NBond
      	DO k = 1, NBond      	         
            atom1 = BondList(1,k)
            atom2 = BondList(2,k)
			
	    !If neither atom belongs to this process, ignore the bond
	    IF (atom1 .lt. NAstart .AND. atom1 .gt. NAend .AND. &
                atom2 .lt. NAstart .AND. atom2 .gt. NAend) CYCLE

	    IF (((atom1.lt.atom2) .and. ((MOD(atom1+atom2, 2).eq.1))) .or. & 
		((atom1.gt.atom2) .and. ((MOD(atom1+atom2, 2).eq.0)))) THEN
		Dx = Xx(atom1) - Xx(atom2)
		Dy = Yy(atom1) - Yy(atom2)

		!Apply the periodic boundary conditions
		Dx = Dx - Box*nint(Dx/Box)
		Dy = Dy - Box*nint(Dy/Box)


		R_square = Dx*Dx + Dy*Dy

		sqR_square = sqrt(R_square)
		Upot = Upot + K_bond*((sqR_square - Rcut_bond)**2.0)  !The bond potential at Rcut = 0 

		IF (sqR_square .gt. 0.0) THEN
		Ff = (-K_bond*2.0*(sqR_square-Rcut_bond))/sqR_square
		ELSE
		Ff = 0.0
                END IF

		Press = Press + Ff
		Fx(atom1) = Fx(atom1) + Ff*Dx
		Fy(atom1) = Fy(atom1) + Ff*Dy

		Fx(atom2) = Fx(atom2) - Ff*Dx
		Fy(atom2) = Fy(atom2) - Ff*Dy
	  ENDIF

      	END DO
      	
      	!WRITE (*,*) 'Soft Potential + Bond', UPot/Natom
      	
      	!Scale The Pressure
      	Press = Press/(2.0d0*Box*Box)
      
      	!write(*,*) K_bond,Rcut_bond
      
      	!DO i=1,Natom
      	!	WRITE (*,*) i,Fx(i),Fy(i)
	    !END DO
 
	END SUBROUTINE Force_soft_neighbor
