SUBROUTINE Force
    
  USE Globals
  IMPLICIT NONE
    
  INTEGER :: I,J,K,Type1, Type2, atom1, atom2
  REAL(KIND=8):: Dx, Dy, Ff, R_square, R_square_i, R_six_i, Rcutsq, atom_cout,&
  sigma_square, R_cut, eps, sigma, sqR_square, Press_sub    

  !Initialize the Forces, Potential Energy and Pressure to 0 
  DO I = 1,Natom
    Fx(I) = 0.0
    Fy(I) = 0.0
  END DO
  Upot = 0.0
  Press = 0.0

  ! new
  atom_cout = 0
  Press_sub = 0.0

  !START LOOP THROUGH ALL ATOM INTERACTIONS
  DO I = 1, (Natom - 1)
    ! clear the press_sub and atom_cout
    Press_sub = 0.0
    atom_cout = 0.0

    DO J = (I + 1), Natom

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
        !Press = Press - Ff*sqR_square
        Press_sub = Press_sub + Ff*sqR_square
        atom_cout = atom_cout + 1

        Ff = Ff*R_square_i

        !Update the total force on atoms
        Fx(I) = Fx(I) + Ff*Dx
        Fy(I) = Fy(I) + Ff*Dy

        Fx(J) = Fx(J) - Ff*Dx
        Fy(J) = Fy(J) - Ff*Dy

      END IF

    END DO  !END J

    IF (atom_cout.NE.0) THEN
      !WRITE (*,*) '11111111 Press LJ',Press_sub, atom_cout
       Press_sub = (Press_sub/atom_cout)
       Press = Press + Press_sub
    ENDIF

  END DO    !END I

  !WRITE (*,*) '11111111 Press LJ',Press

  !Harmonic Bonds             
  !Run through the bond list and grab the interacting atoms
  Press_sub = 0.0
  atom_cout = 0
  DO k = 1, NBond
    atom1 = BondList(1,k)
    atom2 = BondList(2,k)

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
    Press_sub = Press_sub + Ff*sqR_square
    atom_cout = atom_cout + 1
    Fx(atom1) = Fx(atom1) + Ff*Dx
    Fy(atom1) = Fy(atom1) + Ff*Dy

    Fx(atom2) = Fx(atom2) - Ff*Dx
    Fy(atom2) = Fy(atom2) - Ff*Dy
  END DO

  Press_sub =  Press_sub/atom_cout
      
  !WRITE (*,*) 'Potential LJ',Upot/Natom

  !Scale The Pressure
  Press = (Press + Press_sub)/(2.0d0*Box*Box)

  !WRITE (*,*) '2222222222 Press LJ',Press
      
END SUBROUTINE Force

