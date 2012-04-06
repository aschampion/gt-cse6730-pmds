!    Fx/Fy       : Forces
!    Xx/Yy       : Positions
!    Xp/Yp       : Previous positions
!    Vx/Vy       : Velocities
!    Box         : Dimensions of the periodic boundary
!    At          : Atom types

!    dT          : Timestep

!    Nstep       : Timestep iteration count
!    Natom       : Number of atoms

!    Ukin        : Kinetic energy
!    Upot        : Potential energy
!    Utot        : Total energy
!    Temp        : Temperature
!    Press       : Pressure

!    R_cut_matrix        : N-D array of cutoff radii between atom type pairs
!    Ecut        : N-D array of cutoff energies between atom type pairs

!    Pi          : Set the value to 3.14159265

     Integer Maxatom, MaxNumtypes

     Parameter (Maxatom = 100000)
     Parameter (MaxNumtypes = 10)

     Double Precision Fx(Maxatom),Fy(Maxatom),Xx(Maxatom),Yy(Maxatom),&
     		      Xp(Maxatom),Yp(Maxatom),Vx(Maxatom),Vy(Maxatom),&
     		      R_cut_matrix(MaxNumtypes, MaxNumtypes),&
     		      E_cut_matrix(MaxNumtypes, MaxNumtypes),&
		      sigma_matrix(MaxNumtypes, MaxNumtypes),&
		      epsilon_matrix(MaxNumtypes, MaxNumtypes),&
     		      Box, Press, Temp, Utot, Upot, Ukin, Pi, Ecut,&
		      Energy_bond, bond_Rcut
    
    Integer At(Maxatom), Natom, Nstep, dT
    
    Common /Syss/    Fx, Fy, Xx, Yy, Xp, Yp, Vx, Vy,&
                     Box, Press, Temp, Utot, Upot, Ukin,&
		     R_cut_matrix, E_cut_matrix,&
		     sigma_matrix, epsilon_matrix,&
                     Natom, Nstep, dT

    Parameter (Pi = 3.14159265)
