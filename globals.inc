!    Fx/Fy       : Forces
!    Xx/Yy       : Positions
!    Xp/Yp       : Previous positions
!    Vx/Vy       : Velocities
!    Box         : Dimensions of the periodic boundary
!    At          : Atom types

!    dT          : Timestep
!	 EStep		 : No. of timesteps after which to rescale velocity to given given temp.

!    Nstep       : Timestep iteration count
!    Natom       : Number of atoms
!    Nbond       : Number of bonds

!    Ukin        : Kinetic energy
!    Upot        : Potential energy
!    Utot        : Total energy
!    Temp        : Temperature
!    Press       : Pressure

!    Rcut        : matrix of cutoff radii between atom type pairs
!    Ecut        : matrix array of cutoff energies between atom type pairs
!    sigma_matrix: matrix of potential sigmas between atom type pairs
!    epsilon_matrix: matrix of potential epsilons between atom type pairs
!    A_soft      : A parameter for soft potential
!    Rcut_soft   : radius cutoff for soft potential
!    K_bond      :
!    Rcut_bond      : radius cutoff for bonds

!	 LL			 : Linked list for binned atoms
!    hoc		 : Head of linked list for each cell
!	 nlist		 : No. of neighbour atoms in each verlet list
!	 list		 : Verlet list for each atom
!	 Rskin		 : Skin radius for limiting Verlet list radius

!    Pi          : Set the value to 3.14159265

     Integer Maxatom, MaxNumtypes, MaxBonds, Maxcell, Maxneigh

     Parameter (Maxatom = 100000)
     Parameter (MaxNumtypes = 10)
     Parameter (MaxBonds = 1000)
     Parameter (Maxcell = 1500)
     Parameter (Maxneigh = 50)

     Double Precision Fx(Maxatom),Fy(Maxatom),Xx(Maxatom),Yy(Maxatom),&
     		      Xp(Maxatom),Yp(Maxatom),Vx(Maxatom),Vy(Maxatom),&
     		      Rcut(MaxNumtypes, MaxNumtypes),&
     		      Ecut(MaxNumtypes, MaxNumtypes),&
                      sigma_matrix(MaxNumtypes, MaxNumtypes),&
                      epsilon_matrix(MaxNumtypes, MaxNumtypes),&
     		      Box, Press, Temp, Utot, Upot, Ukin, Pi,Rskin,&
                      K_bond, Rcut_bond, A_soft, Rcut_soft, Temp_target
    
    Integer  At(Maxatom), Natom, Nbond, Nstep,&
    		LL(Maxatom), hoc(Maxcell),nlist(Maxatom),list(Maxatom,Maxneigh),&
                bondlist(2, MaxBonds), dT
    
    Common /Syss/    Fx, Fy, Xx, Yy, Xp, Yp, Vx, Vy,&
                     Rcut, Ecut, sigma_matrix, epsilon_matrix,&
                     Box, Press, Temp, Utot, Upot, Ukin, Rskin,&
                     K_bond, Rcut_bond, A_soft, Rcut_soft, Temp_target,&
                     At, Natom, Nbond, Nstep, LL, hoc, nlist, list, bondlist, dt

    Parameter (Pi = 3.14159265)
