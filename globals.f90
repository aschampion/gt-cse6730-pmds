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

	MODULE Globals

		IMPLICIT NONE
		SAVE

    	INTEGER,PARAMETER :: Maxatom=100000
    	INTEGER,PARAMETER :: MaxNumtypes=10
    	INTEGER,PARAMETER :: MaxBonds=1000
    	INTEGER,PARAMETER :: Maxcell=1500
    	INTEGER,PARAMETER :: Maxneigh=50
    	REAL(KIND=8),PARAMETER :: Pi = 3.14159265
    
    	REAL(KIND=8) :: Fx(Maxatom),Fy(Maxatom),Xx(Maxatom),Yy(Maxatom),&
     		      Xp(Maxatom),Yp(Maxatom),Vx(Maxatom),Vy(Maxatom),&
     		      Rcut(MaxNumtypes, MaxNumtypes),&
     		      Ecut(MaxNumtypes, MaxNumtypes),&
                      sigma_matrix(MaxNumtypes, MaxNumtypes),&
                      epsilon_matrix(MaxNumtypes, MaxNumtypes),&
     		      Box, Press, Temp, Utot, Upot, Ukin,Rskin, Mvel,&
                      K_bond, Rcut_bond, A_soft, Rcut_soft, Temp_target, dT, MMov
   		
   		INTEGER	::	At(Maxatom), Natom, Nbond, Nstep,&
    		LL(Maxatom), hoc(Maxcell,Maxcell),nlist(Maxatom),list(Maxatom,Maxneigh),&
                bondlist(2, MaxBonds),Estep,NAstart,NAend
	END MODULE Globals
