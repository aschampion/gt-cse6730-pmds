!    Fx/Fy       : Forces
!    Xx/Yy       : Positions
!    Xp/Yp       : Previous positions
!    Vx/Vy       : Velocities
!    Box         : Dimensions of the periodic boundary
!    At          : Atom types

!    dT          : Timestep
C	 EStep		 : No. of timesteps after which to rescale velocity to given given temp.

!    Nstep       : Timestep iteration count
!    Natom       : Number of atoms

!    Ukin        : Kinetic energy
!    Upot        : Potential energy
!    Utot        : Total energy
!    Temp        : Temperature
!    Press       : Pressure

!    Rcut        : N-D array of cutoff radii between atom type pairs
!    Ecut        : N-D array of cutoff energies between atom type pairs

!	 LL			 : Linked list for binned atoms
!    hoc		 : Head of linked list for each cell
!	 nlist		 : No. of neighbour atoms in each verlet list
!	 list		 : Verlet list for each atom
!	 Rskin		 : Skin radius for limiting Verlet list radius

!    Pi          : Set the value to 3.14159265

     Integer Maxatom, MaxNumtypes

     Parameter (Maxatom = 100000)
     Parameter (MaxNumtypes = 10)
     Parameter (Maxcell = 1500)
     Parameter (Maxneigh = 50)

     Double Precision Fx(Maxatom),Fy(Maxatom),Xx(Maxatom),Yy(Maxatom),&
     		      Xp(Maxatom),Yp(Maxatom),Vx(Maxatom),Vy(Maxatom),&
     		      Rcut(MaxNumtypes, MaxNumtypes),&
     		      Ecut(MaxNumtypes, MaxNumtypes),&
     		      Box, Press, Temp, Utot, Upot, Ukin, Pi,Rskin
    
    Integer At(Maxatom), Natom, Nstep, dT,&
    		LL(Maxatom), hoc(Maxcell),nlist(Maxatom),list(Maxatom,Maxneigh)
    
    Common /Syss/    Fx, Fy, Xx, Yy, Xp, Yp, Vx, Vy,&
                     Box, Press, Temp, Utot, Upot, Ukin,&
                     Natom, Nstep, dT

    Parameter (Pi = 3.14159265)
