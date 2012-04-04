C    Fx/Fy       : Forces
C    Xx/Yy       : Positions
C    Xp/Yp       : Previous positions
C    Vx/Vy       : Velocities
C    Box         : Dimensions of the periodic boundary

C    dT          : Timestep

C    Nstep       : Timestep iteration count
C    Natom       : Number of atoms

C    Ukin        : Kinetic energy
C    Upot        : Potential energy
C    Utot        : Total energy
C    Temp        : Temperature
C    Press       : Pressure

C    Rcut        : N-D array of cutoff radii between atom type pairs
C    Ecut        : N-D array of cutoff energies between atom type pairs

C    Pi          : Set the value to 3.14159265

     Integer Maxatom, MaxNumtypes

     Parameter (Maxatom = 100000)
     Parameter (MaxNumtypes = 10)

     Double Precision Fx(Maxatom),Fy(Maxatom),Xx(Maxatom),Yy(Maxatom)
    $		      Xp(Maxatom),Yp(Maxatom),Vx(Maxatom),Vy(Maxatom)
    $		      Rcut(MaxNumtypes, MaxNumtypes)
    $		      Ecut(MaxNumtypes, MaxNumtypes)
    $		      Box, Press, Temp, Utot, Upot, Ukin,Pi 
    
    Integer Natom, Nstep, dT
    
    Common /Syss/    Fx, Fy, Xx, Yy, Xp, Yp, Vx, Vy,
   &                 Box, Press, Temp, Utot, Upot, Ukin
   &                 Natom, Nstep, dT

   Pi = 3.14159265
