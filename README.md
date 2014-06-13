Parallel molecular dynamics simulation of micelle formation
===============

This is a Fortran 90/MPI implementation of a 2D PMDS with partial LAMMPS descriptor compatibility, originally for Georgia Tech's CSE 6730 (Spring 2012). 

The original authors are
* Andrew Champion (myself)
* Mary C. Benage
* Mohan Rajendran
* Zhenjiang Dong

To try it yourself with an MPI-enabled `gfortran`:

```sh
git clone git@github.com:aschampion/gt-cse6730-pmds.git
cd gt-cse6730-pmds
make all
./pmds # Note that this will take some time and produce a ~200MB dump file
```

With the simulation results in `out.dump`, you can now open the Processing visualization in `pmds_visualizer/pmds_visualizer.pde`. Note that the visualizer requires the [controlP5](http://www.sojamo.de/libraries/controlP5/) Processing library, and the paths in `pmds_visualizer/visualizer.properties` should resolve to the micelle simulation description file in the repository and `out.dump` file you generated above.

## License

This code is provided under the MIT license, though additional restrictions may apply to LAMMPS descriptor files or specific parallel algorithms implemented.