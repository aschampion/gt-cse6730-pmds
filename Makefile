# FC     = gfortran-mp-4.8
FC = gfortran
FFLAGS =-O
#FFLAGS =-fno-align-comments -O

OBJS   = globals.o force.o force_soft.o initial.o integrate.o neighbour.o parser.o norm_vel.o

globals.o : globals.f90
	$(FC) $(FFLAGS) -o $@ -c $<

force.o : force.f90
	$(FC) $(FFLAGS) -o $@ -c $<

force_soft.o : force_soft.f90
	$(FC) $(FFLAGS) -o $@ -c $<

initial.o : initial.f90
	$(FC) $(FFLAGS) -o $@ -c $<

integrate.o : integrate.f90
	$(FC) $(FFLAGS) -o $@ -c $<

neighbour.o : neighbour.f90
	$(FC) $(FFLAGS) -o $@ -c $<

parser.o : parser.f90
	$(FC) $(FFLAGS) -o $@ -c $<
	
norm_vel.o : normalize_vel.f90
	$(FC) $(FFLAGS) -o $@ -c $<

pmds : pmds.f90 $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

all : pmds

clean :
	rm -f *.o
	rm -f pmds
