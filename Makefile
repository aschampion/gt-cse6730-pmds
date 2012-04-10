# FC     = gfortran-mp-4.8
FC = gfortran
FFLAGS = -fno-align-commons -O

OBJS   = force.o force_soft.o initial.o integrate.o neighbour.o parser.o

force.o : force.f90
	$(FC) $(FFLAGS) -o $@ -c $<

force_soft.o : force_soft.f90
	$(FC) $(FFLAGS) -o $@ -c $<

initial.o : initial.f90
	$(FC) $(FFLAGS) -o $@ -c $<

integrate.o : integrate.f95
	$(FC) $(FFLAGS) -o $@ -c $<

neighbour.o : neighbour.f90
	$(FC) $(FFLAGS) -o $@ -c $<

parser.o : parser.f95
	$(FC) $(FFLAGS) -o $@ -c $<

pmds : pmds.f95 $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

all : pmds

clean :
	rm -f *.o
	rm -f pmds
