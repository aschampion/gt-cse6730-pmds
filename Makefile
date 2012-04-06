# FC     = gfortran-mp-4.8
FC = gfortran
FFLAGS = -O

OBJS   = initial.o force.o force_soft.o parser.o

initial.o : initial.f90
	$(FC) $(FFLAGS) -o $@ -c $<

force.o : force.f90
	$(FC) $(FFLAGS) -o $@ -c $<

force_soft.o : force_soft.f90
	$(FC) $(FFLAGS) -o $@ -c $<

parser.o : parser.f95
	$(FC) $(FFLAGS) -o $@ -c $<

pmds : pmds.f95 $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

all : pmds

clean :
	rm -f *.o
	rm -f pmds
