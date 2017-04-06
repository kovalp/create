FFLAGS = -C
FFLAGS = -O2
FC = pgf90

OBJECTS = gauss_create.o gauss_subs.o gauss_density.o 

gauss_create.x: $(OBJECTS)
	$(FC)  $(OBJECTS) -o gauss_create.x $(LFLAGS)

