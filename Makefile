CC=
CFLAGS= 
F90=gfortran
#F90FLAGS=-g -fbounds-check -cpp
F90FLAGS=-g -cpp -O3
PARAFLAGS=-fopenmp

SOURCES= variables.f90 \
         library.f90 \
         defect.f90 \
         prepare_main.f90


OBJ=$(addsuffix .o, $(basename $(SOURCES)))


.SUFFIXES :.c .f90

.f90.o:
	$(F90) $(F90FLAGS) $(PARAFLAGS) -c  $< 

.c.o:
	$(CC) $(CFLAGS) -c $<


fss: $(OBJ)
	@echo "Finding first shell structure ... "
	@echo "Current objects: $(OBJ)"
	$(F90) $(F90FLAGS) $(PARAFLAGS) $(OBJ) -o $@

clean:	
	@echo "Cleaning Directory ... "
	rm -f $(OBJ) 
 
veryclean: 
	@echo "Cleaning Directory AND Executables ..."
	rm -f $(OBJ) Hbond *.mod 
