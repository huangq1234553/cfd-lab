INCLUDES = \
	-I. \
	-I$(PRECICE_ROOT)/src/precice

LIBS = \
	-L$(PRECICE_ROOT)/build/last -lprecice\
	-lm

CC = gcc
CFLAGS = -Wall -pedantic -Werror $(INCLUDES)
.c.o:  ; $(CC) -c $(CFLAGS) $<

OBJ = 	helper.o\
      	init.o\
      	boundary_val.o\
      	uvp.o\
      	sor.o\
      	main.o\
      	visual.o\
      	logger.o\
      	boundary_configurator.o\
      	timing.o\
      	precice_adapter.o


all:  $(OBJ)
	$(CC) $(CFLAGS) -o sim $(OBJ) $(LIBS)

%.o : %.c
	$(CC) -c $(CFLAGS) $*.c -o $*.o $(LIBS)

clean:
	rm $(OBJ)

runsolidplate:
	./sim Cases/heated_plate/heated-plate.dat -q -o Solid_plate/Out

runsolidconvection:
	./sim Cases/natural_convection/convection.dat -q -o Solid_convection/Out

runsolidexchangeF1:
	./sim Cases/heat_exchange/F1-heat-exchange.dat -q -o Solid_exchange/OutF1

runsolidexchangeF2:
	./sim Cases/heat_exchange/F2-heat-exchange.dat -q -o Solid_exchange/OutF2

helper.o            : helper.h logger.h
init.o              : helper.h init.h boundary_configurator.h logger.h
boundary_val.o      : helper.h boundary_val.h logger.h
uvp.o               : helper.h uvp.h logger.h
visual.o            : helper.h logger.h
logger.o            : logger.h timing.h
precice_adapter.o 	: precice_adapter.h

main.o        : helper.h init.h boundary_val.h uvp.h visual.h sor.h logger.h boundary_configurator.h timing.h precice_adapter.h
