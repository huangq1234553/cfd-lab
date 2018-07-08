INCLUDES = \
	-I. \
	-I$(PRECICE_ROOT)/src/precice

LIBS = \
	-L$(PRECICE_ROOT)/build/last -lprecice\
	-lm

CC = gcc
CFLAGS = -O2 -Wall -pedantic -Werror -Wno-unused-value -Wno-unused-variable -Wno-unused-result -g

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

tests:
	./sim Tests/EmptyCavity/EmptyCavity.dat --compact -q --notemp
	./sim Tests/DrivenCavity/DrivenCavity.dat --compact -q --notemp
	./sim Tests/FlowOverStep/FlowOverStep.dat --compact -q --notemp
	./sim Tests/FluidTrap1/FluidTrap1.dat --compact -q
	./sim Tests/FluidTrap2/FluidTrap2.dat --compact -q
	./sim Tests/KarmanVortex/KarmanVortex.dat --compact -q --notemp
	./sim Tests/NaturalConvection/NaturalConvection.dat --compact -q
	./sim Tests/RayleighBernardConvection/RayleighBenardConvection1.dat --compact -q

debugtest:
	./sim DebugTest/F1HeatExchangerPressureBoundaries/heat_exchanger_Pboundaries.dat --debug

helper.o            : helper.h logger.h
init.o              : helper.h init.h boundary_configurator.h logger.h
boundary_val.o      : helper.h boundary_val.h logger.h
uvp.o               : helper.h uvp.h logger.h
visual.o            : helper.h logger.h
logger.o            : logger.h timing.h
precice_adapter.o 	: precice_adapter.h

main.o        : helper.h init.h boundary_val.h uvp.h visual.h sor.h logger.h boundary_configurator.h timing.h precice_adapter.h
