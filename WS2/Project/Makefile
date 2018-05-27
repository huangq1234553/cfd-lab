CC = mpicxx
CFLAGS = -Wall -pedantic -Werror -Wno-write-strings -ggdb
.c.o:  ; $(CC) -c $(CFLAGS) $<

OBJ = 	helper.o\
      	init.o\
      	boundary_val.o\
      	uvp.o\
      	sor.o\
      	visual.o\
      	logger.o\
      	timing.o\
        parallel.o

SIM_OBJ = main.o
TEST_OBJ = 	 Tests/tests.o\
			 Tests/testing.o\
			 Tests/getProcessCoordinates.test.o\
			 Tests/getProcessNeighbours.test.o\
			 Tests/uvComm.test.o\
			 Tests/pressureComm.test.o
ALL_OBJ = $(OBJ) $(SIM_OBJ) $(TEST_OBJ)

TEST_OUT_DIR = Tests/Out

.PHONY: sim tests clean-core clean-sim clean-tests

all: sim tests

sim:  $(OBJ) $(SIM_OBJ)
	$(CC) $(CFLAGS) -o sim $(OBJ) $(SIM_OBJ) -lm

tests: $(OBJ) $(TEST_OBJ)
	$(CC) $(CFLAGS) -o Tests/tests $(OBJ) $(TEST_OBJ) -lm

runsim:
	mkdir -p Viz && cp problem.dat Viz/ && mpirun -np 4 ./sim Viz/problem

runsimseq:
	mkdir -p Viz && cp problem_seq.dat Viz/ && mpirun -np 1 ./sim Viz/problem_seq

runtests:
	$(MAKE) -C Tests/ all

%.o : %.c
	$(CC) -c $(CFLAGS) $*.c -o $*.o

clean-core:
	rm $(OBJ)

clean-sim:
	rm $(SIM_OBJ)

clean-viz:
	rm Viz/*.vtk

clean-tests:
	rm $(TEST_OBJ)

clean-tests-output:
	rm $(TEST_OUT_DIR)/*
	rmdir $(TEST_OUT_DIR)

clean: clean-core clean-sim clean-tests

helper.o      : helper.h logger.h
init.o        : helper.h init.h logger.h
boundary_val.o: helper.h boundary_val.h logger.h
uvp.o         : helper.h uvp.h logger.h
visual.o      : helper.h logger.h visual.h
logger.o      : helper.h logger.h timing.h
timing.o      : timing.h
sor.o         : sor.h
parallel.o    : parallel.h

main.o        : helper.h init.h boundary_val.h uvp.h visual.h sor.h logger.h parallel.h
tests.o       : helper.h init.h boundary_val.h uvp.h visual.h sor.h timing.h logger.h parallel.h

# eof

