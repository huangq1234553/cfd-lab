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
TEST_OBJ = Tests/tests.o Tests/testing.o Tests/uvComm.test.o
ALL_OBJ = $(OBJ) $(SIM_OBJ) $(TEST_OBJ)

all:  $(OBJ) $(SIM_OBJ)
	$(CC) $(CFLAGS) -o sim $(OBJ) $(SIM_OBJ) -lm

tests: $(OBJ) $(TEST_OBJ)
	$(CC) $(CFLAGS) -o Tests/tests $(OBJ) $(TEST_OBJ) -lm

%.o : %.c
	$(CC) -c $(CFLAGS) $*.c -o $*.o

clean:
	rm $(ALL_OBJ)

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

