# This is just to execute the tests automatically

.PHONY: all compile

all:
	mpirun -np 4 ./tests && cat Out/*.log

compile:
	$(MAKE) -C ../ tests

# eof

