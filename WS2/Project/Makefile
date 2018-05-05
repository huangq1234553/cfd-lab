CC = gcc
CFLAGS = -Wall -pedantic -Werror
.c.o:  ; $(CC) -c $(CFLAGS) $<

OBJ = 	helper.o\
      	init.o\
      	boundary_val.o\
      	uvp.o\
      	sor.o\
      	main.o\
      	visual.o\
      	logger.o


all:  $(OBJ)
	$(CC) $(CFLAGS) -o sim $(OBJ)  -lm

%.o : %.c
	$(CC) -c $(CFLAGS) $*.c -o $*.o

clean:
	rm $(OBJ)

image:
    convert -posterize 2 -depth 1 -compress none -negate testGeometry.jpg testGeometry.pgm; sed -i s/"255"/"1"/g testGeometry.pgm

helper.o      : helper.h logger.h
init.o        : helper.h init.h 
boundary_val.o: helper.h boundary_val.h 
uvp.o         : helper.h uvp.h
visual.o      : helper.h

main.o        : helper.h init.h boundary_val.h uvp.h visual.h sor.h logger.h
