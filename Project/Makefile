CC = gcc
CFLAGS = -Wall -pedantic -Werror -std=gnu99 -g   
.c.o:  ; $(CC) -c $(CFLAGS) $<

OBJ = 	helper.o\
      	init.o\
      	boundary_val.o\
		boundary.o\
      	uvp.o\
      	main.o\
      	visual.o\
		sor.o\
                FSI.o\

all:  $(OBJ)
	$(CC) $(CFLAGS) -o sim $(OBJ)  -lm

%.o : %.c
	$(CC) -c $(CFLAGS) $*.c -o $*.o

clean:
	rm $(OBJ)

sor.o	      : sor.h
helper.o      : helper.h 
init.o        : helper.h init.h 
boundary_val.o: helper.h boundary_val.h init.h
boundary.o    : boundary.h init.h
uvp.o         : helper.h uvp.h
visual.o      : helper.h
FSI.o         : FSI.h init.h
main.o        : helper.h init.h boundary_val.h uvp.h visual.h sor.h boundary.h FSI.h
