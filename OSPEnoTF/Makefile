SHELL = /bin/bash

CC = gcc
CFLAGS = 
ALL_CFLAGS = -Wall -Wno-unused -O2 $(CFLAGS)

OBJS =  HTIRestor.o octupole.o bsstep.o mmid.o nrutil.o odeint.o pzextr.o rkck.o rkqs.o 
##OBJS =  HTI.o octupole.o bsstep.o mmid.o nrutil.o odeint.o pzextr.o rkck.o rkqs.o 
LIBS = -lm
EXE = octupole

# pattern rule to compile object files from C files
# might not work with make programs other than GNU make
%.o : %.c Makefile nrutil.h HTI.h
	$(CC) $(ALL_CFLAGS) -c $< -o $@

all: $(EXE)

$(EXE): $(OBJS) Makefile
	$(CC) $(ALL_CFLAGS) $(OBJS) -o $(EXE) $(LIBS)

.PHONY : clean
clean:
	rm -f $(OBJS) $(EXE)
