#
#
# DECs
#CFLAGS = -O2
CFLAGS = -ipo -O3
LDFLAGS = -lm
#
CC = icc
#CC = /usr/local/bin/gcc
#CC = cc

SRCS = \
main.c array_init.c bounds.c diag.c dump.c image.c init.c moc_ct.c \
ranc.c step_drive.c step_ie.c step_pg.c step_trans.c step_visc.c stepvar.c \
sweepx.c sweepy.c sweepym.c sweepz.c timestep.c remap.c dfour1.c drealft.c

OBJS = \
main.o array_init.o bounds.o diag.o dump.o image.o init.o moc_ct.o \
ranc.o step_drive.o step_ie.o step_pg.o step_trans.o step_visc.o stepvar.o \
sweepx.o sweepy.o sweepym.o sweepz.o timestep.o remap.o dfour1.o drealft.o

mthreed: $(OBJS) $(SRCS) makefile
	$(CC) $(CFLAGS) -o mthreed $(OBJS) $(LDFLAGS)

# dependencies
$(OBJS) : defs.h decs.h
 
