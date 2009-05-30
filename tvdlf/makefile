#
#
#
#CFLAGS = -D_INTRINSICS -D_INLINE_INTRINSICS -D_FASTMATH -float_const \
#	-fp_reorder -readonly_strings -O3
#CFLAGS = -O0 -Wall
#CFLAGS = -funroll-all-loops -O2 
CFLAGS = -O3
#CFLAGS = -pg
#CFLAGS = -g
#CC = ccc
CC = icc

SRCS = \
bounds.c diag.c dump.c image.c init.c main.c nrutil.c set_arrays.c \
step_ch.c tensor.c phys.c
 
OBJS = \
bounds.o diag.o dump.o image.o init.o main.o nrutil.o set_arrays.o \
step_ch.o tensor.o phys.o

twod: $(OBJS) makefile
	$(CC) $(CFLAGS) -o twod $(OBJS) -lm 

# dependencies
$(OBJS) : defs.h decs.h nrutil.h makefile

 
