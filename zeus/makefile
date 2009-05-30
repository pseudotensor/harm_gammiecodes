#
#
CFLAGS = -fast -arch ev67
LDFLAGS = -lm
#CC = /usr/local/bin/gcc
CC = ccc

SRCS = \
main.c aupdate.c bounds.c diag.c dump.c fft_mayer.c image.c init.c init_spec.c \
moc_ct.c nrutil.c potsolv.c ranc.c set_arrays.c \
step_bz.c step_drive.c step_grid.c step_ie.c step_pb.c step_pg.c \
step_res.c step_trans.c step_visc.c stepvar.c sweepx.c \
sweepy.c tfft.c timestep.c \
spec_diag.c zero_arrays.c

OBJS = \
main.o aupdate.o bounds.o diag.o dump.o fft_mayer.o image.o init.o init_spec.o \
moc_ct.o nrutil.o potsolv.o ranc.o set_arrays.o \
step_bz.o step_drive.o step_grid.o step_ie.o step_pb.o step_pg.o \
step_res.o step_trans.o step_visc.o stepvar.o sweepx.o \
sweepy.o tfft.o timestep.o \
spec_diag.o zero_arrays.o

twod: $(OBJS) $(SRCS) makefile
	$(CC) $(CFLAGS) -o twod $(OBJS) $(LDFLAGS)

# dependencies
$(OBJS) : defs.h decs.h
 
