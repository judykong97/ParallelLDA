CC=g++

MPI=-DMPI
MPICC = mpic++

DEBUG=0
CFLAGS= -g -O3 -Wall -DDEBUG=$(DEBUG)
LDFLAGS= -lm
# DDIR = ./data

CFILES = main.cpp lda.cpp lda_utils.cpp # crun.c graph.c simutil.c sim.c rutil.c cycletimer.c
HFILES = lda.h lda_utils.h # HFILES = crun.h rutil.h cycletimer.h

# GFILES = gengraph.py grun.py rutil.py sim.py viz.py  regress.py benchmark.py grade.py


# DFILES = $(DDIR)/g-t3600.gph $(DDIR)/g-t32400.gph $(DDIR)/g-t4.gph $(DDIR)/g-t400.gph \
	 $(DDIR)/g-u3600.gph $(DDIR)/g-u32400.gph $(DDIR)/g-u4.gph $(DDIR)/g-u400.gph 

# RFILES  = $(DDIR)/r-3600-d10.rats  $(DDIR)/r-3600-u10.rats \
	$(DDIR)/r-32400-d32.rats  $(DDIR)/r-32400-u32.rats \
	$(DDIR)/r-4-d1.rats  $(DDIR)/r-4-u1.rats \
	$(DDIR)/r-400-d10.rats $(DDIR)/r-400-u10.rats 

all: lda lda-mpi
	
lda: $(CFILES) $(HFILES)
	$(CC) -std=c++11 -lstdc++ -o lda $(CFILES) $(LDFLAGS)

lda-mpi: $(CFILES) $(HFILES) $(XCFILES) $(XHFILES)
	$(MPICC) -std=c++11 -lstdc++ $(MPI) -DOMPI_SKIP_MPICXX -o lda-mpi $(CFILES) $(XCFILES) $(LDFLAGS)

# all: crun crun-mpi

# crun: crun-seq
# 	cp -p crun-seq crun

# crun-seq: $(CFILES) $(HFILES) 
# 	$(CC) $(CFLAGS) -o crun-seq $(CFILES) $(LDFLAGS)

# crun-mpi: $(CFILES) $(XCFILES) $(HFILES) $(XHFILES)
# 	$(MPICC) $(CFLAGS) $(MPI) -o crun-mpi $(CFILES) $(XCFILES) $(LDFLAGS)

clean:
	rm -f lda 
	rm -f lda-mpi

