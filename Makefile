 
    CC     = gcc -g
    CFLAGS = -DLIKWID_PERFMON -I${LIKWID_INCLUDE} -O3 -march=native -mavx
    LFLAGS = -lm -llikwid -L${LIKWID_LIB} -O3 -march=native -mavx

      PROG = cgSolver
      OBJS = sisDiag.o \
             Metodos.o \
             utils.o \
             $(PROG).o

.PHONY: limpa faxina clean distclean purge all

%.o: %.c %.h utils.h sislin.h
	$(CC) -c $(CFLAGS) $<

all: $(PROG)

unroll:         CFLAGS += -DUNROLL
unroll:         $(PROG)

$(PROG):  $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS)

clean:
	@rm -f *~ *.bak

purge:   clean
	@rm -f *.o core a.out
	@rm -f $(PROG)

