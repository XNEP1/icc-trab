    CC     = gcc -g
    CFLAGS = -DLIKWID_PERFMON -I${LIKWID_INCLUDE} -O3 -mavx -march=native
    LFLAGS = -lm -llikwid -L${LIKWID_LIB} -O3 -mavx -march=native

      PROG = cgSolver
      OBJS = sisDiag.o \
             Metodos.o \
             utils.o \
             $(PROG).o

.PHONY: limpa faxina clean distclean purge all

%.o: %.c %.h utils.h sislin.h
	$(CC) -c $(CFLAGS) $<

$(PROG):  $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS)

clean:
	@rm -f *~ *.bak

purge:   clean
	@rm -f *.o core a.out
	@rm -f $(PROG)

