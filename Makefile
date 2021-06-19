
CC     = gcc 
CFLAGS = -Wall
LFLAGS = -lm

PROG = edoedp
OBJS = edoedp.o

.PHONY: limpa faxina clean purge all

%.o: %.c 
	$(CC) -c $(CFLAGS) $<

$(PROG) : % :  $(OBJS) 
	$(CC) -o $@ $^ $(LFLAGS)

limpa clean:
	@rm -f *~ *.bak

faxina purge:   limpa
	@rm -f *.o core a.out *.txt
	@rm -f $(PROG)
