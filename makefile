objdir = object
srcdir = source
incdir = include

#headers = $(wildcard $(incdir)/*.h)
#objects = $(subst $(incdir), $(objdir), $(headers:.h=.o))
objects = $(objdir)/windows.o $(objdir)/signal_processing.o $(objdir)/brent.o $(objdir)/frequency.o $(objdir)/fft.o 


cc = gcc
#cflags = --shared -fPIC -O2 -std=c99 -Wall -I$(INCDIR)
#cflags = -I$(incdir)
cflags = -O2 -std=c99 -Wall -fPIC -I$(incdir)
#cflags = -O2 -std=c99 -Wall -fPIC -I$(incdir)
ldflags = -lm

all: nafflib.so

nafflib.so: $(objects)
	@#$(cc) $(cflags) $^ -o $@  
	@$(cc) --shared -fPIC $(cflags) $^ -o $@  

#$(objdir)/toy1.o: $(srcdir)/toy1.c $(incdir)/toy1.h 
#	@$(cc) -c $(cflags) $< -o $@ 

$(objdir)/windows.o: $(srcdir)/windows.c $(incdir)/windows.h 
	@$(cc) -c $(cflags) $< -o $@ 

$(objdir)/signal_processing.o: $(srcdir)/signal_processing.c $(incdir)/signal_processing.h 
	@$(cc) -c $(cflags) $< -o $@ 

$(objdir)/brent.o: $(srcdir)/brent.c $(incdir)/brent.h 
	@$(cc) -c $(cflags) $< -o $@ 

$(objdir)/frequency.o: $(srcdir)/frequency.c $(incdir)/frequency.h 
	@$(cc) -c $(cflags) $< -o $@ 

$(objdir)/fft.o: $(srcdir)/fft.c $(incdir)/fft.h 
	@$(cc) -c $(cflags) $< -o $@ 

clean:
	rm $(objdir)/*
