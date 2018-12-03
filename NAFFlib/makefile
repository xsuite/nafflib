objdir = object
srcdir = source
incdir = include

python3inc = `python3 -c "from distutils import sysconfig; print(sysconfig.get_python_inc())"`
numpy3inc = `python3 -c "import numpy; print(numpy.get_include())"`
python2inc = `python2 -c "from distutils import sysconfig; print(sysconfig.get_python_inc())"`
numpy2inc = `python2 -c "import numpy; print(numpy.get_include())"`

objects = $(objdir)/windows.o $(objdir)/signal_processing.o $(objdir)/brent.o $(objdir)/frequency.o $(objdir)/fft.o  

cc = gcc
cflags = -O3 -std=c99 -Wall -fPIC -I$(incdir) 
ldflags =  -lm #-lpython

all: NAFFlib_c.so NAFFlib2_c.so
py3: NAFFlib_c.so
py2: NAFFlib2_c.so

NAFFlib_c.so: $(objects) $(objdir)/pynafflib3.o
	@$(cc) --shared -fPIC $(cflags) $^ -o $@ $(ldflags) 

NAFFlib2_c.so: $(objects) $(objdir)/pynafflib2.o
	@$(cc) --shared -fPIC $(cflags) $^ -o $@ $(ldflags) 

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

$(objdir)/pynafflib3.o: $(srcdir)/pynafflib.c $(incdir)/pynafflib.h 
	@echo python3 include path: $(python3inc)
	@echo numpy \(for python3\) include path: $(numpy3inc)
	@$(cc) -c $(cflags) -I$(python3inc) -I$(numpy3inc) $< -o $@ 

$(objdir)/pynafflib2.o: $(srcdir)/pynafflib.c $(incdir)/pynafflib.h 
	@echo python2 include path: $(python2inc)
	@echo numpy \(for python2\) include path: $(numpy2inc)
	@$(cc) -c $(cflags) -I$(python2inc) -I$(numpy2inc) $< -o $@ 

clean:
	rm $(objdir)/*
	rm NAFFlib_c.so
	rm NAFFlib2_c.so
