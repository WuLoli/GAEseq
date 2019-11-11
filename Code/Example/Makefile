
CFLAGS  =	-O3 -DHAVE_TRUNC -fPIC -std=gnu++0x

EXOBJS	=	ExtractMatrix.o


all:	ExtractMatrix

clean:
	rm $(EXOBJS) ExtractMatrix

ExtractMatrix:	$(EXOBJS)
	g++ $(CFLAGS)  -o $@  $(EXOBJS) 
ExtractMatrix.o:	ExtractMatrix.cpp 
	g++   $(CFLAGS)  -c -o $@  $<
