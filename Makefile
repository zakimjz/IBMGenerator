HFILES1 = glob.h dist.h
HFILES2 = gen.h

CFILES1 =  ran0.cc ran1.cc expdev.cc gammln.cc gasdev.cc poidev.cc dist.cc
CFILES2 =  gen.cc main.cc command.cc

OBJECTS1 = ran0.o ran1.o expdev.o gammln.o gasdev.o poidev.o dist.o
OBJECTS2 = gen.o main.o command.o

LIBES = -lm
#CC = cxx
CC = g++
CFLAGS = -O2 -Wno-write-strings
.SUFFIXES: .o .cc

gen:	$(OBJECTS1) $(OBJECTS2)
	$(CC) $(CFLAGS) $(OBJECTS1) $(OBJECTS2) $(LIBES) -o gen

test:	$(OBJECTS1) test.o
	$(CC) $(CFLAGS) $(OBJECTS1) test.o $(LIBES) -o test

clean:
	rm *~ *.o gen

.cc.o:
	$(CC) $(CFLAGS) -c $<

$(OBJECTS1): $(HFILES1)
$(OBJECTS2): $(HFILES1) $(HFILES2)
ran0.o: ran0.cc
