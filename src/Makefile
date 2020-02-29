CC=gcc
STANDARD=-O3 -pipe -march=native -fomit-frame-pointer -fdevirtualize-at-ltrans 
WARNINGS=-Wall -Wextra
LTO=-flto=4
GRAPHITE=-fgraphite-identity -floop-nest-optimize -floop-interchange -ftree-loop-distribution -floop-strip-mine -floop-block
IPA=-fipa-pta
SECURITY=-no-pie -fno-stack-protector -fno-stack-clash-protection
PGO=-fprofile-use=profile/
PGO-GEN=-fprofile-generate=profile/

IDIR=../include
_DEPS = emmisionfunctions.h
DEPS= $(patsubst %,$(IDIR)/%,$(_DEPS))

LIBS =-lm

OBJ = addativeMC.o nfcaMC.o multiGluonMC.o 

CFLAGS =$(STANDARD) $(LTO) $(GRAPHITE) $(IPA) $(SECURITY) -I$(IDIR) 

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) 

../montecarlo: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f *.o ../montecarlo 

