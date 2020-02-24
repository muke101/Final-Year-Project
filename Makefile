CC=gcc
STANDARD=-O2 -pipe -march=native -fomit-frame-pointer -fdevirtualize-at-ltrans -lm
LTO=-flto=4
GRAPHITE=-fgraphite-identity -floop-nest-optimize -floop-interchange -ftree-loop-distribution -floop-strip-mine -floop-block
IPA=-fipa-pta
PGO=-profile-use=profile/
CFLAGS=$(STANDARD) $(LTO) $(GRAPHITE) $(IPA)
DEPS=montecarlos.h

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) 

make: nfcaMC.o multiGluonMC.o addativeMC.o
	$(CC) -o montecarlo nfcaMC.o multiGluonMC.o addativeMC.o $(CFLAGS) 
	
