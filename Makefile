CC=gcc
STANDARD=-O2 -pipe -march=native -fomit-frame-pointer -fdevirtualize-at-ltrans -lm
LTO=-flto=4
GRAPHITE=-fgraphite-identity -floop-nest-optimize -floop-interchange -ftree-loop-distribution -floop-strip-mine -floop-block
IPA=-fipa-pta
PGO=-profile-use=profile/
PGO-GEN=-profile-generate=profile/
CFLAGS=$(STANDARD) $(LTO) $(GRAPHITE) $(IPA)

make: addativeMC.o
	$(CC) -o montecarlo addativeMC.o $(CFLAGS) 

