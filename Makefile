CC = clang
CFLAGS= -g -Wall -O3
LIBS = -lgmp -lgsl -lgslcblas -lssl -lcrypto

OREFOLDER = ore
PHEFOLDER = libpaillier-0.8
BUILDFOLDER = build

MAIN_SRC = sim.c
SRC = encoding.c matrix_ops.c cmp_lists.c
OBJ = $(patsubst %.c,$(BUILDFOLDER)/%.o, $(SRC))

SIMNAME = sim
SIM = $(BUILDFOLDER)/$(SIMNAME)
OREDEP = $(OREFOLDER)/build/ore_blk.o $(OREFOLDER)/build/crypto.o
PHEDEP = $(PHEFOLDER)/paillier.o

.PHONY: clean all

all: $(SIM)

$(SIM): $(MAIN_SRC) $(OBJ) $(OREDEP) $(PHEDEP)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

$(BUILDFOLDER)/%.o: %.c %.h
	$(CC) -c -o $@ $< $(CFLAGS)

$(BUILDFOLDER):
	mkdir -p $(BUILDFOLDER)

$(OREDEP):
	$(MAKE) -C $(OREFOLDER)

$(PHEDEP):
	$(MAKE) -C $(PHEFOLDER)

clean:
	rm -f $(BUILDFOLDER)/*


