CC = mpicc
CFLAGS= -g -Wall -O3
LIBS = -lgmp -lgsl -lgslcblas -lssl -lcrypto

PHEFOLDER = libpaillier-0.8
BUILDFOLDER = build

MAIN_SRC = sim.c
SRC = key_setup.c time_hash.c sensor.c navigator.c encoding.c matrix_ops.c
OBJ = $(patsubst %.c,$(BUILDFOLDER)/%.o, $(SRC))

SIMNAME = sim
SIM = $(BUILDFOLDER)/$(SIMNAME)
PHEDEP = $(PHEFOLDER)/paillier.o

.PHONY: clean all

all: $(SIM)

$(SIM): $(MAIN_SRC) $(OBJ) $(PHEDEP)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

$(BUILDFOLDER)/%.o: %.c %.h
	$(CC) -c -o $@ $< $(CFLAGS)

$(BUILDFOLDER):
	mkdir -p $(BUILDFOLDER)

$(PHEDEP):
	$(MAKE) -C $(PHEFOLDER)

clean:
	rm -f $(BUILDFOLDER)/*


