CC = mpicc
CFLAGS= -g -Wall -O2
LIBS = -lgmp -lgsl -lgslcblas -lm -lcrypto

PHEFOLDER = libpaillier-0.8
BUILDFOLDER = build

MAIN_SRC = sim.c
SRC = encoding.c encoded_paillier_agg.c key_distribution.c enc_matrix.c sensor.c navigator.c
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


