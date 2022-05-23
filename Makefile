CC = mpicc
CFLAGS = -std=c99 -g -O3
LIBS = -lm

BIN = malaria_sim

all: $(BIN)

malaria_sim: malaria_sim.c
	$(CC) -fsanitize=address $(CFLAGS) -o $@ $< $(LIBS)
	
clean:
	$(RM) $(BIN)
