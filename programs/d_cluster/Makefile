
INCLUDE = ../../include
SRC = ../../src

CC = g++
CFLAGS = -Ofast -ffast-math -fno-trapping-math -ftree-vectorize -std=c++11 -Wall -g -Werror -I$(INCLUDE)
LDFLAGS = -lm

OBJS = cluster.o $(SRC)/point_functions.o $(SRC)/cluster_methods.o $(SRC)/hash_table.o $(SRC)/cube_table.o $(SRC)/hash_functions.o $(SRC)/algorithms.o $(SRC)/grid_table.o $(SRC)/frechet.o $(SRC)/curve.o $(SRC)/config.o $(SRC)/point.o $(SRC)/interval.o $(SRC)/simplification.o
EXEC = cluster.out


# -update <Mean_Frechet or Mean_Vector>
# -assignment <Classic or LSH or Hypercube or LSH_Frechet>

ARGS = -i ../../nasd_input.csv -c ../../cluster.conf -o ../../outputFile -update Mean_Frechet -assignment LSH_Frechet -complete

$(EXEC): $(OBJS) $(LIBS)
	$(CC) $(CFLAGS) $(OBJS) -o $(EXEC) $(LDFLAGS)

clean:
	rm -f $(OBJS) $(EXEC)

run: $(EXEC)
	./$(EXEC) $(ARGS)

valgrind: $(EXEC)
	valgrind --tool=memcheck --leak-check=yes --leak-check=full --show-leak-kinds=all ./$(EXEC) $(ARGS)
