
INCLUDE = ../../include
SRC = ../../src

CC = g++
CFLAGS = -Ofast -ffast-math -fno-trapping-math -ftree-vectorize -std=c++11 -Wall -g -Werror -I$(INCLUDE)
LDFLAGS = -lm

OBJS = search.o $(SRC)/point_functions.o $(SRC)/hash_table.o $(SRC)/cube_table.o $(SRC)/hash_functions.o $(SRC)/algorithms.o $(SRC)/search_methods.o $(SRC)/grid_table.o $(SRC)/frechet.o $(SRC)/curve.o $(SRC)/config.o $(SRC)/point.o $(SRC)/interval.o $(SRC)/simplification.o
EXEC = search.out

#----------------------------------------------------------------------------#
    # Arguments for make run / make valgrind
#---------------------------------------#
    # Ai) LSH
# ARGS = -i ../../nasd_input.csv -q ../../nasd_query.csv -k 4 -L 4 -M 100 -probes 10 -o ../../outputFile -algorithm LSH
#---------------------------------------#
    # Ai) Hypercube
# ARGS = -i ../../nasd_input.csv -q ../../nasd_query.csv -k 4 -L 4 -M 100 -probes 10 -o ../../outputFile -algorithm Hypercube
#---------------------------------------#
    # Aii) DISCRETE FRECHET
ARGS = -i ../../nasd_input.csv -q ../../nasd_query.csv -k 4 -L 4 -M 100 -probes 10 -o ../../outputFile -algorithm Frechet -metric discrete -delta 1.0
#---------------------------------------#
    # Aiii) CONTINUOUS FRECHET
# ARGS = -i ../../nasd_input.csv -q ../../nasd_query.csv -k 4 -L 1 -M 100 -probes 10 -o ../../outputFile -algorithm Frechet -metric continuous -delta 1.0 -nobrute
#----------------------------------------------------------------------------#

$(EXEC): $(OBJS) $(LIBS)
	$(CC) $(CFLAGS) $(OBJS) -o $(EXEC) $(LDFLAGS)

clean:
	rm -f $(OBJS) $(EXEC)

run: $(EXEC)
	./$(EXEC) $(ARGS)

valgrind: $(EXEC)
	valgrind --tool=memcheck --leak-check=yes --leak-check=full --show-leak-kinds=all ./$(EXEC) $(ARGS)
