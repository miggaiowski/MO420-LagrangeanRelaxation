CC=/opt/local/bin/g++-mp-4.4
LDFLAGS= -lm
CFLAGS= -Wall -pedantic -O3 -fopenmp

.PHONY: all clean cleanall

all: relaxlag

relaxlag: main.cpp
	$(CC) $(LDFLAGS) $(CFLAGS) main.cpp -o relaxlag

# main.o: main.cpp #graph.h graph.cpp
# 	$(CC) $(CFLAGS) -c main.cpp

# graph.o: graph.cpp graph.h
# 	$(CC) $(CFLAGS) -c graph.cpp

# asserting.o: asserting.cpp asserting.h
# 	$(CC) $(CFLAGS) -c asserting.cpp

clean:
	rm -v relaxlag *.o 

cleanall:
	rm -v relaxlag *.o *~

