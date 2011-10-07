CC=g++
LDFLAGS= -lm
CFLAGS= -Wall -pedantic 

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

