CC	= g++ -std=c++11
CPPFLAGS= -g -Wno-deprecated -O3 -c -I $(INC1) -I $(INC2)
INC1=./local_fibheap/Fibonacci
INC2=./regHeap
LDFLAGS	= -O3 
SOURCES	= main.cpp graph.cpp bitfield.cpp graphDualCriteria.cpp regHeap/regHeap.cpp local_fibheap/Fibonacci/fibheap.cpp
OBJECTS	= $(SOURCES:.cpp=.o)
EXECUTABLE=XuantemporalGraph

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE) : $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o : 
	$(CC) $(CPPFLAGS) $< -o $@

clean:
	rm -f *.o $(INC1)/*.o $(INC2)/*.o XuantemporalGraph
