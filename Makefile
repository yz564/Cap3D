CC := g++
CFLAGS := -Wall -Wextra -pedantic -g -fopenmp 
LFLAGS := -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
TARGET := cap3d
SRCS := main.cpp doConfig.cpp doMesh.cpp doMatrix.cpp setQuadrature.cpp
OBJS := $(SRCS:.cpp=.o)
DEPS := doConfig.h doMesh.h doMatrix.h

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(OBJS) -o $@ $(CFLAGS) $(LFLAGS)

%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

