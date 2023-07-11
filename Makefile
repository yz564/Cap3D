CC := g++
CFLAGS := -Wall -Wextra -pedantic -g -fopenmp
TARGET := cap3d
SRCS := main.cpp doConfig.cpp doMesh.cpp doMatrix.cpp
OBJS := $(SRCS:.cpp=.o)
DEPS := doConfig.h doMesh.h doMatrix.h

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@

%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

