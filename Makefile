CC := g++
CFLAGS := -Wall -Wextra -pedantic -g
TARGET := cap3d
SRCS := main.cpp doConfig.cpp doMesh.cpp
OBJS := $(SRCS:.cpp=.o)
DEPS := doConfig.cpp doMesh.h 

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@

%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

