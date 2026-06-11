CXX = g++

CXXFLAGS = -O3 -march=native -std=c++23 -fopenmp -flto=4 
LDFLAGS  = -fopenmp -flto=4 -lboost_filesystem -lboost_iostreams

TARGET = main2

SRC = main2.cpp ElectricField.cpp distrib2d.cpp
OBJ = $(SRC:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(OBJ) $(LDFLAGS) -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)

rebuild: clean all

.PHONY: all clean rebuild