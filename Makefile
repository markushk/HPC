CXX = mpic++
CXXFLAGS = -std=c++11 -Wall -O3
INCLUDE_DIR = include
SRC_DIR = src
LIB_DIR = lib

TARGET = Ex1
MAIN_TARGET = Ex2

$(TARGET): $(SRC_DIR)/Ex1.cpp
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC_DIR)/Ex1.cpp

$(MAIN_TARGET): $(SRC_DIR)/main.cpp $(LIB_DIR)/GameOfLife.cpp
	$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) -o $(MAIN_TARGET) $(SRC_DIR)/main.cpp $(LIB_DIR)/GameOfLife.cpp


clean:
	rm -f $(TARGET) $(MAIN_TARGET)
