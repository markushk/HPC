CXX = mpic++
CXXFLAGS = -std=c++11 -Wall

TARGET = Ex1

$(TARGET): Ex1.cpp
	$(CXX) $(CXXFLAGS) -o $(TARGET) Ex1.cpp

clean:
	rm -f $(TARGET)