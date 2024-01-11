CXX = mpic++
CXXFLAGS = -std=c++11 -Wall

TARGET = Ex1
MAIN_TARGET = test

$(TARGET): Ex1.cpp
	$(CXX) $(CXXFLAGS) -o $(TARGET) Ex1.cpp

$(MAIN_TARGET): test.cpp
	$(CXX) $(CXXFLAGS) -o $(MAIN_TARGET) test.cpp

clean:
	rm -f $(TARGET) $(MAIN_TARGET)
