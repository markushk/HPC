CXX = mpic++
CXXFLAGS = -std=c++11 -Wall

TARGET = test1
MAIN_TARGET = test2

$(TARGET): test.cpp
	$(CXX) $(CXXFLAGS) -o $(TARGET) test.cpp

$(MAIN_TARGET): test_copy.cpp
	$(CXX) $(CXXFLAGS) -o $(MAIN_TARGET) test_copy.cpp

clean:
	rm -f $(TARGET) $(MAIN_TARGET)
