FLAGS = -std=c++11#-Wall -Wextra -Werror 

TARGET = l

HDRS = \
		includes

SRCS = \
	   src/*.cpp


.PHONY: all clean

all: $(SRCS)
	$(CXX) $(FLAGS) -I $(HDRS) -o $(TARGET) $(CXXFLAGS) $(SRCS) 
	./$(TARGET)
clean:
	rm -rf $(TARGET)