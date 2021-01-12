CXXFLAGS =	-O3 -g -pg -Wall -fmessage-length=0 -I/home/tdinh/igraph/include 

OBJS =		main.o

#LIBS =		-L/home/tdinh/igraph/lib -ligraph

LIBS = 

TARGET =	main

$(TARGET):	$(OBJS)
	$(CXX) -pg -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
