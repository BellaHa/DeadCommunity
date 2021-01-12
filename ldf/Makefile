PLATFORM = linux64
INC_GRAPH= graph
CPPFLAG  = -m64 -g -O3 -w 
INC_DIR  = -I$(INC_GRAPH)
LIBS	 = 
OBJS	= ldf.o clustering.o

ldf: $(OBJS)
	 g++ $(CPPFLAG) -o ldf $(OBJS) $(INC_DIR) $(LIBS) -lpthread -lm

ldf.o: ldf.cpp
	g++ $(CPPFLAG) -c ldf.cpp $(INC_DIR)

clustering.o: clustering.cpp clustering.h
	g++ $(CPPFLAG) -c clustering.cpp $(INC_DIR)
clean:
	rm -f *.o *_c *_c++ *.class *.log *.rlp *.lp *.bas *.ilp
