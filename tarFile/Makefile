# Makefile for breadth-first search
# Sharlee Climer
# October, 2007


CC	= g++
CFLAGS 	= -g
TARGET	= bfs
OBJS	= bfsNet.o network.o

$(TARGET):	$(OBJS)
		$(CC) -o $(TARGET) $(OBJS)

bfsNet.o:	bfsNet.cpp bfsNet.h timer.h
		$(CC) $(CFLAGS) -c bfsNet.cpp

network.o:	network.cpp network.h
		$(CC) $(CFLAGS) -c network.cpp

clean:
		/bin/rm -f *.o $(TARGET)
