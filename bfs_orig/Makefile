# Makefile for breadth-first search
# Sharlee Climer
# October, 2007

# Aditya Karnam
# October, 2017
# Added cleaning *.out objects
# Added cleaning *.log files

CC	= g++
CFLAGS 	= -g -std=c++11
TARGET	= bfs
OBJS	= bfsNet.o network.o cd_spnet.o

$(TARGET):	$(OBJS)
		$(CC) -o $(TARGET) $(OBJS)

bfsNet.o:	bfsNet.cpp bfsNet.h timer.h
		$(CC) $(CFLAGS) -c bfsNet.cpp

network.o:	network.cpp network.h
		$(CC) $(CFLAGS) -c network.cpp

cd_spnet.o:	cd_spnet.cpp cd_spnet.h network.h
		$(CC) $(CFLAGS) -c cd_spnet.cpp

clean:
		/bin/rm -f *.o $(TARGET) *.out

cleanlog:
		/bin/rm -f *.log
