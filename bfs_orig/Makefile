# Makefile for breadth-first search
# Sharlee Climer
# October, 2007

# Aditya Karnam
# October, 2017
# Added cleaning *.out objects
# Added cleaning *.log files

CC	= g++
CFLAGS 	= -g -std=c++11
TARGET	= erbga
OBJS	= erbga_driver.o network.o erbga.o

$(TARGET):	$(OBJS)
		$(CC) -o $(TARGET) $(OBJS)

bfsNet.o:	erbga_driver.cpp erbga_driver.h timer.h
		$(CC) $(CFLAGS) -c erbga_driver.cpp

network.o:	network.cpp network.h
		$(CC) $(CFLAGS) -c network.cpp

cd_spnet.o:	erbga.cpp erbga.h network.h
		$(CC) $(CFLAGS) -c erbga.cpp

clean:
		/bin/rm -f *.o $(TARGET) *.out

cleanlog:
		/bin/rm -f *.log
