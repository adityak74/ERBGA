.SUFFIXES: .c .cc .h

SOURCES = gml_scanner.c gml_parser.c gml_xgmml.c
OBJECTS = $(SOURCES:.c=.o)
CC = gcc
CFLAGS = -g

all : gml_xgmml

gml_xgmml: $(OBJECTS)
	$(CC) -o gml_xgmml.exe $(OBJECTS)

.PHONY: clean

clean: 
	-rm -f $(OBJECTS) gml_xgmml.exe



