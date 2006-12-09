#
# PYTHON_INC points to the Pythn include tree 
# 
#
PYTHON_INC = -I/usr/include/python2.3/ -I/usr/include/python2.3/Numeric/
PY_LIB     = $(DEST)/lib/python

#
# modify this if you want to use other compilers or flags
#
CPP    = g++
CC     = gcc
INC    = -I./ -I../api $(PYTHON_INC) -DTIMESTAMP=$(DATE)


OPT    = -c -O3 
CFLAGS = $(OPT) 

COBJS = \
	Algorithms.o

DATE    = "\"`date`\""

#
# don't change anything below
#

all :	$(COBJS) $(OBJS)
	$(CPP)  $(COBJS) -shared -lstdc++ -w -o Algorithms.so

clean:	
	rm -f $(OBJS) $(COBJS)
	rm -f *.so

#dependecies
vpath %.hh  ./api/
vpath %.cpp ./api/

SRC = ./

$(COBJS): %.o: %.cpp
	$(CPP) $(CFLAGS) $(INC) $(SRC)$<
