LIBS     =  -lm -lgmp -lcddgmp  -DGMPRATIONAL

CFLAGS   =  -O3 -DGMPRATIONAL

CC       =  g++

OBJECTS  =  main.o Geometry.o McmcState.o Msg.o Plane.o Polyhedron.o Probability.o RandomVariable.o RateMatrix.o Vector.o Vertex.o VertexFactory.o

PROGS    = rjtest

all:		$(PROGS)

rjtest:		$(OBJECTS)
		$(CC) $(CFLAGS) $(OBJECTS) $(LIBS) -o rjtest
		
main.o:	main.cpp
		$(CC) $(CFLAGS) -c main.cpp

Geometry.o:	Geometry.cpp
		$(CC) $(CFLAGS) -c Geometry.cpp

McmcState.o:	McmcState.cpp
		$(CC) $(CFLAGS) -c McmcState.cpp

Msg.o:	Msg.cpp
		$(CC) $(CFLAGS) -c Msg.cpp

Plane.o:	Plane.cpp
		$(CC) $(CFLAGS) -c Plane.cpp

Polyhedron.o:	Polyhedron.cpp
		$(CC) $(CFLAGS) -c Polyhedron.cpp

Probability.o:	Probability.cpp
		$(CC) $(CFLAGS) -c Probability.cpp

RandomVariable.o:	RandomVariable.cpp
		$(CC) $(CFLAGS) -c RandomVariable.cpp

RateMatrix.o:	RateMatrix.cpp
		$(CC) $(CFLAGS) -c RateMatrix.cpp

Vector.o:	Vector.cpp
		$(CC) $(CFLAGS) -c Vector.cpp

Vertex.o:	Vertex.cpp
		$(CC) $(CFLAGS) -c Vertex.cpp

VertexFactory.o:	VertexFactory.cpp
		$(CC) $(CFLAGS) -c VertexFactory.cpp

clean:		
		rm -f *.o
