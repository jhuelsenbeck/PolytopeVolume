#ifndef Edge_H
#define Edge_H

#include "Line.hpp"
#include "Vector.hpp"
#include "Vertex.hpp"
class Vector;


class Edge : public Line {

    public:
                        Edge(void) = delete;
                        Edge(Vertex* v1, Vertex* v2);
        double          getDistance(void) { return d; }
        Vertex*         isIntersected(Line& line);

    private:
        class Vector*   v1;
        class Vector*   v2;
        double          d;
};

#endif
