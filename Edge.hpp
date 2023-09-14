#ifndef Edge_H
#define Edge_H

#include "Vector.hpp"
#include "Vertex.hpp"


class Edge : public Line {

    public:
                        Edge(void) = delete;
                        Edge(Vertex* v1, Vertex* v2);
        mpq_class&      getDistanceSquared(void) { return dd; }

    private:
        class Vector*   v1;
        class Vector*   v2;
        mpq_class       dd;
};

#endif
