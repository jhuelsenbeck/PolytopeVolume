#include "Edge.hpp"
#include "Geometry.hpp"
#include "VertexFactory.hpp"



Edge::Edge(Vertex* _v1, Vertex* _v2) : Line() {

    v1 = _v1;
    v2 = _v2;
    
    setPoint(*v1);
    Vector temp;
    temp.x = v2->x - v1->x;
    temp.y = v2->y - v1->y;
    temp.z = v2->z - v1->z;
    setDirection( temp );
    dd = Geometry::distanceSquared(*v1, *v2);
    //dd = v1->distanceSquared(*v2);
}

