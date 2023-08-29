#include "Edge.hpp"
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
    d = v1->distance(*v2);
}

Vertex* Edge::isIntersected(Line& line) {

    // check to see if it intersects the line from the edge
    Vector v = this->direction.cross(line.getDirection());
    //std::cout << "   Is intersected internals: " << v << std::endl;
    if(v.x == 0 && v.y == 0 && v.z == 0)
        return NULL;
    
    // check that the intersection point is inbetween the two points
    Vector i = this->intersect(line);
    double d1 = i.distance(*v1);
    double d2 = i.distance(*v2);
    //std::cout << "   " << d1 << " + " << d2 << " = " << d ;
    if ( fabs(d1 + d2 - d) > 10e-8)
        {
        //std::cout << "  (fail)" << std::endl;
        return NULL;
        }
    //std::cout << "  (success)" << std::endl;
    VertexFactory& vFactory = VertexFactory::vertexFactoryInstance();
    Vertex* newVertex = vFactory.getVertex(i);
    return newVertex;
}
