#ifndef Geometry_hpp
#define Geometry_hpp

#include <gmpxx.h>
#include <iostream>
#include <string>
#include "Plane.hpp"
#include "Vector.hpp"
#include "Vertex.hpp"



namespace Geometry {

    double      angle(Vector& a, Vector& b, Vector& c, Vector& n);
    Vector      centroid(std::vector<Vector*>& vecs);
    Vector      centroid(std::vector<Vertex*>& vecs);
    void        crossProduct(Vector& v1, Vector& v2, Vector& res);
    mpq_class   distanceSquared(Vector& v1, Vector& v2);
    mpq_class   dotProduct(Vector& v1, Vector& v2);
    bool        intersect(Plane& plane1, Plane& plane2, Plane& plane3, Vector& intersection);
}

#endif
