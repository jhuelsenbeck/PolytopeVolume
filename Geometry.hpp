#ifndef Geometry_hpp
#define Geometry_hpp

#include <gmpxx.h>
#include <iostream>
#include <string>
#include "Line.hpp"
#include "Plane.hpp"
#include "Vector.hpp"



namespace Geometry {

    bool intersect(Plane& plane1, Plane& plane2, Plane& plane3, Vector& intersection);
    bool intersect(Plane& plane1, Plane& plane2, Line& intersection);
}

#endif
