#include <cmath>
#include <iostream>
#include "Line.hpp"
#include "Msg.hpp"
#include "Vertex.hpp"
#include "Vector.hpp"


Line::Line(void) {

    direction = Vector();
    point = Vector();
}

bool Line::operator==(Line& rhs) {

    if ( this->point != rhs.point )
        return false;
    if ( this->direction != rhs.direction )
        return false;
    return true;
}

mpq_class Line::distanceToPoint(Vertex& p) {
    
    mpq_class x1 = point.x - p.x;
    mpq_class y1 = point.y - p.y;
    mpq_class z1 = point.z - p.z;
    Vector v(x1, y1, z1);
    Vector v2 = v.cross(this->direction);

    mpq_t x2, y2, z2, dx, dy, dz, ns1, ns2, ds1, ds2, r1;
    mpq_inits(x2, y2, z2, dx, dy, dz, ns1, ns2, ds1, ds2, r1, NULL);
    mpq_mul(x2, v2.x.get_mpq_t(), v2.x.get_mpq_t());
    mpq_mul(y2, v2.y.get_mpq_t(), v2.y.get_mpq_t());
    mpq_mul(z2, v2.z.get_mpq_t(), v2.z.get_mpq_t());
    mpq_mul(dx, direction.x.get_mpq_t(), direction.x.get_mpq_t());
    mpq_mul(dy, direction.y.get_mpq_t(), direction.y.get_mpq_t());
    mpq_mul(dz, direction.z.get_mpq_t(), direction.z.get_mpq_t());
    mpq_add(ns1, x2, y2);
    mpq_add(ns2, ns1, z2);
    mpq_add(ds1, dx, dy);
    mpq_add(ds2, ds1, dz);
    mpq_div(r1, ns2, ds2);
    
    mpq_class d(r1);
    
    return d;
}
