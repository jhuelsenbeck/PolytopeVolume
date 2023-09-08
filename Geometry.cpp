#include "Geometry.hpp"



double Geometry::angle(Vector& a, Vector& b, Vector& c) {

    mpq_class x = a.getX() - b.getX();
    mpq_class y = a.getY() - b.getY();
    mpq_class z = a.getZ() - b.getZ();
    Vector ab(x, y, z);
    x = a.getX() - c.getX();
    y = a.getY() - c.getY();
    z = a.getZ() - c.getZ();
    Vector ac(x, y, z);
    
    // return angle between [0, 180]
    double l1 = ab.length();
    double l2 = ac.length();
    mpq_class dQ = Geometry::dotProduct(ab, ac);
    double dd = dQ.get_d();
    double angle = acosf(dd / (l1 * l2)) / 3.141592f * 180.0f;
    return angle;
}

Vector Geometry::centroid(std::vector<Vector*>& vecs) {

    mpq_class sumX = 0;
    mpq_class sumY = 0;
    mpq_class sumZ = 0;
    for (Vector* v : vecs)
        {
        sumX += v->x;
        sumY += v->y;
        sumZ += v->z;
        }
    mpq_class n = (int)vecs.size();
    mpq_class x = sumX / n;
    mpq_class y = sumY / n;
    mpq_class z = sumZ / n;
    return Vector(x, y, z);
}

Vector Geometry::centroid(std::vector<Vertex*>& vecs) {

    mpq_class sumX = 0;
    mpq_class sumY = 0;
    mpq_class sumZ = 0;
    for (Vector* v : vecs)
        {
        sumX += v->x;
        sumY += v->y;
        sumZ += v->z;
        }
    mpq_class n = (int)vecs.size();
    mpq_class x = sumX / n;
    mpq_class y = sumY / n;
    mpq_class z = sumZ / n;
    return Vector(x, y, z);
}

mpq_class Geometry::distanceSquared(Vector& v1, Vector& v2) {

    mpq_class xDiff = v1.getX() - v2.getX();
    mpq_class yDiff = v1.getY() - v2.getY();
    mpq_class zDiff = v1.getZ() - v2.getZ();
    mpq_class dd = xDiff * xDiff + yDiff * yDiff + zDiff * zDiff;
    return dd;
}

bool Geometry::intersect(Plane& plane1, Plane& plane2, Plane& plane3, Vector& intersection) {

    mpq_class& a1 = plane1.getA();
    mpq_class& b1 = plane1.getB();
    mpq_class& c1 = plane1.getC();
    mpq_class& d1 = plane1.getD();
    mpq_class& a2 = plane2.getA();
    mpq_class& b2 = plane2.getB();
    mpq_class& c2 = plane2.getC();
    mpq_class& d2 = plane2.getD();
    mpq_class& a3 = plane3.getA();
    mpq_class& b3 = plane3.getB();
    mpq_class& c3 = plane3.getC();
    mpq_class& d3 = plane3.getD();

    mpq_class detA  = a1 * (b2 * c3 - c2 * b3) + b1 * (c2 * a3 - a2 * c3) + c1 * (a2 * b3 - b2 * a3);
    if (detA == 0)
        return false;
    mpq_class detAx = -d1 * (b2 * c3 - c2 * b3) - d2 * (b3 * c1 - c3 * b1) - d3 * (b1 * c2 - c1 * b2);
    mpq_class detAy = -d1 * (c2 * a3 - a2 * c3) - d2 * (c3 * a1 - a3 * c1) - d3 * (c1 * a2 - a1 * c2);
    mpq_class detAz = -d1 * (a2 * b3 - b2 * a3) - d2 * (a3 * b1 - b3 * a1) - d3 * (a1 * b2 - b1 * a2);
    
#   if 0
    std::cout << "detA  = " << detA << std::endl;
    std::cout << "detAx = " << detAx << std::endl;
    std::cout << "detAy = " << detAy << std::endl;
    std::cout << "detAz = " << detAz << std::endl;
#   endif
    
    mpq_class x = detAx / detA;
    mpq_class y = detAy / detA;
    mpq_class z = detAz / detA;
    
    intersection.setX(x);
    intersection.setY(y);
    intersection.setZ(z);
    
    return true;
}

bool Geometry::intersect(Plane& plane1, Plane& plane2, Line& intersection) {

    // find direction vector of the intersection line
    mpq_class& a1 = plane1.getA();
    mpq_class& b1 = plane1.getB();
    mpq_class& c1 = plane1.getC();
    mpq_class& d1 = plane1.getD();
    mpq_class& a2 = plane2.getA();
    mpq_class& b2 = plane2.getB();
    mpq_class& c2 = plane2.getC();
    mpq_class& d2 = plane2.getD();
    mpq_class x = b1 * c2 - c1 * b2;
    mpq_class y = a1 * c2 - c1 * a2;
    mpq_class z = a1 * b2 - b1 * a2;
    if (x == 0 && y == 0 && z == 0)
        return false;
    Vector v(x, y, z);
        
    mpq_class dot = x * x + y * y + z * z;  // dot product
    mpq_class xTemp = a1 * d2;
    mpq_class yTemp = b1 * d2;
    mpq_class zTemp = c1 * d2;
    Vector u1(xTemp, yTemp, zTemp);          // d2 * n1
    mpq_class negD1 = -d1;
    xTemp = a2 * negD1;
    yTemp = b2 * negD1;
    zTemp = c2 * negD1;
    Vector u2(xTemp, yTemp, zTemp);           //-d1 * n2
    xTemp = u1.getX() + u2.getX();
    yTemp = u1.getY() + u2.getY();
    zTemp = u1.getZ() + u2.getZ();
    Vector sum(xTemp, yTemp, zTemp);
    xTemp = sum.getY() * v.getZ() - sum.getZ() * v.getY();
    yTemp = sum.getZ() * v.getX() - sum.getX() * v.getZ();
    zTemp = sum.getX() * v.getY() - sum.getY() * v.getX();
    Vector crs(xTemp, yTemp, zTemp);
    xTemp = crs.getX() / dot;
    yTemp = crs.getY() / dot;
    zTemp = crs.getZ() / dot;
    Vector p(xTemp, yTemp, zTemp);       // (d2*N1-d1*N2) X V / V dot V
    
    intersection.setDirection(v);
    intersection.setPoint(p);

    return true;
}

bool Geometry::intersect(Line& line1, Line& line2, Vector& intersection) {

    Vector& p = line1.getPoint();           // P1
    Vector& v = line1.getDirection();       // v
    Vector& q = line2.getPoint();           // Q1
    Vector& u = line2.getDirection();       // u

    // find a = v x u, Vector a = v.cross(u);
    Vector a;
    Geometry::crossProduct(v, u, a);

    // find dot product = (v x u) . (v x u)
    mpq_class dot = Geometry::dotProduct(a, a);              // inner product

    // if v and u are parallel (v x u = 0), then no intersection, return NaN point
    if(dot == 0)
        return false;

    // find b = (Q1-P1) x u, Vector b = (q - p).cross(u);
    mpq_class diffX = q.getX() - p.getX();
    mpq_class diffY = q.getY() - p.getY();
    mpq_class diffZ = q.getZ() - p.getZ();
    Vector diffPQ(diffX, diffY, diffZ);
    Vector b;
    Geometry::crossProduct(diffPQ, u, b);

    // find t = (b.a)/(a.a) = ((Q1-P1) x u) .(v x u) / (v x u) . (v x u), float t = b.dot(a) / dot;
    mpq_class t = Geometry::dotProduct(b, a) / dot;
    
    // find intersection point by substituting t to the line1 eq, Vector point = p + (t * v);
    mpq_class x = p.x + v.getX() * t;
    mpq_class y = p.y + v.getY() * t;
    mpq_class z = p.z + v.getZ() * t;
    intersection.setX( x );
    intersection.setY( y );
    intersection.setZ( z );
    return true;
}

bool Geometry::intersect(Plane& plane, Line& line, Vector& intersection) {

    // from line = p + t * v
    Vector p = line.getPoint();        // (x0, y0, z0)
    Vector v = line.getDirection();    // (x,  y,  z)
    Vector normal;
    plane.normal(normal);

    // dot products
    mpq_class dot1 = Geometry::dotProduct(normal, p); // a*x0 + b*y0 + c*z0
    mpq_class dot2 = Geometry::dotProduct(normal, v); // a*x + b*y + c*z

    // if denominator=0, no intersect
    if (dot2 == 0)
        return false;

    // find t = -(a*x0 + b*y0 + c*z0 + d) / (a*x + b*y + c*z)
    mpq_class t = -(dot1 + plane.getD()) / dot2;

    // find intersection point
    mpq_class x = p.getX();
    mpq_class y = p.getY();
    mpq_class z = p.getZ();
    x += (v.getX() * t);
    y += (v.getY() * t);
    z += (v.getZ() * t);
    intersection.setX(x);
    intersection.setY(y);
    intersection.setZ(z);
    //return p + (v * t);

    return true;
}

bool Geometry::isIntersected(Plane& p, Line& line) {

    // direction vector of line
    Vector& v = line.getDirection();

    // normal of the plane
    Vector normal;
    p.normal(normal);

    // dot product with normal of the plane
    mpq_class dot = Geometry::dotProduct(normal, v);

    if (dot == 0)
        return false;
    return true;
}

void Geometry::crossProduct(Vector& v1, Vector& v2, Vector& res) {

    mpq_class& v1X = v1.getX();
    mpq_class& v1Y = v1.getY();
    mpq_class& v1Z = v1.getZ();
    mpq_class& v2X = v2.getX();
    mpq_class& v2Y = v2.getY();
    mpq_class& v2Z = v2.getZ();
    mpq_class x = v1Y * v2Z - v1Z * v2Y;
    mpq_class y = v1Z * v2X - v1X * v2Z;
    mpq_class z = v1X * v2Y - v1Y * v2X;
    res.setX(x);
    res.setY(y);
    res.setZ(z);
}

mpq_class Geometry::dotProduct(Vector& v1, Vector& v2) {

    mpq_class d = v1.getX() * v2.getX() + v1.getY() * v2.getY() + v1.getZ() * v2.getZ();
    return d;
}
