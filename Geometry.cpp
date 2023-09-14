#include "Geometry.hpp"



double Geometry::angle(Vector& a, Vector& b, Vector& c, Vector& n) {

    mpq_class x = b.getX() - a.getX();
    mpq_class y = b.getY() - a.getY();
    mpq_class z = b.getZ() - a.getZ();
    Vector ab(x, y, z);
    x = c.getX() - a.getX();
    y = c.getY() - a.getY();
    z = c.getZ() - a.getZ();
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
