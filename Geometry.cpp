#include "Geometry.hpp"



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
