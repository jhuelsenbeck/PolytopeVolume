#include "Geometry.hpp"



Point::Point(void) {

    this->x = 0;
    this->y = 0;
    this->z = 0;
}

Point::Point(int xI, int yI, int zI) {

    this->x = xI;
    this->y = yI;
    this->z = zI;
}

Point::Point(double xD, double yD, double zD) {

    this->x = xD;
    this->y = yD;
    this->z = zD;
}

Point::Point(mpq_class& xQ, mpq_class yQ, mpq_class zQ) {

    this->x = xQ;
    this->y = yQ;
    this->z = zQ;
}

std::string Point::getStr(void) {

    std::string str = "(";
    str += std::to_string(x.get_d()) + ", " + std::to_string(y.get_d()) + ", " + std::to_string(z.get_d());
    str += ")";
    return str;
}

Plane::Plane(Point pt1, Point pt2, Point pt3) {

    mpq_class& x1 = pt1.getX();
    mpq_class& y1 = pt1.getY();
    mpq_class& z1 = pt1.getZ();
    mpq_class& x2 = pt2.getX();
    mpq_class& y2 = pt2.getY();
    mpq_class& z2 = pt2.getZ();
    mpq_class& x3 = pt3.getX();
    mpq_class& y3 = pt3.getY();
    mpq_class& z3 = pt3.getZ();

    mpq_class a1 = x2 - x1;
    mpq_class b1 = y2 - y1;
    mpq_class c1 = z2 - z1;
    mpq_class a2 = x3 - x1;
    mpq_class b2 = y3 - y1;
    mpq_class c2 = z3 - z1;
    this->a = b1 * c2 - b2 * c1;
    this->b = a2 * c1 - a1 * c2;
    this->c = a1 * b2 - b1 * a2;
    this->d = (- a * x1 - b * y1 - c * z1);
}

Plane::Plane(Point& pt1, Point& pt2, Point& pt3) {

    mpq_class& x1 = pt1.getX();
    mpq_class& y1 = pt1.getY();
    mpq_class& z1 = pt1.getZ();
    mpq_class& x2 = pt2.getX();
    mpq_class& y2 = pt2.getY();
    mpq_class& z2 = pt2.getZ();
    mpq_class& x3 = pt3.getX();
    mpq_class& y3 = pt3.getY();
    mpq_class& z3 = pt3.getZ();

    mpq_class a1 = x2 - x1;
    mpq_class b1 = y2 - y1;
    mpq_class c1 = z2 - z1;
    mpq_class a2 = x3 - x1;
    mpq_class b2 = y3 - y1;
    mpq_class c2 = z3 - z1;
    this->a = b1 * c2 - b2 * c1;
    this->b = a2 * c1 - a1 * c2;
    this->c = a1 * b2 - b1 * a2;
    this->d = (- a * x1 - b * y1 - c * z1);
}

bool Geometry::intersect(Plane& plane1, Plane& plane2, Plane& plane3, Point& intersection) {

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
