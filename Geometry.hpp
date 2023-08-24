#ifndef Geometry_hpp
#define Geometry_hpp

#include <gmpxx.h>
#include <iostream>
#include <string>



class Point {

    public:
                        Point(void);
                        Point(int xI, int yI, int zI);
                        Point(double xD, double yD, double zD);
                        Point(mpq_class& xQ, mpq_class yQ, mpq_class zQ);
        mpq_class&      getX(void) { return x; }
        mpq_class&      getY(void) { return y; }
        mpq_class&      getZ(void) { return z; }
        std::string     getStr(void);
        void            setX(mpq_class& xQ) { x = xQ; }
        void            setY(mpq_class& yQ) { y = yQ; }
        void            setZ(mpq_class& zQ) { z = zQ; }
        
    private:
        mpq_class       x;
        mpq_class       y;
        mpq_class       z;

    friend std::ostream& operator<<(std::ostream& os, const Point& pt);
};

inline std::ostream& operator<<(std::ostream& os, const Point& pt) {

    os << "{" << pt.x << ", " << pt.y << ", " << pt.z << "}";
    return os;
}

class Plane {

    public:
                        Plane(void) = delete;
                        Plane(Point pt1, Point pt2, Point pt3);
                        Plane(Point& pt1, Point& pt2, Point& pt3);
        mpq_class&      getA(void) { return a; }
        mpq_class&      getB(void) { return b; }
        mpq_class&      getC(void) { return c; }
        mpq_class&      getD(void) { return d; }
    
    private:
        mpq_class       a;
        mpq_class       b;
        mpq_class       c;
        mpq_class       d;

    friend std::ostream& operator<<(std::ostream& os, const Plane& pl);
};

inline std::ostream& operator<<(std::ostream& os, const Plane& pl) {

    os << "{" << pl.a << ", " << pl.b << ", " << pl.c << ", " << pl.d << "}";
    return os;
}

namespace Geometry {

    bool intersect(Plane& plane1, Plane& plane2, Plane& plane3, Point& intersection);
}

#endif
