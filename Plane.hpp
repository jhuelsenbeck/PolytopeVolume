#ifndef PLANE_H
#define PLANE_H

#include <gmpxx.h>
#include "Vector.hpp"

class Plane {

    public:
                        Plane(void);
                        Plane(Vector pt1, Vector pt2, Vector pt3);
                        Plane(Vector& pt1, Vector& pt2, Vector& pt3);
        bool            operator==(const Plane& rhs) const;
        bool            operator<(const Plane& rhs) const;
        mpq_class&      getA(void) { return a; }
        mpq_class&      getB(void) { return b; }
        mpq_class&      getC(void) { return c; }
        mpq_class&      getD(void) { return d; }
        mpf_class       getDistance(Vector& point);
        std::string     getStr(void);
        void            normal(Vector& n);
    
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

#endif
