#ifndef Vector_hpp
#define Vector_hpp

#include <gmpxx.h>
#include <iostream>

class Vector {

    public:
                        Vector(void);
                        Vector(int xI, int yI, int zI);
                        Vector(double xD, double yD, double zD);
                        Vector(mpq_class& xQ, mpq_class yQ, mpq_class zQ);
        bool            operator==(Vector& rhs);
        bool            operator!=(Vector& rhs);
        Vector          cross(Vector& rhs);
        mpq_class       distanceSquared(const Vector& vec) const;
        mpq_class&      getX(void) { return x; }
        mpq_class&      getY(void) { return y; }
        mpq_class&      getZ(void) { return z; }
        std::string     getStr(void);
        double          length(void);
        void            normalize(void);
        void            setX(mpq_class& xQ) { x = xQ; }
        void            setY(mpq_class& yQ) { y = yQ; }
        void            setZ(mpq_class& zQ) { z = zQ; }
        mpq_class       x;
        mpq_class       y;
        mpq_class       z;

    friend std::ostream& operator<<(std::ostream& os, const Vector& pt);
};

inline std::ostream& operator<<(std::ostream& os, const Vector& pt) {

    os << "{" << pt.x << ", " << pt.y << ", " << pt.z << "}";
    return os;
}

#endif
