#ifndef LINE_H
#define LINE_H

#include <cmath>
#include <iostream>
#include "Vector.hpp"
#include "Vertex.hpp"


class Line {

    public:
                        Line(void);
                        Line(const Vector& v, const Vector& p) : direction(v), point(p) { }
        bool            operator==(Line& rhs);
        mpq_class       distanceToPoint(Vertex& p);
        Vector&         getDirection(void) { return direction; }
        Vector&         getPoint(void) { return point; }
        void            setDirection(Vector& p) { direction = p; }
        void            setPoint(Vector& p) { point = p; }
    
    private:
        Vector           direction;
        Vector           point;

    friend std::ostream& operator<<(std::ostream& os, const Line& lne);
};

inline std::ostream& operator<<(std::ostream& os, const Line& lne) {

    os << "{" << lne.point << ", " << lne.direction << "}";
    return os;
}


#endif

