#ifndef Polyhedral_hpp
#define Polyhedral_hpp

#include <vector>
#include <gmpxx.h>
#include "Geometry.hpp"
class Vector;



class Polyhedral {

    public:
                    Polyhedral(void) = delete;
                    Polyhedral(std::vector<mpq_class>& W);
    
    private:
        void        initializePlanes(void);
        bool        isLineInList(Line& x, std::vector<Line>& lines);
        bool        isValid(Vector& pt);
        mpq_class   wAC;
        mpq_class   wAG;
        mpq_class   wAT;
        mpq_class   wCG;
        mpq_class   wCT;
        mpq_class   wGT;
        mpq_class   minX;
        mpq_class   minY;
        mpq_class   minZ;
        mpq_class   maxX;
        mpq_class   maxY;
        mpq_class   maxZ;
};

#endif
