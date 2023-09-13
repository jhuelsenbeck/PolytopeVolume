#ifndef Polyhedral_hpp
#define Polyhedral_hpp

#include <map>
#include <vector>
#include <gmpxx.h>
#include "Geometry.hpp"
class Facet;
class Vector;
class Vertex;


typedef std::map<Plane*,std::vector<Vertex*>> plain_vertex_map;
typedef std::map< std::pair<Plane*,Plane*>, std::vector<Vertex*>> line_vertex_map;



class Polyhedral {


    public:
                    Polyhedral(void) = delete;
                    Polyhedral(std::vector<mpq_class>& W);
    
    private:
        mpq_class   facetArea(Vertex* first, Plane* pln);
        Vertex*     findOtherVertex(Vertex* from, line_vertex_map& linesMap, Vertex* v, Plane* pln);
        void        initializeFacets(plain_vertex_map& verticesMap, line_vertex_map& linesMap);
        void        initializePlanes(void);
        void        insertVertex(line_vertex_map& linesMap, Plane* p1, Plane* p2, Vertex& v);
        bool        isLineInList(Line& x, std::vector<Line>& lines);
        bool        isValid(Vector& pt);
        double      monteCarloVolume(int numberReplicates);
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
        mpf_class   volume;
        std::vector<Facet*> facets;
};

#endif
