#ifndef Polyhedral_hpp
#define Polyhedral_hpp

#include <map>
#include <vector>
#include <gmpxx.h>
#include "Geometry.hpp"
#include "Plane.hpp"
class Facet;
class Vector;
class Vertex;


typedef std::map<Plane*,std::vector<Vertex*>> plain_vertex_map;
typedef std::map< std::pair<Plane*,Plane*>, std::vector<Vertex*>> line_vertex_map;



class Polyhedral {


    public:
                            Polyhedral(void);
                            Polyhedral(const Polyhedral& p) = delete;
        double              monteCarloVolume(int numberReplicates);
        mpf_class           volume(std::vector<mpq_class>& W);
        mpf_class           volume(std::vector<mpq_class>& W, Vector& pt);
    
    private:
        mpq_class           facetArea(Vertex* first, Plane* pln);
        Vertex*             findOtherVertex(Vertex* from, line_vertex_map& linesMap, Vertex* v, Plane* pln);
        void                initializeFacets(plain_vertex_map& verticesMap, line_vertex_map& linesMap);
        void                initializePlanes(void);
        void                insertVertex(line_vertex_map& linesMap, Plane* p1, Plane* p2, Vertex* v);
        bool                isValid(Vector& pt);
        void                samplePolytope(Vector& pt);
        void                setWeights(std::vector<mpq_class>& W);
        
        mpq_class           wAC;
        mpq_class           wAG;
        mpq_class           wAT;
        mpq_class           wCG;
        mpq_class           wCT;
        mpq_class           wGT;
        
        mpq_class           minX;
        mpq_class           minY;
        mpq_class           minZ;
        mpq_class           maxX;
        mpq_class           maxY;
        mpq_class           maxZ;
        mpq_class           diffX;
        mpq_class           diffY;
        mpq_class           diffZ;
        
        mpq_class           xzMaxA;
        mpq_class           xzMaxB;
        mpq_class           xzMinA;
        mpq_class           xzMinB;
        mpq_class           xyMaxA;
        mpq_class           xyMaxB;
        mpq_class           xyMinA;
        mpq_class           xyMinB;
        mpq_class           yzMaxA;
        mpq_class           yzMaxB;
        mpq_class           yzMinA;
        mpq_class           yzMinB;

        mpq_class           zeroQ;
        mpq_class           oneQ;
        mpq_class           oneHalfQ;
        mpq_class           twoQ;
        
        Vector              center;
        Vector              xzMinA_Zero_Zero;
        Vector              xzMaxA_Zero_Zero;
        Vector              xzMinA_One_Zero;
        Vector              xzMaxA_One_Zero;
        Vector              xzMinB_Zero_One;
        Vector              xzMaxB_Zero_One;
        Vector              zero_xyMaxA_Zero;
        Vector              zero_xyMinA_Zero;
        Vector              one_xyMaxB_Zero;
        Vector              one_xyMinB_Zero;
        Vector              zero_xyMaxA_One;
        Vector              zero_xyMinA_One;
        Vector              zero_Zero_yzMaxA;
        Vector              zero_Zero_yzMinA;
        Vector              zero_One_yzMaxB;
        Vector              zero_One_yzMinB;
        Vector              one_Zero_yzMaxA;
        Vector              one_Zero_yzMinA;

        Plane               front;
        Plane               back;
        Plane               top;
        Plane               bottom;
        Plane               left;
        Plane               right;
        Plane               xz1;
        Plane               xz2;
        Plane               xy1;
        Plane               xy2;
        Plane               yz1;
        Plane               yz2;

        std::vector<Plane*> planes;
        plain_vertex_map    verticesMap;
        line_vertex_map     linesMap;

        bool                computeExtrema;
        mpf_class           polytopeVolume;
};

#endif
