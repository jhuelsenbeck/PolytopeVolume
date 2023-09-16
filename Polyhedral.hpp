#ifndef Polyhedral_hpp
#define Polyhedral_hpp

#include <map>
#include <vector>
#include <gmpxx.h>
#include <string>
#include "Geometry.hpp"
#include "Plane.hpp"
class Facet;
class Vector;
class Vertex;

typedef std::map<Plane*,std::vector<Vertex*>> plane_vertex_map;
typedef std::map< std::pair<Plane*,Plane*>, std::vector<Vertex*>> line_vertex_map;
typedef std::map<Vector*,mpf_class> vector_volume_map;



class Polyhedral {

    public:
                            Polyhedral(void);
                            Polyhedral(const Polyhedral& p) = delete;
        double              monteCarloVolume(int numberReplicates);
        mpf_class           volume(std::vector<mpq_class>& W);
        mpf_class           volume(std::vector<mpq_class>& W, Vector& pt);
    
    private:
        void                calculateTetrahedronVolume(Vector* v1, Vector* v2, Vector* v3, mpf_class& d, mpf_class& vol);
        mpq_class           facetArea(Vertex* first, Plane* pln);
        Vertex*             findOtherVertex(Vertex* from, Vertex* v, Plane* pln);
        void                initializeFacets(void);
        void                initializePlanes(void);
        void                insertVertex(Plane* p1, Plane* p2, Vertex* v);
        bool                isValid(Vector& pt);
        void                sampleTetrahedra(std::vector<Vertex*>& vertices, mpf_class& d);
        void                sampleTetrahedron(Vector* center, Vector* v1, Vector* v2, Vector* v3, Vector& pt);
        void                setWeights(std::vector<mpq_class>& W);
        
        mpq_class           wAC;
        mpq_class           wAG;
        mpq_class           wAT;
        mpq_class           wCG;
        mpq_class           wCT;
        mpq_class           wGT;
                
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
        plane_vertex_map    verticesMap;
        line_vertex_map     linesMap;

        bool                randomlySample;
        vector_volume_map   tetrahedra;
        Vector              randomPoint;
        
        mpf_class           polytopeVolume;
};

#endif
