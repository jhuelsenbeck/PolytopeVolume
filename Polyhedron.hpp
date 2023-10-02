#ifndef Polyhedron_hpp
#define Polyhedron_hpp

#include <map>
#include <vector>
#include <gmpxx.h>
#include <string>
#include "Plane.hpp"
class Vector;
class Vertex;

typedef std::map<Plane*,std::vector<Vertex*> > plane_vertex_map;
typedef std::map< std::pair<Plane*,Plane*>, std::vector<Vertex*> > line_vertex_map;
typedef std::map<Vector*,mpq_class> vector_volume_map;



class Polyhedron {

    public:
                            Polyhedron(void);
                            Polyhedron(const Polyhedron& p) = delete;
                           ~Polyhedron(void);
        double              monteCarloVolume(int numberReplicates);
        void                print(std::vector<mpq_class>& W);
        mpq_class           volume(std::vector<mpq_class>& W);
        mpq_class           volume(std::vector<mpq_class>& W, Vector& pt);
        mpq_class           volume(std::vector<mpq_class>& W, Vector& pt, double fac);
    
    private:
        void                calculateTetrahedronVolume(Vector* v1, Vector* v2, Vector* v3, mpq_class& vol);
        void                clearTetrahedraMap(void);
        void                facetVolume(std::vector<Vertex*>& vertices, mpq_class& vol);
        Vertex*             findOtherVertex(Vertex* from, Vertex* v, Plane* pln);
        void                initializeFacets(void);
        void                initializePlanes(void);
        void                insertVertex(Plane* p1, Plane* p2, Vertex* v);
        bool                intersect(Plane& plane1, Plane& plane2, Plane& plane3, Vector& intersection);
        bool                isValid(Vector& pt);
        std::string         mathematicaPolyhedronOutput(void);
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

        mpq_class           aQ;       // used in calculateTetrahedronVolume
        mpq_class           bQ;
        mpq_class           cQ;
        mpq_class           dQ;
        mpq_class           eQ;
        mpq_class           fQ;
        mpq_class           gQ;
        mpq_class           hQ;
        mpq_class           iQ;

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
        
        bool                mathematicaPolyhedron;
        
        mpq_class           polytopeVolume;
};

#endif
