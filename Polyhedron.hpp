#ifndef Polyhedron_hpp
#define Polyhedron_hpp

#include <map>
#include <vector>
#include <gmpxx.h>
#include <string>
#include "Plane.hpp"
class RateMatrix;
class Vector;
class Vertex;

class MpqMatrix {

    public:
        mpq_class&          operator()(size_t r, size_t c) { return this->m[r * 4 + c]; }
        const mpq_class&    operator()(size_t r, size_t c) const { return this->m[r * 4 + c]; }
        void                print(void);
        void                set(int idx, Vector* v);
        
    private:
        mpq_class           m[16];
};

struct VectorInfo {
    
    mpq_class               volume;
    double                  alphaC;
};

typedef std::map<Plane*,std::vector<Vertex*> > plane_vertex_map;
typedef std::map< std::pair<Plane*,Plane*>, std::vector<Vertex*> > line_vertex_map;
typedef std::map<Vector*,VectorInfo> vector_volume_map;



class Polyhedron {

    /**
     * This class is used to calculate the probability of proposing a non-reversible model from
     * a time reversible model for reversible jump MCMC. Proposing a non-reversible model from
     * a time reversible one involves the generation of three uniform(0,1) random variables,
     * u1, u2, and u3. The probability density would appear to be 1 X 1 X 1 = 1, but this is not
     * the case because there are constraints on the values that result in a valid non-reversible
     * model (e.g., one with all q_{ij} >= 0 (i != j). The constraints are (in Latex form):
     *
     * 0 \leq w_{AC} + w_{CG} (2 u_1 - 1) + w_{CT}(2 u_2 - 1) \leq 2 w_{AC}
     * 0 \leq w_{AG} - w_{CG} (2 u_1 - 1) + w_{GT}(2 u_3 - 1) \leq 2 w_{AG}
     * 0 \leq w_{AT} - w_{CT} (2 u_2 - 1) - w_{GT}(2 u_3 - 1) \leq 2 w_{AT}
     *
     * where w_{ij} = \pi_i q_{ij} is the average rate of changing from nucleotide i to j. The
     * class implements the constraints as planes (in 3D space) with the coordinates being (u1,u2,u3).
     * Six of the planes represent the bounds of the uniform(0,1) random variables (i.e., they
     * represent the unit cube). The other six planes represent the six constraints, above. The
     * planes for the unit cube never change and are initialized on instantiation of a Polyhedron
     * object. The other six planes (from the constraints, above) are initialized each time the
     * lnProbability function is called. In fact, the lnProbability functions take as a parameter
     * a vector of weights: w_{AC}, w_{AG}, w_{AT}, w_{CG}, w_{CT}, w_{GT}.
     *
     * The constraints form a polyhedron, with the details depending on the weights. The goal
     * is to calculate the volume of this polyhedron because the inverse of the volume is
     * the probability of proposing the non-reversible model from the time reversible model,
     * with the time reversible model being the point in the very center of the cube with
     * coordinates (1/2, 1/2, 1/2).
     * 
     * This mixture can be considered as a multinomial distribution. We specify a vector of probabilities
     * and a vector of values. Then, a value drawn from this distribution takes each value corresponding to
     * its probability.
     * The values are already of the correct mixture type. You may want to apply a mixture allocation move
     * to change between the current value. The values themselves change automatically when the input parameters change.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (John Huelsenbeck)
     * @since 2014-11-18, version 1.0
     */

    public:
                            Polyhedron(void);
                            Polyhedron(const Polyhedron& p) = delete;
                           ~Polyhedron(void);
        void                setAlphaT(double x) { alphaT = x; }
        double              lnProbabilityForward(std::vector<mpq_class>& W, Vector& pt);
        double              lnProbabilityReverse(std::vector<mpq_class>& W, Vector& pt);
        double              monteCarloVolume(int numberReplicates);
        void                pointFromNonreversibleWeights(RateMatrix& Q, std::vector<mpq_class> pi, Vector& pt);
        void                print(std::vector<mpq_class>& W);
    
    private:
        void                calculateTetrahedronVolume(Vector* v1, Vector* v2, Vector* v3, mpq_class& vol);
        void                clearTetrahedraMap(void);
        void                computeLandU(MpqMatrix& aMat, MpqMatrix& lMat, MpqMatrix& uMat);
        mpq_class           det(MpqMatrix& m);
        void                facetVolume(Plane* pln, std::vector<Vertex*>& vertices, mpq_class& vol);
        Vertex*             findOtherVertex(Vertex* from, Vertex* v, Plane* pln);
        void                initializeFacets(void);
        void                initializePlanes(void);
        void                insertVertex(Plane* p1, Plane* p2, Vertex* v);
        bool                intersect(Plane& plane1, Plane& plane2, Plane& plane3, Vector& intersection);
        bool                isInTetrahedron(Vector* pt, Vector* center, Vector* v1, Vector* v2, Vector* v3,
                            mpq_class& b1, mpq_class& b2, mpq_class& b3, mpq_class& b4);
        bool                isValid(Vector& pt);
        std::string         mathematicaPolyhedronOutput(void);
        void                sampleTetrahedron(Plane* pln, Vector* center, Vector* v1, Vector* v2, Vector* v3, Vector& pt, VectorInfo& info);
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
        bool                pointFoundInPolyhedron;
        vector_volume_map   tetrahedra;
        Vector              randomPoint;
        double              alphaT;
        double              alphaC;
        
        bool                mathematicaPolyhedron;
        
        mpq_class           sumJacobians;
};

#endif
