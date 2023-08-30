#ifndef Facet_H
#define Facet_H

#include <set>
#include <string>
#include <vector>
#include "Plane.hpp"
#include "Vertex.hpp"



class Facet : public Plane {

    public:
                                Facet(void);
                                Facet(const Plane& p);
                                Facet(const Plane& p, std::vector<Vertex*>& verts);
                               ~Facet(void);
        void                    addVertex(Vertex* p);
        void                    addVertex(std::vector<Vertex*> verts);
        mpq_class               area(void);
        Vector                  centroid(void);
        Vertex*                 getFirstVertex(void) { return firstVertex; }
        std::string             getString(void);
        std::vector<Vertex*>    getVertices(void);
        int                     numVertices(void) { return (int)vertices.size(); }
        void                    print(void);
        void                    print(int idx);
        void                    removeVerticesBehindPlane(Plane& p);
        mpf_class               volume(Vertex p);

    private:
        double                  findArea(int n, double* x, double* y);
        Vertex                  findNormal(int n, double* x, double* y, double* z);
        void                    removeVertex(Vertex* v);
        std::vector<Vertex*>    vertices;
        Vertex*                 firstVertex;

        double                  facetArea;
        Vertex                  facetNormal;
        bool                    areaCalculated;
        bool                    normalCalculated;
};

#endif
