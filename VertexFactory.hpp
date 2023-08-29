#ifndef VertexFactory_H
#define VertexFactory_H

#include <set>
#include <vector>
#include "Vector.hpp"

class VertexFactory {

    public:
        static VertexFactory&   vertexFactoryInstance(void)
                                    {
                                    static VertexFactory singleNodeFactory;
                                    return singleNodeFactory;
                                    }
        void                    drainPool(void);
        Vertex*                 getVertex(void);
        Vertex*                 getVertex(Vertex& v);
        Vertex*                 getVertex(Vector& v);
        void                    returnToPool(Vertex* nde);

    private:
                                VertexFactory(void);
                                VertexFactory(const VertexFactory&);
                                VertexFactory& operator=(const VertexFactory&);
                               ~VertexFactory(void);
        std::vector<Vertex*>    vertexPool;
        std::set<Vertex*>       allocatedVertices;
};

#endif
