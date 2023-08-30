#include <cmath>
#include <map>
#include "Edge.hpp"
#include "Facet.hpp"
#include "Geometry.hpp"
#include "VertexFactory.hpp"



Facet::Facet(void) : Plane() {

    firstVertex = NULL;
    areaCalculated = false;
}

Facet::Facet(const Plane& p) : Plane(p) {

    areaCalculated = false;
}

Facet::Facet(const Plane& p, std::vector<Vertex*>& verts) : Plane(p) {

    vertices = verts;
    firstVertex = vertices[0];
    areaCalculated = false;
}

Facet::~Facet(void) {

    VertexFactory& vFactory = VertexFactory::vertexFactoryInstance();
    for (Vertex* v : vertices)
        vFactory.returnToPool(v);
    vertices.clear();
}

void Facet::addVertex(Vertex* p) {

    vertices.push_back(p);
}

void Facet::addVertex(std::vector<Vertex*> verts) {

    for (Vertex* v : verts)
        vertices.push_back(v);
}

mpq_class Facet::area(void) {

    if (areaCalculated == true)
        return facetArea;

    Vector unitNormal;
    this->normal(unitNormal);
    unitNormal.normalize();

    Vertex* f = firstVertex;
    Vertex* p = f;
    Vector sumCrossProducts;
    mpq_class x = sumCrossProducts.getX();
    mpq_class y = sumCrossProducts.getY();
    mpq_class z = sumCrossProducts.getZ();
    do
        {
        Vector crs = p->cross(*(p->getTo()));
        x += crs.getX();
        y += crs.getY();
        z += crs.getZ();
        //sumCrossProducts += crs;

        p = p->getTo();
        } while (p != f);
    mpq_class dot = Geometry::dotProduct(unitNormal, sumCrossProducts);
    if (dot < 0)
        dot = -dot;
    mpq_class facetArea = dot / 2;
    areaCalculated = true;
    
    return facetArea;
}

Vector Facet::centroid(void) {

    mpq_class sumX = 0;
    mpq_class sumY = 0;
    mpq_class sumZ = 0;
    for (Vertex* v : vertices)
        {
        sumX += v->x;
        sumY += v->y;
        sumZ += v->z;
        }
    int n = (int)vertices.size();
    sumX /= n;
    sumY /= n;
    sumZ /= n;
    return Vector( sumX, sumY, sumZ );
}

double Facet::findArea(int n, double* x, double* y) {

#   if 0
    // guarantee the first two vertices are also at array end
    x[n] = x[0];
    y[n] = y[0];
    x[n+1] = x[1];
    y[n+1] = y[1];
    double sum = 0.0;
    double *xptr = x+1, *ylow = y, *yhigh = y+2;
    for (int i=1; i <= n; i++)
        {
        sum += (*xptr++) * ( (*yhigh++) - (*ylow++) );
        }
    return (sum / 2.0);
#   else
    return 0.0;
#   endif
}

Vertex Facet::findNormal(int n, double* x, double* y, double* z) {

    if (normalCalculated == true)
        return facetNormal;
    
    // get the Newell normal
    double nwx = findArea(n, y, z);
    double nwy = findArea(n, z, x);
    double nwz = findArea(n, x, y);
    
    // get length of the Newell normal
    double nlen = sqrt( nwx*nwx + nwy*nwy + nwz*nwz );
    
    // compute the unit normal
    facetNormal.x = nwx / nlen;
    facetNormal.y = nwy / nlen;
    facetNormal.z = nwz / nlen;
    normalCalculated = true;
    
    return facetNormal;
}

std::string Facet::getString(void) {

    std::string str = "";

    Vertex* f = firstVertex;
    Vertex* p = f;
    str += "{";
    do
        {
        if (p != f)
            str += ",";
        std::string strX = p->x.get_str();
        std::string strY = p->y.get_str();
        std::string strZ = p->z.get_str();
        str += "{" + strX + ", " + strY + ", " + strZ + "}";
        p = p->getTo();
        } while (p != f);
    str += "}";

    return str;
}

std::vector<Vertex*> Facet::getVertices(void) {

    std::vector<Vertex*> sortedVertices;
    
    if (firstVertex == NULL)
        return sortedVertices;
    
    Vertex* v = firstVertex;
    do
        {
        sortedVertices.push_back(v);
        v = v->getTo();
        } while (v != firstVertex);

    return sortedVertices;
}

void Facet::print(void) {

    Vertex* f = firstVertex;
    Vertex* p = f;
    std::cout << "{";
    do
        {
        if (p != f)
            std::cout << ",";
        std::cout << *p;
        p = p->getTo();
        } while (p != f);
    std::cout << "}" << std::endl;;
}

void Facet::print(int idx) {

    std::cout << "f" << idx << " = ";
    print();
}

void Facet::removeVerticesBehindPlane(Plane& p) {

    mpq_class oneHalf = 1/2;
    Vertex midPoint(oneHalf, oneHalf, oneHalf);
    
    std::vector<Vertex*> verticesToRemove;
    for (Vertex* v : vertices)
        {
        Edge edge(&midPoint,v);

        //if ( p.isIntersected(edge) == true)
        if ( Geometry::isIntersected(p, edge) == true)
            {
            Vector i;
            Geometry::intersect(p, edge, i);
            
            mpq_class dMiSq = Geometry::distanceSquared(i, midPoint);
            mpq_class dViSq = Geometry::distanceSquared(i, *v);
            mpq_class& edSq = edge.getDistanceSquared();
            
            mpf_class dMi(dMiSq, 1000);
            mpf_class dVi(dViSq, 1000);
            mpf_class ed(edSq, 1000);
            dMi = sqrt(dMi);
            dVi = sqrt(dVi);
            ed  = sqrt(ed);
            mpf_class diff = dMi + dVi - ed;
            if (diff < 0.0)
                diff = -diff;
            if (abs(diff) < 10e-8)
                verticesToRemove.push_back(v);
            }
        }
    
    for (int i=0; i<verticesToRemove.size(); i++)
        removeVertex(verticesToRemove[i]);
}

void Facet::removeVertex(Vertex* v) {

    VertexFactory& vFactory = VertexFactory::vertexFactoryInstance();
    std::remove(vertices.begin(), vertices.end(), v);
    vFactory.returnToPool(v);
}

mpf_class Facet::volume(Vertex p) {

    mpf_class d = getDistance(p);
    mpf_class a = area();
    mpf_class v = d * a / 3.0;
    return v;
}

