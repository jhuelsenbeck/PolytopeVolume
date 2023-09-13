#include <map>
#include "Facet.hpp"
#include "Geometry.hpp"
#include "Msg.hpp"
#include "Polyhedral.hpp"
#include "RandomVariable.hpp"
#include "VertexFactory.hpp"


Polyhedral::Polyhedral(std::vector<mpq_class>& W) {

    // extract symbols from W
    wAC = W[0];
    wAG = W[1];
    wAT = W[2];
    wCG = W[3];
    wCT = W[4];
    wGT = W[5];
    
    initializePlanes();
}

mpq_class Polyhedral::facetArea(Vertex* first, Plane* pln) {

    Vector unitNormal;
    pln->normal(unitNormal);
    unitNormal.normalize();

    mpq_class normalLengthRational = (unitNormal.getX() * unitNormal.getX()) + (unitNormal.getY() * unitNormal.getY()) + (unitNormal.getZ() * unitNormal.getZ());
    mpf_class normalLengthReal = normalLengthRational;
    normalLengthReal = sqrt(normalLengthReal);
    normalLengthRational = normalLengthReal;
    mpq_class normalLengthInv = 1 / normalLengthRational;
    mpq_class& x = unitNormal.getX();
    mpq_class& y = unitNormal.getY();
    mpq_class& z = unitNormal.getZ();
    x *= normalLengthInv;
    y *= normalLengthInv;
    z *= normalLengthInv;

    Vertex* v = first;
    Vector sumCrossProducts;
    Vector cross;
    do
        {
        std::cout << v << " " << v->getTo() << std::endl;
        Geometry::crossProduct(*v, *(v->getTo()), cross);
        sumCrossProducts += cross;
        v = v->getTo();
        } while (v != first);
    mpq_class dot = Geometry::dotProduct(unitNormal, sumCrossProducts);
    if (dot < 0)
        dot = -dot;
    mpq_class facetArea = dot / 2;

    return facetArea;
}

Vertex* Polyhedral::findOtherVertex(Vertex* from, std::map< std::pair<Plane*,Plane*>, std::vector<Vertex*> >& linesMap, Vertex* v, Plane* pln) {

    for (auto lne : linesMap)
        {
        if (lne.first.first == pln || lne.first.second == pln)
            {
            bool inList = false;
            for (int i=0, n=(int)lne.second.size(); i<n; i++)
                {
                if (lne.second[i] == v)
                    {
                    inList = true;
                    break;
                    }
                }
                
            if (inList == true)
                {
                if (lne.second[0] != from && lne.second[1] != from)
                    {
                    if (v == lne.second[0])
                        return lne.second[1];
                    else
                        return lne.second[0];
                    }
                }
            
            }
        }
    return nullptr;
}

void Polyhedral::initializeFacets(std::map<Plane*,std::vector<Vertex*>>& verticesMap, std::map< std::pair<Plane*,Plane*>, std::vector<Vertex*> >& linesMap) {

    mpq_class oneHalf = 1;
    oneHalf /= 2;
    Vector center(oneHalf, oneHalf, oneHalf);
    volume = 0.0;
    for (auto pln : verticesMap)
        {
        mpf_class distance = pln.first->getDistance(center);

        // clean vertices
        for (int i=0, n=(int)pln.second.size(); i<n; i++)
            {
            pln.second[i]->setTo(nullptr);
            pln.second[i]->setFrom(nullptr);
            }
            
        // order the vertices
        Vertex* first = pln.second[0];
        Vertex* v = first;
        do {
            Vertex* nextV = findOtherVertex(v->getFrom(), linesMap, v, pln.first);
            v->setTo(nextV);
            nextV->setFrom(v);
            v = nextV;
            } while (v != first);
        
        // get the facet area
        mpq_class facetBaseArea = facetArea(first, pln.first);
        std::cout << pln.first << " -- ";
        v = first;
        do
            {
            std::cout << v << " -> ";
            v = v->getTo();
            } while (v != first);
            
        std::cout << " (" << facetBaseArea.get_d() << ", " << distance << ", ";
        // calculate the pyramid volume
        facetBaseArea /= 3;
        mpf_class facetBaseAreaF = facetBaseArea;
        volume += facetBaseAreaF * distance;
        
        std::cout << volume << ")" << std::endl;
        }
}

void Polyhedral::initializePlanes(void) {

    // make planes that will slice up the cube
    // u1 -> (-wAG + wCG - wGT + 2.0 * u3 * wGT) / (2.0 * wCG)   max in x,z
    // u1 -> ( wAG + wCG - wGT + 2.0 * u3 * wGT) / (2.0 * wCG)   min in x,z
    // u2 -> ( wAC + wCG - 2.0 * u1 * wCG + wCT) / (2.0 * wCT)   max in x,y
    // u2 -> (-wAC + wCG - 2.0 * u1 * wCG + wCT) / (2.0 * wCT)   min in x,y
    // u3 -> (-wAT + wCT - 2.0 * u2 * wCT + wGT) / (2.0 * wGT)   max in y,z
    // u3 -> ( wAT + wCT - 2.0 * u2 * wCT + wGT) / (2.0 * wGT)   min in y,z
    
    mpq_class zeroQ = 0;
    mpq_class oneQ  = 1;
    mpq_class twoQ  = 2;
    
    mpq_class xzMaxA = (-wAG + wCG - wGT + twoQ * zeroQ * wGT) / (twoQ * wCG);
    mpq_class xzMaxB = (-wAG + wCG - wGT + twoQ * oneQ  * wGT) / (twoQ * wCG);
    mpq_class xzMinA = ( wAG + wCG - wGT + twoQ * zeroQ * wGT) / (twoQ * wCG);
    mpq_class xzMinB = ( wAG + wCG - wGT + twoQ * oneQ  * wGT) / (twoQ * wCG);
    
    mpq_class xyMaxA = (-wAC + wCG - twoQ * zeroQ * wCG + wCT) / (twoQ * wCT);
    mpq_class xyMaxB = (-wAC + wCG - twoQ * oneQ  * wCG + wCT) / (twoQ * wCT);
    mpq_class xyMinA = ( wAC + wCG - twoQ * zeroQ * wCG + wCT) / (twoQ * wCT);
    mpq_class xyMinB = ( wAC + wCG - twoQ * oneQ  * wCG + wCT) / (twoQ * wCT);
    
    mpq_class yzMaxA = (-wAT + wCT - twoQ * zeroQ * wCT + wGT) / (twoQ * wGT);
    mpq_class yzMaxB = (-wAT + wCT - twoQ * oneQ  * wCT + wGT) / (twoQ * wGT);
    mpq_class yzMinA = ( wAT + wCT - twoQ * zeroQ * wCT + wGT) / (twoQ * wGT);
    mpq_class yzMinB = ( wAT + wCT - twoQ * oneQ  * wCT + wGT) / (twoQ * wGT);
            
    Plane xz1 = Plane( Vector(xzMaxA, zeroQ, zeroQ), Vector(xzMaxB, zeroQ, oneQ), Vector(xzMaxA, oneQ, zeroQ) );
    Plane xz2 = Plane( Vector(xzMinA, zeroQ, zeroQ), Vector(xzMinB, zeroQ, oneQ), Vector(xzMinA, oneQ, zeroQ) );
    Plane xy1 = Plane( Vector(zeroQ, xyMaxA, zeroQ), Vector(oneQ, xyMaxB, zeroQ), Vector(zeroQ, xyMaxA, oneQ) );
    Plane xy2 = Plane( Vector(zeroQ, xyMinA, zeroQ), Vector(oneQ, xyMinB, zeroQ), Vector(zeroQ, xyMinA, oneQ) );
    Plane yz1 = Plane( Vector(zeroQ, zeroQ, yzMaxA), Vector(zeroQ, oneQ, yzMaxB), Vector(oneQ, zeroQ, yzMaxA) );
    Plane yz2 = Plane( Vector(zeroQ, zeroQ, yzMinA), Vector(zeroQ, oneQ, yzMinB), Vector(oneQ, zeroQ, yzMinA) );

    Plane front  = Plane( Vector(zeroQ, zeroQ, zeroQ), Vector(zeroQ,  oneQ, zeroQ), Vector( oneQ,  oneQ, zeroQ) );
    Plane back   = Plane( Vector(zeroQ, zeroQ,  oneQ), Vector(zeroQ,  oneQ,  oneQ), Vector( oneQ,  oneQ,  oneQ) );
    Plane top    = Plane( Vector(zeroQ,  oneQ, zeroQ), Vector(zeroQ,  oneQ,  oneQ), Vector( oneQ,  oneQ,  oneQ) );
    Plane bottom = Plane( Vector(zeroQ, zeroQ, zeroQ), Vector(zeroQ, zeroQ,  oneQ), Vector( oneQ, zeroQ,  oneQ) );
    Plane left   = Plane( Vector(zeroQ, zeroQ, zeroQ), Vector(zeroQ,  oneQ, zeroQ), Vector(zeroQ,  oneQ,  oneQ) );
    Plane right  = Plane( Vector( oneQ, zeroQ, zeroQ), Vector( oneQ,  oneQ, zeroQ), Vector( oneQ,  oneQ,  oneQ) );

    std::vector<Plane> planes;
    planes.push_back(front);
    planes.push_back(back);
    planes.push_back(top);
    planes.push_back(bottom);
    planes.push_back(left);
    planes.push_back(right);
    planes.push_back(xz1);
    planes.push_back(xz2);
    planes.push_back(xy1);
    planes.push_back(xy2);
    planes.push_back(yz1);
    planes.push_back(yz2);
    
    VertexFactory& vf = VertexFactory::vertexFactoryInstance();
    plain_vertex_map verticesMap;
    line_vertex_map linesMap;
    for (int i=0, n1 = (int)planes.size(); i<n1; i++)
        {
        for (int j=i+1, n2 = (int)planes.size(); j<n2; j++)
            {
            for (int k=j+1, n3 = (int)planes.size(); k<n3; k++)
                {
                if (i != j && i != k && j != k)
                    {
                    Vertex* intersectionPoint = vf.getVertex();
                    bool planesIntersect = Geometry::intersect(planes[i], planes[j], planes[k], *intersectionPoint);
                    if (planesIntersect == true && isValid(*intersectionPoint) == true)
                        {
                        // add intersection Vector to planes map
                        plain_vertex_map::iterator it = verticesMap.find(&planes[i]);
                        if (it == verticesMap.end())
                            {
                            std::vector<Vertex*> vec;
                            vec.push_back(intersectionPoint);
                            verticesMap.insert( std::make_pair(&planes[i],vec) );
                            }
                        else
                            it->second.push_back(intersectionPoint);

                        it = verticesMap.find(&planes[j]);
                        if (it == verticesMap.end())
                            {
                            std::vector<Vertex*> vec;
                            vec.push_back(intersectionPoint);
                            verticesMap.insert( std::make_pair(&planes[j],vec) );
                            }
                        else
                            it->second.push_back(intersectionPoint);

                        it = verticesMap.find(&planes[k]);
                        if (it == verticesMap.end())
                            {
                            std::vector<Vertex*> vec;
                            vec.push_back(intersectionPoint);
                            verticesMap.insert( std::make_pair(&planes[k],vec) );
                            }
                        else
                            it->second.push_back(intersectionPoint);
                            
                            
                        insertVertex(linesMap, &planes[i], &planes[j], *intersectionPoint);
                        insertVertex(linesMap, &planes[i], &planes[k], *intersectionPoint);
                        insertVertex(linesMap, &planes[j], &planes[k], *intersectionPoint);
                        }
                        
                        
                    }
                }
            }
        }

    std::cout << "Planes" << std::endl;
    for (auto p : verticesMap)
        {
        std::cout << p.first << " -- ";
        for (int i=0; i<p.second.size(); i++)
            std::cout << p.second[i] << " ";
        std::cout << std::endl;
        }
        
    std::cout << "Planes count (" << verticesMap.size() << "):" << std::endl;
    for (auto p : verticesMap)
        {
        std::cout << p.first->getStr() << " -- " << p.second.size() << std::endl;
        }

    std::cout << "Lines" << std::endl;
    for (auto p : linesMap)
        {
        std::cout << p.first.first << " " << p.first.second << " -- ";
        for (int i=0; i<p.second.size(); i++)
            std::cout << p.second[i] << " ";
        std::cout << std::endl;
        }
        
    // set up facets
    initializeFacets(verticesMap, linesMap);
    double mcVolume = monteCarloVolume(100000);
    std::cout << mcVolume << std::endl;

}

void Polyhedral::insertVertex(std::map< std::pair<Plane*,Plane*>, std::vector<Vertex*> >& linesMap, Plane* p1, Plane* p2, Vertex& v) {

    std::pair<Plane*,Plane*> key(p1, p2);
    if (p2 < p1)
        key = std::make_pair(p2, p1);
        
    std::map< std::pair<Plane*,Plane*>, std::vector<Vertex*> >::iterator it = linesMap.find(key);
    if (it == linesMap.end())
        {
        std::vector<Vertex*> vec;
        vec.push_back(&v);
        linesMap.insert( std::make_pair(key,vec) );
        }
    else
        {
        it->second.push_back(&v);
        }
}

bool Polyhedral::isLineInList(Line& x, std::vector<Line>& lines) {

    for (int i=0; i<lines.size(); i++)
        {
        if (x == lines[i])
            return true;
        }
    return false;
}

bool Polyhedral::isValid(Vector& pt) {

    mpq_class& x = pt.getX();
    mpq_class& y = pt.getY();
    mpq_class& z = pt.getZ();

    // first, check that the points are part of or inside of the unit cube
    if (x < 0 || x > 1)
        return false;
    if (y < 0 || y > 1)
        return false;
    if (z < 0 || z > 1)
        return false;

    // second, check that the points satisfy the constraints:
    // 0 \leq w_{AC} + w_{CG} (2 u_1 - 1) + w_{CT}(2 u_2 - 1) \leq 2 w_{AC}
    // 0 \leq w_{AG} - w_{CG} (2 u_1 - 1) + w_{GT}(2 u_3 - 1) \leq 2 w_{AG}
    // 0 \leq w_{AT} - w_{CT} (2 u_2 - 1) - w_{GT}(2 u_3 - 1) \leq 2 w_{AT}
    mpq_class v1 = wAC + wCG * (2 * x - 1) + wCT * (2 * y - 1);
    mpq_class v2 = wAG - wCG * (2 * x - 1) + wGT * (2 * z - 1);
    mpq_class v3 = wAT - wCT * (2 * y - 1) - wGT * (2 * z - 1);
    mpq_class v1Max = 2 * wAC;
    mpq_class v2Max = 2 * wAG;
    mpq_class v3Max = 2 * wAT;
    if ( (v1 >= 0) && (v1 <= v1Max) )
        ;
    else
        return false;
    if ( (v2 >= 0) && (v2 <= v2Max) )
        ;
    else
        return false;
    if ( (v3 >= 0) && (v3 <= v3Max) )
        ;
    else
        return false;
        
    return true;
}

double Polyhedral::monteCarloVolume(int numberReplicates) {

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    int numberInPolytope = 0;
    for (int i=0; i<numberReplicates; i++)
        {
        double u1 = rng.uniformRv();
        double u2 = rng.uniformRv();
        double u3 = rng.uniformRv();
        Vector pt(u1, u2, u3);
        //bool validVertex = checkVertex(pt, wAC, wAG, wAT, wCG, wCT, wGT);
        bool validVertex = isValid(pt);
        if (validVertex == true)
            numberInPolytope++;
        }
    return (double)numberInPolytope / numberReplicates;
}
