#include <iomanip>
#include <map>
#include "Facet.hpp"
#include "Geometry.hpp"
#include "Msg.hpp"
#include "Polyhedron.hpp"
#include "RandomVariable.hpp"
#include "VertexFactory.hpp"



Polyhedron::Polyhedron(void) {
    
    zeroQ = 0;
    oneQ = 1;
    oneHalfQ = 1;
    oneHalfQ /= 2;
    twoQ = 2;
    
    randomlySample = false;

    center.setX(oneHalfQ);
    center.setY(oneHalfQ);
    center.setZ(oneHalfQ);

    // set up fixed planes of cube
    front.set( Vector(zeroQ, zeroQ, zeroQ), Vector(zeroQ,  oneQ, zeroQ), Vector( oneQ,  oneQ, zeroQ) );
    back.set( Vector(zeroQ, zeroQ,  oneQ), Vector(zeroQ,  oneQ,  oneQ), Vector( oneQ,  oneQ,  oneQ) );
    top.set( Vector(zeroQ,  oneQ, zeroQ), Vector(zeroQ,  oneQ,  oneQ), Vector( oneQ,  oneQ,  oneQ) );
    bottom.set( Vector(zeroQ, zeroQ, zeroQ), Vector(zeroQ, zeroQ,  oneQ), Vector( oneQ, zeroQ,  oneQ) );
    left.set( Vector(zeroQ, zeroQ, zeroQ), Vector(zeroQ,  oneQ, zeroQ), Vector(zeroQ,  oneQ,  oneQ) );
    right.set( Vector( oneQ, zeroQ, zeroQ), Vector( oneQ,  oneQ, zeroQ), Vector( oneQ,  oneQ,  oneQ) );

    planes.push_back(&front);
    planes.push_back(&back);
    planes.push_back(&top);
    planes.push_back(&bottom);
    planes.push_back(&left);
    planes.push_back(&right);
    planes.push_back(&xz1);
    planes.push_back(&xz2);
    planes.push_back(&xy1);
    planes.push_back(&xy2);
    planes.push_back(&yz1);
    planes.push_back(&yz2);
    
    polytopeVolume = 0.0;
}

void Polyhedron::calculateTetrahedronVolume(Vector* v1, Vector* v2, Vector* v3, mpf_class& d, mpf_class& vol) {

    mpq_class x1 = v2->getX() - v1->getX();
    mpq_class y1 = v2->getY() - v1->getY();
    mpq_class z1 = v2->getZ() - v1->getZ();
    mpq_class x2 = v3->getX() - v1->getX();
    mpq_class y2 = v3->getY() - v1->getY();
    mpq_class z2 = v3->getZ() - v1->getZ();
    
    mpq_class x3 = y1 * z2 - y2 * z1;
    mpq_class y3 = x1 * z2 - x2 * z1;
    mpq_class z3 = x1 * y2 - x2 * y1;
    
    mpq_class u = (x3 * x3) + (y3 * y3) + (z3 * z3);
    u /= 4;
    mpf_class uF = u;
    mpf_class area = sqrt(uF);
    vol = area * d / 3.0;
}

mpq_class Polyhedron::facetArea(Vertex* first, Plane* pln) {

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

Vertex* Polyhedron::findOtherVertex(Vertex* from, Vertex* v, Plane* pln) {

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

void Polyhedron::initializeFacets(void) {
    
    tetrahedra.clear();
        
    polytopeVolume = 0.0;
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
            Vertex* nextV = findOtherVertex(v->getFrom(), v, pln.first);
            v->setTo(nextV);
            nextV->setFrom(v);
            v = nextV;
            } while (v != first);
        
        // get the facet area
        mpq_class facetBaseArea = facetArea(first, pln.first);
                
        // calculate the pyramid volume
        facetBaseArea /= 3;
        mpf_class facetBaseAreaF = facetBaseArea;
        polytopeVolume += facetBaseAreaF * distance;
        
        // triangluate and sample from each tetrahedron
        if (randomlySample == true)
            sampleTetrahedra(pln.second, distance);
        }
        
    if (randomlySample == true)
        {
        mpf_class u = RandomVariable::randomVariableInstance().uniformRv();
        u *= polytopeVolume;
        mpf_class sumVol;
        for (auto tet : tetrahedra)
            {
            sumVol += tet.second;
            if (u < sumVol)
                {
                randomPoint.set(tet.first->getX(), tet.first->getY(), tet.first->getZ());
                break;
                }
            }
            
        for (auto tet : tetrahedra)
            delete tet.first;
        tetrahedra.clear();
        }
}

void Polyhedron::initializePlanes(void) {

    // make planes that will slice up the cube
    // u1 -> (-wAG + wCG - wGT + 2.0 * u3 * wGT) / (2.0 * wCG)   max in x,z
    // u1 -> ( wAG + wCG - wGT + 2.0 * u3 * wGT) / (2.0 * wCG)   min in x,z
    // u2 -> ( wAC + wCG - 2.0 * u1 * wCG + wCT) / (2.0 * wCT)   max in x,y
    // u2 -> (-wAC + wCG - 2.0 * u1 * wCG + wCT) / (2.0 * wCT)   min in x,y
    // u3 -> (-wAT + wCT - 2.0 * u2 * wCT + wGT) / (2.0 * wGT)   max in y,z
    // u3 -> ( wAT + wCT - 2.0 * u2 * wCT + wGT) / (2.0 * wGT)   min in y,z
        
    xzMaxA = (-wAG + wCG - wGT + twoQ * zeroQ * wGT) / (twoQ * wCG);
    xzMaxB = (-wAG + wCG - wGT + twoQ * oneQ  * wGT) / (twoQ * wCG);
    xzMinA = ( wAG + wCG - wGT + twoQ * zeroQ * wGT) / (twoQ * wCG);
    xzMinB = ( wAG + wCG - wGT + twoQ * oneQ  * wGT) / (twoQ * wCG);

    xyMaxA = (-wAC + wCG - twoQ * zeroQ * wCG + wCT) / (twoQ * wCT);
    xyMaxB = (-wAC + wCG - twoQ * oneQ  * wCG + wCT) / (twoQ * wCT);
    xyMinA = ( wAC + wCG - twoQ * zeroQ * wCG + wCT) / (twoQ * wCT);
    xyMinB = ( wAC + wCG - twoQ * oneQ  * wCG + wCT) / (twoQ * wCT);

    yzMaxA = (-wAT + wCT - twoQ * zeroQ * wCT + wGT) / (twoQ * wGT);
    yzMaxB = (-wAT + wCT - twoQ * oneQ  * wCT + wGT) / (twoQ * wGT);
    yzMinA = ( wAT + wCT - twoQ * zeroQ * wCT + wGT) / (twoQ * wGT);
    yzMinB = ( wAT + wCT - twoQ * oneQ  * wCT + wGT) / (twoQ * wGT);

    xzMaxA_Zero_Zero.set(xzMaxA, zeroQ, zeroQ); // Vector(xzMaxA, zeroQ, zeroQ)
    xzMinA_Zero_Zero.set(xzMinA, zeroQ, zeroQ); // Vector(xzMinA, zeroQ, zeroQ)
    xzMaxA_One_Zero.set(xzMaxA, oneQ, zeroQ);   // Vector(xzMaxA,  oneQ, zeroQ)
    xzMinA_One_Zero.set(xzMinA, oneQ, zeroQ);   // Vector(xzMinA,  oneQ, zeroQ)
    xzMaxB_Zero_One.set(xzMaxB, zeroQ, oneQ);   // Vector(xzMaxB, zeroQ,  oneQ)
    xzMinB_Zero_One.set(xzMinB, zeroQ, oneQ);   // Vector(xzMinB, zeroQ,  oneQ)
    zero_xyMaxA_Zero.set(zeroQ, xyMaxA, zeroQ); // Vector(zeroQ, xyMaxA, zeroQ)
    zero_xyMinA_Zero.set(zeroQ, xyMinA, zeroQ); // Vector(zeroQ, xyMinA, zeroQ)
    one_xyMaxB_Zero.set(oneQ, xyMaxB, zeroQ);   // Vector(oneQ, xyMaxB, zeroQ)
    one_xyMinB_Zero.set(oneQ, xyMinB, zeroQ);   // Vector(oneQ, xyMinB, zeroQ)
    zero_xyMaxA_One.set(zeroQ, xyMaxA, oneQ);   // Vector(zeroQ, xyMaxA, oneQ)
    zero_xyMinA_One.set(zeroQ, xyMinA, oneQ);   // Vector(zeroQ, xyMinA, oneQ)
    zero_Zero_yzMaxA.set(zeroQ, zeroQ, yzMaxA); // Vector(zeroQ, zeroQ, yzMaxA)
    zero_Zero_yzMinA.set(zeroQ, zeroQ, yzMinA); // Vector(zeroQ, zeroQ, yzMinA)
    zero_One_yzMaxB.set(zeroQ, oneQ, yzMaxB);   // Vector(zeroQ, oneQ, yzMaxB)
    zero_One_yzMinB.set(zeroQ, oneQ, yzMinB);   // Vector(zeroQ, oneQ, yzMinB)
    one_Zero_yzMaxA.set(oneQ, zeroQ, yzMaxA);   // Vector(oneQ, zeroQ, yzMaxA)
    one_Zero_yzMinA.set(oneQ, zeroQ, yzMinA);   // Vector(oneQ, zeroQ, yzMinA)

    // set up non-fixed planes
    xz1.set( xzMaxA_Zero_Zero, xzMaxB_Zero_One, xzMaxA_One_Zero );
    xz2.set( xzMinA_Zero_Zero, xzMinB_Zero_One, xzMinA_One_Zero );
    xy1.set( zero_xyMaxA_Zero, one_xyMaxB_Zero, zero_xyMaxA_One );
    xy2.set( zero_xyMinA_Zero, one_xyMinB_Zero, zero_xyMinA_One );
    yz1.set( zero_Zero_yzMaxA, zero_One_yzMaxB, one_Zero_yzMaxA );
    yz2.set( zero_Zero_yzMinA, zero_One_yzMinB, one_Zero_yzMinA );
    
    // empty out the tetrahedra map in preparation for finding random points in
    // in a triangulation of the polyhedron
    if (randomlySample == true)
        tetrahedra.clear();
            
    // note that this checks all 12 choose 3 combinations of planes for intersection even though
    // six pairs of the planes are parallel to one another!
    VertexFactory& vf = VertexFactory::vertexFactoryInstance();
    verticesMap.clear();
    linesMap.clear();
    for (int i=0, n1 = (int)planes.size(); i<n1; i++)
        {
        for (int j=i+1, n2 = (int)planes.size(); j<n2; j++)
            {
            for (int k=j+1, n3 = (int)planes.size(); k<n3; k++)
                {
                if (i != j && i != k && j != k)
                    {
                    Vertex* intersectionPoint = vf.getVertex();
                    bool planesIntersect = Geometry::intersect(*planes[i], *planes[j], *planes[k], *intersectionPoint);
                    if (planesIntersect == true && isValid(*intersectionPoint) == true)
                        {
                        // add intersection Vector to planes map
                        plane_vertex_map::iterator it = verticesMap.find(planes[i]);
                        if (it == verticesMap.end())
                            {
                            std::vector<Vertex*> vec;
                            vec.push_back(intersectionPoint);
                            verticesMap.insert( std::make_pair(planes[i],vec) );
                            }
                        else
                            it->second.push_back(intersectionPoint);

                        it = verticesMap.find(planes[j]);
                        if (it == verticesMap.end())
                            {
                            std::vector<Vertex*> vec;
                            vec.push_back(intersectionPoint);
                            verticesMap.insert( std::make_pair(planes[j],vec) );
                            }
                        else
                            it->second.push_back(intersectionPoint);

                        it = verticesMap.find(planes[k]);
                        if (it == verticesMap.end())
                            {
                            std::vector<Vertex*> vec;
                            vec.push_back(intersectionPoint);
                            verticesMap.insert( std::make_pair(planes[k],vec) );
                            }
                        else
                            it->second.push_back(intersectionPoint);
                            
                            
                        insertVertex(planes[i], planes[j], intersectionPoint);
                        insertVertex(planes[i], planes[k], intersectionPoint);
                        insertVertex(planes[j], planes[k], intersectionPoint);
                        }
                        
                    }
                }
            }
        }
        
    // set up facets
    initializeFacets();

    // clean up
    vf.recallAllVertices();
}

void Polyhedron::insertVertex(Plane* p1, Plane* p2, Vertex* v) {

    std::pair<Plane*,Plane*> key(p1, p2);
    if (p2 < p1)
        key = std::make_pair(p2, p1);
        
    line_vertex_map::iterator it = linesMap.find(key);
    if (it == linesMap.end())
        {
        std::vector<Vertex*> vec;
        vec.push_back(v);
        linesMap.insert( std::make_pair(key,vec) );
        }
    else
        {
        it->second.push_back(v);
        }
}

bool Polyhedron::isValid(Vector& pt) {

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

double Polyhedron::monteCarloVolume(int numberReplicates) {

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    int numberInPolytope = 0;
    for (int i=0; i<numberReplicates; i++)
        {
        double u1 = rng.uniformRv();
        double u2 = rng.uniformRv();
        double u3 = rng.uniformRv();
        Vector pt(u1, u2, u3);
        bool validVertex = isValid(pt);
        if (validVertex == true)
            numberInPolytope++;
        }
    return (double)numberInPolytope / numberReplicates;
}

void Polyhedron::sampleTetrahedra(std::vector<Vertex*>& vertices, mpf_class& d) {

    Vector pt;
    Vertex* f = vertices[0];
    Vertex* p = f->getTo();
    Vertex* v1 = f; // point common to all triangulations of this facet
    do
        {
        // form triangulation
        Vertex* v2 = p;
        Vertex* v3 = p->getTo();
        
        // calculate volume
        mpf_class tetrahedronVolume;
        calculateTetrahedronVolume(v1, v2, v3, d, tetrahedronVolume);
        
        // add a random point from tetrahedron to tetrahedra map for latter use
        Vector* newV = new Vector(pt);
        sampleTetrahedron(&center, v1, v2, v3, *newV);
        tetrahedra.insert( std::make_pair(newV,tetrahedronVolume) );
        
        // on to the next triangulation
        p = v3;
        } while (p->getTo() != f);
}

void Polyhedron::sampleTetrahedron(Vector* center, Vector* v1, Vector* v2, Vector* v3, Vector& pt) {

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    
    mpq_class s = rng.uniformRv();
    mpq_class t = rng.uniformRv();
    mpq_class u = rng.uniformRv();

    if (s + t > 1)
        {
        // cut'n fold the cube into a prism
        s = 1 - s;
        t = 1 - t;
        }
    if (t + u > 1)
        {
        // cut'n fold the prism into a tetrahedron
        mpq_class tmp = u;
        u = 1 - s - t;
        t = 1 - tmp;
        }
    else if (s + t + u > 1)
        {
        mpq_class tmp = u;
        u = s + t + u - 1;
        s = 1 - t - tmp;
        }
    mpq_class a = 1 - s - t - u; // a,s,t,u are the barycentric coordinates of the random point.
    mpq_class x = center->getX() * a + v1->getX() * s + v2->getX() * t + v3->getX() * u;
    mpq_class y = center->getY() * a + v1->getY() * s + v2->getY() * t + v3->getY() * u;
    mpq_class z = center->getZ() * a + v1->getZ() * s + v2->getZ() * t + v3->getZ() * u;
    pt.set(x, y, z); // Vector pt = v0*a + v1*s + v2*t + v3*u;
}

void Polyhedron::sampleTetrahedron(Vector* center, Vector* v1, Vector* v2, Vector* v3, Vector& pt, double shrinkageFactor) {

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    
    mpq_class factor = shrinkageFactor;
    mpq_class oneHalf = 1;
    oneHalf /= 2;
    
    Vector newV1(*v1);
    Vector newV2(*v2);
    Vector newV3(*v3);
    
    mpq_class s = rng.uniformRv();
    mpq_class t = rng.uniformRv();
    mpq_class u = rng.uniformRv();

    if (s + t > 1)
        {
        // cut'n fold the cube into a prism
        s = 1 - s;
        t = 1 - t;
        }
    if (t + u > 1)
        {
        // cut'n fold the prism into a tetrahedron
        mpq_class tmp = u;
        u = 1 - s - t;
        t = 1 - tmp;
        }
    else if (s + t + u > 1)
        {
        mpq_class tmp = u;
        u = s + t + u - 1;
        s = 1 - t - tmp;
        }
    mpq_class a = 1 - s - t - u; // a,s,t,u are the barycentric coordinates of the random point.
    mpq_class x = center->getX() * a + newV1.getX() * s + newV2.getX() * t + newV3.getX() * u;
    mpq_class y = center->getY() * a + newV1.getY() * s + newV2.getY() * t + newV3.getY() * u;
    mpq_class z = center->getZ() * a + newV1.getZ() * s + newV2.getZ() * t + newV3.getZ() * u;
    pt.set(x, y, z); // Vector pt = v0*a + v1*s + v2*t + v3*u;
}

void Polyhedron::setWeights(std::vector<mpq_class>& W) {

    // extract symbols from W
    wAC = W[0];
    wAG = W[1];
    wAT = W[2];
    wCG = W[3];
    wCT = W[4];
    wGT = W[5];
    
    initializePlanes();
}

void Polyhedron::setWeights(std::vector<mpq_class>& W, mpq_class& f) {

    // extract symbols from W
    wAC = W[0];
    wAG = W[1];
    wAT = W[2];
    wCG = W[3];
    wCT = W[4];
    wGT = W[5];
    
    initializePlanes();
}

mpf_class Polyhedron::volume(std::vector<mpq_class>& W) {

    setWeights(W);
    return polytopeVolume;
}

mpf_class Polyhedron::volume(std::vector<mpq_class>& W, Vector& pt) {

    randomlySample = true;
    setWeights(W);
    pt.set(randomPoint.getX(), randomPoint.getY(), randomPoint.getZ());
    if (isValid(pt) == false)
        Msg::error("Random point is not in polyhedron");
    randomlySample = false;
    return polytopeVolume;
}

mpf_class Polyhedron::volume(std::vector<mpq_class>& W, double fac) {

    mpq_class fQ = fac;
    setWeights(W, fQ);
    return polytopeVolume;
}
