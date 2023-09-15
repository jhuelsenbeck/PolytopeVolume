#include <map>
#include "Facet.hpp"
#include "Geometry.hpp"
#include "Msg.hpp"
#include "Polyhedral.hpp"
#include "RandomVariable.hpp"
#include "VertexFactory.hpp"


Polyhedral::Polyhedral(void) {

    zeroQ = 0;
    oneQ = 1;
    oneHalfQ = 1;
    oneHalfQ /= 2;
    twoQ = 2;
    
    computeExtrema = false;

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

void Polyhedral::initializeFacets(plain_vertex_map& verticesMap, line_vertex_map& linesMap) {

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
            Vertex* nextV = findOtherVertex(v->getFrom(), linesMap, v, pln.first);
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
    
    // initialize min and max values
    if (computeExtrema == true)
        {
        minX = 1;
        minY = 1;
        minZ = 1;
        maxX = 0;
        maxY = 0;
        maxZ = 0;
        }
    
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
                        plain_vertex_map::iterator it = verticesMap.find(planes[i]);
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
                            
                            
                        insertVertex(linesMap, planes[i], planes[j], intersectionPoint);
                        insertVertex(linesMap, planes[i], planes[k], intersectionPoint);
                        insertVertex(linesMap, planes[j], planes[k], intersectionPoint);
                        
                        // check extrema
                        if (computeExtrema == true)
                            {
                            if (intersectionPoint->getX() < minX)
                                minX = intersectionPoint->getX();
                            if (intersectionPoint->getY() < minY)
                                minY = intersectionPoint->getY();
                            if (intersectionPoint->getZ() < minZ)
                                minZ = intersectionPoint->getZ();
                            if (intersectionPoint->getX() > maxX)
                                maxX = intersectionPoint->getX();
                            if (intersectionPoint->getY() > maxY)
                                maxY = intersectionPoint->getY();
                            if (intersectionPoint->getZ() > maxZ)
                                maxZ = intersectionPoint->getZ();
                            }
                        }
                        
                        
                    }
                }
            }
        }
        
    // set up facets
    initializeFacets(verticesMap, linesMap);
    
    // clean up
    vf.recallAllVertices();
    
    if (computeExtrema == true)
        {
        diffX = maxX - minX;
        diffY = maxY - minY;
        diffZ = maxZ - minZ;
        }
}

void Polyhedral::insertVertex(line_vertex_map& linesMap, Plane* p1, Plane* p2, Vertex* v) {

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

void Polyhedral::samplePolytope(Vector& pt) {

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    mpq_class x;
    mpq_class y;
    mpq_class z;
    bool isValidPoint = false;
    do {
        x = minX + (diffX * rng.uniformRv());
        y = minY + (diffY * rng.uniformRv());
        z = minZ + (diffZ * rng.uniformRv());
        pt.setX(x);
        pt.setY(y);
        pt.setZ(z);
        if (isValid(pt) == true)
            isValidPoint = true;
        } while(isValidPoint == false);
        
#   if 1
    std::cout << "(" << minX.get_d() << " " << x.get_d() << " " << maxX.get_d() << ") ";
    std::cout << "(" << minY.get_d() << " " << y.get_d() << " " << maxY.get_d() << ") ";
    std::cout << "(" << minZ.get_d() << " " << z.get_d() << " " << maxZ.get_d() << ") " << std::endl;
#   endif
}

void Polyhedral::setWeights(std::vector<mpq_class>& W) {

    // extract symbols from W
    wAC = W[0];
    wAG = W[1];
    wAT = W[2];
    wCG = W[3];
    wCT = W[4];
    wGT = W[5];
    
    initializePlanes();
}

mpf_class Polyhedral::volume(std::vector<mpq_class>& W) {

    setWeights(W);
    return polytopeVolume;
}

mpf_class Polyhedral::volume(std::vector<mpq_class>& W, Vector& pt) {

    computeExtrema = true;
    setWeights(W);
    samplePolytope(pt);
    computeExtrema = false;
    return polytopeVolume;
}
