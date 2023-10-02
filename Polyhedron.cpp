#include <iomanip>
#include <map>
#include "Msg.hpp"
#include "Polyhedron.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "Vertex.hpp"
#include "VertexFactory.hpp"



Polyhedron::Polyhedron(void) {
    
    zeroQ     = 0;
    oneQ      = 1;
    oneHalfQ  = 1;
    oneHalfQ /= 2;
    twoQ      = 2;
    
    randomlySample = false;
    mathematicaPolyhedron = false;

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

    // planes persist and live in the planes vector
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
    
    // initialize the volume
    polytopeVolume = 0.0;
}

Polyhedron::~Polyhedron(void) {

}

void Polyhedron::calculateTetrahedronVolume(Vector* v1, Vector* v2, Vector* v3, mpq_class& vol) {

    aQ = v1->getX() - oneHalfQ;
    bQ = v2->getX() - oneHalfQ;
    cQ = v3->getX() - oneHalfQ;
    dQ = v1->getY() - oneHalfQ;
    eQ = v2->getY() - oneHalfQ;
    fQ = v3->getY() - oneHalfQ;
    gQ = v1->getZ() - oneHalfQ;
    hQ = v2->getZ() - oneHalfQ;
    iQ = v3->getZ() - oneHalfQ;
    
    // volume is 1/6 of the determinant
    vol = (aQ * eQ * iQ) - (aQ * fQ * hQ) - (bQ * dQ * iQ) + (bQ * fQ * gQ) + (cQ * dQ * hQ) - (cQ * eQ * gQ);
    vol /= 6;
    if (vol < 0)
        vol = -vol;
}

void Polyhedron::clearTetrahedraMap(void) {

    for (vector_volume_map::iterator it = tetrahedra.begin(); it != tetrahedra.end(); it++)
        delete it->first;
    tetrahedra.clear();
}

void Polyhedron::facetVolume(std::vector<Vertex*>& vertices, mpq_class& vol) {

    Vector pt;
    Vertex* v1 = vertices[0];
    Vertex* p = v1->getTo();
    do
        {
        // form triangulation, with v1 common to all triangulations of this facet
        Vertex* v2 = p;
        Vertex* v3 = p->getTo();
        
        // calculate volume
        mpq_class tetrahedronVolume;
        calculateTetrahedronVolume(v1, v2, v3, tetrahedronVolume);
        vol += tetrahedronVolume;
        
        // add a random point from tetrahedron to tetrahedra map for later use
        if (randomlySample == true)
            {
            Vector* newV = new Vector(pt);
            sampleTetrahedron(&center, v1, v2, v3, *newV);
            tetrahedra.insert( std::make_pair(newV, tetrahedronVolume) );
            }
        
        // on to the next triangulation
        p = v3;
        } while (p->getTo() != v1);
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
    
    clearTetrahedraMap();

    mpq_class vol;
    polytopeVolume = 0.0;
    // loop over all planes/facets of the polyhedron
    for (auto pln : verticesMap)
        {
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
            
        // calculate the facet volume and (potentially) randomly sample
        facetVolume(pln.second, vol);
        }
    polytopeVolume = vol;
    
    if (randomlySample == true)
        {
        mpq_class u = RandomVariable::randomVariableInstance().uniformRv();
        u *= polytopeVolume;
        mpq_class sumVol;
        for (auto tet : tetrahedra)
            {
            sumVol += tet.second;
            if (u < sumVol)
                {
                randomPoint.set(tet.first->getX(), tet.first->getY(), tet.first->getZ());
                break;
                }
            }
        clearTetrahedraMap();
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
    
    // empty out the tetrahedra map in preparation for finding random points
    // in triangulations of the polyhedron
    if (randomlySample == true)
        clearTetrahedraMap();
        
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
                    bool planesIntersect = intersect(*planes[i], *planes[j], *planes[k], *intersectionPoint);
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
        
    if (mathematicaPolyhedron == true)
        {
        std::string mathStr = mathematicaPolyhedronOutput();
        std::cout << mathStr << std::endl;
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

bool Polyhedron::intersect(Plane& plane1, Plane& plane2, Plane& plane3, Vector& intersection) {

    mpq_class& a1 = plane1.getA();
    mpq_class& b1 = plane1.getB();
    mpq_class& c1 = plane1.getC();
    mpq_class& d1 = plane1.getD();
    mpq_class& a2 = plane2.getA();
    mpq_class& b2 = plane2.getB();
    mpq_class& c2 = plane2.getC();
    mpq_class& d2 = plane2.getD();
    mpq_class& a3 = plane3.getA();
    mpq_class& b3 = plane3.getB();
    mpq_class& c3 = plane3.getC();
    mpq_class& d3 = plane3.getD();

    mpq_class detA  = a1 * (b2 * c3 - c2 * b3) + b1 * (c2 * a3 - a2 * c3) + c1 * (a2 * b3 - b2 * a3);
    if (detA == 0)
        return false;
    mpq_class detAx = -d1 * (b2 * c3 - c2 * b3) - d2 * (b3 * c1 - c3 * b1) - d3 * (b1 * c2 - c1 * b2);
    mpq_class detAy = -d1 * (c2 * a3 - a2 * c3) - d2 * (c3 * a1 - a3 * c1) - d3 * (c1 * a2 - a1 * c2);
    mpq_class detAz = -d1 * (a2 * b3 - b2 * a3) - d2 * (a3 * b1 - b3 * a1) - d3 * (a1 * b2 - b1 * a2);
    
#   if 0
    std::cout << "detA  = " << detA << std::endl;
    std::cout << "detAx = " << detAx << std::endl;
    std::cout << "detAy = " << detAy << std::endl;
    std::cout << "detAz = " << detAz << std::endl;
#   endif
    
    mpq_class x = detAx / detA;
    mpq_class y = detAy / detA;
    mpq_class z = detAz / detA;
    
    intersection.setX(x);
    intersection.setY(y);
    intersection.setZ(z);
    
    return true;
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

std::string Polyhedron::mathematicaPolyhedronOutput(void) {
 
    std::cout << std::fixed << std::setprecision(5);
    std::string str = "Graphics3D[{";
    clearTetrahedraMap();
    mpq_class vol;
    polytopeVolume = 0.0;
    int cnt = 0;
    for (auto pln : verticesMap)
        {
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
            
        // add the facet string
        str += "Polygon[{";
        v = first;
        do
            {
            str += v->getStr();
            v = v->getTo();
            if (v != first)
                str += ",";
            } while (v != first);
        str += "}]";
        
        if (cnt+1 != verticesMap.size())
            str += ",";
            
        cnt++;
        }
                
    str += "}, PlotRange->{{0,1},{0,1},{0,1}},BoxStyle->Dashed]";

    return str;
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

void Polyhedron::print(std::vector<mpq_class>& W) {

    // print the polyhedron
    mathematicaPolyhedron = true;
    setWeights(W);
    mathematicaPolyhedron = false;
    
    int numPoints = 2000;
    Vector pt;
    randomlySample = true;
    
    // uniform at random Beta(inf,0)
    std::cout << "ListPointPlot3D[{";
    for (int i=0; i<numPoints; i++)
        {
        setWeights(W);
        mpq_class& x = randomPoint.getX();
        mpq_class& y = randomPoint.getY();
        mpq_class& z = randomPoint.getZ();
        
        x -= oneHalfQ;
        y -= oneHalfQ;
        z -= oneHalfQ;
        mpq_class factor = 1.0;
        //mpq_class factor = Probability::Beta(1.0,1.0);
        x *= factor;
        y *= factor;
        z *= factor;
        x += oneHalfQ;
        y += oneHalfQ;
        z += oneHalfQ;
        pt.set(randomPoint.getX(), randomPoint.getY(), randomPoint.getZ());
        if (isValid(pt) == false)
            Msg::error("Random point is not in polyhedron");
            
        std::cout << pt.getStr();
        if (i+1 != numPoints)
            std::cout << ",";
        }
    std::cout << "}, PlotRange -> {{0, 1}, {0, 1}, {0, 1}}, BoxStyle -> Dashed, PlotStyle -> {Blue, PointSize[0.003]},BoxRatios -> 1, LabelStyle -> Opacity[0], Ticks -> None]" << std::endl;
    
    randomlySample = false;
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

void Polyhedron::setWeights(std::vector<mpq_class>& W) {

    // extract symbols from W
    this->wAC = W[0];
    this->wAG = W[1];
    this->wAT = W[2];
    this->wCG = W[3];
    this->wCT = W[4];
    this->wGT = W[5];
    
    // construct polyhedron
    initializePlanes();
}

mpq_class Polyhedron::volume(std::vector<mpq_class>& W) {

    setWeights(W);
    return polytopeVolume;
}

mpq_class Polyhedron::volume(std::vector<mpq_class>& W, Vector& pt) {

    randomlySample = true;
    setWeights(W);
    
    pt.set(randomPoint.getX(), randomPoint.getY(), randomPoint.getZ());
    if (isValid(pt) == false)
        Msg::error("Random point is not in polyhedron");
    randomlySample = false;
    
    return polytopeVolume;
}

mpq_class Polyhedron::volume(std::vector<mpq_class>& W, Vector& pt, double fac) {

    randomlySample = true;
    setWeights(W);
    
    mpq_class& x = randomPoint.getX();
    mpq_class& y = randomPoint.getY();
    mpq_class& z = randomPoint.getZ();
    x -= oneHalfQ;
    y -= oneHalfQ;
    z -= oneHalfQ;
    mpq_class factor = fac;
    x *= factor;
    y *= factor;
    z *= factor;
    x += oneHalfQ;
    y += oneHalfQ;
    z += oneHalfQ;
    pt.set(randomPoint.getX(), randomPoint.getY(), randomPoint.getZ());
    if (isValid(pt) == false)
        Msg::error("Random point is not in polyhedron");
    randomlySample = false;
    
    return polytopeVolume;
}
