#include <map>
#include "Facet.hpp"
#include "Geometry.hpp"
#include "Msg.hpp"
#include "Polyhedral.hpp"
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
    
    std::map<Plane*,std::vector<Vector>> verticesMap;
    std::map<Plane*,std::vector<Line> > lines;
    for (int i=0, n1 = (int)planes.size(); i<n1; i++)
        {
        for (int j=i+1, n2 = (int)planes.size(); j<n2; j++)
            {
            for (int k=j+1, n3 = (int)planes.size(); k<n3; k++)
                {
                if (i != j && i != k && j != k)
                    {
                    
                    Vector intersectionPoint;
                    bool planesIntersect = Geometry::intersect(planes[i], planes[j], planes[k], intersectionPoint);
                    if (planesIntersect == true && isValid(intersectionPoint) == true)
                        {
                        // add intersection Vector to planes map
                        std::map<Plane*,std::vector<Vector> >::iterator it = verticesMap.find(&planes[i]);
                        if (it == verticesMap.end())
                            {
                            std::vector<Vector> vec;
                            vec.push_back(intersectionPoint);
                            verticesMap.insert( std::make_pair(&planes[i],vec) );
                            }
                        else
                            it->second.push_back(intersectionPoint);

                        it = verticesMap.find(&planes[j]);
                        if (it == verticesMap.end())
                            {
                            std::vector<Vector> vec;
                            vec.push_back(intersectionPoint);
                            verticesMap.insert( std::make_pair(&planes[j],vec) );
                            }
                        else
                            it->second.push_back(intersectionPoint);

                        it = verticesMap.find(&planes[k]);
                        if (it == verticesMap.end())
                            {
                            std::vector<Vector> vec;
                            vec.push_back(intersectionPoint);
                            verticesMap.insert( std::make_pair(&planes[k],vec) );
                            }
                        else
                            it->second.push_back(intersectionPoint);

                        // add line of intersection to lines map
                        Line lineIJ;
                        Line lineIK;
                        Line lineJK;
                        bool intersectsIJ = Geometry::intersect(planes[i], planes[j], lineIJ);
                        bool intersectsIK = Geometry::intersect(planes[i], planes[k], lineIK);
                        bool intersectsJK = Geometry::intersect(planes[j], planes[k], lineJK);
                        if (intersectsIJ == false || intersectsIK == false || intersectsJK == false)
                            Msg::error("Planes should intersect forming line");
                        
                        std::map<Plane*,std::vector<Line> >::iterator lit = lines.find(&planes[i]);
                        if (lit != lines.end())
                            {
                            if (isLineInList(lineIJ, lit->second) == false)
                                lit->second.push_back(lineIJ);
                            if (isLineInList(lineIK, lit->second) == false)
                                lit->second.push_back(lineIK);
                            }
                        else
                            {
                            std::vector<Line> v;
                            v.push_back(lineIJ);
                            v.push_back(lineIK);
                            lines.insert( std::make_pair(&planes[i], v) );
                            }
                        lit = lines.find(&planes[j]);
                        if (lit != lines.end())
                            {
                            if (isLineInList(lineIJ, lit->second) == false)
                                lit->second.push_back(lineIJ);
                            if (isLineInList(lineJK, lit->second) == false)
                                lit->second.push_back(lineJK);
                            }
                        else
                            {
                            std::vector<Line> v;
                            v.push_back(lineIJ);
                            v.push_back(lineJK);
                            lines.insert( std::make_pair(&planes[j], v) );
                            }
                        lit = lines.find(&planes[k]);
                        if (lit != lines.end())
                            {
                            if (isLineInList(lineIK, lit->second) == false)
                                lit->second.push_back(lineIK);
                            if (isLineInList(lineJK, lit->second) == false)
                                lit->second.push_back(lineJK);
                            }
                        else
                            {
                            std::vector<Line> v;
                            v.push_back(lineIK);
                            v.push_back(lineJK);
                            lines.insert( std::make_pair(&planes[k], v) );
                            }


                        }
                        
                    }
                }
            }
        }

    std::cout << "Planes" << std::endl;
    for (auto p : verticesMap)
        {
        std::cout << *p.first << " -- ";
        for (int i=0; i<p.second.size(); i++)
            std::cout << p.second[i].getStr()    << " ";
        std::cout << std::endl;
        }
    std::cout << "Lines" << std::endl;
    for (auto l : lines)
        {
        std::cout << *l.first << " -- ";
        for (int i=0; i<l.second.size(); i++)
            std::cout << l.second[i]    << " ";
        std::cout << std::endl;
        }
        
    // set up facets
    VertexFactory& vFactory = VertexFactory::vertexFactoryInstance();
    for (std::map<Plane*,std::vector<Vector> >::iterator it = verticesMap.begin(); it != verticesMap.end(); it++)
        {
        // get the information for the plane
        Plane* p = it->first;
        std::vector<Vertex*> vertices;
        for (size_t i=0, n=it->second.size(); i<n; i++)
            vertices.push_back(vFactory.getVertex(it->second[i]));
        std::vector<Line> planeLines;
        std::map<Plane*,std::vector<Line> >::iterator lit = lines.find(p);
        if (lit != lines.end())
            planeLines = lit->second;
        else
            Msg::error("Couldn't find lines associated with plane");
        
        // configure vertices for each plane
        for (Line x : planeLines)
            {
            // find all of the vertices for this line
            std::vector<Vertex*> lineVertices;
            for (int i=0; i<vertices.size(); i++)
                {
                mpq_class d = x.distanceToPoint(*vertices[i]);
                if ( d == 0 )
                    lineVertices.push_back(vertices[i]);
                }
            if (lineVertices.size() != 2)
                Msg::error("There only be two vertices on this line, but we found " + std::to_string(lineVertices.size()));
                
            // link up the two vertices that were found
            if (lineVertices[0]->getFrom() == NULL && lineVertices[1]->getTo() == NULL)
                {
                lineVertices[0]->setFrom(lineVertices[1]);
                lineVertices[1]->setTo(lineVertices[0]);
                }
            else if (lineVertices[0]->getTo() == NULL && lineVertices[1]->getFrom() == NULL)
                {
                lineVertices[0]->setTo(lineVertices[1]);
                lineVertices[1]->setFrom(lineVertices[0]);
                }
            else
                Msg::error("Cannot find empty slots to hook up vertices");
            }
            
        // check zero values on vertices
        /*for (int i=0; i<vertices.size(); i++)
            {
            if (PolytopeMath::isZero(vertices[i]->x, 10e-10) == true)
                vertices[i]->x = 0.0;
            if (PolytopeMath::isZero(vertices[i]->y, 10e-10) == true)
                vertices[i]->y = 0.0;
            if (PolytopeMath::isZero(vertices[i]->z, 10e-10) == true)
                vertices[i]->z = 0.0;
            }*/ // TEMP: JPH 8/21/23 Is this necessary?
            
        // get polytope bounds
        minX = vertices[0]->x;
        minY = vertices[0]->y;
        minZ = vertices[0]->z;
        maxX = minX;
        maxY = minY;
        maxZ = minZ;
        for (int i=1; i<vertices.size(); i++)
            {
            Vertex v = *vertices[i];
            if (v.x < minX)
                minX = v.x;
            if (v.x > maxX)
                maxX = v.x;
            if (v.y < minY)
                minY = v.y;
            if (v.y > maxY)
                maxY = v.y;
            if (v.z < minZ)
                minZ = v.z;
            if (v.z > maxZ)
                maxZ = v.z;
            }
        
        // allocate the facet and add the vertices
        Facet* f = new Facet(*p, vertices);
        facets.push_back(f);
        }
        
        
        
        
        
        
        
        
#   if 0

    // check all (12 choose 3) combinations of planes, getting vertices for any that intersect
    std::map<Plane*,std::vector<Vector> > verticesMap;
    std::map<Plane*,std::vector<Line> > lines;
    std::Vector<Vector> uniquePoints;
    int n = 0;
    for (int i=0; i<planes.size(); i++)
        {
        for (int j=i+1; j<planes.size(); j++)
            {
            for (int k=j+1; k<planes.size(); k++)
                {
                if (i != j && i != k && j != k)
                    {
                    n++;
                   bool doPlanesIntersect = false;
                    Vector pt = planes[i].intersect(planes[j], planes[k], doPlanesIntersect);
                    if (doPlanesIntersect == true)
                        {
                        bool validVertex = checkVertex(pt, wAC, wAG, wAT, wCG, wCT, wGT);
                        if (validVertex == true)
                            {
                            uniquePoints.push_back(pt);
                            // add the vertex to all three planes
                            std::map<Plane*,std::Vector<Vector> >::iterator it = verticesMap.find(&planes[i]);
                            if (it != verticesMap.end())
                                {
                                it->second.push_back(pt);
                                }
                            else
                                {
                                std::Vector<Vector> v;
                                v.push_back(pt);
                                verticesMap.insert( std::make_pair(&planes[i], v) );
                                }
                            it = verticesMap.find(&planes[j]);
                            if (it != verticesMap.end())
                                {
                                it->second.push_back(pt);
                                }
                            else
                                {
                                std::Vector<Vector> v;
                                v.push_back(pt);
                                verticesMap.insert( std::make_pair(&planes[j], v) );
                                }
                            it = verticesMap.find(&planes[k]);
                            if (it != verticesMap.end())
                                {
                                it->second.push_back(pt);
                                }
                            else
                                {
                                std::Vector<Vector> v;
                                v.push_back(pt);
                                verticesMap.insert( std::make_pair(&planes[k], v) );
                                }
                                
                            Line lineIJ = planes[i].intersect(planes[j]);
                            Line lineIK = planes[i].intersect(planes[k]);
                            Line lineJK = planes[j].intersect(planes[k]);
                            std::map<Plane*,std::Vector<Line> >::iterator lit = lines.find(&planes[i]);
                            if (lit != lines.end())
                                {
                                if (isLineInList(lineIJ, lit->second) == false)
                                    lit->second.push_back(lineIJ);
                                if (isLineInList(lineIK, lit->second) == false)
                                    lit->second.push_back(lineIK);
                                }
                            else
                                {
                                std::Vector<Line> v;
                                v.push_back(lineIJ);
                                v.push_back(lineIK);
                                lines.insert( std::make_pair(&planes[i], v) );
                                }
                            lit = lines.find(&planes[j]);
                            if (lit != lines.end())
                                {
                                if (isLineInList(lineIJ, lit->second) == false)
                                    lit->second.push_back(lineIJ);
                                if (isLineInList(lineJK, lit->second) == false)
                                    lit->second.push_back(lineJK);
                                }
                            else
                                {
                                std::Vector<Line> v;
                                v.push_back(lineIJ);
                                v.push_back(lineJK);
                                lines.insert( std::make_pair(&planes[j], v) );
                                }
                            lit = lines.find(&planes[k]);
                            if (lit != lines.end())
                                {
                                if (isLineInList(lineIK, lit->second) == false)
                                    lit->second.push_back(lineIK);
                                if (isLineInList(lineJK, lit->second) == false)
                                    lit->second.push_back(lineJK);
                                }
                            else
                                {
                                std::Vector<Line> v;
                                v.push_back(lineIK);
                                v.push_back(lineJK);
                                lines.insert( std::make_pair(&planes[k], v) );
                                }
                            }
                        std::cout << n << " " << std::setw(3) << i << std::setw(3) << j << std::setw(3) << k << " -- " << pt << " " << validVertex << std::endl;
                        }
                    }
                }
            }
        }


    // add the vertices to the Polyhedral
    for (int i=0; i<uniquePoints.size(); i++)
        vertices.push_back( new Vector(uniquePoints[i]) );

    // set up facets
    VertexFactory& vFactory = VertexFactory::vertexFactoryInstance();
    for (std::map<Plane*,std::Vector<Vector> >::iterator it = verticesMap.begin(); it != verticesMap.end(); it++)
        {
        // get the information for the plane
        Plane* p = it->first;
        std::Vector<Vertex*> vertices;
        for (int i=0; i<it->second.size(); i++)
            vertices.push_back(vFactory.getVertex(it->second[i]));
        std::Vector<Line> planeLines;
        std::map<Plane*,std::Vector<Line> >::iterator lit = lines.find(p);
        if (lit != lines.end())
            planeLines = lit->second;
        else
            Msg::error("Couldn't find lines associated with plane");
        
        // configure vertices for each plane
        for (Line x : planeLines)
            {
            // find all of the vertices for this line
            std::Vector<Vertex*> lineVertices;
            for (int i=0; i<vertices.size(); i++)
                {
                double d = x.distanceToPoint(*vertices[i]);
                if ( d < 10e-8 )
                    lineVertices.push_back(vertices[i]);
                }
            if (lineVertices.size() != 2)
                Msg::error("There only be two vertices on this line, but we found " + std::to_string(lineVertices.size()));
                
            // link up the two vertices that were found
            if (lineVertices[0]->getFrom() == NULL && lineVertices[1]->getTo() == NULL)
                {
                lineVertices[0]->setFrom(lineVertices[1]);
                lineVertices[1]->setTo(lineVertices[0]);
                }
            else if (lineVertices[0]->getTo() == NULL && lineVertices[1]->getFrom() == NULL)
                {
                lineVertices[0]->setTo(lineVertices[1]);
                lineVertices[1]->setFrom(lineVertices[0]);
                }
            else
                Msg::error("Cannot find empty slots to hook up vertices");
            }
            
        // check zero values on vertices
        /*for (int i=0; i<vertices.size(); i++)
            {
            if (PolytopeMath::isZero(vertices[i]->x, 10e-10) == true)
                vertices[i]->x = 0.0;
            if (PolytopeMath::isZero(vertices[i]->y, 10e-10) == true)
                vertices[i]->y = 0.0;
            if (PolytopeMath::isZero(vertices[i]->z, 10e-10) == true)
                vertices[i]->z = 0.0;
            }*/ // TEMP: JPH 8/21/23 Is this necessary?
            
        // get polytope bounds
        minX = vertices[0]->x;
        minY = vertices[0]->y;
        minZ = vertices[0]->z;
        maxX = minX;
        maxY = minY;
        maxZ = minZ;
        for (int i=1; i<vertices.size(); i++)
            {
            Vertex v = *vertices[i];
            if (v.x < minX)
                minX = v.x;
            if (v.x > maxX)
                maxX = v.x;
            if (v.y < minY)
                minY = v.y;
            if (v.y > maxY)
                maxY = v.y;
            if (v.z < minZ)
                minZ = v.z;
            if (v.z > maxZ)
                maxZ = v.z;
            }
        
        // allocate the facet and add the vertices
        Facet* f = new Facet(*p, vertices);
        facets.push_back(f);
        }
#   endif
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
