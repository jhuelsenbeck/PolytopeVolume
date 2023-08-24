#include <map>
#include "Geometry.hpp"
#include "Polyhedral.hpp"


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
            
    Plane xz1 = Plane( Point(xzMaxA, zeroQ, zeroQ), Point(xzMaxB, zeroQ, oneQ), Point(xzMaxA, oneQ, zeroQ) );
    Plane xz2 = Plane( Point(xzMinA, zeroQ, zeroQ), Point(xzMinB, zeroQ, oneQ), Point(xzMinA, oneQ, zeroQ) );
    Plane xy1 = Plane( Point(zeroQ, xyMaxA, zeroQ), Point(oneQ, xyMaxB, zeroQ), Point(zeroQ, xyMaxA, oneQ) );
    Plane xy2 = Plane( Point(zeroQ, xyMinA, zeroQ), Point(oneQ, xyMinB, zeroQ), Point(zeroQ, xyMinA, oneQ) );
    Plane yz1 = Plane( Point(zeroQ, zeroQ, yzMaxA), Point(zeroQ, oneQ, yzMaxB), Point(oneQ, zeroQ, yzMaxA) );
    Plane yz2 = Plane( Point(zeroQ, zeroQ, yzMinA), Point(zeroQ, oneQ, yzMinB), Point(oneQ, zeroQ, yzMinA) );

    Plane front  = Plane( Point(zeroQ, zeroQ, zeroQ), Point(zeroQ,  oneQ, zeroQ), Point( oneQ,  oneQ, zeroQ) );
    Plane back   = Plane( Point(zeroQ, zeroQ,  oneQ), Point(zeroQ,  oneQ,  oneQ), Point( oneQ,  oneQ,  oneQ) );
    Plane top    = Plane( Point(zeroQ,  oneQ, zeroQ), Point(zeroQ,  oneQ,  oneQ), Point( oneQ,  oneQ,  oneQ) );
    Plane bottom = Plane( Point(zeroQ, zeroQ, zeroQ), Point(zeroQ, zeroQ,  oneQ), Point( oneQ, zeroQ,  oneQ) );
    Plane left   = Plane( Point(zeroQ, zeroQ, zeroQ), Point(zeroQ,  oneQ, zeroQ), Point(zeroQ,  oneQ,  oneQ) );
    Plane right  = Plane( Point( oneQ, zeroQ, zeroQ), Point( oneQ,  oneQ, zeroQ), Point( oneQ,  oneQ,  oneQ) );
    

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
    
    std::map<Plane*,std::vector<Point>> planeIntersects;
    for (int i=0, n1 = (int)planes.size(); i<n1; i++)
        {
        for (int j=i+1, n2 = (int)planes.size(); j<n2; j++)
            {
            for (int k=j+1, n3 = (int)planes.size(); k<n3; k++)
                {
                if (i != j && i != k && j != k)
                    {
                    
                    Point intersectionPoint;
                    bool planesIntersect = Geometry::intersect(planes[i], planes[j], planes[k], intersectionPoint);
                    if (planesIntersect == true && isValid(intersectionPoint) == true)
                        {
                       // std::cout << intersectionPoint << std::endl;
                        
                        std::map<Plane*,std::vector<Point> >::iterator it = planeIntersects.find(&planes[i]);
                        if (it == planeIntersects.end())
                            {
                            std::vector<Point> vec;
                            vec.push_back(intersectionPoint);
                            planeIntersects.insert( std::make_pair(&planes[i],vec) );
                            }
                        else
                            it->second.push_back(intersectionPoint);

                        it = planeIntersects.find(&planes[j]);
                        if (it == planeIntersects.end())
                            {
                            std::vector<Point> vec;
                            vec.push_back(intersectionPoint);
                            planeIntersects.insert( std::make_pair(&planes[j],vec) );
                            }
                        else
                            it->second.push_back(intersectionPoint);

                        it = planeIntersects.find(&planes[k]);
                        if (it == planeIntersects.end())
                            {
                            std::vector<Point> vec;
                            vec.push_back(intersectionPoint);
                            planeIntersects.insert( std::make_pair(&planes[k],vec) );
                            }
                        else
                            it->second.push_back(intersectionPoint);

                        }
                        
                    }
                }
            }
        }

    for (auto p : planeIntersects)
        {
        std::cout << *p.first << " -- ";
        for (int i=0; i<p.second.size(); i++)
            std::cout << p.second[i].getStr()    << " ";
        std::cout << std::endl;
        }
        
#   if 0

    // check all (12 choose 3) combinations of planes, getting vertices for any that intersect
    std::map<Plane*,std::vector<Point> > verticesMap;
    std::map<Plane*,std::vector<Line> > lines;
    std::Point<Point> uniquePoints;
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
                    Point pt = planes[i].intersect(planes[j], planes[k], doPlanesIntersect);
                    if (doPlanesIntersect == true)
                        {
                        bool validVertex = checkVertex(pt, wAC, wAG, wAT, wCG, wCT, wGT);
                        if (validVertex == true)
                            {
                            uniquePoints.push_back(pt);
                            // add the vertex to all three planes
                            std::map<Plane*,std::Point<Point> >::iterator it = verticesMap.find(&planes[i]);
                            if (it != verticesMap.end())
                                {
                                it->second.push_back(pt);
                                }
                            else
                                {
                                std::Point<Point> v;
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
                                std::Point<Point> v;
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
                                std::Point<Point> v;
                                v.push_back(pt);
                                verticesMap.insert( std::make_pair(&planes[k], v) );
                                }
                                
                            Line lineIJ = planes[i].intersect(planes[j]);
                            Line lineIK = planes[i].intersect(planes[k]);
                            Line lineJK = planes[j].intersect(planes[k]);
                            std::map<Plane*,std::Point<Line> >::iterator lit = lines.find(&planes[i]);
                            if (lit != lines.end())
                                {
                                if (isLineInList(lineIJ, lit->second) == false)
                                    lit->second.push_back(lineIJ);
                                if (isLineInList(lineIK, lit->second) == false)
                                    lit->second.push_back(lineIK);
                                }
                            else
                                {
                                std::Point<Line> v;
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
                                std::Point<Line> v;
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
                                std::Point<Line> v;
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
        vertices.push_back( new Point(uniquePoints[i]) );

    // set up facets
    VertexFactory& vFactory = VertexFactory::vertexFactoryInstance();
    for (std::map<Plane*,std::Point<Point> >::iterator it = verticesMap.begin(); it != verticesMap.end(); it++)
        {
        // get the information for the plane
        Plane* p = it->first;
        std::Point<Vertex*> vertices;
        for (int i=0; i<it->second.size(); i++)
            vertices.push_back(vFactory.getVertex(it->second[i]));
        std::Point<Line> planeLines;
        std::map<Plane*,std::Point<Line> >::iterator lit = lines.find(p);
        if (lit != lines.end())
            planeLines = lit->second;
        else
            Msg::error("Couldn't find lines associated with plane");
        
        // configure vertices for each plane
        for (Line x : planeLines)
            {
            // find all of the vertices for this line
            std::Point<Vertex*> lineVertices;
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

bool Polyhedral::isValid(Point& pt) {

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
