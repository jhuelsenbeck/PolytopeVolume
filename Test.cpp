#include "Geometry.hpp"
#include "Line.hpp"
#include "Msg.hpp"
#include "Test.hpp"
#include "Vector.hpp"
#include "Vertex.hpp"



bool CodeTest::testAll(void) {

    std::vector<std::string> errors;
    if (testVector(errors) == false)
        Msg::warning("Vector test failed");
    if (testVertex(errors) == false)
        Msg::warning("Vertex test failed");
    if (testLine(errors) == false)
        Msg::warning("Line test failed");
        
    if (errors.size() > 0)
        std::cout << "WARNING: The following errors were found during testing:" << std::endl;
    for (int i=0, n=(int)errors.size(); i<n; i++)
        std::cout << errors[i] << std::endl;
        
    if (errors.size() > 0)
        return false;
    return true;
}

bool CodeTest::testLine(std::vector<std::string>& errors) {

    Vector v1(0, 1, 0);
    Vector v2(0, 1, 0);
    Vector v3(1, 0, 0);
    Line l1(v1, v2);
    Line l2(v1, v3);
    Vector i;
    bool inter = Geometry::intersect(l1, l2, i);
    std::cout << "intercect: " << inter << std::endl;
    if (inter == true)
        std::cout << i << std::endl;

    std::cout << l1 << std::endl;
    std::cout << l2 << std::endl;
    
    return true;
}

bool CodeTest::testVector(std::vector<std::string>& errors) {

    int numInitialErrors = (int)errors.size();

    if (numInitialErrors != errors.size())
        return false;
    return true;
}

bool CodeTest::testVertex(std::vector<std::string>& errors) {

    return true;
}
