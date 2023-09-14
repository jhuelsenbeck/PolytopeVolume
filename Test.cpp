#include "Geometry.hpp"
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
        
    if (errors.size() > 0)
        std::cout << "WARNING: The following errors were found during testing:" << std::endl;
    for (int i=0, n=(int)errors.size(); i<n; i++)
        std::cout << errors[i] << std::endl;
        
    if (errors.size() > 0)
        return false;
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
