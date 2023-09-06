#ifndef Test_hpp
#define Test_hpp

#include <string>
#include <vector>



namespace CodeTest {

    bool testAll(void);
    bool testLine(std::vector<std::string>& errors);
    bool testVector(std::vector<std::string>& errors);
    bool testVertex(std::vector<std::string>& errors);
};

#endif
