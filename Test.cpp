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

    Line l1(1.0, 1.0);
    Line l2(-1.0, 1.0);
    std::cout << "intercect: " << l1.isIntersected(l2) << std::endl;
    
    if (l1.isIntersected(l2) == true)
        {
        Vector v = l1.intersect(l2);
        std::cout << v << std::endl;
        }

    std::cout << l1 << std::endl;
    std::cout << l2 << std::endl;
    
    return true;
}

bool CodeTest::testVector(std::vector<std::string>& errors) {

    int numInitialErrors = (int)errors.size();

    // test equality and inequality operators
    Vector v1;
    Vector v2(v1);
    if (v1 == v2)
        ;
    else
        errors.push_back("Vector: equality operator (==) fail");
    if (v1 != v2)
        errors.push_back("Vector: inequality operator (!=) fail");
        
    // test addition and subtraction operators
    Vector v3(1L, 2L, 3L);
    Vector v4(1.0f, 2.0f, 3.0f);
    Vector v5(2L, 4L, 6L);
    Vector v6 = v3 + v4;
    if (v6 != v5)
        errors.push_back("Vector: addition operator (+) fail");
    v3 += v4;
    if (v3 != v6)
        errors.push_back("Vector: addition operator (+=) fail");
    Vector v7 = v6 - v4;
    if (v7 != v4)
        errors.push_back("Vector: subtraction operator (-) fail");
    v6 -= v4;
    if (v6 != v4)
        errors.push_back("Vector: subtraction operator (-=) fail");
        
    // test set and getter functions
    v1.set(10.0, 20.0, 30.0);
    mpq_t a, b, c, d;
    mpq_inits(a, b, c, d, NULL);
    mpq_set_si(a, 10, 1);
    mpq_set_si(b, 20, 1);
    mpq_set_si(c, 30, 1);
    mpq_canonicalize(a);
    mpq_canonicalize(b);
    mpq_canonicalize(c);
    v2.set(a, b, c);
    if (v1 != v2)
        errors.push_back("Vector: set functions fail");
    mpq_class& valRef = v1.getX();
    if ((int)mpq_get_d(valRef.get_mpq_t()) != 10)
        errors.push_back("Vector: getX function fail");
    valRef = v1.getY();
    if ((int)mpq_get_d(valRef.get_mpq_t()) != 20)
        errors.push_back("Vector: getY function fail");
    valRef = v1.getZ();
    if ((int)mpq_get_d(valRef.get_mpq_t()) != 30)
        errors.push_back("Vector: getZ function fail");
        
    // test negate operator
    v1 = -v2;
    v2 = -v2;
    if (v1 != v2)
        errors.push_back("Vector: negate operator (-) fail");
    
    // test multiplication operators
    Vector v8(100L, 400L, 900L);
    v3 = v1 * v2;
    if (v3 != v8)
        errors.push_back("Vector: multiplication operator (*) fail");
    v1 *= v2;
    if (v1 != v8)
        errors.push_back("Vector: multiplication operator (*=) fail");
        
    // test scale operators
    Vector v9(10L, 20L, 30L);
    double scale = 10.0;
    mpq_set_d(a, scale);
    v1 = v6 * scale;
    v2 = v6 * a;
    if (v1 != v9)
        errors.push_back("Vector: native scale operator (*) fail");
    if (v2 != v9)
        errors.push_back("Vector: gmp scale operator (*) fail");
    v6 *= scale;
    v7 *= a;
    if (v6 != v9)
        errors.push_back("Vector: native scale operator (*=) fail");
    if (v7 != v9)
        errors.push_back("Vector: gmp scale operator (*=) fail");
        
    // test inverse scale operators
    mpq_set_si(a, 10L, 7L);
    mpq_set_si(b, 20L, 7L);
    mpq_set_si(c, 30L, 7L);
    Vector v10(a, b, c);
    scale = 7;
    mpq_set_d(d, scale);
    v1 = v6 / scale;
    v2 = v6 / d;
    if (v1 != v10)
        errors.push_back("Vector: native inverse scale operator (/) fail");
    if (v2 != v10)
        errors.push_back("Vector: gmp inverse scale operator (/) fail");
    v6 /= scale;
    v7 /= d;
    if (v6 != v10)
        errors.push_back("Vector: native inverse scale operator (/=) fail");
    if (v7 != v10)
        errors.push_back("Vector: gmp inverse scale operator (/=) fail");
        
    // test dot product
    v1.set(1, 3, -5);
    v2.set(4, -2, -1);
    int dot = v1.dot(v2);
    if (dot != 3)
        errors.push_back("Vector: dot product function fail");
    
    // test angle
    v1.set(0, -1, 0);
    v2.set(0, 1, 0);
    if ((int)v1.angle(v2) != 180)
        errors.push_back("Vector:angle function fail");
        
    // test equal function
    v1.set(0, 1, 0);
    v2.set(0, 1, 0);
    v3.set(0, 1, 2);
    if (v1.equal(v2) == false || v1.equal(v3) == true)
        errors.push_back("Vector: equal function fail");
        
    // test cross product function
    v1.set(1, 2, 3);
    v2.set(4, 5, 6);
    v3.set(-3, 6, -3);
    v4 = v1.cross(v2);
    if (v4 != v3)
        errors.push_back("Vector: cross product function fail");

    // test length
    v1.set(1, 0, 0);
    v2.set(1, 1, 0);
    if ((int)v1.length() != 1)
        errors.push_back("Vector: length function fail");
    if (fabs(v2.length() - sqrt(2.0)) > 1.0e-10)
        errors.push_back("Vector: length function fail");
        
    // test length squared
    v1.lengthSquared(a);
    mpq_set_si(b, 1, 1);
    mpq_canonicalize(b);
    if (mpq_equal(a, b) == 0)
        errors.push_back("Vector: lengthSquared function fail");
    v2.lengthSquared(a);
    mpq_set_si(b, 2, 1);
    mpq_canonicalize(b);
    if (mpq_equal(a, b) == 0)
        errors.push_back("Vector: lengthSquared function fail");
    
    // test distance function
    v1.set(1, 0, 0);
    v2.set(0, 1, 0);
    if (fabs(v1.distance(v2) - sqrt(2.0)) > 1.0e-10)
        errors.push_back("Vector: distance function fail");
        
    // test normalize function
    v1.set(10, 10, 0);
    v1.normalize();
    if (fabs(v1.length() - 1.0) > 1.0e-10)
        errors.push_back("Vector: normalize function fail");


    std::cout << "v1" << v1 << std::endl;
    std::cout << "v2" << v2 << std::endl;
    std::cout << "v3" << v3 << std::endl;
    std::cout << "v4" << v4 << std::endl;
    std::cout << "v5" << v5 << std::endl;
    std::cout << "v6" << v6 << std::endl;
    std::cout << "v7" << v7 << std::endl;
    std::cout << "v8" << v8 << std::endl;
    std::cout << "v9" << v9 << std::endl;
    std::cout << "v10" << v10 << std::endl;

    mpq_clears(a, b, c, d, NULL);

    if (numInitialErrors != errors.size())
        return false;
    return true;
}

bool CodeTest::testVertex(std::vector<std::string>& errors) {

    return true;
}
