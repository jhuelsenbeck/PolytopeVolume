#ifndef RateMatrix_hpp
#define RateMatrix_hpp

#include <gmpxx.h>
#include <iomanip>
#include <iostream>
#include <vector>



class RateMatrix {

    public:
                            RateMatrix(void);
                            RateMatrix(const RateMatrix& m);
                           ~RateMatrix(void);
        mpq_class&          operator()(size_t r, size_t c) { return this->q[r * 4 + c]; }
        const mpq_class&    operator()(size_t r, size_t c) const { return this->q[r * 4 + c]; }
        RateMatrix&         operator=(const RateMatrix& rhs);
        mpq_class*          begin(void) { return q; }
        mpq_class*          begin(void) const { return q; }
        void                calculateStationaryFrequencies(std::vector<mpq_class>& f);
        mpq_class*          end(void) { return endBuffer; }
        mpq_class*          end(void) const { return endBuffer; }
    
    private:
        void                computeLandU(RateMatrix& aMat, RateMatrix& lMat, RateMatrix& uMat);
        void                transposeMatrix(const RateMatrix& a, RateMatrix& t);
        mpq_class*          q;
        mpq_class*          endBuffer;

    friend std::ostream& operator<<(std::ostream& os, RateMatrix& m);
};

inline std::ostream& operator<<(std::ostream& os, RateMatrix& m) {

    std::vector<mpq_class> bf(4);
    m.calculateStationaryFrequencies(bf);
    os << std::fixed << std::setprecision(10);
    mpq_class averageRate;
    for (int i=0; i<4; i++)
        {
        mpq_class sum;
        for (int j=0; j<4; j++)
            {
            if (m(i,j) > 0)
                os << " ";
            os << m(i,j) << " ";
            sum += m(i,j);
            if (i != j)
                averageRate += bf[i] * m(i,j);
            }
        os << "(" << sum << ")";
        os << std::endl;
        }
    os << "Average Rate = " << averageRate << std::endl;
    return os;
}

#endif
