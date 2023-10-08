#ifndef MpqRateMatrix_hpp
#define MpqRateMatrix_hpp

#include <gmpxx.h>
#include <iomanip>
#include <iostream>
#include <vector>
class RandomVariable;



class MpqRateMatrix {

    public:
                                MpqRateMatrix(void);
                                MpqRateMatrix(const MpqRateMatrix& m);
                               ~MpqRateMatrix(void);
        mpq_class&              operator()(size_t r, size_t c) { return this->q[r * 4 + c]; }
        const mpq_class&        operator()(size_t r, size_t c) const { return this->q[r * 4 + c]; }
        MpqRateMatrix&          operator=(const MpqRateMatrix& rhs);
        void                    calculateStationaryFrequencies(std::vector<mpq_class>& f);
        std::vector<mpq_class>& getExchangeabilityRates(void) { return r; }
        bool                    getIsReversible(void) { return isReversible; }
        std::vector<mpq_class>& getPi(void) { return pi; }
        void                    initializeTimeReversibleModel(RandomVariable* rng);
        void                    nonreversibilize(mpq_class& u1, mpq_class& u2, mpq_class& u3);
        void                    print(void);
        void                    reversibilize(void);
        void                    setExchangeabilityRates(void);
    
    private:
        void                    computeLandU(MpqRateMatrix& aMat, MpqRateMatrix& lMat, MpqRateMatrix& uMat);
        void                    transposeMatrix(const MpqRateMatrix& a, MpqRateMatrix& t);
        mpq_class*              q;
        mpq_class*              endBuffer;
        bool                    isReversible;
        std::vector<mpq_class>  pi;
        std::vector<mpq_class>  r;

    friend std::ostream& operator<<(std::ostream& os, MpqRateMatrix& m);
};

inline std::ostream& operator<<(std::ostream& os, MpqRateMatrix& m) {

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
