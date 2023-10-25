#ifndef MpqRateMatrix_hpp
#define MpqRateMatrix_hpp

#include <gmpxx.h>
#include <iomanip>
#include <iostream>
#include <vector>
class RandomVariable;



class MpqRateMatrix {

    /**
     * This class is used to hold a 4 X 4 rate matrix for a continuous-time Markov model. The
     * implementation of this class relies on GMP rationals to hold rates and to calculate
     * stationary frequencies or exchangeability parameters. This class also has functionality
     * to propose a reversible model from a non-reversible one and also to do the reverse.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (John Huelsenbeck)
     * @since 2014-11-18, version 1.0
     */

    public:
                                MpqRateMatrix(void);
                                MpqRateMatrix(const MpqRateMatrix& m);
                               ~MpqRateMatrix(void);
        mpq_class&              operator()(size_t r, size_t c) { return this->q[r * 4 + c]; }
        const mpq_class&        operator()(size_t r, size_t c) const { return this->q[r * 4 + c]; }
        MpqRateMatrix&          operator=(const MpqRateMatrix& rhs);
        void                    adjust(void);
        void                    calculateAverageRate(mpq_class& ave);
        void                    calculateStationaryFrequencies(std::vector<mpq_class>& f);
        void                    calculateWeights(std::vector<mpq_class>& wts);
        bool                    check(void);
        std::vector<mpq_class>& getExchangeabilityRates(void) { return r; }
        bool                    getIsReversible(void) { return isReversible; }
        std::vector<mpq_class>& getPi(void) { return pi; }
        void                    initializeTimeReversibleModel(RandomVariable* rng);
        void                    nonreversibilize(mpq_class& u1, mpq_class& u2, mpq_class& u3);
        void                    print(void);
        void                    reversibilize(void);
        void                    setExchangeabilityRates(void);
        void                    setIsReversible(bool tf) { isReversible = tf; }
        void                    setPi(std::vector<mpq_class>& f);
        double                  updateNonReversibleRates(RandomVariable* rng, double alpha0);
        double                  updateExchangeabilityRates(RandomVariable* rng, double alpha0);
        double                  updateStationaryFrequencies(RandomVariable* rng, double alpha0);
    
    private:
        void                    computeLandU(MpqRateMatrix& aMat, MpqRateMatrix& lMat, MpqRateMatrix& uMat);
        void                    transposeMatrix(const MpqRateMatrix& a, MpqRateMatrix& t);
        mpq_class*              q;                     // elements of the rate matrix
        mpq_class*              endBuffer;             // memory one past the end of the rate matrix array
        bool                    isReversible;          // flag indicating whether or not this rate matrix is time reversible
        std::vector<mpq_class>  pi;                    // the stationary frequencies of the rate matrix
        std::vector<mpq_class>  r;                     // the exchangeability parameters, if time reversible

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
