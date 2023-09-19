#ifndef Mcmc_hpp
#define Mcmc_hpp

#include <gmpxx.h>
#include "RateMatrix.hpp"
class RandomVariable;



class Mcmc {

    public:
                                Mcmc(void);
        std::vector<mpq_class>& getWeights(void) { return W; }
        double                  lnPriorProbability(std::vector<mpq_class>& f);
        void                    update(void);
    
    private:
        void                    calculateStationaryFrequencies(void);
        void                    initializeTimeReversibleRateMatrix(void);
        void                    seqQ(std::vector<mpq_class>& er, std::vector<mpq_class>& bf, RateMatrix& m);
        double                  updateExchangabilityRates(std::vector<mpq_class>& er);
        double                  updateStationaryFrequencies(std::vector<mpq_class>& bf);
        double                  updateRates(void);
        bool                    isTimeReversible;
        std::vector<double>     stationaryFrequenciesAlpha;
        RandomVariable*         rng;
        std::vector<mpq_class>  W;
        RateMatrix              Q;
        std::vector<mpq_class>  pi;
        std::vector<mpq_class>  r;
        RateMatrix              storedQ;
        std::vector<mpq_class>  storedPi;
        std::vector<mpq_class>  storedR;
};

#endif
