#ifndef McmcState_hpp
#define McmcState_hpp

#include <gmpxx.h>
#include <string>
#include "Polyhedron.hpp"
#include "RateMatrix.hpp"
class RandomVariable;



class McmcState {

    public:
                                McmcState(void);
        std::vector<mpq_class>& getWeights(void) { return W; }
        bool                    accept(double lnP);
        bool                    getIsTimeReversible(void) { return isTimeReversible; }
        double                  lnPriorProbability();
        std::string             stateString(void);
        double                  update(std::string& updateType);
        void                    updateForAcceptance(void);
        void                    updateForRejection(void);
    
    private:
        void                    calculateStationaryFrequencies(void);
        bool                    checkRateMatrix(void);
        void                    initializeTimeReversibleRateMatrix(void);
        void                    normalize(std::vector<mpq_class>& vec);
        void                    setReversibleRateMatrix(void);
        void                    setNonReversibleRateMatrix(void);
        double                  updateExchangabilityRates(void);
        double                  updateNonreversibleRates(void);
        double                  updateStationaryFrequencies(void);
        double                  updateToNonReversible(void);
        double                  updateToReversible(void);
        Polyhedron              poly;
        std::vector<double>     stationaryFrequenciesAlpha;
        RandomVariable*         rng;
        std::vector<mpq_class>  W;
        bool                    isTimeReversible;
        RateMatrix              Q;
        std::vector<mpq_class>  pi;
        std::vector<mpq_class>  r;
        bool                    storedIsTimeReversible;
        RateMatrix              storedQ;
        std::vector<mpq_class>  storedPi;
        std::vector<mpq_class>  storedR;
};

#endif
