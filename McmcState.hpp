#ifndef McmcState_hpp
#define McmcState_hpp

#include <gmpxx.h>
#include <string>
#include "MpqRateMatrix.hpp"
#include "Polyhedron.hpp"
class RandomVariable;

enum UpdateType { ToNR, ToR, Pi, ExchR, NRrates };


class McmcState {

    public:
                                McmcState(void);
        std::vector<mpq_class>& getWeights(void) { return W; }
        bool                    accept(double lnP);
        void                    adjust(void);
        bool                    getIsTimeReversible(void);
        double                  lnPriorProbability(void);
        void                    setAlphaT(double x);
        std::string             stateString(void);
        double                  update(std::string& updateType);
        void                    updateForAcceptance(void);
        void                    updateForRejection(void);
    
    private:
        bool                    checkRateMatrix(void);
        void                    normalize(std::vector<mpq_class>& vec);
        double                  updateExchangabilityRates(void);
        double                  updateNonreversibleRates(void);
        double                  updateStationaryFrequencies(void);
        double                  updateToNonReversible(void);
        double                  updateToReversible(void);
        Polyhedron              poly;
        std::vector<double>     stationaryFrequenciesAlpha;
        RandomVariable*         rng;
        MpqRateMatrix           Q;
        MpqRateMatrix           storedQ;
        std::vector<mpq_class>  W;
        UpdateType              updateType;
};

#endif
