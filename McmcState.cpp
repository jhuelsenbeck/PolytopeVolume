#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include "McmcState.hpp"
#include "Msg.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"

#define MIN_FREQ    10e-4
#define A           0
#define C           1
#define G           2
#define T           3



McmcState::McmcState(void) {
    
    rng = &RandomVariable::getInstance();
    W.resize(6);
    stationaryFrequenciesAlpha.resize(4);
    for (int i=0; i<4; i++)
        stationaryFrequenciesAlpha[i] = 1.0;
    Q.initializeTimeReversibleModel(rng);
    storedQ = Q;
}

bool McmcState::accept(double lnP) {

    if (log(rng->uniformRv()) < lnP)
        return true;
    return false;
}

void McmcState::adjust(void) {

    return;
    Q.adjust();
    storedQ = Q;
}

bool McmcState::checkRateMatrix(void) {

    return Q.check();
}

bool McmcState::getIsTimeReversible(void) {

    return Q.getIsReversible();
}

double McmcState::lnPriorProbability() {

    std::vector<mpq_class>& pi = Q.getPi();
    std::vector<double> bf(4);
    for (int i=0; i<4; i++)
        bf[i] = pi[i].get_d();
    double lnProb = Probability::Dirichlet::lnPdf(stationaryFrequenciesAlpha, bf);
    
    lnProb += 3.0 * (log(bf[0]) + log(bf[1]) + log(bf[2]) + log(bf[3]));
    if (Q.getIsReversible() == true)
        {
        lnProb += log(120.0);
        lnProb += log(4.0) + 0.5 * log(3.0);
        }
    else
        {
        lnProb += log(40320.0);
        lnProb -= 0.5 * log(6.0);
        }
    return lnProb;
}

void McmcState::normalize(std::vector<mpq_class>& vec) {

    mpq_class sum;
    for (int i=0, n=(int)vec.size(); i<n; i++)
        sum += vec[i];
    for (int i=0, n=(int)vec.size(); i<n; i++)
        vec[i] /= sum;
}

void McmcState::setAlphaT(double x) {

    poly.setAlphaT(x);
}

std::string McmcState::stateString(void) {

    if (Q.getIsReversible() == true)
        {
        std::vector<mpq_class>& pi = Q.getPi();
        std::vector<mpq_class>& r = Q.getExchangeabilityRates();
        std::stringstream ss;
        ss << std::fixed << std::setprecision(6);
        ss << "TR: ";
        for (int i=0; i<4; i++)
            ss << pi[i].get_d() << " ";
        ss << "-- ";
        for (int i=0; i<6; i++)
            ss << r[i].get_d() << " ";
        return ss.str();
        }
    else
        {
        std::vector<mpq_class>& pi = Q.getPi();
        std::stringstream ss;
        ss << std::fixed << std::setprecision(6);
        ss << "NR: ";
        for (int i=0; i<4; i++)
            ss << pi[i].get_d() << " ";
        return ss.str();
        }
}

double McmcState::update(std::string& updateType) {

    double lnProb = 0.0;
    if (Q.getIsReversible() == true)
        {
        // update GTR model
        double u = rng->uniformRv();
        if (u < 0.55)
            {
            lnProb = updateStationaryFrequencies();
            updateType = "BF  ";
            }
        else if (u >= 0.55 && u < 0.90)
            {
            lnProb = updateExchangabilityRates();
            updateType = "ER  ";
            }
        else
            {
            lnProb = updateToNonReversible();
            updateType = "ToNR";
            }
        }
    else
        {
        // update non-reversible model
        double u = rng->uniformRv();
        if (u < 0.90)
            {
            lnProb = updateNonreversibleRates();
            updateType = "Rts ";
            }
        else
            {
            lnProb = updateToReversible();
            updateType = "ToTR";
            }
        }
    if (checkRateMatrix() == false)
        {
        std::vector<mpq_class>& pi = Q.getPi();
        std::cout << "Problem rate matrix in " << updateType << std::endl;
        std::cout << Q << std::endl;
        for (int i=0; i<4; i++)
            std::cout << pi[i] << " ";
        std::cout << std::endl;
        mpq_class sum;
        for (int i=0; i<4; i++)
            {
            for (int j=0; j<4; j++)
                {
                if (i != j)
                    sum += pi[i] * Q(i,j);
                }
            }
        std::cout << "Average rate = " << sum.get_d() << std::endl;
        Msg::error("Problem with rate matrix");
        }
    return lnProb;
}

double McmcState::updateExchangabilityRates(void) {

    if (Q.getIsReversible() == false)
        Msg::error("Can only update the exchangability rates for time reversible models");
        
    // register the update type
    updateType = ExchR;
    
    // store old values
    storedQ = Q;

    // update and return log probability of Hastings ratio
    return Q.updateExchangeabilityRates(rng, 100.0);
}

void McmcState::updateForAcceptance(void) {

    storedQ = Q;
}

void McmcState::updateForRejection(void) {
        
    Q = storedQ;
}

double McmcState::updateNonreversibleRates(void) {

    if (Q.getIsReversible() == true)
        Msg::error("Can only update the 12 rates (directly) for non-reversible models");

    // register the update type
    updateType = NRrates;

    // store old values
    storedQ = Q;

    // update and return log probability of Hastings ratio
    return Q.updateNonReversibleRates(rng, 500.0);
}

double McmcState::updateToNonReversible(void) {

    if (Q.getIsReversible() == false)
        Msg::error("Can only move to a non-reversible model from a reversible one");

    // register the update type
    updateType = ToNR;

    // store old values
    storedQ = Q;

    // calculate the Jacobian
    std::vector<mpq_class>& pi = Q.getPi();
    double piC = pi[C].get_d();
    double piG = pi[G].get_d();
    double wCG = piC * storedQ(C,G).get_d();
    double wCT = piC * storedQ(C,T).get_d();
    double wGT = piG * storedQ(G,T).get_d();
    double lnJacobian = (log(64.0) + log(wCG) + log(wCT) + log(wGT));
    lnJacobian -= log(piC) + 2.0 * log(piG);

    // polyhedron density for forward move
    Q.calculateWeights(W);
    Vector randomPoint;
    double lnRv = -poly.lnProbabilityForward(W, randomPoint);
    
    // update rates (off diagonal components)
    mpq_class u1 = randomPoint.getX();
    mpq_class u2 = randomPoint.getY();
    mpq_class u3 = randomPoint.getZ();
#   if 1
    Q.nonreversibilize(u1, u2, u3);
#   else
    Q.setIsReversible(false);
    Q(A,C) = storedQ(A,C) + (pi[C]/pi[A]) * storedQ(C,G) * (2 * u1 - 1) + (pi[C]/pi[A]) * storedQ(C,T) * (2 * u2 - 1);
    Q(A,G) = storedQ(A,G) - (pi[C]/pi[A]) * storedQ(C,G) * (2 * u1 - 1) + (pi[G]/pi[A]) * storedQ(G,T) * (2 * u3 - 1);
    Q(A,T) = storedQ(A,T) - (pi[C]/pi[A]) * storedQ(C,T) * (2 * u2 - 1) - (pi[G]/pi[A]) * storedQ(G,T) * (2 * u3 - 1);
    Q(C,G) = 2 * storedQ(C,G) * u1;
    Q(C,T) = 2 * storedQ(C,T) * u2;
    Q(G,T) = 2 * storedQ(G,T) * u3;
    Q(C,A) = (pi[A]/pi[C]) * storedQ(A,C) - storedQ(C,G) * (2 * u1 - 1) - storedQ(C,T) * (2 * u2 - 1);
    Q(G,A) = (pi[A]/pi[G]) * storedQ(A,G) + (pi[C]/pi[G]) * storedQ(C,G) * (2 * u1 - 1) - storedQ(G,T) * (2 * u3 - 1);
    Q(T,A) = (pi[A]/pi[T]) * storedQ(A,T) + (pi[C]/pi[T]) * storedQ(C,T) * (2 * u2 - 1) + (pi[G]/pi[T]) * storedQ(G,T) * (2 * u3 - 1);
    Q(G,C) = 2 * (pi[C]/pi[G]) * storedQ(C,G) * (1 - u1);
    Q(T,C) = 2 * (pi[C]/pi[T]) * storedQ(C,T) * (1 - u2);
    Q(T,G) = 2 * (pi[G]/pi[T]) * storedQ(G,T) * (1 - u3);
    
    // update the diagonal components of the rate matrix
    mpq_class sum;
    for (int i=0; i<4; i++)
        {
        sum = 0;
        for (int j=0; j<4; j++)
            {
            if (i != j)
                sum += Q(i,j);
            }
        Q(i,i) = -sum;
        }
#   endif
    
    return lnRv + lnJacobian;
}

double McmcState::updateToReversible(void) {

    if (Q.getIsReversible() == true)
        Msg::error("Can only move to a time reversible model from a non-reversible one");

    // register the update type
    updateType = ToR;

    // store old values
    storedQ = Q;

    // average rates with same stationary frequencies
#   if 1
    Q.reversibilize();
#   else
    Q.setIsReversible(true);
    std::vector<mpq_class>& pi = Q.getPi();
    mpq_class averageRate;
    mpq_class sum;
    for (int i=0; i<4; i++)
        {
        sum = 0;
        for (int j=0; j<4; j++)
            {
            if (i != j)
                {
                Q(i,j) = (pi[i] * storedQ(i,j) + pi[j] * storedQ(j,i)) / (2 * pi[i]);
                sum += Q(i,j);
                averageRate += pi[i] * Q(i,j);
                }
            }
        Q(i,i) = -sum;
        }
        
    Q.setExchangeabilityRates();
    
    // make certain average rate is one
    if (averageRate != 1)
        {
        mpq_class factor = 1 / averageRate;
        std::cout << "Warning: the average rate should be one in updateToReversible" << std::endl;
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                Q(i,j) *= factor;
        }
#   endif

    // jacobian
    std::vector<mpq_class>& pi = Q.getPi();
    double piC = pi[C].get_d();
    double piG = pi[G].get_d();
    double piT = pi[T].get_d();
    double wCG = piC * storedQ(C,G).get_d();
    double wGC = piG * storedQ(G,C).get_d();
    double wCT = piC * storedQ(C,T).get_d();
    double wTC = piT * storedQ(T,C).get_d();
    double wGT = piG * storedQ(G,T).get_d();
    double wTG = piT * storedQ(T,G).get_d();
    //double lnJacobian = (piC * piG * piG) / (8.0 * (wCG + wGC) * (wCT + wTC) * (wGT + wTG));
    double lnJacobian = log(piC) + 2.0 * log(piG);
    lnJacobian -= (log(8.0) + log(wCG + wGC) + log(wCT + wTC) + log(wGT + wTG));
    
    // polyhedron parameters
    Q.calculateWeights(W);
    
    // figure out where in the (u1,u2,u3) space the non-reversible model lives
    mpq_class u1 = storedQ(C,G) / (2 * Q(C,G));
    mpq_class u2 = storedQ(C,T) / (2 * Q(C,T));
    mpq_class u3 = storedQ(G,T) / (2 * Q(G,T));
    Vector pt(u1, u2, u3);
    double lnRv = poly.lnProbabilityReverse(W, pt);

    return lnRv + lnJacobian;
}

double McmcState::updateStationaryFrequencies(void) {

    if (Q.getIsReversible() == false)
        Msg::error("Can only update the stationary frequencies (directly) for time reversible models");

    // register the update type
    updateType = Pi;

    // store old values
    storedQ = Q;

    // update and return log probability of Hastings ratio
    return Q.updateStationaryFrequencies(rng, 100.0);
}
