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
    
    rng = &RandomVariable::randomVariableInstance();
    r.resize(6);
    storedR.resize(6);
    pi.resize(4);
    storedPi.resize(4);
    W.resize(6);
    isTimeReversible = true;
    stationaryFrequenciesAlpha.resize(4);
    for (int i=0; i<4; i++)
        stationaryFrequenciesAlpha[i] = 1.0;
    initializeTimeReversibleRateMatrix();
    storedQ = Q;
}

bool McmcState::accept(double lnP) {

    if (log(rng->uniformRv()) < lnP)
        return true;
    return false;
}

void McmcState::calculateStationaryFrequencies(void) {

    Q.calculateStationaryFrequencies(pi);
}

bool McmcState::checkRateMatrix(void) {

    mpq_class averageRate;
    for (int i=0; i<4; i++)
        {
        mpq_class sum;
        for (int j=0; j<4; j++)
            {
            if (i != j)
                sum += pi[i] * Q(i,j);
            }
        if (pi[i] * Q(i,i) != -sum)
            return false;
        averageRate += sum;
        }
    if (averageRate != 1)
        return false;
    return true;
}

void McmcState::initializeTimeReversibleRateMatrix(void) {

    mpq_class oneHalf;
    oneHalf = 1;
    oneHalf /= 2;
    mpq_class one = 1;

    // randomly initialize parameters of the GTR model
    
    // first, choose the stationary frequency
    std::vector<double> alpha4(4, 1.0);
    std::vector<double> bf(4);
    Probability::Dirichlet::rv(rng, alpha4, bf);
    for (int i=0; i<4; i++)
        {
        if (bf[i] < MIN_FREQ)
            bf[i] = MIN_FREQ;
        }
    double sumD = 0.0;
    for (int i=0; i<4; i++)
        sumD += bf[i];
    for (int i=0; i<4; i++)
        bf[i] /= sumD;
    pi[0] = bf[0];
    pi[1] = bf[1];
    pi[2] = bf[2];
    mpq_class s = pi[0] + pi[1] + pi[2];
    pi[3] = one - s;
    for (int i=0; i<4; i++)
        std::cout << "pi[" << i << "] = " << pi[i] << std::endl;
    
    // also, choose the exchangability rates
    std::vector<double> alpha6(6, 1.0);
    std::vector<double> er(6);
    Probability::Dirichlet::rv(rng, alpha6, er);
    for (int i=0; i<6; i++)
        {
        if (er[i] < MIN_FREQ)
            er[i] = MIN_FREQ;
        }
    sumD = 0.0;
    for (int i=0; i<6; i++)
        sumD += er[i];
    for (int i=0; i<6; i++)
        er[i] /= sumD;
    r[0] = er[0];
    r[1] = er[1];
    r[2] = er[2];
    r[3] = er[3];
    r[4] = er[4];
    s = r[0] + r[1] + r[2] + r[3] + r[4];
    r[5] = one - s;

    // set the rate matrix
    for (int i=0, k=0; i<4; i++)
        {
        for (int j=i+1; j<4; j++)
            {
            Q(i,j) = r[k] * pi[j];
            Q(j,i) = r[k] * pi[i];
            k++;
            }
        }
    mpq_class averageRate = 0;
    for (int i=0; i<4; i++)
        {
        mpq_class sum = 0;
        for (int j=0; j<4; j++)
            {
            if (i != j)
                sum += Q(i,j);
            }
        Q(i,i) = -sum;
        averageRate += pi[i] * sum;
        }
    mpq_class factor = 1 / averageRate;
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            Q(i,j) *= factor;
            
    // calculate the weights
    std::vector<mpq_class> w(6);
    for (int i=0, k=0; i<4; i++)
        {
        for (int j=i+1; j<4; j++)
            {
            w[k] = pi[i] * Q(i,j);
            k++;
            }
        }
                
    // move the weights to GMP rational numbers
    W[0] = w[0];
    W[1] = w[1];
    W[2] = w[2];
    W[3] = w[3];
    W[4] = w[4];
    s = W[0] + W[1] + W[2] + W[3] + W[4];
    W[5] = oneHalf - s;
    
    storedIsTimeReversible = isTimeReversible;
    for (int i=0; i<4; i++)
        storedPi[i] = pi[i];
    for (int i=0; i<6; i++)
        storedR[i] = r[i];
    storedQ = Q;
        

#   if 1
    std::cout << "W[AC] = " << W[0] << " " << W[0].get_d() << " " << w[0] << std::endl;
    std::cout << "W[AG] = " << W[1] << " " << W[1].get_d() << " " << w[1] << std::endl;
    std::cout << "W[AT] = " << W[2] << " " << W[2].get_d() << " " << w[2] << std::endl;
    std::cout << "W[CG] = " << W[3] << " " << W[3].get_d() << " " << w[3] << std::endl;
    std::cout << "W[CT] = " << W[4] << " " << W[4].get_d() << " " << w[4] << std::endl;
    std::cout << "W[GT] = " << W[5] << " " << W[5].get_d() << " " << w[5] << std::endl;
    std::cout << "Sum = " << W[0] + W[1] + W[2] + W[3] + W[4] + W[5] << std::endl;
#   endif
}

double McmcState::lnPriorProbability() {

    std::vector<double> bf(4);
    for (int i=0; i<4; i++)
        bf[i] = pi[i].get_d();
    double lnProb = Probability::Dirichlet::lnPdf(stationaryFrequenciesAlpha, bf);
    
    lnProb += 3.0 * (log(bf[0]) + log(bf[1]) + log(bf[2]) + log(bf[3]));
    if (isTimeReversible == true)
        {
        lnProb += log(120.0);
        lnProb += log(4.0) + 0.5 * log(3.0); // original
        //lnProb += log(4.0);  // temporary: this works...we're off by a factor of the square root of 3
        }
    else
        {
        lnProb += log(40320.0);
        lnProb -= 0.5 * log(6.0); // original
        
        //lnProb -= 0.5 * log(2.0); // temporary: this works as an alternative to the original factor, -= 0.5 * log(6.0)
        //lnProb += 0.5 * log(3.0); // temporary: gets very close to 50/50 if added
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

void McmcState::setNonReversibleRateMatrix(void) {

    // update pi for new rate matrix
    calculateStationaryFrequencies();
    
    // rescale rate matrix so average rate is one
    mpq_class averageRate;
    for (int i=0; i<4; i++)
        averageRate += -Q(i,i) * pi[i];
    mpq_class factor = 1 / averageRate;
    for (mpq_class* p=Q.begin(); p != Q.end(); p++)
        (*p) *= factor;
        
    if (checkRateMatrix() == false)
        Msg::error("Unexpected bad matrix in setNonReversibleRateMatrix");
}

void McmcState::setReversibleRateMatrix(void) {

    // store old value of rate matrix
    mpq_class* p1 = storedQ.begin();
    for (mpq_class* p2=Q.begin(); p2 != Q.end(); p2++)
        {
        *p1 = *p2;
        p1++;
        }
    
    // set the rate matrix
    for (int i=0, k=0; i<4; i++)
        {
        for (int j=i+1; j<4; j++)
            {
            Q(i,j) = r[k] * pi[j];
            Q(j,i) = r[k] * pi[i];
            k++;
            }
        }
    mpq_class averageRate = 0;
    for (int i=0; i<4; i++)
        {
        mpq_class sum = 0;
        for (int j=0; j<4; j++)
            {
            if (i != j)
                sum += Q(i,j);
            }
        Q(i,i) = -sum;
        averageRate += pi[i] * sum;
        }
    mpq_class factor = 1 / averageRate;
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            Q(i,j) *= factor;
}

std::string McmcState::stateString(void) {

    if (isTimeReversible == true)
        {
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
        std::stringstream ss;
        ss << std::fixed << std::setprecision(6);
        ss << "NR: ";
        for (int i=0; i<4; i++)
            ss << pi[i].get_d() << " ";
        /*for (int i=0; i<4; i++)
            {
            for (int j=0; j<4; j++)
                {
                if (i != j)
                    {
                    ss << Q(i,j).get_d() << " ";
                    }
                }
            }*/
        return ss.str();
        }
}

double McmcState::update(std::string& updateType) {

    double lnProb = 0.0;
    if (isTimeReversible == true)
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

    if (isTimeReversible == false)
        Msg::error("Can only update the exchangability rates for time reversible models");
        
    // store old values
    for (int i=0; i<6; i++)
        storedR[i] = r[i];
        
#   if 1
        
    std::vector<double> oldRates(6);
    std::vector<double> newRates(6);
    std::vector<double> alphaForward(6);
    std::vector<double> alphaReverse(6);
    double alpha0 = 100.0;
    for (int i=0; i<6; i++)
        oldRates[i] = r[i].get_d();
    for (int i=0; i<6; i++)
        alphaForward[i] = oldRates[i] * alpha0;
    Probability::Dirichlet::rv(rng, alphaForward, newRates);
    Probability::Helper::normalize(newRates, MIN_FREQ);

    for (int i=0; i<6; i++)
        alphaReverse[i] = newRates[i] * alpha0;
        
    r[0] = newRates[0];
    r[1] = newRates[1];
    r[2] = newRates[2];
    r[3] = newRates[3];
    r[4] = newRates[4];
    mpq_class sum = r[0] + r[1] + r[2] + r[3] + r[4];
    r[5] = 1 - sum;
    setReversibleRateMatrix();
    
    double lnProposalProb = Probability::Dirichlet::lnPdf(alphaReverse, oldRates) -
                            Probability::Dirichlet::lnPdf(alphaForward, newRates);
    return lnProposalProb;

#   else

    for (int i=0; i<6; i++)
        r[i] = exp(rng->uniformRv());
    normalize(r);
    setReversibleRateMatrix();
    return 0.0;
    
#   endif
}

void McmcState::updateForAcceptance(void) {

    storedIsTimeReversible = isTimeReversible;
}

void McmcState::updateForRejection(void) {

    isTimeReversible = storedIsTimeReversible;
    
    for (int i=0; i<4; i++)
        pi[i] = storedPi[i];
        
    for (int i=0; i<6; i++)
        r[i] = storedR[i];
        
    mpq_class* p1 = storedQ.begin();
    for (mpq_class* p2=Q.begin(); p2 != Q.end(); p2++)
        {
        *p2 = *p1;
        p1++;
        }
}

double McmcState::updateNonreversibleRates(void) {

    if (isTimeReversible == true)
        Msg::error("Can only update the 12 rates (directly) for non-reversible models");

    // store old values
    storedQ = Q;
    
#   if 1
    std::vector<double> oldRates(12);
    std::vector<double> newRates(12);
    std::vector<double> alphaForward(12);
    std::vector<double> alphaReverse(12);
    double alpha0 = 100.0;

    double sum = 0.0;
    for (int i=0, k=0; i<4; i++)
        {
        for (int j=0; j<4; j++)
            {
            if (i != j)
                {
                oldRates[k] = Q(i,j).get_d();
                sum += oldRates[k];
                k++;
                }
            }
        }
    for (int i=0; i<12; i++)
        oldRates[i] /= sum;
    for (int i=0; i<12; i++)
        alphaForward[i] = oldRates[i] * alpha0;
        
    Probability::Dirichlet::rv(rng, alphaForward, newRates);
    Probability::Helper::normalize(newRates, MIN_FREQ);
    
    for (int i=0; i<12; i++)
        alphaReverse[i] = newRates[i] * alpha0;
        
    for (int i=0, k=0; i<4; i++)
        {
        mpq_class sumQ;
        for (int j=0; j<4; j++)
            {
            if (i != j)
                {
                Q(i,j) = newRates[k];
                sumQ += Q(i,j);
                k++;
                }
            }
        Q(i,i) = -sumQ;
        }
    setNonReversibleRateMatrix();

    double lnProposalProb = Probability::Dirichlet::lnPdf(alphaReverse, oldRates)-
                            Probability::Dirichlet::lnPdf(alphaForward, newRates);
    return lnProposalProb;
    
#   else

    for (int i=0; i<4; i++)
        {
        mpq_class sum;
        for (int j=0; j<4; j++)
            {
            if (i != j)
                {
                Q(i,j) = exp(rng->uniformRv());
                sum += Q(i,j);
                }
            }
        Q(i,i) = -sum;
        }
    setNonReversibleRateMatrix();
    return 0.0;
    
#   endif
}

double McmcState::updateToNonReversible(void) {

    if (isTimeReversible == false)
        Msg::error("Can only move to a non-reversible model from a reversible one");

    // store old values
    storedQ = Q;
    storedIsTimeReversible = isTimeReversible;

    // calculate the Jacobian
    double piC = pi[C].get_d();
    double piG = pi[G].get_d();
    //double piT = pi[T].get_d();
    double wCG = piC * storedQ(C,G).get_d();
    double wCT = piC * storedQ(C,T).get_d();
    double wGT = piG * storedQ(G,T).get_d();
    //double lnJacobian = (64.0 * wCG * wCT * wGT) / (piC * piG * piG) ;
    double lnJacobian = (log(64.0) + log(wCG) + log(wCT) + log(wGT));
    lnJacobian -= log(piC) + 2.0 * log(piG);

    // polyhedron density for reverse move
    W[0] = (pi[A] * Q(A,C) + pi[C] * Q(C,A)) / 2;
    W[1] = (pi[A] * Q(A,G) + pi[G] * Q(G,A)) / 2;
    W[2] = (pi[A] * Q(A,T) + pi[T] * Q(T,A)) / 2;
    W[3] = (pi[C] * Q(C,G) + pi[G] * Q(G,C)) / 2;
    W[4] = (pi[C] * Q(C,T) + pi[T] * Q(T,C)) / 2;
    W[5] = (pi[G] * Q(G,T) + pi[T] * Q(T,G)) / 2;
    Vector randomPoint;
    mpf_class v = poly.volume(W,randomPoint);
    //v = 1 / v;
    double lnRv = log(v.get_d());
    
    // update rates
    mpq_class u1 = randomPoint.getX();
    mpq_class u2 = randomPoint.getY();
    mpq_class u3 = randomPoint.getZ();
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
    for (int i=0; i<4; i++)
        {
        mpq_class sum;
        for (int j=0; j<4; j++)
            {
            if (i != j)
                sum += Q(i,j);
            }
        Q(i,i) = -sum;
        }
        
    calculateStationaryFrequencies();
    
    mpq_class averageRate;
    for (int i=0; i<4; i++)
        averageRate += -pi[i] * Q(i,i);
    mpq_class factor = 1 / averageRate;
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            Q(i,j) *= factor;
    
    isTimeReversible = false;

    return lnRv + lnJacobian;
}

double McmcState::updateToReversible(void) {

    if (isTimeReversible == true)
        Msg::error("Can only move to a time reversible model from a non-reversible one");

    // store old values
    storedQ = Q;
    storedIsTimeReversible = isTimeReversible;
    
    // average rates with same stationary frequencies
    mpq_class averageRate;
    for (int i=0; i<4; i++)
        {
        mpq_class sum;
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
        
    // make certain average rate is one
    mpq_class factor = 1 / averageRate;
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            Q(i,j) *= factor;
        
    // jacobian
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
    
    // polyhedron density for reverse move
    W[0] = pi[A] * Q(A,C);
    W[1] = pi[A] * Q(A,G);
    W[2] = pi[A] * Q(A,T);
    W[3] = pi[C] * Q(C,G);
    W[4] = pi[C] * Q(C,T);
    W[5] = pi[G] * Q(G,T);
    mpf_class v = poly.volume(W);
    v = 1 / v;
    double lnRv = log(v.get_d());

    // figure out exchangability rates
    r[0] = Q(A,C) / pi[A];
    r[1] = Q(A,G) / pi[G];
    r[2] = Q(A,T) / pi[T];
    r[3] = Q(C,G) / pi[C];
    r[4] = Q(C,T) / pi[C];
    r[5] = Q(G,T) / pi[G];
    mpq_class sum;
    for (int i=0; i<6; i++)
        sum += r[i];
    for (int i=0; i<6; i++)
        r[i] /= sum;
    
    isTimeReversible = true;
    
    return lnRv + lnJacobian;
}

double McmcState::updateStationaryFrequencies(void) {

    if (isTimeReversible == false)
        Msg::error("Can only update the stationary frequencies (directly) for time reversible models");
        
    // store old values
    for (int i=0; i<4; i++)
        storedPi[i] = pi[i];
        
#   if 1
        
    std::vector<double> oldFreqs(4);
    std::vector<double> newFreqs(4);
    std::vector<double> alphaForward(4);
    std::vector<double> alphaReverse(4);
    double alpha0 = 100.0;
    for (int i=0; i<4; i++)
        oldFreqs[i] = pi[i].get_d();
    for (int i=0; i<4; i++)
        alphaForward[i] = oldFreqs[i] * alpha0;
    Probability::Dirichlet::rv(rng, alphaForward, newFreqs);
    Probability::Helper::normalize(newFreqs, MIN_FREQ);

    for (int i=0; i<4; i++)
        alphaReverse[i] = newFreqs[i] * alpha0;

    pi[0] = newFreqs[0];
    pi[1] = newFreqs[1];
    pi[2] = newFreqs[2];
    mpq_class s = pi[0] + pi[1] + pi[2];
    pi[3] = 1 - s;
    setReversibleRateMatrix();

    double lnProposalProb = Probability::Dirichlet::lnPdf(alphaReverse, oldFreqs)-
                            Probability::Dirichlet::lnPdf(alphaForward, newFreqs);
    return lnProposalProb;
    
#   else
    for (int i=0; i<4; i++)
        pi[i] = exp(rng->uniformRv());
    normalize(pi);
    setReversibleRateMatrix();
    return 0.0;
#   endif
}
