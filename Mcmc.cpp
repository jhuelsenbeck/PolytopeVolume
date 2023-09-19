#include <cmath>
#include <iostream>
#include <vector>
#include "Mcmc.hpp"
#include "Msg.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"

#define MIN_FREQ    10e-10



Mcmc::Mcmc(void) {
    
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
    update();
    std::cout << Q << std::endl;
    std::cout << storedQ << std::endl;
    
    std::vector<mpq_class> f(4);
    Q.calculateStationaryFrequencies(f);
    for (int i=0; i<4; i++)
        std::cout << "f[" << i << "] = " << f[i] << std::endl;

}

void Mcmc::calculateStationaryFrequencies(void) {


}

void Mcmc::initializeTimeReversibleRateMatrix(void) {

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

double Mcmc::lnPriorProbability(std::vector<mpq_class>& f) {

    std::vector<double> bf(4);
    for (int i=0; i<4; i++)
        bf[i] = f[i].get_d();
    double lnProb = Probability::Dirichlet::lnPdf(stationaryFrequenciesAlpha, bf);
    
    lnProb += 3.0 * (log(bf[0]) + log(bf[1]) + log(bf[2]) + log(bf[3]));
    if (isTimeReversible == true)
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

void Mcmc::seqQ(std::vector<mpq_class>& er, std::vector<mpq_class>& bf, RateMatrix& m) {

    // set the rate matrix
    for (int i=0, k=0; i<4; i++)
        {
        for (int j=i+1; j<4; j++)
            {
            m(i,j) = er[k] * bf[j];
            m(j,i) = er[k] * bf[i];
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
                sum += m(i,j);
            }
        m(i,i) = -sum;
        averageRate += bf[i] * sum;
        }
    mpq_class factor = 1 / averageRate;
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            m(i,j) *= factor;
}

void Mcmc::update(void) {

    if (isTimeReversible == true)
        {
        double u = rng->uniformRv();
        if (u < 0.5)
            updateStationaryFrequencies(pi);
        else
            updateExchangabilityRates(r);
        seqQ(r, pi, Q);
        }
}

double Mcmc::updateExchangabilityRates(std::vector<mpq_class>& er) {

    if (isTimeReversible == false)
        Msg::error("Can only update the exchangability rates for time reversible models");
        
    std::vector<double> oldRates(6);
    std::vector<double> newRates(6);
    for (int i=0; i<6; i++)
        oldRates[i] = r[i].get_d();
    double alpha0 = 100.0;
    std::vector<double> alphaForward(6);
    for (int i=0; i<6; i++)
        alphaForward[i] = oldRates[i] * alpha0;
    Probability::Dirichlet::rv(rng, alphaForward, newRates);

    std::vector<double> alphaReverse(6);
    for (int i=0; i<6; i++)
        alphaReverse[i] = newRates[i] * alpha0;
        
    er[0] = newRates[0];
    er[1] = newRates[1];
    er[2] = newRates[2];
    er[3] = newRates[3];
    er[4] = newRates[4];
    mpq_class sum = er[0] + er[1] + er[2] + er[3] + er[4];
    er[5] = 1 - sum;

    double lnProposalProb = Probability::Dirichlet::lnPdf(alphaReverse, oldRates);
    lnProposalProb -= Probability::Dirichlet::lnPdf(alphaForward, newRates);
    return lnProposalProb;
}

double Mcmc::updateStationaryFrequencies(std::vector<mpq_class>& bf) {

    if (isTimeReversible == false)
        Msg::error("Can only update the stationary frequencies (directly) for time reversible models");
        
    std::vector<double> oldFreqs(6);
    std::vector<double> newFreqs(6);
    for (int i=0; i<4; i++)
        oldFreqs[i] = pi[i].get_d();
    double alpha0 = 100.0;
    std::vector<double> alphaForward(4);
    for (int i=0; i<4; i++)
        alphaForward[i] = oldFreqs[i] * alpha0;
    Probability::Dirichlet::rv(rng, alphaForward, newFreqs);

    std::vector<double> alphaReverse(4);
    for (int i=0; i<4; i++)
        alphaReverse[i] = newFreqs[i] * alpha0;

    bf[0] = newFreqs[0];
    bf[1] = newFreqs[1];
    bf[2] = newFreqs[2];
    mpq_class s = bf[0] + bf[1] + bf[2];
    bf[3] = 1 - s;

    double lnProposalProb = Probability::Dirichlet::lnPdf(alphaReverse, oldFreqs);
    lnProposalProb -= Probability::Dirichlet::lnPdf(alphaForward, newFreqs);
    return lnProposalProb;
}

double Mcmc::updateRates(void) {

    if (isTimeReversible == true)
        Msg::error("Can only update the rates for a non-reversible models");

    return 0.0;
}
