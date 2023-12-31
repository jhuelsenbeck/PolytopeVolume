#include <gmpxx.h>
#include <iomanip>
#include <iostream>
#include <vector>
#include "McmcState.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"

#define MIN_FREQ    10e-10

std::vector<mpq_class> initializeRateMatrix(void);



int main(int argc, const char* argv[]) {

    // interface, such as it is
    int numCycles = 1000000;

    McmcState state;
    state.setAlphaT(1.0);
    double lnLikelihoodRatio = 0.0;
    double curLnPrior = state.lnPriorProbability();
    int reversibleCount = 0;
    for (int n=1; n<=numCycles; n++)
        {
        std::string updateType;
        double lnProposalRatio = state.update(updateType);
        double newLnPrior = state.lnPriorProbability();
        double lnPriorRatio = newLnPrior - curLnPrior;

        std::cout << std::fixed << std::setprecision(6);
        std::cout << std::setw(5) << n << " -- " << std::setw(9) << lnPriorRatio << " " << std::setw(9) << lnProposalRatio << " " << updateType << " -- ";
        if (state.accept(lnLikelihoodRatio + lnPriorRatio + lnProposalRatio) == true)
            {
            state.updateForAcceptance();
            curLnPrior = newLnPrior;
            std::cout << "Accepted -- ";
            }
        else
            {
            state.updateForRejection();
            std::cout << "Rejected -- ";
            }
        std::cout << std::fixed << std::setprecision(3) << (double)reversibleCount / n << " ";
        std::cout << state.stateString();
        std::cout << std::endl;
        
        if (state.getIsTimeReversible() == true)
            reversibleCount++;
            
        if (n % 1000 == 0)
            state.adjust();
        }
    std::cout << "Pr[Reversible] = " << (double)reversibleCount / numCycles << std::endl;
    
    
    
    
    
#if 0
    Vector randomPoint;
    Polyhedron poly;
    std::vector<mpq_class> W = initializeRateMatrix();
    poly.print(W);
#endif
        
        
    
    return 0;
}

std::vector<mpq_class> initializeRateMatrix(void) {

    // randomly initialize parameters of the GTR model
    std::vector<double> f(4);
    std::vector<double> r(6);
    std::vector<double> alpha4(4, 1.0);
    std::vector<double> alpha6(6, 1.0);
    RandomVariable& rng = RandomVariable::getInstance();
    Probability::Dirichlet::rv(&rng, alpha4, f);
    Probability::Dirichlet::rv(&rng, alpha6, r);
    
    // check bounds
    for (int i=0; i<4; i++)
        {
        if (f[i] < MIN_FREQ)
            f[i] = MIN_FREQ;
        }
    for (int i=0; i<6; i++)
        {
        if (r[i] < MIN_FREQ)
            r[i] = MIN_FREQ;
        }
    double sumD = 0.0;
    for (int i=0; i<4; i++)
        sumD += f[i];
    for (int i=0; i<4; i++)
        f[i] /= sumD;
    sumD = 0.0;
    for (int i=0; i<6; i++)
        sumD += r[i];
    for (int i=0; i<6; i++)
        r[i] /= sumD;
    
    // set the rate matrix
    double Q[4][4];
    for (int i=0, k=0; i<4; i++)
        {
        for (int j=i+1; j<4; j++)
            {
            Q[i][j] = r[k] * f[j];
            Q[j][i] = r[k] * f[i];
            k++;
            }
        }
    double averageRate = 0.0;
    for (int i=0; i<4; i++)
        {
        double sum = 0.0;
        for (int j=0; j<4; j++)
            {
            if (i != j)
                sum += Q[i][j];
            }
        Q[i][i] = -sum;
        averageRate += f[i] * sum;
        }
    double factor = 1.0 / averageRate;
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            Q[i][j] *= factor;
            
    // calculate the weights
    std::vector<double> w(6);
    for (int i=0, k=0; i<4; i++)
        {
        for (int j=i+1; j<4; j++)
            {
            w[k] = f[i] * Q[i][j];
            k++;
            }
        }
                
    // move the weights to GMP rational numbers
    mpq_class s;
    mpq_class oneHalf;
    oneHalf = 1;
    oneHalf /= 2;
    std::vector<mpq_class> W;
    for (int i=0; i<5; i++)
        {
        mpq_class x(w[i]);
        W.push_back(x);
        s += x;
        }
        
    mpq_class wGT;
    wGT = oneHalf - s;
    W.push_back(wGT);

#   if 0
    mpq_class sum;
    for (int i=0; i<6; i++)
        {
        std::cout << "W[" << i << "] = " << W[i] << " " << W[i].get_d() << " " << w[i] << std::endl;
        sum += W[i];
        }
    std::cout << "W[ ] = " << sum << std::endl;
#   endif

    return W;
}

