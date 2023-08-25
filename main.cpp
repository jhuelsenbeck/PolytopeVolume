#include <gmpxx.h>
#include <iostream>
#include <vector>
#include "Polyhedral.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"

std::vector<mpq_class> initializeRateMatrix(void);



int main(int argc, const char* argv[]) {

    // get the weights from the prior
    std::vector<mpq_class> W = initializeRateMatrix();

    Polyhedral poly(W);
    
    return 0;
}

std::vector<mpq_class> initializeRateMatrix(void) {

    // randomly initialize parameters of the GTR model
    std::vector<double> f(4);
    std::vector<double> r(6);
    std::vector<double> alpha4(4, 1.0);
    std::vector<double> alpha6(6, 1.0);
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    Probability::Dirichlet::rv(&rng, alpha4, f);
    Probability::Dirichlet::rv(&rng, alpha6, r);
    
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
    mpq_t s, oneHalf;
    mpq_inits(s, oneHalf, NULL);
    mpq_set_d(s, 0.0);
    mpq_set_d(oneHalf, 0.5);
    std::vector<mpq_class> W;
    for (int i=0; i<5; i++)
        {
        mpq_class x(w[i]);
        x.canonicalize();
        W.push_back(x);
        
        mpq_t x2;
        mpq_init(x2);
        mpq_set_d(x2, w[i]);
        mpq_add(s, x2, s);
        }
    mpq_t wGT;
    mpq_init(wGT);
    mpq_sub(wGT, oneHalf, s);
    mpq_add(s, wGT, s);
    mpq_class x(wGT);
    W.push_back(x);

    mpq_clears(s, oneHalf, wGT, NULL);

#   if 1
    mpq_t sum;
    mpq_init(sum);
    mpq_set_d(sum, 0.0);
    for (int i=0; i<6; i++)
        {
        std::cout << "W[" << i << "] = " << W[i] << " " << W[i].get_d() << std::endl;
        mpq_add(sum, W[i].get_mpq_t(), sum);
        }
    std::cout << "W[ ] = " << sum << std::endl;
    mpq_clear(sum);
#   endif

    return W;
}

