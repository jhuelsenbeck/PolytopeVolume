#include "MpqRateMatrix.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "RbException.h"

#define MIN_FREQ    10e-4



MpqRateMatrix::MpqRateMatrix(void) {

    q = new mpq_class[16];
    endBuffer = q + 16;
    pi.resize(4);
    r.resize(6);
    isReversible = false;
}

MpqRateMatrix::MpqRateMatrix(const MpqRateMatrix& m) {
    
    this->isReversible = m.isReversible;
    
    q = new mpq_class[16];
    endBuffer = q + 16;
    pi.resize(4);
    r.resize(6);
    
    mpq_class* r = q;
    for (mpq_class* p=m.q; p != m.endBuffer; p++)
        {
        *r = *p;
        r++;
        }
    for (int i=0; i<4; i++)
        this->pi[i] = m.pi[i];
    for (int i=0; i<6; i++)
        this->r[i] = m.r[i];
}

MpqRateMatrix::~MpqRateMatrix(void) {

    delete [] q;
}

MpqRateMatrix& MpqRateMatrix::operator=(const MpqRateMatrix& rhs) {

    if (this != &rhs)
        {
        this->isReversible = rhs.isReversible;
        mpq_class* r = q;
        for (mpq_class* p=rhs.q; p != rhs.endBuffer; p++)
            {
            *r = *p;
            r++;
            }
        for (int i=0; i<4; i++)
            this->pi[i] = rhs.pi[i];
        for (int i=0; i<6; i++)
            this->r[i] = rhs.r[i];
        }
    return *this;
}

void MpqRateMatrix::adjust(void) {

    if (isReversible == true)
        {
        std::vector<double> adjR(6);
        for (int i=0; i<6; i++)
            adjR[i] = r[i].get_d();
        for (int i=0; i<5; i++)
            r[i] = adjR[i];
        r[5] = 1 - (r[0] + r[1] + r[2] + r[3] + r[4]);
        std::vector<double> adjPi(4);
        for (int i=0; i<4; i++)
            adjPi[i] = pi[i].get_d();
        for (int i=0; i<3; i++)
            pi[i] = adjPi[i];
        pi[3] = 1 - (pi[0] + pi[1] + pi[2]);
            
        for (int i=0, k=0; i<4; i++)
            {
            for (int j=i+1; j<4; j++)
                {
                (*this)(i,j) = r[k] * pi[j];
                (*this)(j,i) = r[k] * pi[i];
                k++;
                }
            }
        mpq_class averageRate;
        for (int i=0; i<4; i++)
            {
            mpq_class sum;
            for (int j=0; j<4; j++)
                {
                if (i != j)
                    sum += (*this)(i,j);
                }
            (*this)(i,i) = -sum;
            averageRate += pi[i] * sum;
            }
        mpq_class factor = 1 / averageRate;
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                (*this)(i,j) *= factor;
                
        }
    else
        {
        mpq_class adjQ[4][4];
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                adjQ[i][j] = (*this)(i,j);
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                (*this)(i,j) = adjQ[i][j];
                
        calculateStationaryFrequencies(this->pi);
                
        mpq_class averageRate;
        for (int i=0; i<4; i++)
            {
            mpq_class sum;
            for (int j=0; j<4; j++)
                {
                if (i != j)
                    sum += (*this)(i,j);
                }
            (*this)(i,i) = -sum;
            averageRate += pi[i] * sum;
            }
        mpq_class factor = 1 / averageRate;
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                (*this)(i,j) *= factor;
        }
}

void MpqRateMatrix::calculateAverageRate(mpq_class& ave) {

    ave = 0;
    for (int i=0; i<4; i++)
        ave += -(pi[i] * (*this)(i,i));
}

void MpqRateMatrix::calculateStationaryFrequencies(std::vector<mpq_class>& f) {

	// transpose the rate matrix (qMatrix) and put into QT
	MpqRateMatrix QT;
	transposeMatrix(*this, QT);

	// compute the LU decomposition of the transposed rate matrix
	MpqRateMatrix L;
	MpqRateMatrix U;
	computeLandU(QT, L, U);
	
	// back substitute into z = 0 to find un-normalized stationary frequencies
	// start with x_n = 1.0
	f[3] = 1;
	for (int i=4-2; i>=0; i--)
		{
		mpq_class dotProduct;
		for (int j=i+1; j<4; j++)
			dotProduct += U(i,j) * f[j];
		f[i] = (0 - dotProduct) / U(i,i);
		}
		
	// normalize the solution vector
	mpq_class sum;
	for (int i=0; i<4; i++)
		sum += f[i];
	for (int i=0; i<4; i++)
		f[i] /= sum;
  
    // make certain to initialize the instance variable
    for (int i=0; i<4; i++)
        this->pi[i] = f[i];
}

void MpqRateMatrix::calculateWeights(std::vector<mpq_class>& wts) {

    if (wts.size() != 6)
        throw(RbException("Weights array must have 6 elements"));
        
    // wts[0] = (pi[0] * (*this)(0,1) + pi[1] * (*this)(1,0)) / 2; // (pi[A] * (*this)(A,C) + pi[C] * (*this)(C,A)) / 2
    // wts[1] = (pi[0] * (*this)(0,2) + pi[2] * (*this)(2,0)) / 2; // (pi[A] * (*this)(A,G) + pi[G] * (*this)(G,A)) / 2
    // wts[2] = (pi[0] * (*this)(0,3) + pi[3] * (*this)(3,0)) / 2; // (pi[A] * (*this)(A,T) + pi[T] * (*this)(T,A)) / 2
    // wts[3] = (pi[1] * (*this)(1,2) + pi[2] * (*this)(2,1)) / 2; // (pi[C] * (*this)(C,G) + pi[G] * (*this)(G,C)) / 2
    // wts[4] = (pi[1] * (*this)(1,3) + pi[3] * (*this)(3,1)) / 2; // (pi[C] * (*this)(C,T) + pi[T] * (*this)(T,C)) / 2
    // wts[5] = (pi[2] * (*this)(2,3) + pi[3] * (*this)(3,2)) / 2; // (pi[G] * (*this)(G,T) + pi[T] * (*this)(T,G)) / 2
    for (int i=0, k=0; i<4; i++)
        {
        for (int j=i+1; j<4; j++)
            wts[k++] = (pi[i] * (*this)(i,j) + pi[j] * (*this)(j,i)) / 2;
        }
}

bool MpqRateMatrix::check(void) {

    mpq_class averageRate;
    for (int i=0; i<4; i++)
        {
        mpq_class sum;
        for (int j=0; j<4; j++)
            {
            if (i != j)
                sum += (*this)(i,j);
            }
        if ((*this)(i,i) != -sum)
            return false;
        averageRate += pi[i] * sum;
        }
    if (averageRate != 1)
        return false;
    return true;
}

void MpqRateMatrix::computeLandU(MpqRateMatrix& aMat, MpqRateMatrix& lMat, MpqRateMatrix& uMat) {

	for (int j=0; j<4; j++)
		{
		for (int k=0; k<j; k++)
			for (int i=k+1; i<j; i++)
				aMat(i,j) = aMat(i,j) - aMat(i,k) * aMat(k,j);

		for (int k=0; k<j; k++)
			for (int i=j; i<4; i++)
				aMat(i,j) = aMat(i,j) - aMat(i,k) * aMat(k,j);

		for (int m=j+1; m<4; m++)
	  		aMat(m,j) /= aMat(j,j);
		}

	for (int row=0; row<4; row++)
		{
		for (int col=0; col<4; col++)
			{
			if ( row <= col )
				{
				uMat(row,col) = aMat(row,col);
				lMat(row,col) = (row == col ? 1 : 0);
				}
			else
				{
				lMat(row,col) = aMat(row,col);
				uMat(row,col) = 0;
				}
			}
		}
}

void MpqRateMatrix::initializeTimeReversibleModel(RandomVariable* rng) {

    isReversible = true;
    
    mpq_class one = 1;
    
    // first, choose the stationary frequency
    std::vector<double> alpha4(4, 1.0);
    std::vector<double> bf(4);
    Probability::Dirichlet::rv(rng, alpha4, bf);
    Probability::Helper::normalize(bf, MIN_FREQ);
    this->pi[0] = bf[0];
    pi[1] = bf[1];
    pi[2] = bf[2];
    mpq_class s = pi[0] + pi[1] + pi[2];
    pi[3] = one - s;
    
    // also, choose the exchangability rates
    std::vector<double> alpha6(6, 1.0);
    std::vector<double> er(6);
    Probability::Dirichlet::rv(rng, alpha6, er);
    Probability::Helper::normalize(er, MIN_FREQ);
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
            (*this)(i,j) = this->r[k] * this->pi[j];
            (*this)(j,i) = this->r[k] * this->pi[i];
            k++;
            }
        }
        
    // set the diagonal and calculate the average rate
    mpq_class averageRate = 0;
    for (int i=0; i<4; i++)
        {
        mpq_class sum = 0;
        for (int j=0; j<4; j++)
            {
            if (i != j)
                sum += (*this)(i,j);
            }
        (*this)(i,i) = -sum;
        averageRate += pi[i] * sum;
        }
        
    // rescale the matrix such that the average rate is one
    mpq_class factor = 1 / averageRate;
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            (*this)(i,j) *= factor;
            
    // set the exchangeability rates
    setExchangeabilityRates();
}

void MpqRateMatrix::nonreversibilize(mpq_class& u1, mpq_class& u2, mpq_class& u3) {

    if (isReversible == false)
        throw(RbException("Cannot make a non-reversibilize model non-reversible (again)"));
        
    // we weren't nonreversible before, but we are now (or will be in a microsecond)
    isReversible = false;
    
    // update rates (off diagonal components)
    mpq_class qAC = (*this)(0,1) + (pi[1]/pi[0]) * (*this)(1,2) * (2 * u1 - 1) + (pi[1]/pi[0]) * (*this)(1,3) * (2 * u2 - 1);
    mpq_class qAG = (*this)(0,2) - (pi[1]/pi[0]) * (*this)(1,2) * (2 * u1 - 1) + (pi[2]/pi[0]) * (*this)(2,3) * (2 * u3 - 1);
    mpq_class qAT = (*this)(0,3) - (pi[1]/pi[0]) * (*this)(1,3) * (2 * u2 - 1) - (pi[2]/pi[0]) * (*this)(2,3) * (2 * u3 - 1);
    mpq_class qCG = 2 * (*this)(1,2) * u1;
    mpq_class qCT = 2 * (*this)(1,3) * u2;
    mpq_class qGT = 2 * (*this)(2,3) * u3;
    mpq_class qCA = (pi[0]/pi[1]) * (*this)(0,1) - (*this)(1,2) * (2 * u1 - 1) - (*this)(1,3) * (2 * u2 - 1);
    mpq_class qGA = (pi[0]/pi[2]) * (*this)(0,2) + (pi[1]/pi[2]) * (*this)(1,2) * (2 * u1 - 1) - (*this)(2,3) * (2 * u3 - 1);
    mpq_class qTA = (pi[0]/pi[3]) * (*this)(0,3) + (pi[1]/pi[3]) * (*this)(1,3) * (2 * u2 - 1) + (pi[2]/pi[3]) * (*this)(2,3) * (2 * u3 - 1);
    mpq_class qGC = 2 * (pi[1]/pi[2]) * (*this)(1,2) * (1 - u1);
    mpq_class qTC = 2 * (pi[1]/pi[3]) * (*this)(1,3) * (1 - u2);
    mpq_class qTG = 2 * (pi[2]/pi[3]) * (*this)(2,3) * (1 - u3);
    (*this)(0,1) = qAC;
    (*this)(0,2) = qAG;
    (*this)(0,3) = qAT;
    (*this)(1,0) = qCA;
    (*this)(1,2) = qCG;
    (*this)(1,3) = qCT;
    (*this)(2,0) = qGA;
    (*this)(2,1) = qGC;
    (*this)(2,3) = qGT;
    (*this)(3,0) = qTA;
    (*this)(3,1) = qTC;
    (*this)(3,2) = qTG;
    
    // update the diagonal components of the rate matrix
    mpq_class sum;
    mpq_class averageRate;
    for (int i=0; i<4; i++)
        {
        sum = 0;
        for (int j=0; j<4; j++)
            {
            if (i != j)
                sum += (*this)(i,j);
            }
        (*this)(i,i) = -sum;
        averageRate += pi[i] * sum;
        }
            
    // the average rate should remain one
    if (averageRate != 1)
        throw(RbException("Average rate should be one when moving to nonreversible model"));
}

void MpqRateMatrix::print(void) {

    std::cout << std::fixed << std::setprecision(8);
    for (int i=0; i<4; i++)
        {
        for (int j=0; j<4; j++)
            {
            if ((*this)(i,j) >= 0)
                std::cout << " ";
            std::cout << (*this)(i,j).get_d() << " ";
            }
        std::cout << std::endl;
        }
}

void MpqRateMatrix::reversibilize(void) {

    if (isReversible == true)
        throw(RbException("Cannot reversibilize a time-reversible rate matrix"));
        
    // we weren't reversible before, but we are now (or will be in a microsecond)
    isReversible = true;

    // set off diagonals
    mpq_class w;
    for (int i=0; i<4; i++)
        {
        for (int j=i+1; j<4; j++)
            {
            w = this->pi[i] * (*this)(i,j) + this->pi[j] * (*this)(j,i);
            (*this)(i,j) = w / (2 * this->pi[i]);
            (*this)(j,i) = w / (2 * this->pi[j]);
            }
        }
        
    // set diagonals
    mpq_class sum;
    mpq_class averageRate;
    for (int i=0; i<4; i++)
        {
        sum = 0;
        for (int j=0; j<4; j++)
            {
            if (i != j)
                sum += (*this)(i,j);
            }
        (*this)(i,i) = -sum;
        averageRate += pi[i] * sum;
        }
        
    setExchangeabilityRates();
        
    // make certain average rate is one
    if (averageRate != 1)
        throw(RbException("Average rate should be one when moving to a reversible model"));
}

void MpqRateMatrix::setExchangeabilityRates(void) {

    if (isReversible == false)
        throw(RbException("Cannot set exchangeability rates for a non-reversible rate matrix"));
        
    // this->r[0] = (*this)(0,1) / pi[1]; // r_AC = Q(A,C) / pi[C]
    // this->r[1] = (*this)(0,2) / pi[2]; // r_AG = Q(A,G) / pi[G]
    // this->r[2] = (*this)(0,3) / pi[3]; // r_AT = Q(A,T) / pi[T]
    // this->r[3] = (*this)(1,2) / pi[2]; // r_CG = Q(C,G) / pi[G]
    // this->r[4] = (*this)(1,3) / pi[3]; // r_CT = Q(C,T) / pi[T]
    // this->r[5] = (*this)(2,3) / pi[3]; // r_GT = Q(G,T) / pi[T]
        
    mpq_class sum;
    for (int i=0, k=0; i<4; i++)
        {
        for (int j=i+1; j<4; j++)
            {
            this->r[k] = (*this)(i,j) / pi[j];
            sum += this->r[k];
            k++;
            }
        }
    for (int i=0; i<6; i++)
        this->r[i] /= sum;
}

void MpqRateMatrix::setPi(std::vector<mpq_class>& f) {

    for (int i=0; i<4; i++)
        this->pi[i] = f[i];
}

void MpqRateMatrix::transposeMatrix(const MpqRateMatrix& a, MpqRateMatrix& t) {
	
	for (int i=0; i<4; i++)
		for (int j=0; j<4; j++)
			t(j,i) = a(i,j);
}

double MpqRateMatrix::updateNonReversibleRates(RandomVariable* rng, double alpha0) {

    if (isReversible == true)
        throw(RbException("Can only update the 12 rates (directly) for non-reversible models"));
    
    // update the rates
    std::vector<double> oldRates(12);
    std::vector<double> newRates(12);
    std::vector<double> alphaForward(12);
    std::vector<double> alphaReverse(12);
    double sum = 0.0;
    for (int i=0, k=0; i<4; i++)
        {
        for (int j=0; j<4; j++)
            {
            if (i != j)
                {
                oldRates[k] = (*this)(i,j).get_d();
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
        
    // update the rate matrix
    for (int i=0, k=0; i<4; i++)
        {
        mpq_class sumQ;
        for (int j=0; j<4; j++)
            {
            if (i != j)
                {
                (*this)(i,j) = newRates[k];
                sumQ += (*this)(i,j);
                k++;
                }
            }
        (*this)(i,i) = -sumQ;
        }
    
    // update pi for new rate matrix
    calculateStationaryFrequencies(pi);
    
    // rescale rate matrix so average rate is one
    mpq_class averageRate;
    calculateAverageRate(averageRate);
    mpq_class factor = 1 / averageRate;
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            (*this)(i,j) *= factor;

    // return the log of the proposal probability
    double lnProposalProb = Probability::Dirichlet::lnPdf(alphaReverse, oldRates)-
                            Probability::Dirichlet::lnPdf(alphaForward, newRates);
    return lnProposalProb;
}

double MpqRateMatrix::updateExchangeabilityRates(RandomVariable* rng, double alpha0) {

    if (isReversible == false)
        throw(RbException("Can only update the exchangability rates for time reversible models"));
            
    // update the rates
    std::vector<double> oldRates(6);
    std::vector<double> newRates(6);
    std::vector<double> alphaForward(6);
    std::vector<double> alphaReverse(6);
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

    // update the off-diagonal components of the rate matrix
    for (int i=0, k=0; i<4; i++)
        {
        for (int j=i+1; j<4; j++)
            {
            (*this)(i,j) = r[k] * pi[j];
            (*this)(j,i) = r[k] * pi[i];
            k++;
            }
        }

    // update the diagonals and calculate the average rate
    mpq_class averageRate;
    for (int i=0; i<4; i++)
        {
        sum = 0;
        for (int j=0; j<4; j++)
            {
            if (i != j)
                sum += (*this)(i,j);
            }
        (*this)(i,i) = -sum;
        averageRate += pi[i] * sum;
        }
        
    // rescale such that the average rate is one
    mpq_class factor = 1 / averageRate;
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            (*this)(i,j) *= factor;

    // return the log of the proposal probability
    double lnProposalProb = Probability::Dirichlet::lnPdf(alphaReverse, oldRates) -
                            Probability::Dirichlet::lnPdf(alphaForward, newRates);
    return lnProposalProb;
}

double MpqRateMatrix::updateStationaryFrequencies(RandomVariable* rng, double alpha0) {

    if (isReversible == false)
        throw(RbException("Can only update the stationary frequencies (directly) for time reversible models"));
                
    // update the frequencies
    std::vector<double> oldFreqs(4);
    std::vector<double> newFreqs(4);
    std::vector<double> alphaForward(4);
    std::vector<double> alphaReverse(4);
    for (int i=0; i<4; i++)
        oldFreqs[i] = pi[i].get_d();
    for (int i=0; i<4; i++)
        alphaForward[i] = oldFreqs[i] * alpha0;
    Probability::Dirichlet::rv(rng, alphaForward, newFreqs);
    Probability::Helper::normalize(newFreqs, MIN_FREQ);
    for (int i=0; i<4; i++)
        alphaReverse[i] = newFreqs[i] * alpha0;

    // change to GMP rationals
    pi[0] = newFreqs[0];
    pi[1] = newFreqs[1];
    pi[2] = newFreqs[2];
    mpq_class sum = pi[0] + pi[1] + pi[2];
    pi[3] = 1 - sum;
    
    // update the off-diagonal components of the rate matrix
    for (int i=0, k=0; i<4; i++)
        {
        for (int j=i+1; j<4; j++)
            {
            (*this)(i,j) = r[k] * pi[j];
            (*this)(j,i) = r[k] * pi[i];
            k++;
            }
        }

    // update the diagonals and calculate the average rate
    mpq_class averageRate;
    for (int i=0; i<4; i++)
        {
        sum = 0;
        for (int j=0; j<4; j++)
            {
            if (i != j)
                sum += (*this)(i,j);
            }
        (*this)(i,i) = -sum;
        averageRate += pi[i] * sum;
        }
        
    // rescale such that the average rate is one
    mpq_class factor = 1 / averageRate;
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            (*this)(i,j) *= factor;

    // return the log of the proposal probability
    double lnProposalProb = Probability::Dirichlet::lnPdf(alphaReverse, oldFreqs)-
                            Probability::Dirichlet::lnPdf(alphaForward, newFreqs);
    return lnProposalProb;
}
