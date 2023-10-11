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
        {
        std::cout << "hey!" << std::endl;
        throw(RbException("Cannot non-reversibilize a non-reversible rate matrix"));
        }
        
    // we weren't nonreversible before, but we are now (or will be in a microsecond)
    isReversible = false;
    
    enum States {A,C,G,T};

    // update rates (off diagonal components)
    mpq_class qAC = (*this)(A,C) + (pi[C]/pi[A]) * (*this)(C,G) * (2 * u1 - 1) + (pi[C]/pi[A]) * (*this)(C,T) * (2 * u2 - 1);
    mpq_class qAG = (*this)(A,G) - (pi[C]/pi[A]) * (*this)(C,G) * (2 * u1 - 1) + (pi[G]/pi[A]) * (*this)(G,T) * (2 * u3 - 1);
    mpq_class qAT = (*this)(A,T) - (pi[C]/pi[A]) * (*this)(C,T) * (2 * u2 - 1) - (pi[G]/pi[A]) * (*this)(G,T) * (2 * u3 - 1);
    mpq_class qCG = 2 * (*this)(C,G) * u1;
    mpq_class qCT = 2 * (*this)(C,T) * u2;
    mpq_class qGT = 2 * (*this)(G,T) * u3;
    mpq_class qCA = (pi[A]/pi[C]) * (*this)(A,C) - (*this)(C,G) * (2 * u1 - 1) - (*this)(C,T) * (2 * u2 - 1);
    mpq_class qGA = (pi[A]/pi[G]) * (*this)(A,G) + (pi[C]/pi[G]) * (*this)(C,G) * (2 * u1 - 1) - (*this)(G,T) * (2 * u3 - 1);
    mpq_class qTA = (pi[A]/pi[T]) * (*this)(A,T) + (pi[C]/pi[T]) * (*this)(C,T) * (2 * u2 - 1) + (pi[G]/pi[T]) * (*this)(G,T) * (2 * u3 - 1);
    mpq_class qGC = 2 * (pi[C]/pi[G]) * (*this)(C,G) * (1 - u1);
    mpq_class qTC = 2 * (pi[C]/pi[T]) * (*this)(C,T) * (1 - u2);
    mpq_class qTG = 2 * (pi[G]/pi[T]) * (*this)(G,T) * (1 - u3);
    (*this)(A,C) = qAC;
    (*this)(A,G) = qAG;
    (*this)(A,T) = qAT;
    (*this)(C,G) = qCG;
    (*this)(C,T) = qCT;
    (*this)(G,T) = qGT;
    (*this)(C,A) = qCA;
    (*this)(G,A) = qGA;
    (*this)(T,A) = qTA;
    (*this)(G,C) = qGC;
    (*this)(T,C) = qTC;
    (*this)(T,G) = qTG;
    
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
            
    // the average rate should remain one. This is just a sanity check
    if (averageRate != 1)
        {
        std::cout << "Warning: the average rate should be one in updateToReversible" << std::endl;
        mpq_class factor = 1 / averageRate;
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                (*this)(i,j) *= factor;
        }
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
        {
        std::cout << "hey!" << std::endl;
        throw(RbException("Cannot reversibilize a time-reversible rate matrix"));
        }
        
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
        
    // make certain average rate is one
    if (averageRate != 1)
        {
        mpq_class factor = 1 / averageRate;
        std::cout << "Warning: the average rate should be one in reversibilize" << std::endl;
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                (*this)(i,j) *= factor;
        }
}

void MpqRateMatrix::setExchangeabilityRates(void) {

    if (isReversible == false)
        throw(RbException("Cannot set exchangeability rates for a non-reversible rate matrix"));
        
    this->r[0] = (*this)(0,1) / pi[1]; // r_AC = Q(A,C) / pi[C]
    this->r[1] = (*this)(0,2) / pi[2]; // r_AG = Q(A,G) / pi[G]
    this->r[2] = (*this)(0,3) / pi[3]; // r_AT = Q(A,T) / pi[T]
    this->r[3] = (*this)(1,2) / pi[2]; // r_CG = Q(C,G) / pi[G]
    this->r[4] = (*this)(1,3) / pi[3]; // r_CT = Q(C,T) / pi[T]
    this->r[5] = (*this)(2,3) / pi[3]; // r_GT = Q(G,T) / pi[T]
}

void MpqRateMatrix::transposeMatrix(const MpqRateMatrix& a, MpqRateMatrix& t) {
	
	for (int i=0; i<4; i++)
		for (int j=0; j<4; j++)
			t(j,i) = a(i,j);
}
