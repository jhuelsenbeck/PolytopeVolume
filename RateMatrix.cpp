#include "RateMatrix.hpp"



RateMatrix::RateMatrix(void) {

    q = new mpq_class[16];
    endBuffer = q + 16;
}

RateMatrix::RateMatrix(const RateMatrix& m) {
    
    q = new mpq_class[16];
    endBuffer = q + 16;
    
    mpq_class* r = q;
    for (mpq_class* p=m.begin(); p != m.end(); p++)
        {
        *r = *p;
        r++;
        }
}

RateMatrix::~RateMatrix(void) {

    delete [] q;
}

RateMatrix& RateMatrix::operator=(const RateMatrix& rhs) {

    if (this != &rhs)
        {
        mpq_class* r = q;
        for (mpq_class* p=rhs.begin(); p != rhs.end(); p++)
            {
            *r = *p;
            r++;
            }
        }
    return *this;
}

void RateMatrix::calculateStationaryFrequencies(std::vector<mpq_class>& f) {

	// transpose the rate matrix (qMatrix) and put into QT
	RateMatrix QT;
	transposeMatrix(*this, QT);

	// compute the LU decomposition of the transposed rate matrix
	RateMatrix L;
	RateMatrix U;
	computeLandU(QT, L, U);
	
	// back substitute into z = 0 to find un-normalized stationary frequencies
	// start with x_n = 1.0
	f[3] = 1;
	for (int i=4-2; i>=0; i--)
		{
		mpq_class dotProduct = 0.0;
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
}

void RateMatrix::computeLandU(RateMatrix& aMat, RateMatrix& lMat, RateMatrix& uMat) {

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

void RateMatrix::transposeMatrix(const RateMatrix& a, RateMatrix& t) {
	
	for (int i=0; i<4; i++)
		for (int j=0; j<4; j++)
			t(j,i) = a(i,j);
}
