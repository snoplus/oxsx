/*********************************************************************************************/
/* A square response Matix for the experiment. Takes a binned pdf and applies the detector   */
/* to produce a new BinnedPhysDist. Inside is a vector of vectors, component fResponse[i][j] =    */
/* R_i_j = fraction of contents in bin j of original pdf -> bin i in the output pdf          */
/* the output bin contents are then x'_j = sum(R_i_j * x_j)                                  */
/* The systematic object is responsible for pdf renormalisation - not here                   */
/*********************************************************************************************/

#ifndef __OXSX_DENSE_MATRIX__
#define __OXSX_DENSE_MATRIX__
#include <AxisCollection.h>
#define ARMA_DONT_USE_CXX11
#include <armadillo>
class BinnedPhysDist;

class DenseMatrix{
 public:
 DenseMatrix() : fNRows(0), fNCols(0) {}
    DenseMatrix(int rows_, int cols_);

    std::vector<double> operator() (const std::vector<double>& input_) const;

    void   SetComponent(size_t row_, size_t column_, double val_);
    double GetComponent(size_t row_, size_t column_) const;
    
    DenseMatrix operator*=(const DenseMatrix& other_);
   
    void   SetZeros();
    void   SetToIdentity();

    void SetSymmetricMatrix(const std::vector<double>& _input);

    void   Print(const std::string&);
    void   PrintSparse(const std::string&);

 private:
    // N x M matrix
    int fNRows;
    int fNCols;
    arma::mat fArmaMat;

};
#endif
