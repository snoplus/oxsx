/*********************************************************************************************/
/* A square response Matix for the experiment. Takes a binned pdf and applies the detector   */
/* to produce a new BinnedPhysDist. Inside is a vector of vectors, component fResponse[i][j] =    */
/* R_i_j = fraction of contents in bin j of original pdf -> bin i in the output pdf          */
/* the output bin contents are then x'_j = sum(R_i_j * x_j)                                  */
/* The systematic object is responsible for pdf renormalisation - not here                   */
/*********************************************************************************************/

#ifndef __OXSX_SPARSE_MATRIX__
#define __OXSX_SPARSE_MATRIX__
#include <AxisCollection.h>
#define ARMA_DONT_USE_CXX11
#include <armadillo>
class BinnedPhysDist;

class SparseMatrix{
 public:
    SparseMatrix() : fNRows(0), fNCols(0) {}
    SparseMatrix(int rows_, int cols_);
    std::vector<double> operator() (const std::vector<double>& input_) const;

    void   SetComponent(size_t row_, size_t column_, double val_);
    double GetComponent(size_t row_, size_t column_) const;

    void   SetComponents(const std::vector<unsigned>& rowIndices_,
                         const std::vector<unsigned>& colIndices_,
                       const std::vector<double>& values_);

    SparseMatrix operator*=(const SparseMatrix& other_);
    SparseMatrix operator*(const SparseMatrix& other_);
    size_t GetNRows() const {return fNRows;}
    size_t GetNCols() const {return fNCols;}
    void   SetZeros();
    void   SetToIdentity();

    void   Print(const std::string&);
    void   PrintDense(const std::string&);

 private:
    arma::sp_mat fArmaMat;
    int fNRows;
    int fNCols;
};
#endif
