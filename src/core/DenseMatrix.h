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
#include <armadillo>
class BinnedPhysDist;

class DenseMatrix
{
public:
   DenseMatrix() : fNRows(0), fNCols(0) {}
   DenseMatrix(size_t rows_, size_t cols_);

   std::vector<double> operator()(const std::vector<double> &input_) const;

   void SetComponent(size_t row_, size_t column_, double val_);
   double GetComponent(size_t row_, size_t column_) const;

   size_t GetNRows() const { return fNRows; }
   size_t GetNCols() const { return fNCols; }

   DenseMatrix operator*=(const DenseMatrix &other_);

   void SetZeros();
   void SetToIdentity();

   void SetSymmetricMatrix(const std::vector<double> &_input);

   void Print(const std::string &);
   void PrintSparse(const std::string &);

   void SetMatrix(const arma::mat& input_matrix);
   arma::mat GetMatrix() const { return fArmaMat; }

private:
   // N x M matrix
   size_t fNRows;
   size_t fNCols;
   arma::mat fArmaMat;
};
#endif
