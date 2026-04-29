/*********************************************************************************************/
/* A (dense) matrix class, a wrapper for the underlying Armadillo (dense) matrix object.     */
/*********************************************************************************************/

#ifndef __OXSX_DENSE_MATRIX__
#define __OXSX_DENSE_MATRIX__
#include <armadillo>

class DenseMatrix
{
public:
   DenseMatrix() : fNRows(0), fNCols(0) {}
   DenseMatrix(size_t rows_, size_t cols_);

   std::vector<double> operator()(const std::vector<double> &input_) const;

   void SetComponent(size_t row_, size_t column_, double val_);
   double GetComponent(size_t row_, size_t column_) const;

   DenseMatrix operator*=(const DenseMatrix &other_);

   size_t GetNRows() const { return fNRows; }
   size_t GetNCols() const { return fNCols; }
   void SetZeros();
   void SetToIdentity();

   void SetSymmetricMatrix(const std::vector<double> &_input);

   void Print(const std::string &prefix_ = "");
   void PrintSparse(const std::string &prefix_ = "");

private:
   // N x M matrix
   size_t fNRows;
   size_t fNCols;
   arma::mat fArmaMat;
};
#endif
