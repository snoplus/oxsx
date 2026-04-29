/*********************************************************************************************/
/* A sparse matrix class, a wrapper for the underlying Armadillo sparse matrix object.       */
/* Its main use in OXO is for holding the response matrix when systematics are applied to    */
/* binned event distributions.                                                               */
/*********************************************************************************************/

#ifndef __OXSX_SPARSE_MATRIX__
#define __OXSX_SPARSE_MATRIX__
#include <armadillo>

class SparseMatrix
{
public:
   SparseMatrix() : fNRows(0), fNCols(0) {}
   SparseMatrix(size_t rows_, size_t cols_);
   std::vector<double> operator()(const std::vector<double> &input_) const;

   void SetComponent(size_t row_, size_t column_, double val_);
   double GetComponent(size_t row_, size_t column_) const;

   void SetComponents(const std::vector<long long unsigned int> &rowIndices_,
                      const std::vector<long long unsigned int> &colIndices_,
                      const std::vector<double> &values_);

   SparseMatrix operator*=(const SparseMatrix &other_);
   SparseMatrix operator*(const SparseMatrix &other_) const;
   size_t GetNRows() const { return fNRows; }
   size_t GetNCols() const { return fNCols; }
   void SetZeros();
   void SetToIdentity();
   void Scale(double);

   void Print(const std::string &prefix_ = "") const;
   void PrintDense(const std::string &prefix_ = "") const;

private:
   arma::sp_mat fArmaMat;
   size_t fNRows;
   size_t fNCols;
};
#endif
