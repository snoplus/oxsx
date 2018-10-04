#include  <SparseMatrix.h>
#include  <iostream>
#include  <Exceptions.h>
#include  <Formatter.hpp>

// Initialise to zeros
SparseMatrix::SparseMatrix(int rows_, int cols_){
    fNCols = cols_;
    fNRows  = rows_; 
    fArmaMat = arma::sp_mat(fNRows, fNCols);
}

void 
SparseMatrix::SetComponent(size_t row_, size_t col_, double val_){
    if (col_ >= fNCols || row_ >= fNRows)
        throw NotFoundError(Formatter() 
                            << "Attempted out of bounds access on  matrix (" 
                            << row_ <<  "," << col_ << ")."
                            << "Matrix is (" << fNRows << "x" << fNCols
                            << ")"
                            );

    fArmaMat(row_,col_) = val_;
}

double 
SparseMatrix::GetComponent(size_t row_, size_t col_) const{
    if (col_ >= fNCols || row_ >= fNRows)
        throw NotFoundError(Formatter() 
                            << "Attempted out of bounds access on  matrix (" 
                            << row_ <<  "," << col_ << ")."
                            << "Matrix is (" << fNRows << "x" << fNCols
                            << ")"
                            );
    return fArmaMat(row_, col_);
}

std::vector<double>
SparseMatrix::operator() (const std::vector<double>& input_) const{
    arma::vec newContents;
    try{
        // convert to armadillo vec
        newContents = fArmaMat * arma::vec(input_);
    }
    catch(const std::logic_error& e_){
        throw DimensionError(Formatter() << "DenseMatrix::opeator() : Input v\
ector ("
                             << input_.size() << ")"
                             << " wrong size for Matrix ("
                             << fNRows << "x" << fNCols
                             << " ) to act on");
    }

    // armadillo function for quick transfer to std::vector double
    return arma::conv_to<std::vector<double> >::from((newContents));
}


SparseMatrix
SparseMatrix::operator*=(const SparseMatrix& other_){
  fArmaMat = fArmaMat * other_.fArmaMat;
  return *this;
}

SparseMatrix
SparseMatrix::operator*(const SparseMatrix& other_){
  fArmaMat = fArmaMat * other_.fArmaMat;
  return *this;
}

void
SparseMatrix::SetZeros(){
    if(!fNRows || !fNCols)
        throw DimensionError(Formatter()<<
                "SparseMatrix:: Can't set elements to zero. (rows,cols) : ("<<
                fNRows<<","<<fNCols<<")"  
                );
    fArmaMat = arma::sp_mat(fNRows, fNCols);
}

void
SparseMatrix::SetToIdentity(){
    if(!fNRows || !fNCols || fNCols!=fNRows )
        throw DimensionError(Formatter()<<
                "SparseMatrix:: Can't set identity as matrix is not square. (rows,cols) : ("<<
                fNRows<<","<<fNCols<<")"  
                );
    fArmaMat.eye();
}

// FIXME: unsigned vs. size_t
void 
SparseMatrix::SetComponents(const std::vector<unsigned>& rowIndices_,
                          const std::vector<unsigned>& colIndices_,
                          const std::vector<double>& values_){
    if(rowIndices_.size() != values_.size() || colIndices_.size() != values_.size())
        throw DimensionError("SparseMatrix::SetComponent() #values != #locations");

    arma::umat locs(2, rowIndices_.size());
    locs.row(0) = arma::urowvec(rowIndices_);
    locs.row(1) = arma::urowvec(colIndices_);

    fArmaMat = arma::sp_mat(locs, arma::vec(values_));
}
