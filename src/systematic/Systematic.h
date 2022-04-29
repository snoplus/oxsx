/***************************************************************************************************/
/* A systematic object is responsible for manipulating the indicies of a pdf to match those of     */
/* its underlying detector response matrix so the response can be applied. And modifying the       */
/* response matrix according to some set of parameters.                                            */
/* The class knows the data representation of the indicies it affects, as well as the rep of       */
/*    the pdfs it acts on, this is nessecary to properly calulate its affect                       */
/* The construct() method is called in fit routines to recalculate the matrix after pdf adjustment */
/* Classes derived from this construct the response matrix according to some physical process.     */
/***************************************************************************************************/

#ifndef __OXSX_SYSTEMATIC__
#define __OXSX_SYSTEMATIC__
#include <SparseMatrix.h>
#include <BinnedED.h>
#include <FitComponent.h>
#include <ObsSet.h>
#include <vector>

class Systematic : public FitComponent{
 public:
    virtual ~Systematic()  {}

    BinnedED 
    operator()(const BinnedED& pdf_, double* norm=nullptr) const;
        
    void SetResponse(const SparseMatrix& responseMatrix_);
    const SparseMatrix& GetResponse() const;
        
    void   SetTransformationObs(const ObsSet&);
    ObsSet GetTransformationObs() const;

    void   SetDistributionObs(const ObsSet&);
    ObsSet GetDistributionObs() const;

    const AxisCollection& GetAxes() const;
    void  SetAxes(const AxisCollection& axes_);

    virtual void Construct() = 0;

 protected:
    SparseMatrix   fResponse;
    AxisCollection fAxes; // The full bin definition of the distributions
    ObsSet fTransObs;     // the observables to tranform
    ObsSet fDistObs;      // the full set of observables for the distriubtion to act on
    // obs of the pdfs it will act on, needs to be at least the lenth of the 
    // systematics representation
    
    // methods below used for index manipilation between pdf and response
    bool BinsCompatible(size_t bin1 , size_t bin2) const;
    bool VectorContains(const std::vector<size_t>&, size_t) const;
};
#endif
