/******************************************************************************************************/
/* A group of Axes that defines the binning for a pdf.                                                */
/* Each bin is assigned a global bin ID, switch between this global ID and the indicies in each       */
/* axis using FlattenIndices and UnpackIndices                                                      */
/******************************************************************************************************/

#ifndef __AXIS_COLLECTION__
#define __AXIS_COLLECTION__

#include <BinAxis.h>
#include <vector>
#include <exception>

class AxisCollection{
    
 public:
    AxisCollection() : fNBins(0), fNDimensions(0) {}
    size_t  FindBin(const std::vector<double>& vals_) const;
    size_t  GetNBins() const;
    void    AddAxis(const BinAxis& axis_);
    void    AddAxes(const std::vector<BinAxis>& axes_);
    const BinAxis& GetAxis(size_t axisIndex_) const;    
    const BinAxis& GetAxis(const std::string& axisName_) const;
    
    size_t GetNDimensions() const;
    size_t FlattenIndices(const std::vector<size_t>& indicies_) const;
    size_t  UnflattenIndex(size_t index_, size_t dim_)  const;
    std::vector<size_t> UnpackIndices(size_t index_) const;

    void GetBinCentres(size_t bin_, std::vector<double>& output_) const;
    void GetBinLowEdges(size_t bin_, std::vector<double>& output_) const;
    void GetBinHighEdges(size_t bin_, std::vector<double>& output_) const;

    double GetBinLowEdge(size_t bin_, size_t dim_)  const;
    double GetBinHighEdge(size_t bin_, size_t dim_) const;
    double GetBinCentre(size_t bin_, size_t dim_)   const;
    double GetBinWidth(size_t bin_, size_t dim_)    const;
    
    // comparison
    bool operator==(const AxisCollection& other_) const;
    bool operator!=(const AxisCollection& other_) const;

    
    std::vector<std::string> GetAxisNames() const;

 private:
    std::vector<BinAxis> fAxes;
    std::vector<std::string> fAxisNames;
    std::vector<size_t> fAxisNbins;
    size_t fNBins;
    size_t fNDimensions;
    
    void    CountBins();
    bool    HasAxis(const std::string& name_);
};

#endif
