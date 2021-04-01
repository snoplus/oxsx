/******************************************************************************************/
/* Generate mixed data samples from a set of other data samples.                          */
/* Either take the expected number or a poisson fluctuation around the underlying rate    */
/* WARNING: The generator only owns pointers to the other data sets to avoid large copies */
/******************************************************************************************/

#ifndef __OXSX_DATA_SET_GENERATOR__
#define __OXSX_DATA_SET_GENERATOR__
#include <cstdlib>
#include <vector>
#include <CutCollection.h>

class OXSXDataSet;
class DataSet;
class Event;

class DataSetGenerator{
 public:
    DataSetGenerator() {}
    void SetDataSets(const std::vector<DataSet*> sets_);
    void SetExpectedRates(const std::vector<double>& rates_);

    void AddDataSet(DataSet* data_, double rates_, bool sequential_);
    
    void SetCuts(const CutCollection& cuts_);
    void AddCut(const Cut& cut_);

    const std::vector<bool>& GetBootstrap() const;
    void SetBootstrap(const std::vector<bool>&);

    OXSXDataSet ExpectedRatesDataSet(std::vector<int>* eventsTaken_ = NULL);
    OXSXDataSet PoissonFluctuatedDataSet(std::vector<int>* eventsTaken_ = NULL);
    OXSXDataSet AllValidEvents(std::vector<int>* eventsTaken_ = NULL);

    OXSXDataSet* AllRemainingEvents(size_t dataSet_, int* eventsTaken_);
    /* std::vector<OXSXDataSet*> AllRemainingEvents(std::vector<int>* eventsTaken_ = NULL); */
    void ClearDataSets();
    void Reset();
    
    void SetSequentialFlags(const std::vector<bool>&);
    const std::vector<bool>& GetSequentialFlags() const;

 private:
    std::vector<bool>        fSequentialFlags;
    std::vector<DataSet*>    fDataSets;
    std::vector<double>      fExpectedRates;
    std::vector<std::vector<size_t> > fEventIndicies;
    std::vector<size_t>      fMaxs;
    void                     RandomDrawsNoReplacement(size_t handleIndex_, int nEvents_, OXSXDataSet& data_);
    void                     SequentialDrawsNoReplacement(size_t handleIndex_, int nEvents_, OXSXDataSet& data_);
    void                     RandomDrawsWithReplacement(size_t handleIndex_, int nEvents_, OXSXDataSet& data_);
    std::vector<bool>        fBootstraps;
    CutCollection            fCuts;
};

#endif
