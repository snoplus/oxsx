/***************************************************************************************************/
/* Writes data sets to disk as hdf5 files and loads them back into oxsx objects. All data sets are */
/* saved as a single data set - a 1D array. The attribute "n_observables" gives the number of      */
/* observables in each event, which is then used to unflatten the 1D array.                        */
/* The attribute "observed_quantities" gives the ordered list of parameters separated by fDelimiter*/
/***************************************************************************************************/

#ifndef __OXSX_IO__
#define __OXSX_IO__
#include <string>
#include <vector>

class OXSXDataSet;
class DataSet;
class Histogram;

class IO{
 public:
    static void SaveDataSet(const DataSet&, const std::string& filename_, const std::string& = "oxsx_saved");    
    static void SaveDataSetH5(const DataSet&, const std::string& filename_);
    static void SaveDataSetROOT(const DataSet&, const std::string& filename_, const std::string& = "oxsx_saved");
      

    static void SaveHistogramH5(const Histogram& , const std::string& filename_);
    static void SaveHistogramROOT(const Histogram& , const std::string& filename_, const std::string& = "oxsx_saved");
    static void SaveHistogram(const Histogram&, const std::string& filename_, const std::string& = "oxsx_saved");

    static OXSXDataSet* LoadDataSet(const std::string& filename_);
    static Histogram    LoadHistogram(const std::string& filename_);

    
 private:
    static std::string GetExt(const std::string&);
    const static char fDelimiter;

    // used for combining and unpacking the observable names into a delimited string for hdf5 attr
    static std::string FlattenStringVector(const std::vector<std::string>&);
    static std::vector<std::string> UnpackString(const std::string&);

};
#endif
