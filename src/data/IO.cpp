#include <IO.h>
#include <DataSet.h>
#include <OXSXDataSet.h>
#include <H5Cpp.h>
#include <iostream>
#include <cassert>
#include <Exceptions.h>
#include <Histogram.h>
#include <Formatter.hpp>
#include <TFile.h>
#include <TNtuple.h>
#include <DistTools.h>
#include <TH1D.h>
#include <TH2D.h>

const char IO::fDelimiter = ':';


std::string
IO::GetExt(const std::string& path_){
    if(path_.find_last_of(".") != std::string::npos)
        return path_.substr(path_.find_last_of(".") + 1);
    else
        throw IOError("Can't interpret extension of " + path_);
}

void
IO::SaveDataSet(const DataSet& dataSet_, const std::string& filename_, const std::string& treename_){
    std::string ext = GetExt(filename_);
    if(ext == "h5")
        SaveDataSetH5(dataSet_, filename_);

    else if(ext == "root")
        SaveDataSetROOT(dataSet_, filename_, treename_);

    else
        throw IOError("IO::SaveDataSet Don't know how to save file as ." + ext);
}


void
IO::SaveDataSetROOT(const DataSet& dataSet_, const std::string& filename_, const std::string& treename_){
    // note if this is already a ROOT tree this is a super silly way to do it, but it doesn't make sense to do it
    // anyway
    if(!dataSet_.GetObservableNames().size()){
        std::cout << "Tried to save a dataset with nothing in it.. nothing done" << std::endl;
        return;
    }

    std::string branchString = dataSet_.GetObservableNames().at(0);
    for(size_t i = 1; i < dataSet_.GetObservableNames().size(); i++){
        branchString += ":";
        branchString += dataSet_.GetObservableNames().at(i);
    }

    TFile f(filename_.c_str(), "RECREATE");
    TNtuple nt(treename_.c_str(), "", branchString.c_str());

    int oneTenth = dataSet_.GetNEntries()/10;
    for(size_t i = 0; i < dataSet_.GetNEntries(); i++){
        if(oneTenth && !(i%oneTenth))
	  std::cout << "(" <<  10 * i/oneTenth << "%)" << std::endl;

        // save these as floats
        Event ev = dataSet_.GetEntry(i);
	const std::vector<double>& datad = ev.GetData();
        std::vector<float> dataf(datad.begin(), datad.end());

        f.cd();
        // and write
        nt.Fill(&dataf[0]);
    }    
    f.cd();
    
    nt.Write();
    f.Close();
    
}


void
IO::SaveDataSetH5(const DataSet& dataSet_, const std::string& filename_){
    // Get the relevant params
    hsize_t nEntries     = dataSet_.GetNEntries();
    hsize_t nObs  = dataSet_.GetNObservables();
    hsize_t nData = nObs * nEntries;

    // create colon separated string from list of observables
    std::vector<std::string> observableNames = dataSet_.GetObservableNames();
    if (observableNames.size() != nObs)
        throw IOError("SaveDataSet::Require one name and one name only for each observable");
    
    // Set up files
    H5::H5File file(filename_, H5F_ACC_TRUNC);
    // Set up the data set
    // 1D, ndata long, called "observations". Saved as native doubles on this computer
    H5::DataSpace dataSpace(1, &nData);  

    // Data Set
    H5::DataSet   theData(file.createDataSet("observations", H5::PredType::NATIVE_DOUBLE, dataSpace));

    //  Set up the attributes - the number of obs per event and the names of the observables
    //  64 chars max in str to save space
    H5::StrType   strType(H5::PredType::C_S1, 64);
    H5::DataSpace attSpace(H5S_SCALAR);
    H5::Attribute obsListAtt = theData.createAttribute("observed_quantities", strType, attSpace);
    obsListAtt.write(strType, FlattenStringVector(observableNames));

    H5::Attribute countAtt = theData.createAttribute("n_observables",
                                                     H5::PredType::NATIVE_INT,
                                                     attSpace);
    countAtt.write(H5::PredType::NATIVE_INT, &nObs);
        
    //  Write the data
    //  Flatten data into 1D array
    //  HDF5 likes c arrays. Here use a vector and pass pointer to first element 
    //  memory guaranteed to be contiguous
    //  for very large datasets copying to flatten will give std::bad_alloc so use hyperslabs
    hsize_t startSlab[1];
    hsize_t count[1];
    
    hsize_t startEntry;
    hsize_t stopEntry;
    hsize_t nSlabEntries;
    
    std::vector<double> flattenedData;
    std::vector<double> eventData;

    // only both to split the data up into chunks of 1e5 to avoid a big copy
    int nSlabs = nData/10000 + 1;
    for(int i = 0; i < nSlabs; i++){
      /// work out the event range
      startEntry =  i * (nEntries / nSlabs);
      stopEntry  =  (i + 1) * (nEntries/ nSlabs);
      if(i == (nSlabs - 1))
        stopEntry = nEntries;
      
      nSlabEntries = stopEntry - startEntry;
      
      // work out the data range
      startSlab[0] = startEntry * nObs;
      count[0] = nSlabEntries * nObs;
      dataSpace.selectHyperslab(H5S_SELECT_SET, count, startSlab);
      H5::DataSpace mspace(1, count);

      // now flatten those events and write to the hyperslab
      flattenedData.clear();
      flattenedData.reserve(count[0]);
      for(hsize_t j = startEntry; j < stopEntry; j++){
        eventData = dataSet_.GetEntry(j).GetData();
        flattenedData.insert(flattenedData.end(), eventData.begin(), eventData.end());
      }
      theData.write(&flattenedData.at(0), H5::PredType::NATIVE_DOUBLE, mspace, dataSpace);
    } // slab loop
}

OXSXDataSet*
IO::LoadDataSet(const std::string& filename_){  
    std::cout << "IO::Loading " << filename_ << std::endl;
    // Get Data Set
    H5::H5File file;
    try{
      file.openFile(filename_, H5F_ACC_RDONLY);   
    }
    
    catch(const H5::FileIException&){
      throw IOError("Failed to open data set file : "  + filename_ + "\n check the path");
    }

    H5::DataSet dataSet = file.openDataSet("observations"); 

    // read meta information
    unsigned nObs = 0;
    H5::Attribute nameAtt  = dataSet.openAttribute("observed_quantities");
    H5::Attribute countAtt  = dataSet.openAttribute("n_observables");
    H5std_string strreadbuf("");
    nameAtt.read(nameAtt.getDataType(), strreadbuf);
    countAtt.read(countAtt.getDataType(), &nObs);    

    // Assemble into an OXSX data set
    OXSXDataSet* oxsxDataSet= new OXSXDataSet;

    // Set the variable names
    oxsxDataSet -> SetObservableNames(UnpackString(strreadbuf));

    // Read data out as 1D array
    hsize_t nData = 0;
    dataSet.getSpace().getSimpleExtentDims(&nData, NULL);
    size_t nEntries = nData/nObs;
    oxsxDataSet->Reserve(nEntries);

    // if the data set is small, just load it up all in one go
    if(nEntries < 1000000){
      std::vector<double> flatData(nData, 0);
      dataSet.read(&flatData.at(0), H5::PredType::NATIVE_DOUBLE);

      std::vector<double> oneEventObs(nObs, 0);
      for(size_t i = 0; i < nEntries; i++){
        for(size_t j = 0; j < nObs; j++)
          oneEventObs[j] = flatData.at(i * nObs + j);
        
        oxsxDataSet -> AddEntry(Event(oneEventObs));
      }
    }
    
    // For larger files read in chunks of 100000 events
    else{
      H5::DataSpace dataSpace = dataSet.getSpace();
      hsize_t extent[1] = {100000 * nObs};
      H5::DataSpace memSpace(1, extent);
      std::vector<double> flatData(100000 * nObs, 0);
      std::vector<double> oneEventObs(nObs, 0);

      hsize_t stride[1] = {1}; // move along to each entry in data set one at a time
      hsize_t count[1]  = {100000 * nObs}; // take 100 events at a time

      for(unsigned i = 0; i < unsigned(nEntries/100000); i++){
        hsize_t offset[1] = {i * 100000 * nObs};
        dataSpace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, NULL);    
        
        dataSet.read(&flatData.at(0), H5::PredType::NATIVE_DOUBLE, memSpace, dataSpace);
      
       for(size_t j = 0; j < 100000; j++){
         for(size_t k = 0; k < nObs; k++){
          oneEventObs[k] = flatData.at(j * nObs + k);
         }
         oxsxDataSet->AddEntry(Event(oneEventObs));
       }
          
      } // loop over chunks

      // Now take care of the left overs
      hsize_t nLeftOver = (nEntries % 100000);
      if(nLeftOver){
          extent[0] = nObs * nLeftOver;
          memSpace = H5::DataSpace(1, extent);
          flatData.resize(extent[0]);
          hsize_t offset[1] = {100000 * (nEntries/100000) * nObs};
          dataSpace.selectHyperslab(H5S_SELECT_SET, extent, offset, stride, NULL);
          dataSet.read(&flatData.at(0), H5::PredType::NATIVE_DOUBLE, memSpace, dataSpace);
          for(size_t i = 0; i < nLeftOver; i++){
              for(size_t j = 0; j < nObs; j++)
                  oneEventObs[j] = flatData.at(i * nObs + j);
              
              oxsxDataSet->AddEntry(Event(oneEventObs));
          }
      }
      
    } // else.. 
    return oxsxDataSet;
}

std::string
IO::FlattenStringVector(const std::vector<std::string>& vec_){
    std::string flatString;
    if(!vec_.size())
        return flatString;

    flatString += vec_.at(0);
    for(size_t i = 1; i < vec_.size(); i++)
        flatString +=  fDelimiter + vec_.at(i);

    return flatString;
}

std::vector<std::string>
IO::UnpackString(const std::string& str_){
    // count instances
    unsigned count = 0;
    for (size_t i = 0; i < str_.size(); i++)
        if(str_.at(i) == fDelimiter)
            count++;

    std::vector<std::string> strs(count + 1); // one more string than delimiter a:b
    size_t currentStr = 0;
    for(size_t i = 0; i < str_.size(); i++){
        if(str_[i] == fDelimiter)
            currentStr++;
        else
            strs[currentStr] += str_[i];
    }
    
    return strs;
}


void
IO::SaveHistogram(const Histogram& hist_, const std::string& filename_, const std::string& histname_){
    std::string ext = GetExt(filename_);
    if(ext == "h5")
        SaveHistogramH5(hist_, filename_);

    else if(ext == "root")
        SaveHistogramROOT(hist_, filename_, histname_);

    else
        throw IOError("IO::SaveDataSet Don't know how to save file as ." + ext);
}

void
IO::SaveHistogramROOT(const Histogram& hist_, const std::string& filename_, const std::string& histname_){
    int dim = hist_.GetNDims();
    if(dim == 1){
        TH1D hist = DistTools::ToTH1D(hist_);
        hist.SetName(histname_.c_str());
        hist.SaveAs(filename_.c_str());
    }
    else if(dim == 2){
        TH2D hist = DistTools::ToTH2D(hist_);
        hist.SetName(histname_.c_str());
        hist.SaveAs(filename_.c_str());
    }
    else
        throw IOError(Formatter() << "IO::SaveHistogramROOT don't know how to save a " << dim
                      << "histogram (just 1D or 2D)");
}

void 
IO::SaveHistogramH5(const Histogram& histo_, const std::string& fileName_){
    hsize_t nDims[] = {histo_.GetNDims()};
    hsize_t nBins[] = {histo_.GetNBins()};
    std::vector<double> contents = histo_.GetBinContents();

    // create the file
    H5::H5File file(fileName_, H5F_ACC_TRUNC);

     // Save nDims, nBins, names and latex names as attributes. Main data is the list of bin contents
    H5::DataSpace dataSpace(1, nBins);
    H5::DataSet   theData(file.createDataSet("histogram", H5::PredType::NATIVE_DOUBLE, dataSpace));

    // grab axis names
    std::vector<std::string> axisNames;
    std::vector<std::string> axisLatexNames;
    for(size_t i = 0; i < histo_.GetNDims(); i++){
        axisNames.push_back(histo_.GetAxes().GetAxis(i).GetName());
        axisLatexNames.push_back(histo_.GetAxes().GetAxis(i).GetLatexName());
    }

    // save axes names and latex names
    H5::StrType   strType(H5::PredType::C_S1, 64);
    H5::DataSpace attSpace(H5S_SCALAR);

    H5::Attribute nameListAtt = theData.createAttribute("axis_names", strType, attSpace);
    nameListAtt.write(strType, FlattenStringVector(axisNames));

    H5::Attribute latexNameListAtt = theData.createAttribute("axis_latex_names", strType, attSpace);
    latexNameListAtt.write(strType, FlattenStringVector(axisLatexNames));

    // save nbins, ndims
    H5::Attribute dimAtt = theData.createAttribute("n_dims",
                                                   H5::PredType::NATIVE_INT,
                                                   attSpace);
    H5::Attribute binAtt = theData.createAttribute("n_bins",
                                                   H5::PredType::NATIVE_INT,
                                                   attSpace);
    binAtt.write(H5::PredType::NATIVE_INT, &nBins);
    dimAtt.write(H5::PredType::NATIVE_INT, &nDims);

    // save the bin boundaries as attributes
    for(size_t i = 0; i < histo_.GetNDims(); i++){
        std::vector<double> lowEdges = histo_.GetAxes().GetAxis(i).GetBinLowEdges();
        std::vector<double> highEdges = histo_.GetAxes().GetAxis(i).GetBinHighEdges();
        hsize_t axisBins = histo_.GetAxes().GetAxis(i).GetNBins();
        H5::DataSpace binSpace(1, &axisBins);
        H5::Attribute lowEdgeAtt = theData.createAttribute(Formatter() << i << "_lows",
                                                           H5::PredType::NATIVE_DOUBLE,
                                                           binSpace);

        H5::Attribute highEdgeAtt = theData.createAttribute(Formatter() << i << "_highs",
                                                            H5::PredType::NATIVE_DOUBLE,
                                                            binSpace);

        lowEdgeAtt.write(H5::PredType::NATIVE_DOUBLE, &lowEdges[0]);
        highEdgeAtt.write(H5::PredType::NATIVE_DOUBLE, &highEdges[0]);
    }

    // now save the bin contents    
    theData.write(&contents[0], H5::PredType::NATIVE_DOUBLE, dataSpace);
    
    return;
}

Histogram
IO::LoadHistogram(const std::string& fileName_){
    std::cout << "IO::Loading " << fileName_ << std::endl;

    // Get Data Set
    H5::H5File file;
    try{
        file.openFile(fileName_, H5F_ACC_RDONLY);
    }

    catch(const H5::FileIException&){
        throw IOError("Failed to open histogram file : "  + fileName_ + "\n check the path");
    }
    
    H5::DataSet dataSet = file.openDataSet("histogram");
    
    // read the attributes
    int nDims;
    int nBins;
    H5std_string axisLatexNamesReadbuf("");
    H5std_string axisNamesReadbuf("");


    /// first the names
    H5::Attribute axisNameAtt  = dataSet.openAttribute("axis_names");
    H5::Attribute axisLatexNameAtt  = dataSet.openAttribute("axis_latex_names");
    axisNameAtt.read(axisNameAtt.getDataType(), axisNamesReadbuf);
    axisLatexNameAtt.read(axisLatexNameAtt.getDataType(), axisLatexNamesReadbuf);
    
    /// now the numbers
    H5::Attribute nBinsAtt  = dataSet.openAttribute("n_bins");
    H5::Attribute nDimsAtt  = dataSet.openAttribute("n_dims");
    nBinsAtt.read(nBinsAtt.getDataType(), &nBins);
    nDimsAtt.read(nDimsAtt.getDataType(), &nDims);
    
    //  now read the bin contents
    hsize_t nData = 0;
    dataSet.getSpace().getSimpleExtentDims(&nData, NULL);
    std::vector<double> binContents(nData, 0);
    dataSet.read(&binContents[0], H5::PredType::NATIVE_DOUBLE);

    // Now read the axis information and use it to build axes objects in memory
    // first build the axes
    AxisCollection axes;
    std::vector<std::string> axisNames = UnpackString(axisNamesReadbuf);
    std::vector<std::string> axisLatexNames = UnpackString(axisLatexNamesReadbuf);
    for(int i = 0; i < nDims; i++){
        // how many bins in this axis
        //        H5::Attribute axisBins = dataSet.openAttribute(Formatter() << i << "_nBins");
        
        // open up the attributes
        H5::Attribute lowEdgeAttribute = dataSet.openAttribute(Formatter() << i << "_lows");
        H5::Attribute highEdgeAttribute = dataSet.openAttribute(Formatter() << i << "_highs");

        /// work out how many bin edges there are
        hsize_t nLowEdges = 0;
        hsize_t nHighEdges = 0;

        lowEdgeAttribute.getSpace().getSimpleExtentDims(&nLowEdges, NULL);        
        highEdgeAttribute.getSpace().getSimpleExtentDims(&nHighEdges, NULL);
        
        // read them
        std::vector<double> lowEdges(nLowEdges, 0);
        std::vector<double> highEdges(nHighEdges, 0);
        H5::DataSpace dataSpace = lowEdgeAttribute.getSpace();
        lowEdgeAttribute.read(H5::PredType::NATIVE_DOUBLE, &lowEdges[0]);
        highEdgeAttribute.read(H5::PredType::NATIVE_DOUBLE, &highEdges[0]);

        axes.AddAxis(BinAxis(axisNames.at(i), lowEdges, highEdges, axisLatexNames.at(i)));
    }

    Histogram loadedHisto(axes);
    loadedHisto.SetBinContents(binContents);

    return loadedHisto;
}
