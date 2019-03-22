#include <BinnedNLLH.h>
#include <math.h>
#include <DataSet.h>
#include <Exceptions.h>
#include <DistFiller.h>
#include <CutLog.h>
#include <Exceptions.h>
#include <Formatter.hpp>
#include <iostream>

double 
BinnedNLLH::Evaluate(){
    if(!fDataSet && !fCalculatedDataDist) 
        throw LogicError("BinnedNNLH function called with no data set and no DataDist! set one of these first");
    
    if (!fCalculatedDataDist)
        BinData();
    
    if(!fAlreadyShrunk){
        fDataDist = fPdfShrinker.ShrinkDist(fDataDist);
        fAlreadyShrunk = true;
    }

    // Construct systematics 
    fSystematicManager.Construct(); 
    // Apply systematics
    fPdfManager.ApplySystematics(fSystematicManager);

    // Apply Shrinking
    fPdfManager.ApplyShrink(fPdfShrinker);

    // loop over bins and calculate the likelihood
    double nLogLH = 0;
    for(size_t i = 0; i < fDataDist.GetNBins(); i++){
        double prob = fPdfManager.BinProbability(i);
        if(!prob)
            throw std::runtime_error("BinnedNLLH::Encountered zero probability bin!");
        nLogLH -= fDataDist.GetBinContent(i) *  log(prob);        
    }

    // Extended LH correction
    const std::vector<double>& normalisations = fPdfManager.GetNormalisations();
    for(size_t i = 0; i < normalisations.size(); i++)
        nLogLH += normalisations.at(i);
            
    // Constraints
    for(std::map<std::string, QuadraticConstraint>::iterator it = fConstraints.begin();
        it != fConstraints.end(); ++it)
        nLogLH += it->second.Evaluate(fComponentManager.GetParameter(it->first));
   
    return nLogLH;
}

void
BinnedNLLH::BinData(){
    fDataDist =  BinnedED(fPdfManager.GetOriginalPdf(0)); // make a copy for same binning and data rep
    fDataDist.Empty();
    CutLog log(fCuts.GetCutNames());
    DistFiller::FillDist(fDataDist, *fDataSet, fCuts, log);
    fCalculatedDataDist = true;    
    fSignalCutLog = log;
}

void
BinnedNLLH::AddDist(const std::vector<BinnedED>& pdfs, const std::vector<std::vector<std::string> >& sys_){
    if (pdfs.size() != sys_.size())
       throw DimensionError(Formatter()<<"BinnedNLLH:: #sys_ != #group_");
    for (int i = 0; i < pdfs.size(); ++i)
        AddDist( pdfs.at(i), sys_.at(i) );
}

void
BinnedNLLH::AddDist(const BinnedED& pdf_, const std::vector<std::string>& syss_){
    fPdfManager.AddPdf(pdf_);
    fSystematicManager.AddDist(pdf_,syss_);
}

void
BinnedNLLH::AddDist(const BinnedED& pdf_){
    fPdfManager.AddPdf(pdf_);
    fSystematicManager.AddDist(pdf_,"");
}
void
BinnedNLLH::AddPdf(const BinnedED& pdf_){
    AddDist(pdf_);
}

void
BinnedNLLH::SetPdfManager(const BinnedEDManager& man_){
    fPdfManager = man_;
}

void
BinnedNLLH::SetSystematicManager(const SystematicManager& man_){
    fSystematicManager = man_;
}

void 
BinnedNLLH::AddSystematic(Systematic* sys_){
    fSystematicManager.Add(sys_);
}

void 
BinnedNLLH::AddSystematic(Systematic* sys_, const std::string&  group_){
    fSystematicManager.Add(sys_, group_);
}

void
BinnedNLLH::SetDataSet(DataSet* dataSet_){
    fDataSet = dataSet_;
    fCalculatedDataDist = false;
}

DataSet*
BinnedNLLH::GetDataSet(){
    return fDataSet;
}

void
BinnedNLLH::SetDataDist(const BinnedED& binnedPdf_){
    fDataDist = binnedPdf_;
    fCalculatedDataDist = true;
}

BinnedED
BinnedNLLH::GetDataDist() const{
    return fDataDist;
}


void
BinnedNLLH::SetBuffer(size_t dim_, unsigned lower_, unsigned upper_){
    fPdfShrinker.SetBuffer(dim_, lower_, upper_);
}

std::pair<unsigned, unsigned>
BinnedNLLH::GetBuffer(size_t dim_) const{
    return fPdfShrinker.GetBuffer(dim_);
}

void
BinnedNLLH::SetBufferAsOverflow(bool b_){
    fPdfShrinker.SetUsingOverflows(b_);
}

bool
BinnedNLLH::GetBufferAsOverflow() const{
    return fPdfShrinker.GetUsingOverflows();
}

void
BinnedNLLH::AddSystematics(const std::vector<Systematic*> systematics_){
    for(size_t i = 0; i < systematics_.size(); i++)
        AddSystematic(systematics_.at(i));
}

void
BinnedNLLH::AddSystematics(const std::vector<Systematic*> sys_, const std::vector<std::string> & groups_){
    if (groups_.size() != sys_.size())
       throw DimensionError(Formatter()<<"BinnedNLLH:: #sys_ != #group_");
    for(size_t i = 0; i <sys_.size(); i++)
        AddSystematic(sys_.at(i), groups_.at(i));
}

void
BinnedNLLH::SetNormalisations(const std::vector<double>& norms_){    
    fPdfManager.SetNormalisations(norms_);
}

std::vector<double>
BinnedNLLH::GetNormalisations() const{
    return fPdfManager.GetNormalisations();
}

void
BinnedNLLH::AddCut(const Cut& cut_){
    fCuts.AddCut(cut_);
}

void 
BinnedNLLH::SetCuts(const CutCollection& cuts_){
    fCuts = cuts_;
}

void 
BinnedNLLH::SetConstraint(const std::string& paramName_, double mean_, double sigma_){
    fConstraints[paramName_] = QuadraticConstraint(mean_, sigma_);
}


double
BinnedNLLH::GetSignalCutEfficiency() const{
    return fSignalCutEfficiency;
}

void
BinnedNLLH::SetSignalCutEfficiency(double eff_){
    fSignalCutEfficiency = eff_;
}

CutLog
BinnedNLLH::GetSignalCutLog() const{
    return fSignalCutLog;
}

void
BinnedNLLH::SetSignalCutLog(const CutLog& lg_){
    fSignalCutLog = lg_;
}

/////////////////////////////////////////////////////////
// Declare which objects should be adjusted by the fit //
/////////////////////////////////////////////////////////
void
BinnedNLLH::RegisterFitComponents(){
    fComponentManager.Clear();
    fComponentManager.AddComponent(&fPdfManager);
    
    //Because the limits are set by name they can be added in any order.
    const std::map<std::string, std::vector<Systematic*> > sys_ = fSystematicManager.GetSystematicsGroup();
    std::vector<std::string> alreadyAdded;
    for (std::map<std::string, std::vector<Systematic*> >::const_iterator group_ = sys_.begin(); group_ !=sys_.end(); ++group_) {
        for (int i = 0; i < group_->second.size(); ++i) {
            if( std::find( alreadyAdded.begin() , alreadyAdded.end() , group_->second.at(i)->GetName() ) == alreadyAdded.end() ){
                fComponentManager.AddComponent( group_->second.at(i) );
                alreadyAdded.push_back( group_->second.at(i)->GetName() );
            }
        }//End of group
    }//End of groups
}


void
BinnedNLLH::SetParameters(const ParameterDict& params_){
    try{
        fComponentManager.SetParameters(params_);
    }
    catch(const ParameterError& e_){
        throw ParameterError(std::string("BinnedNLLH::") + e_.what());
    }
}
                                             
ParameterDict
BinnedNLLH::GetParameters() const{
    return fComponentManager.GetParameters();
}

int
BinnedNLLH::GetParameterCount() const{
    return fComponentManager.GetTotalParameterCount();
}

std::set<std::string>
BinnedNLLH::GetParameterNames() const{
    return fComponentManager.GetParameterNames();
}
