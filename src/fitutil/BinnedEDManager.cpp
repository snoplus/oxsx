#include <BinnedEDManager.h>
#include <SystematicManager.h>
#include <BinnedEDShrinker.h>
#include <BinnedED.h>
#include <Exceptions.h>
#include <sstream>

unsigned 
BinnedEDManager::GetNPdfs() const{
    return fOriginalPdfs.size();
}

size_t
BinnedEDManager::GetNDims() const{
    return fNDims;
}

double 
BinnedEDManager::Probability(const Event& data_) const{
    double sum = 0;

    for(size_t i = 0; i < fWorkingPdfs.size(); i++){
        sum += fNormalisations.at(i) * fWorkingPdfs[i].Probability(data_);
    }

    return sum;
}

double
BinnedEDManager::BinProbability(size_t bin_) const{
    double sum = 0;
    try{
        for(size_t i = 0; i < fWorkingPdfs.size(); i++){
            sum += fNormalisations.at(i) * fWorkingPdfs.at(i).GetBinContent(bin_);
        }
    }
    catch(const std::out_of_range&){
        throw LogicError("BinnedEDManager:: Normalisation vector doesn't match pdf vector - are the normalisations set?");
    }
    return sum;
}


void
BinnedEDManager::SetNormalisations(const std::vector<double>& normalisations_){
    if (normalisations_.size() != fOriginalPdfs.size())
        throw LogicError("BinnedEDManager: number of norms doesn't match #pdfs");
    fNormalisations = normalisations_;
}

void
BinnedEDManager::ApplySystematics(const SystematicManager& sysMan_){
    // If there are no systematics dont do anything
    //  ( working pdfs = original pdfs from initialisation)

    if(!sysMan_.GetSystematicsGroup().size())
        return;

    sysMan_.DistortEDs(fOriginalPdfs,fWorkingPdfs);
}

const BinnedED&
BinnedEDManager::GetOriginalPdf(size_t index_) const{
    return fOriginalPdfs.at(index_);
}

void
BinnedEDManager::AddPdf(const BinnedED& pdf_){
    fOriginalPdfs.push_back(pdf_);
    fWorkingPdfs.push_back(pdf_);
    fNPdfs++;
    fNormalisations.resize(fOriginalPdfs.size(), 0);
    RegisterParameters();
}

void 
BinnedEDManager::AddPdfs(const std::vector<BinnedED>& pdfs_){
    for(size_t i = 0; i < pdfs_.size(); i++){
        AddPdf(pdfs_.at(i));
    }
    RegisterParameters();
}

const std::vector<double>&
BinnedEDManager::GetNormalisations() const{
    return fNormalisations;
}

void
BinnedEDManager::ApplyShrink(const BinnedEDShrinker& shrinker_){
    if (!shrinker_.GetBuffers().size())
        return;
        
    // only shrink if not already shrunk! FIXME: more obvious behaviour
    if (!fWorkingPdfs.size() || fWorkingPdfs.at(0).GetNBins() != fOriginalPdfs.at(0).GetNBins())
        return;

    for (size_t i = 0; i < fWorkingPdfs.size(); i++){
        fWorkingPdfs[i] = shrinker_.ShrinkDist(fWorkingPdfs.at(i));
        fWorkingPdfs[i].Normalise();
    }
    
}

////////////////////////////////
// Make this a fit component! //
////////////////////////////////

std::string
BinnedEDManager::GetName() const{
    return fName;
}

void
BinnedEDManager::SetName(const std::string& n_){
    fName = n_;
}

void
BinnedEDManager::RenameParameter(const std::string& old_, const std::string& new_){
    fParameterManager.RenameParameter(old_, new_);
}

void
BinnedEDManager::SetParameter(const std::string& name_, double value_){
    fParameterManager.SetParameter(name_, value_);
}

double
BinnedEDManager::GetParameter(const std::string& name_) const{
    return fParameterManager.GetParameter(name_);
}

void
BinnedEDManager::SetParameters(const ParameterDict& ps_){
    fParameterManager.SetParameters(ps_);
}

ParameterDict
BinnedEDManager::GetParameters() const{
    return fParameterManager.GetParameters();
}

size_t
BinnedEDManager::GetParameterCount() const{
    return fParameterManager.GetParameterCount();
}

std::set<std::string>
BinnedEDManager::GetParameterNames() const{
    return fParameterManager.GetParameterNames();
}

void
BinnedEDManager::RegisterParameters(){
    fParameterManager.Clear();
    std::vector<std::string> parameterNames;
    for(size_t i = 0; i < fOriginalPdfs.size(); i++)
        parameterNames.push_back(fOriginalPdfs.at(i).GetName() + "_norm");
    fParameterManager.AddContainer(fNormalisations, parameterNames);
}    
