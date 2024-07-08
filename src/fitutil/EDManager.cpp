#include <EDManager.h>
#include <EventDistribution.h>
#include <Event.h>
#include <Exceptions.h>

EDManager::~EDManager(){
    for(size_t i = 0; i < fDists.size(); i++)
        delete fDists[i];
}

void 
EDManager::AddDist(EventDistribution * pdf_){
    if (!fDists.size())
        fNDims = pdf_->GetNDims();

    else if(pdf_->GetNDims() != fNDims)
        throw DimensionError("EDManager::AddDist", fNDims, pdf_->GetNDims(),
                             " dimensions in added pdf");

    fDists.push_back(pdf_->Clone());    
    fNDists++;
    fNormalisations.resize(fNDists, 0);
    RegisterParameters();
}

void
EDManager::AddDists(const std::vector<EventDistribution*>& pdfs_){
    for(size_t i = 0; i < pdfs_.size(); i++)
        AddDist(pdfs_.at(i));
    RegisterParameters();
}

double
EDManager::Probability(const Event& event_) const{
    double prob = 0;
    for(size_t i = 0; i < fDists.size(); i++)
        prob += fNormalisations.at(i) * fDists.at(i)->Probability(event_);
    
    return prob;
}

const std::vector<double>&
EDManager::GetNormalisations() const{
    return fNormalisations;
}

void
EDManager::SetNormalisations(const std::vector<double>& norms_){
    if (norms_.size() != fNDists)
        throw DimensionError("EDManager::SetNormalisations", fNDists, 
                             norms_.size());
    fNormalisations = norms_;
}

size_t 
EDManager::GetNDims() const{
    return fNDims;
}

size_t 
EDManager::GetNDists() const{
    return fNDists;
}

// Make a fittable component - i.e. rescale pdfs inside to fit

void
EDManager::RegisterParameters(){
    fParameterManager.Clear();
    std::vector<std::string> parameterNames;
    for(size_t i = 0; i < fDists.size(); i++)
        parameterNames.push_back(fDists.at(i)->GetName());
    
    fParameterManager.AddContainer(fNormalisations, parameterNames);
}    


std::string
EDManager::GetName() const{
    return fName;
}

void
EDManager::SetName(const std::string& n_){
    fName = n_;
}

void
EDManager::RenameParameter(const std::string& old_, const std::string& new_){
    fParameterManager.RenameParameter(old_, new_);
}

void
EDManager::SetParameter(const std::string& name_, double value_){
    fParameterManager.SetParameter(name_, value_);
}

double
EDManager::GetParameter(const std::string& name_) const{
    return fParameterManager.GetParameter(name_);
}

void
EDManager::SetParameters(const ParameterDict& ps_){
    fParameterManager.SetParameters(ps_);
}

ParameterDict
EDManager::GetParameters() const{
    return fParameterManager.GetParameters();
}

size_t
EDManager::GetParameterCount() const{
    return fParameterManager.GetParameterCount();
}

std::set<std::string>
EDManager::GetParameterNames() const{
    return fParameterManager.GetParameterNames();
}

