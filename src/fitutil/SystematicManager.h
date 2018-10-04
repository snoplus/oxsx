/******************************************************************************************************/
/* Manages a set of pdfs, recieves a list of parameters and passes them out to each of the systematics*/
/* and triggeres their reconstruction. Systematics inside are passed to a set of pdfs to change       */
/* them                                                                                               */
/******************************************************************************************************/

#ifndef __SYSTEMATIC_MANAGER__
#define __SYSTEMATIC_MANAGER__
#include <vector>
#include <set>
#include <Systematic.h>
#include <SparseMatrix.h>

class SystematicManager{
 public:
    SystematicManager(){
      std::vector<Systematic*> holder;
      fGroups[""]=holder;
    }
    ~SystematicManager() {}

    void Add(Systematic*);

    void Add(Systematic*,const std::string& groupName_ );

    const std::map<std::string,std::vector<Systematic*> >& GetSystematicsGroup() const;

    const std::vector<Systematic*>& GetSystematicsInGroup(const std::string & name) const;

    const std::set<std::string> GetGroupNames() const;
    
    const size_t GetNSystematics() const;
    const size_t GetNGroups() const;

    const size_t GetNSystematicsInGroup(const std::string& name_) const;

    const std::vector<std::string> GetSystematicsNamesInGroup(const std::string& name) const;

    const std::vector<std::string> GetGroup(const std::string& name) const;

    void AddDist(const BinnedED& pdf, const std::vector<std::string>& syss_);
    void AddDist(const BinnedED& pdf, const std::string& syss_);

    const SparseMatrix& GetTotalResponse(const std::string& groupName_ = "" ) const;

    void DistortEDs(std::vector<BinnedED>& OrigEDs,std::vector<BinnedED>& WorkingEDs) const;

    void Construct();
    
 private:
    size_t fNGroups;
    std::map<std::string,SparseMatrix> fTotalReponses;
    std::map<std::string,std::vector<Systematic*> > fGroups;

    std::map<std::string,std::vector<std::string> > fEDGroups;
    void UniqueSystematics(const std::vector<std::string>&);
    void checkAllOtherGroups(const Systematic* syss_);
    const size_t CountNSystematics() const;
};
#endif
