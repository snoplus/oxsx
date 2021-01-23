/******************************************************************************************************/
/* Manages a set of pdfs, recieves a list of parameters and passes them out to each of the systematics*/
/* and triggeres their reconstruction. Systematics inside are passed to a set of pdfs to change       */
/* them                                                                                               */
/*                                                                                                    */
/* A group logic is employed, where systematics can be associated with a group of pdfs. A global      */
/* group ("") is created upon initialisation, which gets applied first to all distributions added.    */
/* If no group is specified when adding a systematics, it is automatically assigned to this global    */
/* group.                                                                                             */
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
    
    size_t GetNSystematics() const;
    size_t GetNGroups() const;

    size_t GetNSystematicsInGroup(const std::string& name_) const;

    const std::vector<std::string> GetSystematicsNamesInGroup(const std::string& name) const;

    const std::vector<std::string> GetGroup(const std::string& name) const;

    const std::vector<Systematic*>& GetGlobalSystematics() const;

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
    size_t CountNSystematics() const;
};
#endif
