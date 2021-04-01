/***************************************************************************************/
/* Light weight class to represent a particular combination of the data 
   for a single event 
*/
/**************************************************************************************/
#ifndef __OBS_SET__
#define __OBS_SET__
#include <vector>
#include <string>
#include <map>

class ObsSet{
 public:    
    ObsSet(){}
    ~ObsSet(){}
    ObsSet(const std::vector<std::string>& names_);
    ObsSet(const std::string& name_);
    void SetNames(const std::vector<std::string>& indices_);
    
    const std::vector<std::string>& GetNames() const; 

    size_t GetNObservables() const;
    
    size_t GetIndex(const std::string& name_) const;

    // Call this something else
    std::vector<size_t> GetRelativeIndices(const ObsSet&) const; 
    
    // comparision
    bool operator==(const ObsSet& other_) const;
    bool operator!=(const ObsSet& other_) const;

 private:
    std::vector<std::string> fNames;
    std::map<std::string, size_t> fInverseMap; //cache this
};
#endif
