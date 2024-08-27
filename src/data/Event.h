/***********************************************************************************************/
/* Represents the observations associated with a single event.                                 */
/* Wraps a vector of the data points. Pdfs can query it with their own representation of the   */
/* data to get the relevant indicies ordered correctly.                                        */
/* For example, the data may come as (energy, radius, time). But the pdf is only interested in */
/* (radius, energy) pairs so it has the representation (1, 0).                                 */
/* Event(4, 2, 1).ToRepresentation(pdfRep) = {2, 4}                                        */
/***********************************************************************************************/

#ifndef __OXSX_EVENT__
#define __OXSX_EVENT__
#include <cstdlib>
#include <vector>
#include <map>
#include <string>

class ObsSet;
class Event
{
public:
   Event(const std::vector<double> &obs_);

   double GetDatum(size_t) const;
   double GetDatum(const std::string &name_) const;

   void SetDatum(size_t, double val_);
   void SetDatum(const std::string &name_, double val_);

   std::vector<double> GetData() const;
   std::map<std::string, double> GetDataMap() const;

   size_t GetNObservables() const { return fNObservables; }

   // used to link this to the data set that knows what the numbers mean
   void SetObservableNames(std::vector<std::string> const *p_);
   const std::vector<std::string> &GetObservableNames() const;

   // the method used by pdfs etc to get the data they need
   std::vector<double> ToObsSet(const ObsSet &rep_) const;

private:
   size_t fNObservables;
   std::vector<double> fObservations;
   std::vector<std::string> const *fObsNames;
};
#endif
