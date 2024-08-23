/***************************************************************************************************/
/* The outer product of two pdfs is also a pdf composite = pdf1 * pdf2                             */
/* This implemented here using recursion, so any level of composition can be used.                 */
/* This class owns its own copies of constituent pdfs                                              */
/* Calls to operator() and Probability() are forwarded to the constituent pdfs, in the latter case */
/* each will know where to look from its data representation                                       */
/***************************************************************************************************/

#ifndef __COMPOSITE_ED__
#define __COMPOSITE_ED__
#include <EventDistribution.h>
#include <vector>
#include <string>

class Event;
class CompositeED : public EventDistribution
{
public:
   CompositeED(const EventDistribution *p1_, const EventDistribution *p2_);
   CompositeED(const std::vector<EventDistribution *> &pdfs_);
   virtual ~CompositeED();

   virtual double Probability(const Event &) const;
   virtual EventDistribution *Clone() const; // required for futher compositions

   virtual double Integral() const;
   virtual void Normalise();
   unsigned GetNDims() const;

   std::string GetName() const;
   void SetName(const std::string &);

private:
   std::string fName;
   std::vector<EventDistribution *> fDistPtrs;
};

CompositeED operator*(const EventDistribution &, const EventDistribution &);

#endif
