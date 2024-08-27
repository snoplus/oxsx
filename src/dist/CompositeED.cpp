#include <CompositeED.h>
#include <Event.h>
#include <Exceptions.h>
#include <iostream>

CompositeED::CompositeED(const EventDistribution *p1_, const EventDistribution *p2_)
{
    fDistPtrs.push_back(p1_->Clone());
    fDistPtrs.push_back(p2_->Clone());
    fName = p1_->GetName() + "*" + p2_->GetName();
}

CompositeED::CompositeED(const std::vector<EventDistribution *> &pdfs_)
{
    // if one of the pdfs is composite itself the copy will happen recursively all the way down
    std::string name;
    for (size_t i = 0; i < pdfs_.size(); i++)
    {
        fDistPtrs.push_back(pdfs_[i]->Clone());
        name += pdfs_[i]->GetName();
        if (i != pdfs_.size() - 1)
            name += "*";
    }
    fName = name;
}
CompositeED::~CompositeED()
{
    for (size_t i = 0; i < fDistPtrs.size(); i++)
        delete fDistPtrs[i];
}

double
CompositeED::Probability(const Event &data_) const
{
    double prob = 1;
    for (size_t i = 0; i < fDistPtrs.size(); i++)
        prob *= fDistPtrs[i]->Probability(data_);
    return prob;
}

void CompositeED::Normalise()
{
    for (size_t i = 0; i < fDistPtrs.size(); i++)
        fDistPtrs[i]->Normalise();
}

double
CompositeED::Integral() const
{
    double integral = 1;
    for (size_t i = 0; i < fDistPtrs.size(); i++)
        integral *= fDistPtrs[i]->Integral();
    return integral;
}

EventDistribution *
CompositeED::Clone() const
{
    EventDistribution *cp = new CompositeED(fDistPtrs);
    return static_cast<EventDistribution *>(cp);
}

CompositeED
operator*(const EventDistribution &pdf1_, const EventDistribution &pdf2_)
{
    return CompositeED(&pdf1_, &pdf2_);
}

unsigned
CompositeED::GetNDims() const
{
    unsigned nDims = 0;
    for (size_t i = 0; i < fDistPtrs.size(); i++)
        nDims += fDistPtrs.at(i)->GetNDims();
    return nDims;
}

std::string
CompositeED::GetName() const
{
    return fName;
}

void CompositeED::SetName(const std::string &name_)
{
    fName = name_;
}
