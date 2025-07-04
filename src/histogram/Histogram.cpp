#include <Histogram.h>
#include <Exceptions.h>
#include <Formatter.hpp>
#include <Combinations.hpp>
#include <ContainerTools.hpp>
#include <iostream>
#include <set>
#include <algorithm>

Histogram::Histogram(const AxisCollection &axes_)
{
    SetAxes(axes_);
}

void Histogram::SetAxes(const AxisCollection &axes_)
{
    fAxes = axes_;
    fNBins = fAxes.GetNBins();
    fNDims = fAxes.GetNDimensions();
    fBinContents.resize(fNBins, 0);
}

const AxisCollection &
Histogram::GetAxes() const
{
    return fAxes;
}

double
Histogram::Integral() const
{
    double sum = 0;
    for (size_t i = 0; i < fNBins; i++)
        sum += fBinContents[i];
    return sum;
}

double
Histogram::Integral(double min, double max)
{
    size_t low = FindBin(min);
    size_t high = FindBin(max);
    double sum = 0;
    for (size_t i = low; i < high; i++)
        sum += fBinContents[i];
    return sum;
  }

void Histogram::Normalise()
{
    double sum = Integral();
    for (size_t i = 0; i < fNBins; i++)
        fBinContents[i] /= sum;
}

void Histogram::Scale(double s_)
{
    for (size_t i = 0; i < fNBins; i++)
        fBinContents[i] *= s_;
}

void Histogram::Fill(const std::vector<double> &vals_, double weight_)
{
    if (vals_.size() != fNDims)
        throw DimensionError("Histogram::Fill", fNDims, vals_.size());

    fBinContents[FindBin(vals_)] += weight_;
}

void Histogram::Fill(const std::map<std::string, double> &vals_, double weight_)
{
    try
    {
        Fill(ContainerTools::GetValues(vals_, GetAxisNames()), weight_);
    }
    catch (const std::out_of_range &e_)
    {
        throw NotFoundError("Tried to fill a histogram with incomplete dictionary!");
    }
}

void Histogram::Fill(double vals_, double weight_)
{
    Fill(std::vector<double>(1, vals_), weight_);
}

size_t
Histogram::FindBin(const std::vector<double> &vals_) const
{
    return fAxes.FindBin(vals_);
}

size_t
Histogram::FindBin(double value)
{
  return fAxes.GetAxis(0).FindBin(value);
}

double
Histogram::GetBinContent(size_t bin_) const
{
    if (bin_ > fNBins)
        throw NotFoundError(Formatter()
                            << "Out of bounds bin access attempted on bin "
                            << bin_ << " !");
    return fBinContents[bin_];
}

void Histogram::AddBinContent(size_t bin_, double content_)
{
    if (bin_ > fNBins)
        throw NotFoundError(Formatter()
                            << "Out of bounds bin access attempted on bin "
                            << bin_ << " !");
    fBinContents[bin_] += content_;
}

void Histogram::SetBinContent(size_t bin_, double content_)
{
    if (bin_ > fNBins)
        throw NotFoundError(Formatter() << "Out of bounds bin access attempted on bin " << bin_ << " !");
    fBinContents[bin_] = content_;
}

size_t
Histogram::GetNBins() const
{
    return fNBins;
}

size_t
Histogram::GetNDims() const
{
    return fNDims;
}

void Histogram::Empty()
{
    for (size_t i = 0; i < fNBins; i++)
        fBinContents[i] = 0;
}

size_t
Histogram::FlattenIndices(const std::vector<size_t> &indices_) const
{
    return fAxes.FlattenIndices(indices_);
}

std::vector<size_t>
Histogram::UnpackIndices(size_t bin_) const
{
    return fAxes.UnpackIndices(bin_);
}

std::vector<double>
Histogram::GetBinContents() const
{
    return fBinContents;
}
void Histogram::SetBinContents(const std::vector<double> &data_)
{
    if (data_.size() != fNBins)
        throw DimensionError("Histogram::SetBinContents", fNBins,
                             data_.size());

    fBinContents = data_;
    fNBins = fBinContents.size();
}

std::vector<double>
Histogram::Means() const
{
    std::vector<double> means(fNDims, 0);
    for (size_t i = 0; i < fNBins; i++)
        for (size_t j = 0; j < fNDims; j++)
            means[j] += fBinContents.at(i) * fAxes.GetBinCentre(i, j);
    return means;
}

std::vector<double>
Histogram::Variances() const
{
    std::vector<double> variances(fNDims, 0);

    for (size_t i = 0; i < fNBins; i++)
        for (size_t j = 0; j < fNDims; j++)
        {
            double binCent = fAxes.GetBinCentre(i, j);
            variances[j] += binCent * binCent * fBinContents.at(i);
        }

    std::vector<double> means = Means();
    for (size_t i = 0; i < fNDims; i++)
        variances[i] -= means.at(i) * means.at(i);

    return variances;
}

Histogram
Histogram::Marginalise(const std::vector<std::string> &axes_) const
{
    std::vector<std::string> allAxisNames = fAxes.GetAxisNames();

    // check the pdf does contain the axes asked for
    for (size_t i = 0; i < axes_.size(); i++)
    {

        if (std::find(allAxisNames.begin(), allAxisNames.end(), axes_.at(i)) == allAxisNames.end())
            throw NotFoundError(Formatter()
                                << "Histogram::Marginalise::Tried "
                                << "to project out non existent axis "
                                << axes_.at(i) << "!");
    }

    // work out which axis number corresponds to each name
    std::vector<size_t> indices;
    for (size_t i = 0; i < axes_.size(); i++)
        indices.push_back(fAxes.GetAxisIndex(axes_.at(i)));

    // Get the axes you are interested in, in the order requested
    AxisCollection newAxes;
    for (size_t i = 0; i < axes_.size(); i++)
        newAxes.AddAxis(fAxes.GetAxis(axes_.at(i)));

    // New histogram
    Histogram marginalised(newAxes);

    std::vector<size_t> oldIndices(fNDims);
    std::vector<size_t> newIndices(axes_.size());
    size_t newBin = -1;

    // Now loop over the bins in old and fill new pdfs
    for (size_t bin = 0; bin < fNBins; bin++)
    {
        for (size_t i = 0; i < fNDims; i++)
            oldIndices[i] = fAxes.UnflattenIndex(bin, i);

        for (size_t i = 0; i < indices.size(); i++)
            newIndices[i] = oldIndices.at(indices.at(i));

        newBin = marginalised.FlattenIndices(newIndices);
        marginalised.AddBinContent(newBin, fBinContents.at(bin));
    }
    return marginalised;
}

Histogram
Histogram::Marginalise(const std::string &index_) const
{
    return Marginalise(std::vector<std::string>(1, index_));
}

void Histogram::Recurse(size_t numFreeIdx, std::vector<size_t> binsEachAxis, std::vector<size_t> coords, std::vector<std::vector<size_t>> &localIdx, std::vector<size_t> free_ax) const
{
    /* Function finds every possible combination of free indexes for GetSlice method by creating N nested for loops.

       INPUTS: numFreeIdx, integer number of free idexes. This is the same as NDims of slice and number of nested
                           loops needed.
               binsEachAxis, integer number of bins. For each nested loop i, loop runs for binsEachAxis[i] times.
               coords, vector of dimension equal to NDims of original Histogram. Holds latest combination of fixed and free
                       axis idxs.
               localIdx, vector of vectors, stores every individual coords vector. To fill this is the ultiamte goal of Recurse
                         function.
               free_ax, vector of integers, specifies the axis location of each free axis, so correct element in coords vector
                        is filled.
    */

    // check if at deepest level
    if (numFreeIdx < 1)
    {
        // save the current coordinate combination
        localIdx.push_back(coords);
    }
    else
    {
        // need to go deeper! Begin this level's loop with loop variable running from 0 -> numBins in this axis
        for (size_t i = 0; i < binsEachAxis.at(numFreeIdx - 1); i++)
        {
            // fill relevant free idx
            coords[free_ax[numFreeIdx - 1]] = i;

            // go DEEPER and note subtracted 1 from numFreeIdx for next loop -> this gets the next freeAxis loop correctly
            Recurse(numFreeIdx - 1, binsEachAxis, coords, localIdx, free_ax);
        }
    }
}

Histogram
Histogram::GetSlice(const std::map<std::string, size_t> &fixedBins_) const
{
    /* Function returns a slice of original N-dimensional histogram. Slices may be N-1 dimensional.
       INPUTS: map, fixedBins_, which specifies the axis name and fixed index defining the slice

       OUTPUT: Histogram, slice, the N-1 dimensional histogram slice.

       Function deals with fixed and "free" indexes to create slice. For each free index, a nested loop
       is performed to obtain every combination of axis indexes. This is achieved using recursive loop function,
       Recurse, defined above. Once every combination is found, these coordinates are transformed to global/flat
       bin index relative to original Histogram. Bin contents are found at that globalBin and added to the
       corresponding slice bin.
    */

    // all the axis names in the initial (pre sliced) histogram
    const std::vector<std::string> allAxisNames = fAxes.GetAxisNames();

    // check the pdf does contain the axes and bins asked for
    for (std::map<std::string, size_t>::const_iterator it = fixedBins_.begin(); it != fixedBins_.end(); it++)
    {
        if (std::find(allAxisNames.begin(), allAxisNames.end(), it->first) == allAxisNames.end())
            throw NotFoundError(Formatter()
                                << "Histogram::GetSlice::Tried "
                                << "to get slice of non existent axis "
                                << it->first << "!");
    }

    // want to support multidimensional slices
    if (fixedBins_.size() >= fNDims)
        throw DimensionError("Histogram::GetSlice", fNDims - 1, fixedBins_.size());

    std::vector<std::string> newAxisNames; // collect names of axis present in slice (the "free" idxs)
    std::vector<size_t> newAxisIndexes;    // collect free axis idxs relative to original Histogram
    std::vector<size_t> sliceIndices;      // collect the fixed axis idxs relative to original Histogram
    std::vector<bool> isFixed;             // for each axis, store boolean defining if fixed or free axis

    // loop compares every axis in Histogram to those specified in fixedBins_
    for (std::vector<std::string>::const_iterator it = allAxisNames.begin(); it != allAxisNames.end(); it++)
    {
        if (fixedBins_.find(*it) == fixedBins_.end())
        {
            // axis does not exist in fixedBins, thus it is a free idx to slice out
            newAxisNames.push_back(*it);
            newAxisIndexes.push_back(fAxes.GetAxisIndex(*it));
            std::cout << "Free Index: " << *it << std::endl;
            isFixed.push_back(false);
        }
        else
        {
            // axis exists in fixedBins, thus it is a fixed idx
            sliceIndices.push_back(fAxes.GetAxisIndex(*it));
            isFixed.push_back(true);
            std::cout << "Fixed Indices recorded in Axis " << *it << " with idx " << fixedBins_.at(*it) << std::endl;
        }
    }

    // CREATING NEW SLICE HERE
    AxisCollection newAxes;
    std::vector<size_t> binsEachAxis;          // vector stores number of bins in each free axis
    std::vector<std::vector<size_t>> localIdx; // vector of vectors to store every combination of free and fixed idxs
    std::vector<size_t> coords(fNDims);        // vector stores latest combination of idxs

    // creating the slice and populating binsEachAxis
    for (size_t i = 0; i < newAxisNames.size(); i++)
    {
        newAxes.AddAxis(fAxes.GetAxis(newAxisIndexes[i]));

        // get num bins in each axis
        BinAxis individualAxis = newAxes.GetAxis(newAxisNames[i]);
        size_t bins = individualAxis.GetNBins();
        binsEachAxis.push_back(bins);
    }

    Histogram slice(newAxes);

    // fill in the current combination vector with fixed bin values
    for (size_t i = 0; i < fNDims; i++)
    {
        if (isFixed.at(i) == true)
        {
            coords[i] = fixedBins_.at(allAxisNames[i]);
        }
    }

    // recursively fill localIdxs vector with N nested for loops, where N is number of free axis
    Recurse(newAxisNames.size(), binsEachAxis, coords, localIdx, newAxisIndexes);
    std::cout << "Number of combinations found: " << localIdx.size() << std::endl;

    // for every combination, obtain the global flat idx in original histogram
    size_t bin = 0;
    for (size_t idx = 0; idx < localIdx.size(); idx++)
    {
        size_t oldGlobalBin = fAxes.FlattenIndices(localIdx.at(idx));
        slice.AddBinContent(bin, fBinContents.at(oldGlobalBin)); // fill the relevant slice bin with contents
        bin += 1;
    }

    return slice;
}

double
Histogram::GetBinLowEdge(size_t bin_, size_t index_) const
{
    return fAxes.GetBinLowEdge(bin_, index_);
}

double
Histogram::GetBinHighEdge(size_t bin_, size_t index_) const
{
    return fAxes.GetBinHighEdge(bin_, index_);
}

double
Histogram::GetBinCentre(size_t bin_, size_t index_) const
{
    return fAxes.GetBinCentre(bin_, index_);
}

void Histogram::Add(const Histogram &other_, double weight_)
{
    if (other_.GetAxes() != GetAxes())
        throw ValueError(Formatter() << "Histogram::Add can't add histograms with different binning definitions!");

    for (size_t i = 0; i < GetNBins(); i++)
        AddBinContent(i, other_.GetBinContent(i) * weight_);
}

void Histogram::Multiply(const Histogram &other_)
{
    if (other_.GetAxes() != GetAxes())
        throw ValueError(Formatter() << "Histogram::Add can't add histograms with different binning definitions!");

    for (size_t i = 0; i < GetNBins(); i++)
        SetBinContent(i, GetBinContent(i) * other_.GetBinContent(i));
}

void Histogram::Divide(const Histogram &other_)
{
    if (other_.GetAxes() != GetAxes())
        throw ValueError(Formatter() << "Histogram::Add can't add histograms with different binning definitions!");

    for (size_t i = 0; i < GetNBins(); i++)
        SetBinContent(i, GetBinContent(i) / other_.GetBinContent(i));
}

std::vector<std::string>
Histogram::GetAxisNames() const
{
    return fAxes.GetAxisNames();
}

size_t
Histogram::GetAxisIndex(const std::string &name_) const
{
    return fAxes.GetAxisIndex(name_);
}

void Histogram::AddPadding(double padding_)
{
    std::vector<double> newBinContents;
    for (size_t i = 0; i < fBinContents.size(); i++)
    {
        if (fBinContents.at(i) == 0)
            newBinContents.push_back(padding_);
        else
            newBinContents.push_back(fBinContents[i]);
    }
    fBinContents = newBinContents;
}
