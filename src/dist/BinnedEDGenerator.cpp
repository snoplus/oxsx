#include <BinnedEDGenerator.h>
#include <Exceptions.h>
#include <stdlib.h>
#include <Rand.h>
#include <cmath>
#include <iostream>
#include <algorithm>

BinnedED
BinnedEDGenerator::ExpectedRatesED(const std::string &name) const
{
    if (!(fPdfs.size() == fRates.size()))
        throw LogicError("BinnedEDGenerator::Need exactly one rate for each ");

    if (!fPdfs.size())
        throw LogicError("BinnedEDGenerator::No source pdfs!!");

    BinnedED fakePdf(fPdfs.at(0));

    if (name == "")
    {
        std::string fname = "Gen_Exp";
        for (size_t i = 0; i < fPdfs.size(); i++)
        {
            fname += fPdfs.at(i).GetName();
            if (i != fPdfs.size() - 1)
                fname += "_";
        }
        fakePdf.SetName(fname);
    }
    else
    {
        fakePdf.SetName(name);
    }

    fakePdf.Empty();
    for (size_t i = 0; i < fPdfs.size(); i++)
    {
        unsigned counts = round(fRates.at(i));
        for (unsigned _ = 0; _ < counts; _++)
            fakePdf.AddBinContent(RandomBin(i), 1);
    }

    return fakePdf;
}

BinnedED
BinnedEDGenerator::PoissonFluctuatedED(const std::string &name) const
{
    if (!(fPdfs.size() == fRates.size()))
        throw LogicError("BinnedEDGenerator::Need exactly one rate for each ");

    if (!fPdfs.size())
        throw LogicError("BinnedEDGenerator::No source pdfs!!");

    BinnedED fakePdf(fPdfs.at(0));

    if (name == "")
    {
        std::string fname = "Gen_Poi";
        for (size_t i = 0; i < fPdfs.size(); i++)
        {
            fname += fPdfs.at(i).GetName();
            if (i != fPdfs.size() - 1)
                fname += "_";
        }
        fakePdf.SetName(fname);
    }
    else
    {
        fakePdf.SetName(name);
    }

    fakePdf.Empty();
    for (size_t i = 0; i < fPdfs.size(); i++)
    {
        unsigned counts = Rand::Poisson(fRates.at(i));
        for (unsigned _ = 0; _ < counts; _++)
            fakePdf.AddBinContent(RandomBin(i), 1);
    }

    return fakePdf;
}

size_t
BinnedEDGenerator::RandomBin(size_t pdfIndex_) const
{
    if (!fPdfs.at(pdfIndex_).GetNDims())
        throw LogicError("BinnedEDGenerator::Pdf of dimension 0!");

    // construct cdf
    std::vector<double> cdf(fPdfs.at(pdfIndex_).GetNBins(), 0);
    double sum = 0;
    for (size_t i = 0; i < cdf.size(); i++)
    {
        sum += fPdfs.at(pdfIndex_).GetBinContent(i);
        cdf[i] = sum;
    }

    for (size_t i = 0; i < cdf.size(); i++)
        cdf[i] /= sum;

    // random val on unit line
    double rnd = Rand::Uniform();
    return std::lower_bound(cdf.begin(), cdf.end(), rnd) - cdf.begin();
}
