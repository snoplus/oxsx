#include <FitResult.h>
#include <Histogram.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <Exceptions.h>
#include <Formatter.hpp>

FitResult::~FitResult()
{
    delete fStatSpace;
}

void FitResult::SetBestFit(const ParameterDict &bestFit_)
{
    fBestFit = bestFit_;
}

const ParameterDict &
FitResult::GetBestFit() const
{
    return fBestFit;
}

void FitResult::SetStatSpace(const Histogram &statSpace_)
{
    fStatSpace = new Histogram(statSpace_);
}

const Histogram &
FitResult::GetStatSpace() const
{
    return *fStatSpace;
}

FitResult
FitResult::operator=(const FitResult &other_)
{
    if (!fStatSpace)
        fStatSpace = NULL;
    else
        fStatSpace = new Histogram(*other_.fStatSpace);
    fIsValid = other_.fIsValid;
    fBestFit = other_.fBestFit;
    fExtremeVal = other_.fExtremeVal;
    fStatSample = other_.fStatSample;
    f1DProjections = other_.f1DProjections;
    f2DProjections = other_.f2DProjections;
    fCovarianceMatrix = other_.fCovarianceMatrix;
    return *this;
}

FitResult::FitResult(const FitResult &other_)
{
    if (!other_.fStatSpace)
        fStatSpace = NULL;

    else
        fStatSpace = new Histogram(*other_.fStatSpace);
    fStatSample = other_.fStatSample;
    fBestFit = other_.fBestFit;
    fIsValid = other_.fIsValid;
    fExtremeVal = other_.fExtremeVal;
    f1DProjections = other_.f1DProjections;
    f2DProjections = other_.f2DProjections;
    fCovarianceMatrix = other_.fCovarianceMatrix;
}

const std::vector<std::vector<double>> &
FitResult::GetStatSample() const
{
    return fStatSample;
}

void FitResult::SetStatSample(const std::vector<std::vector<double>> &statSample_)
{
    fStatSample = statSample_;
}

bool FitResult::GetValid() const
{
    return fIsValid;
}

void FitResult::SetValid(bool b_)
{
    fIsValid = b_;
}

void FitResult::Print() const
{
    std::cout << AsString() << std::endl;
}

void FitResult::SaveAs(const std::string &fileName_) const
{
    std::ofstream fs;
    fs.open(fileName_.c_str());
    fs << AsString();
    fs.close();
}

std::string
FitResult::AsString() const
{
    time_t t = time(0);
    struct tm *now = localtime(&t);

    std::stringstream ss;
    ss << "-----------------------------------"
       << "-----------------------------------" << std::endl
       << "OXSX Fit Result: "
       << now->tm_year + 1900 << '-'
       << now->tm_mon + 1 << '-'
       << now->tm_mday
       << "  "
       << now->tm_hour
       << ":"
       << now->tm_min
       << ":"
       << now->tm_sec
       << std::endl
       << "-----------------------------------"
       << "-----------------------------------" << std::endl
       << std::endl;

    ss << "Best Fit Values: " << std::endl
       << std::endl;
    for (ParameterDict::const_iterator it = fBestFit.begin(); it != fBestFit.end(); ++it)
    {
        ss << std::setw(25)
           << it->first << "\t\t"
           << std::setw(10)
           << std::setprecision(printPrecision)
           << it->second
           << std::endl;
    }
    return ss.str();
}

void FitResult::SetPrintPrecision(const size_t &prec_)
{
    printPrecision = prec_;
}

size_t
FitResult::GetPrintPrecision() const
{
    return printPrecision;
}
void FitResult::Set1DProjections(const HistMap &hists_)
{
    f1DProjections = hists_;
}

void FitResult::Set2DProjections(const HistMap &hists_)
{
    f2DProjections = hists_;
}

HistMap
FitResult::Get1DProjections() const
{
    return f1DProjections;
}

HistMap
FitResult::Get2DProjections() const
{
    return f2DProjections;
}

double
FitResult::GetExtremeVal() const
{
    return fExtremeVal;
}

void FitResult::SetExtremeVal(double val_)
{
    fExtremeVal = val_;
}

void FitResult::SetCovarianceMatrix(const DenseMatrix covarianceMatrix_)
{
    fCovarianceMatrix = covarianceMatrix_;
}

const DenseMatrix
FitResult::GetCovarianceMatrix() const
{
    return fCovarianceMatrix;
}
