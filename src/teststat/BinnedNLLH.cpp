#include <BinnedNLLH.h>
#include <math.h>
#include <DataSet.h>
#include <Exceptions.h>
#include <DistFiller.h>
#include <CutLog.h>
#include <Exceptions.h>
#include <Formatter.hpp>
#include <iostream>
#include <Formatter.hpp>

double
BinnedNLLH::Evaluate()
{
    if (!GetDataSet() && !GetCalculatedDataDist())
        throw LogicError("BinnedNNLH function called with no data set and no DataDist! set one of these first");

    if (!GetCalculatedDataDist())
        BinData();
    // Make sure that normalisation params in PDF Manager are setup correctly
    fPdfManager.ReassertNorms();
    if (!fAlreadyShrunk)
    {
        fPdfShrinker.SetBinMap(fDataDist);
        fDataDist = fPdfShrinker.ShrinkDist(fDataDist);
        fAlreadyShrunk = true;
    }

    GetSystematicManager().Construct();
    // Apply systematics
    fPdfManager.ApplySystematics(GetSystematicManager());
    // Marginalise pdfs back to observable dimensions only, if necessary
    fPdfManager.AssertDimensions(fDataDist.GetObservables());
    // Apply Shrinking
    fPdfManager.ApplyShrink(fPdfShrinker);

    const std::vector<double> &normalisations = fPdfManager.GetNormalisations();

    // loop over bins and calculate the likelihood
    if (GetDebugMode())
    {
        std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    }
    double nLogLH = 0;

    for (size_t i = 0; i < fDataDist.GetNBins(); i++)
    {

        double prob = fPdfManager.BinProbability(i);
        if (!prob)
            throw std::runtime_error(Formatter() << "BinnedNLLH::Encountered zero probability bin! #" << i);

        double newProb = prob;
        double penalty = 0;

        if (fUseBarlowBeeston)
        {
            // Calculate MC Statistical uncertainty using Beeston Barlow: sqrt(sum(weights^2))/prob
            // Essentially solving Eqn 11 in https://arxiv.org/abs/1103.0354 for Beta

            // sigma = Sum[root(weight_int)] (sum is over interaction types)
            // weight_int = bincontent_int / generated_int
            // bin_content_int = prob_int_bin * normalisation_int
            // sigma = Sum[root(prob_int_bin * normalisation_int / generated_int_bin)]
            // we care about sq of fractional bin uncertainty (total bin content is prob)
            // sigma2 = (sigma/bin_content_tot)^2
            // sigma2 = (Sum[root(prob_int_bin * normalisation_int / generated_int_bin)] / prob)^2
            // generated_int = GenRate_int * prob_int_bin
            // sigma2 = (Sum[normalisation_int / generated_int] / prob)^2

            double sigma2 = 0;
            for (unsigned int j = 0; j < normalisations.size(); j++)
                sigma2 += (normalisations.at(j)) / ((double)fGenRates.at(j) * prob * prob);

            // Eqn 11 in https://arxiv.org/abs/1103.0354 for Beta: beta^2 + (mu x sigma^2 - 1)beta - n x sigma^2   (prob is mu, the mc events in bin. n is data events in bin)

            // b in quadratic
            double b = prob * sigma2 - 1;
            // b^2-4ac in quadratic
            double det = b * b + 4 * fDataDist.GetBinContent(i) * sigma2;

            if (det < 0.)
            { // imaginary roots
                throw OXSXException(Formatter() << "BinnedNLLH:: Negative determinant in beta calculation for Barlow-Beeston. Something has gone terribly wrong.");
            }

            double sign = (b >= 0) ? 1. : -1.;
            double q = -0.5 * (b + sign * sqrt(det));

            double x1 = q;
            double x2 = fDataDist.GetBinContent(i) * sigma2 / q;

            double beta = (x1 > x2) ? x1 : x2;

            // and update mc bin content
            newProb = beta * prob;
            // get beta penalty term
            penalty = (beta - 1) * (beta - 1) / (2 * sigma2);
        }

        if (GetDebugMode())
        {
            std::cout << "Bin " << i << ", MC bin probability: " << prob << ", data bin probability: ";
            std::cout << fDataDist.GetBinContent(i) << std::endl;
        }

        // LogL = mu_i - data * log(update MC) + Beta Penalty + Norm Correction
        nLogLH = nLogLH - fDataDist.GetBinContent(i) * log(newProb) + penalty + newProb;
    }

    if (GetDebugMode())
    {
        std::cout << "NLLH after summing over bins only: " << nLogLH << std::endl;
    }

    // Constraints
    nLogLH += GetConstraintManager().Evaluate(GetParameters());

    if (GetDebugMode())
    {
        std::cout << "\nTotal NLLH: " << nLogLH << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    }
    return nLogLH;
}

void BinnedNLLH::BinData()
{
    fDataDist = BinnedED(fPdfManager.GetOriginalPdf(0)); // make a copy for same binning and data rep
    fDataDist.Empty();
    CutLog log(GetCuts().GetCutNames());
    DistFiller::FillDist(fDataDist, *GetDataSet(), GetCuts(), log);
    SetCalculatedDataDist(true);
    SetSignalCutLog(log);
}

void BinnedNLLH::SetNormalisations(const std::vector<double> &norms_)
{
    fPdfManager.SetNormalisations(norms_);
}

std::vector<double>
BinnedNLLH::GetNormalisations() const
{
    return fPdfManager.GetNormalisations();
}

void BinnedNLLH::AddPdfs(const std::vector<BinnedED> &pdfs,
                         const std::vector<std::vector<std::string>> &sys_,
                         const std::vector<NormFittingStatus> *norm_fitting_statuses)
{

    if (fUseBarlowBeeston)
        throw OXSXException(Formatter() << "BinnedNLLH:: Must set generated rates if using Barlow-Beeston");

    if (pdfs.size() != sys_.size())
        throw DimensionError(Formatter() << "BinnedNLLH:: #sys_ != #group_");
    if (norm_fitting_statuses != nullptr && pdfs.size() != norm_fitting_statuses->size())
    {
        throw DimensionError("BinnedNLLH: number of norm_fittable bools doesn't the number of pdfs");
    }
    for (size_t i = 0; i < pdfs.size(); ++i)
    {
        if (norm_fitting_statuses == nullptr)
        {
            AddPdf(pdfs.at(i), sys_.at(i));
        }
        else
        {
            AddPdf(pdfs.at(i), sys_.at(i), norm_fitting_statuses->at(i));
        }
    }
}

void BinnedNLLH::AddPdfs(const std::vector<BinnedED> &pdfs,
                         const std::vector<int> &genrates_,
                         const std::vector<NormFittingStatus> *norm_fitting_statuses)
{

    if (norm_fitting_statuses != nullptr && pdfs.size() != norm_fitting_statuses->size())
    {
        throw DimensionError("BinnedNLLH: number of norm_fittable bools doesn't the number of pdfs");
    }
    for (size_t i = 0; i < pdfs.size(); i++)
    {
        if (norm_fitting_statuses != nullptr)
        {
            AddPdf(pdfs.at(i), genrates_.at(i), norm_fitting_statuses->at(i));
        }
        else
        {
            AddPdf(pdfs.at(i), genrates_.at(i));
        }
    }
}

void BinnedNLLH::AddPdfs(const std::vector<BinnedED> &pdfs,
                         const std::vector<NormFittingStatus> *norm_fitting_statuses)
{
    if (fUseBarlowBeeston)
        throw OXSXException(Formatter() << "BinnedNLLH:: Must set generated rates if using Barlow-Beeston");
    if (norm_fitting_statuses != nullptr && pdfs.size() != norm_fitting_statuses->size())
    {
        throw DimensionError("BinnedNLLH: number of norm_fittable bools doesn't the number of pdfs");
    }
    for (size_t i = 0; i < pdfs.size(); i++)
    {
        if (norm_fitting_statuses != nullptr)
        {
            AddPdf(pdfs.at(i), norm_fitting_statuses->at(i));
        }
        else
        {
            AddPdf(pdfs.at(i));
        }
    }
}

void BinnedNLLH::AddPdfs(const std::vector<BinnedED> &pdfs,
                         const std::vector<std::vector<std::string>> &sys_,
                         const std::vector<int> &genrates_,
                         const std::vector<NormFittingStatus> *norm_fitting_statuses)
{
    if (pdfs.size() != sys_.size())
        throw DimensionError(Formatter() << "BinnedNLLH:: #sys_ != #group_");
    for (size_t i = 0; i < pdfs.size(); i++)
    {
        if (norm_fitting_statuses != nullptr)
        {
            AddPdf(pdfs.at(i), sys_.at(i), genrates_.at(i), norm_fitting_statuses->at(i));
        }
        else
        {
            AddPdf(pdfs.at(i), sys_.at(i), genrates_.at(i));
        }
    }
}

void BinnedNLLH::AddPdf(const BinnedED &pdf_,
                        const std::vector<std::string> &syss_,
                        const NormFittingStatus norm_fitting_status)
{
    if (fUseBarlowBeeston)
        throw OXSXException(Formatter() << "BinnedNLLH:: Must set generated rates if using Barlow-Beeston");
    fPdfManager.AddPdf(pdf_, norm_fitting_status);
    GetSystematicManager().AddDist(pdf_, syss_);
}

void BinnedNLLH::AddPdf(const BinnedED &pdf_,
                        const int &genrate_,
                        const NormFittingStatus norm_fitting_status)
{
    fPdfManager.AddPdf(pdf_, norm_fitting_status);
    GetSystematicManager().AddDist(pdf_, "");
    fGenRates.push_back(genrate_);
}

void BinnedNLLH::AddPdf(const BinnedED &pdf_,
                        const NormFittingStatus norm_fitting_status)
{
    if (fUseBarlowBeeston)
        throw OXSXException(Formatter() << "BinnedNLLH:: Must set generated rates if using Barlow-Beeston");
    fPdfManager.AddPdf(pdf_, norm_fitting_status);
    GetSystematicManager().AddDist(pdf_, "");
}

void BinnedNLLH::AddPdf(const BinnedED &pdf_,
                        const std::vector<std::string> &syss_,
                        const int &genrate_,
                        const NormFittingStatus norm_fitting_status)
{
    fPdfManager.AddPdf(pdf_, norm_fitting_status);
    GetSystematicManager().AddDist(pdf_, syss_);
    fGenRates.push_back(genrate_);
}

void BinnedNLLH::SetDataDist(const BinnedED &binnedPdf_)
{
    fDataDist = binnedPdf_;
    SetCalculatedDataDist(true);
}

BinnedED
BinnedNLLH::GetDataDist() const
{
    return fDataDist;
}

void BinnedNLLH::SetBarlowBeeston(const bool useBarlowBeeston_)
{
    fUseBarlowBeeston = useBarlowBeeston_;
}

void BinnedNLLH::SetBuffer(const std::string &dim_, unsigned lower_, unsigned upper_)
{
    fPdfShrinker.SetBuffer(dim_, lower_, upper_);
}

std::pair<unsigned, unsigned>
BinnedNLLH::GetBuffer(const std::string &dim_) const
{
    return fPdfShrinker.GetBuffer(dim_);
}

void BinnedNLLH::SetBufferAsOverflow(bool b_)
{
    fPdfShrinker.SetUsingOverflows(b_);
}

bool BinnedNLLH::GetBufferAsOverflow() const
{
    return fPdfShrinker.GetUsingOverflows();
}

/////////////////////////////////////////////////////////
// Declare which objects should be adjusted by the fit //
/////////////////////////////////////////////////////////
void BinnedNLLH::RegisterFitComponents()
{
    GetComponentManager().Clear();
    GetComponentManager().AddComponent(&fPdfManager);

    // Because the limits are set by name they can be added in any order.
    const std::map<std::string, std::vector<Systematic *>> sys_ = GetSystematicManager().GetSystematicsGroup();
    std::vector<std::string> alreadyAdded;
    for (const auto &group_ : sys_)
    {
        for (const auto &item : group_.second)
        {
            if (std::find(alreadyAdded.begin(), alreadyAdded.end(), item->GetName()) == alreadyAdded.end())
            {
                GetComponentManager().AddComponent(item);
                alreadyAdded.push_back(item->GetName());
            }
        } // End of group
    } // End of groups
}

void BinnedNLLH::SetParameters(const ParameterDict &params_)
{
    try
    {
        GetComponentManager().SetParameters(params_);
    }
    catch (const ParameterError &e_)
    {
        throw ParameterError(std::string("BinnedNLLH::") + e_.what());
    }
}

ParameterDict
BinnedNLLH::GetParameters() const
{
    return GetComponentManager().GetParameters();
}

size_t
BinnedNLLH::GetParameterCount() const
{
    return GetComponentManager().GetTotalParameterCount();
}

std::set<std::string>
BinnedNLLH::GetParameterNames() const
{
    return GetComponentManager().GetParameterNames();
}

void BinnedNLLH::SetPdfManager(const BinnedEDManager &man_)
{
    fPdfManager = man_;
}