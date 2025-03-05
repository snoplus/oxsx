#include <BinnedEDManager.h>
#include <SystematicManager.h>
#include <BinnedEDShrinker.h>
#include <BinnedED.h>
#include <Exceptions.h>
#include <sstream>

unsigned
BinnedEDManager::GetNPdfs() const
{
    return fOriginalPdfs.size();
}

size_t
BinnedEDManager::GetNDims() const
{
    return fNDims;
}

double
BinnedEDManager::Probability(const Event &data_)
{
    ReassertNorms(true);
    double sum = 0;
    size_t j = 0;
    for (size_t i = 0; i < fWorkingPdfs.size(); i++)
    {
        if (fAllowNormsFittable.at(i) == DIRECT)
        {
            sum += fNormalisations.at(i) * fWorkingPdfs[i].Probability(data_);
            j++;
        }
        else if (fAllowNormsFittable.at(i) == INDIRECT)
        {
            sum += fFittableNorms.at(j) * fWorkingPdfs.at(i).Probability(data_);
            j++;
        }
        else
        {
            // In case where we don't auto-normalise the pdf,
            // the normalisation is held within the pdf itself!
            sum += fWorkingPdfs.at(i).Probability(data_);
        }
    }

    return sum;
}

double
BinnedEDManager::BinProbability(size_t bin_)
{
    ReassertNorms(true);
    double sum = 0;
    size_t j = 0;
    try
    {
        for (size_t i = 0; i < fWorkingPdfs.size(); i++)
        {
            if (fAllowNormsFittable.at(i) == DIRECT)
            {
                sum += fNormalisations.at(i) * fWorkingPdfs.at(i).GetBinContent(bin_);
                j++;
            }
            else if (fAllowNormsFittable.at(i) == INDIRECT)
            {
                sum += fFittableNorms.at(j) * fWorkingPdfs.at(i).GetBinContent(bin_);
                j++;
            }
            else
            {
                // In case where we don't auto-normalise the pdf,
                // the normalisation is held within the pdf itself!
                sum += fWorkingPdfs.at(i).GetBinContent(bin_);
            }
        }
    }
    catch (const std::out_of_range &)
    {
        throw LogicError("BinnedEDManager:: Normalisation vector doesn't match pdf vector - are the normalisations set?");
    }
    return sum;
}

void BinnedEDManager::SetNormalisations(const std::vector<double> &normalisations_)
{
    if (normalisations_.size() != fOriginalPdfs.size())
        throw LogicError("BinnedEDManager: number of norms doesn't match #pdfs");
    fNormalisations = normalisations_;
}

void BinnedEDManager::ApplySystematics(const SystematicManager &sysMan_)
{
    // If there are no systematics dont do anything
    //  ( working pdfs = original pdfs from initialisation)

    if (!sysMan_.GetNSystematics())
        return;
    if (fAllNormsDirFittable)
    {
        sysMan_.DistortEDs(fOriginalPdfs, fWorkingPdfs);
    }
    else
    {
        std::vector<double> norms(fNormalisations.size(), 0);
        sysMan_.DistortEDs(fOriginalPdfs, fWorkingPdfs, &norms);
        for (size_t i = 0; i < norms.size(); i++)
        {
            if (fAllowNormsFittable.at(i) == FALSE)
            {
                fNormalisations[i] = norms.at(i);
            }
            else if (fAllowNormsFittable.at(i) == INDIRECT)
            {
                fNormalisations[i] *= norms.at(i);
            }
        }
    }
}

const BinnedED &
BinnedEDManager::GetOriginalPdf(size_t index_) const
{
    return fOriginalPdfs.at(index_);
}

void BinnedEDManager::AddPdf(const BinnedED &pdf_, const NormFittingStatus norm_fitting_status)
{
    fOriginalPdfs.push_back(pdf_);
    fWorkingPdfs.push_back(pdf_);
    fNPdfs++;
    fNormalisations.resize(fOriginalPdfs.size(), 0);
    fAllowNormsFittable.push_back(norm_fitting_status);
    if (norm_fitting_status)
    {
        fFittableNorms.push_back(0);
    }
    if (norm_fitting_status != DIRECT)
    {
        fAllNormsDirFittable = false;
    }
    RegisterParameters();
}

void BinnedEDManager::AddPdfs(const std::vector<BinnedED> &pdfs_,
                              const std::vector<NormFittingStatus> *norm_fitting_statuses)
{
    if (norm_fitting_statuses != nullptr && pdfs_.size() != norm_fitting_statuses->size())
    {
        throw DimensionError("BinnedEDManager: number of norm_fittable bools doesn't the number of pdfs");
    }

    for (size_t i = 0; i < pdfs_.size(); i++)
    {
        if (norm_fitting_statuses != nullptr)
        {
            AddPdf(pdfs_.at(i), norm_fitting_statuses->at(i));
        }
        else
        {
            AddPdf(pdfs_.at(i));
        }
    }
    RegisterParameters();
}

const std::vector<double> &
BinnedEDManager::GetNormalisations() const
{
    return fNormalisations;
}

void BinnedEDManager::ApplyShrink(const BinnedEDShrinker &shrinker_)
{

    // If there's no buffer, there's no need to apply a shrink!
    if (!shrinker_.GetBuffers().size())
        return;

    // If the number of bins in the working PDFs =/= the number of bins in the corresponding original PDF, it's likely the
    // shrinking has already been applied. This probably only happens if a buffer is being used without systematics. So
    // we check if the bins equals the original number minus the shrinking, and if it has we assume the shrinking has already
    // been applied successfully. If it doesn't equal that or the original number of bins, something bad has happened! For these
    // tests we assume the buffer and number of bins are the same for all PDFs so just look at the 0th element
    int bufferBins = 0;
    std::map<std::string, std::pair<unsigned, unsigned>> bufferMap = shrinker_.GetBuffers();
    for (std::map<std::string, std::pair<unsigned, unsigned>>::iterator it = bufferMap.begin(); it != bufferMap.end(); ++it)
    {
        bufferBins += it->second.first;
        bufferBins += it->second.second;
    }
    if (fWorkingPdfs.at(0).GetNBins() == fOriginalPdfs.at(0).GetNBins() - bufferBins)
        return;
    else if (fWorkingPdfs.at(0).GetNBins() != fOriginalPdfs.at(0).GetNBins() )
        throw DimensionError(Formatter() << "BinnedEDManager:: #bins in Original PDF != #bins in Working PDF");

    for (size_t i = 0; i < fWorkingPdfs.size(); i++)
    {

        const double integral_before = fWorkingPdfs[i].Integral();
        fWorkingPdfs[i] = shrinker_.ShrinkDist(fWorkingPdfs.at(i));
        const double integral_after = fWorkingPdfs[i].Integral();

        // Normalise if normalisation is a fittable param, but if indirect then track any change
        if (fAllowNormsFittable.at(i) == DIRECT)
        {
            fWorkingPdfs[i].Normalise();
        }
        else if (fAllowNormsFittable.at(i) == FALSE)
        {
            fNormalisations[i] = integral_after;
        }
        else
        {
            if (integral_before == 0. && integral_after == 0.)
            {
                fNormalisations[i] = 0.;
            }
            else
            {
                fNormalisations[i] *= integral_after / integral_before;
            }
        }
    }
}

////////////////////////////////
// Make this a fit component! //
////////////////////////////////

std::string
BinnedEDManager::GetName() const
{
    return fName;
}

void BinnedEDManager::SetName(const std::string &n_)
{
    fName = n_;
}

void BinnedEDManager::RenameParameter(const std::string &old_, const std::string &new_)
{
    fParameterManager.RenameParameter(old_, new_);
}

void BinnedEDManager::SetParameter(const std::string &name_, double value_)
{
    fParameterManager.SetParameter(name_, value_);
}

double
BinnedEDManager::GetParameter(const std::string &name_) const
{
    return fParameterManager.GetParameter(name_);
}

void BinnedEDManager::SetParameters(const ParameterDict &ps_)
{
    fParameterManager.SetParameters(ps_);
}

ParameterDict
BinnedEDManager::GetParameters() const
{
    return fParameterManager.GetParameters();
}

size_t
BinnedEDManager::GetParameterCount() const
{
    return fParameterManager.GetParameterCount();
}

std::set<std::string>
BinnedEDManager::GetParameterNames() const
{
    return fParameterManager.GetParameterNames();
}

void BinnedEDManager::RegisterParameters()
{
    fParameterManager.Clear();
    std::vector<std::string> parameterNames;
    if (fAllNormsDirFittable)
    {
        for (size_t i = 0; i < fOriginalPdfs.size(); i++)
        {
            parameterNames.push_back(fOriginalPdfs.at(i).GetName());
        }
        fParameterManager.AddContainer(fNormalisations, parameterNames);
    }
    else
    {
        for (size_t i = 0; i < fOriginalPdfs.size(); i++)
        {
            if (fAllowNormsFittable.at(i))
            {
                parameterNames.push_back(fOriginalPdfs.at(i).GetName());
            }
        }
        fParameterManager.AddContainer(fFittableNorms, parameterNames);
    }
}

void BinnedEDManager::ReassertNorms(bool calcing_binprob)
{
    /*
     * Normalisations that are fittable might have been changed
     * in fFittableNorms, but not in fNormalisations!
     * This method corrects any changes that may have occured,
     * so that both vectors are consistent.
     */
    if (!fAllNormsDirFittable)
    {
        size_t j = 0;
        for (size_t i = 0; i < fNormalisations.size(); i++)
        {
            if (fAllowNormsFittable.at(i) == ((calcing_binprob) ? DIRECT : INDIRECT))
            {
                fNormalisations[i] = fFittableNorms.at(j);
                j++;
            }
        }
    }
}

void BinnedEDManager::AssertDimensions(const std::vector<std::string> &observables)
{
    /*
     * Check whether each PDF in the manager actually has dimensions
     * matching that of the input vector.
     * Whenever this is not the case, the given PDF is marginalised onto
     * the input vector's axes.
     */
    for (auto &working_pdf : fWorkingPdfs)
    {
        if (working_pdf.GetObservables() != observables)
        {
            working_pdf = working_pdf.Marginalise(observables);
        }
    }
}
