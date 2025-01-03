#include <Convolution.h>
#include <PDF.h>
#include <JumpPDF.h>
#include <Exceptions.h>
#include <string>
#include <algorithm>

void Convolution::SetFunction(PDF *function_)
{
    // wrap this up if position independent kernel of the form P(x | x2) = P(x - x2)
    delete fDist;
    fDist = static_cast<ConditionalPDF *>(new JumpPDF("kernel", function_));
}

void Convolution::SetConditionalPDF(ConditionalPDF *c_)
{
    delete fDist;
    fDist = c_->Clone();
}

Convolution::~Convolution()
{
    delete fDist;
}

void Convolution::ConstructSubmatrix(std::vector<long long unsigned int> &column_indices, std::vector<long long unsigned int> &row_indices,
                                     std::vector<double> &vals) const
{
    /*
     * Construct the sub-matrix associated with the convolution in the transformed axes.
     * Returns this information, by reference, in the form of column & row indices for non-zero entries of this submatrix,
     * and their associated values. This form is what Armadillo prefers for sparse matrices!
     */
    // variables storing the axes bin information
    std::vector<double> binCentres(fSubMapAxes.GetNDimensions());
    std::vector<double> lowEdges(fSubMapAxes.GetNDimensions());
    std::vector<double> highEdges(fSubMapAxes.GetNDimensions());
    // Pre-allocate memory for the sub-matrix data
    // likely to need only a fraction of this memory!
    column_indices.reserve(fSubMapAxes.GetNBins() * fSubMapAxes.GetNBins());
    row_indices.reserve(fSubMapAxes.GetNBins() * fSubMapAxes.GetNBins());
    vals.reserve(fSubMapAxes.GetNBins() * fSubMapAxes.GetNBins());
    // Loop over all entries of the sub-matrix to determine their values
    for (long long unsigned int origBin = 0; origBin < fSubMapAxes.GetNBins(); origBin++)
    {
        // get the centre of the bin. Need to offset by this for a convolution
        fSubMapAxes.GetBinCentres(origBin, binCentres);

        // loop over the bins it can be smeared into
        for (long long unsigned int destBin = 0; destBin < fSubMapAxes.GetNBins(); destBin++)
        {
            fSubMapAxes.GetBinLowEdges(destBin, lowEdges);
            fSubMapAxes.GetBinHighEdges(destBin, highEdges);
            const double integral = fDist->Integral(lowEdges, highEdges, binCentres);
            // only bother adding to matrix if non-zero!
            if (integral > 0.)
            {
                column_indices.push_back(origBin);
                row_indices.push_back(destBin);
                vals.push_back(integral);
            }
        }
    }
}

void Convolution::MakeFullMatrix(const std::vector<long long unsigned int> &column_indices, const std::vector<long long unsigned int> &row_indices,
                                 const std::vector<double> &vals, SparseMatrix &response_blocked) const
{
    /*
     * Given the response submatrix, make the full response matrix!
     * The actual element values have already been calculated, we're just trying to put them in the right places now.
     */
    // Set up variables first; pre-allocate memory for the large vectors associated with the full matrix's info
    const size_t N = fAxes.GetNBins();
    const size_t n_sub = fSubMapAxes.GetNBins();
    const size_t n_blocks = N / n_sub;
    const size_t size_block = vals.size();
    std::vector<long long unsigned int> column_indices_bl;
    std::vector<long long unsigned int> row_indices_bl;
    std::vector<double> vals_bl;
    column_indices_bl.reserve(size_block * n_blocks);
    row_indices_bl.reserve(size_block * n_blocks);
    vals_bl.reserve(size_block * n_blocks);
    // Loop over the submatrix "blocks" in the full response matrix
    for (size_t block_idx = 0; block_idx < n_blocks; block_idx++)
    {
        // Just add a copy of the values from vals into vals_bl, using fancy STL insert
        vals_bl.insert(vals_bl.end(), std::begin(vals), std::end(vals));
        // Now loop over the numbers in column_indices and row_indices
        for (size_t idx = 0; idx < size_block; idx++)
        {
            // Use the cached mapping from blocked index -> true bin index for both columns and rows!
            column_indices_bl.push_back(fColumnPerms.at(column_indices.at(idx) + block_idx * n_sub));
            row_indices_bl.push_back(fColumnPerms.at(row_indices.at(idx) + block_idx * n_sub));
        }
    }
    // Fill the matrix object with these values
    // Surprisngly, this ends up typically being the slowest part of the whole Construct() method!!
    // This is because Armadillo has to first sort the vectors we're giving it by column & row indices.
    response_blocked.SetComponents(row_indices_bl, column_indices_bl, vals_bl);
}

void Convolution::Construct()
{
    /*
     * Method that constructs the response matrix associated with this convlution systematic.
     * Attempts to be clever and speedy about building this, because this method can get computationally-
     * expensive, fast, when dealing with BinnedED objects with numerous dimensions and many bins!
     *
     * The key insights are that:
     *   (a) Each subspace of the BinnedED object that is being transformed by this convolution
     *       will have an identical response "sub"-matrix, S, acting only in that subspace.
     *   (b) If you were to index the bins of the BinnedED object along the transformed axes first,
     *       then because of (a) the full response matrix, R, will be of block-diagonal form:
     *
     *       [S      ]
     *       [ S     ]
     *       [  S    ]
     *  R =  [   .   ],   where S is the response sub-matrix for the transforming subspace.
     *       [    .  ]
     *       [     . ]
     *       [      S]
     *
     * Annoyingly, (b) doesn't always line up with reality - often, the indexing scheme doesn't go
     * transforming-axis first. So, we permute the rows and columns of the elements such that this works
     * right.
     */
    // First - got to make sure we have axes and a kernel to convolve with!
    const size_t N = fAxes.GetNBins();
    if (!fDist || !N)
        throw LogicError("Convolution::Construct() : Tried to construct convolution without axes or function/distribution, or both!!");
    // Only work out what the mapping from blocked index to true bin index if not cached already
    if (!fCachedIndexPermutations)
        CacheIndexPermutations();

    // Work out the transition probabilitites within this sub set of the bins
    std::vector<long long unsigned int> column_indices;
    std::vector<long long unsigned int> row_indices;
    std::vector<double> vals;
    ConstructSubmatrix(column_indices, row_indices, vals);
    // Expand to full size matrix
    SparseMatrix response_blocked_sm(N, N);
    MakeFullMatrix(column_indices, row_indices, vals, response_blocked_sm);
    // Actually set response matrix to construction!
    fResponse = response_blocked_sm;
}

AxisCollection Convolution::DetermineAxisSubCollection(const std::vector<size_t> &rel_indices) const
{
    /*
     * Given a vector of axis indices, return an AxisCollection with those axes from fAxes.
     */
    AxisCollection ax;
    for (size_t i = 0; i < rel_indices.size(); i++)
        ax.AddAxis(fAxes.GetAxis(rel_indices.at(i)));

    return ax;
}

size_t Convolution::BlockedBinningIndex(size_t bin_index, const std::vector<size_t> &relativeIndices) const
{
    /* Given a bin index, using the existing binning order,
     * convert to the bin number that would be used if we asserted that
     * we counted in the transformed axes first.
     */
    if (relativeIndices.size() == fAxes.GetNDimensions())
    {
        return bin_index; // case where all axes are being transformed is trivial!
    }
    // First - unpack indices for each axis
    const std::vector<size_t> bin_idxs = fAxes.UnpackIndices(bin_index);
    // Now split these indices into those within the transforing sub-collection, and those which aren't
    std::vector<size_t> bin_idxs_sub;
    std::vector<size_t> bin_idxs_notsub;
    for (size_t ax_idx = 0; ax_idx < fAxes.GetNDimensions(); ax_idx++)
    {
        if (std::find(relativeIndices.begin(), relativeIndices.end(), ax_idx) != std::end(relativeIndices))
        {
            bin_idxs_sub.push_back(bin_idxs.at(ax_idx));
        }
        else
        {
            bin_idxs_notsub.push_back(bin_idxs.at(ax_idx));
        }
    }
    // Re-pack each individually, the combine to get bin index if going through transformed
    // sub-collection first
    const size_t idx_sub = fSubMapAxes.FlattenIndices(bin_idxs_sub);
    const size_t idx_notsub = fNotSubMapAxes.FlattenIndices(bin_idxs_notsub);
    const size_t idx_blocked = idx_sub + idx_notsub * fSubMapAxes.GetNBins();

    return idx_blocked;
}

void Convolution::CacheIndexPermutations()
{
    /*
     * During construction of the response matrix, we first create it in block diagonal form,
     * which would be correct if we count the bins along the transformed axes first.
     * However, if this is not the case, we must permute the elements of the block matrix
     * to get the matrix entry values in the correct place.
     * This method creates a mapping from blocked matrix bin index to true bin index.
     */
    //  get the axes that this systematic will act on
    std::vector<size_t> relativeIndices = fTransObs.GetRelativeIndices(fDistObs);
    fSubMapAxes = DetermineAxisSubCollection(relativeIndices);
    // get the axes that this systematic will /not/ act on
    std::vector<size_t> relativeIndicesNot;
    for (size_t idx = 0; idx < fDistObs.GetNObservables(); idx++)
    {
        if (std::find(relativeIndices.begin(), relativeIndices.end(), idx) == std::end(relativeIndices))
        {
            relativeIndicesNot.push_back(idx);
        }
    }
    fNotSubMapAxes = DetermineAxisSubCollection(relativeIndicesNot);
    // Fill a vector of the mapping from bin index when matrix is in its blocked form,
    // to the index when it is in its default state.
    const size_t N = fAxes.GetNBins();
    fColumnPerms.resize(N);
    for (size_t i = 0; i < N; i++)
    {
        const size_t idx_blocked = BlockedBinningIndex(i, relativeIndices);
        fColumnPerms.at(idx_blocked) = i;
    }

    fCachedIndexPermutations = true;
}

///////////////////////////////
// Make this object fittable //
///////////////////////////////

// Fitting this dist to data means adjusting the underlying function

void Convolution::RenameParameter(const std::string &old_, const std::string &new_)
{
    fDist->RenameParameter(old_, new_);
}

void Convolution::SetParameter(const std::string &name_, double value_)
{
    fDist->SetParameter(name_, value_);
}

double
Convolution::GetParameter(const std::string &name_) const
{
    return fDist->GetParameter(name_);
}

void Convolution::SetParameters(const ParameterDict &ps_)
{
    try
    {
        fDist->SetParameters(ps_);
    }
    catch (const ParameterError &e_)
    {
        throw ParameterError("Convolution internal function: " + std::string(e_.what()));
    }
}

ParameterDict
Convolution::GetParameters() const
{
    return fDist->GetParameters();
}

size_t
Convolution::GetParameterCount() const
{
    return fDist->GetParameterCount();
}

std::set<std::string>
Convolution::GetParameterNames() const
{
    return fDist->GetParameterNames();
}

std::string
Convolution::GetName() const
{
    return fName;
}
void Convolution::SetName(const std::string &n_)
{
    fName = n_;
}
