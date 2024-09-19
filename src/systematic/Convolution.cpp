#include <Convolution.h>
#include <PDF.h>
#include <JumpPDF.h>
#include <DenseMatrix.h>
#include <Exceptions.h>
#include <string>
#include <armadillo>
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

void Convolution::ConstructSubmatrix(std::vector<long long unsigned int>& column_indices, std::vector<long long unsigned int>& row_indices,
                                     std::vector<double>& vals) const
{
    std::vector<double> binCentres(fSubMapAxes.GetNDimensions());
    std::vector<double> lowEdges(fSubMapAxes.GetNDimensions());
    std::vector<double> highEdges(fSubMapAxes.GetNDimensions());

    column_indices.reserve(fSubMapAxes.GetNBins()*fSubMapAxes.GetNBins()); // likely to need only a fraction of this memory!
    row_indices.reserve(fSubMapAxes.GetNBins()*fSubMapAxes.GetNBins());
    vals.reserve(fSubMapAxes.GetNBins()*fSubMapAxes.GetNBins());

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
            if (integral > 1.0e-20) { // only bother adding to matrix if non-zero!
                column_indices.push_back(origBin);
                row_indices.push_back(destBin);
                vals.push_back(integral);
            }
            // subMap.SetComponent(destBin, origBin, fDist->Integral(lowEdges, highEdges, binCentres));
        }
    }
}


void Convolution::MakeBlockMatrix(const std::vector<long long unsigned int>& column_indices, const std::vector<long long unsigned int>& row_indices,
                                  const std::vector<double>& vals, SparseMatrix& response_blocked) {
    const size_t N = fAxes.GetNBins();
    const size_t n_sub = fSubMapAxes.GetNBins();
    const size_t n_blocks = N/n_sub;
    const size_t size_block = vals.size();
    std::vector<long long unsigned int> column_indices_bl;
    std::vector<long long unsigned int> row_indices_bl;
    std::vector<double> vals_bl;
    column_indices_bl.reserve(size_block*n_blocks);
    row_indices_bl.reserve(size_block*n_blocks);
    vals_bl.reserve(size_block*n_blocks);
    for (size_t block_idx=0; block_idx<n_blocks; block_idx++) {
        vals_bl.insert(vals_bl.end(), std::begin(vals), std::end(vals));
        for (size_t idx=0; idx<size_block; idx++) {
            column_indices_bl.push_back(fColumnPerms.at(column_indices.at(idx)+block_idx*n_sub));
            row_indices_bl.push_back(fRowPerms.at(row_indices.at(idx)+block_idx*n_sub));
        }
    }
    response_blocked.SetComponents(row_indices_bl, column_indices_bl, vals_bl);
}


void Convolution::Construct()
{
    const size_t N = fAxes.GetNBins();
    if (!fDist || !N)
        throw LogicError("Convolution::Construct() : Tried to construct convolution without axes or function/distribution, or both!!");

    if (!fCachedPermutationMatrix)
        CachePermutationMatrix();

    // Work out the transition probabilitites within this sub set of the bins
    std::vector<long long unsigned int> column_indices;
    std::vector<long long unsigned int> row_indices;
    std::vector<double> vals;
    ConstructSubmatrix(column_indices, row_indices, vals);
    // Expand to full size matrix by first generating a block diagonal matrix
    SparseMatrix response_blocked_sm(N, N);
    MakeBlockMatrix(column_indices, row_indices, vals, response_blocked_sm);

    // arma::sp_mat response_blocked(N, N);
    // const size_t n_blocks = N/n_sub;

    // for (size_t block_idx=0; block_idx<n_blocks; block_idx++) {
    //     response_blocked.submat(block_idx*n_sub,block_idx*n_sub, (block_idx+1)*n_sub-1,(block_idx+1)*n_sub-1) = subMap.GetMatrix();
    // }
    // SparseMatrix response_blocked_sm(N, N);
    // response_blocked_sm.SetMatrix(response_blocked);
    // Pre- and post-multiply by permutation matrices to go from and to default PDF indexing -> transformation axis-first indexing,
    // So that block diagonal matrix can be used nicely
    fResponse = response_blocked_sm;
    // if (fPermMatrixIdentity) {
    //     fResponse = response_blocked_sm;
    // } else {
    //     fResponse = fBinningPermT * response_blocked_sm * fBinningPerm;
    // }

    // // Now expand to the full size matrix. Elements are zero by default
    // // compatible bins are cached, values must match the smaller matrix above
    // size_t destBin = -1;
    // std::vector<long long unsigned int> nonZeroRowIndices;
    // std::vector<long long unsigned int> nonZeroColIndices;
    // std::vector<double> values;
    // nonZeroRowIndices.reserve(fCompatibleBins.at(0).size());
    // nonZeroColIndices.reserve(fCompatibleBins.at(0).size());

    // for (size_t origBin = 0; origBin < fAxes.GetNBins(); origBin++)
    // {
    //     for (size_t i = 0; i < fCompatibleBins.at(origBin).size(); i++)
    //     {
    //         destBin = fCompatibleBins.at(origBin).at(i);
    //         nonZeroRowIndices.push_back(origBin);
    //         nonZeroColIndices.push_back(destBin);
    //         values.push_back(subMap.GetComponent(fSysBins.at(origBin),
    //                                              fSysBins.at(destBin)));
    //     }
    // }

    // fResponse.SetComponents(nonZeroRowIndices, nonZeroColIndices, values);
}

AxisCollection Convolution::DetermineAxisSubCollection(const std::vector<size_t>& rel_indices)
{
    AxisCollection ax;
    for (size_t i = 0; i < rel_indices.size(); i++)
        ax.AddAxis(fAxes.GetAxis(rel_indices.at(i)));

    return ax;
}


size_t Convolution::BlockedBinningIndex(size_t bin_index, const std::vector<size_t>& relativeIndices)
{
    /* Given a bin index, using the existing binning order,
     * convert to the bin number that would be used if we asserted that
     * we counted in the transformed axes first.
     */
    // First - unpack indices for each axis
    const std::vector<size_t> bin_idxs = fAxes.UnpackIndices(bin_index);
    // Now split these indices into those within the transforing sub-collection, and those which aren't
    std::vector<size_t> bin_idxs_sub;
    std::vector<size_t> bin_idxs_notsub;
    for (size_t ax_idx=0; ax_idx<fAxes.GetNDimensions(); ax_idx++) {
        if (std::find(relativeIndices.begin(), relativeIndices.end(), ax_idx) != std::end(relativeIndices)) {
            bin_idxs_sub.push_back(bin_idxs.at(ax_idx));
        } else {
            bin_idxs_notsub.push_back(bin_idxs.at(ax_idx));
        }
    }
    // Re-pack each individually, the combine to get bin index if going through transformed
    // sub-collection first
    const size_t idx_sub = fSubMapAxes.FlattenIndices(bin_idxs_sub);
    const size_t idx_notsub = fNotSubMapAxes.FlattenIndices(bin_idxs_notsub);
    const size_t idx_blocked = idx_sub + idx_notsub*fSubMapAxes.GetNBins();

    return idx_blocked;
}


void Convolution::CachePermutationMatrix()
{
    //  get the axes that this systematic will act on
    std::vector<size_t> relativeIndices = fTransObs.GetRelativeIndices(fDistObs);
    fSubMapAxes = DetermineAxisSubCollection(relativeIndices);
    // get the axes that this systematic will /not/ act on
    std::vector<size_t> relativeIndicesNot;
    for (size_t idx=0; idx<fDistObs.GetNObservables(); idx++) {
        if (std::find(relativeIndices.begin(), relativeIndices.end(), idx) == std::end(relativeIndices)) {
            relativeIndicesNot.push_back(idx);
        }
    }
    fNotSubMapAxes = DetermineAxisSubCollection(relativeIndicesNot);
    //
    const size_t N = fAxes.GetNBins();
    fColumnPerms.resize(N);
    fRowPerms.reserve(N);
    for (size_t i = 0; i < N; i++) {
        const size_t idx_blocked = BlockedBinningIndex(i, relativeIndices);
        fColumnPerms.at(idx_blocked) = i;
        fRowPerms.push_back(idx_blocked);
    }

    //
    // std::vector<long long unsigned int> axesBinningIndices;
    // std::vector<long long unsigned int> blockedBinningIndices;
    // std::vector<double> values(fAxes.GetNBins(), 1.0);
    // fPermMatrixIdentity = true;
    // for (size_t i = 0; i < fAxes.GetNBins(); i++) {
    //     axesBinningIndices.push_back(i);
    //     const size_t idx_blocked = BlockedBinningIndex(i, relativeIndices);
    //     if (idx_blocked != i) { fPermMatrixIdentity = false; }
    //     blockedBinningIndices.push_back(idx_blocked);
    // }
    // // std::cout << blockedBinningIndices.size() << " " << axesBinningIndices.size() << " " << values.size() << std::endl;
    // if (!fPermMatrixIdentity) {
    //     fBinningPerm.SetComponents(blockedBinningIndices, axesBinningIndices, values);
    //     fBinningPermT.SetComponents(axesBinningIndices, blockedBinningIndices, values);
    // }
    // // std::cout << "fPermMatrixIdentity: " << fPermMatrixIdentity << std::endl;
    // fCachedPermutationMatrix = true;
}

// void Convolution::CacheCompatibleBins()
// {
//     fCompatibleBins.resize(fAxes.GetNBins());
//     // only need to look at one side of the matrix, its symmetric
//     for (size_t i = 0; i < fAxes.GetNBins(); i++)
//     {
//         fCompatibleBins.at(i).push_back(i); // always true
//         for (size_t j = i + 1; j < fAxes.GetNBins(); j++)
//         {
//             if (BinsCompatible(i, j))
//             {
//                 fCompatibleBins.at(i).push_back(j);
//                 fCompatibleBins.at(j).push_back(i);
//             }
//         }
//     }

//     std::vector<size_t> relativeIndices = fTransObs.GetRelativeIndices(fDistObs);
//     const AxisCollection &axes = fAxes;

//     //  get the axes that this systematic will act on
//     fSubMapAxes = AxisCollection();
//     for (size_t i = 0; i < relativeIndices.size(); i++)
//         fSubMapAxes.AddAxis(axes.GetAxis(relativeIndices.at(i)));

//     // cache the equivilent index in the binning system of the systematic
//     fSysBins.resize(fAxes.GetNBins());
//     std::vector<size_t> sysIndices(relativeIndices.size(), 0);
//     for (size_t i = 0; i < axes.GetNBins(); i++)
//     {
//         for (size_t dim = 0; dim < relativeIndices.size(); dim++)
//             sysIndices[dim] = axes.UnflattenIndex(i, relativeIndices.at(dim));

//         fSysBins[i] = fSubMapAxes.FlattenIndices(sysIndices);
//     }
//     fCachedCompatibleBins = true;
// }

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
