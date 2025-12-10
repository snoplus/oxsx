#include "ShapeInterpConstraint.h"
#include <Exceptions.h>
#include <string>
#include <bitset>


std::vector<double> ShapeInterpConstraint::ExtractConstraintParams(const ParameterDict& params) const
{
    /*
    We only need the parameter values relevant to this particular constraint! Extract them
    from the full parameter dictionary object.
    */
    const std::vector<std::string> axis_names = fHist.GetAxisNames();
    std::vector<double> con_params;

    for (const auto &name : axis_names)
    {
        if (params.count(name) == 0)
        {
            throw NotFoundError("ShapeInterpConstraint::ExtractConstraintParams(): constraint param " + name + " not found in parameter dictionary given.");
        }
        con_params.push_back(params.at(name));
    }
    return con_params;
}


double ShapeInterpConstraint::NDLinearInterp(const std::vector<double>& grid_vals, const std::vector<std::vector<double>>& grid_locs, const std::vector<double>& con_params) const
{
    /*
     * Perform n-dimensional linear interpolation, given a sample point in nD space (con_params),
     * the locations of the bin centres which surround that point in the shape histogram (grid_locs),
     * and the values at those "lattice points" (grid_locs).
     * Uses a generalised version of the algorithm implied by:
     * https://en.wikipedia.org/wiki/Trilinear_interpolation#Formulation
     * We take advantage of the strong symmetry in this iterpolation alogrithm.
     */
    // First: get coordinates of sample point in "lattice" coordinates
    std::vector<double> lattice_coords;
    const size_t N_DIMS = con_params.size();
    for (size_t dim=0; dim<N_DIMS; dim++)
    {
        const double sample_val = con_params.at(dim);
        const double x0 = grid_locs.at(0).at(dim);
        const double x1 = grid_locs.at(1 << dim).at(dim);

        const double xd = (sample_val-x0)/(x1-x0);
        lattice_coords.push_back(xd);
    }
    // Now calculate interpolated value, via binary representation of lattice point index!
    const size_t N_LATT_POINTS = 1 << N_DIMS; // N = 2^N_DIMS
    double interp_val = 0;
    for (size_t latt_idx=0; latt_idx<N_LATT_POINTS; latt_idx++)
    {
        const std::bitset<20> latt_idx_bits(latt_idx); // (Arbitrary) max dimension of 20, so we can use bitset!
        double val = grid_vals.at(latt_idx);
        for (size_t dim=0; dim<N_DIMS; dim++)
        {
            val *= (latt_idx_bits.test(dim)) ? lattice_coords.at(dim) : (1.0 - lattice_coords.at(dim));
        }
        interp_val += val;
    }

    return interp_val;
}


double ShapeInterpConstraint::Evaluate(const ParameterDict& params) const
{
    // Extract relevant params from dict which constraint cares about
    const std::vector<double> con_params = ExtractConstraintParams(params);
    // Get grid point values relevant for this sample point
    std::vector<std::vector<double>> grid_locs;
    const std::vector<double> grid_vals = fHist.GetLatticePointValues(con_params, grid_locs);
    // Run n-D linear interpolation
    return NDLinearInterp(grid_vals, grid_locs, con_params);
}