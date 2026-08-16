#include "TimeIntegrator.hpp"
#include "Solver.hpp"

// Explicit:
void ExplicitRK::advance(Solver &s)
{
    const int nc = s.numCells();
    const std::vector<Cell> &cells = s.mesh().cells();
    std::vector<double> &cphi     = s.cphi();
    std::vector<double> &cphi_new = s.cphiNew();

    flux_.assign(nc, 0.0);
    dphiexp_.assign(nc, 0.0);

    for (int count = stages_; count >= 1; --count)
    {
        s.assembleFlux(flux_);
        for (int i = 0; i < nc; ++i)
        {
            double d = (1.0 / cells[i].volume) * flux_[i];
            d = -(s.dt()[i] / static_cast<double>(count)) * d;
            cphi_new[i] = cphi[i] + d;
        }
    }
}

// Implicit:
void ImplicitSGS::advance(Solver &s)
{
    const int nc = s.numCells();
    const std::vector<Cell> &cells = s.mesh().cells();
    std::vector<double> &cphi     = s.cphi();
    std::vector<double> &cphi_new = s.cphiNew();

    flux_.assign(nc, 0.0);
    dphiexp_.assign(nc, 0.0);
    dphimp_.assign(nc, 0.0);

    s.assembleFlux(flux_);
    for (int i = 0; i < nc; ++i)
        dphiexp_[i] = -(s.dt()[i]) * (1.0 / cells[i].volume) * flux_[i];

    s.sgs(dphiexp_, dphimp_);
    for (int i = 0; i < nc; ++i)
        cphi_new[i] = cphi[i] + dphimp_[i];
}
