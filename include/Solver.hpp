#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <vector>
#include "Mesh.hpp"
#include "Input.hpp"
#include "Convection.hpp"

struct ResidualRecord { int iter; double value; };

class Solver
{
public:
    Solver(const Mesh &mesh, const Input &input);
    ~Solver();

    void run();

    // --- pieces used by the time integrators -------------------------------
    void assembleFlux(std::vector<double> &flux);        // the two face loops
    void sgs(const std::vector<double> &duexp,
             std::vector<double> &duimp) const;          // implicit smoother

    std::vector<double>       &cphi()     { return cphi_; }
    std::vector<double>       &cphiNew()  { return cphi_new_; }
    const std::vector<double> &dt()       const { return dt_; }
    const Mesh                &mesh()      const { return mesh_; }
    const Input               &input()     const { return input_; }
    int                        numCells() const { return mesh_.numCells(); }

    // --- results ------------------------------------------------------------
    double exactCell(int i) const { return phi_exact_[i]; }
    const std::vector<double> &nodeValues() const { return vphi_; }
    const std::vector<double> &nodeExact()  const { return vphi_exact_; }
    const std::vector<ResidualRecord> &history() const { return history_; }

private:
    void   initialise();
    void   exactSolution();
    void   timeStep(int iteration);
    void   vertexInterpolation();          // fills vphi_ (area or pseudo-Laplacian)
    void   areaWeighted();
    void   pseudoLaplacian();
    void   outletLeastSquares();
    void   gradientGG(const double *px, const double *py, const double *p,
                      int cc, double &dpx, double &dpy) const;
    double convectiveFlux(int face, double nx, double ny) const;
    double residue(int iteration);
    double errorCell() const;

    const Mesh  &mesh_;
    const Input &input_;
    double pi_;

    std::vector<double> cphi_, cphi_new_, vphi_;
    std::vector<double> phi_exact_, vphi_exact_;
    std::vector<double> dt_;
    std::vector<Vec2>   dphi_;      // cell gradients

    double cfl_;                    // current CFL (ramped for implicit)
    double res1_;                   // first-iteration residual

    ConvectionScheme *scheme_;
    class TimeIntegrator *integrator_;
    std::vector<ResidualRecord> history_;
};

#endif
