#include "Solver.hpp"
#include "TimeIntegrator.hpp"

#include <cmath>
#include <iostream>
#include <iomanip>

Solver::Solver(const Mesh &mesh, const Input &input)
    : mesh_(mesh), input_(input),
      cphi_(mesh.numCells(), 0.0),
      cphi_new_(mesh.numCells(), 0.0),
      vphi_(mesh.numNodes(), 0.0),
      phi_exact_(mesh.numCells(), 0.0),
      vphi_exact_(mesh.numNodes(), 0.0),
      dt_(mesh.numCells(), 0.0),
      dphi_(mesh.numCells(), Vec2()),
      cfl_(input.cfl), res1_(0.0),
      scheme_(0), integrator_(0)
{
    pi_ = std::atan(1.0) * 4.0;
    scheme_     = (input.order == 2) ? static_cast<ConvectionScheme *>(new SecondOrder())
                                     : static_cast<ConvectionScheme *>(new FirstOrder());
    integrator_ = (input.mode_impli == 1) ? static_cast<TimeIntegrator *>(new ImplicitSGS())
                                          : static_cast<TimeIntegrator *>(new ExplicitRK(input.rkstage));
}

Solver::~Solver()
{
    delete scheme_;
    delete integrator_;
}

void Solver::exactSolution()
{
    const double Pe = input_.Pe, L = input_.L;
    const double m1 = Pe / 2.0;
    const double m2 = m1 * m1 + pi_ * pi_;
    const double r1 = m1 + std::sqrt(m2);
    const double r2 = m1 - std::sqrt(m2);
    const double deno = r2 * std::exp(r2 * L) - r1 * std::exp(r1 * L);

    for (int i = 0; i < mesh_.numCells(); ++i)
    {
        const Vec2 &c = mesh_.cells()[i].centroid;
        const double num = r2 * std::exp(r1 * c.x + r2 * L) - r1 * std::exp(r1 * L + r2 * c.x);
        phi_exact_[i] = std::sin(pi_ * c.y) * (num / deno);
    }
    for (int i = 0; i < mesh_.numNodes(); ++i)
    {
        const Vec2 &p = mesh_.nodes()[i].coord;
        const double num = r2 * std::exp(r1 * p.x + r2 * L) - r1 * std::exp(r1 * L + r2 * p.x);
        vphi_exact_[i] = std::sin(pi_ * p.y) * (num / deno);
    }
}

// initialise_variables : zero the field, set inlet node values.
void Solver::initialise()
{
    std::fill(cphi_.begin(), cphi_.end(), 0.0);
    std::fill(cphi_new_.begin(), cphi_new_.end(), 0.0);
    std::fill(vphi_.begin(), vphi_.end(), 0.0);
    for (int i = 0; i < mesh_.numNodes(); ++i)
        if (mesh_.nodes()[i].bc == BC_INLET)
            vphi_[i] = std::sin(pi_ * mesh_.nodes()[i].coord.y);
}

// time_step_calculation : dt = cfl * sqrt(volume) / |velocity|.
void Solver::timeStep(int iteration)
{
    if (input_.mode_impli != 0)                 
    {
        if (input_.cfl_option == 1)
            cfl_ = input_.cflmulti;
        else if (input_.cfl_option == 2)
            cfl_ = input_.cflmulti * static_cast<double>(iteration);
        else if (input_.cfl_option == 3)
        {
            cfl_ = input_.cflmulti * static_cast<double>(iteration);
            if (cfl_ >= input_.maxcfl) cfl_ = input_.maxcfl;
        }
    }

    const double q = std::sqrt(input_.uvel * input_.uvel + input_.vvel * input_.vvel);
    for (int i = 0; i < mesh_.numCells(); ++i)
    {
        const double sl = std::sqrt(mesh_.cells()[i].volume);
        dt_[i] = (cfl_ * sl) / q;
    }

    if (input_.timestep == 1)                   // global: use the smallest dt
    {
        double mintime = dt_[0];
        for (int i = 1; i < mesh_.numCells(); ++i)
            if (dt_[i] < mintime) mintime = dt_[i];
        std::fill(dt_.begin(), dt_.end(), mintime);
    }
}

// gradient_gg
void Solver::gradientGG(const double *px, const double *py, const double *p,
                        int cc, double &dpx, double &dpy) const
{
    dpx = 0.0; dpy = 0.0;
    for (int k = 0; k < cc; ++k)
    {
        const int k2 = (k == cc - 1) ? 0 : k + 1;
        const double norm_x =  (py[k2] - py[k]);
        const double norm_y = -(px[k2] - px[k]);
        dpx += 0.5 * (p[k] + p[k2]) * norm_x;
        dpy += 0.5 * (p[k] + p[k2]) * norm_y;
    }
}

// vertex_value_interpolation
void Solver::vertexInterpolation()
{
    if (input_.nodal_value_option == 2) pseudoLaplacian();
    else                                areaWeighted();
}

// area_weighted_average : vphi = sum(volume*cphi)/sum(volume), then boundaries.
void Solver::areaWeighted()
{
    const int np = mesh_.numNodes();
    for (int i = 0; i < np; ++i)
    {
        double sum = 0.0, sw = 0.0;
        const std::vector<int> &cellsAt = mesh_.support()[i].cells;
        for (std::size_t j = 0; j < cellsAt.size(); ++j)
        {
            const double w = mesh_.cells()[cellsAt[j]].volume;
            sum += cphi_new_[cellsAt[j]] * w;
            sw  += w;
        }
        vphi_[i] = (sw > 0.0) ? sum / sw : 0.0;
    }
    for (int i = 0; i < np; ++i)
    {
        if (mesh_.nodes()[i].bc == BC_INLET)
            vphi_[i] = std::sin(pi_ * mesh_.nodes()[i].coord.y);
        else if (mesh_.nodes()[i].bc == BC_WALL)
            vphi_[i] = 0.0;
    }
    outletLeastSquares();
}

// pseudo_laplacian
void Solver::pseudoLaplacian()
{
    const int np = mesh_.numNodes();
    for (int i = 0; i < np; ++i)
    {
        const Node &nd = mesh_.nodes()[i];
        if (nd.bc == BC_INTERIOR)
        {
            double Ixx = 0, Ixy = 0, Iyy = 0, Rx = 0, Ry = 0;
            const std::vector<int> &cellsAt = mesh_.support()[i].cells;
            for (std::size_t j = 0; j < cellsAt.size(); ++j)
            {
                const Vec2 &c = mesh_.cells()[cellsAt[j]].centroid;
                const double dxj = c.x - nd.coord.x;
                const double dyj = c.y - nd.coord.y;
                Ixx += dxj * dxj; Ixy += dxj * dyj; Iyy += dyj * dyj;
                Rx += dxj; Ry += dyj;
            }
            const double det = Ixx * Iyy - Ixy * Ixy;

            double sum = 0.0, sw = 0.0;
            for (std::size_t j = 0; j < cellsAt.size(); ++j)
            {
                const Vec2 &c = mesh_.cells()[cellsAt[j]].centroid;
                const double delx = c.x - nd.coord.x;
                const double dely = c.y - nd.coord.y;
                double w;
                if (det == 0.0) { w = 1.0; }
                else
                {
                    const double lx = (Ixy * Ry - Iyy * Rx) / det;
                    const double ly = (Ixy * Rx - Ixx * Ry) / det;
                    w = std::fabs(1.0 + lx * delx + ly * dely);
                    if (w == 0.0)      w = 1.0;
                    else if (w > 2.0)  w = 2.0;
                }
                sum += w * cphi_new_[cellsAt[j]];
                sw  += w;
            }
            vphi_[i] = sum / sw;
        }
        else if (nd.bc == BC_INLET)
            vphi_[i] = std::sin(pi_ * nd.coord.y);
        else if (nd.bc == BC_WALL)
            vphi_[i] = 0.0;
    }
    outletLeastSquares();
}

void Solver::outletLeastSquares()
{
    const int np = mesh_.numNodes();
    for (int i = 0; i < np; ++i)
    {
        if (mesh_.nodes()[i].bc != BC_OUTLET) continue;
        const Vec2 &pi = mesh_.nodes()[i].coord;

        double sumdxdy = 0, sumdysq = 0, sumdx = 0, sumdy = 0;
        double sumdphidx = 0, sumdphidy = 0;
        const std::vector<int> &nb = mesh_.support()[i].nodes;
        for (std::size_t j = 0; j < nb.size(); ++j)
        {
            const int k = nb[j];
            if (k < 0) continue;
            const double dx = mesh_.nodes()[k].coord.x - pi.x;
            const double dy = mesh_.nodes()[k].coord.y - pi.y;
            sumdxdy += dx * dy; sumdysq += dy * dy; sumdx += dx; sumdy += dy;
            sumdphidy += vphi_[k] * dy;
            sumdphidx += vphi_[k] * dx;
        }
        const double deno = sumdxdy * sumdy - sumdysq * sumdx;
        vphi_[i] = (sumdxdy * sumdphidy - sumdysq * sumdphidx) / deno;
    }
}

// update_variable :
// flux(i) = residual R(U) for each cell. Cell gradients are formed on the way.
void Solver::assembleFlux(std::vector<double> &flux)
{
    const int nc = mesh_.numCells();
    const std::vector<Face> &faces = mesh_.faces();
    const std::vector<Cell> &cells = mesh_.cells();
    const std::vector<Node> &nodes = mesh_.nodes();
    const double invPe = 1.0 / input_.Pe;

    vertexInterpolation();

    for (int i = 0; i < nc; ++i) dphi_[i] = Vec2();
    std::vector<double> vol_sum(nc, 0.0);
    std::fill(flux.begin(), flux.end(), 0.0);

    // ---- viscous (diffusion) loop -----------------------------------------
    for (std::size_t i = 0; i < faces.size(); ++i)
    {
        const Face &f = faces[i];
        const int in = f.incell, out = f.outcell, n1 = f.node1, n2 = f.node2;
        double dphix = 0.0, dphiy = 0.0;

        if (f.bc == BC_INTERIOR)
        {
            const double px[4] = {cells[in].centroid.x, nodes[n1].coord.x,
                                  cells[out].centroid.x, nodes[n2].coord.x};
            const double py[4] = {cells[in].centroid.y, nodes[n1].coord.y,
                                  cells[out].centroid.y, nodes[n2].coord.y};
            const double p[4]  = {cphi_new_[in], vphi_[n1], cphi_new_[out], vphi_[n2]};
            gradientGG(px, py, p, 4, dphix, dphiy);

            dphi_[in].x += dphix; dphi_[in].y += dphiy; vol_sum[in] += f.covolume;
            dphi_[out].x += dphix; dphi_[out].y += dphiy; vol_sum[out] += f.covolume;

            dphix /= f.covolume; dphiy /= f.covolume;
            const double rf = invPe * (dphix * f.sx + dphiy * f.sy);
            flux[out] += rf;
            flux[in]  -= rf;
        }
        else
        {
            if (in >= 0)
            {
                const double px[3] = {nodes[n1].coord.x, nodes[n2].coord.x, cells[in].centroid.x};
                const double py[3] = {nodes[n1].coord.y, nodes[n2].coord.y, cells[in].centroid.y};
                const double p[3]  = {vphi_[n1], vphi_[n2], cphi_new_[in]};
                gradientGG(px, py, p, 3, dphix, dphiy);

                if (f.bc == BC_OUTLET)
                {
                    dphix = 0.0;
                    dphiy = ((vphi_[n2] - vphi_[n1]) / (nodes[n2].coord.y - nodes[n1].coord.y))
                            * f.covolume;
                }
                dphi_[in].x += dphix; dphi_[in].y += dphiy; vol_sum[in] += f.covolume;

                dphix /= f.covolume; dphiy /= f.covolume;
                const double rf = invPe * (dphix * f.sx + dphiy * f.sy);
                flux[in] -= rf;
            }
            else if (out >= 0)
            {
                const double px[3] = {nodes[n1].coord.x, nodes[n2].coord.x, cells[out].centroid.x};
                const double py[3] = {nodes[n1].coord.y, nodes[n2].coord.y, cells[out].centroid.y};
                const double p[3]  = {vphi_[n1], vphi_[n2], cphi_new_[out]};
                gradientGG(px, py, p, 3, dphix, dphiy);

                if (f.bc == BC_OUTLET)
                {
                    dphix = 0.0;
                    dphiy = ((vphi_[n2] - vphi_[n1]) / (nodes[n2].coord.y - nodes[n1].coord.y))
                            * f.covolume;
                }
                dphi_[out].x += dphix; dphi_[out].y += dphiy; vol_sum[out] += f.covolume;

                dphix /= f.covolume; dphiy /= f.covolume;
                const double rf = invPe * (dphix * f.sx + dphiy * f.sy);
                flux[out] += rf;
            }
        }
    }

    for (int i = 0; i < nc; ++i)
    {
        dphi_[i].x /= vol_sum[i];
        dphi_[i].y /= vol_sum[i];
    }

    // ---- convective loop --------------------------------------------------
    for (std::size_t i = 0; i < faces.size(); ++i)
    {
        const Face &f = faces[i];
        const double area = std::sqrt(f.sx * f.sx + f.sy * f.sy);
        const double nx = f.sx / area;
        const double ny = f.sy / area;

        const double rf = convectiveFlux(static_cast<int>(i), nx, ny);

        if (f.incell < 0)       flux[f.outcell] -= rf;
        else if (f.outcell < 0) flux[f.incell]  += rf;
        else { flux[f.outcell] -= rf; flux[f.incell] += rf; }
    }
}

// convective_flux : upwind flux with first/second order reconstruction and the
// inlet/wall/outlet boundary states.
double Solver::convectiveFlux(int face, double nx, double ny) const
{
    const Face &f = mesh_.faces()[face];
    const std::vector<Cell> &cells = mesh_.cells();
    const std::vector<Node> &nodes = mesh_.nodes();

    const double xcf = 0.5 * (nodes[f.node1].coord.x + nodes[f.node2].coord.x);
    const double ycf = 0.5 * (nodes[f.node1].coord.y + nodes[f.node2].coord.y);
    const double area = std::sqrt(f.sx * f.sx + f.sy * f.sy);

    double phil = 0.0, phir = 0.0;
    if (f.incell >= 0)
    {
        const Vec2 d(xcf - cells[f.incell].centroid.x, ycf - cells[f.incell].centroid.y);
        phil = scheme_->state(cphi_new_[f.incell], dphi_[f.incell], d);
    }
    if (f.outcell >= 0)
    {
        const Vec2 d(xcf - cells[f.outcell].centroid.x, ycf - cells[f.outcell].centroid.y);
        phir = scheme_->state(cphi_new_[f.outcell], dphi_[f.outcell], d);
    }

    if (f.bc == BC_INLET)
    {
        if (f.incell < 0)       phil = std::sin(pi_ * ycf);
        else if (f.outcell < 0) phir = std::sin(pi_ * ycf);
    }
    else if (f.bc == BC_OUTLET || f.bc == BC_WALL)
    {
        if (f.incell < 0)       phil = 0.0;
        else if (f.outcell < 0) phir = 0.0;
    }

    const double un = input_.uvel * nx + input_.vvel * ny;
    double flux = 0.5 * (un + std::fabs(un)) * phil
                + 0.5 * (un - std::fabs(un)) * phir;
    return flux * area;
}

// S_G_S : symmetric Gauss-Seide
void Solver::sgs(const std::vector<double> &duexp, std::vector<double> &duimp) const
{
    const int nc = mesh_.numCells();
    const std::vector<Face> &faces = mesh_.faces();
    const std::vector<Cell> &cells = mesh_.cells();
    const std::vector<Node> &nodes = mesh_.nodes();
    const double alpha = 1.0 / input_.Pe;

    std::fill(duimp.begin(), duimp.end(), 0.0);

    for (int sweep = 1; sweep <= input_.no_of_sweeps; ++sweep)
    {
        int istart, iend, incr;
        if (sweep % 2 != 0) { istart = 0;      iend = nc;  incr = 1; }
        else                { istart = nc - 1; iend = -1;  incr = -1; }

        for (int i = istart; i != iend; i += incr)
        {
            double lhs = 0.0, rhs = 0.0;
            const std::vector<int> &cf = cells[i].face;
            for (std::size_t j = 0; j < cf.size(); ++j)
            {
                const int jf = cf[j];
                const Face &f = faces[jf];
                const double area = std::sqrt(f.sx * f.sx + f.sy * f.sy);
                double nx = f.sx / area, ny = f.sy / area;
                const double xm = 0.5 * (nodes[f.node1].coord.x + nodes[f.node2].coord.x);
                const double ym = 0.5 * (nodes[f.node1].coord.y + nodes[f.node2].coord.y);

                if (f.bc == BC_INTERIOR)
                {
                    int in, out;
                    if (i == f.incell) { in = i; out = f.outcell; }
                    else { nx = -nx; ny = -ny; in = i; out = f.incell; }

                    const double nor_vel = input_.uvel * nx + input_.vvel * ny;
                    const double dist_x = cells[out].centroid.x - cells[in].centroid.x;
                    const double dist_y = cells[out].centroid.y - cells[in].centroid.y;
                    const double proj = dist_x * nx + dist_y * ny;
                    const double diff_add = alpha / proj;

                    const double lhs_add = std::fabs(nor_vel) / 2.0 + diff_add;
                    lhs += lhs_add * area;
                    const double rhs_add = (nor_vel - std::fabs(nor_vel)) / 2.0 - diff_add;
                    rhs += rhs_add * duimp[out] * area;
                }
                else
                {
                    int in;
                    if (i == f.incell) { in = i; }
                    else { nx = -nx; ny = -ny; in = i; }

                    const double nor_vel = input_.uvel * nx + input_.vvel * ny;
                    const double dist_x = xm - cells[in].centroid.x;
                    const double dist_y = ym - cells[in].centroid.y;
                    const double proj = dist_x * nx + dist_y * ny;
                    const double diff_add = alpha / proj;

                    const double lhs_add = std::fabs(nor_vel) / 2.0 + diff_add;
                    lhs += lhs_add * area;
                }
            }
            const double factor = dt_[i] / cells[i].volume;
            lhs = 1.0 + factor * lhs;
            rhs = factor * rhs;
            duimp[i] = (duexp[i] - rhs) / lhs;
        }
    }

    for (int i = 0; i < nc; ++i) duimp[i] *= input_.relax;
}

// residue_calculation :
double Solver::residue(int iteration)
{
    double resid = 0.0;
    for (int i = 0; i < mesh_.numCells(); ++i)
    {
        const double diff = std::fabs((cphi_new_[i] - cphi_[i]) / dt_[i]);
        resid += diff * diff;
    }
    resid = std::sqrt(resid / mesh_.numCells());
    if (iteration == 1) res1_ = resid;
    return resid / res1_;
}

double Solver::errorCell() const
{
    double error = 0.0, sumarea = 0.0;
    for (int i = 0; i < mesh_.numCells(); ++i)
    {
        const double err = phi_exact_[i] - cphi_[i];
        error   += err * err * mesh_.cells()[i].volume;
        sumarea += mesh_.cells()[i].volume;
    }
    return std::sqrt(error / sumarea);
}

void Solver::run()
{
    exactSolution();
    initialise();

    for (int iter = 1; iter <= input_.itermax; ++iter)
    {
        timeStep(iter);
        integrator_->advance(*this);

        const double rel = residue(iter);
        history_.push_back((ResidualRecord){iter, rel});

            std::cout << "  iter " << std::setw(6) << iter
                      << "   residue " << std::scientific << std::setprecision(4)
                      << rel << "   cfl " << cfl_ << "\n";

        for (int i = 0; i < mesh_.numCells(); ++i) cphi_[i] = cphi_new_[i];

        if (rel <= input_.maximum_residue || iter == input_.itermax)
        {
            std::cout << "  iter " << std::setw(6) << iter
                      << "   residue " << std::scientific << rel << "   stop\n";
            areaWeighted();     // final node values for output
            std::cout << "\n  global L2 error = " << std::scientific
                      << std::setprecision(6) << errorCell() << "\n";
            break;
        }
    }
}
