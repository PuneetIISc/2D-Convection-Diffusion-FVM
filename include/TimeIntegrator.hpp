#ifndef TIMEINTEGRATOR_HPP
#define TIMEINTEGRATOR_HPP

#include <vector>

class Solver;

class TimeIntegrator
{
public:
    virtual ~TimeIntegrator() {}
    virtual void advance(Solver &s) = 0;   // computes cphi_new from cphi
};

// Explicit multistage Runge-Kutta: stage denominators are count = rkstage..1
// (dppar = 1,2,3,4 in the F90), so the last stage is a full step.
class ExplicitRK : public TimeIntegrator
{
public:
    explicit ExplicitRK(int stages) : stages_(stages) {}
    void advance(Solver &s);
private:
    int stages_;
    std::vector<double> flux_, dphiexp_;
};

// Implicit symmetric Gauss-Seidel.
class ImplicitSGS : public TimeIntegrator
{
public:
    ImplicitSGS() {}
    void advance(Solver &s);
private:
    std::vector<double> flux_, dphiexp_, dphimp_;
};

#endif
