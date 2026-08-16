#ifndef INPUT_HPP
#define INPUT_HPP

#include <string>

// Run settings, read from userchoice.inp
struct Input
{
    double Pe = 100.0, L = 1.0, uvel = 1.0, vvel = 0.0;
    int order = 1, nodal_value_option = 1, timestep = 2;
    int mode_impli = 0, rkstage = 3;
    double cfl = 0.4;
    int cfl_option = 1;
    int no_of_sweeps = 4;
    double cflmulti = 1.0, maxcfl = 0.4, relax = 1.0;
    int itermax = 50000;
    double maximum_residue = 1.0e-8;

    void read(const std::string &file);
    void echo() const;
};

// file.input : number of files, then the three mesh file names.
struct MeshFiles
{
    std::string geometry, nodesupport, cell;
    void read(const std::string &file);
};

#endif
