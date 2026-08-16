#ifndef WRITER_HPP
#define WRITER_HPP

#include <string>
#include "Solver.hpp"

namespace Writer
{
    void solutionDat(const std::string &file, const Solver &s);  // Tecplot
    void vtk        (const std::string &file, const Solver &s);  // ParaView
    void residue    (const std::string &file, const Solver &s);
}

#endif
