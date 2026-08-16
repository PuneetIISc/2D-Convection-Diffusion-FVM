#include "Input.hpp"

#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>

namespace
{
    std::istringstream nextLine(std::ifstream &in)
    {
        std::string line;
        while (std::getline(in, line))
        {
            const std::size_t p = line.find_first_not_of(" \t\r\n");
            if (p == std::string::npos) continue;
            if (line[p] == '!' || line[p] == '#') continue;
            return std::istringstream(line);
        }
        return std::istringstream("");
    }
}

void Input::read(const std::string &file)
{
    std::ifstream in(file.c_str());
    if (!in)
    {
        std::cerr << "Error : cannot open '" << file << "'" << std::endl;
        std::exit(1);
    }
    nextLine(in) >> Pe >> L >> uvel >> vvel;
    nextLine(in) >> order >> nodal_value_option >> timestep;
    nextLine(in) >> mode_impli >> rkstage >> cfl >> cfl_option;
    nextLine(in) >> no_of_sweeps >> cflmulti >> maxcfl >> relax;
    nextLine(in) >> itermax >> maximum_residue;
}

void Input::echo() const
{
    std::cout << "Pe = " << Pe << "   alpha = " << uvel * L / Pe << "\n";
    std::cout << "order = " << order
              << "   nodal_value_option = " << nodal_value_option
              << "   timestep = " << timestep << "\n";
    std::cout << "mode = " << (mode_impli == 0 ? "explicit RK" : "implicit SGS");
    if (mode_impli == 0) std::cout << " (" << rkstage << " stages)";
    else                 std::cout << " (" << no_of_sweeps << " sweeps)";
    std::cout << "\n\n";
}

void MeshFiles::read(const std::string &file)
{
    std::ifstream in(file.c_str());
    if (!in)
    {
        std::cerr << "Mesh Files details file : file.input - not found " << std::endl;
        std::exit(1);
    }
    int n = 0;
    in >> n;
    in >> geometry >> nodesupport >> cell;
}
