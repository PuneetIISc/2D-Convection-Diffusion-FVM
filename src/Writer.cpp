#include "Writer.hpp"
#include <fstream>
#include <iomanip>

namespace Writer
{

void solutionDat(const std::string &file, const Solver &s)
{
    const Mesh &mesh = s.mesh();
    const std::vector<double> &vphi = s.nodeValues();

    std::ofstream out(file.c_str());
    out << std::setprecision(9);
    out << "TITLE = \"FVM\"\n";
    out << "VARIABLES = \"X\", \"Y\", \"Phi\"\n";
    out << "ZONE N = " << mesh.numNodes() << "  E = " << mesh.numCells()
        << " , F=FEPOINT  ET=QUADRILATERAL\n";

    for (int i = 0; i < mesh.numNodes(); ++i)
    {
        const Vec2 &p = mesh.nodes()[i].coord;
        out << p.x << " " << p.y << " " << vphi[i] << "\n";
    }
    for (int i = 0; i < mesh.numCells(); ++i)
    {
        const std::vector<int> &nd = mesh.cells()[i].node;
        out << nd[0] + 1 << " " << nd[1] + 1 << " " << nd[2] + 1 << " " << nd[2] + 1 << "\n";
    }
}

void vtk(const std::string &file, const Solver &s)
{
    const Mesh &mesh = s.mesh();
    const std::vector<double> &vphi = s.nodeValues();

    std::ofstream out(file.c_str());
    out << std::setprecision(9);
    out << "# vtk DataFile Version 3.0\nFVM\nASCII\nDATASET UNSTRUCTURED_GRID\n";

    out << "POINTS " << mesh.numNodes() << " double\n";
    for (int i = 0; i < mesh.numNodes(); ++i)
    {
        const Vec2 &p = mesh.nodes()[i].coord;
        out << p.x << " " << p.y << " 0\n";
    }
    out << "CELLS " << mesh.numCells() << " " << mesh.numCells() * 4 << "\n";
    for (int i = 0; i < mesh.numCells(); ++i)
    {
        const std::vector<int> &nd = mesh.cells()[i].node;
        out << "3 " << nd[0] << " " << nd[1] << " " << nd[2] << "\n";
    }
    out << "CELL_TYPES " << mesh.numCells() << "\n";
    for (int i = 0; i < mesh.numCells(); ++i) out << "5\n";

    out << "POINT_DATA " << mesh.numNodes() << "\n";
    out << "SCALARS phi double 1\nLOOKUP_TABLE default\n";
    for (int i = 0; i < mesh.numNodes(); ++i) out << vphi[i] << "\n";
    out << "SCALARS phi_exact double 1\nLOOKUP_TABLE default\n";
    const std::vector<double> &vex = s.nodeExact();
    for (int i = 0; i < mesh.numNodes(); ++i) out << vex[i] << "\n";
}

void residue(const std::string &file, const Solver &s)
{
    std::ofstream out(file.c_str());
    const std::vector<ResidualRecord> &h = s.history();
    for (std::size_t i = 0; i < h.size(); ++i)
        out << h[i].iter << " " << std::scientific << std::setprecision(6)
            << h[i].value << "\n";
}

}
