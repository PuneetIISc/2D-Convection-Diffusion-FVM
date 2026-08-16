/************************************************************************
 *            2-D CONVECTION DIFFUSION SOLVER                           *                
 *            using FVM Framework for unstructured mesh                 *
 *        with Explicit and Implicit Time Integration                   *                                  
 *                              -                                       *
 *  This Code is written as Project Assignment for FVM Course at IISc   *
 *              under Guidance and Support of Bala Sir                  *
 *                    by Puneet Pushkar, M.Tech (AE)                    *
 ************************************************************************/


#include <iostream>

#include "Input.hpp"
#include "Mesh.hpp"
#include "Solver.hpp"
#include "Writer.hpp"

int main(int argc, char **argv)
{
    const std::string choiceFile = (argc > 1) ? argv[1] : "userchoice.inp";
    const std::string listFile   = (argc > 2) ? argv[2] : "file.input";

    Input input;
    std::cout << "Reading Choice Input File" << "\n";
    input.read(choiceFile);
    input.echo();

    MeshFiles files;
    std::cout << "Reading Mesh File" << "\n";
    files.read(listFile);

    Mesh mesh;
    mesh.read(files.geometry, files.nodesupport, files.cell);
    std::cout << "Number of points = " << mesh.numNodes() << "\n"
              << "Number of cells  = " << mesh.numCells() << "\n"
              << "Number of faces  = " << mesh.numFaces() << "\n\n";

    Solver solver(mesh, input);
    std::cout << "Update starts ..." << "\n";
    solver.run();

    Writer::solutionDat("Solution.dat", solver);
    Writer::vtk        ("Solution.vtk", solver);
    Writer::residue    ("residue.dat",  solver);
    std::cout << "  wrote Solution.dat, Solution.vtk, residue.dat\n";

    return 0;
}
