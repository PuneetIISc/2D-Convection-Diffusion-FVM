#include "Mesh.hpp"

#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <cmath>

void Mesh::read(const std::string &geometryFile,
                const std::string &nodeSupportFile,
                const std::string &cellFile)
{
    readGeometry(geometryFile);
    readNodeSupport(nodeSupportFile);
    readCells(cellFile);
    computeCovolumes();
}

void Mesh::readGeometry(const std::string &file)
{
    std::ifstream in(file.c_str());
    if (!in) { std::cerr << "Error : Can't open " << file << std::endl; std::exit(1); }

    int np, nc, nf;
    in >> np >> nc >> nf;

    nodes_.resize(np);
    for (int i = 0; i < np; ++i)
        in >> nodes_[i].coord.x >> nodes_[i].coord.y >> nodes_[i].bc;

    cells_.resize(nc);
    for (int i = 0; i < nc; ++i)
        in >> cells_[i].centroid.x >> cells_[i].centroid.y >> cells_[i].volume;

    faces_.resize(nf);
    for (int i = 0; i < nf; ++i)
    {
        int n1, n2, in_, out_, bc;
        double sx, sy;
        in >> n1 >> n2 >> in_ >> out_ >> sx >> sy >> bc;
        Face &f = faces_[i];
        f.node1   = n1 - 1;
        f.node2   = n2 - 1;
        f.incell  = (in_  == 0) ? -1 : in_  - 1;
        f.outcell = (out_ == 0) ? -1 : out_ - 1;
        f.sx = sx; f.sy = sy; f.bc = bc;
    }
}

void Mesh::readNodeSupport(const std::string &file)
{
    std::ifstream in(file.c_str());
    if (!in) { std::cerr << "Error : Can't open " << file << std::endl; std::exit(1); }

    const int np = numNodes();
    support_.resize(np);

    std::string line;
    for (int i = 0; i < np; ++i)               
    {
        std::getline(in, line);
        if (line.empty()) { --i; continue; }
        std::istringstream ss(line);
        int count; ss >> count;
        for (int k = 0; k < count; ++k) { int c; ss >> c; support_[i].cells.push_back(c - 1); }
    }
    for (int i = 0; i < np; ++i)                 
    {
        std::getline(in, line);
        if (line.empty()) { --i; continue; }
        std::istringstream ss(line);
        int count; ss >> count;
        for (int k = 0; k < count; ++k) { int m; ss >> m; support_[i].nodes.push_back(m - 1); }
    }
}

void Mesh::readCells(const std::string &file)
{
    std::ifstream in(file.c_str());
    if (!in) { std::cerr << "Error : Can't open " << file << std::endl; std::exit(1); }

    const int nc = numCells();
    for (int i = 0; i < nc; ++i)                
    {
        int nn; in >> nn;
        cells_[i].node.resize(nn);
        for (int k = 0; k < nn; ++k) { int v; in >> v; cells_[i].node[k] = v - 1; }
    }
    for (int i = 0; i < nc; ++i)                
    {
        int nf; in >> nf;
        cells_[i].face.resize(nf);
        for (int k = 0; k < nf; ++k) { int v; in >> v; cells_[i].face[k] = v - 1; }
    }
}

// co_volume(i) = area of triangle (n1,n2,incell) + area of triangle (n1,n2,outcell),
void Mesh::computeCovolumes()
{
    for (std::size_t i = 0; i < faces_.size(); ++i)
    {
        Face &f = faces_[i];
        const Vec2 &A = nodes_[f.node1].coord;
        const Vec2 &B = nodes_[f.node2].coord;
        double cov = 0.0;

        if (f.incell >= 0)
        {
            const Vec2 &C = cells_[f.incell].centroid;
            const double vol = std::fabs((A.x * B.y + B.x * C.y + C.x * A.y)
                                       - (A.y * B.x + B.y * C.x + C.y * A.x));
            cov += 0.5 * vol;
        }
        if (f.outcell >= 0)
        {
            const Vec2 &C = cells_[f.outcell].centroid;
            const double vol = std::fabs((A.x * B.y + B.x * C.y + C.x * A.y)
                                       - (A.y * B.x + B.y * C.x + C.y * A.x));
            cov += 0.5 * vol;
        }
        f.covolume = cov;
    }
}
