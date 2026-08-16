#ifndef MESH_HPP
#define MESH_HPP

#include <string>
#include <vector>
#include "Vec2.hpp"

// Boundary codes as written in geometry.out.
const int BC_INTERIOR = 0;
const int BC_INLET    = 2;
const int BC_WALL     = 21;
const int BC_OUTLET   = 3;

struct Node
{
    Vec2 coord;
    int  bc = BC_INTERIOR;
};

struct Cell
{
    Vec2             centroid;
    double           volume = 0.0;
    std::vector<int> node;   // node_num  (cell -> nodes)
    std::vector<int> face;   // face_num  (cell -> faces), used by the SGS
};

struct Face
{
    int    node1 = -1, node2 = -1;
    int    incell = -1, outcell = -1;   // -1 means "no cell" (boundary side)
    double sx = 0.0, sy = 0.0;          // area-normal, incell -> outcell
    int    bc = BC_INTERIOR;
    double covolume = 0.0;              // diamond area (co_volume)
};

// Cells and nodes surrounding each node (from nodesupport.out).
struct NodeSupport
{
    std::vector<int> cells;
    std::vector<int> nodes;
};

// Reads geometry.out / nodesupport.out / cell.out and computes the co-volumes.
class Mesh
{
public:
    void read(const std::string &geometryFile,
              const std::string &nodeSupportFile,
              const std::string &cellFile);

    int numNodes() const { return static_cast<int>(nodes_.size()); }
    int numCells() const { return static_cast<int>(cells_.size()); }
    int numFaces() const { return static_cast<int>(faces_.size()); }

    const std::vector<Node>        &nodes()   const { return nodes_; }
    const std::vector<Cell>        &cells()   const { return cells_; }
    const std::vector<Face>        &faces()   const { return faces_; }
    const std::vector<NodeSupport> &support() const { return support_; }

private:
    void readGeometry(const std::string &file);
    void readNodeSupport(const std::string &file);
    void readCells(const std::string &file);
    void computeCovolumes();

    std::vector<Node>        nodes_;
    std::vector<Cell>        cells_;
    std::vector<Face>        faces_;
    std::vector<NodeSupport> support_;
};

#endif
