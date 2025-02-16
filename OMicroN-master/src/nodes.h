#ifndef NODES_H
#define NODES_H

#include <eigen3/Eigen/Dense>

struct Node {
    int id;                        // Unique node ID
    Eigen::Vector3f position;      // Node position (x, y, z)
    Eigen::Vector3f displacement;  // Displacement (dx, dy, dz)
    bool isFixed[3];               // Boundary condition flags (x, y, z)

    Node(int nodeId, const Eigen::Vector3f& pos)
        : id(nodeId), position(pos), displacement(Eigen::Vector3f::Zero()) {
        isFixed[0] = isFixed[1] = isFixed[2] = false;
    }

    void setFixed(int direction) { isFixed[direction] = true; }
};

#endif // NODE_H
