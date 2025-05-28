// bh_tree.cpp
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <vector>
#include <cmath>

namespace py = pybind11;

struct BHNode
{
    double center[3];
    double half_size;
    std::vector<int> indices;
    BHNode *children[8] = {nullptr};
    double mass = 0.0;
    double mass_center[3] = {0, 0, 0};
    bool is_leaf = true;

    BHNode(double cx, double cy, double cz, double hs, const std::vector<int> &idx)
        : half_size(hs), indices(idx)
    {
        center[0] = cx;
        center[1] = cy;
        center[2] = cz;
    }

    ~BHNode()
    {
        for (auto &c : children)
            if (c)
                delete c;
    }
};

void subdivide(BHNode *node, const double *pos, const double *masses, size_t n, double min_size = 1e-2)
{
    if (node->indices.size() <= 1 || node->half_size < min_size)
    {
        double msum = 0;
        double mcx = 0, mcy = 0, mcz = 0;
        for (auto idx : node->indices)
        {
            double m = masses[idx];
            msum += m;
            mcx += pos[3 * idx] * m;
            mcy += pos[3 * idx + 1] * m;
            mcz += pos[3 * idx + 2] * m;
        }
        node->mass = msum;
        if (msum > 0)
        {
            node->mass_center[0] = mcx / msum;
            node->mass_center[1] = mcy / msum;
            node->mass_center[2] = mcz / msum;
        }
        else
        {
            node->mass_center[0] = node->center[0];
            node->mass_center[1] = node->center[1];
            node->mass_center[2] = node->center[2];
        }
        return;
    }

    node->is_leaf = false;
    std::vector<int> children_indices[8];
    double hs2 = node->half_size / 2;

    for (auto idx : node->indices)
    {
        const double *p = pos + 3 * idx;
        int oct = 0;
        if (p[0] > node->center[0])
            oct |= 4;
        if (p[1] > node->center[1])
            oct |= 2;
        if (p[2] > node->center[2])
            oct |= 1;
        children_indices[oct].push_back(idx);
    }

    for (int i = 0; i < 8; i++)
    {
        double offset[3] = {((i & 4) ? 0.5 : -0.5) * node->half_size * 2,
                            ((i & 2) ? 0.5 : -0.5) * node->half_size * 2,
                            ((i & 1) ? 0.5 : -0.5) * node->half_size * 2};
        node->children[i] = new BHNode(
            node->center[0] + offset[0],
            node->center[1] + offset[1],
            node->center[2] + offset[2],
            hs2,
            children_indices[i]);
        subdivide(node->children[i], pos, masses, n, min_size);
    }

    node->mass = 0;
    double mcx = 0, mcy = 0, mcz = 0;
    for (int i = 0; i < 8; i++)
    {
        BHNode *c = node->children[i];
        node->mass += c->mass;
        mcx += c->mass * c->mass_center[0];
        mcy += c->mass * c->mass_center[1];
        mcz += c->mass * c->mass_center[2];
    }
    if (node->mass > 0)
    {
        node->mass_center[0] = mcx / node->mass;
        node->mass_center[1] = mcy / node->mass;
        node->mass_center[2] = mcz / node->mass;
    }
    else
    {
        node->mass_center[0] = node->center[0];
        node->mass_center[1] = node->center[1];
        node->mass_center[2] = node->center[2];
    }
}

void acc_from_node(BHNode *node, const double *pos, int i, double theta, double G, double softening, double *acc)
{
    if (node->mass == 0)
        return;

    const double *p = pos + 3 * i;
    if (node->is_leaf && node->indices.size() == 1 && node->indices[0] == i)
        return;

    double rx = node->mass_center[0] - p[0];
    double ry = node->mass_center[1] - p[1];
    double rz = node->mass_center[2] - p[2];
    double dist = std::sqrt(rx * rx + ry * ry + rz * rz) + 1e-10;

    if (node->is_leaf || (node->half_size / dist < theta))
    {
        double denom = dist * dist * dist + softening * softening * softening;
        double f = G * node->mass / denom;
        acc[0] += f * rx;
        acc[1] += f * ry;
        acc[2] += f * rz;
    }
    else
    {
        for (int c = 0; c < 8; c++)
        {
            if (node->children[c])
            {
                acc_from_node(node->children[c], pos, i, theta, G, softening, acc);
            }
        }
    }
}

py::array_t<double> compute_acceleration(py::array_t<double> pos_, py::array_t<double> masses_, double theta = 0.5, double G = 4.302e-3, double softening = 10.0)
{
    auto pos = pos_.unchecked<2>();
    auto masses = masses_.unchecked<1>();
    size_t n = pos.shape(0);

    // Compute bounding box center and half size
    double minx = pos(0, 0), maxx = pos(0, 0);
    double miny = pos(0, 1), maxy = pos(0, 1);
    double minz = pos(0, 2), maxz = pos(0, 2);
    for (size_t i = 1; i < n; i++)
    {
        if (pos(i, 0) < minx)
            minx = pos(i, 0);
        if (pos(i, 0) > maxx)
            maxx = pos(i, 0);
        if (pos(i, 1) < miny)
            miny = pos(i, 1);
        if (pos(i, 1) > maxy)
            maxy = pos(i, 1);
        if (pos(i, 2) < minz)
            minz = pos(i, 2);
        if (pos(i, 2) > maxz)
            maxz = pos(i, 2);
    }
    double cx = (minx + maxx) / 2;
    double cy = (miny + maxy) / 2;
    double cz = (minz + maxz) / 2;
    double max_range = std::max({maxx - minx, maxy - miny, maxz - minz}) / 2;

    std::vector<int> all_indices(n);
    for (size_t i = 0; i < n; i++)
        all_indices[i] = i;

    BHNode *root = new BHNode(cx, cy, cz, max_range, all_indices);
    subdivide(root, pos.data(0), masses.data(0), n);

    // Output acceleration array
    py::array_t<double> acc({n, 3});
    auto acc_mut = acc.mutable_unchecked<2>();
    for (size_t i = 0; i < n; i++)
    {
        double a[3] = {0, 0, 0};
        acc_from_node(root, pos.data(0), i, theta, G, softening, a);
        acc_mut(i, 0) = a[0];
        acc_mut(i, 1) = a[1];
        acc_mut(i, 2) = a[2];
    }

    delete root;
    return acc;
}

PYBIND11_MODULE(bh_tree, m)
{
    m.def("compute_acceleration", &compute_acceleration,
          py::arg("pos"), py::arg("masses"), py::arg("theta") = 0.5, py::arg("G") = 4.302e-3, py::arg("softening") = 10.0,
          "Compute acceleration using Barnes-Hut tree");
}
