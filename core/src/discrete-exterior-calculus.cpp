// PLEASE READ:
//
// This file additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because we are
// "inside" the class, we no longer have to call
//
//          geometry->inputVertexPositions[v], etc.
//
// We can just call
//
//          this->inputVertexPositions[v], etc.
//
// or simply
//
//          inputVertexPositions[v], etc.
//
// In addition, we no longer access the corresponding surface mesh via
//
//          mesh->vertices(), etc.
//
// but instead <mesh> is not a pointer anymore, so we use
//
//          mesh.vertices(), etc.
//
// Functions in this file can be called from other projects simply by using geometry->buildHodgeStar0Form(), etc. where
// "geometry" is a pointer to a VertexPositionGeometry. This avoids having to declare a GeometryRoutines object in every
// project, and also mimics the way that geometry routines are normally called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.

#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral {
namespace surface {


/*
 * Build Hodge operator on 0-forms.
 * By convention, the area of a vertex is 1.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar0Form() const {
    Vector<double> result = Vector<double>::Zero(mesh.nVertices());
    for (Vertex v : mesh.vertices()) {
        result[v.getIndex()] = barycentricDualArea(v);
    }
    return SparseMatrix<double>{result.asDiagonal()};
}

/*
 * Build Hodge operator on 1-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar1Form() const {
    Vector<double> result = Vector<double>::Zero(mesh.nEdges());
    for (Edge e : mesh.edges()) {
        Halfedge he = e.halfedge();
        result[e.getIndex()] = 0.5 * (cotan(he) + cotan(he.twin()));
    }
    return SparseMatrix<double>{result.asDiagonal()};
}

/*
 * Build Hodge operator on 2-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 2-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar2Form() const {
    Vector<double> result = Vector<double>::Zero(mesh.nFaces());
    for (Face f : mesh.faces()) {
        Halfedge ha = f.halfedge();
        Halfedge hb = ha.next();

        Vector3 a = halfedgeVector(ha);
        Vector3 b = halfedgeVector(hb);

        result[f.getIndex()] = 2.0 / cross(a, b).norm();
    }
    return SparseMatrix<double>{result.asDiagonal()};
}

/*
 * Build exterior derivative on 0-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative0Form() const {
    auto rows = mesh.nEdges();
    auto columns = mesh.nVertices();

    std::vector<Eigen::Triplet<double>> triplets;
    // Each row has 2 non-zero entries, because each edge is adjacent to 2 vertices.
    triplets.reserve(2 * rows);

    for (Edge e : mesh.edges()) {
        auto i = e.getIndex();
        Vertex v1 = e.firstVertex();
        Vertex v2 = e.secondVertex();
        triplets.push_back(Eigen::Triplet<double>(i, v1.getIndex(), -1.0));
        triplets.push_back(Eigen::Triplet<double>(i, v2.getIndex(),  1.0));
    }

    SparseMatrix<double> matrix(rows, columns);
    matrix.setFromTriplets(triplets.begin(), triplets.end());
    return matrix;
}

/*
 * Build exterior derivative on 1-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative1Form() const {
    auto rows = mesh.nFaces();
    auto columns = mesh.nEdges();

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(3 * rows);

    for (Face f : mesh.faces()) {
        auto i = f.getIndex();
        for (Edge e : f.adjacentEdges()) {
            auto j = e.getIndex();
            bool matchesOrientation = (e.halfedge().face() == f);
            triplets.push_back(Eigen::Triplet<double>(i, j, matchesOrientation ? 1.0 : -1.0));
        }
    }

    SparseMatrix<double> matrix(rows, columns);
    matrix.setFromTriplets(triplets.begin(), triplets.end());    
    return matrix;
}

} // namespace surface
} // namespace geometrycentral