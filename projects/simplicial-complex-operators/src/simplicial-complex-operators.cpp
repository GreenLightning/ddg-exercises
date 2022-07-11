// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:
    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        idx = geometry->vertexIndices[v];
    }

    for (Edge e : mesh->edges()) {
        idx = geometry->edgeIndices[e];
    }

    for (Face f : mesh->faces()) {
        idx = geometry->faceIndices[f];
    }

    // You can more easily get the indices of mesh elements using the function getIndex(), albeit less efficiently and
    // technically less safe (although you don't need to worry about it), like so:
    //
    //      v.getIndex()
    //
    // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

    for (Vertex v : mesh->vertices()) {
        idx = v.getIndex(); // == geometry->vertexIndices[v])
    }

    // Geometry Central already sets the indices for us, though, so this function is just here for demonstration.
    // You don't have to do anything :)
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {
    auto rows = mesh->nEdges();
    auto columns = mesh->nVertices();

    std::vector<Eigen::Triplet<size_t>> triplets;
    // Each row has 2 non-zero entries, because each edge is adjacent to 2 vertices.
    triplets.reserve(2 * rows);

    for (Edge e : mesh->edges()) {
        auto i = e.getIndex();
        Vertex v1 = e.firstVertex();
        Vertex v2 = e.secondVertex();
        triplets.push_back(Eigen::Triplet<size_t>(i, v1.getIndex(), 1));
        triplets.push_back(Eigen::Triplet<size_t>(i, v2.getIndex(), 1));
    }

    SparseMatrix<size_t> matrix(rows, columns);
    matrix.setFromTriplets(triplets.begin(), triplets.end());
    return matrix;
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {
    auto rows = mesh->nFaces();
    auto columns = mesh->nEdges();

    std::vector<Eigen::Triplet<size_t>> triplets;
    triplets.reserve(3 * rows);

    for (Face f : mesh->faces()) {
        auto i = f.getIndex();
        for (Edge e : f.adjacentEdges()) {
            auto j = e.getIndex();
            triplets.push_back(Eigen::Triplet<size_t>(i, j, 1));
        }
    }

    SparseMatrix<size_t> matrix(rows, columns);
    matrix.setFromTriplets(triplets.begin(), triplets.end());
    return matrix;
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {
    Vector<size_t> result = Vector<size_t>::Zero(mesh->nVertices());
    for (size_t i : subset.vertices) {
        result[i] = 1;
    }
    return result;
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {
    Vector<size_t> result = Vector<size_t>::Zero(mesh->nEdges());
    for (size_t i : subset.edges) {
        result[i] = 1;
    }
    return result;
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {
    Vector<size_t> result = Vector<size_t>::Zero(mesh->nFaces());
    for (size_t i : subset.faces) {
        result[i] = 1;
    }
    return result;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {
    // The star always includes the original subset.
    MeshSubset result = subset;

    Vector<size_t> vertices = buildVertexVector(subset);
    Vector<size_t> connectedEdges = A0 * vertices;

    for (size_t i = 0; i < connectedEdges.size(); i++) {
        if (connectedEdges[i] != 0) {
            result.addEdge(i);
        }
    }

    // Add edges from the original subset, so that we get the faces connected
    // to either a vertex (via one of the edges already in connectedEdges) or
    // an edge in the original subset.
    connectedEdges += buildEdgeVector(subset);
    Vector<size_t> connectedFaces = A1 * connectedEdges;

    for (size_t i = 0; i < connectedFaces.size(); i++) {
        if (connectedFaces[i] != 0) {
            result.addFace(i);
        }
    }

    return result;
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {
    // The closure always includes the original subset as well.
    MeshSubset result = subset;

    Vector<size_t> faces = buildFaceVector(subset);
    Vector<size_t> edges = A1.transpose() * faces;

    for (size_t i = 0; i < edges.size(); i++) {
        if (edges[i] != 0) {
            result.addEdge(i);
        }
    }

    edges += buildEdgeVector(subset);
    Vector<size_t> vertices = A0.transpose() * edges;

    for (size_t i = 0; i < vertices.size(); i++) {
        if (vertices[i] != 0) {
            result.addVertex(i);
        }
    }

    return result;
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {
    MeshSubset result = closure(star(subset));
    result.deleteSubset(star(closure(subset)));
    return result;
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {
    MeshSubset cl = closure(subset);
    return subset.equals(cl);
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {
    // We can skip the computation for the sets below, if there are no
    // corresponding simplices (in which case the empty set is correct).
    bool hasFaces = !subset.faces.empty();
    bool hasEdges = !subset.edges.empty();

    std::set<size_t> faceEdges;
    if (hasFaces) {
        Vector<size_t> faces = buildFaceVector(subset);
        Vector<size_t> edges = A1.transpose() * faces;
        for (size_t i = 0; i < edges.size(); i++) {
            if (edges[i] != 0) {
                faceEdges.insert(i);
            }
        }
    }

    std::set<size_t> edgeVertices;
    if (hasEdges) {
        Vector<size_t> edges = buildEdgeVector(subset);
        Vector<size_t> vertices = A0.transpose() * edges;
        for (size_t i = 0; i < vertices.size(); i++) {
            if (vertices[i] != 0) {
                edgeVertices.insert(i);
            }
        }
    }

    if (hasFaces) {
        // All edges in the subset must belong to a face and all vertices must belong to an edge (and transitively to a face).
        return (subset.edges == faceEdges && subset.vertices == edgeVertices) ? 2 : -1;
    } else if (hasEdges) {
        // All vertices must belong to an edge.
        return (subset.vertices == edgeVertices) ? 1 : -1;
    } else {
        // A set of vertices always is a pure 0-complex.
        // However, if the subset is completely empty, it is no longer a 0-complex,
        // which requires at least one simplex of degree 0 (i.e. at least one vertex).
        return (!subset.vertices.empty()) ? 0 : -1;
    }
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {
    int degree = isPureComplex(subset);

    // The boundary is defined as the closure of the set of all simplices that
    // are proper faces of exactly one simplex of the subset. In a pure
    // k-complex, this set will contain only (k-1)-simplices, so we switch on
    // the degree of the given complex.

    switch (degree) {
        case 0: {
            // A 0-complex does not have a boundary.
            return MeshSubset{};
        }

        case 1: {
            // For a 1-complex we have to find the vertices that are adjacent
            // to only one edge from the subset. The vertices vector below
            // will contain the number of adjacent edges. We do not have to
            // anything to get the closure, as the closure of a vertex is
            // always itself.
            MeshSubset result;

            Vector<size_t> edges = buildEdgeVector(subset);
            Vector<size_t> vertices = A0.transpose() * edges;

            for (size_t i = 0; i < vertices.size(); i++) {
                if (vertices[i] == 1) {
                    result.addVertex(i);
                }
            }

            return result;
        }

        case 2: {
            // The idea for 2-complexes is the same as for 1-complexes, except
            // now we are looking for edges that are adjacent to exactly one
            // triangle.
            MeshSubset result;

            Vector<size_t> faces = buildFaceVector(subset);
            Vector<size_t> edges = A1.transpose() * faces;

            for (size_t i = 0; i < edges.size(); i++) {
                if (edges[i] == 1) {
                    result.addEdge(i);
                } else {
                    edges[i] = 0;
                }
            }

            // Now, we also have to find the vertices corresponding to these
            // edges to get the closure of the set.
            Vector<size_t> vertices = A0.transpose() * edges;
            for (size_t i = 0; i < vertices.size(); i++) {
                if (vertices[i] != 0) {
                    result.addVertex(i);
                }
            }

            return result;

        }

        default:
            // The boundary is not defined if the subset is not a pure simplicial complex.
            // We return the original subset, so as not to clear the user's selection.
            return subset;
    }
}