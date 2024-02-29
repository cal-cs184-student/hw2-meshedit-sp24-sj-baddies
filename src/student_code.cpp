#include "student_code.h"
#include "CGL/vector2D.h"
#include "CGL/vector3D.h"
#include "halfEdgeMesh.h"
#include "mutablePriorityQueue.h"
#include <vector>

using namespace std;

namespace CGL
{

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (class member).
   *
   * @param points A vector of points in 2D
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector2D> BezierCurve::evaluateStep(std::vector<Vector2D> const &points)
  { 
    std::vector<Vector2D> intermediate_points;
    for(int i = 0; i < points.size() - 1; i++){
      Vector2D p_i = ((1 - t) * points.at(i)) + t * points.at(i + 1);
      intermediate_points.push_back(p_i);
    }
    return intermediate_points;
  }

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (function parameter).
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector3D> BezierPatch::evaluateStep(std::vector<Vector3D> const &points, double t) const
  {
    std::vector<Vector3D> intermediate_points;
    for(int i = 0; i < points.size() - 1; i++){
      Vector3D p_i = ((1 - t) * points.at(i)) + t * points.at(i + 1);
      intermediate_points.push_back(p_i);
    }
    return intermediate_points;
  }

  /**
   * Fully evaluates de Casteljau's algorithm for a vector of points at scalar parameter t
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> const &points, double t) const
  {
    if(points.size() == 1){
      return points.at(0);
    }
    return evaluate1D(evaluateStep(points, t), t);
  }

  /**
   * Evaluates the Bezier patch at parameter (u, v)
   *
   * @param u         Scalar interpolation parameter
   * @param v         Scalar interpolation parameter (along the other axis)
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate(double u, double v) const 
  {  
    std::vector<Vector3D> intermediate_points;
    for(int i = 0; i < controlPoints.size(); i++){
      intermediate_points.push_back(evaluate1D(controlPoints[i], u));
    }
    return evaluate1D(intermediate_points, v);
  }

  Vector3D Vertex::normal( void ) const
  {
    // Returns an approximate unit normal at this vertex, computed by
    // taking the area-weighted average of the normals of neighboring
    // triangles, then normalizing.

    Vector3D weighted = Vector3D(0, 0, 0);
    
    HalfedgeCIter h = this->halfedge();
    FaceCIter f = h->face();
    // 
    do {
      HalfedgeCIter h_twin = h->twin();
      VertexCIter v = h_twin->vertex();

      weighted += h->face()->normal();

      h = h_twin->next();

    } while(h != this->halfedge());

    return weighted.unit();
  }

  EdgeIter HalfedgeMesh::flipEdge( EdgeIter e0 )
  {
    // This method should flip the given edge and return an iterator to the flipped edge.

    if (e0->isBoundary()) {
        return e0;
    }

    // PHASE 1
    HalfedgeIter h0 = e0->halfedge();
    HalfedgeIter h1 = h0->next();
    HalfedgeIter h2 = h1->next();
    HalfedgeIter h3 = h0->twin();
    HalfedgeIter h4 = h3->next();
    HalfedgeIter h5 = h4->next();
    HalfedgeIter h6 = h1->twin();
    HalfedgeIter h7 = h2->twin();
    HalfedgeIter h8 = h4->twin();
    HalfedgeIter h9 = h5->twin();

    VertexIter v0 = h0->vertex();
    VertexIter v1 = h3->vertex();
    VertexIter v2 = h6->vertex();
    VertexIter v3 = h8->vertex();

    EdgeIter e1 = h1->edge();
    EdgeIter e2 = h2->edge();
    EdgeIter e3 = h4->edge();
    EdgeIter e4 = h5->edge();

    FaceIter f0 = h0->face();
    FaceIter f1 = h3->face();

    h0->setNeighbors(h1, h3, v3, e0, f0);
    h1->setNeighbors(h2, h7, v2, e2, f0);
    h2->setNeighbors(h0, h8, v0, e3, f0);
    h3->setNeighbors(h4, h0, v2, e0, f1);
    h4->setNeighbors(h5, h9, v3, e4, f1);
    h5->setNeighbors(h3, h6, v1, e1, f1);
    h6->setNeighbors(h6->next(), h5, v2, e1, h6->face());
    h7->setNeighbors(h7->next(), h1, v0, e2, h7->face());
    h8->setNeighbors(h8->next(), h2, v3, e3, h8->face());
    h9->setNeighbors(h9->next(), h4, v1, e4, h9->face());

    v0->halfedge() = h2;
    v1->halfedge() = h5;
    v2->halfedge() = h3;
    v3->halfedge() = h0;

    e0->halfedge() = h0;
    e1->halfedge() = h6;
    e2->halfedge() = h7;
    e3->halfedge() = h8;
    e4->halfedge() = h9;

    f0->halfedge() = h0;   
    f1->halfedge() = h3;
    return e0;
  }

  VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
  {
    // This method should split the given edge and return an iterator to the newly inserted vertex.
    // The halfedge of this vertex should point along the edge that was split, rather than the new edges.

    if (e0->isBoundary()) {
        return e0->halfedge()->vertex();
    }

    // COLLECT =========================
    HalfedgeIter h0 = e0->halfedge();
    HalfedgeIter h1 = h0->next();
    HalfedgeIter h2 = h1->next();
    HalfedgeIter h3 = h0->twin();
    HalfedgeIter h4 = h3->next();
    HalfedgeIter h5 = h4->next();
    HalfedgeIter h6 = h1->twin();
    HalfedgeIter h7 = h2->twin();
    HalfedgeIter h8 = h4->twin();
    HalfedgeIter h9 = h5->twin();

    VertexIter v0 = h0->vertex();
    VertexIter v1 = h3->vertex();
    VertexIter v2 = h6->vertex();
    VertexIter v3 = h8->vertex();

    EdgeIter e1 = h1->edge();
    EdgeIter e2 = h2->edge();
    EdgeIter e3 = h4->edge();
    EdgeIter e4 = h5->edge();

    FaceIter f0 = h0->face();
    FaceIter f1 = h3->face();
    // COLLECT =========================

    // ALLOCATE NEW HALFEDGES
    HalfedgeIter h10 = newHalfedge();
    HalfedgeIter h11 = newHalfedge();
    HalfedgeIter h12 = newHalfedge();
    HalfedgeIter h13 = newHalfedge();
    HalfedgeIter h14 = newHalfedge();
    HalfedgeIter h15 = newHalfedge();
    // ALLOCATE NEW EDGES
    EdgeIter e5 = newEdge();
    EdgeIter e6 = newEdge();
    EdgeIter e7 = newEdge();
    // ALLOCATE NEW VERTEX
    VertexIter v4 = newVertex();
    // ALLOCATE NEW FACES
    FaceIter f2 = newFace();
    FaceIter f3 = newFace();

    // REASSIGN ===========================
    // next, twin, vertex, edge, face
    h0->setNeighbors(h13, h3, v0, e0, f2);
    h1->setNeighbors(h12, h6, v1, e1, f1);
    h2->setNeighbors(h0, h7, v2, e2, f2);
    h3->setNeighbors(h4, h0, v4, e0, f3);
    h4->setNeighbors(h11, h8, v0, e3, f3);
    h5->setNeighbors(h14, h9, v3, e4, f0);

    // h6->setNeighbors(h9, h1, v2, e1, f1);
    // h7->setNeighbors(h6, h2, v0, e2, f2);
    // h8->setNeighbors(h7, h4, v3, e3, f3);
    // h9->setNeighbors(h8, h5, v1, e4, f0);

    h6->setNeighbors(h6->next(), h1, v2, e1, h6->face());
    h7->setNeighbors(h7->next(), h2, v0, e2, h7->face());
    h8->setNeighbors(h8->next(), h4, v3, e3, h8->face());
    h9->setNeighbors(h9->next(), h5, v1, e4, h9->face());
    
    h10->setNeighbors(h5, h11, v4, e5, f0);
    h11->setNeighbors(h3, h10, v3, e5, f3);
    h12->setNeighbors(h15, h13, v2, e6, f1);
    h13->setNeighbors(h2, h12, v4, e6, f2);
    h14->setNeighbors(h10, h15, v1, e7, f0);
    h15->setNeighbors(h1, h14, v4, e7, f1);

    v0->halfedge() = h0;
    v1->halfedge() = h14;
    v2->halfedge() = h12;
    v3->halfedge() = h11;
    v4->halfedge() = h3;

    e0->halfedge() = h0;
    e1->halfedge() = h6;
    e2->halfedge() = h7;
    e3->halfedge() = h8;
    e4->halfedge() = h9;
    e5->halfedge() = h10;
    e6->halfedge() = h12;
    e7->halfedge() = h15;

    f0->halfedge() = h5;   
    f1->halfedge() = h1;
    f2->halfedge() = h2;
    f3->halfedge() = h4;
    // REASSIGN ===========================
    v4->position = 0.5 * (v0->position + v1->position);
    return v4;
  }

  void MeshResampler::upsample( HalfedgeMesh& mesh )
  {
    // TODO Part 6.
    // This routine should increase the number of triangles in the mesh using Loop subdivision.
    // One possible solution is to break up the method as listed below.

    // 1. Compute new positions for all the vertices in the input mesh, using the Loop subdivision rule,
    // and store them in Vertex::newPosition. At this point, we also want to mark each vertex as being
    // a vertex of the original mesh.

    for (VertexIter vertex = mesh.verticesBegin(); vertex != mesh.verticesEnd(); vertex++){
      
      Vector3D original_neighbor_position_sum = Vector3D();
      HalfedgeCIter halfedge = vertex->halfedge();

      do { // loop neighbors to sum neighbor positions
        HalfedgeCIter twin = halfedge->twin();
        VertexCIter neighbor = twin->vertex();
        original_neighbor_position_sum += neighbor->position;
        halfedge = twin->next();
      } while(halfedge != vertex->halfedge()); 

      float n = vertex->degree();
      float u;
      if (n ==3) {
        u = 3 / 16.0;
      }
      else {
        u = 3 / (8 * n);
      }
      // (1 - n * u) * original_position + u * original_neighbor_position_sum
      vertex->newPosition = (1 - n * u) * vertex->position + u * original_neighbor_position_sum;
      vertex->isNew = false;
    }
    
    // 2. Compute the updated vertex positions associated with edges, and store it in Edge::newPosition.
    for (EdgeIter edge = mesh.edgesBegin(); edge != mesh.edgesEnd(); edge++) {

      // 3/8 * (A + B) + 1/8 * (C + D)
      HalfedgeIter halfedge = edge->halfedge();
      HalfedgeIter twin = halfedge->twin();

      Vector3D vertex_A = halfedge->vertex()->position;
      Vector3D vertex_B = twin->vertex()->position;
      Vector3D vertex_C = halfedge->next()->next()->vertex()->position;
      Vector3D vertex_D = twin->next()->next()->vertex()->position;

      edge->newPosition = (3/8.0) * (vertex_A + vertex_B) + (1/8.0) * (vertex_C + vertex_D);
      edge->isNew = false;
    }
    
    // 3. Split every edge in the mesh, in any order. For future reference, we're also going to store some
    // information about which subdivide edges come from splitting an edge in the original mesh, and which edges
    // are new, by setting the flat Edge::isNew. Note that in this loop, we only want to iterate over edges of
    // the original mesh---otherwise, we'll end up splitting edges that we just split (and the loop will never end!)
    for (EdgeIter edge = mesh.edgesBegin(); edge != mesh.edgesEnd(); edge++) {
      if(!edge->halfedge()->vertex()->isNew && !edge->halfedge()->twin()->vertex()->isNew){ // only split old edges
        VertexIter split = mesh.splitEdge(edge);
        HalfedgeIter halfedge = split->halfedge();
        split->newPosition = halfedge->edge()->newPosition;
        halfedge->next()->next()->edge()->isNew = true;
        halfedge->twin()->next()->edge()->isNew = true;
        split->isNew = true;
      }
    }
    
    // 4. Flip any new edge that connects an old and new vertex.
    for (EdgeIter edge = mesh.edgesBegin(); edge != mesh.edgesEnd(); edge++) {
      if (edge->isNew) {
        bool split_new = edge->halfedge()->vertex()->isNew;
        bool split_twin_new = edge->halfedge()->twin()->vertex()->isNew;
        if ((split_new && !split_twin_new) || (!split_new && split_twin_new))
          edge = mesh.flipEdge(edge);
      }
    }

    // 5. Copy the new vertex positions into final Vertex::position.
    for (VertexIter vertex = mesh.verticesBegin(); vertex != mesh.verticesEnd(); vertex++) {
      vertex->position = vertex->newPosition;
    }
  }
}
