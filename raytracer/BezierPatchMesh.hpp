#ifndef _BEZIERPATCHMESH_HPP_HPP_INCLUDE_ONCE
#define _BEZIERPATCHMESH_HPP_HPP_INCLUDE_ONCE

#include "raytracerConfig.hpp"

#include "BVHIndexedTriangleMesh.hpp"

namespace rt
{

// Bezier surface stores the set of control points and allows
// creation of a triangle mesh by sampling the parametric surface
class BezierPatchMesh : public BVHIndexedTriangleMesh
{
public: 
  /**
   * @brief BezierMeshPatch constructor
   * @param m    number of control points in u direction
   * @param n    number of control points in v direction
   * @param resu triangle resolution in u direction
   * @param resv triangle resolution in v direction
   */
  RAYTRACER_EXPORTS BezierPatchMesh(size_t m,    size_t n,
                  size_t resu, size_t resv);
  float BezierPatchMesh::factorial(float n)const;

  // Must be called before rendering and after control point manipulation
  // Creates the set of triangles.
  RAYTRACER_EXPORTS void initialize();

  RAYTRACER_EXPORTS void setControlPoint(size_t i, size_t j, const Vec3d& p)
  {
    mControlPoints[mM*j + i] = p;
  }
  RAYTRACER_EXPORTS const Vec3d& controlPoint(size_t i, size_t j) const
  {
    return mControlPoints[mM*j + i];
  }

  RAYTRACER_EXPORTS struct BezierPatchSample
  {
    Vec3d position;
    Vec3d normal;
    Vec2d uv;
  };
  
  /** Samples the Bezier surface at the given point.
   * @brief sample
   * @param u coordinate of domain point
   * @param v coordinate of domain point
   * @return BezierPatchSample sample
   */
   RAYTRACER_EXPORTS BezierPatchSample sample(double u, double v) const;

protected:
    // Override this method to recompute the bounding box of this object.
    RAYTRACER_EXPORTS BoundingBox computeBoundingBox() const override;

private:

  size_t mM, mN;                    //!< patch control point dimensions
  size_t mResU, mResV;              //!< triangle resolution in both parameter directions
  std::vector<Vec3d> mControlPoints; //!< patch control points

};
} //namespace rt

#endif // _BEZIERPATCHMESH_HPP_HPP_INCLUDE_ONCE
