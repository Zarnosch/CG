#ifndef _TRIANGLEMESH_HPP_HPP_INCLUDE_ONCE
#define _TRIANGLEMESH_HPP_HPP_INCLUDE_ONCE

#include "raytracerConfig.hpp"

#include "Renderable.hpp"
#include "Triangle.hpp"

namespace rt
{

/// Stores a collection of triangles (or rather TriangleElements)
/// Supports the intersection test of a ray with any of the stored
/// triangles ( linear time w.r.t. the number of triangles).
class TriangleMesh : public Renderable
{
private:

  /// Lightweight storage container for triangle data:
  /// vertex position, vertex normal and uv parameters
  RAYTRACER_EXPORTS struct TriangleElement
  {
    /// Create a triangle with constant normals over the triangle and
    /// and no parametrization.
    TriangleElement(const Vec3d &_v0,const Vec3d &_v1,const Vec3d &_v2) :
      v0(_v0),v1(_v1),v2(_v2)
    {
      n0 = n1 = n2 = cross(v1-v0,v2-v0).normalize();
      uvw0 = uvw1 = uvw2 = Vec3d(0,0,0);
    }

    /// Create a triangle with varying normals over the triangle and
    /// and no parametrization.
    TriangleElement(const Vec3d &_v0,const Vec3d &_v1,const Vec3d &_v2,
                    const Vec3d &_n0,const Vec3d &_n1,const Vec3d &_n2) :
    v0(_v0),v1(_v1),v2(_v2),n0(_n0),n1(_n1),n2(_n2)
    {
      uvw0 = uvw1 = uvw2 = Vec3d(0,0,0);
    }

    /// Create a triangle with varying normals over the triangle and
    /// and a uv parametrization.
    TriangleElement(const Vec3d &_v0,const Vec3d &_v1,const Vec3d &_v2,
                    const Vec3d &_n0,const Vec3d &_n1,const Vec3d &_n2,
                    const Vec3d &_uvw0,const Vec3d &_uvw1,const Vec3d &_uvw2) :
      v0(_v0),v1(_v1),v2(_v2),
      n0(_n0),n1(_n1),n2(_n2),
      uvw0(_uvw0),uvw1(_uvw1),uvw2(_uvw2)
    {
      
    }
    Vec3d v0;
    Vec3d v1;
    Vec3d v2;
    Vec3d n0;
    Vec3d n1;
    Vec3d n2;
    Vec3d uvw0;
    Vec3d uvw1;
    Vec3d uvw2;
  };

public:
  
  /// Implements the intersection computation between ray and any stored triangle.
  RAYTRACER_EXPORTS bool
    closestIntersectionModel(const Ray &ray, double maxLambda, RayIntersection& intersection) const override;

  RAYTRACER_EXPORTS bool anyIntersectionModel(const Ray &ray, double maxLambda) const override;

  RAYTRACER_EXPORTS void addTriangle(const Vec3d &v0,const Vec3d &v1,const Vec3d &v2)
  {
    mTriangles.push_back(TriangleElement(v0,v1,v2));
  }
  RAYTRACER_EXPORTS void addTriangle(const Vec3d &v0,const Vec3d &v1,const Vec3d &v2,
                   const Vec3d &n0,const Vec3d &n1,const Vec3d &n2)
  {
    mTriangles.push_back(TriangleElement(v0,v1,v2,n0,n1,n2));
  }
  RAYTRACER_EXPORTS void addTriangle(const Vec3d &v0,const Vec3d &v1,const Vec3d &v2,
                   const Vec3d &n0,const Vec3d &n1,const Vec3d &n2,
                   const Vec3d &uvw0,const Vec3d &uvw1,const Vec3d &uvw2)
  {
    mTriangles.push_back(TriangleElement(v0,v1,v2,n0,n1,n2,uvw0,uvw1,uvw2));
  }

protected:

  // Override this method to recompute the bounding box of this object.
  RAYTRACER_EXPORTS BoundingBox computeBoundingBox() const override;

private:
  std::vector<TriangleElement> mTriangles;

};
} //namespace rt

#endif //_TRIANGLEMESH_HPP_HPP_INCLUDE_ONCE