#ifndef _TRIANGLE_HPP_INCLUDE_ONCE
#define _TRIANGLE_HPP_INCLUDE_ONCE

#include "raytracerConfig.hpp"

#include "Renderable.hpp"
#include "Ray.hpp"
#include "Material.hpp"

namespace rt
{

class TriangleRayIntersection;

/// Triangle defined by position and radius.
class Triangle : public Renderable
{
public:
  RAYTRACER_EXPORTS Triangle( const Vec3d &v0, const Vec3d &v1, const Vec3d &v2,
            const Vec3d &uvw0=Vec3d(0,0,0), const Vec3d &uvw1=Vec3d(0,0,0),
            const Vec3d &uvw2=Vec3d(0,0,0));

  /// Implements the intersection computation between ray and Triangle.
  RAYTRACER_EXPORTS bool
  closestIntersectionModel(const Ray &ray, double maxLambda, RayIntersection& intersection) const override;

protected:

  // Override this method to recompute the bounding box of this object.
  RAYTRACER_EXPORTS BoundingBox computeBoundingBox() const override;

private:

  // Vertex positions
  Vec3d mVertices[3];
  Vec3d mUVW[3];
};
} //namespace rt

#endif //_Triangle_HPP_INCLUDE_ONCE
