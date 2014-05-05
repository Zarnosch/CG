#ifndef PLANE_HPP_INCLUDE_ONCE
#define PLANE_HPP_INCLUDE_ONCE

#include "raytracerConfig.hpp"

#include <memory>

#include "Renderable.hpp"
#include "Ray.hpp"

namespace rt
{

/// Infinite plane defined by point on plane and plane normal.
class Plane : public Renderable
{
public:
  RAYTRACER_EXPORTS Plane(const Vec3d& normal = Vec3d(0.0,0.0,1.0));

  /// Implements the intersection computation between ray and plane
  RAYTRACER_EXPORTS virtual bool
  closestIntersectionModel(const Ray &ray, double maxLambda, RayIntersection& intersection) const override;

  RAYTRACER_EXPORTS const Vec3d& normal() const { return mNormal; }

  RAYTRACER_EXPORTS void setNormal(const Vec3d &normal ) { mNormal=normal; mNormal.normalize(); }

protected:

  // Override this method to recompute the bounding box of this object.
  RAYTRACER_EXPORTS BoundingBox computeBoundingBox() const override;

private:
  Vec3d mNormal, mTangent, mBitangent;
};

}
#endif //PLANE_HPP_INCLUDE_ONCE
