#ifndef _SPHERE_HPP_INCLUDE_ONCE
#define _SPHERE_HPP_INCLUDE_ONCE

#include "raytracerConfig.hpp"

#include "Renderable.hpp"
#include "Ray.hpp"
#include "Material.hpp"

namespace rt
{

class SphereRayIntersection;

//Unit sphere; i.e. sphere centered at origin with radius = 1
class Sphere : public Renderable
{
public:
  /// Implements the intersection computation between ray and sphere.
  RAYTRACER_EXPORTS bool
  closestIntersectionModel(const Ray &ray, double maxLambda, RayIntersection& intersection) const override;

  // Override this method to recompute the bounding box of this object.
  RAYTRACER_EXPORTS BoundingBox computeBoundingBox() const override;
};
} //namespace rt

#endif //_SPHERE_HPP_INCLUDE_ONCE
