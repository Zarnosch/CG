#ifndef DIFFUSEMATERIAL_HPP_INCLUDE_ONCE
#define DIFFUSEMATERIAL_HPP_INCLUDE_ONCE

#include "raytracerConfig.hpp"

#include "Material.hpp"

namespace rt
{

class RayIntersection;
class Light;

/// Material based on Lambert's cosine law
/// see http://en.wikipedia.org/wiki/Lambert%27s_cosine_law
class DiffuseMaterial : public Material
{
public:
  RAYTRACER_EXPORTS DiffuseMaterial(const Vec3d& color = Vec3d(0,0.4,0.8));

  RAYTRACER_EXPORTS Vec4d shade(const RayIntersection& intersection, 
             const Light& light) const override;
};

} //namespace rt

#endif // DIFFUSEMATERIAL_HPP_INCLUDE_ONCE
