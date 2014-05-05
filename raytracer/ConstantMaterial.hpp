#ifndef CONSTANTMATERIAL_HPP_INCLUDE_ONCE
#define CONSTANTMATERIAL_HPP_INCLUDE_ONCE

#include "raytracerConfig.hpp"

#include "Material.hpp"

namespace rt
{

class RayIntersection;
class Light;

/// Simple material that has no shading.
class ConstantMaterial : public Material
{
public:
  RAYTRACER_EXPORTS ConstantMaterial(const Vec3d& color = Vec3d(0,0.4,0.8));

  RAYTRACER_EXPORTS Vec4d shade(const RayIntersection& intersection, 
             const Light& light) const override;
};

} //namespace rt

#endif // CONSTANTMATERIAL_HPP_INCLUDE_ONCE
