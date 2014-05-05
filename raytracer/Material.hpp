#ifndef MATERIAL_HPP_INCLUDE_ONCE
#define MATERIAL_HPP_INCLUDE_ONCE

#include "raytracerConfig.hpp"

#include "Math.hpp"
#include <memory>

namespace rt
{
class RayIntersection;
class Light;

/// Abstract material with shade and BRDF functions.
class Material
{
public:
  RAYTRACER_EXPORTS Material(const Vec3d &color = Vec3d(0.5,0.5,0.5), double reflectance=0.0) :
    mColor(color), mReflectance(reflectance) {}

  RAYTRACER_EXPORTS virtual ~Material() {}

  // Returns the RGBA color for a point viewed from a direction 
  // lit by a light.
  RAYTRACER_EXPORTS virtual Vec4d shade(const RayIntersection& intersection, 
                     const Light& light) const = 0;

  RAYTRACER_EXPORTS const Vec3d& color() const {return mColor;}

  RAYTRACER_EXPORTS double reflectance() const {
    return mReflectance;
  }

  /// Valid RGB color components have range [0,1].
  RAYTRACER_EXPORTS void setColor(const Vec3d& color) { mColor=color; }

private:
Vec3d mColor;
double mReflectance;
};

}

#endif //MATERIAL_HPP_INCLUDE_ONCE
