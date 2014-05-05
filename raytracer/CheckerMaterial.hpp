#ifndef _CHECKERMATERIAL_HPP_INCLUDE_ONCE
#define _CHECKERMATERIAL_HPP_INCLUDE_ONCE

#include "raytracerConfig.hpp"

#include "Material.hpp"

namespace rt
{

/// Material that creates a checkerboard based on the uv parametrization
/// at the intersection point. It calls the shade function of either 
/// material1 o material2
class CheckerMaterial : public Material
{
public:
  RAYTRACER_EXPORTS CheckerMaterial(std::shared_ptr<Material> material1,
                  std::shared_ptr<Material> material2,
                  const Vec2d &tiles = Vec2d(1,1));

  RAYTRACER_EXPORTS Vec4d shade(const RayIntersection& intersection, 
    const Light& light) const override;

private:
  std::shared_ptr<Material> mMaterial1;
  std::shared_ptr<Material> mMaterial2;
  Vec2d mTiles; ///< number of tiles per uv in [0,1]x[0,1]
};
} //namespace rt

#endif //_CHECKERMATERIAL_HPP_INCLUDE_ONCE