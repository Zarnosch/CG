#include "CheckerMaterial.hpp"
#include "Intersection.hpp"

namespace rt
{

CheckerMaterial::CheckerMaterial(std::shared_ptr<Material> material1,
                                 std::shared_ptr<Material> material2,
                                 const Vec2d &tiles) : 
  Material(), mMaterial1(material1), mMaterial2(material2), mTiles(tiles)
{}

Vec4d CheckerMaterial::shade(const RayIntersection& intersection, 
           const Light& light) const
{
  const Vec3d &uvw = intersection.uvw();

  const bool left  = (fmod(fabs(uvw[0]), 1/mTiles[0])
    < (1/mTiles[0] / double(2)))
    ^ (uvw[0] < double(0));
  const bool lower = (fmod(fabs(uvw[1]), 1/mTiles[1])
    < (1/mTiles[1] / double(2)))
    ^ (uvw[1] < double(0));

  if(left ^ lower)
    return mMaterial1->shade(intersection, light);
  else
    return mMaterial2->shade(intersection, light);
}

} //namespace rt