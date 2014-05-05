#include "ConstantMaterial.hpp"

namespace rt {

ConstantMaterial::ConstantMaterial(const Vec3d& color) :  
  Material(color)
{

}

Vec4d ConstantMaterial::shade(const RayIntersection& intersection,
                             const Light& light) const 
{
   return Vec4d(this->color(),1.0);
}

} //namespace rt
