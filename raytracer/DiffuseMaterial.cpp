#include "DiffuseMaterial.hpp"
#include "Light.hpp"
#include "Math.hpp"
#include "Ray.hpp"

namespace rt {

DiffuseMaterial::DiffuseMaterial(const Vec3d& color) :  
  Material(color)
{

}

Vec4d DiffuseMaterial::shade(const RayIntersection& intersection,
                             const Light& light) const 
{
  Vec3d N = intersection.normal();
  Vec3d L = (light.position() - intersection.position()).normalize();

  double cosNL = std::max(dot(N,L),double(0));

  return Vec4d(this->color()*cosNL,1.0);
}

} //namespace rt
