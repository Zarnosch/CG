#include "Plane.hpp"
#include "Ray.hpp"
#include "Math.hpp"

namespace rt
{

Plane::Plane(const Vec3d& normal) :
  mNormal(normal)
{
  mNormal.normalize();

  mTangent = Vec3d(1,0,0);

  mBitangent = cross(mNormal,mTangent);

  const double bitangentNorm = mBitangent.length();

  if(bitangentNorm > Math::safetyEps())
  {
    mBitangent /= bitangentNorm;
    mTangent = cross(mBitangent,mNormal).normalize();
  } else
  {
    mTangent = Vec3d(0,1,0);
    mBitangent = Vec3d(0,0,1);
  }
}

bool
Plane::closestIntersectionModel(const Ray &ray, double maxLambda, RayIntersection& intersection) const
{
  // Solve linear equation for lambda
  // see also http://en.wikipedia.org/wiki/Line-plane_intersection
  double a=dot((-ray.origin()),(mNormal));
  double d=dot(ray.direction(),(mNormal));

  // No intersection if ray is (almost) parallel to plane
  if (fabs(d)<Math::safetyEps())
    return false;

  double lambda = a / d;

  // Only intersections in [0,1] range are valid.
  if (lambda<0.0 || lambda>maxLambda)
    return false;

  const Vec3d p = ray.pointOnRay(lambda);
  const Vec3d uvw(dot(p,mTangent),dot(p,mBitangent), double(0));

  intersection=RayIntersection(ray,shared_from_this(),lambda,mNormal,uvw);
  return true;
}

BoundingBox Plane::computeBoundingBox() const
{
  return BoundingBox(Vec3d(-std::numeric_limits<double>::infinity(),
                          -std::numeric_limits<double>::infinity(),
                          -std::numeric_limits<double>::infinity()),
                     Vec3d( std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::infinity()));
}

} //namespace rt
