#include "Sphere.hpp"
#include "Ray.hpp"
#include <memory>

namespace rt
{

bool
Sphere::closestIntersectionModel(const Ray &ray, double maxLambda, RayIntersection& intersection) const
{
  // Implicit sphere: (x-c)^2-r^2 = 0
  // r = 1, c=0
  // x^2 = 1
  // Ray: x = o+t*d
  
  // Solve: 0=(o+t*d)^2 -1
  // 0 = t^2*(d*d) + t*(2*o*d) + (o*o-1)

  // Then, these are the 2 possible solutions
  // t0 = 0.5*(-b-sqrt(b*b - 4*a))

  // simplify
  // t01 = -o.d +- sqrt( (o.d)^2 - o.o + 1)

  const Vec3d &d = ray.direction();
  const Vec3d &o = ray.origin();

  double b = dot(o,d);
  double c = o.lengthSquared() - 1;

  double ds = b*b - c;

  //discriminant is negative, no intersection
  if(ds < 0)
    return false;

  double dssqr = sqrt(ds);

  double t0 = -b-dssqr;
  double t1 = -b+dssqr;

  double lambda = t0;

  //ray starts inside sphere, discard first hit
  if(t0 < 0)
    lambda = t1;

  if(lambda < 0 || lambda > maxLambda)
    return false;

  // Compute intersection point
  const Vec3d p = ray.pointOnRay(lambda);

  // Compute parameterization of intersection point
  const double theta = std::atan2(p[1], p[0]);
  const double phi   = std::acos (p[2]);
  const Vec3d uvw(theta,phi,double(0));

  intersection = RayIntersection(ray,shared_from_this(),lambda,ray.pointOnRay(lambda),uvw);
  return true;
}

BoundingBox Sphere::computeBoundingBox() const
{
  BoundingBox box;
  box.setMin(Vec3d(-1,-1,-1));
  box.setMax(Vec3d(1,1,1));
  return box;
}

} //namespace rt
