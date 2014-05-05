#include "Triangle.hpp"
#include "Ray.hpp"
#include <memory>
#include "Intersection.hpp"

namespace rt
{
Triangle::Triangle(const Vec3d &v0, const Vec3d &v1, const Vec3d &v2,
                   const Vec3d &uvw0, const Vec3d &uvw1, const Vec3d &uvw2)
{
  mVertices[0] = v0;
  mVertices[1] = v1;
  mVertices[2] = v2;
  mUVW[0] = uvw0;
  mUVW[1] = uvw1;
  mUVW[2] = uvw2;
}

bool
Triangle::closestIntersectionModel(const Ray &ray, double maxLambda, RayIntersection& intersection) const
{
  Vec3d bary;
  double lambda;

  if(!Intersection::lineTriangle(ray,mVertices[0],mVertices[1],mVertices[2],bary,lambda))
    return false;

  // Intersection is inside triangle if 0<=u,v,w<=1
  if(lambda<0 || lambda>maxLambda)
    return false;

  const Vec3d normal = cross(mVertices[1]-mVertices[0],mVertices[2]-mVertices[0]).normalize();
  const Vec3d uvw = mUVW[0]*bary[0]+mUVW[1]*bary[1]+mUVW[2]*bary[2];

  intersection = RayIntersection(ray,shared_from_this(),lambda,normal,uvw);
  return true;                      
}

BoundingBox Triangle::computeBoundingBox() const
{
  BoundingBox bbox;
  bbox.expandByPoint(mVertices[0]);
  bbox.expandByPoint(mVertices[1]);
  bbox.expandByPoint(mVertices[2]);
  return bbox;
}

} //namespace rt
