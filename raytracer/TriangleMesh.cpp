#include "TriangleMesh.hpp"
#include "Intersection.hpp"

namespace rt
{


bool TriangleMesh::closestIntersectionModel(const Ray &ray, double maxLambda, RayIntersection& intersection) const
{
  double closestLambda = maxLambda;
  Vec3d closestbary;
  int  closestTri = -1;

  Vec3d bary;
  double lambda;

  // Loop over all stored triangles to find a possible intersection between ray 
  // and any triangle
  for (size_t i=0;i<mTriangles.size();++i)
  {
    const TriangleElement &tri = mTriangles[i];

    if (Intersection::lineTriangle(ray,tri.v0,tri.v1,tri.v2,bary,lambda) &&
      lambda > 0 && lambda < closestLambda)
    {
      closestLambda = lambda;
      closestbary = bary;
      closestTri = (int)i;
    }
  }

  // If an intersection occurred, compute the normal and the uv coordinates based
  // on the barycentric coordinates of the hit point
  if (closestTri >= 0)
  {
    const TriangleElement &tri = mTriangles[closestTri];
    Vec3d n = (tri.n0*closestbary[0]+tri.n1*closestbary[1]+
              tri.n2*closestbary[2]).normalize();
    Vec3d uvw = tri.uvw0*closestbary[0]+tri.uvw1*closestbary[1]+
               tri.uvw2*closestbary[2];

    intersection = RayIntersection(ray,shared_from_this(),closestLambda,n,uvw);
    return true;
  }
  return false;
}

bool TriangleMesh::anyIntersectionModel(const Ray &ray, double maxLambda) const
{
  Vec3d uvw;
  double lambda;
  for (size_t i=0;i<mTriangles.size();++i)
  {
    const TriangleElement &tri = mTriangles[i];
    if (Intersection::lineTriangle(ray,tri.v0,tri.v1,tri.v2,uvw,lambda) &&
      lambda > 0 && lambda < maxLambda)
      return true;
  }
  return false;
}

BoundingBox TriangleMesh::computeBoundingBox() const
{
  BoundingBox bbox;
  for (size_t i=0;i<mTriangles.size();++i)
  {
    const TriangleElement& tri = mTriangles[i];
    bbox.expandByPoint(tri.v0);
    bbox.expandByPoint(tri.v1);
    bbox.expandByPoint(tri.v2);
  }
  return bbox;
}

} //namespace rt