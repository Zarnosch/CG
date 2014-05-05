#ifndef BVHINDEXEDTRIANGLEMESH_HPP_INCLUDE_ONCE
#define BVHINDEXEDTRIANGLEMESH_HPP_INCLUDE_ONCE

#include "raytracerConfig.hpp"

#include "IndexedTriangleMesh.hpp"
#include "BVTree.hpp"

namespace rt
{
class BVHIndexedTriangleMesh : public IndexedTriangleMesh
{
public:
  RAYTRACER_EXPORTS BVHIndexedTriangleMesh();

  RAYTRACER_EXPORTS void initialize() override;

  RAYTRACER_EXPORTS bool closestIntersectionModel(const Ray &ray, double maxLambda, RayIntersection& intersection) const override;

  RAYTRACER_EXPORTS bool anyIntersectionModel(const Ray &ray, double maxLambda) const override;

private:
  BVTree mTree;
};
} //namespace rt

#endif //BVHINDEXEDTRIANGLEMESH_HPP_INCLUDE_ONCE