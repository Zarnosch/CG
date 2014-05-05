#ifndef RENDERABLE_HPP_INCLUDE_ONCE
#define RENDERABLE_HPP_INCLUDE_ONCE

#include "raytracerConfig.hpp"

#include <memory>
#include <vector>
#include "BoundingBox.hpp"
#include "Math.hpp"

namespace rt
{
class Material;
class Ray;
class RayIntersection;

/// Abstract class for visible geometry.
class Renderable : public std::enable_shared_from_this<Renderable>
{
public:
  // Constructor.
  RAYTRACER_EXPORTS Renderable();
  // Destructor.
  RAYTRACER_EXPORTS virtual ~Renderable();

  // Computes the point of intersection between ray and object.
  RAYTRACER_EXPORTS bool
    closestIntersection(const Ray &ray, double maxLambda, RayIntersection& intersection) const;

  // This is used for so-called 'any hit' rays (returns true if there is at
  // least one intersection.)
  RAYTRACER_EXPORTS bool anyIntersection(const Ray &ray, double maxLambda) const;

  // Gets the transformation.
  RAYTRACER_EXPORTS Mat4x4d& transform()
  {
    mTransformClean = false;
    return mTransform;
  }

  // Gets the material.
  RAYTRACER_EXPORTS std::shared_ptr<const Material> material() const { return mMaterial; }
  // Sets the material.
  RAYTRACER_EXPORTS void setMaterial(std::shared_ptr<Material> material) { mMaterial = material; }

  // Recomputes the bounding box.
  RAYTRACER_EXPORTS void updateBoundingBox() { mBoundingBox = this->computeBoundingBox();}

  // Override this method for pre-render initialization
  RAYTRACER_EXPORTS virtual void initialize() {} 

  //Update the internal cached matrices based on the transformation.
  void updateTransforms();

protected:

  // This function does the ray intersection test in the local model coordinate
  // system of the object. Override this function for each class inheriting
  // from Renderable.
  RAYTRACER_EXPORTS virtual bool
    closestIntersectionModel(const Ray &ray, double maxLambda,RayIntersection& intersection) const = 0;

  // This function does an any hit ray intersection test in the local model
  // coordinate system of the object. By default, closestIntersectionLocal
  // is called. For most geometric primitives, there are faster methods to
  // test if there is at least one intersection rather than computing the
  // closest. If so, override this method and perform the faster test.
  RAYTRACER_EXPORTS virtual bool anyIntersectionModel(const Ray &ray, double maxLambda) const;

  // Override this method to recompute the bounding box of this object.
  RAYTRACER_EXPORTS virtual BoundingBox computeBoundingBox() const = 0;

private:
  Mat4x4d mTransform;
  std::shared_ptr<Material> mMaterial;

  bool mTransformClean;
  Mat4x4d mTransformInv, mTransformInvTransp;
  Ray  transformRayWorldToModel(const Ray &ray) const;
  double transformRayLambdaWorldToModel(const Ray &ray, const double lambda) const;
  BoundingBox mBoundingBox;
};

} //namespace rt

#endif //RENDERABLE_HPP_INCLUDE_ONCE
