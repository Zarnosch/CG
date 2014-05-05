#include "Renderable.hpp"
#include "Ray.hpp"
namespace rt
{

Renderable::Renderable() : mTransformClean(true)
{

}

Renderable::~Renderable()
{

}

bool Renderable::anyIntersectionModel(const Ray &ray, double maxLambda) const
{
  RayIntersection intersection;
  return this->closestIntersectionModel(ray,maxLambda,intersection);
}

bool
Renderable::closestIntersection(const Ray &ray, double maxLambda, RayIntersection& intersection) const
{
  //adapt maximal lambda value relative to transformation properties
  //(e.g., scaling)
  Ray  modelRay       = transformRayWorldToModel(ray);
  maxLambda = transformRayLambdaWorldToModel(ray, maxLambda);

  // early out test - do we hit the bounding box of the object?
  if (!mBoundingBox.anyIntersection(modelRay, maxLambda))
    return false;

  //transform ray from world to model coordinate system
  RayIntersection isect_local;
  if (!this->closestIntersectionModel(modelRay,maxLambda,isect_local))
    return false;

  //transform intersection from model to world coordinate system
  isect_local.transform(mTransform, mTransformInvTransp);
  intersection=isect_local;
  return true;
}

bool Renderable::anyIntersection(const Ray &ray, double maxLambda) const
{
  Ray  modelRay       = transformRayWorldToModel(ray);
  maxLambda = transformRayLambdaWorldToModel(ray, maxLambda);
  // early out test - do we hit the bounding box of the object?
  if (!mBoundingBox.anyIntersection(modelRay,maxLambda))
    return false;

  return this->anyIntersectionModel(modelRay,maxLambda);
}

void Renderable::updateTransforms() 
{
  if(!mTransformClean)
  {
    mTransformInv = mTransform;

    // The matrix should always be invertible. If this wasn't the case
    // an error should be thrown.
    double det;
    mTransformInv.invert(&det);
    if(det==0)
    {
      std::cerr<<"Renderable::transformRayLocal: Error: tranformation not invertible"<<std::endl;
    }
    mTransformInvTransp = mTransformInv.getTransposed();
    mTransformClean = true;
  }
}

Ray Renderable::transformRayWorldToModel(const Ray &ray) const
{
  return ray.transformed(mTransformInv);
}

double Renderable::transformRayLambdaWorldToModel(const Ray &ray, const double lambda) const
{
  Vec3d model_direction = mTransformInv.as3x3()*ray.direction();

  return lambda * model_direction.length();
}

}
