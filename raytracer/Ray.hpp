#ifndef RAY_HPP_INCLUDE_ONCE
#define RAY_HPP_INCLUDE_ONCE

#include "raytracerConfig.hpp"

#include <memory>

#include "Math.hpp"

namespace rt
{

class Renderable;

/// Ray consists of point and direction
class Ray
{
public:
  RAYTRACER_EXPORTS Ray() : mOrigin(0.0,0.0,0.0), mDirection(0.0,0.0,0.0) {}
  RAYTRACER_EXPORTS Ray(const Vec3d &origin, const Vec3d &direction) :
    mOrigin(origin), mDirection(Vec3d(direction).normalize())  {}

  RAYTRACER_EXPORTS Vec3d pointOnRay(double lambda) const { return mOrigin + mDirection*lambda; }

  RAYTRACER_EXPORTS const Vec3d& origin()    const { return mOrigin; }
  RAYTRACER_EXPORTS const Vec3d& direction() const { return mDirection; }

  RAYTRACER_EXPORTS void setOrigin(const Vec3d& origin)       { mOrigin=origin; }
  RAYTRACER_EXPORTS void setDirection(const Vec3d& direction) { mDirection=Vec3d(direction).normalize(); }

  RAYTRACER_EXPORTS void transform(const Mat4x4d &transform)
  {
    mOrigin    = transform*mOrigin;
    mDirection = (transform.as3x3()*mDirection).normalize();
  }

  RAYTRACER_EXPORTS Ray transformed(const Mat4x4d &transform) const
  {
    return Ray(transform*mOrigin,
               (transform.as3x3()*mDirection).normalize());
  }

private:
  Vec3d mOrigin;
  Vec3d mDirection;
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/// Container class for intersection between a renderable and a ray
class RayIntersection
{
public:

  RAYTRACER_EXPORTS RayIntersection(){}

  RAYTRACER_EXPORTS RayIntersection(const Ray &ray,
                  std::shared_ptr<const Renderable> renderable,
                  const double lambda, const Vec3d &normal, const Vec3d &uvw) :
    mRay(ray), mRenderable(renderable), mLambda(lambda), mNormal(normal), mUVW(uvw)
  {
    mPosition=ray.pointOnRay(mLambda);
  }

  RAYTRACER_EXPORTS virtual ~RayIntersection() {}

  RAYTRACER_EXPORTS const Ray& ray()                               const { return mRay; }
  RAYTRACER_EXPORTS std::shared_ptr<const Renderable> renderable() const { return mRenderable; }
  RAYTRACER_EXPORTS double lambda()                                  const { return mLambda; }
  RAYTRACER_EXPORTS const Vec3d& position()                         const { return mPosition; }
  RAYTRACER_EXPORTS const Vec3d& normal()                           const { return mNormal; }
  RAYTRACER_EXPORTS const Vec3d& uvw()                              const { return mUVW; }

  RAYTRACER_EXPORTS virtual void transform(const Mat4x4d &transform,
                         const Mat4x4d &transformInvTransp)
  {

    //transform local intersection information to global system
    mPosition = transform*mPosition;
    mNormal   = (transformInvTransp.as3x3()*mNormal).normalize();
    mRay.transform(transform);
    mLambda   = (mPosition - mRay.origin()).length();
  }

protected:
  Ray mRay;
  std::shared_ptr<const Renderable> mRenderable;
  double mLambda;
  Vec3d mPosition;
  Vec3d mNormal;
  Vec3d mUVW;
};

} //namespace rt

#endif //RAY_HPP_INCLUDE_ONCE
