#ifndef CAMERA_HPP_INCLUDE_ONCE
#define CAMERA_HPP_INCLUDE_ONCE

#include "raytracerConfig.hpp"

#include "Ray.hpp"

namespace rt {

/// Abstract base class for camera.
class Camera
{
public:
  RAYTRACER_EXPORTS Camera();
  RAYTRACER_EXPORTS virtual ~Camera();

  /// Compute the primary ray passing through pixel x,y.
  RAYTRACER_EXPORTS virtual Ray ray(size_t x, size_t y) const = 0;

  // Constant accessors.
  RAYTRACER_EXPORTS const Vec3d& position()     const { return mPosition; }
  RAYTRACER_EXPORTS const Vec3d& lookAt()       const { return mLookAt; }
  RAYTRACER_EXPORTS const Vec3d& up()           const { return mUp; }
  RAYTRACER_EXPORTS double horizontalFOV() const { return mHorizontalFOV; }
  RAYTRACER_EXPORTS double verticalFOV()   const { return mVerticalFOV; }
  RAYTRACER_EXPORTS size_t xResolution()       const { return mXResolution; }
  RAYTRACER_EXPORTS size_t yResolution()       const { return mYResolution; }

  // Mutator methods requiring reinitialization.

  RAYTRACER_EXPORTS void setPosition(const Vec3d &position)
  {
    mPosition=position;
    this->init();
  }

  /// Camera target point must not coincide with camera eye position.
  RAYTRACER_EXPORTS void setLookAt(const Vec3d &lookAt)
  {
    mLookAt=lookAt;
    this->init();
  }

  /// Camera target point must not coincide with direction from camera
  /// target to eye.
  RAYTRACER_EXPORTS void setUp(const Vec3d &up)
  {
    mUp=up;
    this->init();
  }

  RAYTRACER_EXPORTS void setFOV(double horizontalFOV, double verticalFOV)
  {
    mHorizontalFOV=horizontalFOV;
    mVerticalFOV=verticalFOV;
    this->init();
  }

  /// Set the camera image resolution. Valid values are > 0.
  RAYTRACER_EXPORTS void setResolution(size_t x, size_t y)
  {
    mXResolution=x;
    mYResolution=y;
    this->init();
  }

protected:

  // Use these accessors in derived classes.
  RAYTRACER_EXPORTS const Vec3d& direction() const { return mDirection; }
  RAYTRACER_EXPORTS const Vec3d& down()      const { return mDown; }
  RAYTRACER_EXPORTS const Vec3d& right()     const { return mRight; }
  RAYTRACER_EXPORTS const Vec3d& topLeft()   const { return mTopLeft; }

private:

  void init();

  Vec3d  mPosition;            //!< World-space position of camera eye.
  Vec3d  mLookAt;              //!< World-space position of camera target.
  Vec3d  mUp;                  //!< World-space up vector.

  size_t mXResolution;        //!< Resolution in x-axis, i.e. width.
  size_t mYResolution;        //!< Resolution in y-axis, i.e. height.

  double  mHorizontalFOV;       //!< Horizontal field of view in degrees.
  double  mVerticalFOV;         //!< Vertical field of view in degrees.

  Vec3d  mDirection;           //!< World-space vector from camera eye to target.
  Vec3d  mDown;                //!< World-space down direction on camera plane.
  Vec3d  mRight;               //!< World-space right direction on camera plane.
  Vec3d  mTopLeft;             //!< World-space position of top
                              //!< left point on camera plane.
};

} //namespace rt

#endif //CAMERA_HPP_INCLUDE_ONCE
