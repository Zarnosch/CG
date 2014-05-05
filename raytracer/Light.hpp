#ifndef LIGHT_HPP_INCLUDE_ONCE
#define LIGHT_HPP_INCLUDE_ONCE

#include "raytracerConfig.hpp"

#include "Math.hpp"

namespace rt
{

/// Point light with position and spectral intensity.
class Light
{
public:
  RAYTRACER_EXPORTS Light(const Vec3d &position=Vec3d(0,0,0),
        const Vec3d &spectralIntensity=Vec3d(1,1,1));

  RAYTRACER_EXPORTS const Vec3d& position()          const { return mPosition; }
  RAYTRACER_EXPORTS const Vec3d& spectralIntensity() const { return mSpectralIntensity; }

  RAYTRACER_EXPORTS void setPosition(const Vec3d &position) { mPosition=position; }
  RAYTRACER_EXPORTS void setSpectralIntensity(const Vec3d &spectralIntensity)
  {
    mSpectralIntensity=spectralIntensity;
  }

private:
  Vec3d mPosition;
  Vec3d mSpectralIntensity;
};

} //namespace rt

#endif //LIGHT_HPP_INCLUDE_ONCE
