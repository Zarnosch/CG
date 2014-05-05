#include "Light.hpp"

namespace rt
{

Light::Light(const Vec3d &position, const Vec3d &spectralIntensity) :
  mPosition(position), mSpectralIntensity(spectralIntensity)
{

}

} //namespace rt
