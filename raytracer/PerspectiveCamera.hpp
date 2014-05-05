#ifndef PERSPECTIVECAMERA_HPP_INCLUDE_ONCE
#define PERSPECTIVECAMERA_HPP_INCLUDE_ONCE

#include "raytracerConfig.hpp"

#include "Camera.hpp"

namespace rt {

class PerspectiveCamera : public Camera
{
public:
  RAYTRACER_EXPORTS Ray ray(size_t x, size_t y) const override;
};

} //namespace rt

#endif //PERSPECTIVECAMERA_HPP_INCLUDE_ONCE
