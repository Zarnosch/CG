#ifndef RAYTRACER_HPP_INCLUDE_ONCE
#define RAYTRACER_HPP_INCLUDE_ONCE

#include "raytracerConfig.hpp"

#include <vector>
#include <string>
#include <memory>

#include "Math.hpp"

namespace rt
{

class Scene;
class Ray;
class RayIntersection;
class Image;

/// Performs recursive raytracing.
class Raytracer
{
public:
  RAYTRACER_EXPORTS Raytracer(size_t maxDepth=10);
  RAYTRACER_EXPORTS virtual ~Raytracer();

  /// The scene contains all Renderables, Lights, and a Camera.
  RAYTRACER_EXPORTS void setScene(std::shared_ptr<Scene> scene) { mScene=scene; }

  /// Writes RGBA values to an image.
  RAYTRACER_EXPORTS void renderToImage(std::shared_ptr<Image> image) const;

protected:

  /// Returns the color of a traced ray.
  RAYTRACER_EXPORTS Vec4d trace(const Ray &ray,
             size_t depth) const;

  /// Determines the color of an intersection point.
  RAYTRACER_EXPORTS Vec4d shade(const RayIntersection& intersection,
             size_t depth) const;

private:
  size_t mMaxDepth;              ///< Maximum number of ray indirections.
  std::shared_ptr<Scene> mScene;
};

} //namespace rt

#endif //RAYTRACER_HPP_INCLUDE_ONCE
