#ifndef SCENE_HPP_INCLUDE_ONCE
#define SCENE_HPP_INCLUDE_ONCE

#include "raytracerConfig.hpp"

#include <memory>
#include <vector>

#include "Math.hpp"
#include "Renderable.hpp"

namespace rt
{
  class Image;
  class Camera;
  class Light;
  class Renderable;
  class Ray;
  class RayIntersection;

class Scene
{
public:
  RAYTRACER_EXPORTS Scene();
  RAYTRACER_EXPORTS virtual ~Scene();

  /// Add geometry to the scene.
  RAYTRACER_EXPORTS void addRenderable(std::shared_ptr<Renderable> renderable) {
    mRenderables.push_back(renderable);
  }

  /// Add a point light to the scene.
  RAYTRACER_EXPORTS void addLight(std::shared_ptr<Light> light) {
    mLights.push_back(light);
  }

  /// Returns a vector containing all lights in the scene.
  RAYTRACER_EXPORTS const std::vector<std::shared_ptr<Light>>& lights() const { return mLights; }

  /// Computes the closest intersection of a ray and any object in scene.
  RAYTRACER_EXPORTS bool
  closestIntersection(const Ray &ray, RayIntersection& intersection,
                      double maxLambda = std::numeric_limits<double>::infinity()) const; 

  /// Checks whether a ray intersects any object in the scene.
  RAYTRACER_EXPORTS bool anyIntersection(const Ray &ray,
                       double maxLambda = std::numeric_limits<double>::infinity()) const;

  RAYTRACER_EXPORTS const Vec4d& backgroundColor() const { return mBackgroundColor; }
  std::shared_ptr<Camera> camera()    {return mCamera;}

  RAYTRACER_EXPORTS void setBackgroundColor(const Vec4d& rgba)      { mBackgroundColor = rgba; }
  RAYTRACER_EXPORTS void setCamera(std::shared_ptr<Camera> camera) {mCamera=camera; }

  //prepare scene for rendering
  RAYTRACER_EXPORTS void prepareScene();

private:
  Vec4d mBackgroundColor;

  std::shared_ptr<Camera> mCamera;

  std::vector<std::shared_ptr<Light>> mLights;
  std::vector<std::shared_ptr<Renderable>> mRenderables;
};

} //namespace rt

#endif //SCENE_HPP_INCLUDE_ONCE
