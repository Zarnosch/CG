#include "Scene.hpp"
#include "Camera.hpp"
#include "Image.hpp"
#include <algorithm>

namespace rt
{

Scene::Scene() : mBackgroundColor(0,0,0,0)
{
}

Scene::~Scene()
{
}

bool
Scene::closestIntersection(const Ray &ray,RayIntersection& intersection, double maxLambda) const
{ 
  double closestLambda = maxLambda;
  RayIntersection tmpIntersection;
  RayIntersection closestIntersection;
  bool hit(false);
  for (size_t i=0;i<mRenderables.size();++i)
  {
    std::shared_ptr<Renderable> r = mRenderables[i];
    if(r->closestIntersection(ray,closestLambda,tmpIntersection))
    {
      if(tmpIntersection.lambda() < closestLambda)
      {
        hit=true;
        closestLambda = tmpIntersection.lambda();
        closestIntersection = tmpIntersection;
      }
    }
  }
  intersection=closestIntersection;
  return hit;
}

bool Scene::anyIntersection(const Ray &ray, double maxLambda) const
{
  for (size_t i=0;i<mRenderables.size();++i)
  {
    std::shared_ptr<Renderable> r = mRenderables[i];
    if(r->anyIntersection(ray,maxLambda))
      return true;
  }
  return false;
}

void Scene::prepareScene()
{
  for(size_t i=0;i<mRenderables.size();++i)
  {
    mRenderables[i]->updateBoundingBox();
    mRenderables[i]->initialize();
    mRenderables[i]->updateTransforms();
  }
}

}
