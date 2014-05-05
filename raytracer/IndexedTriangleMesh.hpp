#ifndef _INDEXEDTRIANGLEMESH_HPP_INCLUDE_ONCE
#define _INDEXEDTRIANGLEMESH_HPP_INCLUDE_ONCE

#include "raytracerConfig.hpp"

#include "Renderable.hpp"
#include "Ray.hpp"

namespace rt
{

class IndexedTriangleMesh : public Renderable
{
public:

  /// Implements the intersection computation between ray and any stored triangle.
  RAYTRACER_EXPORTS bool
    closestIntersectionModel(const Ray &ray, double maxLambda, RayIntersection& intersection) const override;

  RAYTRACER_EXPORTS bool anyIntersectionModel(const Ray &ray, double maxLambda) const override;

  RAYTRACER_EXPORTS bool loadFromOBJ(const std::string &filePath);
  RAYTRACER_EXPORTS bool saveToOBJ(const std::string &filePath, bool textureCoordinates=true, bool normals=true) const;

  RAYTRACER_EXPORTS int addVertex(const Vec3d &v, const Vec3d &n, const Vec3d &uvw)
  {
    mVertexPosition.push_back(v);
    mVertexNormal.push_back(n);
    mVertexTextureCoordinate.push_back(uvw);
    return int(mVertexPosition.size()-1);
  }
  RAYTRACER_EXPORTS void addTriangle(const int i0,const int i1,const int i2)
  {
    mIndices.push_back(i0);
    mIndices.push_back(i1);
    mIndices.push_back(i2);
  }

  RAYTRACER_EXPORTS const std::vector<Vec3d>& vertexPositions()          const {return mVertexPosition;}
  RAYTRACER_EXPORTS const std::vector<Vec3d>& vertexTextureCoordinates() const {return mVertexTextureCoordinate;}
  RAYTRACER_EXPORTS const std::vector<Vec3d>& vertexNormals()            const {return mVertexNormal;}
  RAYTRACER_EXPORTS const std::vector<int>&  triangleIndices()          const {return mIndices;}

protected:

  // Override this method to recompute the bounding box of this object.
  RAYTRACER_EXPORTS BoundingBox computeBoundingBox() const override;

private:
  std::vector<Vec3d>                   mVertexPosition;
  std::vector<Vec3d>                   mVertexTextureCoordinate;
  std::vector<Vec3d>                   mVertexNormal;
  std::vector<int>                    mIndices;
};
} //namespace rt


#endif //_INDEXEDTRIANGLEMESH_HPP_INCLUDE_ONCE