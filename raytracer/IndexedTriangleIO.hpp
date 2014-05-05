#ifndef _INDEXEDTRIANGLEIO_INCLUDE_ONCE
#define _INDEXEDTRIANGLEIO_INCLUDE_ONCE

#include "raytracerConfig.hpp"

#include "Math.hpp"
#include <map>
#include <vector>

namespace rt
{

class IndexedTriangleIO
{
public:
  RAYTRACER_EXPORTS void clear();
  RAYTRACER_EXPORTS bool loadFromOBJ(const std::string &filePath);
  RAYTRACER_EXPORTS bool saveToOBJ(const std::string &filePath, bool textureCoordinates=true, bool normals=true) const;

  RAYTRACER_EXPORTS const std::vector<Vec3d>& vertexPositions()          const {return mVertexPosition;}
  RAYTRACER_EXPORTS const std::vector<Vec3d>& vertexTextureCoordinates() const {return mVertexTextureCoordinate;}
  RAYTRACER_EXPORTS const std::vector<Vec3d>& vertexNormals()            const {return mVertexNormal;}
  RAYTRACER_EXPORTS const std::vector<int>&  triangleIndices()          const {return mIndices;}

  RAYTRACER_EXPORTS void setVertexPositions(const std::vector<Vec3d>& v) {mVertexPosition=std::vector<Vec3d>(v);}
  RAYTRACER_EXPORTS void setVertexTextureCoordinates(const std::vector<Vec3d>& v) {mVertexTextureCoordinate=std::vector<Vec3d>(v);}
  RAYTRACER_EXPORTS void setVertexNormals(const std::vector<Vec3d>& v) {mVertexNormal=std::vector<Vec3d>(v);}
  RAYTRACER_EXPORTS void setTriangleIndices(const std::vector<int>& v) {mIndices=std::vector<int>(v);}

private:

  // Container for printing warnings on Wavefront OBJ format specs
  // that are unsupported by this loader. 
  struct UnsupportedFeatures
  {
    UnsupportedFeatures() : hasLine(false),
      hasGroup(false),
      hasSmoothingGroup(false),
      hasMaterial(false) 
    { }

    bool hasLine;
    bool hasGroup;
    bool hasSmoothingGroup;
    bool hasMaterial;
  };
  
  // Vertex flattening or duplication needs strict weak ordering
  struct CompareVec3i
  {
    bool operator()(const Vec3i& lhs, const Vec3i& rhs) const
    {
      for (int i=0;i<3;++i)
      {
        if (lhs[i] < rhs[i]) return true;
        if (lhs[i] > rhs[i]) return false;
      }
      return false;
    }
  };

  void insertFlattenedTriangleV  (const Vec3i &v);
  void insertFlattenedTriangleVT (const Vec3i &v,const Vec3i &t);
  void insertFlattenedTriangleVN (const Vec3i &v,const Vec3i &n);
  void insertFlattenedTriangleVTN(const Vec3i &v,const Vec3i &t,const Vec3i &n);
  int insertFlattenedVertexVTN(const int v, const int t, const  int n);
  void clearCaches();

  std::vector<Vec3d>                   mVertexPosition;
  std::vector<Vec3d>                   mVertexTextureCoordinate;
  std::vector<Vec3d>                   mVertexNormal;
  std::vector<int>                    mIndices;

  std::vector<Vec3d>                   mVertexPositionCache;
  std::vector<Vec3d>                   mVertexTextureCoordinateCache;
  std::vector<Vec3d>                   mVertexNormalCache;
  std::map<Vec3i,size_t,CompareVec3i> mIndicesCache;
};

} //namespace rt


#endif //_INDEXEDTRIANGLEIO_INCLUDE_ONCE