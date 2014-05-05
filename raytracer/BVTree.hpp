#ifndef BVTREE_HPP_INCLUDE_ONCE
#define BVTREE_HPP_INCLUDE_ONCE

#include "raytracerConfig.hpp"

#include <vector>
#include <stack>
#include "Math.hpp"
#include "BoundingBox.hpp"

namespace rt
{

class BVTree
{
public:
  RAYTRACER_EXPORTS BVTree();


  //build from indexed triangle set
  RAYTRACER_EXPORTS void build(const std::vector<Vec3d> &vertexPositions,const std::vector<Vec3i> &triangleIndices);

  //returns a set of triangle indices as candidates for ray-triangle intersection
  RAYTRACER_EXPORTS const std::vector<int>& intersectBoundingBoxes(const Ray &ray, const double maxLambda) const;
private:

  struct Node
  {
    Node() {left=0;right=0;}
    int left;
    int right;
    BoundingBox bbox;
  };
  void sortTriangles();
  void createNodes(const std::vector<Vec3d> &vertexPositions,
    const std::vector<Vec3i> &triangleIndices);

  void printSortedIndicesStatus(unsigned int offset, unsigned int numTriangles);
  void buildHierarchy(unsigned int rootNodeIndex, unsigned int numTriangles,unsigned int offset);

  void computeBoundingBoxAreas(unsigned int offset, unsigned int numTriangles);

  std::vector<Node> mNodes;

  std::vector<bool>        mTempMarker;
  std::vector<BoundingBox> mTempTriangleBoxes;
  std::vector<Vec3i>       mTempSortedTriangleIndices;
  std::vector<Vec3i>       mTempBufferTriangleIndices;
  std::vector<Vec3d>        mTempAreasLeft;
  std::vector<Vec3d>        mTempAreasRight;
  mutable std::vector<std::vector<int>>         mTempCandidates;   
  mutable std::vector<std::stack<int>>          mTempTraversalJobs;

};
}

#endif //BVTREE_HPP_INCLUDE_ONCE
