#include "BVTree.hpp"
#if defined(_OPENMP)
#include <omp.h>
#endif

namespace rt
{


BVTree::BVTree()
{
#if defined(_OPENMP)
  mTempCandidates.resize(omp_get_max_threads());
  mTempTraversalJobs.resize(omp_get_max_threads());
#else
  mTempCandidates.resize(1);
  mTempTraversalJobs.resize(1);
#endif
}

const std::vector<int>& BVTree::intersectBoundingBoxes(const Ray &ray, const double maxLambda) const
{
#if defined(_OPENMP)
  std::vector<int>& tempCandidates=mTempCandidates[omp_get_thread_num()];
  std::stack<int>& tempTraversalJobs = mTempTraversalJobs[omp_get_thread_num()];
#else
  std::vector<int>& tempCandidates=mTempCandidates[0];
  std::stack<int>& tempTraversalJobs = mTempTraversalJobs[0];
#endif

 // int* jobs = (int*)alloca(sizeof(int)*100); //yields 25% better performance

  tempCandidates.clear();
  tempTraversalJobs.push(0);

  while(!tempTraversalJobs.empty())
  {
    //take current job
    int node = tempTraversalJobs.top();
    tempTraversalJobs.pop();

    //test ray vs. bounding box of node
    if(mNodes[node].bbox.anyIntersection(ray,maxLambda))
    {
      if(mNodes[node].left <= 0 && mNodes[node].right == -1) // is a leaf node
        tempCandidates.push_back(-mNodes[node].left);
      else//is not a leaf
      {
        tempTraversalJobs.push(mNodes[node].left);
        tempTraversalJobs.push(mNodes[node].right);
      }
    }
  }
  return tempCandidates;
}

void BVTree::createNodes(const std::vector<Vec3d> &vertexPositions,
    const std::vector<Vec3i> &triangleIndices)
{
  size_t n=triangleIndices.size();
  mTempTriangleBoxes.resize(n);

  for(unsigned int i=0;i<n;++i)
  {
    const Vec3d &v0 = vertexPositions[triangleIndices[i][0]];
    const Vec3d &v1 = vertexPositions[triangleIndices[i][1]];
    const Vec3d &v2 = vertexPositions[triangleIndices[i][2]];

    Vec3d newMin,newMax;
    for(unsigned int j=0;j<3;++j)
    {
      newMin[j]=std::min(v0[j],std::min(v1[j],v2[j]));
      newMax[j]=std::max(v0[j],std::max(v1[j],v2[j]));
    }
    mTempTriangleBoxes[i].setMin(newMin);
    mTempTriangleBoxes[i].setMax(newMax);
  }
}

void BVTree::sortTriangles()
{
  const std::vector<BoundingBox> & triangleNodes = mTempTriangleBoxes;
  size_t n=triangleNodes.size();

  mTempSortedTriangleIndices.resize(n);

  std::vector<int> permutation(n);

  //init permutation
  for(unsigned int i=0;i<n;++i)
    permutation[i]=int(i);

  for(unsigned int dim=0;dim<3;++dim)
  {
    //sort according to dim-coordinate of bounding box centers
    std::sort(permutation.begin(), permutation.end(), [&triangleNodes,dim](int a, int b) -> bool
    {
      return (triangleNodes[a].max()[dim]+triangleNodes[a].min()[dim]) <
             (triangleNodes[b].max()[dim]+triangleNodes[b].min()[dim]);
    });

    //copy permutation indices to result
    for(unsigned int i=0;i<n;++i)
      mTempSortedTriangleIndices[i][dim]=permutation[i];
  }
}

void BVTree::build(const std::vector<Vec3d> &vertexPositions,const std::vector<Vec3i> &triangleIndices)
{
  //create bounding boxes for all triangles
  this->createNodes(vertexPositions,triangleIndices);
  unsigned n = unsigned(triangleIndices.size());
  //sort triangles by x,y and z
  //sortedTriangles now holds sorting permutations for orderings w.r.t x,y,z in the three components
  this->sortTriangles();

  mTempMarker.resize(n);
  mTempBufferTriangleIndices.resize(n);
  mTempAreasLeft.resize(n);
  mTempAreasRight.resize(n);

  //Create root node
  {
    Node root;
    BoundingBox nodeBox;
    for(unsigned int i=0;i<mTempTriangleBoxes.size();++i)
      nodeBox.merge(mTempTriangleBoxes[i]);
    root.bbox=nodeBox;
    mNodes.push_back(root);
  }

  this->buildHierarchy(0,0,n);

  //clear temporary storage
  std::vector<bool>().swap(mTempMarker);
  std::vector<BoundingBox>().swap(mTempTriangleBoxes);
  std::vector<Vec3i>().swap(mTempSortedTriangleIndices);
  std::vector<Vec3i>().swap(mTempBufferTriangleIndices);
  std::vector<Vec3d>().swap(mTempAreasLeft);
  std::vector<Vec3d>().swap(mTempAreasRight);


//   std::cerr<<"created bvh tree with "<<mNodes.size()<<" nodes"<<std::endl;
//
//   //print leaf node indices
//   for(size_t i=0;i<mNodes.size();++i)
//   {
//     const Node &node=mNodes[i];
//     if(node.left<1) //either root or leaf
//     {
//       std::cerr<<"triangle index: "<<abs(node.left)<<std::endl;
//     }
//   }
}
void BVTree::printSortedIndicesStatus(unsigned int offset, unsigned int numTriangles)
{
  std::cerr<<"Indices Status\nX:";
  for(size_t i=offset;i<numTriangles+offset;++i)
    std::cerr<<int(mTempMarker[mTempSortedTriangleIndices[i][0]]);
  std::cerr<<"\nY:";
  for(size_t i=offset;i<numTriangles+offset;++i)
    std::cerr<<int(mTempMarker[mTempSortedTriangleIndices[i][1]]);
  std::cerr<<"\nZ:";
  for(size_t i=offset;i<numTriangles+offset;++i)
    std::cerr<<int(mTempMarker[mTempSortedTriangleIndices[i][2]]);
  std::cerr<<std::endl;

  std::cerr<<"Indices Order\nX:";
  for(size_t i=offset;i<numTriangles+offset;++i)
    std::cerr<<int(mTempSortedTriangleIndices[i][0])<<",";
  std::cerr<<"\nY:";
  for(size_t i=offset;i<numTriangles+offset;++i)
    std::cerr<<int(mTempSortedTriangleIndices[i][1])<<",";
  std::cerr<<"\nZ:";
  for(size_t i=offset;i<numTriangles+offset;++i)
    std::cerr<<int(mTempSortedTriangleIndices[i][2])<<",";
  std::cerr<<std::endl;
}
void BVTree::buildHierarchy(unsigned int rootNodeIndex,unsigned int offset, unsigned int numTriangles)
{
//  std::cerr<<std::endl;
//  std::cerr<<"Num Triangles: "<<numTriangles<<std::endl;

  //1. Find split axis: can be either x,y or z axis
  unsigned int splitDimension=0;
  unsigned int splitIndex=0;
  double minCosts=std::numeric_limits<double>::max();

  //compute total area
  double totalArea=mNodes[rootNodeIndex].bbox.computeArea();

  // The split heuristics requires the computation of costs
  this->computeBoundingBoxAreas(offset,numTriangles);

  //now evaluate costs for each split axis
  for(unsigned int dim=0;dim<3;++dim)
  {
    //evaluate costs
    for(unsigned int i=0;i<numTriangles-1;++i)
    {
      double leftArea = mTempAreasLeft[i+offset][dim];
      double rightArea = mTempAreasRight[i+offset][dim];
      double splitCosts = (leftArea*(i+1)+rightArea*(numTriangles-1-i))/totalArea;
      if(splitCosts < minCosts)
      {
        minCosts=splitCosts;
        splitDimension=dim;
        splitIndex=i;
      }
    }
  }

//  std::cerr<<"Split Costs: "<<minCosts<<std::endl;
//  std::cerr<<"Split Dimension: "<<splitDimension<<std::endl;
//  std::cerr<<"Split Index: "<<splitIndex<<std::endl;

  //now mark the first i triangles as left, the remaining as right
  for(unsigned int i=0;i<=splitIndex;++i)
    mTempMarker[mTempSortedTriangleIndices[i+offset][splitDimension]]=true; //TODO debug here  , maybe splitIndex is not compatible with offset?!
  for(unsigned int i=splitIndex+1;i<numTriangles;++i)
    mTempMarker[mTempSortedTriangleIndices[i+offset][splitDimension]]=false;
//  this->printSortedIndicesStatus(offset,numTriangles);

  //the indices in the remaining two dimensions must be put into left, right order
  // Copy to temporary buffer
  unsigned int leftCount,rightCount;
  for(unsigned int dim=0;dim<3;++dim)
  {
    if(dim==splitDimension)
      continue;
    leftCount=0;
    rightCount=0;
    for(unsigned int i=0;i<numTriangles;++i)
    {
      if(mTempMarker[mTempSortedTriangleIndices[i+offset][dim]])
        mTempBufferTriangleIndices[offset+leftCount++][dim]=mTempSortedTriangleIndices[i+offset][dim];
      else
        mTempBufferTriangleIndices[offset+numTriangles-(++rightCount)][dim]=mTempSortedTriangleIndices[i+offset][dim];
    }
  }

  // Copy from temp buffer into workspace array
  for(unsigned int dim=0;dim<3;++dim)
  {
    if(dim==splitDimension)
      continue;
    for(unsigned int i=0;i<leftCount;++i)
      mTempSortedTriangleIndices[i+offset][dim]=mTempBufferTriangleIndices[offset+i][dim];
    for(unsigned int i=0;i<rightCount;++i)
      mTempSortedTriangleIndices[offset+leftCount+i][dim]=mTempBufferTriangleIndices[offset+numTriangles-1-i][dim];
  }

//  this->printSortedIndicesStatus(offset,numTriangles);

  //create nodes and prepare for recursive call
  Node nodeLeft;
  BoundingBox leftBox;
  for(unsigned int i=0;i<=splitIndex;++i)
    leftBox.merge(mTempTriangleBoxes[mTempSortedTriangleIndices[i+offset][splitDimension]]);
  nodeLeft.bbox = leftBox;
  mNodes.push_back(nodeLeft);
  unsigned idxLeft=unsigned(mNodes.size())-1;
  mNodes[rootNodeIndex].left=int(idxLeft);


  Node nodeRight;
  BoundingBox rightBox;
  for(unsigned int i=splitIndex+1;i<numTriangles;++i)
    rightBox.merge(mTempTriangleBoxes[mTempSortedTriangleIndices[i+offset][splitDimension]]);
  nodeRight.bbox = rightBox;
  mNodes.push_back(nodeRight);
  unsigned idxRight=unsigned(mNodes.size())-1;
  mNodes[rootNodeIndex].right=int(idxRight);

  if(splitIndex==0) //left count == 1  -> left child is a leaf
  {
    mNodes[idxLeft].left=-mTempSortedTriangleIndices[0+offset][splitDimension];
    mNodes[idxLeft].right=-1;
  }
  else
    this->buildHierarchy(idxLeft,offset,splitIndex+1);

  if(splitIndex==numTriangles-2) //right count == 1 -> right child is a leaf
  {
    mNodes[idxRight].left=-mTempSortedTriangleIndices[numTriangles+offset-1][splitDimension];
    mNodes[idxRight].right=-1;
  }
  else
  {
    unsigned int offsetNew=offset+splitIndex+1;
    this->buildHierarchy(idxRight,offsetNew,numTriangles-1-splitIndex);
  }
}

void BVTree::computeBoundingBoxAreas(unsigned int offset, unsigned int numTriangles)
{
  //Left areas
  for(unsigned int dim=0;dim<3;++dim)
  {
    BoundingBox boxLeft;
    for(unsigned int i=offset;i<numTriangles+offset;++i)
    {
      boxLeft.merge(mTempTriangleBoxes[mTempSortedTriangleIndices[i][dim]]); //expand from left to right
      mTempAreasLeft[i][dim]=boxLeft.computeArea();
    }
  }

  //Right areas
  for(unsigned int dim=0;dim<3;++dim)
  {
    BoundingBox boxRight;
    for(unsigned int i=numTriangles+offset-1;i>offset;--i)
    {
      boxRight.merge(mTempTriangleBoxes[mTempSortedTriangleIndices[i][dim]]); //expand from left to right
      mTempAreasRight[i-1][dim]=boxRight.computeArea();
    }
  }
}
} //namespace rt
