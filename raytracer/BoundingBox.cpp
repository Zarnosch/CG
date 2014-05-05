#include "BoundingBox.hpp"

namespace rt
{

void BoundingBox::merge(const BoundingBox& bb)
{
  mMin = Vec3d(std::min(bb.min()[0],mMin[0]),std::min(bb.min()[1],mMin[1]),std::min(bb.min()[2],mMin[2]));
  mMax = Vec3d(std::max(bb.max()[0],mMax[0]),std::max(bb.max()[1],mMax[1]),std::max(bb.max()[2],mMax[2]));
}

bool BoundingBox::anyIntersection(const Ray &ray, double maxLambda) const
{
  double tnear = -std::numeric_limits<double>::max(); 
  double tfar = std::numeric_limits<double>::max();
//   double tnear = 0; 
//   double tfar = maxLambda;
  double t1;
  double t2;

  Vec3d direction = ray.direction();
  Vec3d origin = ray.origin();

  for (int i=0; i<3; i++)
  {
    if (direction[i]) // this is numerically unstable for small values!
    {
      t1 = (mMin[i]-origin[i])/direction[i];
      t2 = (mMax[i]-origin[i])/direction[i];
      if (t1 > t2){
        double tmp = t2;
        t2 = t1;
        t1 = tmp;
      }
      if (t1 > tnear)
        tnear = t1;
      if (t2 < tfar)
        tfar = t2;
      if (tnear > tfar)
        return false;
      if (tfar < 0)
        return false;
    }else{
      if (origin[i]<mMin[i] || origin[i]>mMax[i]) //is this right? we might not reject some trivial cases
        return false;
    }
  }

  return true;
}

void BoundingBox::expandByPoint(const Vec3d& p)
{
  for (int i=0;i<3;++i)
  {
    mMin[i]=std::min(mMin[i],p[i]);
    mMax[i]=std::max(mMax[i],p[i]);
  }
}

double BoundingBox::computeArea() const
{
  Vec3d d = mMax-mMin;

  return 2*(d[0]*d[1]+d[0]*d[2]+d[1]*d[2]);
}

}