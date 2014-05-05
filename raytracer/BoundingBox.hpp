#ifndef BOUNDINGBOX_HPP_INCLUDE_ONCE
#define BOUNDINGBOX_HPP_INCLUDE_ONCE

#include "raytracerConfig.hpp"

#include "Math.hpp"
#include "Ray.hpp"

namespace rt
{

/**
 * Represents an axis-aligned bounding box.
 */
class BoundingBox
{
public:
  RAYTRACER_EXPORTS BoundingBox() : 
    mMin( std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()),
    mMax(-std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity())
  { } 
  RAYTRACER_EXPORTS BoundingBox(const Vec3d &min, const Vec3d &max) : mMin(min), mMax(max) { }
        
  // Returns true in case of any intersection
  RAYTRACER_EXPORTS bool anyIntersection(const Ray &ray, double maxLambda) const;

  RAYTRACER_EXPORTS void merge(const BoundingBox& box); // merge with a given bounding box
        
  RAYTRACER_EXPORTS const Vec3d& min() const { return mMin; }
  RAYTRACER_EXPORTS const Vec3d& max() const { return mMax; }

  RAYTRACER_EXPORTS void setMin(const Vec3d &min) { mMin = min; }
  RAYTRACER_EXPORTS void setMax(const Vec3d &max) { mMax = max; }

  RAYTRACER_EXPORTS void expandByPoint(const Vec3d& point);

  RAYTRACER_EXPORTS double computeArea() const;

private:
    Vec3d mMin; // Coordinate of min corner
    Vec3d mMax; // Coordinate of max corner
};

}

#endif

