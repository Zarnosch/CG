#include "BezierPatchMesh.hpp"
#include <functional>

namespace rt
{

  BezierPatchMesh::BezierPatchMesh(size_t m,    size_t n,
    size_t resu, size_t resv) : BVHIndexedTriangleMesh(),
  mM(m), mN(n), mResU(resu), mResV(resv)
  {
    // Allocate memory for Bezier points
    mControlPoints.resize(m*n);
  }

  void BezierPatchMesh::initialize()
  {
    //this function samples the underlying continuous patch and tessellates it
    //regularly with triangles

    //sample at triangle vertices at uniform uv parameters
    std::vector<BezierPatchSample> samples; samples.reserve(mResU*mResV);

    for(size_t j=0; j<mResV; ++j)
      for(size_t i=0; i<mResU; ++i)
        samples.push_back(this->sample(double(i) / (mResU-1), double(j) / (mResV-1)));

    //construct triangle strips based on the computed samples
    for(size_t j=0; j<mResV-1; ++j)
    {
      for(size_t i=0; i<mResU-1; ++i)
      {
        const size_t i00 = mResU* j    +  i,
          i10 = mResU* j    + (i+1),
          i01 = mResU*(j+1) +  i,
          i11 = mResU*(j+1) + (i+1);

        {
          //construct lower triangle
          const Vec3d &v0 = samples[i00].position,
            &v1 = samples[i10].position,
            &v2 = samples[i01].position;
          const Vec2d &p0 = samples[i00].uv,
            &p1 = samples[i10].uv,
            &p2 = samples[i01].uv;
          Vec3d  n0 = samples[i00].normal,
            n1 = samples[i10].normal,
            n2 = samples[i01].normal;

          if(!n0.lengthSquared() || !n1.lengthSquared() || !n2.lengthSquared())
          {
            //fall back to triangle normals if normals are not defined
            const Vec3d normal = cross((v1-v0),(v2-v0)).normalize();
            n0 = normal; n1 = normal; n2 = normal;
          }

          int i0=this->addVertex(v0,n0,Vec3d(p0,0));
          int i1=this->addVertex(v1,n1,Vec3d(p1,0));
          int i2=this->addVertex(v2,n2,Vec3d(p2,0));
          this->addTriangle(i0,i1,i2);
        }

        {
          //construct upper triangle
          const Vec3d &v0 = samples[i10].position,
            &v1 = samples[i11].position,
            &v2 = samples[i01].position;
          const Vec2d &p0 = samples[i10].uv,
            &p1 = samples[i11].uv,
            &p2 = samples[i01].uv;
          Vec3d n0 = samples[i10].normal,
            n1 = samples[i11].normal,
            n2 = samples[i01].normal;

          if(!n0.lengthSquared() || !n1.lengthSquared() || !n2.lengthSquared())
          {
            //fall back to triangle normals if normals are not defined
            const Vec3d normal = cross((v1-v0),(v2-v0)).normalize();
            n0 = normal; n1 = normal; n2 = normal;
          }
          int i0=this->addVertex(v0,n0,Vec3d(p0,0));
          int i1=this->addVertex(v1,n1,Vec3d(p1,0));
          int i2=this->addVertex(v2,n2,Vec3d(p2,0));
          this->addTriangle(i0,i1,i2);
        }
      }
    }
    BVHIndexedTriangleMesh::initialize();
  }

  BoundingBox BezierPatchMesh::computeBoundingBox() const
  {
    BoundingBox bbox;
    for (size_t i=0;i<mControlPoints.size();++i)
      bbox.expandByPoint(mControlPoints[i]);
    return bbox;
  }

  float BezierPatchMesh::factorial(float n) const {
	  float f = 1;
	  if (n == 0) return 1;
	  for (unsigned int p = 1; p <= n; p++) 
	  {
		  f *= p;
	  }
	  return f;
  }
  /*
  BezierPatchMesh::BezierPatchSample BezierPatchMesh::sample(double u, double v) const
  {
    BezierPatchSample ret;
    ret.uv = Vec2d(u,v);
	Vec3d x (0,0,0);
    // Programming TASK 2: implement this method
    // You need to compute ret.position and ret.normal!

    // This data will be used within the initialize() function for the 
    // triangle mesh construction.
	//Tensorprodukt
	//std::cout << "::::::::::::::::::::::::::::::::::::::::::::::::::" << std::endl;
	//std::cout << "claculating for u:" << u << " v: " << v << std::endl;
	for (int i = 0; i < mM ; ++i)
	{
		for (int j = 0; j < mN ; ++j)
		{
			double Bernsteinpolynom1;
			double Bernsteinpolynom2;
			double Binomialkoefizient1 = 1;
			double Binomialkoefizient2 = 1;
			if (i <= mM)
			{
				Binomialkoefizient1 = factorial(mM) / (factorial(mM - i) * factorial(i));
			}
			else Binomialkoefizient1 = 0;
			if (j <= mN)
			{
				Binomialkoefizient2 = factorial(mN) / (factorial(mN - j) * factorial(j));
			}
			else Binomialkoefizient2 = 0;
			//if (mM - 1 <= Math::safetyEps() && mM - 1 >= -Math::safetyEps())	Bernsteinpolynom1 = Binomialkoefizient1 * pow(u, i) * 1;
			Bernsteinpolynom1 = Binomialkoefizient1 * pow(u, i) * pow((1 - u), (mM - i));
			//if (mN - 1 <= Math::safetyEps() && mN - 1 >= -Math::safetyEps())	Bernsteinpolynom2 = Binomialkoefizient2 * pow(v, j) * 1;
			Bernsteinpolynom2 = Binomialkoefizient2 * pow(v, j) * pow((1 - v), (mN - j));
			x += (controlPoint(i, j) * Bernsteinpolynom1 * Bernsteinpolynom2);
			//std::cout << "i: " << i << "j: " << j << ", Bp1: " << Bernsteinpolynom1 << ", Bp2: " << Bernsteinpolynom2 << ", cP: " << controlPoint(i, j) << std::endl;
		}
	}
    ret.position = x;
    ret.normal = Vec3d();
    return ret;
  }
  */
  Vec3d BezierPatchMesh::deCasteljau(Vec3d b, double t, size_t r, size_t i) const
  {
	  if (r == 0) {
		  return controlPoint(r, i);
	  }
	  else {
		  return ((1 - t) * (deCasteljau(b, t, (r - 1), i) + (t * deCasteljau(b, t, (r - 1), (i + 1)))));
	  }
  }

  Vec3d BezierPatchMesh::deCasteljauQ(Vec3d b, double t, size_t r, size_t i, std::vector<Vec3d> q) const
  {
	  if (r == 0) {
		  return q.at(i);
	  }
	  else {
		  return ((1 - t) * (deCasteljauQ(b, t, (r - 1), i, q) + (t * deCasteljauQ(b, t, (r - 1), (i + 1), q))));
	  }
  }

  BezierPatchMesh::BezierPatchSample BezierPatchMesh::sample(double u, double v) const
  {
	  BezierPatchSample ret;
	  ret.uv = Vec2d(u, v);

	  //std::vector<Vec3d> P = mControlPoints;

	  size_t m = mM;
	  size_t n = mN;

	  //controlPoints(i,j)
	  std::vector<Vec3d> q;

	  Vec3d p;

	  for (int i = 0; i < m; ++i) {

		  for (size_t r = 1; r <= n; ++r) {
			  for (size_t j = 0; j <= (n - r); ++j) { // 3 - 1 = 2, 3 - 2 = 1, 3 - 3 = 0
				  q[i] = deCasteljau(controlPoint(j, 0), v, i + 1, j);
			  }
		  }


		  for (int j = 0; j < n; ++j) {

		  }

		  deCasteljauQ(q[i], u, , , q)

	  }






	  for (size_t r = 1; r <= n; ++r) {
		  for (size_t i = 0; i <= (n + r); ++i) {
			  p = deCasteljauQ(temp[i], u, r, i, temp);
		  }
	  }

	  // Programming TASK 2: implement this method
	  // You need to compute ret.position and ret.normal!

	  // This data will be used within the initialize() function for the 
	  // triangle mesh construction.

	  ret.position = p;
	  ret.normal = Vec3d();
	  return ret;
  }

} //namespace rt
