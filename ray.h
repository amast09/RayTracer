#ifndef RAY_H
#define RAY_H

#include <cstdlib>
//#include "photon.h"
#include "kdtree.h"

// forward declarations
template <typename T, typename P, typename C> 
class kdtree_t;
class photon_t;
class photon_c;


#define MAX_DIST 100
class model_t;
class ray_t
{
  public:
  // constructors (overloaded)
  ray_t() : \
	dis(0.0), \
	pos(0.0,0.0,0.0), \
	dir(0.0,0.0,0.0) \
	{ };

  ray_t(const vec_t& o, const vec_t& d, double r=0.0) : \
	dis(r), \
	pos(o), \
	dir(d) \
	{ };

  // copy constructor
  ray_t(const ray_t& rhs) : \
	dis(rhs.dis), \
	pos(rhs.pos), \
	dir(rhs.dir) \
	{ };

  // operators (incl. assignment operator)
  const ray_t& operator=(const ray_t& rhs)
	{
	  if(this != &rhs) {
	    dis = rhs.dis;
	    pos = rhs.pos;
	    dir = rhs.dir;
	  }
          return *this;
	}

  // methods
  void trace(model_t&,rgb_t<double>&, int bounce, kdtree_t<photon_t, photon_t*, photon_c> kdtree);
  void trace(model_t&,rgb_t<double>&, int bounce);

  // destructors (default ok, no 'new' in constructor)
  ~ray_t()
	{ };

  protected:
  double   dis;	// distance
  vec_t    pos;	// position
  vec_t    dir;	// direction
};
#endif
