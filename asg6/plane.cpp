#include <iostream>
#include <string>
#include <cassert>
#include <cmath>

#include "vector.h"
#include "pixel.h"
#include "material.h"
#include "object.h"
#include "plane.h"

std::ostream& plane_t::put(std::ostream& s) const
{
//std::cerr << "in plane_t::operator<<\n";

  // get parent to print itself
  object_t::put(s);

  s << "  normal " << normal << std::endl;
  s << "  point  "<< point << std::endl;

  s << "}" << std::endl << std::endl;

  return s;
}

std::istream& plane_t::get(std::istream& s)
{
	char		c;
	std::string	attrname;

//std::cerr << "in plane_t::operator>>\n";

  // get parent to read itself
  object_t::get(s);

  // loop until we hit '}'
  while((c = s.peek()) != '}') {

    // read in attribute name and trailing whitespace
    s >> attrname >> std::ws;

    if(attrname == "normal")
      s >> normal >> std::ws;

    if(attrname == "point")
      s >> point >> std::ws;
  }

  // eat '}' character
  while(s.good() && s.get(c) && (c != '}'));

  // init
  ndotq = point.dot(normal.norm());

  return s;
}

double	plane_t::hits(const vec_t& pos, const vec_t& dir, vec_t& hit, vec_t& N)
{
	double	ndotd,t,ndotb,precision=0.000001;
	vec_t	last_hit;

  assert(cookie == OBJ_COOKIE);

  // get angle between view direction and plane normal
  ndotd = dir.dot(normal.norm());

  // if ndotd == 0, ray is parallel to the plane, no intersection
  if(fabs(ndotd) < precision) return(-1);

  // find distance of ray pos to plane: N . b + ndotq
  ndotb = pos.dot(normal.norm());

  // here we need to temporarily negate N to find distance along ray to plane
  t = (ndotq - ndotb) / ndotd;

  // if t < 0 then intersection "behind" where ray started from
  if(t <= 0) return(-1);

  // hit point is obtained by moving along dir by t amount & adding to ray pos
  last_hit = pos + t * dir;

  // if z > 0 then hit point is in front of image plane; invalid intersection
  if(last_hit[2] > 0.0) return(-1);
  else                  { hit = last_hit; N = normal.norm(); }

  return(t);
}
