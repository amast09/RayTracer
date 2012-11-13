#include <iostream>
#include <string>
#include <cassert>
#include <cmath>

#include "vector.h"
#include "pixel.h"
#include "material.h"
#include "object.h"
#include "sphere.h"

std::ostream& sphere_t::put(std::ostream& s) const
{
//std::cerr << "in sphere_t::operator<<\n";

  // get parent to print itself
  object_t::put(s);

  s << "  center " << center << std::endl;
  s << "  radius  "<< radius << std::endl;

  s << "}" << std::endl << std::endl;

  return s;
}

std::istream& sphere_t::get(std::istream& s)
{
	char		c;
	std::string	attrname;

//std::cerr << "in sphere_t::operator>>\n";

  // get parent to read itself
  object_t::get(s);

  // loop until we hit '}'
  while((c = s.peek()) != '}') {

    // read in attribute name and trailing whitespace
    s >> attrname >> std::ws;

    if(attrname == "center")
      s >> center >> std::ws;

    if(attrname == "radius")
      s >> radius >> std::ws;
  }

  // eat '}' character
  while(s.good() && s.get(c) && (c != '}'));

  return s;
}

double	sphere_t::hits(const vec_t& pos, const vec_t& dir, vec_t& hit, vec_t& N)
{
	vec_t	pc;
	double	a=1.0,b,c;
	double	d, t0, t1, tf;

  assert(cookie == OBJ_COOKIE);

  // get vector from sphere center to ray base (pos - center)
  pc = pos - center;

  // compute coeffs for quadratic formula, a should be 1.0 if dir normalized
  a = dir.dot(dir);
  b = 2.0 * pc.dot(dir);
  c = pc.dot(pc) - radius * radius;

  // determine the discriminant from the quadratic formula
  d = b*b - 4.0*a*c;

  // if discriminant positive, solve for t
  if(d > 0) {
    // t is the distance from ray's base to hit on sphere, always take the
    // smaller of the two roots as we want the entry wound and not the exit
    t0 = (-b - sqrt(d))/(2.0*a);
    t1 = (-b + sqrt(d))/(2.0*a);

    if(t0 > .00001 && t0 < t1)
    {
	tf = t0;
    }
	else
	{
		tf = t1;
	}

    hit = pos + tf * dir;

    // the final step is to compute the normal at the surface---for a sphere
    // the normal is simply a vector pointing from the center to the hit point
    N = 1.0/radius * (hit - center);
  }
     return(tf);
}


