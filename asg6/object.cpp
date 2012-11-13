#include <iostream>
#include <string>

#include "vector.h"
#include "pixel.h"
#include "material.h"
#include "object.h"

std::ostream& object_t::put(std::ostream& s) const
{
//std::cerr << "in object_t::operator<<\n";

  s << type << " " << name << std::endl;
  s << "{" << std::endl;
  s << "  material " << material << std::endl;

  return s;
}

std::istream& object_t::get(std::istream& s)
{
	char		c;
	std::string	attrname;

//std::cerr << "in object_t::operator>>\n";

  s >> name;

  // consume all chars until we get to '{'
  while(s.good() && s.get(c) && (c != '{'));

  // read in material and trailing whitespace (first attribute MUST be material)
  s >> attrname >> attrname >> std::ws;

  // read in attribute and consume whitespace at EOL
  material = attrname;

  return s;
}

double	object_t::hits(const vec_t& pos,const vec_t& dir,vec_t& hit,vec_t& N)
{
  std::cerr << "object_t::hits: shouldn't be called" << std::endl;
}
