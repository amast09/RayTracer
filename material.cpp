#include <iostream>
#include <string>

#include "pixel.h"
#include "material.h"

std::ostream& operator<<(std::ostream& s, const material_t& rhs)
{
  s << "material " << rhs.name.c_str() << std::endl;
  s << "{" << std::endl;
  if(!rhs.ambient.iszero()) s << "   ambient " << rhs.ambient << std::endl;
  if(!rhs.diffuse.iszero()) s << "   diffuse " << rhs.diffuse << std::endl;
  if(!rhs.specular.iszero()) s << "  specular " << rhs.specular << std::endl;
  if(rhs.alpha != 0.0) s << "  alpha " << rhs.specular << std::endl;
  if(rhs.ior != 0.0) s << "  ior " << rhs.specular << std::endl;
  s << "}" << std::endl << std::endl;

  return s;
}

std::istream& operator>>(std::istream& s, material_t& rhs)
{
	char		c;
	std::string	attrname;

  s >> rhs.name;

  // consume all chars until we get to '{'
  while(s.good() && s.get(c) && (c != '{'));

  // loop until we hit '}'
  while((c = s.peek()) != '}') {

    // read in attribute name
    s >> attrname;

    // read in attribute and consume whitespace at EOL
    if(attrname == "ambient") s >> rhs.ambient >> std::ws;
    if(attrname == "diffuse") s >> rhs.diffuse >> std::ws;
    if(attrname == "specular") s >> rhs.specular >> std::ws;
    if(attrname == "alpha") s >> rhs.alpha >> std::ws;
    if(attrname == "ior") s >> rhs.ior >> std::ws;
  }

  // eat '}' character
  while(s.good() && s.get(c) && (c != '}'));

  return s;
}
