#include <iostream>
#include <string>
#include <cassert>

#include "vector.h"
#include "pixel.h"
#include "light.h"

std::ostream& operator<<(std::ostream& s, const light_t& rhs)
{
  assert(rhs.cookie == LGT_COOKIE);

  // print out 'light' token and light name
  s << "light " << rhs.name.c_str() << std::endl;

  s << "{" << std::endl;
  s << "  location " << rhs.location << std::endl;
  s << "  emissivity " << rhs.color << std::endl;
  s << "}" << std::endl << std::endl;

  return s;
}

std::istream& operator>>(std::istream& s, light_t& rhs)
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
    if(attrname == "location") s >> rhs.location >> std::ws;
    if(attrname == "emissivity") s >> rhs.color >> std::ws;
  }

  // eat '}' character
  while(s.good() && s.get(c) && (c != '}'));

  return s;
}
