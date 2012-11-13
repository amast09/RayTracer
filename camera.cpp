#include <iostream>
#include <string>
#include <cassert>

#include "vector.h"
#include "pixel.h"
#include "camera.h"

std::ostream& operator<<(std::ostream& s, const camera_t& rhs)
{
  assert(rhs.cookie == CAM_COOKIE);

  // print out 'camera' token and camera name
  s << "camera " << rhs.name.c_str() << std::endl;

  s << "{" << std::endl;
  s << "   pixeldim " << rhs.pixel_dim[0] << " " << rhs.pixel_dim[1] << std::endl;
  s << "   worldim " << rhs.world_dim[0] << " " << rhs.world_dim[1] << std::endl;
  s << "  viewpoint " << rhs.view_point << std::endl;
  s << "}" << std::endl << std::endl;

  return s;
}

std::istream& operator>>(std::istream& s, camera_t& rhs)
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
    if(attrname == "pixeldim") s >> rhs.pixel_dim[0] >> rhs.pixel_dim[1]  >> std::ws;
    if(attrname == "worlddim") s >> rhs.world_dim[0] >> rhs.world_dim[1]  >> std::ws;
    if(attrname == "viewpoint") s >> rhs.view_point >> std::ws;
  }

  // eat '}' character
  while(s.good() && s.get(c) && (c != '}'));

  return s;
}
