#ifndef LIGHT_H
#define LIGHT_H

#define LGT_COOKIE 30492344

class light_t
{
  public:
  // constructors (overloaded)
  light_t() : \
	cookie(LGT_COOKIE), \
	location(0.0,0.0,0.0), \
	color(0.0,0.0,0.0)
	{ };

  // copy constructor
  light_t(const light_t& rhs) : \
	cookie(rhs.cookie), \
	name(rhs.name), \
	location(rhs.location), \
	color(rhs.color) \
	{ };

  // destructors (default ok, no 'new' in constructor)
  ~light_t()
	{ };

  // operators (incl. assignment operator)
  const light_t& operator=(const light_t& rhs)
	{
	  if(this != &rhs) {
	    cookie = rhs.cookie;
	    name = rhs.name;
	    location = rhs.location;
	    color = rhs.color;
	  }
          return *this;
	}

  // friends
  friend std::ostream& operator<<(std::ostream& s, const light_t& rhs);
  friend std::ostream& operator<<(std::ostream& s, light_t *rhs)
		{ return(s << (*rhs)); }

  friend std::istream& operator>>(std::istream& s, light_t& rhs);
  friend std::istream& operator>>(std::istream& s, light_t *rhs)
		{ return(s >> (*rhs)); }

  // methods
  int		getcookie()	{ return cookie; }
  std::string	getname()	{ return name; }
  vec_t         getlocation()	{ return location; }
  rgb_t<double>	getcolor()	{ return color; }

  private:
  int		cookie;		// magic number
  std::string	name;		// name
  vec_t		location;	// location in world coords
  rgb_t<double>	color;		// light color
};

#endif
