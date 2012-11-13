#ifndef SPHERE_H
#define SPHERE_H

class sphere_t : public object_t
{
  public:
  // constructors (overloaded)
  sphere_t(std::string token) : \
        object_t(token), \
	radius(0.0) \
	{ };

  // copy constructor
  sphere_t(const sphere_t& rhs) : \
        object_t(rhs), \
	center(rhs.center), \
	radius(rhs.radius)
	{ };

  // destructors (default not ok)
  ~sphere_t() \
	{ };

  // operators (incl. assignment operator)
  const sphere_t& operator=(const sphere_t& rhs)
	{
	  if(this != &rhs) {
	    cookie = rhs.cookie;
	    type = rhs.type;
	    name = rhs.name;
	    material = rhs.material;
	    center = rhs.center;
	    radius = rhs.radius;
	  }
          return *this;
	}

  // friends
  std::ostream& put(std::ostream& s) const;
  std::istream& get(std::istream& s);

  // methods
  double	hits(const vec_t&, const vec_t&, vec_t&, vec_t&);
  vec_t		getcenter()	{ return center; }

  private:
  vec_t	 center;
  double radius;
};

#endif
