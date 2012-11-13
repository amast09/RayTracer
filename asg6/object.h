#ifndef OBJECT_H
#define OBJECT_H

#define OBJ_COOKIE 12345678

class object_t
{
  public:
  // constructors (overloaded)
  object_t(std::string itype = std::string("unknown")) : \
	type(itype), \
	cookie(OBJ_COOKIE) \
	{ };

  // copy constructor
  object_t(const object_t& rhs) : \
	cookie(rhs.cookie), \
	type(rhs.type), \
	name(rhs.name), \
	material(rhs.material) \
	{ };

  // destructors (default ok)
  ~object_t()
	{ };

  // operators (incl. assignment operator)
  const object_t& operator=(const object_t& rhs)
	{
	  if(this != &rhs) {
	    cookie = rhs.cookie;
	    type = rhs.type;
	    name = rhs.name;
	    material = rhs.material;
	  }
          return *this;
	}

  // friends
  friend std::ostream& operator<<(std::ostream& s, const object_t& rhs)
	{ return rhs.put(s); }
  friend std::ostream& operator<<(std::ostream& s, object_t *rhs)
	{ return(s << (*rhs)); }

  friend std::istream& operator>>(std::istream& s, object_t& rhs)
	{ return rhs.get(s); }
  friend std::istream& operator>>(std::istream& s, object_t *rhs)
	{ return(s >> (*rhs)); }

  // methods
  int		getcookie()	{ return cookie; }
  std::string	gettype()	{ return type; }
  std::string	getname()	{ return name; }
  std::string	getmaterial()	{ return material; }

  virtual std::ostream& put(std::ostream& s) const;
  virtual std::istream& get(std::istream& s);

  // methods (virtual void, dreived class must define)
  virtual double	hits(const vec_t&,const vec_t&,vec_t&,vec_t&);

  protected:
  int		cookie;		// magic number
  std::string	type;		// e.g., plane, sphere, etc.
  std::string	name;		// e.g., left_wall, center_sphere, etc.
  std::string	material;	// material
};

#endif
