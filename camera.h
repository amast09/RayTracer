#ifndef CAMERA_H
#define CAMERA_H

#define CAM_COOKIE 49495923

class camera_t
{
  public:
  // constructors (overloaded)
  camera_t() : \
	cookie(CAM_COOKIE), \
	view_point(0.0,0.0,0.0) \
	{
	  for(int i=0; i<2; i++) pixel_dim[i] = 0;
	  for(int i=0; i<2; i++) world_dim[i] = 0.0;
	}

  // copy constructor
  camera_t(const camera_t& rhs) : \
	cookie(rhs.cookie), \
	name(rhs.name), \
	view_point(rhs.view_point) \
	{
	  for(int i=0; i<2; i++) pixel_dim[i] = rhs.pixel_dim[i];
	  for(int i=0; i<2; i++) world_dim[i] = rhs.world_dim[i];
	}

  // destructors (default ok, no 'new' in constructor)
  ~camera_t()
	{ };

  // operators (incl. assignment operator)
  const camera_t& operator=(const camera_t& rhs)
	{
	  if(this != &rhs) {
	    cookie = rhs.cookie;
	    name = rhs.name;
	    view_point = rhs.view_point;
	    for(int i=0; i<2; i++) pixel_dim[i] = rhs.pixel_dim[i];
	    for(int i=0; i<2; i++) world_dim[i] = rhs.world_dim[i];
	  }
          return *this;
	}

  // friends
  friend std::ostream& operator<<(std::ostream& s, const camera_t& rhs);
  friend std::ostream& operator<<(std::ostream& s, camera_t *rhs)
		{ return(s << (*rhs)); }

  friend std::istream& operator>>(std::istream& s, camera_t& rhs);
  friend std::istream& operator>>(std::istream& s, camera_t *rhs)
		{ return(s >> (*rhs)); }

  // methods
  int		getcookie()	{ return cookie; }
  std::string	getname()	{ return name; }
  vec_t		getview_point()	{ return view_point; }
  int           getw()		{ return pixel_dim[0]; }
  int           geth()		{ return pixel_dim[1]; }
  double        getww()		{ return world_dim[0]; }
  double        getwh()		{ return world_dim[1]; }

  private:
  int		cookie;		// magic number
  std::string	name;		// name
  int		pixel_dim[2];	// projection screen size in pixels
  double	world_dim[2];	// screen size in world coords
  vec_t		view_point;	// viewpoint location in world coords
};

#endif
