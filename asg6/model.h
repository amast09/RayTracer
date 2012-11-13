#ifndef MODEL_H
#define MODEL_H

#include <vector>

// forward declarations
class ray_t;
class photon_t;

class model_t
{
  public:
  // constructors (overloaded)
  model_t() : \
	lgts(), \
        cam(), \
	mats(), \
	objs() \
	{ }

  // copy constructor
  model_t(const model_t& rhs)
	{
	  lgts = rhs.lgts;
	  cam = rhs.cam;
	  mats = rhs.mats;
	  objs = rhs.objs;
	}

  // destructors (default ok)
  ~model_t()
	{ };

  // operators (incl. assignment operator)
  const model_t& operator=(const model_t& rhs)
	{
          return *this;
	}

  // friends
  friend std::ostream& operator<<(std::ostream& s, model_t& rhs);
  friend std::ostream& operator<<(std::ostream& s, model_t *rhs)
	{ return(s << (*rhs)); }

  friend std::istream& operator>>(std::istream& s, model_t& rhs);
  friend std::istream& operator>>(std::istream& s, model_t *rhs)
	{ return(s >> (*rhs)); }

  friend class ray_t;


  // members
  int		getpixel_w()	{ return cam.getw(); }
  int		getpixel_h()	{ return cam.geth(); }
  double	getworld_w()	{ return cam.getww(); }
  double	getworld_h()	{ return cam.getwh(); }
  vec_t		getviewpoint()	{ return cam.getview_point(); }

  void		shoot(std::vector<photon_t* >& photons);
  object_t*	find_closest(vec_t&,vec_t&,double&,vec_t&,vec_t&);
  material_t*	getmaterial(std::string name);

  // data members
  private:
  list_t<light_t* >		lgts;
  camera_t			cam;
  list_t<material_t* >		mats;
  list_t<object_t* >		objs;
};

#endif
