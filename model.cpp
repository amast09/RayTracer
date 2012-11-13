#include <omp.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cassert>
#include <cstring>
#include <cmath>

#include "vector.h"
#include "pixel.h"
#include "camera.h"
#include "light.h"
#include "material.h"
#include "object.h"
#include "list.h"
#include "plane.h"
#include "sphere.h"
#include "model.h"
#include "ray.h"
#include "photon.h"
#include "timer.h"
#include "kdtree.h"

#ifndef INFINITY
#define INFINITY HUGE
#endif

std::istream& operator>>(std::istream& s, model_t& rhs)
{
	light_t		*lgt;
	material_t	*mat;
	object_t	*obj;
	std::string	token,name;

  while(!s.eof()) {

    if((s >> token).good()) {

      if(token == "light") {
        s >> (lgt = new light_t());
        std::cerr << "loaded " << lgt->getname() << std::endl;
        rhs.lgts.push_back(lgt);
      }

      if(token == "camera") {
        s >> rhs.cam;
        std::cerr << "loaded " << rhs.cam.getname() << std::endl;
      }

      if(token == "material") {
        s >> (mat = new material_t());
        std::cerr << "loaded " << mat->getname() << std::endl;
        rhs.mats.push_back(mat);
      }

      if(token == "plane") {
        s >> (obj = new plane_t(token));
        std::cerr << "loaded " << obj->getname() << std::endl;
        rhs.objs.push_back(obj);
      }

      if(token == "sphere") {
        s >> (obj = new sphere_t(token));
        std::cerr << "loaded " << obj->getname() << std::endl;
        rhs.objs.push_back(obj);
      }
    }
  }

  return s;
}

std::ostream& operator<<(std::ostream& s, model_t& rhs)
{
	light_t				*lgt;
	object_t			*obj;
	material_t			*mat;

	list_t<light_t* >::iterator	litr;
	list_t<material_t* >::iterator	mitr;
	list_t<object_t* >::iterator	oitr;

  // print out camera 
  s << rhs.cam;

  // print out lights, materials, objects 
  for(litr = rhs.lgts.begin(); litr != rhs.lgts.end(); litr++) s << *litr;
  for(mitr = rhs.mats.begin(); mitr != rhs.mats.end(); mitr++) s << *mitr;
  for(oitr = rhs.objs.begin(); oitr != rhs.objs.end(); oitr++) s << *oitr;

  return s;
}

object_t* model_t::find_closest(vec_t& pos, vec_t& dir,
                                double& dist, vec_t& hit, vec_t& N)
{
	// candidate object
	double				c_dist;
	vec_t				c_hit,c_N;
	object_t			*c_obj=NULL;

	// closest object thus far
	double				closest_dist=INFINITY;
	vec_t				closest_hit,closest_N;
	object_t			*closest_obj=NULL;

	list_t<object_t* >::iterator	oitr;

  // check each object; basic list iteration
  for(oitr = objs.begin(); oitr != objs.end(); oitr++) {
    c_obj = (object_t *)*oitr;
    // if intersection behind us, ignore
    if((c_dist = c_obj->hits(pos,dir,c_hit,c_N))<0)
      continue;
    // if distance to hit is smaller, set closest
    else if((0.00001 < c_dist) && (c_dist < closest_dist)) {
      closest_dist = c_dist;
      closest_obj = c_obj;
      closest_hit = c_hit;
      closest_N = c_N;
    } else
      continue;
  }

  if(closest_obj) {
    hit = closest_hit;
    N = closest_N;
    dist += closest_dist;
  }

  return(closest_obj);
}

material_t* model_t::getmaterial(std::string name)
{
	material_t			*mat;
	list_t<material_t* >::iterator	mitr;

  for(mitr = mats.begin(); mitr != mats.end(); mitr++) {
    mat = (material_t *)*mitr;
    assert(mat->getcookie() == MAT_COOKIE);
    if(mat->getname() == name) return mat;
  }
  return NULL;
}

void	model_t::shoot(std::vector<photon_t* >& photons)
{
  list_t<light_t* >::iterator	litr;

  for(litr = lgts.begin(); litr != lgts.end(); litr++)
  {
    int i;

    for(i = 0; i < 2000; i++)
    {
      photon_t* phot = new photon_t((*litr)->getlocation());

      if(phot -> caustic(*this, 0))
        { photons.push_back(phot); }
    }

    for(i = 0; i < 2000; i++)
    {
      photon_t* phot = new photon_t((*litr)->getlocation());

      if(phot -> global(*this, 0))
       { photons.push_back(phot); }
    }
  }
}

