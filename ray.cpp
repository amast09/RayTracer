#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cassert>
#include <cmath>

#include "vector.h"
#include "pixel.h"
#include "camera.h"
#include "light.h"
#include "material.h"
#include "object.h"
#include "list.h"
#include "plane.h"
#include "model.h"
#include "photon.h"
#include "ray.h"
#include "kdtree.h"

void ray_t::trace(model_t&		model,
                  rgb_t<double>&	color,
		  int bounce)
{
	// (static) object properties
	object_t			*obj=NULL;
	rgb_t<double>			ambient, diffuse, specular;
	material_t			*mat;
	vec_t				hit,N;	   // hit point and normal
	float				alpha, ior;

	// light and illumination model variables
	light_t				*lgt=NULL;
	vec_t				L,V,R,H;
	rgb_t<double>			I_d, I_s;
	double				r,ndotl=0.0,n=32.0;
	list_t<light_t* >::iterator	litr;

// prevent infinite loops
  if(bounce > 5) return;

  // get closest object, if any
  if(!(obj = model.find_closest(pos,dir,dis,hit,N)) || dis > MAX_DIST)
    return;

  // if hit distance valid, compute color at surface
  if(dis > 0) {

    // get object material properties
    if((mat = model.getmaterial(obj->getmaterial())) != NULL) {
      ambient = mat->getamb();
      diffuse = mat->getdiff();
      specular = mat->getspec();
      alpha = mat->getalpha();
      ior = mat->getior();
    }

    // ambient color
    color += 1.0/dis * ambient;			// ambient scaled by ray dist

    // clamp resultant color
    color.clamp(0.0,1.0);

    // view direction (direction from hit point to camera)
    V = -dir;

    // diffuse component from each light...
    if(!diffuse.iszero()) {
      for(litr = model.lgts.begin(); litr != model.lgts.end(); litr++) {
        // pointer to light
        lgt = (light_t *)*litr;
  
        // light direction and distance
        L = lgt->getlocation() - hit;
        r = L.len();				// distance to light
        L = L.norm();				// dir to light
  
        // angle with light
        ndotl = N.dot(L);
  
        // check visibility wrt light
        if(0.0 < ndotl && ndotl < 1.0) {
          // light color scaled by N . L
          I_d = 1.0/r * ndotl * lgt->getcolor();
  
          // add in diffuse contribution from light scaled by material property
          color += (I_d * diffuse);
        }
      }
  
      // clamp resultant color
      color.clamp(0.0,1.0);
    }

    // reflection
    if(!specular.iszero()) {
	rgb_t<double>	refcolor;
	vec_t		r = dir.reflect(N);
        ray_t		*reflection = new ray_t(hit,r,dis);

      // trace the reflection
      reflection->trace(model,refcolor, (bounce + 1));
      delete reflection;

      // composite surface color: blending baesd on surface properties
      // (note that diffuse + specular should add up to 1)
      color  = (diffuse * color) + (specular * refcolor);

      // clamp resultant color
      color.clamp(0.0,1.0);
    }

    // transmission
    if(alpha > 0.0) {
	rgb_t<double>	transmitted_color;
	vec_t		t;

        // decide which N about which to refract
	if(dir.dot(N) < .00001)
		{ t = dir.refract(N, ior); }
	else
		{ t = dir.refract( -N, (1/ior) ); }; 

       ray_t		*refraction = new ray_t(hit,t,dis);

      // trace the refraction
      refraction->trace(model,transmitted_color, (bounce + 1));
      delete refraction;

      // composite surface color: blending baesd on surface properties
      // (note that diffuse + specular should add up to 1)
      color  = (color * (1.0 - alpha)) + (transmitted_color * alpha);

      // clamp resultant color
      color.clamp(0.0,1.0);
    }

    // specular highlights
    if(!specular.iszero()) {

      // add in specular highlight from each light...
      for(litr = model.lgts.begin(); litr != model.lgts.end(); litr++) {
        // pointer to light
        lgt = (light_t *)*litr;

        // light direction and distance
        L = lgt->getlocation() - hit;
        r = L.len();				// distance to light
        L = L.norm();				// dir to light

        // angle with light
        ndotl = N.dot(L);

        // check visibility wrt light
        if(0.0 < ndotl && ndotl < 1.0) {
          // specular reflection direction (and bisector)
          R = L.reflect(N);
          // bisector (not used)
          H = (0.5 * (L + V)).norm();

          // light color scaled by (R . V)^n
          I_s = 1.0/r * pow(R.dot(V),n) * lgt->getcolor();
  
          // add in specular contribution from light scaled by material property
          color += (I_s * specular);
        }

        // clamp resultant color
        color.clamp(0.0,1.0);
      }
    }
  }
}


void ray_t::trace(model_t&		model,
                  rgb_t<double>&	color,
		  int bounce,
		  kdtree_t<photon_t, photon_t*, photon_c> kdtree)
{
	// (static) object properties
	object_t			*obj=NULL;
	rgb_t<double>			ambient, diffuse, specular;
	material_t			*mat;
	vec_t				hit,N;	   // hit point and normal
	float				alpha, ior;

	// light and illumination model variables
	light_t				*lgt=NULL;
	vec_t				L,V,R,H;
	rgb_t<double>			I_d, I_s;
	double				r,ndotl=0.0,n=32.0;
	list_t<light_t* >::iterator	litr;

	double				radius(0.0);
	std::vector<photon_t* >		knearest;

// prevent infinite loops
  if(bounce > 5) return;

  // get closest object, if any
  if(!(obj = model.find_closest(pos,dir,dis,hit,N)) || dis > MAX_DIST)
    return;

  // if hit distance valid, compute color at surface
  if(dis > 0) {

    // get object material properties
    if((mat = model.getmaterial(obj->getmaterial())) != NULL) {
      ambient = mat->getamb();
      diffuse = mat->getdiff();
      specular = mat->getspec();
      alpha = mat->getalpha();
      ior = mat->getior();
    }

    // ambient color
    color += 1.0/dis * ambient;			// ambient scaled by ray dist

    // clamp resultant color
    color.clamp(0.0,1.0);

    // view direction (direction from hit point to camera)
    V = -dir;

    // diffuse component from each light...
    if(!diffuse.iszero()) {
      for(litr = model.lgts.begin(); litr != model.lgts.end(); litr++) {
        // pointer to light
        lgt = (light_t *)*litr;
  
        // light direction and distance
        L = lgt->getlocation() - hit;
        r = L.len();				// distance to light
        L = L.norm();				// dir to light
  
        // angle with light
        ndotl = N.dot(L);
  
        // check visibility wrt light
        if(0.0 < ndotl && ndotl < 1.0) {
          // light color scaled by N . L
          I_d = 1.0/r * ndotl * lgt->getcolor();
  
          // add in diffuse contribution from light scaled by material property
          color += (I_d * diffuse);
        }
      }
  
      // clamp resultant color
      color.clamp(0.0,1.0);
    }

    // reflection
    if(!specular.iszero()) {
	rgb_t<double>	refcolor;
	vec_t		r = dir.reflect(N);
        ray_t		*reflection = new ray_t(hit,r,dis);

      // trace the reflection
      reflection->trace(model,refcolor, (bounce + 1), kdtree);
      delete reflection;

      // composite surface color: blending baesd on surface properties
      // (note that diffuse + specular should add up to 1)
      color  = (diffuse * color) + (specular * refcolor);

      // clamp resultant color
      color.clamp(0.0,1.0);
    }

    // transmission
    if(alpha > 0.0) {
	rgb_t<double>	transmitted_color;
	vec_t		t;

        // decide which N about which to refract
	if(dir.dot(N) < .00001)
		{ t = dir.refract(N, ior); }
	else
		{ t = dir.refract( -N, (1/ior) ); }; 

       ray_t		*refraction = new ray_t(hit,t,dis);

      // trace the refraction
      refraction->trace(model,transmitted_color, (bounce + 1), kdtree);
      delete refraction;

      // composite surface color: blending baesd on surface properties
      // (note that diffuse + specular should add up to 1)
      color  = (color * (1.0 - alpha)) + (transmitted_color * alpha);

      // clamp resultant color
      color.clamp(0.0,1.0);
    }

    // specular highlights
    if(!specular.iszero()) {

      // add in specular highlight from each light...
      for(litr = model.lgts.begin(); litr != model.lgts.end(); litr++) {
        // pointer to light
        lgt = (light_t *)*litr;

        // light direction and distance
        L = lgt->getlocation() - hit;
        r = L.len();				// distance to light
        L = L.norm();				// dir to light

        // angle with light
        ndotl = N.dot(L);

        // check visibility wrt light
        if(0.0 < ndotl && ndotl < 1.0) {
          // specular reflection direction (and bisector)
          R = L.reflect(N);
          // bisector (not used)
          H = (0.5 * (L + V)).norm();

          // light color scaled by (R . V)^n
          I_s = 1.0/r * pow(R.dot(V),n) * lgt->getcolor();
  
          // add in specular contribution from light scaled by material property
          color += (I_s * specular);
        }

        // clamp resultant color
        color.clamp(0.0,1.0);
      }
    }
  }
/*
	// flux computation
        photon_t *query = new photon_t(hit, vec_t(0.0,0.0,0.0), 0.0);
	kdtree.knn(*query, knearest, radius, 20);
	
	rgb_t<double> flux(0.0,0.0,0.0);
	rgb_t<double> tempcolor(0.0,0.0,0.0);

	vec_t	temp_pow;
	vec_t	temp_dir;


	for(int i = 0; i < 20; i++)
	{
		temp_pow = knearest[i]->get_power();
		temp_dir = knearest[i]->get_dir();
		temp_dir = temp_dir.norm();

		temp_dir *= vec_t(-1,-1,-1);

		if(temp_dir.dot(N) > 0)
		{	
			tempcolor[0] = temp_pow[0];
			tempcolor[1] = temp_pow[1];
			tempcolor[2] = temp_pow[2];
		}
		else
			{ tempcolor[0] = 0; tempcolor[1] = 0; tempcolor[2] = 0; }

		flux = flux + tempcolor;
	}
	
	flux[0] = flux[0] * (1 / (3.14 * pow(radius,2)));
	flux[1] = flux[1] * (1 / (3.14 * pow(radius,2)));
	flux[2] = flux[2] * (1 / (3.14 * pow(radius,2)));

	color = color + flux;
	color.clamp(0.0,1.0);
*/
}
