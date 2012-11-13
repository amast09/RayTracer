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
#include "ray.h"
#include "photon.h"


vec_t photon_t::genrand_hemisphere()
{
	double  azimuth = genrand(0.0, 2.0 * M_PI);
	double  elevation = genrand(0.0, 2.0 * M_PI);

	double  sinA = sin(azimuth), sinE = sin(elevation);
	double  cosA = cos(azimuth), cosE = cos(elevation);
	vec_t   dir, vup;

	dir[0] = -sinA * cosE;  vup[0] =  sinA * sinE;
	dir[1] =  sinE;       vup[1] =  cosE;
	dir[2] =  cosA * cosE;  vup[2] = -cosA * sinE;

	return(dir + vup);
}

bool photon_t::caustic(model_t& model, int bounce)
{
	double				rdis = 0.0;
	double				alpha, ior;
	object_t			*obj = NULL;
	rgb_t<double>			ambient, diffuse, specular;
	material_t			*mat;
	vec_t				N, hit, normal;

	// prevent infinite loops
	if(bounce > 5)
		{  return(false);  }

	// get closest object, if any
	if(!(obj = model.find_closest(pos, dir, rdis, hit, normal)) || dis > MAX_DIST)
		{  return(false);  }

	// if hit distance valid
	if(rdis > 0)
	{
		// accumulate distance travelled by ray
		dis += rdis;

		// get object material properties
		if((mat = model.getmaterial(obj -> getmaterial())) != NULL)
		{
			ambient = mat -> getamb();
			diffuse = mat -> getdiff();
			specular = mat -> getspec();
			alpha = mat -> getalpha();
			ior = mat -> getior();
		}

		N = normal.norm();		// surf the normal

		// Refraction
		if(alpha > 0.0)
		{
			vec_t	out;

			out = (dir.dot(N)<0)?dir.defract(N, ior).norm():dir.defract(N*-1.0, 1/ior).norm();
			pos = hit;
			dir = out;

			return(caustic(model, (bounce + 1))); 
		}
		else
		{
			pos = hit;
			return(true);
		}
	}
}

bool photon_t::global(model_t& model, int bounce)
{
	double				rdis = 0.0;
	object_t			*obj = NULL;
	rgb_t<double>			ambient, diffuse, specular;
	material_t			*mat;
	vec_t				N, hit, normal;

	// prevent infinite loops
	if(bounce > 5)
		{  return(false);  }

	// get closest object, if any
	if(!(obj = model.find_closest(pos,dir,rdis,hit,normal)) || dis > MAX_DIST)
		{ return(false); }

	// if hit distance is valid
	if(rdis > 0)
	{

		// accumulate distance travelled by ray
		dis += rdis;

		// get object material properties
		if((mat = model.getmaterial(obj->getmaterial())) != NULL)
		{
			ambient = mat -> getamb();
			diffuse = mat -> getdiff();
			specular = mat -> getspec();
		}

		N = normal.norm();		// surf normal

		double roulette = genrand(0.0, 1.0);
		double Kd = (diffuse[0] + diffuse[1] + diffuse[2])/3.0;
		double Ks = (specular[0] + specular[1] + specular[2])/3.0;

		// diffuse reflection
		if( (0.0 < roulette) && (roulette < Kd))
		{
			// set pos and dir to hit and genrand_hemisphere().norm()
			pos = hit;
			dir = genrand_hemisphere().norm();
			// reflect the photon
			return global(model, (bounce + 1));
		}

		// specular reflection
		else if( (Kd < roulette) && (roulette < Kd + Ks))
		{
			// set pos and dir to hit and dir.reflect(N);
			pos = hit;
			dir = dir.reflect(N);
			// reflect the photon
			return global(model, (bounce + 1));
		}

		// absorbtion ("stick" photon)
		else if( (Kd + Ks < roulette) && (roulette < 1.0))
		{
			pos = hit;
			return(true);
		}
	}
}
