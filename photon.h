#ifndef PHOTON_H
#define PHOTON_H

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "ray.h"
#include "vector.h"
#include "model.h"

using namespace std;

class photon_t : public ray_t
{
	public:

	// Constructor with no parameters
	photon_t() 
	: ray_t(vec_t(0.0, 0.0, 0.0), genrand_hemisphere(), 0.0)
	{
		power = vec_t(100.0, 100.0, 100.0);
	}

	// Constructors given a power vector
	photon_t(const vec_t& p)
	:ray_t(p, genrand_hemisphere(), 0.0)
	{
		power = vec_t(100.0, 100.0, 100.0);
	}

	// Copy constructor
	photon_t(const photon_t& rhs)
	{
		power = rhs.power;
	}

	photon_t(const vec_t& x, const vec_t& y, double z = 0.0)
	:ray_t(x, y, z)
	{
		power = vec_t(100.0, 100.0, 100.0);
	}
	
	// Assignment operator
	const photon_t& operator=(const photon_t& rhs)
	{
		power = rhs.power;
		dis = rhs.dis;
		dir = rhs.dir;
		pos = rhs.pos;
		return *this;
	}


///////// Friend Functions ////////////
	
	// output stream
	friend ostream& operator<<(ostream& s, photon_t& rhs)
	{
		return(s << rhs.pos[0] << " " << rhs.pos[1] << " " << rhs.pos[2] << " "<< rhs.power[0] << " " << rhs.power[1] << " " << rhs.power[2]); 
	}
	// input stream
	friend istream& operator>>(istream& s, photon_t& rhs)
	{
		return s >> rhs.power[0] >> rhs.power[1] >> rhs.power[2];
	}

	// operator [] accessor
	const double& operator[](int i) const  { return pos[i]; }
	double& operator[](int i)              { return pos[i]; }

	// less than operator used to find the min
	friend bool operator<(const photon_t& lhs, const photon_t& rhs)
	{
		for(int i = 0; i < 3; i++)
		{
			if(lhs.pos[i] > rhs.pos[i]) { return(false); }
		}
		return(true);
	}

	// greater than operator used to find the max
	friend bool operator>(const photon_t& lhs, const photon_t& rhs)
	{
		for(int i = 0; i < 3; i++)
		{
			if(lhs.pos[i] < rhs.pos[i]) { return(false); }
		}
		return(true);
	}

///////// Member functions /////////////

	// Generate random number
	double   genrand(double lo, double hi)
	{  return( (double)(((double)rand()/(double)RAND_MAX)*hi + lo) );  }

	// Set the direction for the photon
	vec_t genrand_hemisphere();

	// takes care of caustic photons
	bool caustic(model_t& model, int bounce);

	// takes care of global photons
	bool global(model_t& model, int bounce);

	// lets others access the power member of photon
	vec_t get_power()	{ return power; }

	vec_t get_pos()		{ return pos; }

	vec_t get_dir()		{ return dir; }

	// updates the power value
	void change_power(vec_t rhs)	{ power = rhs; }

	// compute the distance between a passed photon and the "this" photon given a vec_t 
	double distance(const vec_t& rhs)
		{ vec_t   diff = pos - rhs; return(sqrt(diff.dot(diff))); }

	// compute the distance between a passed photon and the "this" photon given a photon_t
	double distance(const photon_t& rhs)
		{ vec_t   diff = pos - rhs.pos; return(sqrt(diff.dot(diff))); }

	double distance(photon_t*& rhs)
		{ vec_t   diff = pos - rhs->pos; return(sqrt(diff.dot(diff))); }

	int dim()	{ return 3; }

	private:

	vec_t power;
};


class photon_c
{
	public:

	// constructor
	photon_c(int inputaxis = 0)  {  axis = inputaxis;  }

	bool operator()(const photon_t& p1, const photon_t& p2) const
		{  return( p1[axis] < p2[axis] );  }

	bool operator()(const photon_t*& p1, const photon_t*& p2) const
		{  return( (*p1)[axis] < (*p2)[axis] ); }

	bool operator()(photon_t * const & p1, photon_t * const & p2) const
		{ return((*p1)[axis] < (*p2)[axis]); }

	private:

	int axis;
};

#endif
