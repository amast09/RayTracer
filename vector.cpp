#include <iostream>
#include <cmath>

#include "vector.h"

// compute the dot product of two vectors
double vec_t::dot(const vec_t& rhs) const
{
	double	s=0.0;

  for(int i=0;i<3;i++) s += vec[i] * rhs[i];

  return(s);
}

// compute the length of the vector v1
double vec_t::len() const
{
  return(sqrt(dot(*this)));
}

// compute the dot product of two vectors
vec_t vec_t::norm() const
{
	vec_t	result(*this);
	double	length = len();

  for(int i=0;i<3;i++) result[i] /= length;

  return(result);
}

// compute the reflection vector: v = u - 2 (u dot n) n
vec_t vec_t::reflect(const vec_t& n) const
{
	vec_t	u(*this);
	vec_t	result;

  // u - 2 (u dot n) n
  result = u - 2.0 * u.dot(n) * n;

  return result;
}

// compute the refraction vector
vec_t vec_t::refract(const vec_t& N, double n2) const
{
	const float n1 = 1.000293;
	double	x;
	double	z;
	vec_t  refraction;
	vec_t  u(*this);

	x = (1 - (pow((n1 / n2), 2) * (1 - pow(u.dot(N), 2))));

	if(x < 0)
		{return(u.reflect(N));}

// refraction = (n1/n2)(u - cosθ1N) - cosθ2N
	refraction = (((n1 / n2) * (u - ((u.dot(N)) * N)))  - ((sqrt(x)) * N));

	return(refraction);
}

// compute the refraction vector
vec_t vec_t::defract(const vec_t& N, double n2) const
{
	const float n1 = 1.000293;
	double	x;
	double	z;
	vec_t  refraction;
	vec_t  u(*this);

	x = (1 - (pow((n1 / n2), 2) * (1 - pow(u.dot(N), 2))));

	if(x < 0)
		{return(u.reflect(N));}

// refraction = (n1/n2)(u - cosθ1N) - cosθ2N
	refraction = (((n1 / n2) * (u + ((u.dot(N)) * N)))  - ((sqrt(x)) * N));

	return(refraction);
}

