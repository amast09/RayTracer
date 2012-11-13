#ifndef VECTOR_H
#define VECTOR_H

class vec_t
{
  public:
  // constructors (overloaded)
  vec_t(double ix=0.0, double iy=0.0, double iz=0.0) \
  { vec[0] = ix; vec[1] = iy; vec[2] = iz; };

  // copy constructor
  vec_t(const vec_t& rhs) \
  { for(int i=0; i<3; i++) vec[i] = rhs[i]; }

  // destructors (but default ok)
  ~vec_t()	{ };

  // assignment operator
  vec_t& operator=(const vec_t& rhs)
  {
    if(this != &rhs) for(int i=0; i<3; i++) vec[i] = rhs[i];

    return *this;
  }

  // compound assignment operator +=
  // style notes from:
  // http://www.cs.caltech.edu/courses/cs11/material/cpp/donnie/cpp-ops.html
  vec_t& operator+=(const vec_t& rhs)
  {
    if(this != &rhs) for(int i=0; i<3; i++) vec[i] += rhs[i];

    return *this;
  }

  vec_t& operator*=(const vec_t& rhs)
  {
    if(this != &rhs) for(int i=0; i<3; i++) vec[i] *= rhs[i];

    return *this;
  }

  // accessors
  const double& operator[](int i) const  { return vec[i]; }
  double& operator[](int i)              { return vec[i]; }

  // unary operator
  vec_t operator-() const
  {
	vec_t	result;

    for(int i=0; i<3; i++) result[i] = -vec[i];

    return result;
  }

  // operators
  vec_t operator-(const vec_t& rhs) const
  {
	vec_t	result;

    for(int i=0; i<3; i++) result[i] = vec[i] - rhs[i];

    return result;
  }

  // add this instance's value to rhs, and return a new instance with result
  // from:
  // http://www.cs.caltech.edu/courses/cs11/material/cpp/donnie/cpp-ops.html
  const vec_t operator+(const vec_t& rhs) const {
    return vec_t(*this) += rhs;
  }
 
  const vec_t operator*(const vec_t& rhs) const {
    return vec_t(*this) *= rhs;
  }

  friend vec_t operator*(const vec_t& lhs, const double& s)
  {
	vec_t	result;

    for(int i=0; i<3; i++) result[i] = s * lhs[i];

    return result;
  }

  friend vec_t operator*(const double& s, const vec_t& rhs)
  {
	vec_t	result;

    for(int i=0; i<3; i++) result[i] = s * rhs[i];

    return result;
  }

  // friends
  friend std::ostream& operator<<(std::ostream& s, const vec_t& rhs)
  { return(s << rhs[0] << " " << rhs[1] << " " << rhs[2]); }
  friend std::ostream& operator<<(std::ostream& s, vec_t *rhs)
  { return(s << (*rhs)); }

  friend std::istream& operator>>(std::istream& s, vec_t& rhs)
  { return(s >> rhs[0] >> rhs[1] >> rhs[2]); }
  friend std::istream& operator>>(std::istream& s, vec_t *rhs)
  { return(s >> (*rhs)); }

  // members
  double dot(const vec_t&) const;
  double len() const;
  vec_t	 norm() const;
  vec_t  reflect(const vec_t&) const;
  vec_t  refract(const vec_t& n, double n_t) const;
  vec_t  defract(const vec_t& n, double n_t) const;

  private:
  double vec[3];
};

#endif
