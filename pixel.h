#ifndef PIXEL_H
#define PIXEL_H

#ifndef uchar
typedef unsigned char uchar;
#endif

// forward declarations
template <typename T> class rgb_t;
template <typename T> std::ostream& operator<<(std::ostream&, const rgb_t<T>&);
template <typename T> std::ostream& operator<<(std::ostream&,       rgb_t<T>*);
template <typename T> std::istream& operator>>(std::istream&,       rgb_t<T>&);
template <typename T> std::istream& operator>>(std::istream&,       rgb_t<T>*);

template <typename T> rgb_t<T> operator*(const rgb_t<T>& rhs, const T& s);
template <typename T> rgb_t<T> operator*(const T& s, const rgb_t<T>& rhs);
template <typename T> rgb_t<T> operator-(const T& s, const rgb_t<T>& rhs);

template <typename T>
class rgb_t
{
  public:
  // constructors (overloaded)
  rgb_t(T ix=0, T iy=0, T iz=0) \
  { pix[0] = ix; pix[1] = iy; pix[2] = iz; };

  // copy constructor
  rgb_t(const rgb_t& rhs) \
  { for(int i=0; i<3; i++) pix[i] = rhs[i]; }

  // destructors (but default ok)
  ~rgb_t()	{ };

  // assignment operator
  rgb_t& operator=(const rgb_t& rhs)
  {
    if(this != &rhs) for(int i=0; i<3; i++) pix[i] = rhs[i];

    return *this;
  }

  // compound assignment operator, from:
  // http://www.cs.caltech.edu/courses/cs11/material/cpp/donnie/cpp-ops.html
  rgb_t& operator*=(const rgb_t& rhs)
  {
    if(this != &rhs) for(int i=0; i<3; i++) pix[i] *= rhs[i];

    return *this;
  }

  rgb_t& operator+=(const rgb_t& rhs)
  {
    if(this != &rhs) for(int i=0; i<3; i++) pix[i] += rhs[i];

    return *this;
  }

  rgb_t& operator-=(const rgb_t& rhs)
  {
    if(this != &rhs) for(int i=0; i<3; i++) pix[i] -= rhs[i];

    return *this;
  }

  rgb_t& operator*=(const T& s)
  {
    for(int i=0; i<3; i++) pix[i] *= s;

    return *this;
  }

  rgb_t& operator+=(const T& s)
  {
    for(int i=0; i<3; i++) pix[i] += s;

    return *this;
  }

  rgb_t& operator-=(const T& s)
  {
    for(int i=0; i<3; i++) pix[i] -= s;

    return *this;
  }

  // accessors
  const T& operator[](int i) const  { return pix[i]; }
  T& operator[](int i)              { return pix[i]; }

  // add this instance's value to rhs, return a new instance with result, from:
  // http://www.cs.caltech.edu/courses/cs11/material/cpp/donnie/cpp-ops.html
  const rgb_t operator*(const rgb_t& rhs) const {
    return rgb_t(*this) *= rhs;
  }

  const rgb_t operator+(const rgb_t& rhs) const {
    return rgb_t(*this) += rhs;
  }
 
  const rgb_t operator*(const T& s) const {
    return rgb_t(*this) *= s;
  }

  const rgb_t operator+(const T& s) const {
    return rgb_t(*this) += s;
  }

  const rgb_t operator-(const T& s) const {
    return rgb_t(*this) -= s;
  }

  // friend operators (can't actually define them here for some reason)
  friend rgb_t (::operator* <>)(const rgb_t& lhs, const T& s);
  friend rgb_t (::operator* <>)(const T& s, const rgb_t& rhs);
  friend rgb_t (::operator- <>)(const T& s, const rgb_t& rhs);

  // friends
  friend std::ostream& operator<< <>(std::ostream& s, const rgb_t& rhs);
  friend std::ostream& operator<< <>(std::ostream& s, rgb_t *rhs);

  friend std::istream& operator>> <>(std::istream& s, rgb_t& rhs);
  friend std::istream& operator>> <>(std::istream& s, rgb_t *rhs);

  // members
  void zero(void) {
    pix[0] = T(0);  pix[1] = T(0);  pix[2] = T(0);
  }
  bool iszero(void) const;
  void clamp(T min, T max)
  {
    for(int i=0;i<3;i++) if(pix[i] < min) pix[i] = min;
    for(int i=0;i<3;i++) if(pix[i] > max) pix[i] = max;
  }

  private:
  T pix[3];
};

#endif
