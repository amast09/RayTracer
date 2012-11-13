#include <iostream>
#include <cmath>

#include "pixel.h"

template <typename T>
bool rgb_t<T>::iszero(void) const
{
	double  precision=0.0000001;

  for(int i=0;i<3;i++) if(fabs((double)pix[i]) > precision) return(false);

  return true;
}

/////////////////////////// friend operators: (const scalar&, ...) //////////// 
template <typename T>
rgb_t<T> operator*(const rgb_t<T>& lhs, const T& s)
{
	rgb_t<T>	result;

    for(int i=0; i<3; i++) result[i] = s * lhs[i];

    return result;
}

template <typename T>
rgb_t<T> operator*(const T& s, const rgb_t<T>& rhs)
{
        rgb_t<T>        result;

  for(int i=0; i<3; i++) result[i] = s * rhs[i];

  return result;
}

template <typename T>
rgb_t<T> operator-(const T& s, const rgb_t<T>& rhs)
{
        rgb_t<T>        result;

  for(int i=0; i<3; i++) result[i] = s - rhs[i];

  return result;
}

template <typename T>
std::ostream& operator<<(std::ostream& s, const rgb_t<T>& rhs)
{
  return(s << rhs[0] << " " << rhs[1] << " " << rhs[2]);
}

template <typename T>
std::ostream& operator<<(std::ostream& s, rgb_t<T> *rhs)
{
  return(s << (*rhs));
}

template <typename T>
std::istream& operator>>(std::istream& s, rgb_t<T>& rhs)
{
  return(s >> rhs[0] >> rhs[1] >> rhs[2]);
}

template <typename T>
std::istream& operator>>(std::istream& s, rgb_t<T> *rhs)
{
  return(s >> (*rhs));
}

/////////////////////////// specializations /////////////////////////////////// 
template class rgb_t<double>;
template std::ostream& operator<<(std::ostream&, const rgb_t<double>&);
template std::ostream& operator<<(std::ostream&,       rgb_t<double>*);
template std::istream& operator>>(std::istream&,       rgb_t<double>&);
template std::istream& operator>>(std::istream&,       rgb_t<double>*);

template rgb_t<double> operator*(const rgb_t<double>& rhs, const double& s);
template rgb_t<double> operator*(const double& s, const rgb_t<double>& rhs);
template rgb_t<double> operator-(const double& s, const rgb_t<double>& rhs);

template class rgb_t<uchar>;
template std::ostream& operator<<(std::ostream&, const rgb_t<uchar>&);
template std::ostream& operator<<(std::ostream&,       rgb_t<uchar>*);
template std::istream& operator>>(std::istream&,       rgb_t<uchar>&);
template std::istream& operator>>(std::istream&,       rgb_t<uchar>*);

template rgb_t<uchar> operator*(const rgb_t<uchar>& rhs, const uchar& s);
template rgb_t<uchar> operator*(const uchar& s, const rgb_t<uchar>& rhs);
template rgb_t<uchar> operator-(const uchar& s, const rgb_t<uchar>& rhs);
