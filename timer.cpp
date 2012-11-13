#include        <iostream>
#include        <string>
#include	<cstdlib>
#include	<sys/time.h>

#include	"timer.h"

namespace atd {

std::ostream& operator<<(std::ostream& s,const timer_t& rhs)
{
  s << "time elapsed: " << rhs.tt << " us" << std::endl;

  return(s);
}


void timer_t::start()
{
  ts = stamp_us();
}


void timer_t::end()
{
  te = stamp_us();
  tt = te - ts;
}


double timer_t::stamp_us()
{
	double		s,us,tod;
	struct timeval	tp;

  // get time of day (tod), return in microseconds

  gettimeofday(&tp,NULL);
  s = static_cast<double>(tp.tv_sec);
  us = static_cast<double>(tp.tv_usec);
  tod = s*1000000.0 + us;
  return(tod);
}

} // namespace atd
