#include <iostream>

#include "vector.h"
#include "pixel.h"
#include "light.h"
#include "material.h"
#include "object.h"
#include "list.h"

template <typename T>
typename list_t<T>::const_iterator list_t<T>::insert(list_t<T>::const_iterator itr,const T& data)
{
	link_t	*p = itr.curp;

  // insert data before itr
  return const_iterator(p->prev = p->prev->next = new link_t(data,p->prev,p));
}

template <typename T>
typename list_t<T>::iterator list_t<T>::erase(list_t<T>::iterator itr)
{
	link_t				*p = itr.curp;
	typename list_t<T>::iterator	ret(p->next);

  // erase item at itr
  p->prev->next = p->next;
  p->next->prev = p->prev;

  delete p;

  return ret;
}

/////////////////////////// specializations /////////////////////////////////// 
template class list_t<light_t *>;
template class list_t<material_t *>;
template class list_t<object_t *>;
