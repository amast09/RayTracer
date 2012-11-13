#ifndef LIST_H
#define LIST_H

// forward declarations
template <typename T> class list_t;

template <typename T>
class list_t
{
  private:
  struct link_t
  {
	T	data;
	link_t	*prev;
	link_t	*next;

	link_t(const T& d=T(), link_t *p = NULL, link_t *n = NULL) : \
		data(d), prev(p), next(n) \
						{ };
  };

  public:
  // returns const reference when de-referenced
  class const_iterator
  {
	public:
	const_iterator() : curp(NULL)	{ };

	const T& operator*() const		// accessor
	{ return retrieve(); }
	const_iterator& operator++()		// prefix ++itr
	{ curp = curp->next; return *this; }
	const_iterator operator++(int)		// postfix itr++
	{ const_iterator old(curp); ++(*this); return old; }
	bool operator==(const const_iterator& rhs) const
	{ return curp == rhs.curp; }
	bool operator!=(const const_iterator& rhs) const
	{ return !(curp == rhs.curp); }

	protected:
	link_t	*curp;

	T& retrieve() const
	{ return curp->data; }

	const_iterator(link_t *p) : curp(p)	{ };

	friend class list_t<T>;
  };

  // returns reference when de-referenced
  class iterator : public const_iterator
  {
	public:
	iterator() 				{ };

	T& operator*()				// mutator
	{ return iterator::retrieve(); }
	const T& operator*() const		// accessor
	{ return const_iterator::operator*(); }
	iterator& operator++()			// prefix ++itr
	{ iterator::curp = iterator::curp->next; return *this; }
	iterator operator++(int)		// postfix itr++
	{ iterator old(*this); ++(*this); return old; }

	protected:
	iterator(link_t *p) : const_iterator(p)	{ };

	friend class list_t<T>;
  };

  public:
  // constructors (overloaded)
  list_t()
  {
    head = new link_t;
    tail = new link_t;
    head->next = tail;
    tail->prev = head;
  }

  // copy constructor
  list_t(const list_t& rhs)
  {
    head = new link_t;
    tail = new link_t;
    head->next = tail;
    tail->prev = head;
    *this = rhs;
  }

  // destructors
  ~list_t()	{ clear(); delete head; delete tail; };

  // operators (incl. assignment operator)
  const list_t& operator=(const list_t& rhs)
  {
    if(this == &rhs) return *this;	// standard alias test

    clear();				// clear this list

    // copy all elements
    for(const_iterator itr = rhs.begin(); itr != rhs.end(); ++itr)
      insert(itr,*itr);

    return *this;
  }

  // iterator functions
  const_iterator	begin() const	{ return const_iterator(head->next); }
  const_iterator	end() const	{ return const_iterator(tail); }
  iterator		begin()		{ return iterator(head->next); }
  iterator		end()		{ return iterator(tail); }

  const_iterator	insert(const_iterator,const T&);
  const_iterator	push_back(const T& data) { return insert(end(),data); }
  iterator		erase(iterator);
  iterator		pop_front()	{ return erase(begin()); }

  // members
  bool   		empty() const	{ return head->next == tail; }

  void			clear()		{ while(!empty()) pop_front(); };

  // data members
  private:
  link_t  *head;		// head node
  link_t  *tail;		// tail node
};

#endif
