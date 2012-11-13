#ifndef KDTREE_H
#define KDTREE_H

#ifndef INFINITY
#define INFINITY MAXFLOAT
#endif


// forward declarations
template <typename T, typename P, typename C> class kdtree_t;
template <typename T, typename P, typename C>
        std::ostream& operator<<(std::ostream&, const kdtree_t<T,P,C>&);

template <typename T, typename P, typename C>

class kdtree_t
{
	private:

////////// Node /////////////////////
	struct kdnode_t
	{
		P        	data;
		T		min, max;
		kdnode_t	*left, *right;
		int		axis;

		kdnode_t(const P& dat = P(), const T& inmin = T(), const T& inmax = T(),
			kdnode_t *l = NULL, kdnode_t *r = NULL, int x = 0) : \
				data(dat), min(inmin), max(inmax), left(l), right(r), axis(x) \
									{ };

		// output stream operators of node
		friend std::ostream& operator<<(std::ostream &s, const kdnode_t& rhs)
		{
			s << rhs.left;
			s << rhs.data;
			s << std::endl;
			s << rhs.right;
			return(s);
		}

		friend std::ostream& operator<<(std::ostream &s, kdnode_t *rhs)
		{ 
			if(rhs) 	{ return(s << (*rhs)); }
			else		{return(s); }
		}
	};
////////////////////////////////////

	public:

	// constructors
	kdtree_t()	{ root = NULL; };

	// destructors
	~kdtree_t()	{ clear(); };

	// assignment operator
	kdtree_t operator=(const kdtree_t& rhs)
	{
		if(this != &rhs)
		{
			clear();		
			root = clone(rhs.root);
		}
		return(*this);
	}

	// output stream operators
	friend std::ostream& operator<< <>(std::ostream& s, const kdtree_t& rhs);
	friend std::ostream& operator<<(std::ostream& s, const kdtree_t *rhs)
		{ return (s << (*rhs)); };

////////// member functions //////////
	
	// checks to see if the tree is empty
	bool		empty() const
				{ return root == NULL ? true : false; }

	// inserts a point into the kd tree
	kdnode_t*	insert(std::vector<P>& pt, const T& min, const T& max)
				{ root = insert(pt, min, max, 0); }

	// finds the nearest neighbor to a given point (pt)
	void		nn(T& center, P& pt, double& radius)
				{ radius = INFINITY; nn(root, center, pt, radius); }

	// finds k nearest neighbors to a given point (pt)
	void		knn(T& center, std::vector<P> pt, double& radius, int k)
				{ radius = INFINITY; knn(root, center, pt, radius, k); }

	// creates a list of data that is within the min and max passed by the function
	void		range(const T& min, const T& max, std::vector<P>& pt)
				{ range(root, min, max, pt); }

	// clears the entire kd tree
	void		clear()		{ clear(root); }

	private:

	kdnode_t	*root;

	kdnode_t*	insert(std::vector<P>&, const T&, const T&, int);
	void		nn(kdnode_t* &, T&, P&, double&);
	void		knn(kdnode_t* &, T&, std::vector<P>&, double&, int);
	void		range(kdnode_t* &, const T&, const T&, std::vector<P>&);
	void		clear(kdnode_t* &);
	kdnode_t*	clone(kdnode_t* ) const;

};

#endif
