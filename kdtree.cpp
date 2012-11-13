#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
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
#include "timer.h"
#include "photon.h"
#include "kdtree.h"

	// outputs the entire tree
	template <typename T, typename P, typename C>
	std::ostream& operator<<(std::ostream& s, const kdtree_t<T, P, C>& rhs)
	{
		s.setf(ios::fixed,ios::floatfield);
		s.precision(2);

		if(rhs.empty())
			{ s << "tree is empty" << std::endl; }
		else
			{ s << rhs.root; }
		return s;
	}

	// deletes the entire tree through recursion, by deleting the left and right
	// nodes of each node before deleteing the actual node, starting from the root
	template <typename T, typename P, typename C>
	void kdtree_t<T,P,C>::clear(kdnode_t* &t)
	{
		if(t != NULL)
		{
			clear(t->left);
			clear(t->right);
			delete t;
		}
		t = NULL;
	}

	// clones tree by recurrisvely cloning the left and right subtrees of the passed node
	template <typename T, typename P, typename C>
	typename kdtree_t<T,P,C>::kdnode_t* kdtree_t<T,P,C>::clone(kdnode_t* nod) const
	{
		if(nod == NULL)	{ return NULL; }
		return( new kdnode_t(nod->data, nod->min, nod->max, clone(nod->left), clone(nod->right), nod->axis) );
	}

	template <typename T, typename P, typename C>
	typename kdtree_t<T,P,C>::kdnode_t* kdtree_t<T,P,C>::insert(std::vector<P>& x, const T& min, const T& max, int d)
	{
		int					axis = x.empty() ? 0 : d % x[0]->dim();
		int					m = 0; // the index for the median
		P					median;
		T					_min, _max; // bounding values of subspaces
		std::vector<P>				left, right;
		typename std::vector<P>::iterator 	itr;

		// return if there are no points in the vector
		if(x.empty()) { return(NULL); }

		// debugging
		std::cerr << "depth: " << d << std::endl;
		std::cerr << "size: " << x.size() << std::endl;

		// sort the vector of points to obtain the median
		sort(x.begin(), x.end(), C(axis));

		// debugging
		for(itr = x.begin(); itr < x.end(); itr++)  { std::cerr << (**itr) << std::endl; }

		// get actual median
		m = x.size() / 2;

		// create left and right subtrees
		for(int i = 0; i < (int)x.size(); i++)
		{
			if(i < m)	{ left.push_back(x[i]); }
			else if(i > m)	{ right.push_back(x[i]); }
			else		{ median = x[m]; }
		}

		// create a new node
		kdnode_t* node = new kdnode_t(median, min, max, NULL, NULL, axis);

		// recursively add the left subtree
		_min = min;
		_max = max;
		_max[axis] = (*median)[axis];
		node->left = insert(left, _min, _max, (d + 1));

		// recursively add the right subtree
		_min = min;
		_max = max;
		_max[axis] = (*median)[axis];
		node->right = insert(right, _min, _max, (d + 1));

		return(node);
		
	}

	// find the nearest point to a point passed in the parameters
	template <typename T, typename P, typename C>
	void kdtree_t<T,P,C>::nn(kdnode_t* &nod, T& center, P& pt, double& radius)
	{
		double	dist;
		int	axis;

		if(nod == NULL)	{ return; }

		// determine the node's distance to the point passed as a parameter
		dist = center.distance(nod->data);

		if(dist < radius)
		{
			radius = dist;
			pt = nod->data;
		}

		// traverse down "closer" side of the tree
		// testing each encountered node as we move
		// down the tree the recursive calls will
		// override the pt and radius arguments
		// if a closer node than the current node is found
		//
		// also need to check to see if the current closest
		// circle defined by the center and radius variables
		// intersects the other side of the tree in which
		// case we must search that subtree
		axis = nod->axis;
		if(center[axis] <= (*nod->data)[axis])
		{
			nn(nod->left, center, pt, radius);
			if((center[axis] + radius) > (*nod->data)[axis])
				{  nn(nod->right, center, pt, radius); }
		}
		else
		{
			nn(nod->right, center, pt, radius);
			if((center[axis] - radius) >= (*nod->data)[axis])
				{ nn(nod->left, center, pt, radius); }		
		}
	}

	// find k nearest points to the point passed as a parameter
	template <typename T, typename P, typename C>
	void kdtree_t<T,P,C>::knn(kdnode_t* &nod, T& center, std::vector<P>& pt, double& radius, int k)
	{
		double					dist;
		int					axis;
		typename std::vector<P>::iterator	pitr;

		if(nod == NULL)	{ return; }

		// determine the node's distance to the point that is passed as a parameter
		dist = center.distance(nod->data);

		// instead of searching a set raduis, search a radius
		// of infinity until k points are found
		if((int)pt.size() < k)
		{
			// if the list is empty or the distance
			// is larger than last node add the node to the end
			if(pt.empty() || (dist > center.distance( pt.back())))
				{ pt.push_back(nod->data); }
			// else iterate through the list to find proper place to insert the point
			else
			{
				for(pitr = pt.begin(); pitr != pt.end(); pitr++)
				{
					if(dist < center.distance(*pitr))
						{ pt.insert(pitr, 1, nod->data); break; }
				}
			}
		}
		// else insert the current node into the list of closest points
		// checking its distance compared to points that are in the list
		// insert the point only if the distance is smaller than the 
		// last point on the list if it isn't just ignore the point
		else
		{
			if(dist < center.distance( pt.back()))
			{
				for(pitr = pt.begin(); pitr != pt.end(); pitr++)
				{
					if(dist < center.distance(*pitr))
						{ pt.insert(pitr, 1, nod->data); break; }
				}
			}
			// since we already have k total points and have added
			// one node we only need to pop one point off of the
			// back of the list of points
			if((int)pt.size() > k)	{ pt.pop_back(); }
		}

		// find the largest distance in the list of points
		radius = center.distance(pt.back());

		// exact same algorithm as nn except recursievly run knn instead of nn
		axis = nod->axis;
		if(center[axis] <= (*nod->data)[axis])
		{
			knn(nod->left, center, pt, radius, k);
			if((center[axis] + radius) > (*nod->data)[axis])
				{  knn(nod->right, center, pt, radius, k); }
		}
		else
		{
			knn(nod->right, center, pt, radius, k);
			if((center[axis] - radius) >= (*nod->data)[axis])
				{ knn(nod->left, center, pt, radius, k); }		
		}

	}

	template <typename T, typename P, typename C>
	void kdtree_t<T,P,C>::range(kdnode_t* &nod, const T& min, const T& max, std::vector<P>& pt)
	{
		bool	truth_value;

		if(nod == NULL)	{ return; }

		// checks to see if the data is in the range
		for(int i = 0; i < (nod->data->dim()); i++)
		{
			if(((*nod->data)[i] < min[i]) || ((*nod->data)[i] > max[i]))
				{ truth_value = false; break; }
		}

		// if data is in range it is added to the list
		if(truth_value == true)		{ pt.push_back(nod->data); }

		// the range is on both sides of the tree
		if(	((*nod->data)[nod->axis] >= min[nod->axis]) &&
			((*nod->data)[nod->axis] <= max[nod->axis])	)
		{
			range(nod->left,min,max,pt);
			range(nod->right,min,max,pt);
		}

		// the range is all smaller and can be placed completely on the left of the tree
		else if((*nod->data)[nod->axis] >= max[nod->axis])
			{ range(nod->left,min,max,pt); }

		// the range is all larger and can be placed completely on the right of the tree
		else if((*nod->data)[nod->axis] <= min[nod->axis])
			{ range(nod->right,min,max,pt); }
	}


	///////// specializations //////////
	template class kdtree_t<photon_t, photon_t*, photon_c>;	
	template std::ostream& operator<<(std::ostream&, const kdtree_t<photon_t, photon_t*, photon_c>&);


