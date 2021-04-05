#ifndef E_INTERVAL_TREE
#define E_INTERVAL_TREE

#include "GBase.h"
#include "GVec.hh"

// This is an interval tree implementation based on red-black-trees
// as described in the book _Introduction_To_Algorithms_ by Cormen, Leisserson, and Rivest.

class GIntervalTreeNode {
  friend class GIntervalTree;
protected:
  GSeg* storedInterval;
  int key;
  int high;
  int maxHigh;
  int red; /* if red=0 then the node is black */
  GIntervalTreeNode* left;
  GIntervalTreeNode* right;
  GIntervalTreeNode* parent;
public:
  void Print(GIntervalTreeNode* nil,
	     GIntervalTreeNode* root) const {
	  printf(", k=%i, h=%i, mH=%i",key,high,maxHigh);
	  printf("  l->key=");
	  if( left == nil) printf("NULL"); else printf("%i",left->key);
	  printf("  r->key=");
	  if( right == nil) printf("NULL"); else printf("%i",right->key);
	  printf("  p->key=");
	  if( parent == root) printf("NULL"); else printf("%i",parent->key);
	  printf("  red=%i\n",red);
  }
  GIntervalTreeNode():storedInterval(NULL), key(0), high(0),maxHigh(0),red(0),
		  left(NULL), right(NULL), parent(NULL) {}
  GIntervalTreeNode(GSeg * newInterval): storedInterval (newInterval),
      key(newInterval->start), high(newInterval->end) ,
      maxHigh(high), red(0), left(NULL), right(NULL), parent(NULL) {  }
  ~GIntervalTreeNode() {}
};

struct G_ITRecursionNode {
public:
  // this structure stores the information needed when we take the
  // right branch in searching for intervals but possibly come back
  // and check the left branch as well.
  GIntervalTreeNode * start_node;
  unsigned int parentIndex;
  int tryRightBranch;
} ;

class GIntervalTree {
private:
  unsigned int recursionNodeStackSize;
  G_ITRecursionNode * recursionNodeStack;
  unsigned int currentParent;
  unsigned int recursionNodeStackTop;
protected:
  // A sentinel is used for root and for nil.  root->left should always
  // point to the node which is the root of the tree.  nil points to a
  // node which should always be black but has arbitrary children and
  // parent and no key or info; These sentinels are used so
  // that the root and nil nodes do not require special treatment in the code
  GIntervalTreeNode* root;
  GIntervalTreeNode* nil;

  //  INPUT:  the node to rotate on
  //  Rotates as described in _Introduction_To_Algorithms by
  //  Cormen, Leiserson, Rivest (Chapter 14).  Basically this
  //  makes the parent of x be to the left of x, x the parent of
  //  its parent before the rotation and fixes other pointers
  //  accordingly. Also updates the maxHigh fields of x and y
  //  after rotation.
  void LeftRotate(GIntervalTreeNode* x) {
	  GIntervalTreeNode* y;
	  //  originally wrote this function to use the sentinel for
	  //  nil to avoid checking for nil.  However this introduces a
	  //  very subtle bug because sometimes this function modifies
	  //  the parent pointer of nil.  This can be a problem if a
	  //  function which calls LeftRotate also uses the nil sentinel
	  //  and expects the nil sentinel's parent pointer to be unchanged
	  //  after calling this function.  For example, when DeleteFixUP
	  //  calls LeftRotate it expects the parent pointer of nil to be
	  //  unchanged.

	  y=x->right;
	  x->right=y->left;

	  if (y->left != nil) y->left->parent=x; // used to use sentinel here
	  // and do an unconditional assignment instead of testing for nil

	  y->parent=x->parent;

	  // instead of checking if x->parent is the root as in the book, we
	  // count on the root sentinel to implicitly take care of this case
	  if( x == x->parent->left) {
	    x->parent->left=y;
	  } else {
	    x->parent->right=y;
	  }
	  y->left=x;
	  x->parent=y;

	  x->maxHigh=GMAX(x->left->maxHigh, GMAX(x->right->maxHigh,x->high));
	  y->maxHigh=GMAX(x->maxHigh,GMAX(y->right->maxHigh,y->high));
  }
  //    make the parent of x be to the left of x, x the parent of
  //    its parent before the rotation and fixes other pointers
  //    accordingly. Also updates the maxHigh fields of x and y
  //    after rotation.
  void RightRotate(GIntervalTreeNode*y) {
	  GIntervalTreeNode* x;

	  x=y->left;
	  y->left=x->right;

	  if (nil != x->right)  x->right->parent=y; //used to use sentinel here
	  // and do an unconditional assignment instead of testing for nil

	  // instead of checking if x->parent is the root as in the book, we
	  // count on the root sentinel to implicitly take care of this case
	  x->parent=y->parent;
	  if( y == y->parent->left) {
	    y->parent->left=x;
	  } else {
	    y->parent->right=x;
	  }
	  x->right=y;
	  y->parent=x;

	  y->maxHigh=GMAX(y->left->maxHigh,GMAX(y->right->maxHigh,y->high));
	  x->maxHigh=GMAX(x->left->maxHigh,GMAX(y->maxHigh,x->high));

	}

  // Inserts z into the tree as if it were a regular binary tree
  //   using the algorithm described in _Introduction_To_Algorithms_
  //   by Cormen et al.  This function is only intended to be called
  //    by the InsertTree function and not by the user
  void TreeInsertHelp(GIntervalTreeNode* z) {
	  //  this should only be called by the Insert method
	  GIntervalTreeNode* x;
	  GIntervalTreeNode* y;

	  z->left=z->right=nil;
	  y=root;
	  x=root->left;
	  while( x != nil) {
	    y=x;
	    if ( x->key > z->key) {
	      x=x->left;
	    } else { // x->key <= z->key
	      x=x->right;
	    }
	  }
	  z->parent=y;
	  if ( (y == root) ||
	       (y->key > z->key) ) {
	    y->left=z;
	  } else {
	    y->right=z;
	  }
	#if defined(DEBUG_ASSERT)
	  Assert(!nil->red,"nil not red in ITTreeInsertHelp");
	  Assert((nil->maxHigh=MIN_INT),
		 "nil->maxHigh != MIN_INT in ITTreeInsertHelp");
	#endif
	}

  void TreePrintHelper(GIntervalTreeNode* x) const {
	  if (x != nil) {
	    TreePrintHelper(x->left);
	    x->Print(nil,root);
	    TreePrintHelper(x->right);
	  }
	}

  //  FUNCTION:  FixUpMaxHigh
  //  INPUTS:  x is the node to start from
  //  EFFECTS:  Travels up to the root fixing the maxHigh fields after
  //            an insertion or deletion
  void FixUpMaxHigh(GIntervalTreeNode* x) {
	  while(x != root) {
	    x->maxHigh=GMAX(x->high,GMAX(x->left->maxHigh,x->right->maxHigh));
	    x=x->parent;
	  }
  }

  //  FUNCTION:  DeleteFixUp
  //    INPUTS:  x is the child of the spliced
  //            out node in DeleteNode.
  //    EFFECT:  Performs rotations and changes colors to restore red-black
  //             properties after a node is deleted
  void DeleteFixUp(GIntervalTreeNode* x) {
	  GIntervalTreeNode * w;
	  GIntervalTreeNode * rootLeft = root->left;

	  while( (!x->red) && (rootLeft != x)) {
	    if (x == x->parent->left) {
	      w=x->parent->right;
	      if (w->red) {
		w->red=0;
		x->parent->red=1;
		LeftRotate(x->parent);
		w=x->parent->right;
	      }
	      if ( (!w->right->red) && (!w->left->red) ) {
		w->red=1;
		x=x->parent;
	      } else {
		if (!w->right->red) {
		  w->left->red=0;
		  w->red=1;
		  RightRotate(w);
		  w=x->parent->right;
		}
		w->red=x->parent->red;
		x->parent->red=0;
		w->right->red=0;
		LeftRotate(x->parent);
		x=rootLeft; // this is to exit while loop
	      }
	    } else { // the code below is has left and right switched from above
	      w=x->parent->left;
	      if (w->red) {
		w->red=0;
		x->parent->red=1;
		RightRotate(x->parent);
		w=x->parent->left;
	      }
	      if ( (!w->right->red) && (!w->left->red) ) {
		w->red=1;
		x=x->parent;
	      } else {
		if (!w->left->red) {
		  w->right->red=0;
		  w->red=1;
		  LeftRotate(w);
		  w=x->parent->left;
		}
		w->red=x->parent->red;
		x->parent->red=0;
		w->left->red=0;
		RightRotate(x->parent);
		x=rootLeft; // this is to exit while loop
	      }
	    }
	  }
	  x->red=0;
	}

  // Make sure the maxHigh fields for everything makes sense.
  void CheckMaxHighFields(GIntervalTreeNode * x) const {
    if (x != nil) {
      CheckMaxHighFields(x->left);
      if(!(CheckMaxHighFieldsHelper(x,x->maxHigh,0) > 0)) {
        GEXIT("Error found in CheckMaxHighFields.\n");
      }
      CheckMaxHighFields(x->right);
    }
  }

  int CheckMaxHighFieldsHelper(GIntervalTreeNode * y,
  				    const int currentHigh,
  					int match) const {
  	if (y != nil) {
  		match = CheckMaxHighFieldsHelper(y->left,currentHigh,match) ?
  				1 : match;
  		GVERIFY(y->high <= currentHigh);
  		if (y->high == currentHigh)
  			match = 1;
  		match = CheckMaxHighFieldsHelper(y->right,currentHigh,match) ?
  				1 : match;
  	}
  	return match;
  }public:
  GIntervalTree():recursionNodeStackSize(128),
	 recursionNodeStack(NULL), currentParent(0), recursionNodeStackTop(1),
	 root(new GIntervalTreeNode), nil(new GIntervalTreeNode) {
	//nil = new IntervalTreeNode;
	nil->left = nil->right = nil->parent = nil;
	nil->red = 0;
	nil->key = nil->high = nil->maxHigh = INT_MIN;
	nil->storedInterval = NULL;

	//root = new IntervalTreeNode;
	root->parent = root->left = root->right = nil;
	root->key = root->high = root->maxHigh = INT_MAX;
	root->red=0;
	root->storedInterval = NULL;

	/* the following are used for the Enumerate function */
	//recursionNodeStackSize = 128;
	GMALLOC(recursionNodeStack, recursionNodeStackSize*sizeof(G_ITRecursionNode));
	//recursionNodeStackTop = 1;
	recursionNodeStack[0].start_node = NULL;
}

~GIntervalTree() {
		GIntervalTreeNode * x = root->left;
		GVec<GIntervalTreeNode *> stuffToFree;
		if (x != nil) {
			if (x->left != nil) {
				stuffToFree.Push(x->left);
			}
			if (x->right != nil) {
				stuffToFree.Push(x->right);
			}
			// delete x->storedInterval;
			delete x;
			while( stuffToFree.Count()>0 ) {
				x = stuffToFree.Pop();
				if (x->left != nil) {
					stuffToFree.Push(x->left);
				}
				if (x->right != nil) {
					stuffToFree.Push(x->right);
				}
				// delete x->storedInterval;
				delete x;
			}
		}
		delete nil;
		delete root;
		GFREE(recursionNodeStack);
  }
  void Print() const { TreePrintHelper(root->left); }


  //  FUNCTION:  DeleteNode
  //
  //    INPUTS:  tree is the tree to delete node z from
  //    OUTPUT:  returns the Interval stored at deleted node
  //    EFFECT:  Deletes z from tree and but don't call destructor
  //             Then calls FixUpMaxHigh to fix maxHigh fields then calls
  //             DeleteFixUp to restore red-black properties
  GSeg* DeleteNode(GIntervalTreeNode* z) {
	  GIntervalTreeNode* y;
	  GIntervalTreeNode* x;
	  GSeg* returnValue = z->storedInterval;

	  y= ((z->left == nil) || (z->right == nil)) ? z : GetSuccessorOf(z);
	  x= (y->left == nil) ? y->right : y->left;
	  if (root == (x->parent = y->parent)) { // assignment of y->p to x->p is intentional
	    root->left=x;
	  } else {
	    if (y == y->parent->left) {
	      y->parent->left=x;
	    } else {
	      y->parent->right=x;
	    }
	  }
	  if (y != z) { // y should not be nil in this case
	#ifdef DEBUG_ASSERT
	    Assert( (y!=nil),"y is nil in DeleteNode \n");
	#endif
	    // y is the node to splice out and x is its child
	    y->maxHigh = INT_MIN;
	    y->left=z->left;
	    y->right=z->right;
	    y->parent=z->parent;
	    z->left->parent=z->right->parent=y;
	    if (z == z->parent->left) {
	      z->parent->left=y;
	    } else {
	      z->parent->right=y;
	    }
	    FixUpMaxHigh(x->parent);
	    if (!(y->red)) {
	      y->red = z->red;
	      DeleteFixUp(x);
	    } else
	      y->red = z->red;
	    delete z;
	  } else {
	    FixUpMaxHigh(x->parent);
	    if (!(y->red)) DeleteFixUp(x);
	    delete y;
	  }
	  return returnValue;
	}

  //  Before calling InsertNode  the node x should have its key set
  //  FUNCTION:  InsertNode
  //  INPUT:  newInterval is the interval to insert
  //  OUTPUT:  This function returns a pointer to the newly inserted node
  //           which is guaranteed to be valid until this node is deleted.
  //           What this means is if another data structure stores this
  //           pointer then the tree does not need to be searched when this
  //           is to be deleted.
  //  EFFECTS:  Creates a node node which contains the appropriate key and
  //           info pointers and inserts it into the tree.
  GIntervalTreeNode * Insert(GSeg* newInterval) {
	  GIntervalTreeNode* y;
	  GIntervalTreeNode* newNode;
	  GIntervalTreeNode* x = new GIntervalTreeNode(newInterval);
	  TreeInsertHelp(x);
	  FixUpMaxHigh(x->parent);
	  newNode = x;
	  x->red=1;
	  while(x->parent->red) { // use sentinel instead of checking for root
	    if (x->parent == x->parent->parent->left) {
	      y=x->parent->parent->right;
	      if (y->red) {
		x->parent->red=0;
		y->red=0;
		x->parent->parent->red=1;
		x=x->parent->parent;
	      } else {
		if (x == x->parent->right) {
		  x=x->parent;
		  LeftRotate(x);
		}
		x->parent->red=0;
		x->parent->parent->red=1;
		RightRotate(x->parent->parent);
	      }
	    } else { // case for x->parent == x->parent->parent->right
	             // this part is just like the section above with
	             // left and right interchanged
	      y=x->parent->parent->left;
	      if (y->red) {
		x->parent->red=0;
		y->red=0;
		x->parent->parent->red=1;
		x=x->parent->parent;
	      } else {
		if (x == x->parent->left) {
		  x=x->parent;
		  RightRotate(x);
		}
		x->parent->red=0;
		x->parent->parent->red=1;
		LeftRotate(x->parent->parent);
	      }
	    }
	  }
	  root->left->red=0;
	  return(newNode);
	}

  //  FUNCTION:  GetSuccessorOf
  //  INPUTS:  x is the node we want the successor of
  //  OUTPUT:  This function returns the successor of x or NULL if no
  //             successor exists.
  GIntervalTreeNode * GetPredecessorOf(GIntervalTreeNode* x) const {
	  GIntervalTreeNode* y;
	  if (nil != (y = x->right)) { // assignment to y is intentional
	    while(y->left != nil) { // returns the minium of the right subtree of x
	      y=y->left;
	    }
	    return(y);
	  } else {
	    y=x->parent;
	    while(x == y->right) { // sentinel used instead of checking for nil
	      x=y;
	      y=y->parent;
	    }
	    if (y == root) return(nil);
	    return(y);
	  }

  }
  //  FUNCTION:  GetPredecessorOf
  //    INPUTS:  x is the node to get predecessor of
  //    OUTPUT:  This function returns the predecessor of x or NULL if no
  //             predecessor exists.
  GIntervalTreeNode * GetSuccessorOf(GIntervalTreeNode* x) const {
	  GIntervalTreeNode* y;
	  if (nil != (y = x->left)) { // assignment to y is intentional
	    while(y->right != nil) { // returns the maximum of the left subtree of x
	      y=y->right;
	    }
	    return(y);
	  } else {
	    y=x->parent;
	    while(x == y->left) {
	      if (y == root) return(nil);
	      x=y;
	      y=y->parent;
	    }
	    return(y);
	  }

  }

  //  FUNCTION:  Enumerate
  //    INPUTS:  tree is the tree to look for intervals overlapping the
  //             closed interval [low,high]
  //    OUTPUT:  stack containing pointers to the nodes overlapping
  //             [low,high]
  //    EFFECT:  Returns a stack containing pointers to nodes containing
  //             intervals which overlap [low,high] in O(max(N,k*log(N)))
  //             where N is the number of intervals in the tree and k is
  //             the number of overlapping intervals

  //    Note:    This basic idea for this function comes from the
  //              _Introduction_To_Algorithms_ book by Cormen et al, but
  //             modifications were made to return all overlapping intervals
  //             instead of just the first overlapping interval as in the
  //             book.  The natural way to do this would require recursive
  //             calls of a basic search function.  I translated the
  //             recursive version into an iterative version with a stack
  //            as described below.

  //  The basic idea for the function below is to take the IntervalSearch
  //  function from the book and modify to find all overlapping intervals
  //  instead of just one.  This means that any time we take the left
  //  branch down the tree we must also check the right branch if and only if
  //  we find an overlapping interval in that left branch.  Note this is a
  //  recursive condition because if we go left at the root then go left
  //  again at the first left child and find an overlap in the left subtree
  //  of the left child of root we must recursively check the right subtree
  //  of the left child of root as well as the right child of root.
  GVec<GSeg*> * Enumerate(int low, int high) {
		GVec<GSeg*> * enumResultStack;
		GIntervalTreeNode* x=root->left;
		int stuffToDo = (x != nil);

		// Possible speed up: add min field to prune right searches

	#ifdef DEBUG_ASSERT
		Assert((recursionNodeStackTop == 1),
				"recursionStack not empty when entering IntervalTree::Enumerate");
	#endif
		currentParent = 0;
		enumResultStack = new GVec<GSeg*>(4);

		while(stuffToDo) {
			//if (Overlap(low,high,x->key,x->high) ) {
			if (low<=x->high && x->key<=high) {
				enumResultStack->Push(x->storedInterval);
				recursionNodeStack[currentParent].tryRightBranch=1;
			}
			if(x->left->maxHigh >= low) { // implies x != nil
				if ( recursionNodeStackTop == recursionNodeStackSize ) {
					recursionNodeStackSize *= 2;
					recursionNodeStack = (G_ITRecursionNode *)
				  realloc(recursionNodeStack,
						  recursionNodeStackSize * sizeof(G_ITRecursionNode));
					if (recursionNodeStack == NULL)
						GEXIT("realloc failed in IntervalTree::Enumerate\n");
				}
				recursionNodeStack[recursionNodeStackTop].start_node = x;
				recursionNodeStack[recursionNodeStackTop].tryRightBranch = 0;
				recursionNodeStack[recursionNodeStackTop].parentIndex = currentParent;
				currentParent = recursionNodeStackTop++;
				x = x->left;
			} else {
				x = x->right;
			}
			stuffToDo = (x != nil);
			while( (!stuffToDo) && (recursionNodeStackTop > 1) ) {
				if(recursionNodeStack[--recursionNodeStackTop].tryRightBranch) {
					x=recursionNodeStack[recursionNodeStackTop].start_node->right;
					currentParent=recursionNodeStack[recursionNodeStackTop].parentIndex;
					recursionNodeStack[currentParent].tryRightBranch=1;
					stuffToDo = ( x != nil);
				}
			}
		}
	#ifdef DEBUG_ASSERT
		Assert((recursionNodeStackTop == 1),
				"recursionStack not empty when exiting IntervalTree::Enumerate");
	#endif
		return(enumResultStack);
	}
};


#endif
