/*
 * usgread.h
 *
 *  Created on: Sep 4, 2023
 *      Author: gpertea
 */

#ifndef USGREAD_H_
#define USGREAD_H_
#include "GVec.hh"
#include "rlink.h"

extern bool useUSG;


typedef int SG_getChrIDFunc(const char* chr);

enum SGNodeType {
    TSTART,
    TEND,
    JSTART,
    JEND
};

struct SGNode; //forward declaration

struct SGJxLink {
    uint pos;
    char strand;
    uint count;
    SGNode* jx; // Link to the node at position pos in this bundle
};

// Structure for representing a node in the splicing graph
struct SGNode {
    SGNodeType type=TSTART; // Type of the node
    uint position=0;    // Position of the node
    double cov_perc=0;
    uint pos_cov=0;
    double avg_cov_next=0;
    double cov_diff=0; // tstart/tend coverage differential
    int bidx=-1; // index in bundle.nodes array
    // Vector to hold junction links (only populated for JSTART and JEND)
    GVec<SGJxLink> jxLinks;
    //comparison operators required for sorting
    bool operator==(SGNode& d){
          return (position==d.position && type==d.type);
      }
    bool operator<(SGNode& d){
         return (position==d.position)? (type<d.type) : (position<d.position);
      }
    const char* catTokens(GDynArray<char*>& tokens) {
    	//in-place restoration of a strsplit result
    	char* r=NULL;
    	if (tokens.Count()==0) return "";
    	r=tokens[0];
    	for (uint i=1;i<tokens.Count();i++) {
    		char* p = tokens[i];
    		p--; p[0]='\t';
    	}
    	return r;
    }
    SGNode() {};
    SGNode(SGNodeType _type, uint _pos):type(_type), position(_pos) {};
    SGNode(GDynArray<char*>& tokens) {
		 if (tokens.Count()<6) {
			 const char* catline=catTokens(tokens);
			 GError("Error: SGNode line must have at least 6 tokens\n%s\n", catline);
		 }
		 this->type=(strcmp(tokens[0], "tstart") == 0) ? TSTART :
				  (strcmp(tokens[0], "tend") == 0) ? TEND :
				  (strcmp(tokens[0], "jstart") == 0) ? JSTART : JEND;
		 this->position = (uint)atoi(tokens[1]);
		 this->cov_perc = atof(tokens[2]);
		 this->pos_cov = (uint)atoi(tokens[3]);
		 this->avg_cov_next = atof(tokens[4]);
		 if (this->type<JSTART) { //tstart/tend
			this->cov_diff=atof(tokens[5]);
		 } else { // jstart/jend - parse list
			for (uint i=5;i<tokens.Count();i++) {
				SGJxLink newLink;
				GDynArray<char*> lnktok(3);
				int numtk=strsplit(tokens[i], lnktok, ':');
				if (numtk!=3) {
					 const char* catline=catTokens(tokens);
					 GError("Error parsing SGNode j-line links:\n%s\n", catline);
				}
				newLink.pos = atoi(lnktok[0]);
				newLink.strand = lnktok[1][0];
				newLink.count = atoi(lnktok[2]);
				this->jxLinks.Add(newLink);
			}
		 }
    }
};

struct SGBundle: public GSeg {
    int id=-1;   // bundle number
    //int chrID=-1; // id of chromosome
    char chr[192]; //contig/chr name
    double avg_cov=0;  // average bundle coverage
    GList<SGNode> nodes; // nodes in this bundle
	GVec<int> tstarts; // indices of tstart nodes
    SGBundle():nodes(true, true) {} // nodes are owned and sorted
    void clear() { nodes.Clear();  tstarts.Clear(); start=end=0; id=-1; chr[0]=0; avg_cov=0; }
    void setup(uint _start, uint _end, int _id, const char* _chr, double cov=0) {
    	start=_start; end=_end; id=_id;avg_cov=cov;
    	int slen=strlen(_chr);
    	if (slen>191) GMessage("Warning: reference name for usg bundle %d is too long (%s). Truncated.\n",
    			_id, _chr);
    	slen=GMIN(slen,191);
    	strncpy(chr, _chr, slen);
    	chr[slen]=0;
    }
};

class SGReader {
     FILE* fh=NULL;
	 char* linebuf=NULL;
	 int linebufcap=5000;
	 bool pushedBack=false;
	 int linelen=0;
	 char* nextLine() { //also sets linelen
		 if (pushedBack) {
			 pushedBack=false;
			 return linebuf;
		 }
		 char* r=fgetline(linebuf, linebufcap, fh, NULL, &linelen);
		 return r;
	 }
	 //SG_getChrIDFunc* getChrId;
  public:
	SGReader(FILE* f, /* SG_getChrIDFunc* _getChrId,*/ uint buflen=0):fh(f) {
	  //getChrId(_getChrId) {
	  if (buflen>0) linebufcap=buflen;
	  GMALLOC(linebuf, linebufcap);
	}
   ~SGReader() { GFREE(linebuf); fclose(fh); }
	bool next(SGBundle& bundle) {
        char* line=nextLine();
        bundle.clear();
        if (!line) return false;
        if (!startsWith(line, "bundle\t"))
        	GError("Error at SGReader::next(), bundle start expected\n%s\n",
        			line);
       GDynArray<char*> tokens;
       int numtok=strsplit(line, tokens, '\t');
       if (numtok!=6) GError("Error parsing bundle line\n%s\n", line);
       //int chrId=getChrId(tokens[1]);
       //if (chrId<0) GError("Error parsing chr Id# from line\n%s\n", line);
       int id=atoi(tokens[2]);
       uint bstart=(uint)atoi(tokens[3]);
       uint bend=(uint)atoi(tokens[4]);
       double bcov=atof(tokens[5]);
       bundle.setup(bstart, bend, id, tokens[1], bcov);
       //now read the bundle content:
       while ((line=nextLine())) {
    	   if (startsWith(line, "bundle\t")) {
    		   pushedBack=true;
    		   break;
    	   }
    	   numtok=strsplit(line, tokens, '\t');
    	   SGNode* node=new SGNode(tokens);
    	   int bidx=bundle.nodes.Add(node);
		   if (node->type==TSTART) bundle.tstarts.Add(bidx);
    	   node->bidx=bidx; // keep track of index in bundle.nodes array
       }
       //2nd pass to update j-node jxLinks pointers
       for (int i = 0; i < bundle.nodes.Count(); ++i) {
               SGNode* node = bundle.nodes[i];
               if (node->type == JSTART || node->type == JEND) {
                   for (int j = 0; j < node->jxLinks.Count(); ++j) {
                       SGJxLink& link = node->jxLinks[j];
                       SGNode findnode(node->type==JSTART ? JEND : JSTART, link.pos);
                       int fidx=-1;
                       if (bundle.nodes.Found(&findnode, fidx)) {
                    	      link.jx=bundle.nodes[fidx];
                    	      //TODO: should we check around for multiple entries of the same type with the same coordinate?
                          } else {
                    		   GMessage("Warning: no matching jx found for link %d %c of jx %d in bundle %d\n",
                    				   link.pos, link.strand, node->position, bundle.id);
                    	  }
                   }
               }
           }
       return true;
	}
};


struct UCov: public GSeg {
	float cov;
	UCov *next;
	UCov(int cstart=0, int cend=0, float _cov=0):GSeg(cstart, cend),
			cov(_cov),next(NULL) {}
};



struct UCNode { // a universal graph "coverage node" (might actually contain more than one node)
	int prevSGnode; // index of previous SGnode in USG bundle
	UCov *covintv[3]; // covered regions in this UCNode on 0:-, 1:., and 2:+
	GPVec<UCNode> parentlink[3]; // links to linked nodes in the same bundle that came before the node -> needed in order to find connected components
	GPVec<UCNode> childlink[3]; // links to coverage nodes in the same bundle
	GVec<float> childcov[3]; // coverage linking this UCNode to child UCNode

	// Constructor for UCNode
	UCNode(int u = 0): prevSGnode(u) {
		// Initialize childlink pointers to nullptr (default is null pointer)
		for (int s = 0; s < 3; s++) {
			covintv[s]=NULL;
			parentlink[s] = GPVec<UCNode>(); // TODO: do I need this?
			childlink[s] = GPVec<UCNode>(); // TODO: do I need this?
			childcov[s]=GVec<float>(); // TODO: do I need this?
		}
	}
};

/*
struct UGroup: public GSeg { // universal group -> links several UNode's together that would belong to the same UGraph
	UGroup *ulk; // link to next ugroup
	GPVec<UNode> unode;
	UGroup(int rstart=0, int rend=0, UGroup *_ulk=NULL): GSeg(rstart, rend), ulk(_ulk),unode(true) { }
};*/

struct BundleData;

int build_usg(BundleData* bdata);


#endif /* USGREAD_H_ */
