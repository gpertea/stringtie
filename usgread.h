/*
 * usgread.h
 *
 *  Created on: Sep 4, 2023
 *      Author: gpertea
 */

#ifndef USGREAD_H_
#define USGREAD_H_
#include "GVec.hh"


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


struct UNode:public GSeg { // a universal graph "node" -> links together several entries in the USG bundle
	SGNode *ustart; // start node in bundle
	SGNode *uend; // end node in bundle
	float cov_sum; // read coverage from sample
	float multi; // proportion of multi-mapped reads
	float neg_prop; // proportion of negative reads assigned to group out of all positives and negatives
	UNode(int rstart=0, int rend=0, SGNode *u_start=NULL, SGNode *u_end=NULL, float _cov_sum=0,
			float _multi=0,float _neg_prop=0): GSeg(rstart, rend), ustart(u_start),uend(u_end),cov_sum(_cov_sum),
			multi(_multi), neg_prop(_neg_prop) { }
};

struct UGroup: public GSeg { // universal group -> links several UNode's together that would belong to the same UGraph
	UGroup *ulk; // link to next ugroup
	GPVec<UNode> unode;
	UGroup(int rstart=0, int rend=0, UGroup *_ulk=NULL): GSeg(rstart, rend), ulk(_ulk),unode(true) { }
};

struct BundleData;

void build_usg(BundleData* bdata, GVec<int> &read2unode);


#endif /* USGREAD_H_ */
