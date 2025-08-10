#ifndef __RLINK_H__
#define __RLINK_H__
#include "GArgs.h"
#include "GStr.h"
#include "GBitVec.h"
#include "time.h"
#include "bundle.h" // includes tablemaker.h
#include "GHashMap.hh"

#define MAX_NODE 1000000
#define KMER 31

#define DROP 0.5
#define ERROR_PERC 0.1
#define DBL_ERROR 0.01

#define CHI_WIN 100
#define CHI_THR 50
#define SMALL_EXON 35 // exons smaller than this have a tendency to be missed by long read data

#define IS_FPKM_FLAG 1
#define IS_TPM_FLAG 2
#define IS_COV_FLAG 4

const double epsilon=0.000001; //-E
const float trthr=1.0;   // transfrag pattern threshold
const float MIN_VAL=-100000.0;
const int MAX_DIST=200; // is 200 too much, or should I set it up to 150?

//const uint largeintron=20000; // don't trust introns longer than this unless there is higher evidence; less than 10% of all human annotated introns are longer than this
//const uint longintron=70000; // don't trust introns longer than this unless there is higher evidence; about 98% of all human annotated introns are shorter than this
const uint longintron=100000; // don't trust introns longer than this unless there is higher evidence; about 99% of all human annotated introns are shorter than this
const uint longintronanchor=25; // I need a higher anchor for long introns -> use a smaller value here? i.e. 20?
//const uint verylongintron=100000; // don't trust introns longer than this unless there is higher evidence; about 99% of all human annotated introns are shorter than this
//const uint verylongintronanchor=35; // I need a higher anchor for very long introns

const float mismatchfrac=0.02;
const float lowcov=1.5;
const float lowisofrac=0.02;

const int max_trf_number=40000; // maximum number of transfrag accepted so that the memory doesn't blow up

extern bool forceBAM; //for stdin alignment data

extern bool verbose;
extern bool debugMode;

struct CBundlenode:public GSeg {
	float cov;
	int bid; // bundle node id in bnode -> to easy retrieve it
	CBundlenode *nextnode; // next node in the same bundle
	CBundlenode(int rstart=0, int rend=0, float _cov=0, int _bid=-1, CBundlenode *_nextnode=NULL):GSeg(rstart, rend),
			cov(_cov),bid(_bid),nextnode(_nextnode) {}
};



struct CPath {
	int node;
	int contnode;
	float abundance;
	CPath(int n1=0,int n2=0,float abund=0):node(n1),contnode(n2),abundance(abund){}
};

struct CTransfrag {
	GVec<int> nodes;
	GBitVec pattern;
	float abundance;
	float srabund; // keeps abundance associated to srfrag
	GVec<CPath> path; // stores all the possible paths that leave from a node to reach next node in a transfrag, and distributes the abundance of the transfrag between all possible continuations
	float usepath;
	int weak; // number of weak links
	bool real:1;
	bool longread:1; // there is at least a longread supporting transfrag
	bool shortread:1; // there is at least one short read supporting transfrag
	int guide;
	uint longstart; // for long reads: min start of all longreads sharing transfrag
	uint longend; // for long reads: max end of all longreads sharing transfrag
	CTransfrag(GVec<int>& _nodes,GBitVec& bit, float abund=0, bool treal=false, int tguide=0,float sr=0):nodes(_nodes),pattern(bit),abundance(abund),srabund(sr),path(),usepath(-1),weak(-1),real(treal),longread(false),shortread(false),guide(tguide),longstart(false),longend(false) {}
	CTransfrag(float abund=0, bool treal=false,int tguide=0):nodes(),pattern(),abundance(abund),srabund(0),path(),usepath(-1),weak(-1),real(treal),longread(false),shortread(false),guide(tguide),longstart(false),longend(false) {
	}
};

struct CMTransfrag { // this is the super-class for transfrag -> to use in case of merging transcripts
	CTransfrag *transfrag;
	GVec<int> read; // all reads' indeces that are connected to this transfrag
	int nf;
	int nl;
	uint len;
	CMTransfrag(CTransfrag *t=NULL):transfrag(t),read(),nf(0),nl(0),len(0) {}
};

struct CGuide {
	CTransfrag *trf;
	//GffObj* t; // this is what I was using before but I need to store the guide index in guides instead
	//CGuide(CTransfrag* _trf=NULL, GffObj* _t=NULL):trf(_trf),t(_t) {}
	int g; // stores guide index in guides instead of the actual pointer
	bool nascent:1;
	CGuide(CTransfrag* _trf=NULL, int _g=-1,bool n=false):trf(_trf),g(_g),nascent(n) {}
};

struct CPartGuide {
	int idx;
	int olen; // overlap with node
	int allolen; // overlap with all nodes
	int glen; // guide length
	bool terminal_in;
	bool terminal_out;
	float gcount;
	float cov; // assigned overall coverage so far
	float ncov; // new coverage in node
	CPartGuide(int _idx=0, int _olen=0, int _aolen=0, int _glen=0):idx(_idx),olen(_olen),allolen(_aolen),glen(_glen),
			terminal_in(false),terminal_out(false),gcount(0),cov(0),ncov(0) {}
};

struct CTrGuidePat { // remember abundances based on node guide pattern
	GBitVec pat;
	float abund;
	bool terminal;
	GVec<int> g;
	CTrGuidePat():pat(),abund(0),terminal(false),g() {} // default constructor
	CTrGuidePat(GBitVec p, float a, bool _t): pat(p),abund(a),terminal(_t),g() {}
};

struct CNodeGuide {
	GVec<CPartGuide> guide;
	GVec<CTrGuidePat> trcount;
	float sumtrcount; // sum of all used in guides transfrag abundances -> this is needed in order to see if a guide would explain all of them
	CNodeGuide():guide(),trcount(),sumtrcount(0) {}
};

struct CGroup:public GSeg {
	int color;
	int grid;
	float cov_sum;
	float multi;
	float neg_prop; // proportion of negative reads assigned to group out of all positives and negatives
	CGroup *next_gr;
	CGroup(int rstart=0, int rend=0, int _color=-1, int _grid=0, float _cov_sum=0,
			float _multi=0,float _neg_prop=0, CGroup *_next_gr=NULL): GSeg(rstart, rend), color(_color), grid(_grid),cov_sum(_cov_sum),
			multi(_multi), neg_prop(_neg_prop),next_gr(_next_gr) { }
};

struct CMerge {
	GStr name;
	GVec<int> fidx; // file indices for the transcripts in the merge
	CMerge(const char* rname=NULL):name(rname),fidx() {}
};


struct CExon{
	int predno;
	int exonno;
	float exoncov;
	CExon(int p=0,int e=0,float c=0):predno(p),exonno(e),exoncov(c) {}
};


struct CPred{
	int predno;
	float cov;
	CPred(int p=0,float c=0):predno(p),cov(c) {}
};

struct CLongTrf{
	int t;
	float cov;
	GVec<int> group;
	CLongTrf(int tno=0,float c=0):t(tno),cov(c),group() {}
};

struct CMPrediction {
	CPrediction *p;
	GVec<int> nodes;
	GBitVec pat; // pattern of nodes and introns in prediction
	GBitVec b; // not retained introns
	CMPrediction(CPrediction* _p=NULL): p(_p),nodes(),pat(),b() {}
	CMPrediction(CPrediction* _p,GVec<int>& _nodes,GBitVec& _pat, GBitVec& _b): p(_p),nodes(_nodes),pat(_pat),b(_b) {}
};

struct CNodeCapacity {
	int id;
	bool left;
	float perc;
	CNodeCapacity(int nid=0,bool leftnode=false,float p=0): id(nid),left(leftnode),perc(p) {}
};

// this class keeps the gene predictions (linked bundle nodes initially)
struct CGene:public GSeg { // I don't necessarily need to make this a GSeg since I can get the start&end from the exons
	char strand;
	char* geneID;
	char* geneName;
        double cov;    // this is the actual gene coverage
        double covsum; // this is a sum of transcripts coverages -> this is what we need for FPKM and TPM estimations
	GVec<GSeg> exons;  // all possible exons in gene (those are bnodes in bundle)
        CGene(int gstart=0, int gend=0, char _strand='.',char *gid=NULL, char *gname=NULL):GSeg(gstart,gend),
                strand(_strand), geneID(gid), geneName(gname), exons() { cov=0; covsum=0;}
	// getGeneID() and getGeneName() functions of gffobj return pointers to this attributes in gffobj so I don't need to clean them up here
};

struct CGraphinfo {
	int ngraph;
	int nodeno;
	CGraphinfo(int ng=-1,int nnode=-1):ngraph(ng),nodeno(nnode){}
};

struct CGJunc {
	int leftnode;
	int rightnode;
	double cov; // ngood
	double goodcov; // ngood_reads
	CGJunc(int n1=0,int n2=0,double _cov=0,double _goodcov=0):leftnode(n1),rightnode(n2),cov(_cov),goodcov(_goodcov){}
};


struct CGNode {
	int id;    // initial id in graphno
	bool last; // if this is last node (to be linked to sink later)
	bool keep; // if I keep it in the final count (true by default)
	bool merge; // if this node needs to be merged to its adjacent node
	bool future;
	CGNode(int _id=0,bool _last=false,bool _keep=true, bool _merge=false, bool _future=false):id(_id),last(_last),keep(_keep),merge(_merge),future(_future){}
};

struct CTreePat {
	int nodeno;
	int childno;
	CTransfrag *tr;
	CTreePat **nextpat;
	CTreePat(int n=0,int cno=0):nodeno(n),childno(cno),tr(NULL),nextpat(NULL){
		if(cno) {
			GCALLOC(nextpat,cno*sizeof(CTreePat *));
			for(int i=0;i<cno;i++) nextpat[i]=NULL;
		}
	}
	void setchilds(int cno) {
		if(cno && !nextpat) {
			GCALLOC(nextpat,cno*sizeof(CTreePat *));
			for(int i=0;i<cno;i++) nextpat[i]=NULL;
		}
		childno=cno;
	}
	void settree(int i, CTreePat *t) {
		if(i<childno) nextpat[i]=t;
	}
	CTreePat *settree(int nextpos,int n,int cno) {
		if(nextpos<childno) {
			if(!nextpat[nextpos]) nextpat[nextpos]=new CTreePat(n,cno);
			return(nextpat[nextpos]);
		}
		return(NULL);
	}
};

struct CTrimPoint { // this can work as a guide keeper too, where pos is the guideidx, abundance is the flow, and start is the included status
	uint pos;
	float abundance;
	bool start:1;
	CTrimPoint(uint _pos=0,float abund=0.0,bool _start=true):pos(_pos),abundance(abund),start(_start) {}
};

struct CInterval {
	uint pos; // interval start position
	float val; // interval value or interval last position depending on use
	CInterval *next; // next interval;
	CInterval(uint _pos=0,float _val=0,CInterval *_next=NULL):pos(_pos),val(_val),next(_next) {}
};


struct CMaxIntv:public GSeg {
	GVec<CExon> node;
	float cov;
	CMaxIntv *next; // next interval;
	CMaxIntv(uint start=0,uint end=0):GSeg(start,end),node(),cov(0.0),next(NULL) {}
	CMaxIntv(GVec<CExon>& _node,uint start,uint end,float _cov=0.0,CMaxIntv *_next=NULL):GSeg(start,end),node(_node),cov(_cov),next(_next) {}
};

struct CNasc{
	CPrediction *pred;
	int exonno;
	float exoncov;
	CNasc(CPrediction *p=NULL,int e=0,float c=0):pred(p),exonno(e),exoncov(c) {}
};

struct CNascIntv:public GSeg {
	GVec<CNasc> node;
	float cov;
	CNascIntv *next; // next interval;
	CNascIntv(uint start=0,uint end=0):GSeg(start,end),node(),cov(0.0),next(NULL) {}
	CNascIntv(GVec<CNasc>& _node,uint start,uint end,float _cov=0.0,CNascIntv *_next=NULL):GSeg(start,end),node(_node),cov(_cov),next(_next) {}
};

struct GInterval {
	uint start;
	uint end;
	GInterval *next;
	GInterval(uint _start, uint _end,GInterval *_next=NULL):start(_start),end(_end),next(_next) {}
};

struct CTrInfo {
	int trno;
	float abundance;
	float penalty;
	CTrInfo(int tr=-1,float _abund=0.0, float _pen=0.0):trno(tr),abundance(_abund),penalty(_pen) {}
};

struct CNetEdge {
	int link;
	float rate;
	bool fake;
	CNetEdge(int lnk=0.0,float r=0.0, bool f=false):link(lnk),rate(r),fake(f){}
};

struct CComponent {
	float size;
	GVec<int> *set;
	CComponent(float _size=0.0,GVec<int> *_set=NULL):size(_size),set(_set) {}
	~CComponent() { if(set) delete set;}
};

struct GEdge { // guide edge
	// if val < endval then this is start; otherwise it is end
	uint val;  // value of the boundary
	uint endval; // value of the other exon boundary shared with val
	int strand;
	bool operator<(const GEdge& o) const {
		return(val<o.val || (val==o.val && strand<o.strand));
	}
	bool operator==(const GEdge& o) const {
		return(val==o.val && strand==o.strand);
	}
	GEdge(uint _val=0,uint _endval=0,int _strand=0):val(_val),endval(_endval),strand(_strand) {}
};

struct CGraphnode:public GSeg {
	int nodeid;
	float cov;
	float abundin; //sum of all transfrags entering node
	float abundout; //sum of all transfrags exiting node
	//float frag; // number of fragments included in node
	GVec<int> child;
	GVec<int> parent;
	GBitVec childpat;
	GBitVec parentpat;
	GVec<int> trf; // transfrags that pass the node
	bool hardstart:1; // verified/strong start
	bool hardend:1;	// verified/strong end
	//CGraphnode(int s=0,int e=0,unsigned int id=MAX_NODE,float nodecov=0,float cap=0,float r=0,float f=0):GSeg(s,e),nodeid(id),cov(nodecov),capacity(cap),rate(r),frag(f),child(),parent(),childpat(),parentpat(),trf(){}
	CGraphnode(int s=0,int e=0,unsigned int id=MAX_NODE,float nodecov=0,float in=0,float out=0):GSeg(s,e),
			nodeid(id),cov(nodecov),abundin(in),abundout(out),child(),parent(),childpat(),parentpat(),trf(),hardstart(false),hardend(false){}
};




void countFragment(BundleData& bdata, GSamRecord& brec, int hi,int nh);

int printResults(BundleData* bundleData, int geneno, GStr& refname);
int printMergeResults(BundleData* bundleData, int geneno, GStr& refname);

int infer_transcripts(BundleData* bundle);

// --- utility functions
void printGff3Header(FILE* f, GArgs& args);

void printTime(FILE* f);


#endif
