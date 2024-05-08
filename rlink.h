#ifndef __RLINK_H__
#define __RLINK_H__
#include "GArgs.h"
#include "GStr.h"
#include "gff.h"
#include "GSam.h"
#include "GBitVec.h"
#include "time.h"
#include "tablemaker.h"
#include "GHashMap.hh"
//#include "cds.h"

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
const int MAX_MAXCOMP=200; // is 200 too much, or should I set it up to 150?

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

extern bool mergeMode;
extern bool forceBAM; //for stdin alignment data

extern bool verbose;
extern bool debugMode;

//collect all refguide transcripts for a single genomic sequence
struct GRefData {
  GList<GffObj> rnas; //all transcripts on this genomic seq
  int gseq_id;
  const char* gseq_name;
  GRefData(int gid=-1):rnas(false,true,false),gseq_id(gid),gseq_name(NULL) {
    gseq_id=gid;
    if (gseq_id>=0)
       gseq_name=GffObj::names->gseqs.getName(gseq_id);
  }

  void add(GffReader* gffr, GffObj* t) {
     if (gseq_id>=0) {
        if (gseq_id!=t->gseq_id)
           GError("Error: invalid call to GRefData::add() - different genomic sequence!\n");
     }
     else { //adding first transcript, initialize storage
        gseq_id=t->gseq_id;
        gseq_name=t->getGSeqName();
        if (gffr->gseqtable[gseq_id]==NULL)
            GError("Error: invalid genomic sequence data (%s)!\n",gseq_name);
        rnas.setCapacity(gffr->gseqtable[gseq_id]->fcount);
     }
     rnas.Add(t);
     t->isUsed(true);
     //setLocus(t); //use the GRefLocus::mexons to quickly find an overlap with existing loci, or create a new one
  }

  bool operator==(GRefData& d){
    return gseq_id==d.gseq_id;
  }
  bool operator<(GRefData& d){
    return (gseq_id<d.gseq_id);
  }
};


struct CBundlenode:public GSeg {
	float cov;
	int bid; // bundle node id in bnode -> to easy retrieve it
	CBundlenode *nextnode; // next node in the same bundle
	CBundlenode(int rstart=0, int rend=0, float _cov=0, int _bid=-1, CBundlenode *_nextnode=NULL):GSeg(rstart, rend),
			cov(_cov),bid(_bid),nextnode(_nextnode) {}
};


enum GPFType {
	GPFT_NONE=0,
	GPFT_TSS,
	GPFT_CPAS
}; //on 4 bits: maximum 15 types

struct GPtFeature { //point feature (single coordinate)
 GPFType ftype : 4;
 int ref_id: 26; //index in a reftable[] with reference names, max 67,108,863
 int strand: 2; //-1=-, 0=unstranded, +1=+
 uint coord; //genomic coordinate for this feature
 GPtFeature(GPFType ft=GPFT_NONE, int rid=-1, int _strand=0, uint loc=0):ftype(ft),
		    ref_id(rid), strand(_strand), coord(loc) {}
 bool operator<(const GPtFeature &o) { return coord<o.coord; }
 bool operator==(const GPtFeature &o) { return coord==o.coord; }
    //-- should really match ftype and strand too,
    //   but we don't care, for bundle inclusion
};


struct GRefPtData {
  int ref_id; //same with GPtFeature::ref_id, also in GffObj::names->gseqs
  GList<GPtFeature> pfs; //all point feature on this genomic seq, sorted
  GRefPtData(int gid=-1):ref_id(gid), pfs(true,false,false) { }

  void add(GPtFeature* t) { //adds a fully formed GPtFeature record
     if (ref_id!=t->ref_id || ref_id<0 || t->ref_id<0)
           GError("Error: invalid call to GRefPtData::add() - cannot add feature with ref id %d to ref data id %d!\n",
        		   t->ref_id, ref_id);
     pfs.Add(t);
  }
  //sorting by unique ref_id
  bool operator==(GRefPtData& d){
    return ref_id==d.ref_id;
  }
  bool operator<(GRefPtData& d){
    return (ref_id<d.ref_id);
  }

};


// bundle data structure, holds all data needed for
// infering transcripts from a bundle
enum BundleStatus {
	BUNDLE_STATUS_CLEAR=0, //available for loading/prepping
	BUNDLE_STATUS_LOADING, //being prepared by the main thread (there can be only one)
	BUNDLE_STATUS_READY //ready to be processed, or being processed
};

struct CBundle {
	int len;
	float cov;
	float multi;
	int startnode;  // id of start node in bundle of same strand
	int lastnodeid; // id of last node added to bundle
	CBundle(int _len=0, float _cov=0, float _multi=0, int _start=-1, int _last=-1):
			len(_len),cov(_cov),multi(_multi), startnode(_start),lastnodeid(_last) {}
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
	CGuide(CTransfrag* _trf=NULL, int _g=-1):trf(_trf),g(_g) {}
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

/*struct TwoFloat{
	float start;
	float end;
	TwoFloat(float v1=0,float v2=0):start(v1),end(v2) {}
};*/

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

struct CPrediction:public GSeg {
	int geneno;
	GffObj* t_eq; //equivalent reference transcript (guide)
	//char *id;
	float cov;
	float longcov;
	char strand;
	//float frag; // counted number of fragments associated with prediction
	int tlen;
	bool flag;
	GVec<GSeg> exons;
	GVec<float> exoncov;
	GStr mergename;
	CPrediction(int _geneno=0, GffObj* guide=NULL, int gstart=0, int gend=0, float _cov=0, char _strand='.',
	int _len=0,bool f=true):GSeg(gstart,gend), geneno(_geneno),t_eq(guide),cov(_cov),longcov(0),strand(_strand),
	//CPrediction(int _geneno=0, char* _id=NULL,int gstart=0, int gend=0, float _cov=0, char _strand='.', float _frag=0,
	//		int _len=0,bool f=true):GSeg(gstart,gend), geneno(_geneno),id(_id),cov(_cov),strand(_strand),frag(_frag),
			tlen(_len),flag(f),exons(),exoncov(),mergename() {}
	void init(int _geneno=0, GffObj* guide=NULL, int gstart=0, int gend=0, float _cov=0, char _strand='.',
	          int _len=0) {
		geneno=_geneno;
		t_eq=guide;
		start=gstart;
		end=gend;
		cov=_cov;
		strand=_strand;
		tlen=_len;
		flag=true;
		exons.Clear();
		exoncov.Clear();
		mergename.clear();
	}

	CPrediction(CPrediction& c):GSeg(c.start, c.end), geneno(c.geneno),
//			id(Gstrdup(c.id)), cov(c.cov), strand(c.strand), frag(c.frag), tlen(c.tlen), flag(c.flag),
			t_eq(c.t_eq), cov(c.cov), longcov(c.longcov),strand(c.strand), tlen(c.tlen), flag(c.flag),
	      exons(c.exons),  exoncov(c.exoncov), mergename(c.mergename) {}
	~CPrediction() { //GFREE(id);
		}
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
	float cov;    // this is the actual gene coverage
	float covsum; // this is a sum of transcripts coverages -> this is what we need for FPKM and TPM estimations
	GVec<GSeg> exons;  // all possible exons in gene (those are bnodes in bundle)
	CGene(int gstart=0, int gend=0, char _strand='.',char *gid=NULL, char *gname=NULL):GSeg(gstart,gend),
		strand(_strand), geneID(gid), geneName(gname), exons() { cov=0; covsum=0;}
	// getGeneID() and getGeneName() functions of gffobj return pointers to this attributes in gffobj so I don't need to clean them up here
};

//holding transcript info for --merge mode
struct TAlnInfo {
	GStr name; //transcript name
	int fileidx; //index of transcript file in the TInputFiles.files array
	double cov;
	double fpkm;
	double tpm;
	int g;
	TAlnInfo(const char* rname=NULL, int fidx=0):name(rname), fileidx(fidx),
			cov(-1),fpkm(-1),tpm(-1),g(-1) { }
};

struct CJunction;

struct CReadAln:public GSeg {
	//DEBUG ONLY:
	// GStr name;
	char strand; // 1, 0 (unkown), -1 (reverse)
	short int nh;
	uint len;
	float read_count;       // keeps count for all reads (including paired and unpaired)
	bool unitig:1;			// set if read come from an unitig
	bool longread:1;	    // set if read comes from long read data
	GVec<float> pair_count;   // keeps count for all paired reads
	GVec<int> pair_idx;     // keeps indeces for all pairs in assembly mode, or all reads that were collapsed in merge mode
	GVec<GSeg> segs; //"exons"
	GPVec<CJunction> juncs;
	union {
		TAlnInfo* tinfo;
		bool in_guide;
	};

	CReadAln(char _strand=0, short int _nh=0,
			int rstart=0, int rend=0, TAlnInfo* tif=NULL): GSeg(rstart, rend), //name(rname),
					strand(_strand),nh(_nh), len(0), read_count(0), unitig(false),longread(false),pair_count(),pair_idx(),
					segs(), juncs(false), tinfo(tif) { }
	CReadAln(CReadAln &rd):GSeg(rd.start,rd.end) { // copy contructor
		strand=rd.strand;
		nh=rd.nh;
		len=rd.len;
		read_count=rd.read_count;
		unitig=rd.unitig;
		longread=rd.longread;
		pair_count=rd.pair_count;
		pair_idx=rd.pair_idx;
		juncs=rd.juncs;
		tinfo=rd.tinfo;
	}
	int overlapSegLen(CReadAln* r) {

		if (r->start>end || start>r->end) return 0;

		int i=0;
		int j=0;
		int len=0;
		while(i<segs.Count()) {
			if(segs[i].end<r->segs[j].start) i++;
			else if(r->segs[j].end<segs[i].start) j++;
			else { // there is overlap
				len+=segs[i].overlapLen(r->segs[j].start,r->segs[j].end);
				if(segs[i].end<r->segs[j].end) i++;
				else j++;
			}
			if(j==r->segs.Count()) break;
		}
		return len;
	}
	~CReadAln() { if(mergeMode) {delete tinfo;} }
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

/*
struct CSegCov:public GSeg {
	bool spliced:1;
	GVec<GStr> rname;
	CSegCov *next; // next interval;
	CSegCov(uint start=0,uint end=0):GSeg(start,end),spliced(false),rname(),next(NULL) {}
};
*/

struct CMaxIntv:public GSeg {
	GVec<CExon> node;
	CMaxIntv *next; // next interval;
	CMaxIntv(uint start=0,uint end=0):GSeg(start,end),node(),next(NULL) {}
	CMaxIntv(GVec<CExon>& _node,uint start,uint end,CMaxIntv *_next=NULL):GSeg(start,end),node(_node),next(_next) {}
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
	float capacity; // sum of all transcripts abundances exiting and through node
	float rate; // conversion rate between in and out transfrags of node
	//float frag; // number of fragments included in node
	GVec<int> child;
	GVec<int> parent;
	GBitVec childpat;
	GBitVec parentpat;
	GVec<int> trf; // transfrags that pass the node
	bool hardstart:1; // verified/strong start
	bool hardend:1;	// verified/strong end
	//CGraphnode(int s=0,int e=0,unsigned int id=MAX_NODE,float nodecov=0,float cap=0,float r=0,float f=0):GSeg(s,e),nodeid(id),cov(nodecov),capacity(cap),rate(r),frag(f),child(),parent(),childpat(),parentpat(),trf(){}
	CGraphnode(int s=0,int e=0,unsigned int id=MAX_NODE,float nodecov=0,float cap=0,float r=0):GSeg(s,e),
			nodeid(id),cov(nodecov),capacity(cap),rate(r),child(),parent(),childpat(),parentpat(),trf(),hardstart(false),hardend(false){}
};

// # 0: strand; 1: start; 2: end; 3: nreads; 4: nreads_good;
struct CJunction:public GSeg {
	char strand; //-1,0,1
	char guide_match; //exact match of a ref intron?
	char consleft; // -1,0,1 -1 is not set up, 0 is non consensus, 1 is consensus
	char consright; // -1,0,1 -1 is not set up, 0 is non consensus, 1 is consensus
	double nreads;
	double nreads_good;
	double leftsupport;
	double rightsupport;
	double nm; // number of reads with a high nm (high mismatch)
	double mm; // number of reads that support a junction with both anchors bigger than longintronanchor
	CJunction(int s=0,int e=0, char _strand=0):GSeg(s,e),
			strand(_strand), guide_match(0), consleft(-1), consright(-1),nreads(0),nreads_good(0),
			leftsupport(0),rightsupport(0),nm(0),mm(0) {}
	bool operator<(CJunction& b) {
		if (start<b.start) return true;
		if (start>b.start) return false;
		if (end<b.end) return true;
		if (end>b.end) return false;
		if (strand>b.strand) return true;
		return false;
	}
	bool operator==(CJunction& b) {
		return (start==b.start && end==b.end && strand==b.strand);
	}
};

struct GReadAlnData {
	GSamRecord* brec;
	char strand; //-1, 0, 1
	int nh;
	int hi;
	GPVec<CJunction> juncs;
	union {
		TAlnInfo* tinfo;
		bool in_guide;
	};
	//GPVec< GVec<RC_ExonOvl> > g_exonovls; //>5bp overlaps with guide exons, for each read "exon"
	GReadAlnData(GSamRecord* bamrec=NULL, char nstrand=0, int num_hits=0,
			int hit_idx=0, TAlnInfo* tif=NULL):brec(bamrec), strand(nstrand),
					nh(num_hits), hi(hit_idx), juncs(true), tinfo(tif) { } //, g_exonovls(true)
	~GReadAlnData() { if(mergeMode) delete tinfo; }
};


/*
struct CTCov { //covered transcript info
	int first_cov_exon;
	int last_cov_exon;
	int numt;
	GffObj* guide;
	bool whole;
	CTCov(GffObj* t, int fex=-1, int lex=0, int ntr=0):first_cov_exon(fex), last_cov_exon(lex),
			   numt(ntr), guide(t), whole(false) {
		whole = (first_cov_exon<0);
	}
	void print(FILE* f) {
		if (whole) { //from get_covered()
			guide->printTranscriptGff(f);
		}
		else { //from get_partial_covered()
			bool partial=true;
			if (last_cov_exon<0) {
				if (guide->exons.Count()==1) partial=false;
				last_cov_exon=first_cov_exon;
			} else {
			 if(last_cov_exon-first_cov_exon+1==guide->exons.Count()) partial=false;
			}
			for(int i=first_cov_exon;i<=last_cov_exon;i++) {
				if(partial) fprintf(f, "%s\tpartial\texon\t%u\t%u\t.\t%c\t.\ttranscript_id \"%s_part%d\";\n",guide->getGSeqName(),
						guide->exons[i]->start,guide->exons[i]->end,guide->strand,guide->getID(), numt);
				else fprintf(f, "%s\tcomplete\texon\t%u\t%u\t.\t%c\t.\ttranscript_id \"%s\";\n",guide->getGSeqName(),
						guide->exons[i]->start,guide->exons[i]->end,guide->strand,guide->getID());
			}
		}
	}
};
*/

// bundle data structure, holds all input data parsed from BAM file
// - r216 regression
struct BundleData {
 BundleStatus status;
 //int64_t bamStart; //start of bundle in BAM file
 int idx; //index in the main bundles array
 int start;
 int end;
 unsigned long numreads; // number of reads in this bundle
 /*
 float wnumreads; // NEW: weighted numreads; a multi-mapped read mapped in 2 places will contribute only 0.5
 double sumreads; // sum of all reads' lengths in bundle
 double sumfrag; // sum of all fragment lengths (this includes the insertion so it is an estimate)
 float num_reads; // number of all reads in bundle that we considered (weighted)
 float num_cov; // how many coverages we added (weighted) to obtain sumcov
 float num_frag; // how many fragments we added to obtain sumfrag
 double num_fragments3;
 double sum_fragments3;
*/
 double num_fragments; //aligned read/pairs
 double frag_len;
 double sum_cov; // sum of all transcripts coverages --> needed to compute TPMs
 char covflags;

 GStr refseq; //reference sequence name
 char* gseq; //actual genomic sequence for the bundle
 GList<CReadAln> readlist;
 GVec<float> bpcov[3];   // this needs to be changed to a more inteligent way of storing the data
 GList<CJunction> junction;
 GPVec<GffObj> keepguides;
 GPVec<GPtFeature> ptfs; //point features for this bundle
 GList<CPrediction> pred;
 RC_BundleData* rc_data;
 BundleData():status(BUNDLE_STATUS_CLEAR), idx(0), start(0), end(0),
		 numreads(0),
		 num_fragments(0), frag_len(0),sum_cov(0),covflags(0),
		 refseq(), gseq(NULL), readlist(false,true), //bpcov(1024),
		 junction(true, true, true),
		 keepguides(false), ptfs(false), pred(false), rc_data(NULL) {
	 for(int i=0;i<3;i++) 	bpcov[i].setCapacity(4096);
 }

 void getReady(int currentstart, int currentend) {
	 //this is only called when the bundle is valid and ready to be processed
	 start=currentstart;
	 end=currentend;
	 //refseq=ref;
	 //tag all these guides
	 for (int i=0;i<this->keepguides.Count();++i) {
		 RC_TData* tdata=(RC_TData*)(keepguides[i]->uptr);
		 tdata->in_bundle=1;
	 }
	 status=BUNDLE_STATUS_READY;
 }

 void rc_init(GffObj* t, GPVec<RC_TData>* rc_tdata,
		 GPVec<RC_Feature>* rc_edata, GPVec<RC_Feature>* rc_idata) {
	  if (rc_data==NULL) {
	  	rc_data = new RC_BundleData(t->start, t->end,
	  			rc_tdata, rc_edata, rc_idata);
	  }
 }
 /* after reference annotation was loaded
 void rc_finalize_refs() {
     if (rc_data==NULL) return;
     //rc_data->setupCov();
	}
	Not needed here, we update the coverage span as each transcript is added
 */
 void keepGuide(GffObj* scaff, GPVec<RC_TData>* rc_tdata=NULL,
		 GPVec<RC_Feature>* rc_edata=NULL, GPVec<RC_Feature>* rc_idata=NULL);

 //bool evalReadAln(GSamRecord& brec, char& strand, int nh); //, int hi);
 bool evalReadAln(GReadAlnData& alndata, char& strand);

 void Clear() {
	keepguides.Clear();
	ptfs.Clear();
	pred.Clear();
	pred.setSorted(false);
	readlist.Clear();
	for(int i=0;i<3;i++) {
		bpcov[i].Clear();
		bpcov[i].setCapacity(1024);
	}
	junction.Clear();
	start=0;
	end=0;
	status=BUNDLE_STATUS_CLEAR;
	numreads=0;
	num_fragments=0;
	frag_len=0;
	sum_cov=0;
	covflags=0;
	delete rc_data;
	GFREE(gseq);
	rc_data=NULL;
 }

 ~BundleData() {
	Clear();
 }
};

void processRead(int currentstart, int currentend, BundleData& bdata,
		 GHash<int>& hashread, GReadAlnData& alndata);
		 //GSamRecord& brec, char strand, int nh, int hi);

void countFragment(BundleData& bdata, GSamRecord& brec, int hi,int nh);

int printResults(BundleData* bundleData, int geneno, GStr& refname);
int printMergeResults(BundleData* bundleData, int geneno, GStr& refname);

int infer_transcripts(BundleData* bundle);

// --- utility functions
void printGff3Header(FILE* f, GArgs& args);

void printTime(FILE* f);


#endif
