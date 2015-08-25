#ifndef __RLINK_H__
#define __RLINK_H__
#include "GArgs.h"
#include "GStr.h"
#include "gff.h"
#include "GBam.h"
#include "GBitVec.h"
#include "time.h"
#include "tablemaker.h"

#define MAX_NODE 100000

#define DROP 0.5

#define CHI_WIN 100
#define CHI_THR 50

const double epsilon=0.00000001; //-E
const float trthr=1.0;   // transfrag pattern threshold
const float MIN_VAL=-100000.0;
const int MAX_MAXCOMP=200; // is 200 too much, or should I set it up to 150?

const int longintron=20000; // don't trust introns longer than this unless there is higher evidence; 93.5% of all annotated introns are shorter than this
const int longintronanchor=25; // I need a higher anchor for long introns

const int max_trf_number=40000; // maximum number of transfrag accepted so that the memory doesn't blow up

extern bool singlePass;

//collect all refguide transcripts for a single genomic sequence
struct GRefData {
  GList<GffObj> rnas; //all transcripts on this genomic seq
  int gseq_id;
  const char* gseq_name;
   //GList<GTData> tdata; //transcript data (uptr holder for all rnas loaded here)
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
        gseq_name=t->getGeneName();
        if (gffr->gseqStats[gseq_id]==NULL)
            GError("Error: invalid genomic sequence data (%s)!\n",gseq_name);
        rnas.setCapacity(gffr->gseqStats[gseq_id]->fcount);
     }
     rnas.Add(t);
     t->isUsed(true);
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
	float nread;
	float multi;
	int startnode;  // id of start node in bundle of same strand
	int lastnodeid; // id of last node added to bundle
	CBundle(int _len=0, float _cov=0, float _nread=0,float _multi=0, int _start=-1, int _last=-1):
		len(_len),cov(_cov),nread(_nread),multi(_multi), startnode(_start),lastnodeid(_last) {}
};

struct CTransfrag {
	GVec<int> nodes;
	GBitVec pattern;
	float abundance;
	bool real;
	CTransfrag(GVec<int>& _nodes,GBitVec& bit, float abund=0, bool treal=true):nodes(_nodes),pattern(bit),abundance(abund),real(treal) {}
	CTransfrag(float abund=0, bool treal=true):nodes(),pattern(),abundance(abund),real(treal) {}
};

struct CGuide {
	CTransfrag *trf;
	//char *id;
	//CGuide(CTransfrag* _trf=NULL,char* _id=NULL):trf(_trf),id(_id) {}
	GffObj* t;
	CGuide(CTransfrag* _trf=NULL, GffObj* _t=NULL):trf(_trf),t(_t) {}
};

struct CGroup:public GSeg {
	int color;
	int grid;
	float cov_sum;
	float nread;
	float multi;
	float neg_prop; // proportion of negative reads assigned to group out of all positives and negatives
	CGroup *next_gr;
	CGroup(int rstart=0, int rend=0, int _color=-1, int _grid=0, float _cov_sum=0,float _nread=0,float _multi=0,float _neg_prop=0,
			CGroup *_next_gr=NULL): GSeg(rstart, rend), color(_color), grid(_grid),
			cov_sum(_cov_sum), nread(_nread),multi(_multi), neg_prop(_neg_prop), next_gr(_next_gr) { }
};

struct CPrediction:public GSeg {
	int geneno;
	GffObj* t_eq; //equivalent reference transcript (guide)
	//char *id;
	float cov;
	char strand;
	float frag; // counted number of fragments associated with prediction
	int tlen;
	bool flag;
	GVec<GSeg> exons;
	GVec<float> exoncov;
	CPrediction(int _geneno=0, GffObj* guide=NULL, int gstart=0, int gend=0, float _cov=0, char _strand='.', float _frag=0,
	int _len=0,bool f=true):GSeg(gstart,gend), geneno(_geneno),t_eq(guide),cov(_cov),strand(_strand),frag(_frag),
	//CPrediction(int _geneno=0, char* _id=NULL,int gstart=0, int gend=0, float _cov=0, char _strand='.', float _frag=0,
	//		int _len=0,bool f=true):GSeg(gstart,gend), geneno(_geneno),id(_id),cov(_cov),strand(_strand),frag(_frag),
			tlen(_len),flag(f),exons(),exoncov() {}
	CPrediction(CPrediction& c):GSeg(c.start, c.end), geneno(c.geneno),
//			id(Gstrdup(c.id)), cov(c.cov), strand(c.strand), frag(c.frag), tlen(c.tlen), flag(c.flag),
			t_eq(c.t_eq), cov(c.cov), strand(c.strand), frag(c.frag), tlen(c.tlen), flag(c.flag),
	      exons(c.exons),  exoncov(c.exoncov) {}
	~CPrediction() { //GFREE(id);
		}
};

struct CJunction;

struct CReadAln:public GSeg {
	//DEBUG ONLY:
	// GStr name;
	// 0: strand; 1: NH; 2: pair's no; 3: coords of read; 4: junctions
	char strand; // 1, 0 (unkown), -1 (reverse)
	short int nh;
	int pair_idx; //mate index alignment in CReadAln list
	GVec<GSeg> segs; //"exons"
	GPVec<CJunction> juncs; //junction index in CJunction list
	//DEBUG ONLY: (discard rname when no debugging needed)
	CReadAln(char _strand=0, short int _nh=0,
			int rstart=0, int rend=0 /*,  const char* rname=NULL */): GSeg(rstart, rend), //name(rname),
					strand(_strand), nh(_nh), pair_idx(-1), segs(), juncs(false) { }
};

struct CGraphinfo {
	int ngraph;
	int nodeno;
	CGraphinfo(int ng=-1,int nnode=-1):ngraph(ng),nodeno(nnode){}
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
	bool start;
	CTrimPoint(uint _pos=0,float abund=0.0,bool _start=true):pos(_pos),abundance(abund),start(_start) {}
};

struct CInterval {
	uint pos; // interval start position
	float val; // interval value
	CInterval *next; // next interval;
	CInterval(uint _pos=0,float _val=0,CInterval *_next=NULL):pos(_pos),val(_val),next(_next) {}
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

struct CGraphnode:public GSeg {
	int nodeid;
	float cov;
	float capacity; // sum of all transcripts abundances exiting and through node
	float rate; // conversion rate between in and out transfrags of node
	float frag; // number of fragments included in node
	GVec<int> child;
	GVec<int> parent;
	GBitVec childpat;
	GBitVec parentpat;
	GVec<int> trf; // transfrags that pass the node
	CGraphnode(int s=0,int e=0,int id=MAX_NODE,float nodecov=0,float cap=0,float r=0,float f=0):GSeg(s,e),nodeid(id),
			cov(nodecov),capacity(cap),rate(r),frag(f),child(),parent(),childpat(),parentpat(),trf(){}
};

// # 0: strand; 1: start; 2: end; 3: nreads; 4: nreads_good;
struct CJunction:public GSeg {
	char strand; //-1,0,1
	char guide_match; //exact match of a ref intron?
	double nreads;
	double nreads_good;
	CJunction(int s=0,int e=0, char _strand=0):GSeg(s,e),
			strand(_strand), guide_match(0), nreads(0), nreads_good(0) {}
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
	GBamRecord* brec;
	char strand; //-1, 0, 1
	int nh;
	int hi;
	GPVec<CJunction> juncs;
	//GPVec< GVec<RC_ExonOvl> > g_exonovls; //>5bp overlaps with guide exons, for each read "exon"
	GReadAlnData(GBamRecord* bamrec=NULL, char nstrand=0, int num_hits=0, int hit_idx=0):brec(bamrec),
			strand(nstrand), nh(num_hits), hi(hit_idx), juncs(true) { } //, g_exonovls(true)
};

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

// bundle data structure, holds all input data parsed from BAM file
// - r216 regression
struct BundleData {
 BundleStatus status;
 //int64_t bamStart; //start of bundle in BAM file
 int idx; //index in the main bundles array
 int start;
 int end;
 bool covSaturated;
 int numreads;
 int num_fragments; //aligned read/pairs
 int frag_len;
 GStr refseq;
 GList<CReadAln> readlist;
 GVec<float> bpcov;
 GList<CJunction> junction;
 GPVec<GffObj> keepguides;
 GPVec<CTCov> covguides;
 GList<CPrediction> pred;
 RC_BundleData* rc_data;
 BundleData():status(BUNDLE_STATUS_CLEAR), idx(0), start(0), end(0),
		 covSaturated(false), numreads(0), num_fragments(0), frag_len(0),refseq(), readlist(false,true),
		 bpcov(1024), junction(true, true, true), keepguides(false), pred(false), rc_data(NULL) { }

 void getReady(int currentstart, int currentend) {
	 start=currentstart;
	 end=currentend;
	 //refseq=ref;
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

 //bool evalReadAln(GBamRecord& brec, char& strand, int nh); //, int hi);
 bool evalReadAln(GReadAlnData& alndata, char& strand);

 void Clear() {
	keepguides.Clear();
	pred.Clear();
	readlist.Clear();
	bpcov.Clear();
	bpcov.setCapacity(1024);
	junction.Clear();
	start=0;
	end=0;
	status=BUNDLE_STATUS_CLEAR;
	covSaturated=false;
	numreads=0;
	num_fragments=0;
	frag_len=0;
	delete rc_data;
	rc_data=NULL;
 }

 ~BundleData() {
	Clear();
 }
};

int processRead(int currentstart, int currentend, BundleData& bdata,
		 GHash<int>& hashread, GReadAlnData& alndata);
		 //GBamRecord& brec, char strand, int nh, int hi);

void countFragment(BundleData& bdata, GBamRecord& brec, int hi);

int printResults(BundleData* bundleData, int ngenes, int geneno, GStr& refname);

int infer_transcripts(BundleData* bundle, bool fast);

// --- utility functions
void printGff3Header(FILE* f, GArgs& args);

void printTime(FILE* f);


#endif
