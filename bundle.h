#include "gff.h"
#include "GSam.h"
#include "GStr.h"

extern bool mergeMode;

#include "tablemaker.h"

extern bool genNascent; // generate nascent synthetic transcripts for each bundle

void setNascent(GffObj* t, byte v=1); // set/clear the synthetic nascent RNA flag
//0 = not nascent, 1 = synthetic nascent, 2 = nascent replacing a guide
byte isNascent(GffObj* t); // check if a transcript is a synthetic nascent RNA

// set/getBundleFlag - should set/get  (RC_TData*)(keepguides[i]->uptr)->in_bundle
//        1 = default, 2 = every intron covered by a read, 3 = stored to be printed

//GffObj* nascentFrom(GffObj* ntx); // get the guide transcript that a nascent transcript was generated from
inline GffObj* nascentFrom(GffObj* ntx) {
    if (ntx && ntx->uptr) {
        return ((RC_TData*)(ntx->uptr))->gen_from;
    }
    return NULL;
}

enum GuideBundleStatus {
  GBST_UNSET = 0,
  GBST_IN_BUNDLE,    // 1: added to bundle
  GBST_ALL_INTR_COV, // 2: all introns covered by at least one read
  GBST_STORED,     // 3: stored to be printed
};
GuideBundleStatus getGuideStatus(GffObj* t);
void setGuideStatus(GffObj* t, GuideBundleStatus status);


//collect all refguide transcripts for a single genomic sequence
struct GRefData {
  GList<GffObj> rnas; //all guides on this genomic seq; they are added in order as sorted in GffReader
  GList<GffObj> synrnas; //just keep track/store of synthetic nascent transcripts generated on this chr
  int gseq_id; //chromosome index in GffObj::names->gseqs
  const char* gseq_name; //chromosome name
                          // rnas list: unsorted, dealloc on free, duplicates allowed
  GRefData(int gid=-1):rnas(false, true, false), synrnas(false, true, false),
                       gseq_id(gid),gseq_name(NULL) {
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
		if (genNascent) synrnas.setCapacity(gffr->gseqtable[gseq_id]->fcount);
     }
     rnas.Add(t);
     t->isUsed(true); //mark as used, to prevent deletion in GffReader destructor by gflst::freeUnused()
  }

  bool operator==(GRefData& d){
    return gseq_id==d.gseq_id;
  }
  bool operator<(GRefData& d){
    return (gseq_id<d.gseq_id);
  }
};


struct Ref_RC_Data {
     GRefData* refdata=nullptr;
	 GPVec<RC_TData>*   rc_tdata=nullptr;
	 GPVec<RC_Feature>* rc_edata =nullptr;
	 GPVec<RC_Feature>* rc_idata = nullptr;
 };

enum GPFType {
	GPFT_NONE=0,
	GPFT_TSS,
	GPFT_CPAS,
	GPFT_JSTART,
	GPFT_JEND
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
	CPrediction* linkpred; // for nascent RNAs prediction of transcript that it is linked to and viceversa
	GVec<GSeg> exons;
	GVec<float> exoncov;
	GStr mergename;
	CPrediction(int _geneno=0, GffObj* guide=NULL, int gstart=0, int gend=0, float _cov=0, char _strand='.',
	int _len=0,bool f=true, CPrediction* lp=NULL):GSeg(gstart,gend), geneno(_geneno),t_eq(guide),cov(_cov),longcov(0),strand(_strand),
	//CPrediction(int _geneno=0, char* _id=NULL,int gstart=0, int gend=0, float _cov=0, char _strand='.', float _frag=0,
	//		int _len=0,bool f=true):GSeg(gstart,gend), geneno(_geneno),id(_id),cov(_cov),strand(_strand),frag(_frag),
			tlen(_len),flag(f),linkpred(lp),exons(),exoncov(),mergename() {}
	void init(int _geneno=0, GffObj* guide=NULL, int gstart=0, int gend=0, float _cov=0, char _strand='.',
	          int _len=0,bool f=true, CPrediction* lp=NULL) {
		geneno=_geneno;
		t_eq=guide;
		start=gstart;
		end=gend;
		cov=_cov;
		strand=_strand;
		tlen=_len;
		flag=f;
		linkpred=lp;
		exons.Clear();
		exoncov.Clear();
		mergename.clear();
	}

	CPrediction(CPrediction& c):GSeg(c.start, c.end), geneno(c.geneno),
//			id(Gstrdup(c.id)), cov(c.cov), strand(c.strand), frag(c.frag), tlen(c.tlen), flag(c.flag),
			t_eq(c.t_eq), cov(c.cov), longcov(c.longcov),strand(c.strand), tlen(c.tlen), flag(c.flag),linkpred(c.linkpred),
	      exons(c.exons),  exoncov(c.exoncov), mergename(c.mergename) {}
	~CPrediction() { //GFREE(id);
		}
};

// bundle data structure, holds all input data parsed from BAM file
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
 GPVec<GffObj> keepguides; //list of guides in this bundle (+ synthetic nascents if genNascent)
 GPVec<GPtFeature> ptfs; //point features for this bundle
 GList<CPrediction> pred;
 int numNascents=0; //number of nascent transcripts generated for this bundle
 RC_BundleData* rc_data; // read count data for this bundle
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
		 //RC_TData* tdata=(RC_TData*)(keepguides[i]->uptr);
		 //tdata->in_bundle=1;
		 setGuideStatus(keepguides[i], GBST_IN_BUNDLE);
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
 void keepGuide(GffObj* scaff, Ref_RC_Data& ref_rc_data);
 void generateAllNascents(int from_guide_idx, Ref_RC_Data& ref_rc); //defined in tablemaker.cpp
 // use this for debug only
 void printBundleGuides() {
	 GStr fname("bundle");
	 fname.appendfmt("_%d_guides.bed",idx);
	 GStr fnasc("bundle");
	 fnasc.appendfmt("_%d_nascents.bed",idx);
	 FILE* f=fopen(fname.chars(),"w");
	 FILE* fn=fopen(fnasc.chars(),"w");
	 for(int i=0;i<keepguides.Count();i++) {
		 GffObj* g=keepguides[i];
		 if (isNascent(g)) {
			g->addAttr("gen_from", nascentFrom(g)->getID());
			g->printBED(fn);
		 } else g->printBED(f);
	 }
	 fclose(f);
	 fclose(fn);
 }


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
	numNascents=0;
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
		 GHash<int>& hashread, GReadAlnData& alndata,bool ovlpguide);
		 //GSamRecord& brec, char strand, int nh, int hi);

