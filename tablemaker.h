// tablemaker.h
 
#ifndef TABLEMAKER_H_
#define TABLEMAKER_H_
#include <vector>
#include <map>
#include <set>
//#include <string>
#include <algorithm>
using namespace std;

#define RC_MIN_EOVL 5 //minimum overlap for exons

extern bool ballgown;

void Ballgown_setupFiles(FILE* &f_tdata, FILE* &f_edata, FILE* &f_idata,
	            FILE* &f_e2t, FILE* &f_i2t);

struct RC_Feature { //read count data for exon or intron of a reference transcript
	uint id; //feature id (>0), +1 to the index either in global guides_RC_exons/introns if ballgown,
	         //                                 or in bundle_RC_exons/introns if not ballgown
	bool nascent_only=true; // if feature is only in nascent transcripts
	GVec<uint> t_ids; //transcripts owning this feature
	 //if -B, this is the index in the global refguides_RC_Data array + 1
	 // otherwise it is the index in the BundleData::keepguides array + 1
	long l=0; 
	long r=0; //genomic coordinates for the feature
	char strand=0;
	mutable uint rcount=0; //# read alignments overlapping this feature (>5bp overlaps for exons;
	                     // exact coord. match for introns)
	mutable uint ucount=0; //# uniquely mapped reads overlapping/matching this ref feature
	mutable double mrcount=0; //multi-map weighted read counts overlapping/matching this feature

	mutable double movlcount=0; //exons only: multi-map weighted sum of overlap lengths

    double avg=0;
    double stdev=0;
    double mavg=0;
    double mstdev=0;
	struct PCompare {
	 bool operator()(const RC_Feature* p1, const RC_Feature* p2) {
	 return (*p1 < *p2);
	 }
	};

	RC_Feature(long l0=0, long r0=0, char s='.', uint fid=0, uint tid=0): id(fid), t_ids(1), l(l0), r(r0),
		                           strand(s) {
	if (l>r) Gswap(l,r);
	if (tid>0) t_ids.Add(tid);
	}
	RC_Feature(const RC_Feature& seg): id(seg.id), t_ids(seg.t_ids), l(seg.l), r(seg.r),
			strand(seg.strand) {
	}

	RC_Feature(const RC_Feature& seg, uint tid): id(seg.id), t_ids(1), l(seg.l), r(seg.r),
		strand(seg.strand) {
	if (l>r) Gswap(l,r);
	if (tid>0) t_ids.Add(tid);
	}

	bool operator<(const RC_Feature& o) const {
	 //if (id == o.id) return false;
	 if (l != o.l) return (l < o.l);
     if (r != o.r) return (r < o.r);
     if (strand == '.' || o.strand == '.') return false;
     if (strand != o.strand) return (strand < o.strand);
     return false;
	 }
	bool operator==(const RC_Feature& o) const {
	 //if (id == o.id) return true;
	 return (l==o.l && r==o.r &&
		 (strand == o.strand || strand == '.' || o.strand == '.'));
	 }
	bool strand_compatible(const RC_Feature& o) const {
		 return (strand == '.' || o.strand == '.' || strand == o.strand);
	}
	//WARNING: the overlap checks IGNORE strand!
	bool overlap(long hl, long hr) const {
	  if (hl>hr) { long t=hl; hl=hr; hr=t; }
      return (l<=hr && r<=hl);
	  }
	bool overlap(long hl, long hr, long minovl) const {
	  if (hl>hr) { long t=hl; hl=hr; hr=t; }
      hl+=minovl;hr-=minovl;
      return (l<=hr && r<=hl);
	  }
	long ovlen(long hl, long hr) const {
     if (hl>hr) Gswap(hl,hr); 
     if (l<hl) {
        if (hl>r) return 0;
        return (hr>r) ? r-hl+1 : hr-hl+1;
        }
       else { //hl<=l
        if (l>hr) return 0;
        return (hr<r)? hr-l+1 : r-l+1;
        }
	 }
};

struct RC_ExonOvl {
	RC_Feature* feature=NULL; //pointer to an item of RC_BundleData::g_exons	
	byte in_guides = 0; // =1 if feature is not only in nascent transcripts
	byte mate_ovl=0; // =1 if the mate of this read overlaps the same exon
	long ovlen=0;
	bool operator<(const RC_ExonOvl& o) const {
		if (in_guides!=o.in_guides)
			return (in_guides>o.in_guides);
		if (mate_ovl!=o.mate_ovl)
			return (mate_ovl>o.mate_ovl);
		if (ovlen!=o.ovlen)
			return (ovlen>o.ovlen);
		if (feature->r-feature->l != o.feature->r-o.feature->l)
			return (feature->r-feature->l > o.feature->r-o.feature->l);
		if (feature->strand != o.feature->strand)
			return (feature->strand<o.feature->strand);
		return (feature->l<o.feature->l);
	} //operator <
	bool operator==(const RC_ExonOvl& o) const {
		return (mate_ovl==o.mate_ovl && ovlen==o.ovlen && feature==o.feature);
	}
	RC_ExonOvl(RC_Feature* f=NULL, long olen=0, byte movl=0):feature(f),
			mate_ovl(movl), ovlen(olen) {
		if (f) in_guides=!f->nascent_only;
	}
};

//typedef set<const RC_Feature*, RC_Feature::PCompare> RC_FeatPtrSet;
typedef set<RC_Feature>::iterator RC_FeatIt;
typedef map<uint, set<uint> > RC_Map2Set;
typedef map<uint, set<uint> >::iterator RC_Map2SetIt;

struct RC_TData { //storing RC data for a transcript (guide/nascent)
    // uptr data block created for each guide, even those not used in a bundle
	// but only fully used with -B (Ballgown Raw/Read Count data collection)
	GffObj* ref_t;
	uint t_id; //index in the global RC_tdata table
	long l;
	long r;
	//char in_bundle; // 1 if used in a bundle (guide added to keepguides, default value)
	                // 2 if all introns are covered by at least one read, 
					// 3 if it is stored to be printed
			//NOTE: superseded by setGuideStatus/getGuideStatus(GffObj* t)
	int eff_len;
	double cov;
	double fpkm;
	GffObj* gen_from; // <-- for nascents, the guide transcript that generated this
	//char strand;
    GPVec<RC_Feature> t_exons;
    GPVec<RC_Feature> t_introns;
	void rc_addFeatures(uint& c_e_id, GList<RC_Feature>& fexons, GPVec<RC_Feature>& edata,
	                      uint& c_i_id, GList<RC_Feature>& fintrons, GPVec<RC_Feature>& idata);
	void addFeature(long fl, long fr, GPVec<RC_Feature>& fvec, uint& f_id,
			          GList<RC_Feature>& fset, GPVec<RC_Feature>& fdata,
					  int& cache_idx);
	RC_TData(GffObj& s, uint id=0):ref_t(&s), t_id(id), l(s.start), r(s.end), //in_bundle(0), 
	        eff_len(s.covlen), cov(0), fpkm(0), gen_from(NULL), //strand(s.strand),
			t_exons(false), t_introns(false) { //, e_idx_cache(-1), i_idx_cache(-1) {
	}
    
    bool operator<(const RC_TData& o) const {
    	if (l != o.l) return (l < o.l);
    	if (r != o.r) return (r < o.r);
    	if (char c=(ref_t->strand - o.ref_t->strand)) return (c<0);
	    return (strcmp(ref_t->getID(), o.ref_t->getID())<0);
    }
    bool operator==(const RC_TData& o) const {
    	if (t_id!=0 && o.t_id!=0) return (t_id==o.t_id);
    	return (l==o.l && r==o.r && ref_t->strand == o.ref_t->strand &&
    			strcmp(ref_t->getID(),o.ref_t->getID())==0);
    }
};

FILE* rc_fwopen(const char* fname);
FILE* rc_frenopen(const char* fname);
void rc_frendel(const char* fname);

struct BundleData;

void rc_writeRC(GPVec<RC_TData>& RC_data,
		GPVec<RC_Feature>& RC_exons,
		GPVec<RC_Feature>& RC_introns,
		FILE* &f_tdata, FILE* &f_edata, FILE* &f_idata,
        FILE* &f_e2t, FILE* &f_i2t);

int rc_cov_inc(int i);

class RC_MultiCovInc {
	float fcov;
  public:
	RC_MultiCovInc(int numhits):fcov(1.0) {
	 if (numhits>1) fcov=1/(float)numhits;
	}
	float operator()(const float& v) {
	  return (v+fcov);
	}
};

struct RC_BundleData {
 long init_lmin;
 long lmin;
 long rmax;
 GPVec<RC_TData> g_tdata; //raw counting data for all transcripts in this bundle
 // RC_TData::t_id-1 = the element index in this array
 GList<RC_Feature> g_exons; //set of guide exons in this bundle, sorted by start coordinate
 GList<RC_Feature> g_introns; //set of guide introns in this bundle, sorted by start coordinate
 //RC_FeatIt xcache; //cache the first exon overlapping xcache_pos to speed up exon-overlap queries (findOvlExons())
 int xcache; //exons index of the first exon overlapping xcache_pos
 long xcache_pos; // left coordinate of last cached exon overlap query (findOvlExons())

 // the following coverage arrays will only used with Ballgown data (-B)
 vector<float> f_mcov; //coverage data, multi-map aware, per strand
 vector<int> f_cov;
 vector<float> r_mcov; //coverage data on the reverse strand
 vector<int> r_cov;

 //-- when no global Ballgown data is to be generated, these are
 // local (bundle) stable order tables of guide features
 GPVec<RC_TData>* bundle_RC_tdata; //pointer to the global RC tdata table
   // RC_Feature::id-1 = the index in these arrays:
 GPVec<RC_Feature>* bundle_RC_exons;  //pointers to global (if ballgown)
 GPVec<RC_Feature>* bundle_RC_introns;// OR local exon/intron RC data
 //local exon/intron ids within the bundle
 uint c_exon_id;
 uint c_intron_id;

 
 //--
 RC_BundleData(long t_l=0, long t_r=0, GPVec<RC_TData>* rc_tdata=NULL,
		 GPVec<RC_Feature>* rc_exons=NULL,GPVec<RC_Feature>* rc_introns=NULL):
			 init_lmin(0), lmin(t_l), rmax(t_r),  g_tdata(false),
	           // features:(sorted, free, unique)
	 g_exons(true, false, true), g_introns(true, false, true),
			 xcache(0), xcache_pos(0),
			 bundle_RC_tdata(rc_tdata), bundle_RC_exons(rc_exons), bundle_RC_introns(rc_introns),
			 c_exon_id(0), c_intron_id(0)
 {
	 if (ballgown) {
		 if (rmax>lmin) updateCovSpan();
	 }else {
		 //g_tdata.setFreeItem(true);
		 //guides_RC_tdata   = &g_tdata;
		 //-- override the passed rc_exons/rc_introns if not ballgown
		 //-- because these are now locally maintained so they'll be deallocated with the bundle
		 bundle_RC_exons  = new GPVec<RC_Feature>(true);
		 bundle_RC_introns= new GPVec<RC_Feature>(true);
	 }
 }

 ~RC_BundleData() {
	 f_cov.clear();
	 f_mcov.clear();
	 r_cov.clear();
	 r_mcov.clear();
	 if (!ballgown) {
		 delete bundle_RC_exons;
		 delete bundle_RC_introns;
	 }
 }
 
 uint addTranscript(GffObj& t) { //should return the guide index in *guides_RC_tdata
	bool boundary_changed=false;
	if (lmin==0 || lmin>(long)t.start) { lmin=t.start; boundary_changed=true; }
	if (rmax==0 || rmax<(long)t.end) { rmax=t.end; boundary_changed=true; }
	GASSERT(t.uptr); //we should always have a RC_TData for each guide
	RC_TData* tdata=(RC_TData*)(t.uptr);
	//tdata->in_bundle=1; //don't tag here, it might be in a read-no-overlap bundle
	g_tdata.Add(tdata);
	if (ballgown) {
	   if (boundary_changed) updateCovSpan();
	   //rc_addFeatures() called before, but we still need to add exons
	   // and introns to the local sets: g_exons, g_introns
	   for (int i=0;i<tdata->t_exons.Count();i++) {
		 g_exons.Add(tdata->t_exons[i]);
	   }
	   for (int i=0;i<tdata->t_introns.Count();i++) {
		   g_introns.Add(tdata->t_introns[i]);
	   }
	}
	else { //no ballgown counting, but we still want to get exon and intron read coverage info
		tdata->rc_addFeatures(c_exon_id, g_exons, *bundle_RC_exons,
		        c_intron_id, g_introns, *bundle_RC_introns);
	}
 return tdata->t_id;
 }

 void updateCovSpan() {
	 //ideally this should be called after all reference transcripts were added
	 // should NEVER be called repeatedly, for the same bundle, with a different lmin !
	 GASSERT(rmax>lmin);
	 int blen=rmax-lmin+1;
	 if (init_lmin==0) init_lmin=lmin;
	 else {
		 if (lmin!=init_lmin) //this should never happen
			 GError("Error setting up Ballgown coverage data (lmin should never change!) !\n");
	 }
	 f_cov.resize(blen, 0);
	 r_cov.resize(blen, 0);
	 f_mcov.resize(blen, 0.0);
	 r_mcov.resize(blen, 0.0);
 }

 void updateCov(char strand, int numhits, long gpos, int glen) {
  if (gpos>rmax || gpos+glen<lmin) return; //no overlap with bundle
  if (gpos<lmin) { //this read maps before the bundle start (left overhang)
	int gadj=lmin-gpos;
	gpos+=gadj;
	glen-=gadj;
  }
  if (gpos+glen>rmax) {
	glen=rmax-gpos;
  }
  if (glen<=0) return; //no overlap (should not happen here)
  int goffs=gpos-lmin;
  if (goffs<0) return; //should not happen
  if (strand=='.' || strand=='+') {
	transform(f_cov.begin()+goffs, f_cov.begin()+goffs+glen,
		        f_cov.begin()+goffs, rc_cov_inc);
    transform(f_mcov.begin()+goffs, f_mcov.begin()+goffs+glen,
    	        f_mcov.begin()+goffs, RC_MultiCovInc(numhits));
  }
  if (strand=='.' || strand=='-') {
	transform(r_cov.begin()+goffs, r_cov.begin()+goffs+glen,
		        r_cov.begin()+goffs, rc_cov_inc);
    transform(r_mcov.begin()+goffs, r_mcov.begin()+goffs+glen,
    	        r_mcov.begin()+goffs, RC_MultiCovInc(numhits));
  }
 }

 bool findOvlExons(GArray<RC_ExonOvl>& exovls, long hl, long hr, char strand='.',
		                                    int mate_pos=0, bool update_cache=true) {
	 //exovls should be clear (unless you want to append to an existing list)
	 bool hasOverlaps=false;
	 if (g_exons.Count()==0) return false;
	 int xstart=0;
	 bool no_cache=(xcache_pos==0 || xcache_pos>hl);
	 if (no_cache) {
		 if (update_cache) {
			 xcache=g_exons.Count()-1;
			 xcache_pos=0;
		 }
	 }
	 else xstart=xcache; //must have a valid value
	 bool upd_cache(update_cache);
	 int last_checked_exon=g_exons.Count()-1;
	 int vlen=hr-hl+1;
	 int min_ovl=GMIN(vlen, RC_MIN_EOVL);
	 for (long p=xstart;p < g_exons.Count();++p) {
		 last_checked_exon=p;
		 long l=g_exons[p]->l;
		 long r=g_exons[p]->r;
		 if (l > hr) break;
		 if (hl > r) continue;
		 //exon overlap here
		 long ovlen=0;
		 if (hl<l) {
			 ovlen = ( hr<r ? hr-l+1 : r-l+1 );
		 }
		 else { // l<=hl
			 ovlen= ( hr<r ? hr-hl+1 : r-hl+1 );
		 }
		 if (upd_cache) {
			 //cache first overlap
			 xcache=p;
			 upd_cache=false;
		 }
		 if (strand!='.' && strand!=g_exons[p]->strand) continue; //non-matching strand
		 byte mate_ovl=0;
		 if (mate_pos && mate_pos+10>l && mate_pos+5<r)
			    mate_ovl=1; //mate read likely overlaps this exon
		 if (mate_ovl || ovlen>=min_ovl) {
			 hasOverlaps=true;
			 RC_ExonOvl fovl(g_exons[p], ovlen, mate_ovl);
			 exovls.Add(fovl);
		 }
	 }
	 if (update_cache) {
		 if (upd_cache) xcache=last_checked_exon; //there was no overlap found
		 xcache_pos=hl;
	 }
	 return hasOverlaps;
 }

 RC_Feature* findIntron(long hl, long hr, char strand) {
   int fidx=0;
   RC_Feature* r=NULL;
   RC_Feature t(hl, hr, strand);
   if (g_introns.Found(&t, fidx))
	   r=g_introns[fidx];
   return r;
 }
}; //struct RC_BundleData

void rc_update_exons(RC_BundleData& rc);
void rc_updateExonCounts(const RC_ExonOvl& exonovl, int nh);
#endif /* TABLEMAKER_H_ */
