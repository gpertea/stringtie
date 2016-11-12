/*
 * tablemaker.h
 *
 *  Created on: Oct 26, 2014
 *      Author: gpertea
 */

#ifndef TABLEMAKER_H_
#define TABLEMAKER_H_
#include <vector>
#include <map>
#include <set>
//#include <string>
#include <algorithm>
using namespace std;

#define RC_MIN_EOVL 5

extern bool ballgown;

void Ballgown_setupFiles(FILE* &f_tdata, FILE* &f_edata, FILE* &f_idata,
	            FILE* &f_e2t, FILE* &f_i2t);

//Bundle raw count data
/*
struct RC_TSeg {
   uint id; //feature id (>0)
   int l; int r; //genomic coordinates
   char strand;
   bool operator<(const RC_TSeg& o) const {
     //if (id == o.id) return false;
     if (l != o.l) return (l < o.l);
     if (r != o.r) return (r < o.r);
     if (strand == '.' || o.strand == '.') return false;
     if (strand != o.strand) return (strand < o.strand);
     return false;
   }

   bool operator==(const RC_TSeg& o) const {
     //if (id == o.id) return true;
     return (l==o.l && r==o.r &&
         (strand == o.strand || strand == '.' || o.strand == '.'));
   }

   RC_TSeg(int fl=0, int fr=0, char s='.', int fid=0) : id(fid),
         l(fl), r(fr), strand(s) { }

};
*/

struct RC_Feature { //exon or intron of a reference transcript
	uint id; //feature id (>0), +1 to the index either in global guides_RC_exons/introns if ballgown,
	         //                                 or in bundle_RC_exons/introns if not ballgown
	GVec<uint> t_ids; //transcripts owning this feature
	 //if -B, this is the index in the global refguides_RC_Data array + 1
	 // otherwise it is the index in the BundleData::keepguides array + 1
	int l; int r; //genomic coordinates for the feature
	char strand;
	mutable uint rcount; //# read alignments overlapping this feature (>5bp overlaps for exons;
	                     // exact coord. match for introns)
	mutable uint ucount; //# uniquely mapped reads overlapping/matching this ref feature
	mutable double mrcount; //multi-map weighted read counts overlapping/matching this feature

	mutable double movlcount; //exons only: multi-map weighted sum of overlap lengths

    double avg;
    double stdev;
    double mavg;
    double mstdev;

	struct PCompare {
	 bool operator()(const RC_Feature* p1, const RC_Feature* p2) {
	 return (*p1 < *p2);
	 }
	};

	RC_Feature(int l0=0, int r0=0, char s='.', uint fid=0, uint tid=0): id(fid), t_ids(1), l(l0), r(r0),
		strand(s), rcount(0),ucount(0),mrcount(0), movlcount(0), avg(0), stdev(0), mavg(0), mstdev(0) {
	if (l>r) { int t=l; l=r; r=t; }
	if (tid>0) t_ids.Add(tid);
	}
	RC_Feature(const RC_Feature& seg): id(seg.id), t_ids(seg.t_ids), l(seg.l), r(seg.r),
			strand(seg.strand), rcount(0),ucount(0),mrcount(0), movlcount(0), avg(0), stdev(0), mavg(0), mstdev(0) {
	}

	RC_Feature(const RC_Feature& seg, uint tid): id(seg.id), t_ids(1), l(seg.l), r(seg.r),
		strand(seg.strand), rcount(0),ucount(0),mrcount(0), movlcount(0), avg(0), stdev(0), mavg(0), mstdev(0) {
	if (l>r) { int t=l; l=r; r=t; }
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
	bool overlap(int hl, int hr) const {
	  if (hl>hr) { int t=hl; hl=hr; hr=t; }
      return (l<=hr && r<=hl);
	  }
	bool overlap(int hl, int hr, int minovl) const {
	  if (hl>hr) { int t=hl; hl=hr; hr=t; }
      hl+=minovl;hr-=minovl;
      return (l<=hr && r<=hl);
	  }
	uint ovlen(int hl, int hr) const {
     if (hl>hr) { int t=hl; hl=hr; hr=t; }
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
/*
//for locus tracking and coverage, keep merged exons/introns in the locus
struct GSegInfo {
	int t_id;  //index of RC_TData in the guides_RC_data table + 1
	int exonum; //exon number
	GSegInfo(int tid=-1, int en=-1):t_id(tid), exonum(en) { }
};
class GMSeg:public GSeg {
  public:
    GVec<GSegInfo> msegs; //keep track of exons contributing to this merged exon
    GMSeg(int l=0, int r=0, int tid=-1, int eno=-1):GSeg(l,r), msegs(tid, eno) {
    }
};

//reference locus - formed by exon-overlapping transcripts
class GRefLocus:public GSeg {
	GArray<GMSeg> mexons; //merged exons in this locus (from any transcript)
	GPVec<GffObj> rnas; //transcripts in this locus
};
*/

struct RC_ExonOvl {
	RC_Feature* feature; //pointer to an item of RC_BundleData::g_exons
	int mate_ovl; // = 1 if the mate of this read overlaps the same exon
	int ovlen;
	bool operator<(const RC_ExonOvl& o) const {
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
	RC_ExonOvl(RC_Feature* f=NULL, int olen=0, int movl=0):feature(f),
			mate_ovl(movl), ovlen(olen) {
	}
};

//typedef set<const RC_Feature*, RC_Feature::PCompare> RC_FeatPtrSet;
typedef set<RC_Feature>::iterator RC_FeatIt;
typedef map<uint, set<uint> > RC_Map2Set;
typedef map<uint, set<uint> >::iterator RC_Map2SetIt;

/*struct RC_Seg { //just a genomic interval holder
	int l;
	int r;
	RC_Seg(int l0=0, int r0=0):l(l0), r(r0) { }
};
*/

struct RC_TData { //storing RC data for a transcript
	//only used with -B (full Ballgown data)
	GffObj* ref_t;
	uint t_id;
	int l;
	int r;
	char in_bundle; // 1 if used by read bundles (present in keepguides), 2 if all introns are covered by at least one read, 3 if it is stored to be printed
	//GRefLocus* locus; //pointer to a locus info
	int eff_len;
	double cov;
	double fpkm;
	//char strand;
    GPVec<RC_Feature> t_exons;
    GPVec<RC_Feature> t_introns;
	void rc_addFeatures(uint& c_e_id, GList<RC_Feature>& fexons, GPVec<RC_Feature>& edata,
	                      uint& c_i_id, GList<RC_Feature>& fintrons, GPVec<RC_Feature>& idata);
	void addFeature(int fl, int fr, GPVec<RC_Feature>& fvec, uint& f_id,
			          GList<RC_Feature>& fset, GPVec<RC_Feature>& fdata,
					  int& cache_idx);
	RC_TData(GffObj& s, uint id=0):ref_t(&s), t_id(id), l(s.start), r(s.end),
			in_bundle(0), eff_len(s.covlen), cov(0), fpkm(0), //strand(s.strand),
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
 int init_lmin;
 int lmin;
 int rmax;
 GPVec<RC_TData> g_tdata; //raw counting data for all transcripts in this bundle
 // RC_TData::t_id-1 = the element index in this array
 GList<RC_Feature> g_exons; //set of guide exons in this bundle, sorted by start coordinate
 GList<RC_Feature> g_introns; //set of guide introns in this bundle, sorted by start coordinate
 //RC_FeatIt xcache; //cache the first exon overlapping xcache_pos to speed up exon-overlap queries (findOvlExons())
 int xcache; //exons index of the first exon overlapping xcache_pos
 int xcache_pos; // left coordinate of last cached exon overlap query (findOvlExons())

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
 RC_BundleData(int t_l=0, int t_r=0, GPVec<RC_TData>* rc_tdata=NULL,
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
	if (lmin==0 || lmin>(int)t.start) { lmin=t.start; boundary_changed=true; }
	if (rmax==0 || rmax<(int)t.end) { rmax=t.end; boundary_changed=true; }
	GASSERT(t.uptr); //we should always have a RC_TData for each guide
	RC_TData* tdata=(RC_TData*)(t.uptr);
	//tdata->in_bundle=1; //don't tag here, it might be in a read-no-overlap bundle
  /*RC_TData* tdata=NULL;
  if (ballgown) {
     tdata=(RC_TData*)(t.uptr);
  }
  else {
    //add RC transcript data locally for the bundle
    tdata=new RC_TData(t, g_tdata.Count()+1);
    t.uptr=tdata;
      //guides_RC_Data.Add(tdata);
      tdata->rc_addFeatures(c_exon_id, g_exons, *guides_RC_exons,
        c_intron_id, g_introns, *guides_RC_introns);
  }*/
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
	else {
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

 void updateCov(char strand, int numhits, int gpos, int glen) {
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

 bool findOvlExons(GArray<RC_ExonOvl>& exovls, int hl, int hr, char strand='.',
		                                    int mate_pos=0, bool update_cache=true) {
	 //exovls should be clear, unless the caller knows what s/he's doing
	 bool hasOverlaps=false;
	 if (g_exons.Count()==0) return false;
	 RC_Feature q(hl, hr);
	 int xstart=0;
	 bool no_cache=(xcache_pos==0 || xcache_pos>hl);
	 if (no_cache) {
		 if (update_cache) {
			 //xcache=exons.end();
			 xcache=g_exons.Count()-1;
			 xcache_pos=0;
		 }
	 }
	 else xstart=xcache; //must have a valid value
	 bool upd_cache(update_cache);
	 int last_checked_exon=g_exons.Count()-1;
	 for (int p=xstart;p < g_exons.Count();++p) {
		 last_checked_exon=p;
		 int l=g_exons[p]->l;
		 int r=g_exons[p]->r;
		 if (l > hr) break;
		 if (hl > r) continue;
		 //exon overlap here
		 int ovlen=0;
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
		 int mate_ovl=0;
		 if (mate_pos && mate_pos+10>l && mate_pos+5<r)
			    mate_ovl=1; //mate read likely overlaps this exon
		 if (mate_ovl || ovlen>=5) {
			 //TODO: check this, arbitrary ovl minimum of 5bp
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
/*
 RC_FeatPtrSet findExons(int hl, int hr, char strand='.', bool update_cache=true) {
   //returns exons overlapping given interval hl-hr
   RC_FeatPtrSet ovlex; //return set
   if (exons.Count()==0) return ovlex;
   RC_Feature q(hl, hr);
   //RC_FeatIt xstart=exons.begin();
   int xstart=0;
   bool no_cache=(xcache_pos==0 || xcache_pos>hl);
   if (no_cache) {
	   if (update_cache) {
	      //xcache=exons.end();
		  xcache=exons.Count()-1;
	      xcache_pos=0;
	   }
   }
   else xstart=xcache; //must have a valid value
   bool upd_cache(update_cache);
   //RC_FeatIt last_checked_exon(exons.end());
   int last_checked_exon=exons.Count()-1;
   //for (RC_FeatIt p=xstart;p != exons.end();++p) {
   for (int p=xstart;p < exons.Count();++p) {
     last_checked_exon=p;
	 if (exons[p]->l > hr) break;
	 if (hl > exons[p]->r) continue;
     //exon overlap
     if (upd_cache) {
         //cache first overlap
		 xcache=p;
		 upd_cache=false;
     }
     if (strand!='.' && strand!=exons[p]->strand) continue;
	 ovlex.insert(exons[p]);
   }
   if (update_cache) {
	 if (upd_cache) xcache=last_checked_exon; //there was no overlap found
	 xcache_pos=hl;
     }
   return ovlex;
  }

 */
 RC_Feature* findIntron(int hl, int hr, char strand) {
   int fidx=0;
   RC_Feature* r=NULL;
   RC_Feature t(hl, hr, strand);
   if (g_introns.Found(&t, fidx))
	   r=g_introns[fidx];
   return r;
 }
}; //struct RC_BundleData

void rc_update_exons(RC_BundleData& rc);

#endif /* TABLEMAKER_H_ */
