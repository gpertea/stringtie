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


void Ballgown_setupFiles(FILE* &f_tdata, FILE* &f_edata, FILE* &f_idata,
	            FILE* &f_e2t, FILE* &f_i2t);

//Bundle raw count data

struct RC_ScaffSeg {
   uint id; //feature id (>0)
   int l; int r; //genomic coordinates
   char strand;
   bool operator<(const RC_ScaffSeg& o) const {
     //if (id == o.id) return false;
     if (l != o.l) return (l < o.l);
     if (r != o.r) return (r < o.r);
     if (strand == '.' || o.strand == '.') return false;
     if (strand != o.strand) return (strand < o.strand);
     return false;
   }

   bool operator==(const RC_ScaffSeg& o) const {
     //if (id == o.id) return true;
     return (l==o.l && r==o.r &&
         (strand == o.strand || strand == '.' || o.strand == '.'));
   }

   RC_ScaffSeg(int fl=0, int fr=0, char s='.', int fid=0) : id(fid),
         l(fl), r(fr), strand(s) { }

};

struct RC_Feature { //exon or intron of a reference transcript
	uint id; //feature id (>0)
	uint t_id; //transcript id;
	int l; int r; //genomic coordinates
	char strand;
	mutable uint rcount; //# reads covering this feature
	mutable uint ucount; //# uniquely mapped reads covering this feature
	mutable double mrcount; //multi-mapping-weighted counts
    double avg;
    double stdev;
    double mavg;
    double mstdev;

    //mutable vector<int> coverage; //per-base exon coverage data
	struct PCompare {
	 bool operator()(const RC_Feature* p1, const RC_Feature* p2) {
	 return (*p1 < *p2);
	 }
	};

	RC_Feature(int l0=0, int r0=0, char s='.', uint fid=0, uint tid=0): id(fid), t_id(tid), l(l0), r(r0),
		strand(s), rcount(0),ucount(0),mrcount(0), avg(0), stdev(0), mavg(0), mstdev(0) {
	if (l>r) { int t=l; l=r; r=t; }
	}

	RC_Feature(RC_ScaffSeg& seg, uint tid=0): id(seg.id), t_id(tid), l(seg.l), r(seg.r),
		strand(seg.strand), rcount(0),ucount(0),mrcount(0), avg(0), stdev(0), mavg(0), mstdev(0) {
	if (l>r) { int t=l; l=r; r=t; }
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


typedef set<const RC_Feature*, RC_Feature::PCompare> RC_FeatPtrSet;
typedef set<RC_Feature>::iterator RC_FeatIt;
typedef map<uint, set<uint> > RC_Map2Set;
typedef map<uint, set<uint> >::iterator RC_Map2SetIt;

struct RC_Seg { //just a genomic interval holder
	int l;
	int r;
	RC_Seg(int l0=0, int r0=0):l(l0), r(r0) { }
};

struct RC_ScaffData { //storing RC data for a transcript
	GffObj* scaff;
	uint t_id;
	GStr t_name; //original GFF ID for the transcript
	int l;
	int r;
	//char strand;
	int num_exons;
	int eff_len;
	double cov;
	double fpkm;
	//other mutable fields here, to be updated by rc_update_scaff()
	char strand;
	//vector<RC_ScaffSeg> exons;
	//vector<RC_ScaffSeg> introns;
    GPVec<RC_Feature> t_exons;
    GPVec<RC_Feature> t_introns;
	//RC_ScaffIds(uint id=0, char s='.'):t_id(id),strand(s) { }
	void rc_addFeatures(uint& c_e_id, set<RC_ScaffSeg>& fexons, GPVec<RC_Feature>& edata,
	                      uint& c_i_id, set<RC_ScaffSeg>& fintrons, GPVec<RC_Feature>& idata);
	void addFeature(int fl, int fr, GPVec<RC_Feature>& fvec, uint& f_id,
			          set<RC_ScaffSeg>& fset, set<RC_ScaffSeg>::iterator& fit, GPVec<RC_Feature>& fdata);
	RC_ScaffData(GffObj& s, uint id=0):scaff(&s), t_id(id), t_name(s.getID()), l(s.start), r(s.end),
			num_exons(s.exons.Count()), eff_len(s.covlen), cov(0), fpkm(0), strand(s.strand),
			t_exons(false), t_introns(false) {
	  /*RC_ScaffIds& sdata = *(scaff->rc_id_data());
	  t_id = sdata.t_id;
	  t_name=scaff->annotated_trans_id();
	  strand=sdata.strand;
	  l=scaff->left();
	  r=scaff->right();
	  num_exons=scaff->exons.Count();
	  strand=scaff->strand;
	  for (size_t i=0;i<exons.size();++i) {
		 RC_ScaffSeg& exon = sdata.exons[i];
		 eff_len+=exon.r-exon.l;
	  }
	  */

	}

    bool operator<(const RC_ScaffData& o) const {
    	if (l != o.l) return (l < o.l);
    	if (r != o.r) return (r < o.r);
    	if (strand != o.strand) return (strand < o.strand);
	    return (t_name < o.t_name);
		return false;
    }
    bool operator==(const RC_ScaffData& o) const {
    	if (t_id!=0 && o.t_id!=0 && t_id!=o.t_id) return false;
    	return (l==o.l && r==o.r && strand == o.strand &&
    			t_name == o.t_name);
    }
};

FILE* rc_fwopen(const char* fname);
FILE* rc_frenopen(const char* fname);
void rc_frendel(const char* fname);

class BundleData;

//void rc_write_counts(const char* refname, BundleData& bundle);

void rc_writeRC(GPVec<RC_ScaffData>& RC_data,
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
 //set<RC_ScaffData> tdata; //all transcripts in this bundle
 //map<uint, set<uint> > e2t; //mapping exon ID to transcript IDs
 //map<uint, set<uint> > i2t; //mapping intron ID to transcript IDs
 //set<RC_Feature> exons; //all exons in this bundle, by their start coordinate
 //set<RC_Feature> introns; //all introns in this bundle, by their start coordinate
 //GList<RC_ScaffData> tdata;
 GPVec<RC_ScaffData> tdata;
 GList<RC_Feature> exons;
 GList<RC_Feature> introns;
 //RC_FeatIt xcache; //cache the first exon overlapping xcache_pos to speed up exon-overlap queries (findExons())
 int xcache; //exons index of the first exon overlapping xcache_pos
 int xcache_pos; // left coordinate of last cached exon overlap query (findExons())
 // -- output files
 /*
 FILE* ftdata; //t_data
 FILE* fedata; //e_data
 FILE* fidata; //i_data
 FILE* fe2t;   //e2t
 FILE* fi2t;   //i2t
 */
 vector<float> f_mcov; //coverage data, multi-map aware, per strand
 vector<int> f_cov;
 vector<float> r_mcov; //coverage data on the reverse strand
 vector<int> r_cov;
 //
 RC_BundleData(int t_l=0, int t_r=0):init_lmin(0), lmin(t_l), rmax(t_r),
	 tdata(false), // e2t(), i2t(), exons(), introns(),
	 exons(true, false, true), introns(true,false,true),
	 xcache(0), xcache_pos(0)
     //, ftdata(NULL), fedata(NULL), fidata(NULL), fe2t(NULL), fi2t(NULL)
	 {
	 if (rmax>lmin) updateCovSpan();
 }

 ~RC_BundleData() {
	 f_cov.clear();
	 f_mcov.clear();
	 r_cov.clear();
	 r_mcov.clear();
 }

/*
 void addBundleFeature(uint t_id, int l, int r, char strand, uint f_id, set<RC_Feature>& fset,
	                         map<uint, set<uint> >& f2t) {
   RC_Feature feat(l, r, strand, f_id);
   fset.insert(feat);
   //pair<RC_FeatIt, bool> in = fset.insert(feat);
   //if (!in.second) { //existing f_id
   // f_id=in.first->id;
   //}
   set<uint> tset;
   tset.insert(t_id);
   pair<RC_Map2SetIt, bool> mapin=f2t.insert(pair<uint, set<uint> >(f_id, tset));
   if (!mapin.second) {
	 //existing f_id
	 (*mapin.first).second.insert(t_id);
   }
  }
*/

 void addTranscript(GffObj& t) {
   //if (!ps.rc_id_data()) return;
   //RC_ScaffIds& sdata = *(ps.rc_id_data());
   GASSERT(t.uptr);
   RC_ScaffData& sdata=*(RC_ScaffData*)(t.uptr);
   //tdata.insert(sdata);
   tdata.Add(&sdata);
   bool boundary_changed=false;
   if (lmin==0 || lmin>(int)t.start) { lmin=t.start; boundary_changed=true; }
   if (rmax==0 || rmax<(int)t.end) { rmax=t.end; boundary_changed=true; }
   if (boundary_changed) updateCovSpan();
   //for (vector<RC_ScaffSeg>::iterator it=sdata.exons.begin();it!=sdata.exons.end();++it) {
   for (int i=0;i<sdata.t_exons.Count();i++) {
	 //addBundleFeature(sdata.t_id, sdata.exons[i], sdata.strand, exons);
	 exons.Add(sdata.t_exons[i]);
   }
   //store introns:
   //for (vector<RC_ScaffSeg>::iterator it=sdata.introns.begin();it!=sdata.introns.end();++it) {
   //   addBundleFeature(sdata.t_id, it->l, it->r, sdata.strand, it->id, introns, i2t);
   for (int i=0;i<sdata.t_introns.Count();i++) {
	   introns.Add(sdata.t_introns[i]);
   }
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

 /*
   RC_FeatIt findIntron(int hl, int hr, char strand) {
   RC_FeatIt ri=introns.find(RC_Feature(hl, hr, strand));
   return ri;
 */
 RC_Feature* findIntron(int hl, int hr, char strand) {
   int fidx=0;
   RC_Feature* r=NULL;
   RC_Feature t(hl, hr, strand);
   if (introns.Found(&t, fidx))
	   r=introns[fidx];
   return r;
 }
}; //struct RC_BundleData

void rc_update_exons(RC_BundleData& rc);

#endif /* TABLEMAKER_H_ */
