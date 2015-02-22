//#include "tablemaker.h"
#include "rlink.h"
#include <numeric>

extern GStr ballgown_dir;

int rc_cov_inc(int i) {
  return ++i;
}

/*
void rc_update_tdata(BundleData& bundle, GffObj& scaff,
	                          double cov, double fpkm) {
	if (bundle.rc_data==NULL) return;
  RC_BundleData& rc = *(bundle.rc_data);
  if (rc.exons.size()==0) return;
  RC_ScaffData& q = (RC_ScaffData*)scaff.uptr;
  set<RC_ScaffData>::iterator tdata = rc.tdata.find(q);
  if (tdata==rc.tdata.end()) {
	fprintf(stderr, "Error: cannot locate bundle ref. transcript %s (%s:%d-%d)!\n",
		scaff.getID(), scaff.getGSeqName(), scaff.start, scaff.end);
	return;
  }
  (*tdata).cov=cov;
  (*tdata).fpkm=fpkm;
}
*/
void BundleData::rc_store_t(GffObj* t) {
	//if (!rc_stage) return;
	if (rc_data==NULL) {
	  rc_init(t);
	}
	rc_data->addTranscript(*t);
 //check this read alignment against ref exons and introns
}

struct COvlSorter {
  bool operator() (pair<int, const RC_Feature*> i,
      pair<int, const RC_Feature*> j) {
    return (i.first>j.first); //sort in decreasing order of overlap length
  }
} OvlSorter;


void rc_updateExonCounts(const RC_Feature* exon, int nh) {
  exon->rcount++;
  exon->mrcount += (nh > 1) ? (1.0/nh) : 1;
  if (nh<=1)  exon->ucount++;
}

bool BundleData::rc_count_hit(GBamRecord& brec, char strand, int nh) { //, int hi) {
 if (rc_data==NULL) return false; //no ref transcripts available for this reads' region
 if (rc_data->tdata.Count()==0) return false; //nothing to do without transcripts

 //check this read alignment against ref exons and introns
 /*
 int gstart=brec.start; //alignment start position on the genome

 We NEVER update read-counting boundaries of the bundle based on reads - they should only be based on reference transcripts
 if (rc_data->f_cov.size()==0 && gstart<rc_data->lmin) {
   fprintf(stderr, "Warning: adjusting lmin coverage bundle from %d to %d !\n", int(rc_data->lmin), (int)gstart);
   rc_data->lmin=gstart;
 }
 */
 if ((int)brec.end<rc_data->lmin || (int)brec.start>rc_data->rmax) {
	 return false; //hit outside coverage area
 }
 /*
 int gpos=brec.start; //current genomic position
 int rlen=0; //read length, obtained here from the cigar string
 int segstart=gstart;
 */
 vector<RC_Seg> rsegs;
 vector<RC_Seg> rintrons;
 for (int i=0;i<brec.exons.Count();i++) {
	 rc_data->updateCov(strand, nh, brec.exons[i].start, brec.exons[i].len());
	 rsegs.push_back(RC_Seg(brec.exons[i].start, brec.exons[i].end) );
	 if (i>0) {
		 //add intron
		 rintrons.push_back(RC_Seg(brec.exons[i-1].end+1, brec.exons[i].start-1));
	 }
 }
 //now check rexons and rintrons with findExons() and findIntron()
 for (size_t i=0;i<rintrons.size();++i) {
   /*RC_FeatIt ri=rc_data->findIntron(rintrons[i].l, rintrons[i].r, strand);
   if (ri!=rc_data->introns.end()) {
     (*ri).rcount++;
     (*ri).mrcount += (nh > 1) ? (1.0/nh) : 1;
     if (nh==1)  (*ri).ucount++;
   }
   */
	RC_Feature* ri=rc_data->findIntron(rintrons[i].l, rintrons[i].r, strand);
	if (ri) {
	    ri->rcount++;
	    ri->mrcount += (nh > 1) ? (1.0/nh) : 1;
	    if (nh==1)  ri->ucount++;
	}
 } //for each intron

 for (size_t i=0;i<rsegs.size();++i) {
     RC_FeatPtrSet ovlex=rc_data->findExons(rsegs[i].l, rsegs[i].r, strand);
     if (ovlex.size()==0) continue;
     if (ovlex.size()>1) {
       vector< pair<int, const RC_Feature*> > xovl; //overlapped exons sorted by decreasing overlap length
       for (RC_FeatPtrSet::iterator ox=ovlex.begin();ox != ovlex.end(); ++ox) {
         int ovlen=(*ox)->ovlen(rsegs[i].l, rsegs[i].r);
         if (ovlen>=5)
           xovl.push_back(pair<int, const RC_Feature*>(ovlen, *ox));
       }
       if (xovl.size()>1) {
          sort(xovl.begin(), xovl.end(), OvlSorter); //larger overlaps first
          //update the counts only for ref exons with max overlap to this segment
          int max_ovl=xovl.begin()->first;
          for (vector<pair<int, const RC_Feature*> >::iterator xo=xovl.begin();xo!=xovl.end();++xo) {
             if (max_ovl - xo->first > 5 ) break; //more than +5 bases coverage for the other exons
             rc_updateExonCounts(xo->second, nh);
          }
       } else if (xovl.size() == 1) {
         rc_updateExonCounts(xovl.begin()->second, nh);
       }
     } else {
       // 1 exon overlap only
       int ovlen=(*ovlex.begin())->ovlen(rsegs[i].l, rsegs[i].r);
       if (ovlen>=5) rc_updateExonCounts(*ovlex.begin(), nh);
     }
   } //for each read "exon"
 return true;
}

FILE* rc_fwopen(const char* fname) {
 if (strcmp(fname,"-")==0) return stdout;
 GStr fpath(ballgown_dir); //always ends with '/'
 //fpath += '.';
 fpath += fname;
 //fpath += "/";
 //fpath += fname;
 fpath += ".ctab";
 FILE* fh=fopen(fpath.chars(), "w");
 if (fh==NULL) {
   fprintf(stderr, "Error: cannot create file %s\n",
					fpath.chars());
   exit(1);
   }
 return fh;
}

FILE* rc_frenopen(const char* fname) {
	//if (strcmp(fname,"-")==0) return stdout;
	GStr fpath(ballgown_dir);
	//fpath += '.';
	fpath += fname;
	//fpath += "/";
	//fpath += fname;
	fpath += ".ctab";
	GStr fren(fpath);
	fren += ".tmp";
	//rename fpath to fren and open fren for reading
	if (rename(fpath.chars(), fren.chars())!=0) {
		 GError("Error: cannot rename %s to %s!\n", fpath.chars(), fren.chars());
	}
	FILE* fh=fopen(fren.chars(), "r");
	if (fh==NULL) {
	  GError("Error: cannot open file %s\n", fren.chars());
	}
	return fh;
}

void rc_frendel(const char* fname) {
	GStr fpath(ballgown_dir);
	//fpath += '.';
	fpath += fname;
	fpath += ".ctab";
	fpath += ".tmp";
	if (remove(fpath.chars())!=0) {
		GMessage("Warning: could not remove file %s!\n",fpath.chars());
	}
}

void rc_write_f2t(FILE* fh, map<uint, set<uint> >& f2t) {
  for (map<uint, set<uint> >::iterator m=f2t.begin(); m!=f2t.end(); ++m) {
    uint f_id=(*m).first;
    set<uint>& tset = (*m).second;
    for (set<uint>::iterator it=tset.begin();it!=tset.end();++it) {
 	 uint t_id = *it;
 	 fprintf(fh, "%u\t%u\n", f_id, t_id);
    }
  }
  fflush(fh);
}


void rc_update_exons(RC_BundleData& rc) {
	//update stdev etc. for all exons in bundle
	for (int f=0;f<rc.exons.Count(); ++f) {
      RC_Feature& exon = *(rc.exons[f]);
      //assert( exon.l >= rc.lmin );
      int L=exon.l-rc.lmin;
      int xlen=exon.r-exon.l+1;
      if (exon.l < rc.lmin) {
    	  //shouldn't be here, bundle read-counting boundaries should be based on exons!
    	  if (exon.r<rc.lmin) continue;
    	  xlen-=(rc.lmin-exon.l+1);
    	  L=0;
      }
      if (rc.rmax<exon.r) {
    	  if (exon.l>rc.rmax) continue; //should never happen
    	  xlen-=(exon.r-rc.rmax+1);
      }
      int R=L+xlen;
      vector<int>::iterator xcov_begin;
      vector<int>::iterator xcov_end;
      vector<float>::iterator xmcov_begin;
      vector<float>::iterator xmcov_end;
      if (exon.strand=='+' || exon.strand=='.') {
         xcov_begin  = rc.f_cov.begin()+L;
         xcov_end = rc.f_cov.begin()+R;
         xmcov_begin = rc.f_mcov.begin()+L;
         xmcov_end = rc.f_mcov.begin()+R;
      } else {
        xcov_begin  = rc.r_cov.begin()+L;
        xcov_end = rc.r_cov.begin()+R;
        xmcov_begin = rc.r_mcov.begin()+L;
        xmcov_end = rc.r_mcov.begin()+R;
      }

      double avg = (double)accumulate(xcov_begin, xcov_end, 0) / xlen;
      vector<double> diff(xlen);
      transform(xcov_begin, xcov_end, diff.begin(),
                     bind2nd( minus<double>(), avg));
      double sq_sum = inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
      double stdev = sqrt(sq_sum / xlen);

      double mavg = (double)accumulate(xmcov_begin, xmcov_end, 0) / xlen;
      vector<double> mdiff(xlen);
      transform(xmcov_begin, xmcov_end, mdiff.begin(),
                     bind2nd( minus<double>(), mavg));
      sq_sum = inner_product(mdiff.begin(), mdiff.end(), mdiff.begin(), 0.0);
      double mstdev = sqrt(sq_sum / xlen);
      exon.avg=avg;
      exon.stdev=stdev;
      exon.mavg=mavg;
      exon.mstdev=mstdev;
	} //for each exon in bundle


}


void rc_write_RCfeature( GPVec<RC_ScaffData>& rcdata, GPVec<RC_Feature>& features, FILE*& fdata, FILE*& f2t,
		               bool is_exon=false) {
  for (int i=0;i<features.Count();++i) {
	RC_Feature& f=*(features[i]);
	const char* ref_name=rcdata[f.t_id-1]->scaff->getGSeqName();
	if (is_exon) {
	  fprintf(fdata, "%u\t%s\t%c\t%d\t%d\t%d\t%d\t%.2f\t%.4f\t%.4f\t%.4f\t%.4f\n",
		  f.id, ref_name, f.strand, f.l, f.r, f.rcount,
		  f.ucount, f.mrcount, f.avg, f.stdev, f.mavg, f.mstdev);
	}
    else { //introns
	  fprintf(fdata,"%u\t%s\t%c\t%d\t%d\t%d\t%d\t%.2f\n",f.id, ref_name,
		  f.strand, f.l, f.r, f.rcount, f.ucount, f.mrcount);
	}
  // f2t -------
	fprintf(f2t, "%u\t%u\n", f.id, f.t_id);
  } //for each feature
  fclose(fdata);
  fclose(f2t);
}


/*void rc_write_counts(const char* refname, BundleData& bundle) {
 RC_BundleData& rc = *bundle.rc_data;
 //if (rc.exons.size()==0) return;
  *
*/
void rc_writeRC(GPVec<RC_ScaffData>& RC_data,
		GPVec<RC_Feature>& RC_exons,
		GPVec<RC_Feature>& RC_introns,
		FILE* &f_tdata, FILE* &f_edata, FILE* &f_idata,
        FILE* &f_e2t, FILE* &f_i2t) {

 for (int t=0;t<RC_data.Count();++t) {
  //File: t_data.ctab
  //t_id tname chr strand start end num_exons gene_id gene_name cufflinks_cov cufflinks_fpkm
   const RC_ScaffData& sd=*RC_data[t];
   const char* refname = sd.scaff->getGSeqName();
   const char* genename= sd.scaff->getGeneName();
   if (genename==NULL) genename=".";
   fprintf(f_tdata, "%u\t%s\t%c\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%f\t%f\n",
	  sd.t_id, refname, sd.strand, sd.l, sd.r, sd.scaff->getID(),
	  sd.num_exons, sd.eff_len, sd.scaff->getGeneID(),
	  genename, sd.cov, sd.fpkm);
 }//for each transcript
 //fflush(f_tdata);
 fclose(f_tdata);

 //File: e_data.ctab
 //e_id chr gstart gend rcount ucount mrcount
 rc_write_RCfeature(RC_data, RC_exons, f_edata, f_e2t, true);
 //File: i_data.ctab
 //i_id chr gstart gend rcount ucount mrcount
 rc_write_RCfeature(RC_data, RC_introns, f_idata, f_i2t);

 // feature-to-transcript link files
 //rc_write_f2t(rc.fe2t,  rc.e2t);
 //rc_write_f2t(rc.fi2t,  rc.i2t);
}

void RC_ScaffData::addFeature(int fl, int fr, GPVec<RC_Feature>& fvec,
                          uint& f_id, set<RC_ScaffSeg>& fset, set<RC_ScaffSeg>::iterator& fit,
						  GPVec<RC_Feature>& fdata) {
  //f_id is the largest f_id inserted so far in fset
  bool add_new = true;
  RC_ScaffSeg newseg(fl,fr,this->strand);
  //RC_Feature* newfeature=NULL;
  if (fset.size()>0) {
    if (fit == fset.end()) --fit;
    if (newseg < (*fit)) {
      bool eq=false;
      while (newseg < (*fit) || (eq = (newseg==(*fit)))) {
        if (eq) {
          add_new = false;
          newseg.id = fit->id;
          break;
        }
        if (fit==fset.begin()) {
          break;
        }
        --fit;
      }
    }
    else { //newseg >= *fit
      bool eq=false;
      while ((*fit) < newseg || (eq = (newseg==(*fit)))) {
        if (eq) {
          add_new = false;
          newseg.id = fit->id;
          break;
        }
        ++fit;
        if (fit==fset.end()) {
          --fit;
          break;
        }
      }
    }
  }
  if (add_new) {
    newseg.id = ++f_id;
    pair< set<RC_ScaffSeg>::iterator, bool> ret = fset.insert(newseg);
    //ret.second = was_inserted (indeed new)
    if (!ret.second) {
      GError("Error: feature %d-%d (%c) already in segment set!\n",
                 newseg.l, newseg.r, newseg.strand);
      //newseg.id = ret.first->id;
    }
    fdata.Add(new RC_Feature(newseg, this->t_id));
#ifdef DEBUG
    if (fdata.Count()!=(int)f_id) {
    	GMessage("Error: fdata.Count=%d, f_id=%d\n", fdata.Count(), f_id);
    }
#endif
    GASSERT((uint)fdata.Count()==f_id);
  }
  //fvec.push_back(newseg);
  GASSERT(fdata[newseg.id-1]->id==newseg.id);
  fvec.Add(fdata[newseg.id-1]);
}

void Ballgown_setupFiles(FILE* &f_tdata, FILE* &f_edata, FILE* &f_idata,
	            FILE* &f_e2t, FILE* &f_i2t) {
  if (f_tdata == NULL) {
	//first call, create the files
	 f_tdata = rc_fwopen("t_data");
	 fprintf(f_tdata, "t_id\tchr\tstrand\tstart\tend\tt_name\tnum_exons\tlength\tgene_id\tgene_name\tcov\tFPKM\n");
	 f_edata = rc_fwopen("e_data");
    fprintf(f_edata, "e_id\tchr\tstrand\tstart\tend\trcount\tucount\tmrcount\tcov\tcov_sd\tmcov\tmcov_sd\n");
	 f_idata = rc_fwopen("i_data");
    fprintf(f_idata, "i_id\tchr\tstrand\tstart\tend\trcount\tucount\tmrcount\n");
    f_e2t = rc_fwopen("e2t");
    fprintf(f_e2t,  "e_id\tt_id\n");
    f_i2t = rc_fwopen("i2t");
    fprintf(f_i2t,  "i_id\tt_id\n");
  }
}

void RC_ScaffData::rc_addFeatures(uint& c_e_id, set<RC_ScaffSeg>& fexons, GPVec<RC_Feature>& edata,
                    uint& c_i_id, set<RC_ScaffSeg>& fintrons, GPVec<RC_Feature>& idata) {
  GASSERT(scaff);
  GffObj& m = *(scaff);
  for (int i = 0; i < m.exons.Count(); ++i)  {
    set<RC_ScaffSeg>::iterator eit=fexons.end();
    set<RC_ScaffSeg>::iterator iit=fintrons.end();
    addFeature(m.exons[i]->start, m.exons[i]->end, t_exons, c_e_id, fexons, eit, edata);
    if (i>0) { //store intron
      addFeature(m.exons[i-1]->end+1, m.exons[i]->start-1, t_introns, c_i_id, fintrons, iit, idata);
    }
  } //for each exon
}
