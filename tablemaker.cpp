#include "rlink.h" // this includes "tablemaker.h", "bundle.h" etc.
#include <numeric>

extern GStr ballgown_dir;

int rc_cov_inc(int i) {
  return ++i;
}

struct COvlSorter {
  bool operator() (pair<int, const RC_Feature*> i,
      pair<int, const RC_Feature*> j) {
    return (i.first>j.first); //sort in decreasing order of overlap length
  }
} OvlSorter;


void rc_updateExonCounts(const RC_ExonOvl& exonovl, int nh) {
  //this only gets read overlaps > 5bp or otherwise filtered in evalReadAln()
  exonovl.feature->rcount++;
  if (nh>1) {
	  exonovl.feature->mrcount += (1.0/nh);
	  exonovl.feature->movlcount += ((double)exonovl.ovlen/nh);
  }
  else { // nh<=1
	  exonovl.feature->mrcount++;
	  exonovl.feature->movlcount += exonovl.ovlen;
	  exonovl.feature->ucount++;
  }
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
	for (int f=0;f<rc.g_exons.Count(); ++f) {
      RC_Feature& exon = *(rc.g_exons[f]);
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


void rc_write_RCfeature( GPVec<RC_TData>& rcdata, GPVec<RC_Feature>& features, FILE*& fdata, FILE*& f2t,
		               bool is_exon=false) {
  for (int i=0;i<features.Count();++i) {
	RC_Feature& f=*(features[i]);
	const char* ref_name=rcdata[f.t_ids[0]-1]->ref_t->getGSeqName();
	if (is_exon) {
	  fprintf(fdata, "%u\t%s\t%c\t%ld\t%ld\t%d\t%d\t%.2f\t%.4f\t%.4f\t%.4f\t%.4f\n",
		  f.id, ref_name, f.strand, f.l, f.r, f.rcount, 
		  f.ucount, f.mrcount, f.avg, f.stdev, f.mavg, f.mstdev);
	}
    else { //introns
	  fprintf(fdata,"%u\t%s\t%c\t%ld\t%ld\t%d\t%d\t%.2f\n",f.id, ref_name,
		  f.strand, f.l, f.r, f.rcount, f.ucount, f.mrcount);
	}
  // f2t -------
	for (int t=0;t<f.t_ids.Count();++t)
	   fprintf(f2t, "%u\t%u\n", f.id, f.t_ids[t]);
  } //for each feature
  fclose(fdata);
  fclose(f2t);
}


void rc_writeRC(GPVec<RC_TData>& RC_data,
		GPVec<RC_Feature>& RC_exons,
		GPVec<RC_Feature>& RC_introns,
		FILE* &f_tdata, FILE* &f_edata, FILE* &f_idata,
        FILE* &f_e2t, FILE* &f_i2t) {

 for (int t=0;t<RC_data.Count();++t) {
  //File: t_data.ctab
  //t_id tname chr strand start end num_exons gene_id gene_name cufflinks_cov cufflinks_fpkm
  RC_TData& sd=*RC_data[t];
   const char* refname = sd.ref_t->getGSeqName();
   const char* genename= sd.ref_t->getGeneName();
   const char* geneID= sd.ref_t->getGeneID();
   if (genename==NULL) genename=".";
   if (geneID==NULL) geneID=".";

   fprintf(f_tdata, "%u\t%s\t%c\t%ld\t%ld\t%s\t%d\t%d\t%s\t%s\t%f\t%f\n",
	  sd.t_id, refname, sd.ref_t->strand, sd.l, sd.r, sd.ref_t->getID(),
	  sd.t_exons.Count(), sd.eff_len, geneID,
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
}
// collect and prepare exons and introns for their raw counts storage etc. 
void RC_TData::rc_addFeatures(uint& c_e_id, GList<RC_Feature>& exonSet, GPVec<RC_Feature>& exonTable,
                    uint& c_i_id, GList<RC_Feature>& intronSet, GPVec<RC_Feature>& intronTable) {
  GASSERT(ref_t);
  GffObj& m = *(ref_t);
  int ecache_idx = exonSet.Count()-1;
  int icache_idx = intronSet.Count()-1;
  //int ecache_idx = e_idx_cache>=0 ? e_idx_cache : exonSet.Count()-1;
  //int icache_idx = i_idx_cache>=0 ? i_idx_cache : intronSet.Count()-1;
  for (int i = 0; i < m.exons.Count(); ++i)  {
    addFeature((long)m.exons[i]->start, (long)m.exons[i]->end, t_exons, c_e_id, exonSet, exonTable, ecache_idx);
    //if (i==0) e_idx_cache=ecache_idx;
    if (i>0) { //store intron
      //if (i==1) i_idx_cache=icache_idx;
      addFeature(long(m.exons[i-1]->end)+1, long(m.exons[i]->start)-1, t_introns, c_i_id,
    		  intronSet, intronTable, icache_idx);
    } //for each intron
  } //for each exon
}

void RC_TData::addFeature(long fl, long fr, GPVec<RC_Feature>& fvec,
                          uint& f_id, GList<RC_Feature>& fset,
						  GPVec<RC_Feature>& fdata, int& cache_idx) {
  //f_id is the largest f_id inserted so far in fset
  bool add_new = true;
  bool notNascent = !isNascent(ref_t);
  RC_Feature* newseg=new RC_Feature(fl, fr, ref_t->strand, 0, this->t_id);
  //RC_Feature* newfeature=NULL;
  int fit=cache_idx<0 ? fset.Count()-1 : cache_idx;
  int fp_id=-1;
  if (fset.Count()>0) {
    if (*newseg < *(fset[fit])) {
      bool eq=false;
      while (*newseg < *(fset[fit]) || (eq = (*newseg==*(fset[fit])))) {
        if (eq) {
          add_new = false;
          fp_id = fset[fit]->id; 
          if (notNascent) fset[fit]->nascent_only = false;
          break;
        }
        //newseg< fset[fit]
        --fit;
        if (fit<0) break; //newseg should be inserted at 0
      } //while newseg<fset[fit]
      if (add_new) ++fit;
         // newseg < fset[fit+1]
         //we'll insert newseg at position fit+1
    }
    else { //newseg >= *fset[fit]
      bool eq=false;
      while (*(fset[fit]) < *newseg || (eq = (*newseg==*(fset[fit])))) {
        if (eq) {
          add_new = false;
          fp_id = fset[fit]->id;
          if (notNascent) fset[fit]->nascent_only = false;
          break;
        }
        ++fit;
        if (fit==fset.Count()) {
        	//newseg should be appended to the list
        	break;
        }
      }
      //if (fit<=fset.Count() && !add_new) {
         //fset[fit-1] < newseg  < fset[fit]
         //we'll insert newseg at position fit
      //}
    }
  } //check existing set
  if (add_new) { //did not see this feature before
    newseg->id = ++f_id;
    if (notNascent) newseg->nascent_only = false;
    if (fit<0) fit = fset.Add(newseg);
    else fset.sortInsert(fit, newseg);
    if (fit<0) {
      GError("Error: feature %d-%d (%c) already in feature set!\n",
                 newseg->l, newseg->r, newseg->strand);
    }
    cache_idx=fit;
    fp_id=fdata.Add(newseg)+1;

#ifdef DEBUG
    if (fdata.Count()!=(int)f_id) {
    	GMessage("Error: fdata.Count=%d, f_id=%d\n", fdata.Count(), f_id);
    }
#endif
    GASSERT((uint)fdata.Count()==f_id);
  }
  else { //feature seen before, only update its parent list
    fdata[fp_id-1]->t_ids.Add(this->t_id);
    delete newseg;
  }
  GASSERT(fdata[fp_id-1]->id==(uint)fp_id);
  // always add the feature to the transcript's t_exons/t_introns vector  
  fvec.Add(fdata[fp_id-1]);
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
