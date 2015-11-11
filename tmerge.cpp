#include "tmerge.h"


void TInputFiles::Add(const char* fn) {
		if (fileExists(fn)<2) {
			    GError("Error: input file %s cannot be found!\n",
			            fn);
		}
		files.Add(new GStr(fn));
	}

GBamRecord* TInputFiles::next() {
	//must free old current record first
	delete crec;
	crec=NULL;
    if (recs.Count()>0) {
    	crec=recs.Pop();//lowest coordinate
    	GBamRecord* rnext=readers[crec->fidx]->next();
    	if (rnext)
    		recs.Add(new TInputRecord(rnext,crec->fidx));
    	return crec->brec;
    }
    else return NULL;
}


int gseqstat_cmpName(const pointer p1, const pointer p2) {
	GSeqStat& g1=*((GffObj*)p1);
	GSeqStat& g2=*((GffObj*)p2);
	return strcmp(g1.gseqname, g2.gseqname);
}


GStr convert2BAM(GStr* gtf) {
  GStr bamfname(tmp_path);
  bamfname+=*gtf;
  GStr samhname(bamfname);
  bamfname+=".bam";
  samhname+=".sam";
  FILE* samh=fopen(samhname.chars(), "w");
  if (samh==NULL) GError("Error creating file: %s\n",samhname.chars());
  fprintf(samh, "@HD\tVN:1.0\tSO:coordinate\n");
  //load GTF as sorted
  GffReader gfr(gtf->chars(), true, true); //transcript only, sorted by location
  gfr.showWarnings(debugMode || verbose);
  gfr.readAll(true, true, true); //keep attributes, merge close exons, no_exon_attributes
  gfr.gseqStats.Sort(gseqstat_cmpName);
  for (int i=0;i<gfr.gseqStats.Count();++i) {
  	fprintf(samh, "@SQ\tSN:%s\tLN:%ul\n", gfr.gseqStats[i]->gseqname,
  			gfr.gseqStats[i]->maxcoord+500);
  }
  fclose(samh);
  GBamWriter bw(bamfname.chars(),samhname.chars());
  for (int i=0;i<gfr.gflst.Count();++i) {
	  GffObj& m = *gfr.gflst[i];
	  int t_id=bw.get_tid(m.getGSeqName());
	  if (t_id<0)
		   GError("Error getting header ID# for gseq %s (file: %s)\n",m.getGSeqName(),gtf->chars());
	  GStr cigar;
	  for (int k=0;k<m.exons.Count();++k) {
		  if (k>0) {
			  cigar+=int(m.exons[k]->start-m.exons[k-1]->end-1); cigar+='N';
		  }
		  cigar+=m.exons[k]->len();
		  cigar+='M';
	  }
	  GBamRecord brec(m.getID(), t_id, m.start, false, "*", cigar.chars());
	  if (m.strand=='-' || m.strand=='+') {
		   GStr tag("XS:A:");
		   tag+=m.strand;
	       brec.add_aux(tag.chars());
	  }
	  //TODO: add special int tag with file index value (qry#)
	  //TODO: add cov, FPKM, TPM as special float tags
  }
}


int TInputFiles::start() {
	GPVec<GStr> bamfiles(this->files);
	if (mergeMode) {//convert to temp BAM files
		bamfiles.Clear();
		for (int i=0;i<files.Count();++i) {
			bamfiles.Add(new GStr(convert2BAM(files[i])));
		}
	}
	//regular stringtie BAM input
	for (int i=0;i<bamfiles.Count();++i) {
		readers.Add(new GBamReader(bamfiles[i]->chars()));
	}
	return readers.Count();
}


