#include "tmerge.h"


void TInputFiles::Add(const char* fn) {
	   GStr sfn(fn);
		if (sfn!="-" && fileExists(fn)<2) {
			    GError("Error: input file %s cannot be found!\n",
			            fn);
		}

		files.Add(sfn);
	}


int gseqstat_cmpName(const pointer p1, const pointer p2) {
	return strcmp(((GSeqStat*)p1)->gseqname, ((GSeqStat*)p2)->gseqname);
}

GStr TInputFiles::convert2BAM(GStr& gtf, int idx) {
  GStr bamfname(tmp_path);
  bamfname.appendfmt("transcripts_s%04d",idx);
  GStr samhname(bamfname);
  bamfname+=".bam";
  samhname+=".sam";
  tmpfiles.Add(bamfname);
  tmpfiles.Add(samhname);
  FILE* samh=fopen(samhname.chars(), "w");
  if (samh==NULL) GError("Error creating file: %s\n",samhname.chars());
  fprintf(samh, "@HD\tVN:1.0\tSO:coordinate\n");
  //load GTF as sorted
  GffReader gfr(gtf.chars(), true, true); //transcript only, sorted by location
  gfr.setRefAlphaSorted(true); //make sure refseq IDs are sorted alphabetically
  gfr.showWarnings(debugMode || verbose);
  gfr.readAll(true, true, true); //keep attributes, merge close exons, no_exon_attributes
  if (gfr.gflst.Count()==0)
	  GError("Error: no transcripts were found in input file %s\n", gtf.chars());
  gfr.gseqStats.Sort(gseqstat_cmpName);
  for (int i=0;i<gfr.gseqStats.Count();++i) {
  	fprintf(samh, "@SQ\tSN:%s\tLN:%u\n", gfr.gseqStats[i]->gseqname,
  			gfr.gseqStats[i]->maxcoord+500);
  }
  fprintf(samh, "@CO\tfn:%s\n",gtf.chars());
  fclose(samh);
  GBamWriter bw(bamfname.chars(),samhname.chars());
  for (int i=0;i<gfr.gflst.Count();++i) {
	  GffObj& m = *gfr.gflst[i];
	  int t_id=bw.get_tid(m.getGSeqName());
	  if (t_id<0)
		   GError("Error getting header ID# for gseq %s (file: %s)\n",m.getGSeqName(),gtf.chars());
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
	  GStr s("ZF:i:");
	  s+=idx;
	  brec.add_aux(s.chars());
	  char *av=m.getAttr("cov");
	  if (av!=NULL) {
		  s="ZS:Z:";
		  s+=av;
		  s+='|';
		  av=m.getAttr("FPKM");
		  if (av) s+=av;
		  av=m.getAttr("TPM");
		  if (av) { s+='|';s+=av; }
		  brec.add_aux(s.chars());
	  }
	  bw.write(&brec);
  } //for each transcript
  return bamfname;
}


int TInputFiles::start() {
	GVec<GStr> bamfiles;
	if (mergeMode && this->files.Count()==1) {
		//special case, if it's only one file it must be a list (usually)
		GStr fname(this->files.First());
		FILE* flst=fopen(this->files.First().chars(),"r");
		if (flst==NULL) GError("Error: could not open input file %s!\n",
				fname.chars());
		GVec<GStr> infiles;
		char* line=NULL;
		int lcap=5000;
		GMALLOC(line, lcap);
		bool firstline=true;
		//bool isalist=true;
		while (fgetline(line,lcap,flst)) {
			GStr s(line);
			s.trim();
			if (s.length()<2 || s[0]=='#') continue; //skip comments/header in the list file, if any
			if (firstline) {
				if (!fileExists(s.chars())) {
					//it must be a GFF
					//this shouldn't normally happen in mergeMode, with one input file
					if (s.count('\t')<8)
						GError("Error: cannot find file '%s' and %s does not look like GFF!\n", s.chars(), fname.chars());
					break;
				}
				firstline=false;
				files.Clear();
			}
			if (!fileExists(s.chars()))
				GError("Error opening transcript file %s !\n",s.chars());
			fname=s;
			files.Add(fname);
		} //for each line in the list file
		GFREE(line);
		fclose(flst);
	}
	if (mergeMode) {//files are GTF/GFF, convert to temp BAM files
		for (int i=0;i<files.Count();++i) {
			//crude way to bypass GTF conversion when resuming/debugging
			if (files[i].endsWith(".bam")) {
				bamfiles.Add(files[i]);
			}
			else {
				GStr s=convert2BAM(files[i], i);
				bamfiles.Add(s);
			}
		}
	}
	else {
		bamfiles=files;
	}
	//stringtie multi-BAM input
	for (int i=0;i<bamfiles.Count();++i) {
		GBamReader* bamreader=new GBamReader(bamfiles[i].chars(),
				(bamfiles[i]=="-" && forceBAM));
		readers.Add(bamreader);
		GBamRecord* brec=bamreader->next();
		if (brec)
		   recs.Add(new TInputRecord(brec, i));
	}
	return readers.Count();
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

void TInputFiles::stop() {
 for (int i=0;i<readers.Count();++i) {
	 readers[i]->bclose();
 }
 if (!keepTempFiles) {
	 for (int i=0;i<tmpfiles.Count();++i) {
		 unlink(tmpfiles[i].chars());
	 }
 }
}

