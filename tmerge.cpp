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

char* convert2BAM(GStr* gtf) {

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


