#ifndef STRINGTIE_MERGE_H_
#define STRINGTIE_MERGE_H_
#include "GStr.h"
#include "GList.hh"
#include "rlink.h"
extern GStr tmp_path;
extern bool keepTempFiles;
/*struct TInFile {
	GStr fpath;
	int ftype; //0=bam file, 1=GTF transfrags
	TInFile(const char* fn=NULL,int ft=0):fpath(fn),ftype(ft) { }
}*/
struct TInputRecord {
	GBamRecord* brec;
	int fidx; //index in files and readers
	bool operator<(TInputRecord& o) {
		 //decreasing location sort
		 GBamRecord& r1=*brec;
		 GBamRecord& r2=*(o.brec);
		 int refcmp=strcmp(r1.refName(),r2.refName());
		 if (refcmp==0) {
		 //higher coords first
			if (r1.start!=r2.start)
				 return (r1.start>r2.start);
			else {
				if (r1.end!=r2.end)
				   return (r1.end>r2.end);
				else if (fidx==o.fidx)
						return (strcmp(r1.name(), r2.name())>0);
					else return fidx>o.fidx;
			}
		 }
		 else { //use lexicographic order of ref seqs
			 return (refcmp>0);
		 }
	}
	bool operator==(TInputRecord& o) {
		 GBamRecord& r1=*brec;
		 GBamRecord& r2=*(o.brec);
		 return ( strcmp(r1.refName(),r2.refName())==0 && r1.start==r2.start && r1.end==r2.end
				 && fidx==o.fidx && strcmp(r1.name(),r2.name())==0);
	}

	TInputRecord(GBamRecord* b=NULL, int i=0):brec(b),fidx(i) {}
	~TInputRecord() {
		delete brec;
	}
};

struct TInputFiles {
 protected:
	TInputRecord* crec;
	GStr convert2BAM(GStr& gtf, int idx);
 public:
	GPVec<GBamReader> readers;
	GVec<GStr> files; //same order
	GVec<GStr> tmpfiles; //all the temp files created by this
	GList<TInputRecord> recs; //next record for each
	TInputFiles():crec(NULL), readers(true), files(), tmpfiles(),
			recs(true, true, true) { }
	void Add(const char* fn);
	int count() { return files.Count(); }
	int start(); //open all files, load 1 record from each
	GBamRecord* next();
	void stop(); //
};


#endif /* STRINGTIE_MERGE_H_ */
