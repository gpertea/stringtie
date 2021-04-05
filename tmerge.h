#ifndef STRINGTIE_MERGE_H_
#define STRINGTIE_MERGE_H_
#include "GStr.h"
#include "GList.hh"
#include "rlink.h"
#include "GIntervalTree.hh"

extern GStr tmp_path;
extern bool keepTempFiles;
extern GStr mfltgff;

struct GSTree {
	GIntervalTree it[3]; //0=unstranded, 1: +strand, 2: -strand
};

extern GHash<GSTree*> map_trees; //map a ref sequence name to its own interval trees (3 per ref seq)

int loadITree(const char* fname); //return number of transcripts loaded from fname into map_trees
bool matchITree(GffObj& t); //return true if t has a match in map_trees

struct TInputRecord {
	GBamRecord* brec;
	int fidx; //index in files and readers
	bool operator<(TInputRecord& o) {
		 //decreasing location sort
		 GBamRecord& r1=*brec;
		 GBamRecord& r2=*(o.brec);
		 //int refcmp=strcmp(r1.refName(),r2.refName());
		 int refcmp=mergeMode ? strcmp(r1.refName(),r2.refName()) : r1.refId()-r2.refId();
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
		 else { //use header order
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
	bool mFlt;
	GPVec<GBamReader> readers;
	GVec<GStr> files; //same order
	GVec<GStr> tmpfiles; //all the temp files created by this
	GList<TInputRecord> recs; //next record for each
	TInputFiles():crec(NULL), mFlt(false), readers(true), files(), tmpfiles(),
			recs(true, true, true) {
		if (!mfltgff.is_empty()) {
			mFlt=(loadITree(mfltgff.chars())>0);
		}
	}
	void Add(const char* fn);
	int count() { return files.Count(); }
	int start(); //open all files, load 1 record from each
	GBamRecord* next();
	void stop(); //
};

#endif /* STRINGTIE_MERGE_H_ */
