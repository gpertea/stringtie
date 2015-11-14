#include "rlink.h"
#include "GBitVec.h"
#include <float.h>

//#define GMEMTRACE 1  //debugging mem allocation

#ifdef GMEMTRACE
#include "proc_mem.h"
#endif

//import globals from main program:

//extern GffNames* gseqNames;
extern FILE *c_out;         // file handle for the input transcripts that are fully covered by reads

extern bool specific;
extern bool trim;
//extern bool partialcov;
extern bool complete;
extern bool eonly;

//extern int maxReadCov;

extern float isofrac;
extern float mcov;
extern int mintranscriptlen; // minimum number for a transcript to be printed
extern int sensitivitylevel;
extern uint junctionsupport; // anchor length for junction to be considered well supported <- consider shorter??
extern int junctionthr; // number of reads needed to support a particular junction
extern float readthr;     // read coverage per bundle bp to accept it; otherwise considered noise
extern uint bundledist;  // reads at what distance should be considered part of separate bundles
                        // <- this is not addressed everywhere, e.g. in infer_transcripts -> look into this

extern bool includesource;
extern bool EM;
extern bool weight;
extern bool geneabundance; // need to comute the gene abundance

extern FILE* f_out;
extern GStr label;


void printTime(FILE* f) {
	time_t ltime; /* calendar time */
	ltime=time(NULL);
	struct tm *t=localtime(&ltime);
	fprintf(f, "[%02d/%02d %02d:%02d:%02d]",t->tm_mon+1, t->tm_mday,
			t->tm_hour, t->tm_min, t->tm_sec);
}

void printBitVec(GBitVec& bv) {
   for (uint i=0;i<bv.size();i++) {
       fprintf(stderr, "%c", bv.test(i)?'1':'0');
   }
}

/*
CJunction* add_junction(int start, int end, int leftsupport, int rightsupport,
		GList<CJunction>& junction, char strand, int nh) {
	int oidx=-1;
	CJunction *nj=NULL;
	//CJunction* nj=junction.AddIfNew(new CJunction(start, end, strand), true, &oidx);
	CJunction jn(start, end, strand);
	if (junction.Found(&jn, oidx)) {
		nj=junction.Get(oidx);
	}
	else {
		nj=new CJunction(start, end, strand);
		junction.Add(nj);
	}
	//if (nh==0) nh=1;
	//nj->nreads+=float(1)/nh;
	if (leftsupport >= junctionsupport && rightsupport >=junctionsupport) {
		nj->nreads_good+=float(1)/nh;
	}
	return nj;
}
*/

CJunction* add_junction(int start, int end, GList<CJunction>& junction, char strand) {
	int oidx=-1;
	CJunction *nj=NULL;
	//CJunction* nj=junction.AddIfNew(new CJunction(start, end, strand), true, &oidx);
	CJunction jn(start, end, strand);
	if (junction.Found(&jn, oidx)) {
		nj=junction.Get(oidx);
	}
	else {
		nj=new CJunction(start, end, strand);
		junction.Add(nj);
	}

	return nj;
}

void cov_add(GVec<float>& bpcov, int i, int j, float v) {
	if (j>=bpcov.Count()) bpcov.Resize(j+1, 0);
	for (int k=i;k<j;k++)
		bpcov[k]+=v;
}

void cov_add(GVec<float>* bpcov, int sno, int i, int j, float v) {
	bool neutral=false;
	if(sno!=1) neutral=true;
	if (j>=bpcov[sno].Count())
		for(int s=0;s<3;s++) bpcov[s].Resize(j+1, 0);
	for (int k=i;k<j;k++) {
		bpcov[sno][k]+=v;
		if(neutral) bpcov[1][k]+=v;
	}
}

float getBCov(GVec<float>& bpcov, int p) {
	//if (p<0) GMessage("Error: invalid bpcov index (%d)!\n", p);
	if (p>=bpcov.Count()) return 0;
	else return bpcov[p];
}

/*
bool maxCovReached(int currentstart, GBamRecord& brec, BundleData& bdata) { // coverage saturation is reached if read spans no portion of saturation less than maxReadCov
	for (int i=0;i<brec.exons.Count();i++) {
		if (getBCov(bdata.bpcov, brec.exons[i].start-currentstart)<=maxReadCov)
			return false;
	}
	if (!bdata.covSaturated) {
		GMessage("Warning: bundle %s:%d-%d(%d) (%djs) reached coverage saturation (%d) starting with read mapped at %d\n",
				bdata.refseq.chars(), bdata.start, bdata.end, bdata.numreads, bdata.junction.Count(),
				maxReadCov, brec.start);
		bdata.covSaturated=true;
	}
	return true;
}
*/

void countFragment(BundleData& bdata, GBamRecord& brec, int nh) {
	static uint32_t BAM_R2SINGLE = BAM_FREAD2 | BAM_FMUNMAP ;


	for (int i=0;i<brec.exons.Count();i++) {
		bdata.frag_len+=float(1)*brec.exons[i].len()/nh;
	}
	if (!brec.isPaired() || ((brec.flags()&BAM_FREAD1)!=0) ||
					((brec.flags()&BAM_R2SINGLE)==BAM_R2SINGLE ) ) {
		bdata.num_fragments+=float(1)/nh;
	}

	/*
	if (hi==0) {
		for (int i=0;i<brec.exons.Count();i++) {
			bdata.frag_len+=brec.exons[i].len();
		}
		if (!brec.isPaired() || ((brec.flags()&BAM_FREAD1)!=0) ||
				((brec.flags()&BAM_R2SINGLE)==BAM_R2SINGLE ) ) {
			bdata.num_fragments++;
		}
	}
	else if (hi==1) {
		for (int i=0;i<brec.exons.Count();i++) {
			bdata.frag_len1+=brec.exons[i].len();
		}
		if (!brec.isPaired() || ((brec.flags()&BAM_FREAD1)!=0) ||
					((brec.flags()&BAM_R2SINGLE)==BAM_R2SINGLE ) ) {
			bdata.num_fragments1++;
		}
	}
	*/
}

bool exonmatch(GVec<GSeg> &prevexons, GVec<GSeg> &exons) {
	if(prevexons.Count() != exons.Count()) return false;
	for(int i=0;i<exons.Count();i++) {
		if(prevexons[i].end!=exons[i].end || prevexons[i].start!=exons[i].start) return false;
	}
	return true;
}

void processRead(int currentstart, int currentend, BundleData& bdata,
		 GHash<int>& hashread,  GReadAlnData& alndata) { // some false positives should be eliminated here in order to break the bundle

	GList<CReadAln>& readlist = bdata.readlist;    // list of reads gathered so far
	GList<CJunction>& junction = bdata.junction;   // junctions added so far
	GBamRecord& brec=*(alndata.brec);			   // bam record
    char strand=alndata.strand;
    int nh=alndata.nh;
    int hi=alndata.hi;
	int readstart=brec.start;
	CReadAln* readaln=NULL;                        // readaln is initialized with NULL
	//bool covSaturated=false;                       // coverage is set to not saturated

	bool match=false;  // true if current read matches a previous read
	int n=readlist.Count()-1;
	while(n>-1 && readlist[n]->start==brec.start) {
		if(strand==readlist[n]->strand) match=exonmatch(readlist[n]->segs,brec.exons);
		if(match) break; // this way I make sure that I keep the n of the matching readlist
		n--;
	}

	if (bdata.end<currentend) {
		bdata.start=currentstart;
		bdata.end=currentend;
	}
	bdata.numreads++;                         // number of reads gets increased no matter what
	//bdata.wnumreads+=float(1)/nh;

	if (!match) { // if this is a new read I am seeing I need to set it up
		readaln=new CReadAln(strand, nh, brec.start, brec.end, alndata.tinfo);
		for (int i=0;i<brec.exons.Count();i++) {
			readaln->len+=brec.exons[i].len();
			if(i) {
				CJunction* nj=add_junction(brec.exons[i-1].end, brec.exons[i].start, junction, strand);
				if (alndata.juncs.Count())
					nj->guide_match=alndata.juncs[i-1]->guide_match;
				if (nj) readaln->juncs.Add(nj);
			}
			readaln->segs.Add(brec.exons[i]);
		}
		n=readlist.Add(readaln);  // reset n for the case there is no match
	}
	else { //redundant read alignment matching a previous alignment
		// keep shortest nh so that I can see for each particular read the multi-hit proportion
		if(nh<readlist[n]->nh) readlist[n]->nh=nh;
		//for mergeMode, we have to free the transcript info:
		if (alndata.tinfo!=NULL) {
			 delete alndata.tinfo;
			 alndata.tinfo=NULL;
		}
	}

	if((int)brec.end>currentend) {
			currentend=brec.end;
	  	bdata.end=currentend;
	}

	float rdcount=float(1)/nh;
	readlist[n]->read_count+=rdcount; // increase single count just in case I don't see the pair

	// now set up the pairing
	if (brec.refId()==brec.mate_refId()) {  //only consider mate pairing data if mates are on the same chromosome/contig
		int pairstart=brec.mate_start();
		if (currentstart<=pairstart) { // if pairstart is in a previous bundle I don't care about it

			GStr readname(brec.name());
			GStr id(readname); // init id with readname
			//id+=':';id+=readstart; // this shouldn't be needed since I have the HI tag
			id+=':';id+=hi;

			if(pairstart<=readstart) { // if I've seen the pair already
				const int* np=hashread[id.chars()];
				if(np) { // the pair was stored --> why wouldn't it be? : only in the case that the pair starts at the same position
					bool notfound=true;
					for(int i=0;i<readlist[*np]->pair_idx.Count();i++)
						if(readlist[*np]->pair_idx[i]==n) {
							readlist[*np]->pair_count[i]+=rdcount;
							notfound=false;
							break;
						}
					if(notfound) { // I didn't see the pairing before
						readlist[*np]->pair_idx.Add(n);
						readlist[*np]->pair_count.Add(rdcount);
					}

					notfound=true;
					for(int i=0;i<readlist[n]->pair_idx.Count();i++)
						if(readlist[n]->pair_idx[i]==*np) {
							readlist[n]->pair_count[i]+=rdcount;
							notfound=false;
							break;
						}
					if(notfound) { // I didn't see the pairing before
						int i=*np;
						readlist[n]->pair_idx.Add(i);
						readlist[n]->pair_count.Add(rdcount);
					}
					hashread.Remove(id.chars());
				}
				/* do not allow pairings starting at the same place --> way too confusing and I don't gain much info from these
				else if(readstart==pairstart) { // I have not seen the pair because it starts at the same place -> need to store it
					hashread.Add(id.chars(), new int(n));
				}
				*/
			}
			else { // I might still see the pair in the future
				hashread.Add(id.chars(), new int(n));
			}
		}
	} //<-- if mate is mapped on the same chromosome
}

int get_min_start(CGroup **currgroup) {
	int nextgr=0;

	if(currgroup[0]!=NULL) {
		if(currgroup[1]!=NULL) {
		    if(currgroup[2]!=NULL) {
		    	int twogr = currgroup[0]->start < currgroup[1]->start ? 0 : 1;

		    	nextgr = currgroup[twogr]->start < currgroup[2]->start ? twogr : 2;
		    	return(nextgr);
		    }
		    else {
		    	nextgr = currgroup[0]->start < currgroup[1]->start ? 0 : 1;
		    	return(nextgr);
		    }
		}
		else {
		    if(currgroup[2]!=NULL) {
		    	nextgr = currgroup[0]->start < currgroup[2]->start ? 0 : 2;
		    	return(nextgr);
		    }
		    else {
		    	return(0);
		    }
		}
	} // end if(currgroup[0]!=NULL)
	else {
		if(currgroup[1]!=NULL) {
		    if(currgroup[2]!=NULL) {
		    	nextgr = currgroup[1]->start < currgroup[2]->start ? 1 : 2;
		    	return(nextgr);
		    }
		    else {
		    	return(1);
		    }
		}
		else {
		    return(2);
		}
	}

	return(nextgr);
}

void set_strandcol(CGroup *prevgroup, CGroup *group, int grcol, GVec<int>& eqcol, GVec<int>& equalcolor){

	int zerocol=eqcol[prevgroup->color];
	if(zerocol>-1) {

		while(eqcol[zerocol]!=-1 && eqcol[zerocol]!=zerocol) {
		    zerocol=eqcol[zerocol];
		}
		int tmpcol=zerocol;
		while(equalcolor[zerocol]!=zerocol) {
			zerocol=equalcolor[zerocol];
		}
		eqcol[prevgroup->color]=zerocol;
		eqcol[tmpcol]=zerocol;

		if(zerocol<grcol) {
		    equalcolor[grcol]=zerocol;
		    group->color=zerocol;
		}
		else if(grcol<zerocol) {
		    equalcolor[zerocol]=grcol;
		    eqcol[prevgroup->color]=grcol;
		}
	} // if(zerocol>-1)
	else {
		eqcol[prevgroup->color]=grcol;
	}

}

void add_group_to_bundle(CGroup *group, CBundle *bundle, GPVec<CBundlenode>& bnode){

	CBundlenode *currlastnode=bnode[bundle->lastnodeid];
	int bid=bnode.Count();

	if(group->start > currlastnode->end) { // group after last bnode
		CBundlenode *currbnode=new CBundlenode(group->start,group->end,group->cov_sum,bid);
		currlastnode->nextnode=currbnode;
		bnode.Add(currbnode);
		bundle->lastnodeid=bid;
		bundle->len+=group->end-group->start+1;
		bundle->cov+=group->cov_sum;
		bundle->nread+=group->nread;
		bundle->multi+=group->multi;
	}
	else { // group overlaps bnode
		if(currlastnode->end < group->end) {
		    bundle->len+= group->end - currlastnode->end;
		    currlastnode->end= group->end;
		}
		bundle->cov+=group->cov_sum;
		bundle->nread+=group->nread;
		bundle->multi+=group->multi;
		currlastnode->cov+=group->cov_sum;
	}
}

int create_bundle(GPVec<CBundle>& bundle,CGroup *group,GPVec<CBundlenode>& bnode) {

	int bid=bnode.Count();
	int bno=bundle.Count();
	CBundlenode *startbnode=new CBundlenode(group->start,group->end,group->cov_sum,bid);

	CBundle *newbundle=new CBundle(group->end-group->start+1,group->cov_sum,group->nread,group->multi,bid,bid);
	bundle.Add(newbundle);
	bnode.Add(startbnode);

	return(bno);
}

int setCmp(const pointer p1, const pointer p2) {
	CTrInfo *a=(CTrInfo*)p1;
	CTrInfo *b=(CTrInfo*)p2;
	if(a->trno<b->trno) return -1;
	if(a->trno>b->trno) return 1;
	return 0;
}

int capCmp(const pointer p1, const pointer p2) {
	CTrInfo *a=(CTrInfo*)p1;
	CTrInfo *b=(CTrInfo*)p2;
	if(a->abundance<b->abundance) return -1;
	if(a->abundance>b->abundance) return 1;
	return 0;
}

int trCmp(const pointer p1, const pointer p2) {
	CTransfrag *a=(CTransfrag*)p1;
	CTransfrag *b=(CTransfrag*)p2;
	if(a->pattern.count()<b->pattern.count()) return 1;
	if(a->pattern.count()>b->pattern.count()) return -1;
	if(!a->real && b->real) return -1;
	if(a->real && !b->real) return 1;
	if(a->abundance<b->abundance) return 1;
	if(a->abundance>b->abundance) return -1;
	return 0;
}

int edgeCmp(const pointer p1, const pointer p2) {
	CNetEdge *a=(CNetEdge*)p1;
	CNetEdge *b=(CNetEdge*)p2;
	if(a->rate<b->rate) return 1;
	if(a->rate>b->rate) return -1;
	return 0;
}

int pointCmp(const pointer p1, const pointer p2) {
	CTrimPoint *a=(CTrimPoint*)p1;
	CTrimPoint *b=(CTrimPoint*)p2;
	if(a->start && !b->start) return -1;
	if(!a->start && b->start) return 1;
	if(a->abundance<b->abundance) return 1;
	if(a->abundance>b->abundance) return -1;
	return 0;
}

int edgeCmpEM(const pointer p1, const pointer p2) {
	CNetEdge *a=(CNetEdge*)p1;
	CNetEdge *b=(CNetEdge*)p2;
	if(a->fake && !b->fake) return 1; // check if this is right -> I want the fake one to come last
	if(!a->fake && b->fake) return -1;
	if(a->rate<b->rate) return 1;
	if(a->rate>b->rate) return -1;
	return 0;
}

int guideabundCmp(const pointer p1, const pointer p2) {
	CTransfrag *a=(CTransfrag*)p1;
	CTransfrag *b=(CTransfrag*)p2;
	if(a->abundance<b->abundance) return 1;
	if(a->abundance>b->abundance) return -1;
	if(a->pattern.count()<b->pattern.count()) return 1;
	if(a->pattern.count()>b->pattern.count()) return -1;
	return 0;
}

int guidedabundCmp(const pointer p1, const pointer p2) {
	CGuide *a=(CGuide*)p1;
	CGuide *b=(CGuide*)p2;
	if(a->trf->abundance<b->trf->abundance) return 1;
	if(a->trf->abundance>b->trf->abundance) return -1;
	if(a->trf->pattern.count()<b->trf->pattern.count()) return 1;
	if(a->trf->pattern.count()>b->trf->pattern.count()) return -1;
	return 0;
}

int juncCmpEnd(const pointer p1, const pointer p2) {
	CJunction* a=(CJunction*)p1;
	CJunction* b=(CJunction*)p2;
	if (a->end<b->end) return -1;
	if (a->end>b->end) return 1;
	if (a->start<b->start) return -1;
	if (a->start>b->start) return 1;
	return 0;
}

void merge_fwd_groups(GPVec<CGroup>& group, CGroup *group1, CGroup *group2, GVec<int>& merge, GVec<int>& eqcol) {

	// get end of group (group1 is assumed to come before group2)
	group1->end=group2->end;

	// get smallest color of group
	while(eqcol[group1->color]!=group1->color) {
		group1->color=eqcol[group1->color];
	}
	while(eqcol[group2->color]!=group2->color) {
		group2->color=eqcol[group2->color];
	}

	if(group1->color<group2->color) {
	   eqcol[group2->color]=group1->color;
	   // print STDERR "in merge set eqcol[",$$group2[2],"]=",$$eqcol{$$group1[2]},"\n";
	}
	else if(group1->color>group2->color) {
		eqcol[group1->color]=group2->color;
		// print STDERR "in merge set eqcol[",$$group1[2],"]=",$$group2[2],"\n";
		group1->color=group2->color;
	}

	group1->cov_sum+=group2->cov_sum;
	group1->next_gr=group2->next_gr; // this is possible because group1->next_gr=group2

	merge[group2->grid]=group1->grid;

	group1->nread+=group2->nread;
	group1->multi+=group2->multi;

	// delete group2
	group.freeItem(group2->grid);
}


int merge_read_to_group(int n,int np, int p, float readcov, int sno,int readcol,GList<CReadAln>& readlist,int color,GPVec<CGroup>& group,CGroup **allcurrgroup,
		//CGroup **startgroup,GVec<int> *readgroup,GVec<int>& eqcol,GVec<int>& merge,float& fraglen,int *usedcol) {
		CGroup **startgroup,GVec<int> *readgroup,GVec<int>& eqcol,GVec<int>& merge,int *usedcol) {

	//fprintf(stderr,"merge readcol=%d\n",readcol);

	CGroup *currgroup=allcurrgroup[sno];

	if(currgroup != NULL) { // this type of group - negative, unknown, or positive - was created before

		//set currgroup first
		CGroup *lastgroup=NULL;
		while(currgroup!=NULL && readlist[n]->start > currgroup->end) { // while read start after the current group's end advance group -> I might have more groups leaving from current group due to splicing
		    lastgroup=currgroup;
		    currgroup=currgroup->next_gr;
		}

		if(currgroup==NULL || readlist[n]->segs[0].end < currgroup->start) // currgroup is null only if we reached end of currgroup list because currgroup is not NULL initially
			currgroup=lastgroup;

		// now process each group of coordinates individually
		CGroup *thisgroup=currgroup;
		int ncoord=readlist[n]->segs.Count(); // number of "exons" in read
		int lastpushedgroup=-1;
		bool added=false;

		for(int i=0;i<ncoord;i++) {

			//fraglen+=readcov*readlist[n]->segs[i].len();  // this might be useful to have if I decide not to use the HI tag anymore

		    // skip groups that are left behind
		    while(thisgroup!=NULL && readlist[n]->segs[i].start > thisgroup->end) { // find the group where "exon" fits
		    	lastgroup=thisgroup;
		    	thisgroup=thisgroup->next_gr;
		    }

		    if(thisgroup && readlist[n]->segs[i].end >= thisgroup->start) { // read overlaps group

		    	// I need to split pairs here if color didn't reach this group: it means there is a gap between these groups and no reads joining them
		    	if(!i && np>-1 && readlist[np]->nh && np<n) { // I only consider first exon here because the rest of the groups need to get the same color

		    		int grouppair=thisgroup->grid;
		    		while( merge[grouppair]!=grouppair) {
		    			grouppair=merge[grouppair];
		    		}
		    		thisgroup->grid=grouppair;

		    		int thiscol=thisgroup->color;
		    		while(eqcol[thiscol]!=thiscol) { // get smallest color
		    			thiscol=eqcol[thiscol];
		    		}
		    		thisgroup->color=thiscol;

		    		if(thiscol!=readcol) { // pair color didn't reach this group

		    			//fprintf(stderr,"Split pairs: %d-%d and %d-%d on strand %d\n",readlist[np]->start,readlist[np]->end,readlist[n]->start,readlist[n]->end,readlist[n]->strand);
		    			readlist[n]->pair_idx[p]=-1;
		    			for(int j=0;j<readlist[np]->pair_idx.Count();j++)
		    				if(readlist[np]->pair_idx[j]==n) {
		    					readlist[np]->pair_idx[j]=-1;
		    					break;
		    				}

		    			readcol=usedcol[sno]; // readcol gets back the new color
		    			if(readcol<0) {
		    				usedcol[sno]=color;
		    				readcol=color;
		    				eqcol.Add(color);
		    				color++;
		    			}
		    		}
		    	}


		    	if(!added) { // read gets added only to first group - why ??
		    		thisgroup->nread+=readcov;
		    		if(readlist[n]->nh>1) thisgroup->multi+=readcov;  // this will probably need to be coded differently if I do super-reads, or collapse more than one read into the same nh
		    		added=true;
		    	}

		    	if(readlist[n]->segs[i].start<thisgroup->start) {
		    		thisgroup->start=readlist[n]->segs[i].start;
		    	}

		    	// find end of new group
		    	CGroup *nextgroup=thisgroup->next_gr;
		    	while(nextgroup!=NULL && readlist[n]->segs[i].end >= nextgroup->start) {
		    		merge_fwd_groups(group,thisgroup,nextgroup,merge,eqcol);
		    		nextgroup=thisgroup->next_gr;
		    	}
		    	if(readlist[n]->segs[i].end > thisgroup->end) {
		    		thisgroup->end=readlist[n]->segs[i].end;
		    	}

				// get smallest color of group
		    	while(eqcol[thisgroup->color]!=thisgroup->color) {
		    		thisgroup->color=eqcol[thisgroup->color];
		    	}

		    	if(readcol!=thisgroup->color) { // read color is different from group color
		    		if(readcol<thisgroup->color) { // set group color to current read color
		    			eqcol[thisgroup->color]=readcol;
		    			thisgroup->color=readcol;
		    		}
		    		else { // read color is bigger than group

		    			eqcol[readcol]=thisgroup->color;
		    			readcol=thisgroup->color;
		    		}
		    	}

		    	if(thisgroup->grid != lastpushedgroup) {
		    		readgroup[n].Add(thisgroup->grid);   // readgroup for read n gets the id of group
		    		lastpushedgroup=thisgroup->grid;
		    	}

		    	thisgroup->cov_sum+=(readlist[n]->segs[i].end-readlist[n]->segs[i].start+1)*readcov; // coverage is different than number of reads
		    } // end if(thisgroup && readlist[n]->segs[i].end >= thisgroup->start)
		    else { // read is at the end of groups, or read is not overlapping other groups -> lastgroup is not null here because currgroup was not null

		    	// I need to split pairs here because I have no overlap => color didn't reach here for sure if I am at the first exon
		    	if(!i && np>-1 && readlist[np]->nh && np<n) {

		    		//fprintf(stderr,"Split pairs: %d-%d and %d-%d on strand %d\n",readlist[np]->start,readlist[np]->end,readlist[n]->start,readlist[n]->end,readlist[n]->strand);
	    			readlist[n]->pair_idx[p]=-1;
	    			for(int j=0;j<readlist[np]->pair_idx.Count();j++)
	    				if(readlist[np]->pair_idx[j]==n) {
	    					readlist[np]->pair_idx[j]=-1;
	    					break;
	    				}

	    			readcol=usedcol[sno]; // readcol gets back the new color
	    			if(readcol<0) {
	    				usedcol[sno]=color;
	    				readcol=color;
	    				eqcol.Add(color);
	    				color++;
	    			}

		    	}

		    	int ngroup=group.Count();
		    	float nread=0;
		    	float multi=0;
		    	if(!added) {
		    		nread=readcov;
		    		if(readlist[n]->nh>1) multi=readcov;
		    		added=true;
		    	}
		    	CGroup *newgroup=new CGroup(readlist[n]->segs[i].start,readlist[n]->segs[i].end,readcol,ngroup,(readlist[n]->segs[i].end-readlist[n]->segs[i].start+1)*readcov,nread,multi);
		    	group.Add(newgroup);
		    	merge.Add(ngroup);
		    	lastgroup->next_gr=newgroup; // can lastgroup be null here -> no from the way I got here
		    	newgroup->next_gr=thisgroup;
		    	readgroup[n].Add(ngroup);
		    	lastpushedgroup=ngroup;
				thisgroup=lastgroup;
		    }
		} // for(int i=0;i<ncoord;i++)
	} // if(currgroup != NULL)
	else { // create new group of this type

		// I need to split pairs here because I have no overlap => color didn't reach here for sure
		if(np>-1 && readlist[np]->nh && np<n) {

			//fprintf(stderr,"Split pairs: %d-%d and %d-%d on strand %d\n",readlist[np]->start,readlist[np]->end,readlist[n]->start,readlist[n]->end,readlist[n]->strand);
			readlist[n]->pair_idx[p]=-1;
			for(int j=0;j<readlist[np]->pair_idx.Count();j++)
				if(readlist[np]->pair_idx[j]==n) {
					readlist[np]->pair_idx[j]=-1;
					break;
				}

			readcol=usedcol[sno]; // readcol gets back the new color
			if(readcol<0) {
				usedcol[sno]=color;
				readcol=color;
				eqcol.Add(color);
				color++;
			}

		}

		int ncoord=readlist[n]->segs.Count();
		CGroup *lastgroup=NULL;
		float nread=readcov;
		float multi=0;
		if(readlist[n]->nh>1) multi=readcov;
		for(int i=0;i<ncoord;i++) {

			//fraglen+=readcov*readlist[n]->segs[i].len();

			int ngroup=group.Count();
			CGroup *newgroup=new CGroup(readlist[n]->segs[i].start,readlist[n]->segs[i].end,readcol,ngroup,(readlist[n]->segs[i].end-readlist[n]->segs[i].start+1)*readcov,nread,multi);
			nread=0;
			multi=0;
			group.Add(newgroup);
			merge.Add(ngroup);
			if(lastgroup!=NULL) {
				lastgroup->next_gr=newgroup;
			}
			else {
				currgroup=newgroup;
			}
			lastgroup=newgroup;
			readgroup[n].Add(ngroup);
		}

	}

	allcurrgroup[sno]=currgroup;

	if(startgroup[sno]==NULL) startgroup[sno]=currgroup;

	return color;
}

int add_read_to_group(int n,GList<CReadAln>& readlist,int color,GPVec<CGroup>& group,CGroup **allcurrgroup,
		//CGroup **startgroup,GVec<int> *readgroup,GVec<int>& eqcol,GVec<int>& merge,float& fraglen,uint& fragno) {
		CGroup **startgroup,GVec<int> *readgroup,GVec<int>& eqcol,GVec<int>& merge) {

	int usedcol[3]={-1,-1,-1};

	float single_count=readlist[n]->read_count; // need to compute read single count

	for(int p=0;p<readlist[n]->pair_idx.Count();p++) {

		// for each pair I have to pretend is independent of the previous one since I grouped together several reads which might be independent
		int sno=readlist[n]->strand+1; // 0: negative strand; 1: zero strand; 2: positive strand (starting from -1,0,1)

		// check if I've seen read's pair and if yes get its readcol; at the least get read's pair strand if available
		int np=readlist[n]->pair_idx[p]; // pair read number

		if(np>-1 && readlist[np]->nh) { // read pair exists and it wasn't deleted

			// see if I have the correct read strand
			char snop=readlist[np]->strand+1;  // snop is the strand of pair read
			if(sno!=snop) { // different strands for read and pair
				if(sno==1) { // read n is on zero (neutral) strand, but pair has strand
					// readlist[n]->strand=readlist[np]->strand; I can not update the read strand anymore -> REMEMBER this later
					sno=snop;  // assign strand of pair to read
				}
				else if(snop!=1) { // conflicting strands -> un-pair reads in the hope that one is right
					readlist[n]->pair_idx[p]=-1;
					for(int j=0;j<readlist[np]->pair_idx.Count();j++)
						if(readlist[np]->pair_idx[j]==n) {
							readlist[np]->pair_idx[j]=-1;
							break;
						}
					np=-1;
				}
			}

			if(np>-1) { // reads didn't get split

				single_count-=readlist[n]->pair_count[p]; // update single count of read here

				int readcol=usedcol[sno];

				if(np<n) { // there is a pair and it came before the current read in sorted order of position
					// first group of pair read is: $$readgroup[$np][0]

					int grouppair=readgroup[np][0];
					while( merge[grouppair]!=grouppair) {
						grouppair=merge[grouppair];
					}
					readgroup[np][0]=grouppair;

					readcol=group[readgroup[np][0]]->color;    // readcol gets assigned the color of the pair's group
					while(eqcol[readcol]!=readcol) { // get smallest color
						readcol=eqcol[readcol];
					}
					//print STDERR "Adjust color of group ",$$readgroup[$np][0]," to $readcol\n";
					group[readgroup[np][0]]->color=readcol;
				}
				else { // it's the first time I see the read in the fragment
					//fragno+=readlist[n]->pair_count[p];
					if(usedcol[sno]<0) { // I didn't use the color yet
						usedcol[sno]=color;
						readcol=color;
						eqcol.Add(color);  // read might keep the new color so we need to add it to the eqcol
						color++;
					}
				}

				color=merge_read_to_group(n,np,p,readlist[n]->pair_count[p],sno,readcol,readlist,color,group,
										allcurrgroup,startgroup,readgroup,eqcol,merge,usedcol);


			} // this ends if(np>-1)
		} // if(np>-1 && readlist[np]->nh) : read pair exists and it wasn't deleted
	} // for(int p=0,...)

	// now I need to deal with the single count
	if(single_count>epsilon) { // my way of controlling for rounding errors
		//fragno+=single_count;
		int readcol=usedcol[readlist[n]->strand+1];
		if(readcol<0) {
			readcol=color;
			usedcol[readlist[n]->strand+1]=color;
			eqcol.Add(color); // read might keep the new color so we need to add it to the eqcol
			color++;
		}
		color=merge_read_to_group(n,-1,-1,single_count,readlist[n]->strand+1,readcol,readlist,color,group,
								allcurrgroup,startgroup,readgroup,eqcol,merge,usedcol);

	}

	return color;
}

CGraphnode *create_graphnode(int s, int g, uint start,uint end,int nodeno,CBundlenode *bundlenode,
		GVec<CGraphinfo> **bundle2graph,GPVec<CGraphnode> **no2gnode) {

	/*
	{ // DEBUG ONLY
		fprintf(stderr,"create_graphnode[%d][%d]:%d-%d\n",s,g,start,end);
	}
	*/

	CGraphnode* gnode=new CGraphnode(start,end,nodeno);
	CGraphinfo ginfo(g,nodeno);
	bundle2graph[s][bundlenode->bid].Add(ginfo);
	no2gnode[s][g].Add(gnode);

	return(gnode);
}

float compute_chi(GArray<float>& winleft, GArray<float>& winright, float sumleft, float sumright) {

	float chi=0;
	for(int j=0;j<CHI_WIN;j++) {
		float mul=(winleft[j]+ winright[j])/(sumleft+sumright);
		float mur=mul;
		mul*=sumleft;
		mur*=sumright;
		chi+= (winleft[j]-mul)/mul+(winright[j]-mur)/mur;
	}
	return(chi);
}

void find_trims(int refstart,int sno,uint start,uint end,GVec<float>* bpcov,uint& sourcestart,float& maxsourceabundance,uint& sinkend,
		float& maxsinkabundance){

	if(end-start<2*CHI_WIN-1) return;

	float sumleft=0;
	float sumright=0;

	//float maxsinkchi=0;
	//float maxsourcechi=0;

	float sinkabundance=0;
	float sourceabundance=0;

	GArray<float> winleft(CHI_WIN,false); // not auto-sort
	GArray<float> winright(CHI_WIN,false); // not auto-sort

	float cov;

	for(uint i=start;i<=end;i++) {

		if(i-start<2*CHI_WIN-1)  { // I have to compute the sumleft and sumright first
			if(i-start<CHI_WIN) {
				cov=bpcov[sno][i-refstart];
				if(bpcov[1][i-refstart]>cov) cov+=(bpcov[1][i-refstart]-bpcov[0][i-refstart]-bpcov[2][i-refstart])/bpcov[1][i-refstart];
				sumleft+=cov;
				winleft.Add(cov);
			}
			else {
				cov=bpcov[sno][i-refstart];
				if(bpcov[1][i-refstart]>cov) cov+=(bpcov[1][i-refstart]-bpcov[0][i-refstart]-bpcov[2][i-refstart])/bpcov[1][i-refstart];
				sumright+=cov;
				winright.Add(cov);
				if(i-start==2*CHI_WIN-2) {
					winleft.setSorted(true);
					winright.setSorted(true);
				}
			}
	    }
	    else { // I can do the actual sumleft, sumright comparision
			cov=bpcov[sno][i-refstart];
			if(bpcov[1][i-refstart]>cov) cov+=(bpcov[1][i-refstart]-bpcov[0][i-refstart]-bpcov[2][i-refstart])/bpcov[1][i-refstart];
	    	sumright+=cov;
			winright.Add(cov);

			float chi=0;
			if(sumleft!=sumright) chi=compute_chi(winleft,winright,sumleft,sumright);

			if(chi>CHI_THR) { // there is a significant difference
				if(sumleft>sumright) { // possible drop (sink cut)
					//if(chi>maxsinkchi) {
					sinkabundance=(sumleft-sumright)/CHI_WIN;
					if(maxsinkabundance<sinkabundance) {
						sinkend=i-CHI_WIN;
						//maxsinkchi=chi;
						maxsinkabundance=sinkabundance;
					}
				}
				else if(sumright>sumleft) { // increase (source cut)
					//if(chi>maxsourcechi) {
					sourceabundance=(sumright-sumleft)/CHI_WIN;
					if(maxsourceabundance<sourceabundance) {
						sourcestart=i-CHI_WIN+1;
						maxsourceabundance=sourceabundance;
						//maxsourcechi=chi;
					}
				}
	    	}

			cov=bpcov[sno][i-refstart-2*CHI_WIN+1];
			if(bpcov[1][i-refstart-2*CHI_WIN+1]>cov)
				cov+=(bpcov[1][i-refstart-2*CHI_WIN+1]-bpcov[0][i-refstart-2*CHI_WIN+1]-bpcov[2][i-refstart-2*CHI_WIN+1])/bpcov[1][i-refstart-2*CHI_WIN+1];
	    	sumleft-=cov;
	    	int idx=winleft.IndexOf(cov);
	    	if(idx>=0) winleft.Delete(idx);

			cov=bpcov[sno][i-refstart-CHI_WIN+1];
			if(bpcov[1][i-refstart-CHI_WIN+1]>cov)
				cov+=(bpcov[1][i-refstart-CHI_WIN+1]-bpcov[0][i-refstart-CHI_WIN+1]-bpcov[2][i-refstart-CHI_WIN+1])/bpcov[1][i-refstart-CHI_WIN+1];
	    	sumleft+=cov;
	    	winleft.Add(cov);
	    	sumright-=cov;
	    	idx=winright.IndexOf(cov);
	    	if(idx>=0) winright.Delete(idx);
	    }
	}

}

CGraphnode *add_trim_to_graph(int s, int g,uint lastpos,CTrimPoint& mytrim,CGraphnode *graphnode,CGraphnode *source,CGraphnode *sink,GVec<float>& futuretr,
		int& graphno,CBundlenode *bundlenode,GVec<CGraphinfo> **bundle2graph,GPVec<CGraphnode> **no2gnode) {

	if(mytrim.start) { // this is a source link
		float tmp=graphno-1;
		CGraphnode *prevnode=NULL;
		if(mytrim.pos>graphnode->start) { // there is place for another node
			uint prevend=graphnode->end;
			graphnode->end=mytrim.pos-1;
			prevnode=graphnode;
			graphnode=create_graphnode(s,g,mytrim.pos,prevend,graphno,bundlenode,bundle2graph,no2gnode);
			graphno++;
		}
		source->child.Add(graphnode->nodeid);  // this node is the child of source
		graphnode->parent.Add(source->nodeid); // this node has source as parent
		if(prevnode) {
			prevnode->child.Add(graphnode->nodeid); // this node is the child of previous node
			graphnode->parent.Add(prevnode->nodeid); // this node has as parent the previous node
		}
		futuretr.Add(tmp);
		futuretr.cAdd(0.0);
		futuretr.Add(mytrim.abundance);
	}
	else { // this is a link to sink
		float tmp=graphno-1;
		if(mytrim.pos<lastpos) { // there is still place for a new node
			uint prevend=graphnode->end;
			graphnode->end=mytrim.pos;
			CGraphnode *prevnode=graphnode;
			graphnode=create_graphnode(s,g,mytrim.pos+1,prevend,graphno,bundlenode,bundle2graph,no2gnode);
			graphno++;
			prevnode->child.Add(graphnode->nodeid); // this node is the child of previous node
			graphnode->parent.Add(prevnode->nodeid); // this node has as parent the previous node
			sink->parent.Add(prevnode->nodeid); // prevnode is the parent of sink
		}
		else {
			sink->parent.Add(graphnode->nodeid); // prevnode is the parent of sink
		}
		futuretr.Add(tmp);
		futuretr.cAdd(1.0);
		futuretr.Add(mytrim.abundance);
	}

	return(graphnode);
}

CGraphnode *trimnode(int s, int g, int refstart,uint newend, CGraphnode *graphnode,CGraphnode *source, CGraphnode *sink, GVec<float>* bpcov,
		GVec<float>& futuretr, int& graphno,CBundlenode *bundlenode,GVec<CGraphinfo> **bundle2graph,GPVec<CGraphnode> **no2gnode, int &edgeno) {

	uint sourcestart=0;
	uint sinkend=0;
	float sinkabundance=0;
	float sourceabundance=0;
	find_trims(refstart,2*s,graphnode->start,newend,bpcov,sourcestart,sourceabundance,sinkend,sinkabundance);

	if(sourcestart < sinkend) { // source trimming comes first

		if(sourcestart) { // there is evidence of graphnode trimming from source
			graphnode->end=sourcestart-1;
			CGraphnode *prevnode=graphnode;
			graphnode=create_graphnode(s,g,sourcestart,newend,graphno,bundlenode,bundle2graph,no2gnode);
			graphno++;
			source->child.Add(graphnode->nodeid);  // this node is the child of source
			graphnode->parent.Add(source->nodeid); // this node has source as parent
			//fprintf(stderr,"edge 0-%d",graphnode->nodeid);
			prevnode->child.Add(graphnode->nodeid); // this node is the child of previous node
			graphnode->parent.Add(prevnode->nodeid); // this node has as parent the previous node
			//fprintf(stderr,"edge %d-%d",prevnode->nodeid,graphnode->nodeid);
			float tmp=graphno-1;
			futuretr.cAdd(0.0);
			futuretr.Add(tmp);
			sourceabundance+=trthr;futuretr.Add(sourceabundance);
			tmp=prevnode->nodeid;futuretr.Add(tmp);
			tmp=graphnode->nodeid;futuretr.Add(tmp);
			tmp=trthr;futuretr.Add(tmp);
			// COUNT 2 EDGES HERE
			edgeno+=2;
		}

		// sinkend is always positive since it's bigger than sourcestart
		float tmp=graphno-1;
		graphnode->end=sinkend;
		CGraphnode *prevnode=graphnode;
		graphnode=create_graphnode(s,g,sinkend+1,newend,graphno,bundlenode,bundle2graph,no2gnode);
		graphno++;
		prevnode->child.Add(graphnode->nodeid); // this node is the child of previous node
		graphnode->parent.Add(prevnode->nodeid); // this node has as parent the previous node
		//fprintf(stderr,"edge %d-%d",prevnode->nodeid,graphnode->nodeid);
		sink->parent.Add(prevnode->nodeid); // prevnode is the parent of sink
		// remember to create transfrag as well -> I don't know the gno yet, so I can not create it here
		//fprintf(stderr,"edge %d-sink",graphnode->nodeid);
		futuretr.Add(tmp);
		futuretr.cAdd(-1.0);
		sinkabundance+=trthr;futuretr.Add(sinkabundance);
		tmp=prevnode->nodeid;futuretr.Add(tmp);
		tmp=graphnode->nodeid;futuretr.Add(tmp);
		tmp=trthr;futuretr.Add(tmp);
		// COUNT 2 EDGES HERE
		edgeno+=2;
	}
	else if(sourcestart > sinkend) { // sink trimming comes first
		if(sinkend) {
			graphnode->end=sinkend;
			CGraphnode *prevnode=graphnode;
			graphnode=create_graphnode(s,g,sinkend+1,newend,graphno,bundlenode,bundle2graph,no2gnode);
			graphno++;
			prevnode->child.Add(graphnode->nodeid); // this node is the child of previous node
			graphnode->parent.Add(prevnode->nodeid); // this node has as parent the previous node
			//fprintf(stderr,"edge %d-%d",prevnode->nodeid,graphnode->nodeid);
			sink->parent.Add(prevnode->nodeid); // prevnode is the parent of sink
			//fprintf(stderr,"edge %d-sink",graphnode->nodeid);
			// remember to create transfrag as well -> I don't know the gno yet, so I can not create it here
			float tmp=graphno-2;
			futuretr.Add(tmp);
			futuretr.cAdd(-1.0);
			sinkabundance+=trthr;futuretr.Add(sinkabundance);
			tmp=prevnode->nodeid;futuretr.Add(tmp);
			tmp=graphnode->nodeid;futuretr.Add(tmp);
			tmp=trthr;futuretr.Add(tmp);
			// COUNT 2 EDGES HERE
			edgeno+=2;
		}

		// sourcestart is positive since it's bigger than sinkend
		graphnode->end=sourcestart-1;
		CGraphnode *prevnode=graphnode;
		graphnode=create_graphnode(s,g,sourcestart,newend,graphno,bundlenode,bundle2graph,no2gnode);
		graphno++;
		source->child.Add(graphnode->nodeid);  // this node is the child of source
		graphnode->parent.Add(source->nodeid); // this node has source as parent
		//fprintf(stderr,"edge 0-%d",graphnode->nodeid);
		prevnode->child.Add(graphnode->nodeid); // this node is the child of previous node
		graphnode->parent.Add(prevnode->nodeid); // this node has as parent the previous node
		//fprintf(stderr,"edge %d-%d",prevnode->nodeid,graphnode->nodeid);
		float tmp=graphno-1;
		futuretr.cAdd(0.0);
		futuretr.Add(tmp);
		sourceabundance+=trthr;futuretr.Add(sourceabundance);
		tmp=prevnode->nodeid;futuretr.Add(tmp);
		tmp=graphnode->nodeid;futuretr.Add(tmp);
		tmp=trthr;futuretr.Add(tmp);
		// COUNT 2 EDGES HERE
		edgeno+=2;
	}
	// else both source and sink trimming are not present

	return(graphnode);
}

inline int edge(int min, int max, int gno) {
	//return((gno-1)*min-min*(min-1)/2+max-min); // this should be changed if source to node edges are also stored
	return((gno-1)*(min+1)-min*(min-1)/2+max-min); // this includes source to node edges
}

GBitVec traverse_dfs(int s,int g,CGraphnode *node,CGraphnode *sink,GBitVec parents,int gno, GVec<bool>& visit,
		GPVec<CGraphnode> **no2gnode,GPVec<CTransfrag> **transfrag, int &edgeno,GIntHash<int> **gpos,int &lastgpos){

	if(visit[node->nodeid]) {
		node->parentpat = node->parentpat | parents;
		for(int n=0;n<gno;n++) {
			if(parents[n]) // add node's children to all parents of node
				no2gnode[s][g][n]->childpat = no2gnode[s][g][n]->childpat | node->childpat;
			else if(node->childpat[n])
				no2gnode[s][g][n]->parentpat = no2gnode[s][g][n]->parentpat | node->parentpat;
		}
	}
	else {
		node->childpat.resize(gno+edgeno);
		node->parentpat.resize(gno+edgeno);
		node->parentpat = node->parentpat | parents;
		visit[node->nodeid]=true;
		parents[node->nodeid]=1; // add the node to the parents

		if(node->parent.Count()==1 && !node->parent[0]) { // node has source only as parent -> add transfrag from source to node
			GBitVec trpat(gno+edgeno);
			trpat[0]=1;
			trpat[node->nodeid]=1;

			int key=edge(0,node->nodeid,gno);
			int *pos=gpos[s][g][key];
			if(pos!=NULL) trpat[*pos]=1;
			else {
				//fprintf(stderr,"s=%d g=%d key=%d lastgpos=%d add edge between %d and %d\n",s,g,key,lastgpos,0,node->nodeid);
				gpos[s][g].Add(key,lastgpos);
				trpat[lastgpos]=1;
				lastgpos++;
			}

			GVec<int> nodes;
			nodes.cAdd(0);
			nodes.Add(node->nodeid);
			CTransfrag *tr=new CTransfrag(nodes,trpat,trthr);
			/*
			{ // DEBUG ONLY
				fprintf(stderr,"Add source transfrag[%d][%d]= %d and pattern",s,g,transfrag[s][g].Count());
				//printBitVec(trpat);
				fprintf(stderr,"\n");
			}
			*/
			transfrag[s][g].Add(tr);
		}

		int n=node->child.Count();
		if(node != sink && !n) {
			node->child.Add(sink->nodeid);  // add sink to the node's children
			sink->parent.Add(node->nodeid); // add node to sink's parents
			// create the transfrag that ends the node
			GBitVec trpat(gno+edgeno);
			trpat[node->nodeid]=1;
			trpat[gno-1]=1;

			int key=edge(node->nodeid,gno-1,gno);
			int *pos=gpos[s][g][key];
			if(pos!=NULL) trpat[*pos]=1;
			else {
				gpos[s][g].Add(key,lastgpos);
				//fprintf(stderr,"s=%d g=%d key=%d lastgpos=%d add edge between %d and %d\n",s,g,key,lastgpos, node->nodeid,gno-1);
				trpat[lastgpos]=1;
				lastgpos++;
			}

			GVec<int> nodes;
			nodes.Add(node->nodeid);
			nodes.Add(sink->nodeid);
			CTransfrag *tr=new CTransfrag(nodes,trpat,trthr);
			/*
			{ // DEBUG ONLY
				fprintf(stderr,"Add sink transfrag[%d][%d]= %d for nodeid=%d and pattern:",s,g,transfrag[s][g].Count(),node->nodeid);
				//printBitVec(trpat);
				fprintf(stderr,"\n");
			}
			*/
			transfrag[s][g].Add(tr);
			n++;
	    }
		//fprintf(stderr,"Add %d children of node %d: ",n,node->nodeid);
		//for(int i=0;i<n;i++) fprintf(stderr," %d",node->child[i]);
		//fprintf(stderr,"\n");
		//edgeno+=n; // this will have to be deleted in the end; now I put it so that I can check equivalence with the one computed when creating the graph

	    for(int i=0; i< n; i++) { // for all children
	    	GBitVec childparents=parents;
	    	int min=node->nodeid; // nodeid is always smaller than child node ?
	    	int max=node->child[i];
	    	if(min>max) {
	    		max=node->nodeid; // nodeid is always smaller than child node ?
	    		min=node->child[i];
	    	}

			int key=edge(min,max,gno);
			int *pos=gpos[s][g][key];
			if(pos!=NULL) {
				childparents[*pos]=1; // add edge from node to child to the set of parents from child
				node->childpat[*pos]=1; // add edge from node to child to the set of node children
			}
			else {
				gpos[s][g].Add(key,lastgpos);
				//fprintf(stderr,"s=%d g=%d key=%d lastgpos=%d add edge between %d and %d child\n",s,g,key,lastgpos,min, max);
				childparents[lastgpos]=1;
				node->childpat[lastgpos]=1;
				lastgpos++;
			}

	    	node->childpat = node->childpat | traverse_dfs(s,g,no2gnode[s][g][node->child[i]],sink,childparents,gno,visit,no2gnode,transfrag,edgeno,gpos,lastgpos);
	    }
	} // end else from if(visit[node->nodeid])

	GBitVec children = node->childpat;
	children[node->nodeid]=1;

	return(children);
}

int create_graph(int refstart,int s,int g,CBundle *bundle,GPVec<CBundlenode>& bnode, GList<CJunction>& junction,GList<CJunction>& ejunction,GVec<CGraphinfo> **bundle2graph,
		GPVec<CGraphnode> **no2gnode,GPVec<CTransfrag> **transfrag,GIntHash<int> **gpos,GVec<float>* bpcov,int &edgeno, int &lastgpos){

	CGraphnode* source=new CGraphnode(0,0,0);
	no2gnode[s][g].Add(source);
	CGraphnode* sink=new CGraphnode();

	int njunctions=junction.Count();

	/*
	{ // DEBUG ONLY
		fprintf(stderr,"Junctions[%d][%d]: ",s,g);
		for(int i=0;i<njunctions;i++) fprintf(stderr," %d-%d:%d",junction[i]->start,junction[i]->end,junction[i]->strand);
		fprintf(stderr,"\n");
	}
	*/

	int njs=0; // index of sorted junction starts
	int nje=0; // index of sorted junction ends

	int graphno=1; // number of nodes in graph
	GHash<GVec<int> > ends; // keeps ids of all nodes ending at a certain position; OR ALL NODES THAT ARE LINKED BY JUNCTIONS TO A CERTAIN POSITION

	CBundlenode *bundlenode=bnode[bundle->startnode];

	GVec<float> futuretr;

	while(bundlenode!=NULL) {

	    uint currentstart=bundlenode->start; // current start is bundlenode's start
	    uint endbundle=bundlenode->end; // initialize end with bundlenode's end for now

	    CGraphnode *graphnode=create_graphnode(s,g,currentstart,endbundle,graphno,bundlenode,bundle2graph,no2gnode); // creates a $graphno graphnode  with start at bundle start, and end at bundle end
	    graphno++;

	    int end=0;
	    while(nje<njunctions && ejunction[nje]->end<=currentstart) { // read all junction ends at or before the current start -> assuming there are any (at this point, smaller junction ends should not be relevant to this bundle/currentstart
	      if(ejunction[nje]->end==currentstart && (ejunction[nje]->strand+1) == 2*s) { // junction ends at current start and is on the same strand and not deleted
	    	  end=1;
	      }
	      nje++;
	    }

	    if(end) { // I might have nodes finishing here; but I have a junction finishing here for sure
	    	GStr cs((int)currentstart);
	    	GVec<int> *e=ends[cs.chars()]; // HOW CAN I HAVE MORE THAN ONE NODE FINISHING HERE???; because this keeps all nodes that are linked by junctions here
	    	if(e) {
	    		for(int i=0;i<e->Count();i++) {
	    			CGraphnode *node=no2gnode[s][g][e->Get(i)];
	    			node->child.Add(graphnode->nodeid);  // this node is the child of previous node
	    			graphnode->parent.Add(node->nodeid); // this node has as parent the previous node
	    			// COUNT EDGE HERE
	    			edgeno++;
	    			//fprintf(stderr,"1 edge %d(%d-%d)-%d(%d-%d)\n",node->nodeid,node->start,node->end,graphnode->nodeid,graphnode->start,graphnode->end);
	    		}
	    	}
	    	else { // I haven't seen nodes before that finish here => link to source
		    	source->child.Add(graphnode->nodeid);  // this node is the child of source
		    	graphnode->parent.Add(source->nodeid); // this node has source as parent
		    	// COUNT EDGE HERE
    			edgeno++;
    			//fprintf(stderr,"2 edge 0-%d(%d-%d)\n",graphnode->nodeid,graphnode->start,graphnode->end);
	    	}
	    }
	    else { // this node comes from source directly
	    	source->child.Add(graphnode->nodeid);  // this node is the child of source
	    	graphnode->parent.Add(source->nodeid); // this node has source as parent
	    	// COUNT EDGE HERE
			edgeno++;
			//fprintf(stderr,"3 edge 0-%d(%d-%d)\n",graphnode->nodeid,graphnode->start,graphnode->end);
	    }


	    bool completed=false;

	    do {

	    	//fprintf(stderr,"\ncompleted=%d njunctions=%d nje=%d njs=%d junc[njs]->start=%d currentsart=%d\n",completed,njunctions,nje,njs,junction[njs]->start,currentstart);

	    	while(nje<njunctions && (((int)ejunction[nje]->strand+1) != 2*s)) nje++; // skip junctions that don't have the same strand
	    	while(njs<njunctions && ((((int)junction[njs]->strand+1)!= 2*s) || (junction[njs]->start<currentstart))) njs++; // junctions that start before the current graphnode and I haven't seen them before are part of a different bundle

	    	//fprintf(stderr,"\ns=%d nje=%d njs=%d junc[njs]->start=%d junc[njs]->strand=%d unc[nje]->start=%d junc[nje]->strand=%d currentsart=%d\n",s,nje,njs,junction[njs]->start,junction[njs]->strand,ejunction[nje]->start,ejunction[nje]->strand,currentstart);

	    	int minjunction = -1; // process next junction -> either a start or an ending whichever has the first position on the genome; if they have same position then process ending first
	    	if((nje<njunctions && (ejunction[nje]->end<endbundle)) || (njs<njunctions && (junction[njs]->start<=endbundle))) {
	    		if(nje<njunctions) { // there are still junctions endings
	    			if(njs<njunctions) { // there are still junctions starting
	    				minjunction = junction[njs]->start >= ejunction[nje]->end ? 1 : 0; // one of them is clearly before the endbundle from the initial if
	    			}
	    			else minjunction = 1;
	    		}
	    		else minjunction = 0;
	    	}

	    	//fprintf(stderr,"minjunction=%d\n",minjunction);

	    	if(minjunction == 0 ) { // found a start junction here

	    		if(trim) graphnode=trimnode(s,g,refstart,junction[njs]->start,graphnode,source,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);// do something to find intermediate nodes; alternatively, I could only do this for end nodes

	    		// if no trimming required just set the end of the node
	    		graphnode->end=junction[njs]->start; // set the end of current graphnode to here; introduce smaller nodes if trimming is activated
	    		uint pos=junction[njs]->start;
	    		while(njs<njunctions && junction[njs]->start==pos ) { // remember ends here
	    			if((junction[njs]->strand+1) == 2*s) {
	    				GStr je((int)junction[njs]->end);
	    				GVec<int> *e=ends[je.chars()];
	    				if(!e) {
	    					e = new GVec<int>();
	    					ends.Add(je.chars(),e);
	    				}
	    				e->Add(graphnode->nodeid);
	    			}
	    			njs++;
	    		}

	    		if(pos<endbundle) { // there is still place for another node in this bundle (I might put a limit of length here for the graphnode -> because otherwise one can assume this is just a pre-mRNA fragment)
	    			//fprintf(stderr,"Create pos+1 node[%d][%d] %d with start %d\n",s,g,no2gnode[s][g].Count(),pos+1);
	    			CGraphnode *nextnode = create_graphnode(s,g,pos+1,endbundle,graphno,bundlenode,bundle2graph,no2gnode);
	    			graphno++;
	    			graphnode->child.Add(nextnode->nodeid); // make nextnode a child of current graphnode
	    			nextnode->parent.Add(graphnode->nodeid);// make graphnode a parent of nextnode
	    			//fprintf(stderr,"4 edge %d(%d-%d)-%d(%d-%d)\n",graphnode->nodeid,graphnode->start,graphnode->end,nextnode->nodeid,nextnode->start,nextnode->end);
	    			graphnode=nextnode;
	    			// COUNT EDGE HERE
	    			edgeno++;
	    		}
	    		else completed=true;
	    	}
	    	else if(minjunction == 1) { // found a junction end here

	    		uint pos=ejunction[nje]->end;
	    		while(nje<njunctions && ejunction[nje]->end==pos) { // read all junction ends at the current start
	    			nje++;
	    		}

	    		if(graphnode->start<pos) { // last created node starts before the position of the new node I want to create

	    			if(trim) graphnode=trimnode(s,g,refstart,pos-1,graphnode,source,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);// do something to find intermediate nodes; alternatively, I could only do this for end nodes

	    			graphnode->end=pos-1; // set end of current graphnode here
	    			CGraphnode *nextnode = create_graphnode(s,g,pos,endbundle,graphno,bundlenode,bundle2graph,no2gnode);
	    			graphno++;
	    			graphnode->child.Add(nextnode->nodeid); // make nextnode a child of current graphnode
	    			nextnode->parent.Add(graphnode->nodeid);// make graphnode a parent of nextnode
	    			//fprintf(stderr,"5 edge %d(%d-%d)-%d(%d-%d)\n",graphnode->nodeid,graphnode->start,graphnode->end,nextnode->nodeid,nextnode->start,nextnode->end);
	    			graphnode=nextnode;
	    			// COUNT EDGE HERE
	    			edgeno++;
	    		}

	    		GStr spos((int)pos);
	    		GVec<int> *e=ends[spos.chars()]; // WHY DOESN'T THIS REPEAT THE SAME THING IN CASE THE START HASN'T BEEN ADJUSTED? because nje is bigger now than the ones that end at the currentstart
	    		if(e) for(int i=0;i<e->Count();i++) {
	    			CGraphnode *node=no2gnode[s][g][e->Get(i)];
	    			node->child.Add(graphnode->nodeid);  // this node is the child of previous node
	    			graphnode->parent.Add(node->nodeid); // this node has as parent the previous node
	    			//fprintf(stderr,"6 edge %d(%d-%d)-%d(%d-%d)\n",node->nodeid,node->start,node->end,graphnode->nodeid,graphnode->start,graphnode->end);
	    			// COUNT EDGE HERE
	    			edgeno++;
	    		}
	    	}

	    } while((nje<njunctions && (ejunction[nje]->end<endbundle)) || (njs<njunctions && (junction[njs]->start<=endbundle)));

	    //fprintf(stderr,"done while completed=%d\n",completed);

	    if(!completed) { // I did not finish node --> this will be an ending node

	      if(trim) // do something to find intermediate nodes; alternatively, I could only do this for end nodes
	    	  graphnode=trimnode(s,g,refstart,endbundle,graphnode,source,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);

	      graphnode->end=endbundle;
	      // COUNT EDGE HERE (this is an edge to sink)
	      edgeno++;
	      //fprintf(stderr,"7 edge %d(%d-%d)-sink\n",graphnode->nodeid,graphnode->start,graphnode->end);
	    }

	    bundlenode=bundlenode->nextnode; // advance to next bundle
	} // end while(bundlenode!=NULL)

	sink->nodeid=graphno;
	no2gnode[s][g].Add(sink);
	graphno++;

	//fprintf(stderr,"This graph has %d nodes and %d edges\n",graphno,edgeno);
	lastgpos=graphno; // nodes are from 0 to graphno-1, so the first "available" position in GBitVec is graphno

	// now I can create the future transfrags because I know graphno
	for(int i=0;i<futuretr.Count();i+=3) {
		// add links between node and sink
		int n1=int(futuretr[i]);
		int n2=int(futuretr[i+1]);
		GBitVec trpat(graphno+edgeno);
		trpat[n1]=1;
		GVec<int> nodes;
		if(n2<0) {
			CGraphnode *node=no2gnode[s][g][n1];
			node->child.Add(sink->nodeid);
			// add node to sink transfrag
			trpat[graphno-1]=1;

			int key=edge(n1,graphno-1,graphno);
			int *pos=gpos[s][g][key];
			if(pos!=NULL) trpat[*pos]=1;
			else {
				gpos[s][g].Add(key,lastgpos);
				trpat[lastgpos]=1;
				lastgpos++;
			}

			nodes.Add(n1);
			nodes.Add(sink->nodeid);
		}
		else {
			trpat[n2]=1;
			int key=edge(n1,n2,graphno);
			int *pos=gpos[s][g][key];
			if(pos!=NULL) trpat[*pos]=1;
			else {
				gpos[s][g].Add(key,lastgpos);
				trpat[lastgpos]=1;
				lastgpos++;
			}

			nodes.cAdd(n1);
			nodes.Add(n2);
		}
		CTransfrag *tr=new CTransfrag(nodes,trpat,futuretr[i+2]);
		/*
		{ // DEBUG ONLY
			fprintf(stderr,"Add future transfrag[%d][%d]= %d with pattern",s,g,transfrag[s][g].Count());
			//printBitVec(trpat);
			fprintf(stderr,"\n");
		}
		*/
		transfrag[s][g].Add(tr);
	}

	// finished reading bundle -> now create the parents' and children's patterns
	GVec<bool> visit;
	visit.Resize(graphno,false);
	GBitVec parents(graphno+edgeno);

	//fprintf(stderr,"traverse graph[%d][%d] now with %d edges and lastgpos=%d....\n",s,g,edgeno,lastgpos);//edgeno=0;
	traverse_dfs(s,g,source,sink,parents,graphno,visit,no2gnode,transfrag,edgeno,gpos,lastgpos);
	//fprintf(stderr,"done traversing with edgeno=%d\n",edgeno);

	// delete variables created here, like e.g. ends; do I need to delete the GVec<int> elements created too?
	ends.Clear();

	return(graphno);

}

/*
// I don't use this one
int create_graph(int s,int g,CBundle *bundle,GPVec<CBundlenode>& bnode, GList<CJunction>& junction,GList<CJunction>& ejunction,GVec<CGraphinfo> **bundle2graph,
		GPVec<CGraphnode> **no2gnode,GPVec<CTransfrag> **transfrag,GVec<CTrimPoint>& trims){

	CGraphnode* source=new CGraphnode(0,0,0);
	no2gnode[s][g].Add(source);
	CGraphnode* sink=new CGraphnode();

	int njunctions=junction.Count();
	int njs=0; // index of sorted junction starts
	int nje=0; // index of sorted junction ends

	int graphno=1; // number of nodes in graph
	GHash<GVec<int> > ends;

	CBundlenode *bundlenode=bnode[bundle->startnode];

	GVec<float> futuretr;
	int nt=0; // index for the current trim

	while(bundlenode!=NULL) {

	    uint currentstart=bundlenode->start; // current start is bundlenode's start
	    uint endbundle=bundlenode->end; // initialize end with bundlenode's end for now

	    CGraphnode *graphnode=create_graphnode(s,g,currentstart,endbundle,graphno,bundlenode,bundle2graph,no2gnode); // creates a $graphno graphnode  with start at bundle start, and end at bundle end
	    graphno++;

	    int end=0;
	    while(nje<njunctions && ejunction[nje]->end<=currentstart) { // read all junction ends at or before the current start -> assuming there are any (at this point, smaller junction ends should not be relevant to this bundle/currentstart
	      if(ejunction[nje]->end==currentstart && (ejunction[nje]->strand+1) == 2*s) { // junction ends at current start and is on the same strand and not deleted
	    	  end=1;
	      }
	      nje++;
	    }

	    if(end) { // I might have nodes finishing here
	    	GStr cs((int)currentstart);
	    	GVec<int> *e=ends[cs.chars()];
	    	if(e) {
	    		for(int i=0;i<e->Count();i++) {
	    			CGraphnode *node=no2gnode[s][g][e->Get(i)];
	    			node->child.Add(graphnode->nodeid);  // this node is the child of previous node
	    			graphnode->parent.Add(node->nodeid); // this node has as parent the previous node
	    		}
	    	}
	    	else { // I haven't seen nodes before that finish here => link to source
		    	source->child.Add(graphnode->nodeid);  // this node is the child of source
		    	graphnode->parent.Add(source->nodeid); // this node has source as parent
	    	}
	    }
	    else { // this node comes from source directly
	    	source->child.Add(graphnode->nodeid);  // this node is the child of source
	    	graphnode->parent.Add(source->nodeid); // this node has source as parent
	    }


	    int completed=0;
	    do {
	    	while(nje<njunctions && ((ejunction[nje]->strand+1) != 2*s)) nje++;
	    	while(njs<njunctions && (((junction[njs]->strand+1)!= 2*s) || (junction[njs]->start<currentstart))) njs++; // junctions that start before the current graphnode and I haven't seen them before are part of a different bundle

	    	int minjunction = -1; // process next junction -> either a start or an ending whichever has the first position on the genome
	    	if((nje<njunctions && (ejunction[nje]->end<endbundle)) || (njs<njunctions && (junction[njs]->start<=endbundle))) {
	    		if(nje<njunctions) { // there are still junctions endings
	    			if(njs<njunctions) { // there are still junctions starting
	    				minjunction = junction[njs]->start > ejunction[nje]->end ? 1 : 0;
	    			}
	    			else minjunction = 1;
	    		}
	    		else minjunction = 0;
	    	}

	    	if(minjunction == 0 ) { // found a start junction here

	    		while(nt < trims.Count() && trims[nt].pos<=junction[njs]->start) {
	    			if(trims[nt].pos>=graphnode->start)
	    				graphnode=add_trim_to_graph(s,g,junction[njs]->start,trims[nt],graphnode,source,sink,futuretr,graphno,bundlenode,bundle2graph,no2gnode);
	    			nt++;
	    		}

	    		// if no trimming required just set the end of the node
	    		graphnode->end=junction[njs]->start; // set the end of current graphnode to here; introduce smaller nodes if trimming is activated

	    		uint pos=junction[njs]->start;
	    		while(njs<njunctions && junction[njs]->start==pos ) { // remember ends here
	    			if((junction[njs]->strand+1) == 2*s) {
	    				GStr je((int)junction[njs]->end);
	    				GVec<int> *e=ends[je.chars()];
	    				if(!e) {
	    					e = new GVec<int>();
	    					ends.Add(je.chars(),e);
	    				}
	    				e->Add(graphnode->nodeid);
	    			}
	    			njs++;
	    		}

	    		if(pos<endbundle) { // there is still place for another node in this bundle (I might put a limit of length here for the graphnode -> because otherwise one can assume this is just a pre-mRNA fragment)
	    			//print STDERR "Create nextnode:",$pos+1,"-",$endbundle,"\n";
	    			CGraphnode *nextnode = create_graphnode(s,g,pos+1,endbundle,graphno,bundlenode,bundle2graph,no2gnode);
	    			graphno++;
	    			graphnode->child.Add(nextnode->nodeid); // make nextnode a child of current graphnode
	    			nextnode->parent.Add(graphnode->nodeid);// make graphnode a parent of nextnode
	    			graphnode=nextnode;
	    		}
	    		else completed=1;
	    	}
	    	else if(minjunction == 1) { // found a junction end here

	    		uint pos=ejunction[nje]->end;
	    		while(nje<njunctions && ejunction[nje]->end==pos) { // read all junction ends at the current start
	    			nje++;
	    		}

	    		if(graphnode->start<pos) { // last created node starts before the position of the new node I want to create

		    		while(nt < trims.Count() && trims[nt].pos<pos) {
		    			if(trims[nt].pos>=graphnode->start)
		    				graphnode=add_trim_to_graph(s,g,pos-1,trims[nt],graphnode,source,sink,futuretr,graphno,bundlenode,bundle2graph,no2gnode);
		    			nt++;
		    		}

	    			graphnode->end=pos-1; // set end of current graphnode here
	    			CGraphnode *nextnode = create_graphnode(s,g,pos,endbundle,graphno,bundlenode,bundle2graph,no2gnode);
	    			graphno++;
	    			graphnode->child.Add(nextnode->nodeid); // make nextnode a child of current graphnode
	    			nextnode->parent.Add(graphnode->nodeid);// make graphnode a parent of nextnode
	    			graphnode=nextnode;
	    		}

	    		// print STDERR "Check how many nodes end here at $pos\n";
	    		GStr spos((int)pos);
	    		GVec<int> *e=ends[spos.chars()];
	    		if(e) for(int i=0;i<e->Count();i++) {
	    			CGraphnode *node=no2gnode[s][g][e->Get(i)];
	    			node->child.Add(graphnode->nodeid);  // this node is the child of previous node
	    			graphnode->parent.Add(node->nodeid); // this node has as parent the previous node
	    		}
	    	}

	    } while((nje<njunctions && (ejunction[nje]->end<endbundle)) || (njs<njunctions && (junction[njs]->start<=endbundle)));


	    if(!completed) {

	    	while(nt < trims.Count() && trims[nt].pos<=endbundle) {
	    		if(trims[nt].pos>=graphnode->start)
	    			graphnode=add_trim_to_graph(s,g,endbundle,trims[nt],graphnode,source,sink,futuretr,graphno,bundlenode,bundle2graph,no2gnode);
    			nt++;
    		}

	      graphnode->end=endbundle;
	    }

	    bundlenode=bundlenode->nextnode; // advance to next bundle
	} // end while(bundlenode!=NULL)

	sink->nodeid=graphno;
	no2gnode[s][g].Add(sink);
	graphno++;

	// now I can create the future transfrags because I know graphno
	for(int i=0;i<futuretr.Count();i+=3) {
		// add links between node and sink
		int n=int(futuretr[i]);
		GBitVec trpat(1+(graphno+1)*graphno/2);
		trpat[n]=1;
		GVec<int> nodes;
		if(futuretr[i+1]) {
			CGraphnode *node=no2gnode[s][g][n];
			node->child.Add(sink->nodeid);
			// add node to sink transfrag
			trpat[graphno-1]=1;
			trpat[edge(n,graphno-1,graphno)]=1;
			nodes.Add(n);
			nodes.Add(sink->nodeid);
		}
		else {
			trpat[0]=1;
			trpat[edge(0,n,graphno)]=1;
			nodes.cAdd(0);
			nodes.Add(n);
		}
		CTransfrag *tr=new CTransfrag(nodes,trpat,futuretr[i+2]);
		transfrag[s][g].Add(tr);
	}

	// finished reading bundle -> now create the parents' and children's patterns
	GVec<bool> visit;
	visit.Resize(graphno,false);
	GBitVec parents(1+(graphno+1)*graphno/2);
	traverse_dfs(s,g,source,sink,parents,graphno,visit,no2gnode,transfrag);

	// delete variables created here, like e.g. ends; do I need to delete the GVec<int> elements created too?
	ends.Clear();

	return(graphno);

}
*/

void get_read_pattern(float readcov,GBitVec& pattern0,GBitVec& pattern1,int *rgno, float *rprop,GVec<int> *rnode,GList<CReadAln>& readlist,int n,
		GVec<int> *readgroup,GVec<int>& merge,GVec<int> *group2bundle,GVec<CGraphinfo> **bundle2graph,GVec<int> *graphno,GVec<int> *edgeno,GIntHash<int> **gpos,
		GPVec<CGraphnode> **no2gnode) {
		//uint readedge,int *rbnode) {

	int lastgnode[2]={-1,-1}; // lastgnode[0] is for - strand; [1] is for + strand -> I need these in order to add the edges to the read pattern; check this: if it's not correct than storage was wrong!
	int ncoord=readlist[n]->segs.Count();

	int k[2]={0,0}; // need to keep track of coordinates already added to coverages of graphnodes
    bool valid[2]={true,true};
    for(int s=0;s<2;s++) if(!rprop[s]) valid[s]=false;

    for(int i=0;i<readgroup[n].Count();i++) // how can a read be associated to multiple groups? ---> I guess if it is spliced
    	if(valid[0] || valid[1]) { // there are still stranded bundles associated with the read
    		int gr=readgroup[n][i];
    		while(merge[gr]!=gr) gr=merge[gr];
    		for(int s=0;s<2;s++)
    			if(valid[s]) {
    				int bnode=group2bundle[2*s][gr];
    				if(bnode>-1 && bundle2graph[s][bnode].Count()) { // group $id has a bundle node associated with it and bundle was processed
    					int nbnode=bundle2graph[s][bnode].Count();
    					int j=0;
    					while(j<nbnode && k[s]<ncoord) {
    						int ngraph=bundle2graph[s][bnode][j].ngraph;
    						int gnode=bundle2graph[s][bnode][j].nodeno;
    						CGraphnode *node=no2gnode[s][ngraph].Get(gnode);
							// compute graphnode coverage by read here (assuming they come in order on the genomic line)
    						bool intersect=false;
    						while(k[s]<ncoord) {
    							int bp = readlist[n]->segs[k[s]].overlapLen(node);
    							if(bp) {
    								intersect=true;
    								//if(readedge>=node->start && readedge<=node->end) rbnode[s]=bnode;
    								/*
    								if(!readlist[n]->strand) { // if read is unstranded then only a certain proportion of it should go to the node coverage
        								rprop[s]=group[gr]->neg_prop; // unspliced read should belong to one group only
        								if(s) rprop[s]=1-rprop[s];
    								}
    								*/
    								//fprintf(stderr,"update cov of node %d for multi=%g and readcov=%g read=%d\n",node->nodeid,multi,bp*readcov,n);
    								node->cov+=rprop[s]*bp*readcov;
    								if(readlist[n]->segs[k[s]].end<=node->end) k[s]++;
    								else break;
				  				}
    							else break;
    						}

    						if(intersect) { // read intersects gnode
    							if(lastgnode[s]>-1 && gnode!=lastgnode[s]) { // this is not the first time I see a gnode for the read
									// need to update the edge pattern
    								int min=lastgnode[s];
    								int max=gnode;
    								if(min>gnode) {
    									max=min;
    									min=gnode;
    								}

    								int *pos=gpos[s][ngraph][edge(min,max,graphno[s][ngraph])];
    								if(pos!=NULL) {
    									if(s) pattern1[*pos]=1; // added edge from previous gnode to current one for the read
    									else pattern0[*pos]=1; // added edge from previous gnode to current one for the read
    								}
    								//else GError("s=%d g=%d edge=%d Found read %d(%d-%d) with edge between nodes %d(%d-%d)-%d(%d-%d) that is not in the graph with %d nodes and %d edges!\n",s,ngraph,edge(min,max,graphno[s][ngraph]),n,readlist[n]->start,readlist[n]->end,min,no2gnode[s][ngraph].Get(min)->start,no2gnode[s][ngraph].Get(min)->end,max,no2gnode[s][ngraph].Get(max)->start,no2gnode[s][ngraph].Get(max)->end,graphno[s][ngraph],edgeno[s][ngraph]);
    								rnode[s].Add(gnode);
    							}
    							else { // first time considering read
    								rgno[s]=ngraph;
    								rnode[s].Add(gnode);
    							}
    							lastgnode[s]=gnode;
    							if(s) {
    								pattern1.resize(graphno[s][ngraph]+edgeno[s][ngraph]);
    								pattern1[gnode]=1; // here I could remember kids as well to speed things up
    							}
    							else {
    								pattern0.resize(graphno[s][ngraph]+edgeno[s][ngraph]);
    								pattern0[gnode]=1; // here I could remember kids as well to speed things up
    							}
    						} // end if(intersect
    						j++;
    					} // end while(j<nbnode && k<ncoord)
    				} // end if(bnode>-1 && bundle2graph[s][bnode].Count())
    				else { // I don't need to process this strand for the read anymore because bundle was not processed or read doesn't belong to that stranded bundle
    					valid[s]=false;
    					continue;
    				}
    			} // endif(valid[s])
    	} // end if(valid[0] || valid[2])

}

void free_treepat(CTreePat *t){
	if(t) {
		for(int i=0;i<t->childno;i++) free_treepat(t->nextpat[i]);
		GFREE(t->nextpat);
		delete t;
	}
}

CTreePat *construct_treepat(int gno, GIntHash<int>& gpos,GPVec<CTransfrag>& transfrag) {

	// create root CTreePat first
	CTreePat *root=new CTreePat(0,gno-1); // if links from source to nodes are desired this should be changed

	// now construct all child CTreePat's
	for(int t=0;t<transfrag.Count();t++)
		if(transfrag[t]->nodes[0]){ // don't include transfrags from source -> not needed
			CTreePat *tree=root;
			int m=0; // previous node in pattern that was set in pattern
			for(int n=1;n<gno;n++){
				if(transfrag[t]->pattern[n]) {
					CTreePat *child;
					if(m) { // there is a node m that was seen before
						int *pos=gpos[edge(m,n,gno)];
						if(pos && transfrag[t]->pattern[*pos]) // there is an edge between m and n
							child=tree->settree(gno-1-m+n-m-1,n,2*(gno-n-1));
						else child=tree->settree(n-m-1,n,2*(gno-n-1));
					}
					else // this is the root tree
						child=tree->settree(n-1,n,2*(gno-n-1));
					tree=child;
					m=n;
				}
			}
			tree->tr=transfrag[t];
		}

	return(root);
}

void print_pattern(CTreePat *tree,GStr& pattern,int gno) {
	if(!tree) return;
	pattern+=tree->nodeno;
	if(tree->tr) {
		fprintf(stderr,"pat: %s %f\n",pattern.chars(),tree->tr->abundance);
	}
	for(int i=0;i<tree->childno;i++)
		if(tree->nextpat[i]) {
			GStr tmppat(pattern.chars());
			if(i<gno-1-tree->nodeno) tmppat+=".";
			else tmppat+="-";
			print_pattern(tree->nextpat[i],tmppat,gno);
		}
}

CTransfrag *findtrf_in_treepat(int gno,GIntHash<int>& gpos,GVec<int>& node,GBitVec& pattern,CTreePat *tr2no) { // doesn't work for patterns including source node

	if(!tr2no) return(NULL);

	CTreePat *tree=tr2no;
	for(int n=0;n<node.Count();n++) {
		if(n) { // not the first node in pattern
			int *pos=gpos[edge(node[n-1],node[n],gno)];
			if(pos && pattern[*pos]) // there is an edge between node[n-1] and node[n]
				tree=tree->nextpat[gno-1-node[n-1]+node[n]-node[n-1]-1];
			else tree=tree->nextpat[node[n]-node[n-1]-1];
		}
		else tree=tree->nextpat[node[n]-1];
		if(!tree) return(NULL);
	}

	return(tree->tr);
}

/*
bool istrf_in_treepat(int gno,GHash<int>& gpos,GVec<int>& node,GBitVec& pattern,CTreePat *tr2no) {

	if(!tr2no) return(false);

	CTreePat *tree=tr2no;
	for(int n=0;n<node.Count();n++) {
		if(n) { // not the first node in pattern
			bool isedge=false;
			GStr spos(edge(node[n-1],node[n],gno));
			const int *pos=gpos[spos.chars()];
			if(pos && pattern[*pos]) isedge=true;
			if(isedge) // there is an edge between node[n-1] and node[n]
				tree=tree->nextpat[gno-1-node[n-1]+node[n]-node[n-1]-1];
			else tree=tree->nextpat[node[n]-node[n-1]-1];
		}
		else tree=tree->nextpat[node[n]-1];
		if(!tree) return(false);
	}

	return(true);
}
*/

void update_abundance(int s,int g,int gno,GIntHash<int>&gpos,GBitVec& pattern,float abundance,GVec<int>& node,GPVec<CTransfrag> **transfrag,
		CTreePat ***tr2no){

	CTransfrag *t=findtrf_in_treepat(gno,gpos,node,pattern,tr2no[s][g]);
	if(!t) { // t is NULL
		t=new CTransfrag(node,pattern,0);

		/*
		{ // DEBUG ONLY
			fprintf(stderr,"Add update transfrag[%d][%d]=%d and pattern",s,g,transfrag[s][g].Count());
			//printBitVec(pattern);
			fprintf(stderr,"\n");
		}
		*/

		transfrag[s][g].Add(t);

		// node.Sort() : nodes should be sorted; if they are not then I should update to sort here
		CTreePat *tree=tr2no[s][g];
		for(int n=0;n<node.Count();n++) {
			CTreePat *child;
			if(n) { // not the first node in pattern
				int *pos=gpos[edge(node[n-1],node[n],gno)];
				if(pos && pattern[*pos]) // there is an edge between node[n-1] and node[n]
					child=tree->settree(gno-1-node[n-1]+node[n]-node[n-1]-1,node[n],2*(gno-node[n]-1));
				else child=tree->settree(node[n]-node[n-1]-1,node[n],2*(gno-node[n]-1));
			}
			else child=tree->settree(node[n]-1,node[n],2*(gno-node[n]-1));
			tree=child;
		}
		tree->tr=t;
	}
	t->abundance+=abundance;

}

void get_fragment_pattern(GList<CReadAln>& readlist,int n, int np,float readcov,GVec<int> *readgroup,GVec<int>& merge,
		GVec<int> *group2bundle,GVec<CGraphinfo> **bundle2graph,GVec<int> *graphno,GVec<int> *edgeno, GIntHash<int> **gpos,GPVec<CGraphnode> **no2gnode,
		GPVec<CTransfrag> **transfrag,CTreePat ***tr2no,GPVec<CGroup> &group) {

	GBitVec rpat[2];
	int rgno[2]={-1,-1};
	GVec<int> rnode[2];
	float rprop[2]={1,1};
	//bool goodfrag=false;

	// compute proportions of read associated to strands
	if(readlist[n]->nh && !readlist[n]->strand && np>-1 && readlist[np]->nh && !readlist[np]->strand) { // both reads are unstranded
		int gr1=readgroup[n][0]; // read is unstranded => it should belong to one group only
		while(merge[gr1]!=gr1) gr1=merge[gr1];
		int gr2=readgroup[n][0]; // read is unstranded => it should belong to one group only
		while(merge[gr2]!=gr2) gr2=merge[gr2];
		rprop[0]=(group[gr1]->neg_prop+group[gr2]->neg_prop)/2;
		rprop[1]=1-rprop[0];
		/*
		if(gr1==gr2) { // both reads are unstranded -> we might be able to determine fragment length
			goodfrag=true;
		}
		*/
	}
	else {
		if(readlist[n]->nh) {
			if(!readlist[n]->strand) { // the paired read is not present otherwise it would have the same strand from add_read_to_group
				int gr=readgroup[n][0]; // read is unstranded => it should belong to one group only
				while(merge[gr]!=gr) gr=merge[gr];
				rprop[0]=group[gr]->neg_prop;
				rprop[1]=1-rprop[0];
			}
			else {
				if(readlist[n]->strand==-1) rprop[1]=0;
				else rprop[0]=0;
			}
		}
		else if(np>-1 && readlist[np]->nh) { // readlist[n] is deleted
			if(!readlist[np]->strand) { // the paired read is not present otherwise it would have the same strand from add_read_to_group
				int gr=readgroup[np][0]; // read is unstranded => it should belong to one group only
				while(merge[gr]!=gr) gr=merge[gr];
				rprop[0]=group[gr]->neg_prop;
				rprop[1]=1-rprop[0];
			}
			else {
				if(readlist[np]->strand==-1) rprop[1]=0;
				else rprop[0]=0;
			}
		}
	}

	//int bnode1[2]={-1,-1};
	if(readlist[n]->nh) {
		get_read_pattern(readcov,rpat[0],rpat[1],rgno,rprop,rnode,readlist,n,readgroup,merge,group2bundle,bundle2graph,graphno,edgeno,gpos,no2gnode);
		//get_read_pattern(readcov,rpat[0],rpat[1],rgno,rprop,rnode,readlist,n,readgroup,merge,group2bundle,bundle2graph,graphno,edgeno,gpos,no2gnode,readlist[n]->end,bnode1);
		//bdata->num_reads+=readcov;
		//bdata->sumreads+=readcov*readlist[n]->len;
	}

	GBitVec ppat[2];
	int pgno[2]={-1,-1};
	GVec<int> pnode[2];
	//int bnode2[2]={-1,-1};
	// get pair pattern if pair exists and it hasn't been deleted
	if(np>-1 && readlist[np]->nh) {
		get_read_pattern(readcov,ppat[0],ppat[1],pgno,rprop,pnode,readlist,np,readgroup,merge,group2bundle,bundle2graph,graphno,edgeno,gpos,no2gnode);
		//get_read_pattern(readcov,ppat[0],ppat[1],pgno,rprop,pnode,readlist,np,readgroup,merge,group2bundle,bundle2graph,graphno,edgeno,gpos,no2gnode,readlist[np]->start,bnode2);

		/* computing some statistics -> might reconsider in the future
		bdata->num_reads+=readcov;
		bdata->sumreads+=readcov*readlist[np]->len;

		if((bnode1[0]!=-1 && bnode1[0]==bnode2[0])||(bnode1[1]!=-1 && bnode1[1]==bnode2[1])) {
			goodfrag=true;
			if(readlist[n]->segs.Count()>1 && readlist[np]->start<readlist[n]->end) goodfrag=false;
		}

		// check if I can get fragment length here
		if(goodfrag) {
			bdata->sumfrag+=(double(readlist[n]->len+readlist[np]->len+readlist[np]->start-1)-(double)readlist[n]->end)*readcov;
			bdata->num_frag+=readcov;
			fprintf(stderr,"fraglen=%d where n=%d np=%d len1=%d len2=%d, end=%d start=%d readcov=%g sumfrag=%g numfrag=%g\n",readlist[n]->len+readlist[np]->len+readlist[np]->start-1-readlist[n]->end,n,np,readlist[n]->len,readlist[np]->len,readlist[n]->end,readlist[np]->start,readcov,bdata->sumfrag,bdata->num_frag);
		}
		*/
	}

	//if(readlist[n]->nh || (np>-1 && readlist[np]->nh)) bdata->num_cov+=readcov;

	for(int s=0;s<2;s++){
		if(rgno[s]>-1) { // read is valid (has pattern) on strand s
			if(pgno[s]>-1) { // pair is also valid => fragment is valid: check if there are conflicts
				// check if there are conflicts between patterns and graph -> could this happen if there are no reads covering  a gap in bundle? --> yes
				bool conflict=true;
				if(rgno[s]==pgno[s]) { // read and pair belong to the same graph
					// check if there is a conflict of patterns
					CGraphnode *gnode=no2gnode[s][rgno[s]][pnode[s][0]];
					GBitVec conflictpattn=gnode->parentpat;
					conflictpattn[pnode[s][0]]=1;

					if((conflictpattn & rpat[s])==rpat[s]) { // there isn't a conflict -> pair parents should contain read pattern
						conflict=false;
						int i=0;
						if(pnode[s][0]==rnode[s].Last()) // read and pair share a node
							i++;
						while(i<pnode[s].Count()) { rnode[s].Add(pnode[s][i]);i++;}
						rpat[s]=rpat[s]|ppat[s];
						update_abundance(s,rgno[s],graphno[s][rgno[s]],gpos[s][rgno[s]],rpat[s],rprop[s]*readcov,rnode[s],transfrag,tr2no);
					}
				}
				if(conflict) { // update both patterns separately
					update_abundance(s,rgno[s],graphno[s][rgno[s]],gpos[s][rgno[s]],rpat[s],rprop[s]*readcov,rnode[s],transfrag,tr2no);
					update_abundance(s,pgno[s],graphno[s][pgno[s]],gpos[s][pgno[s]],ppat[s],rprop[s]*readcov,pnode[s],transfrag,tr2no);
				}
			}
			else { // pair has no valid pattern
				update_abundance(s,rgno[s],graphno[s][rgno[s]],gpos[s][rgno[s]],rpat[s],rprop[s]*readcov,rnode[s],transfrag,tr2no);
			}
		}
		else // read has no valid pattern but pair might
			if(pgno[s]>-1) {
				update_abundance(s,pgno[s],graphno[s][pgno[s]],gpos[s][pgno[s]],ppat[s],rprop[s]*readcov,pnode[s],transfrag,tr2no);
			}
	}

}

void settrf_in_treepat(CTransfrag *t,int gno,GIntHash<int>& gpos,GVec<int>& node,GBitVec& pattern,CTreePat *tr2no) {

	if(!tr2no) return;

	CTreePat *tree=tr2no;
	for(int n=0;n<node.Count();n++) {
		if(n) { // not the first node in pattern
			int *pos=gpos[edge(node[n-1],node[n],gno)];
			if(pos && pattern[*pos]) // there is an edge between node[n-1] and node[n]
				tree=tree->nextpat[gno-1-node[n-1]+node[n]-node[n-1]-1];
			else tree=tree->nextpat[node[n]-node[n-1]-1];
		}
		else tree=tree->nextpat[node[n]-1];
		if(!tree) return;
	}
	tree->tr=t;
}

void eliminate_transfrags_under_thr(int gno,GIntHash<int>& gpos,GPVec<CTransfrag>& transfrag, CTreePat *tr2no,float threshold) {

	for(int t=transfrag.Count()-1;t>=0;t--)
		if(transfrag[t]->abundance<threshold && transfrag[t]->nodes[0] && transfrag[t]->nodes.Last()<gno-1) { // need to delete transfrag that doesn't come from source or ends at sink
			settrf_in_treepat(NULL,gno,gpos,transfrag[t]->nodes,transfrag[t]->pattern,tr2no); // this should be eliminated if I want to store transcripts from 0 node
			transfrag.Exchange(t,transfrag.Count()-1);
			transfrag.Delete(transfrag.Count()-1);
		}

	while(transfrag.Count()>max_trf_number) {
		threshold++;
		for(int t=transfrag.Count()-1;t>=0;t--)
			if(transfrag[t]->abundance<threshold && transfrag[t]->nodes[0] && transfrag[t]->nodes.Last()<gno-1) { // need to delete transfrag that doesn't come from source
				settrf_in_treepat(NULL,gno,gpos,transfrag[t]->nodes,transfrag[t]->pattern,tr2no); // this should be eliminated if I want to store transcripts from 0 node
				transfrag.Exchange(t,transfrag.Count()-1);
				transfrag.Delete(transfrag.Count()-1);
			}
	}

}

bool conflict(int &i,int node,GVec<int>& trnode,int n,GPVec<CGraphnode>& no2gnode,GBitVec& trpat,int gno, GIntHash<int>& gpos) {

  while(i<n && node>trnode[i]) i++;

  if(i<n && node==trnode[i]) return(false);

  int node1=trnode[i-1];

  if(i==n) { // transcript is at the end; $node > $node1
    if(no2gnode[node1]->childpat[node]) // node is among children of node1 -> should be equivalent to testing that node1 is among parents of node
    	return(false);
    return(true);
  }

  int node2=trnode[i]; // node1 < node < node2

  if(no2gnode[node1]->childpat[node] && no2gnode[node]->childpat[node2]) {
	  int *pos=gpos[edge(node1,node2,gno)];
	  if(pos && trpat[*pos]) return(true);
	  else return(false);
  }

  return(true);
}

void binary_insert_link(int n2, float capacity,GVec<int>& link, GVec<float>& val, int first, int last) {

	while(first<=last) {
		int mid=(first + last) / 2;  // compute mid point.
		if(capacity>val[mid]) first=mid+1;
		else if(capacity<val[mid]) last=mid-1;
		else { // found a good position to inset
			first=mid;
			break;
		}
	}

	int pos=first;
	val.idxInsert(pos,capacity);
	link.idxInsert(pos,n2);

}

void binary_insert_trf_to_node(int t, GVec<int>& trf,int first,int last) {

	while(first<=last) {
		int mid=(first + last) / 2;  // compute mid point.
		if(t>trf[mid]) first=mid+1;
		else if(t<trf[mid]) last=mid-1;
		else return; // transcript already inserted
	}

	int pos=first;

	trf.idxInsert(pos,t);

}

void assign_incomplete_trf_to_nodes(int t,int n1, int n2,GPVec<CGraphnode>& no2gnode){

	// a better version of this function could take into account also the fragments inner length

	for(int i=n1+1;i<n2;i++) {
		if(no2gnode[n1]->childpat[i] && no2gnode[i]->childpat[n2]) { // i is a child of n1, and n2 is a child of i
			binary_insert_trf_to_node(t,no2gnode[i]->trf,0,no2gnode[i]->trf.Count()-1);
		}
	}

	/* this version was taking a very long time
	for(int i=0;i<no2gnode[n1]->child.Count();i++) { // for all children of n1
		int child=no2gnode[n1]->child[i];
		if(child!=n2 && no2gnode[child]->childpat[n2]) { // n2 is child of n1's child
			binary_insert_trf_to_node(t,no2gnode[child]->trf,0,no2gnode[child]->trf.Count()-1);
			assign_incomplete_trf_to_nodes(t,child,n2,no2gnode);
		}
	}
	*/

}

inline int comptbl_pos(int t1,int t2,int n){
	return(t2+t1*(2*n-t1-1)/2);
}

void process_transfrags(int gno,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,CTreePat *tr2no,
		GBitVec& compatible,int edgeno,GIntHash<int> &gpos) {

	/*
	{ // DEBUG ONLY
		printTime(stderr);
		fprintf(stderr,"There are %d transfrags before clean up\n",transfrag.Count());
	}
	*/

	// eliminate transfrags below threshold (they represent noise) if they don't come from source
	eliminate_transfrags_under_thr(gno,gpos,transfrag,tr2no,trthr);

	// introduce "fake" transcripts as holders to adjust the abundances, if the EM algorithm is not used -> I don't need to to this for source and sink

	//if(!EM) {
	if(0) {
		GVec<float> f_init(gno, float(0));
		GVec< GVec<float> > n(gno, &f_init, false);
		//GVec<float> n[gno]; // abundances of all transfrags passing through node and leaving node
		GVec< GVec<float> > f(gno, &f_init, false);
		//GVec<float> f[gno]; // abundances of all transfrags leaving node
		//GVec< GVec<int> > u(gno, new GVec<int>(gno,0) );
		//GVec<int> u[gno]; // transfrags that are uncommited to a specific child; this is a more complicated issue -> ignore for now
		/*
		for(int i=1;i<gno-1;i++) {
			n[i].Resize(gno,0);
			f[i].Resize(gno,0);
		}
        */
		for(int t=0;t<transfrag.Count();t++)
			if(transfrag[t]->nodes[0] && transfrag[t]->nodes.Count()) { // transfrag has more than one node, and doesn't start at source
				for(int i=0;i<transfrag[t]->nodes.Count()-1;i++) { // for all nodes in transfrags that are not last
					n[transfrag[t]->nodes[i]][transfrag[t]->nodes[i+1]]+=transfrag[t]->abundance;
					if(!i) f[transfrag[t]->nodes[i]][transfrag[t]->nodes[i+1]]+=transfrag[t]->abundance;
				}
			}

		/* if I were able to know the uncommited u for each node than I could distribute them
		for(int t=0;t<u.Count();t++) { // process the uncommited
				GVec<int> c;
				float sum=0;
				for(int k=0;k<nchild;k++) {
					if(no2gnode[transfrag[u[t]]->nodes[1]]->parentpat[no2gnode[i]->child[k]]) { // child k is parent of second node in transfrag
						c.Add(k);
						sum+=n[k];
					}
				}
				// if I am worried about rounding errors:
				if(c.Count()>1) {
					if(sum) for(int k=0;k<c.Count();k++)
						f[c[k]]+=(n[c[k]]*transfrag[u[t]]->abundance)/sum;
				}
				else f[c[0]]+=transfrag[u[t]]->abundance;
				// otherwise:
				if(sum) for(int k=0;k<c.Count();k++) {
					f[c[k]]+=(n[c[k]]*transfrag[u[t]]->abundance)/sum;
				}
			}
		*/

		for(int i=1;i<gno-1;i++) { // for all nodes except source and sink
			int nchild=no2gnode[i]->child.Count();
			if(nchild==1) { // if an only child I don't need to solve the linear system of equations
				int c=no2gnode[i]->child[0];
				if(!f[i][c] && n[i][c]) { // one child that has no out links but there are transcripts going through it
					GBitVec trpat(gno+edgeno);
					trpat[i]=1;
					trpat[c]=1;
					int *pos=gpos[edge(i,c,gno)];
					if(!pos) // this should not happen -> this is an error
						GError("Link between nodes %d and %d was not properly created\n",i,c);
					trpat[*pos]=1;
					GVec<int> nodes;
					nodes.Add(i);
					nodes.Add(c);
					CTransfrag *tr=new CTransfrag(nodes,trpat,n[i][c],false);
					transfrag.Add(tr);
				}
			}
			else {
				// check if I should continue processing -> don't process if some n[i]'s are 0, or if they are in the same proportions as f[i]'s
				bool contprocess=true;
				bool diff=false;
				int j=0;
				int minj=0;
				float minf=f[i][no2gnode[i]->child[0]];
				while(contprocess && j<nchild) {
					if(n[i][no2gnode[i]->child[j]]==0) contprocess=false;
					else {
						if(!diff && n[i][no2gnode[i]->child[j]]!=f[i][no2gnode[i]->child[j]]) diff=true;
						if(f[i][no2gnode[i]->child[j]]<minf) {
							minj=j;
							minf=f[i][no2gnode[i]->child[j]];
						}
						j++;
					}
				}
				if(contprocess && diff) { // abundances need to be adjusted
					bool strict=true;
					if(f[i][no2gnode[i]->child[minj]]>0) strict=false;
					// first find xj
					float xj=0;
					for(int k=0;k<nchild;k++) if(k!=minj) {
						float max=n[i][no2gnode[i]->child[minj]]*f[i][no2gnode[i]->child[k]]/n[i][no2gnode[i]->child[k]]-
								f[i][no2gnode[i]->child[minj]];
						if(max>xj) {
							xj=max;
							if(f[i][no2gnode[i]->child[k]]) strict=false;
							else strict=true;
						}
						else if(max==xj && !strict && !f[i][no2gnode[i]->child[k]]) strict=true;
					}
					if(strict) xj+=trthr;
					for(int k=0;k<nchild;k++) {
						float xk=0;
						if(k!=minj) xk=n[i][no2gnode[i]->child[k]]*(f[i][no2gnode[i]->child[minj]]+xj)/n[i][no2gnode[i]->child[minj]]-f[i][no2gnode[i]->child[k]];
						else xk=xj;
						if(xk) {
							GBitVec trpat(gno+edgeno);
							trpat[i]=1;
							trpat[no2gnode[i]->child[k]]=1;
							int *pos=gpos[edge(i,no2gnode[i]->child[k],gno)];
							if(!pos) { // this should not happen -> this is an error
								GError("Link between nodes %d and %d was not properly created\n",i,no2gnode[i]->child[k]);
							}
							trpat[*pos]=1;
							GVec<int> nodes;
							nodes.Add(i);
							nodes.Add(no2gnode[i]->child[k]);
							CTransfrag *tr=new CTransfrag(nodes,trpat,xk,false);
							transfrag.Add(tr);
						}
					}
				}
			}
		}
	}

	// sort transfrag with smallest being the one that has the most nodes, and ties are decided by the abundance (largest abundance first); last transfrags all have 1 node
	transfrag.Sort(trCmp);

	// in perl I was setting one node and two node edge' patterns here

	/*
	{ // DEBUG ONLY
		printTime(stderr);
		fprintf(stderr,"There are %d transfrags that remained\n",transfrag.Count());
	}
	*/

	// create compatibilities
	for(int t1=0;t1<transfrag.Count();t1++) { // transfrags are processed in increasing order -> important for the later considerations

		// update nodes
		int n1=transfrag[t1]->nodes.Count();

		if(n1>1) { // add transfrag to nodes' in and out; if a transfrag only has one node then it is node added to a node; I might want to change this for the computation of fpkm's
			for(int n=0;n<n1;n++) { // for all nodes in transfrag

				if(n && n<transfrag[t1]->nodes.Count()-1) {// not first or last node
					// add t1 to in and out of node
					no2gnode[transfrag[t1]->nodes[n]]->trf.Add(t1);

					// check if transfrag t1 is incomplete between node[n-1] and node [n]
					int *pos=gpos[edge(transfrag[t1]->nodes[n-1],transfrag[t1]->nodes[n],gno)];
					if(!pos || !transfrag[t1]->pattern[*pos]) // there is no edge between node[n-1] and node[n]
						assign_incomplete_trf_to_nodes(t1,transfrag[t1]->nodes[n-1],transfrag[t1]->nodes[n],no2gnode); 	// this ensures that I still have compatibilities when going through a certain path: !!! THIS IS NEW COMPARED TO PERL
				}
				else if(n) { // last but not first node
					// add t1 to in of node
					no2gnode[transfrag[t1]->nodes[n]]->trf.Add(t1);

					// check if transfrag t1 is incomplete between node[n-1] and node [n]
					int *pos=gpos[edge(transfrag[t1]->nodes[n-1],transfrag[t1]->nodes[n],gno)];
					if(!pos || !transfrag[t1]->pattern[*pos]) // there is no edge between node[n-1] and node[n]
						assign_incomplete_trf_to_nodes(t1,transfrag[t1]->nodes[n-1],transfrag[t1]->nodes[n],no2gnode); 	// this ensures that I still have compatibilities when going through a certain path: !!! THIS IS NEW COMPARED TO PERL
				}
				else { // first node -> only add transfrag to out of node
					no2gnode[transfrag[t1]->nodes[n]]->trf.Add(t1);
				}
			}
		}
		/*
		else { // this transcript is included completely in node
			no2gnode[transfrag[t1]->nodes[0]]->frag+=transfrag[t1]->abundance;
		}
		*/

		// add t1 to t1 compatibility
		bool comp=true;
		compatible[comptbl_pos(t1,t1,transfrag.Count())]=1;
		for(int t2=t1+1;t2<transfrag.Count();t2++) {
			// here check compatibility between t1 and t2;
			int n2=transfrag[t2]->nodes.Count();
			int i1=0;
			int i2=0;
			comp=true;
			while(i1<n1 && i2<n2) {
				if(transfrag[t1]->nodes[i1]==transfrag[t2]->nodes[i2]) {
					if(i1==n1-1 || i2==n2-1) { // one transcript finishes -> no need to check anymore
						i1=n1;i2=n2;
					}
					else { // advance the smallest one
						if(transfrag[t1]->nodes[i1+1]<transfrag[t2]->nodes[i2+1]) i1++;
						else i2++;
					}
				}
				else if(transfrag[t1]->nodes[i1]<transfrag[t2]->nodes[i2]) {
					i1++;
					if(conflict(i1,transfrag[t2]->nodes[i2],transfrag[t1]->nodes,n1,no2gnode,transfrag[t1]->pattern,gno,gpos)) {
						comp=false;
						break;
					}
				}
				else {
					i2++;
					if(conflict(i2,transfrag[t1]->nodes[i1],transfrag[t2]->nodes,n2,no2gnode,transfrag[t2]->pattern,gno,gpos)) {
						comp=false;
						break;
					}
				}
			}
			if(comp) compatible[comptbl_pos(t1,t2,transfrag.Count())]=1;
		} // end for(int t2=t1+1;t2<transfrag.Count();t2++)
	} // end for(int t1=0;t1<transfrag.Count();t1++)

	// set source-to-child transfrag abundances: optional in order not to keep these abundances too low:
	// update the abundances of the transfrags coming in from source and going to a node that doesn't have other parents than source
	CGraphnode *source=no2gnode[0];
	for(int i=0;i<source->child.Count();i++) {
		float abundance=0;
		int t0=-1;
		if(no2gnode[source->child[i]]->parent.Count()==1 && !no2gnode[source->child[i]]->parent[0]) { // source is the only parent of node
			for(int j=0;j<no2gnode[source->child[i]]->trf.Count();j++) {
				int t=no2gnode[source->child[i]]->trf[j];
				if(transfrag[t]->nodes.Last()==source->child[i]) t0=t;
				else abundance+=transfrag[t]->abundance;
			}
			if(t0>-1) { // found transfrag from source to node
				transfrag[t0]->abundance=abundance;
			}
		}
	}


}

bool is_compatible(int t1,int t2, int n,GBitVec& compatible) {
	if(t1<t2) return(compatible[comptbl_pos(t1,t2,n)]);
	else return(compatible[comptbl_pos(t2,t1,n)]);
}

/*
// I don't use this one
GVec<int> *max_compon_size(int trnumber,float &maxsize,GVec<CTrInfo>& set,GVec<bool>& compatible, GHash<CComponent>& computed) {

	// this max_compon presumes the set is always sorted according to the set.trno

	GVec<int> *result=NULL;

	for(int i=0;i<set.Count();i++) {
		float size=set[i].abundance;
		float maxagreesize=0;
		GVec<CTrInfo> agreeset;
		GStr s;
		for(int j=i+1;j<set.Count();j++) {
			if(compatible[comptbl_pos(set[i].trno,set[j].trno,trnumber)]) { // make sure that the transcripts in set are sorted to speed up things
				agreeset.Add(set[j]);
				s+=set[j].trno;
				maxagreesize+=set[j].abundance;
			}
		}

		GVec<int> *agreeresult=NULL;

		if(size>MIN_VAL && agreeset.Count() && (size+maxagreesize>maxsize)) {
			CComponent *agreecomp=computed[s.chars()];
			if(!agreecomp) {
				float agreesize=MIN_VAL;
				agreeresult=max_compon_size(trnumber,agreesize,agreeset,compatible,computed);
				agreecomp=new CComponent(agreesize,agreeresult);
				computed.Add(s.chars(),agreecomp);
			}
			else agreeresult=agreecomp->set;
			size+=agreecomp->size;
		}

		if(size>maxsize) {
			if(result) {
				delete result;
			}
			result=new GVec<int>();
			result->Add(set[i].trno);
			if(agreeresult) result->Add(*agreeresult);
			maxsize=size;
		}

	}

	return(result);
}
*/

GVec<int> *max_compon_size_with_penalty(int trnumber,float &maxsize,GVec<CTrInfo>& set,GBitVec& compatible,
		GBitVec& mark,GBitVec& removable, GHash<CComponent>& computed) {

	// this max_compon presumes the set is always sorted according to the set.trno

	GVec<int> *result=NULL;
	float penalty=0;

	for(int i=0;i<set.Count();i++) {
		float size=set[i].abundance-penalty;
		float maxagreesize=0;
		GVec<CTrInfo> agreeset;
		GStr s;
		for(int j=i+1;j<set.Count();j++) {
			if(compatible[comptbl_pos(set[i].trno,set[j].trno,trnumber)]) { // make sure that the transcripts in set are sorted to speed up things
				agreeset.Add(set[j]);
				s+=set[j].trno;
				if(set[j].abundance) s+=' '; else s+='.';
				if(set[j].penalty) s+=' ';else s+='.';
				maxagreesize+=set[j].abundance;
			}
			else if(mark[set[j].trno]) { // set[j] is a special interest transfrag that I prefer keeping
				if(removable[set[j].trno]) { // if I can remove set[j]
					size-=set[j].penalty; // I need to pay a penalty for removing j
				}
				else { //# j is not removable -> I can't continue considering i
					size=MIN_VAL;
					break;
				}
			}
		}

		GVec<int> *agreeresult=NULL;

		if(size>MIN_VAL && agreeset.Count() && (size+maxagreesize>maxsize)) {
			CComponent *agreecomp=computed[s.chars()];
			if(!agreecomp) {
				float agreesize=MIN_VAL;
				agreeresult=max_compon_size_with_penalty(trnumber,agreesize,agreeset,compatible,mark,removable,computed);
				agreecomp=new CComponent(agreesize,agreeresult);
				computed.Add(s.chars(),agreecomp);
			}
			else agreeresult=agreecomp->set;
			size+=agreecomp->size;
		}

		if(size>maxsize) {
			if(result) {
				delete result;
			}
			result=new GVec<int>();
			result->Add(set[i].trno);
			if(agreeresult) result->Add(*agreeresult);
			maxsize=size;
		}

		if(mark[set[i].trno]) { // set[i] is a special interest transfrag that I prefer keeping
			if(removable[set[i].trno]) { // if I can remove set[i],
				penalty+=set[i].penalty; // I need to pay a penalty for removing j
			}
			else break; //# i is not removable -> I can't continue considering the rest of the transcripts
		}

	}

	return(result);
}

bool onpath(GBitVec& trpattern,GVec<int>& trnode,GBitVec& pathpattern,int mini,int maxi,GPVec<CGraphnode>& no2gnode,int gno,GIntHash<int>& gpos) {

	if(trnode[0]<mini) // mini can be reached through transcript
	    if(!no2gnode[mini]->parentpat[trnode[0]])	return false;


	if(trnode.Last()>maxi) // from maxi I can reach end of transcript
	    if(!no2gnode[maxi]->childpat[trnode.Last()]) return false;

	int first=1;

	for(int i=0;i<trnode.Count();i++) {
		if(trnode[i]>=mini && trnode[i]<=maxi) {
	      if(!pathpattern[trnode[i]]) return false;
    	  int *pos;
    	  if(i) pos=gpos[edge(trnode[i-1],trnode[i],gno)];
	      if(first) {
	    	  first=0;
	    	  if(i && trnode[i]>mini && pos && trpattern[*pos])
	    		  return false;
	    	  if(i && !no2gnode[mini]->parentpat[trnode[i-1]]) return false; // I can not reach mini from previous node in transfrag
	      }
	      else if(i && pos && trpattern[*pos] && !pathpattern[*pos]) return false;
	    }
	    if(trnode[i]>maxi) {
	    	int *pos=gpos[edge(trnode[i-1],trnode[i],gno)];
	    	if(i && trnode[i-1]<maxi && pos && trpattern[*pos])
	    		return false;
	    	if(!no2gnode[maxi]->childpat[trnode[i]]) return false; // I can not reach this node from maxi
	    	break;
	    }
	}

	return true;
}

void update_transcript_to_path_back(float abundance,GVec<int>& trnode,GVec<int>& path,GVec<float>& pathincov,
		GVec<float>& pathoutcov) {

	int lastnode=trnode.Last();

	if(path[0]<trnode.Last()) // transcript ends at path[0]
		lastnode=path[0];

	if(trnode.Last()!=path.Last()) pathoutcov.Last()+=abundance;

	for(int i=path.Count()-2;i>=0;i--) {
		if(path[i]==lastnode) { // transcript ends at this node
			pathincov[i]+=abundance;
			break;
		}
		/* this case should never be true because the transcript can not start in the middle of the path
		 * else if(path[i]==trnode[0]) { // transcript starts at this node
			pathoutcov[i]+=abundance;
		}*/
		else { // middle of the path
			pathincov[i]+=abundance;
			pathoutcov[i]+=abundance;
		}
	}
}

void update_transcript_to_path_fwd(float abundance,GVec<int>& trnode,GVec<int>& path,GVec<float>& pathincov,
		GVec<float>& pathoutcov) {

	if(trnode[0]!=path.Last()) pathincov.Last()+=abundance;

	for(int i=path.Count()-2;i>=0;i--) {
		if(path[i]==trnode[0]) { // transcript ends at this node
			pathoutcov[i]+=abundance;
			break;
		}
		/* this case should never be true because the transcript can not end in the middle of the path
		 * else if(path[i]==trnode.Last()) { // transcript starts at this node
			pathincov[i]+=abundance;
		}*/
		else { // middle of the path
			pathincov[i]+=abundance;
			pathoutcov[i]+=abundance;
		}
	}
}

/*
// I don't use this one
bool can_be_removed_back_all(float abundance,int lasttrnode,GVec<int>& path,GVec<float>& pathincov,GVec<float>& pathoutcov) { // this assumes all transcripts go through all the nodes

	if(pathoutcov.Last()-abundance<epsilon) return(false); // transcript leaving the last added node to path is the only one

	int lastnode=lasttrnode;
	if(path[0]<lasttrnode) // added transcript ends at path[0]
		lastnode=path[0];
	if(path.Last()!=lasttrnode) { // transcript doesn't end at last node added to path
		for(int i=path.Count()-2;i>=0;i--) {
			if(path[i]==lastnode) { // transcript ends at this node
				if(pathincov[i]-abundance<epsilon) return(false);
				break;
			}
			else // middle of the path
				if((pathincov[i]-abundance<epsilon) || (pathoutcov[i]-abundance<epsilon)) return(false);
		}
	}

  return(true);
}
*/

bool can_be_removed_back(float abundance,GBitVec& pattern,float& penalty,int lasttrnode,GVec<int>& path,GVec<float>& pathincov,GVec<float>& pathoutcov) {  // this assumes the incomplete transcripts only add to the nodes they pass through

	if(pattern[path.Last()] && pathoutcov.Last()-abundance<epsilon) { // transcript leaving the last added node to path is the only one
		penalty=abundance;
		return(false);
	}

	int lastnode=lasttrnode;
	if(path[0]<lasttrnode) // added transcript ends at path[0]
		lastnode=path[0];
	if(path.Last()!=lasttrnode) { // transcript doesn't end at last node added to path
		for(int i=path.Count()-2;i>=0;i--) {
			if(path[i]==lastnode) { // transcript ends at this node
				if(pattern[path[i]]) {
					penalty=abundance;
					if(pathincov[i]-abundance<epsilon) return(false);
				}
				break;
			}
			else // middle of the path
				if(pattern[path[i]]) {
					penalty=abundance;
					if((pathincov[i]-abundance<epsilon) || (pathoutcov[i]-abundance<epsilon)) return(false);
				}
		}
	}

  return(true);
}

bool can_be_removed_fwd(float abundance,GBitVec& pattern, float& penalty, int firstnode,GVec<int>& path,GVec<float>& pathincov,
		GVec<float>& pathoutcov) { // incomplete transcripts are all considered to go through their nodes

	//fprintf(stderr,"can be removed abundance=%f\n",abundance);
	if(pattern[path.Last()] && pathincov.Last()-abundance<epsilon) {
		penalty=abundance;
		//fprintf(stderr,"1: ret false: pathincov=%f abundance=%f\n",pathincov.Last(),abundance);
		return(false); // transcript entering the last added node to path is the only one
	}

	if(path.Last()!=firstnode) { // transcript doesn't start at last node added to path -> this shouldn't happen since the transcript has already been added
		for(int i=path.Count()-2;i>=0;i--) {
			//fprintf(stderr,"pathincov=%f pathoutcov=%f firstnode=%d path[%d]=%d\n",pathincov[i],pathoutcov[i],firstnode,i,path[i]);
			//if(pattern[path[i]]) fprintf(stderr,"pat!!\n");
			if(pattern[path[i]] && path[i]==firstnode) { // transcript starts at this node
				penalty=abundance;
				if(pathoutcov[i]-abundance<epsilon) {
					//fprintf(stderr,"2: ret false: pathoutcov[%d]=%f abundance=%f\n",i,pathoutcov[i],abundance);
					return(false);
				}
				break;
			}
			else // middle of the path
				if(pattern[path[i]]) {
					penalty=abundance;
					if((pathincov[i]-abundance<epsilon) || (pathoutcov[i]-abundance<epsilon)) {
						//fprintf(stderr,"3: ret false: pathincov[%d]=%f pathoutcov[%d]=%f abundance=%f\n",i,pathincov[i],i,pathoutcov[i],abundance);
						return(false);
					}
				}
		}
	}

  return(true);
}

/*
// I don't use this one
bool can_be_removed_fwd_all(float abundance,int firstnode,GVec<int>& path,GVec<float>& pathincov,
		GVec<float>& pathoutcov) {  // all transcripts are treated the same

	if(pathincov.Last()-abundance<epsilon) return(false); // transcript entering the last added node to path is the only one

	if(path.Last()!=firstnode) { // transcript doesn't start at last node added to path
		for(int i=path.Count()-2;i>=0;i--) {
			if(path[i]==firstnode) { // transcript starts at this node
				if(pathoutcov[i]-abundance<epsilon) return(false);
				break;
			}
			else // middle of the path
				if((pathincov[i]-abundance<epsilon) || (pathoutcov[i]-abundance<epsilon)) return(false);
		}
	}

  return(true);
}
*/

bool fwd_to_sink_fast(int i,GVec<int>& path,GBitVec& pathpat,GPVec<CTransfrag>& transfrag,GPVec<CGraphnode>& no2gnode,
		GVec<float>& nodecov,int gno,GIntHash<int>& gpos){

	// find all parents -> if parent is source then go back
	CGraphnode *inode=no2gnode[i];

	int nchildren=inode->child.Count(); // number of children

	if(!nchildren) return true; // node is sink
	int maxc=-1;

	/* I can not do this because then I don't check for continuity
	if(nchildren==1) {
		maxc=inode->child[0]; // only one child to consider
	}
	else {
	*/
	float maxcov=0;
	//int maxchild=inode->child[0];
	//float maxchildcov=-1;
	bool exclude=false;
	for(int c=0;c<nchildren;c++) {
		float childcov=0;
		CGraphnode *cnode=no2gnode[inode->child[c]];
		/*
			if(nodecov[inode->child[c]]>maxchildcov) {
				maxchildcov=nodecov[inode->child[c]];
				maxchild=inode->child[c];
			}
		*/
		if(sensitivitylevel && inode->child[c]==i+1 && i<gno-2 && inode->end+1==no2gnode[i+1]->start && cnode->len()
				&& nodecov[i+1]/cnode->len() <1000 && nodecov[i]*DROP>nodecov[i+1])  { // adjacent to child
			exclude=true;
		}
		else {
			pathpat[inode->child[c]]=1;
	    	int *pos=gpos[edge(i,inode->child[c],gno)];
	    	if(pos) pathpat[*pos]=1;
	    	else GError("Found parent-child: %d-%d not linked by edge!\n",i,inode->child[c]);
			for(int j=0;j<cnode->trf.Count();j++) { // for all transfrags going through child
				int t=cnode->trf[j];
				if(transfrag[t]->abundance<epsilon) { // this transfrag was used before -> needs to be deleted
					cnode->trf.Delete(j);
					j--;
				}
				else if(transfrag[t]->real && transfrag[t]->nodes[0]<=i && transfrag[t]->nodes.Last()>=inode->child[c] &&   // transfrag goes from i to c
						(transfrag[t]->pattern[inode->child[c]] || transfrag[t]->pattern[i]) &&    // transfrag is not incomplete through these nodes
						onpath(transfrag[t]->pattern,transfrag[t]->nodes,pathpat,path[0],inode->child[c],no2gnode,gno,gpos)) { // transfrag is compatible with path
					//fprintf(stderr,"add transfr[%d]->abund=%f\n",t,transfrag[t]->abundance);
					childcov+=transfrag[t]->abundance;
				}
			}

			if(childcov>maxcov) {
				maxcov=childcov;
				maxc=inode->child[c];
			}

			pathpat[inode->child[c]]=0;
			if(pos) pathpat[*pos]=0;
		}
	}
	if(maxc==-1) {
		if(exclude && nodecov[i+1]) {
			CGraphnode *cnode=no2gnode[i+1];
			float childcov=0;
			pathpat[i+1]=1;
	    	int *pos=gpos[edge(i,i+1,gno)];
	    	if(pos) pathpat[*pos]=1;
	    	else GError("Found parent-child: %d-%d not linked by edge\n",i,i+1);
			for(int j=0;j<cnode->trf.Count();j++) { // for all transfrags going through child
				int t=cnode->trf[j];
				if(transfrag[t]->abundance<epsilon) { // this transfrag was used before -> needs to be deleted
					cnode->trf.Delete(j);
					j--;
				}
				else if(transfrag[t]->real && transfrag[t]->nodes[0]<=i && transfrag[t]->nodes.Last()>=i+1 &&   // transfrag goes from i to c
						(transfrag[t]->pattern[i+1] || transfrag[t]->pattern[i]) &&    // transfrag is not incomplete through these nodes
						onpath(transfrag[t]->pattern,transfrag[t]->nodes,pathpat,path[0],i+1,no2gnode,gno,gpos)) { // transfrag is compatible with path
					childcov+=transfrag[t]->abundance;
				}
			}
			pathpat[i+1]=0; // these were 1 in the previous code, but it doesn't make sense -> I hope it was wrong
	    	if(pos) pathpat[*pos]=0;

			if(childcov) {
				maxc=i+1;
			}
			//else return false;
			else {
				//fprintf(stderr,"fwd: return false\n");
				return false;
			}
		}
		//else return false; //maxc=maxchild;
		else {
			//fprintf(stderr,"fwd: return false maxc=maxchild=%c\n",maxc);
			return false; //maxc=maxchild;
		}
	}

	/*
	{ // DEBUG ONLY
		printTime(stderr);
		fprintf(stderr,"maxc=%d\n",maxc);
	}
	*/

	// add maxp to path
	path.Add(maxc);
	pathpat[maxc]=1;

	int *pos=gpos[edge(i,maxc,gno)];
	if(pos) pathpat[*pos]=1;
	else GError("Found parent-child %d-%d not linked by edge\n",i,maxc);

	return fwd_to_sink_fast(maxc,path,pathpat,transfrag,no2gnode,nodecov,gno,gpos);
}

bool back_to_source_fast(int i,GVec<int>& path,GBitVec& pathpat,GPVec<CTransfrag>& transfrag,GPVec<CGraphnode>& no2gnode,
		GVec<float>& nodecov,int gno,GIntHash<int>& gpos){

	// find all parents -> if parent is source then go back
	CGraphnode *inode=no2gnode[i];

	int nparents=inode->parent.Count(); // number of parents

	if(!nparents) return true; // node is source
	int maxp=-1;

	/* I can not do this because then I don't check for continuity
	if(nparents==1) {
		maxp=inode->parent[0];
	}
	else {
	*/
	float maxcov=0;
	//int maxparent=inode->parent[0];
	//float maxparentcov=-1;
	bool exclude=false;
	for(int p=0;p<nparents;p++) {
		float parentcov=0;
		CGraphnode *pnode=no2gnode[inode->parent[p]];
		/*
			if(nodecov[inode->parent[p]]>maxparentcov) {
				maxparentcov=nodecov[inode->parent[p]];
				maxparent=inode->parent[p];
			}
		*/
		if(sensitivitylevel && inode->parent[p]==i-1 && i>1 && inode->start==no2gnode[i-1]->end+1 && pnode->len() &&
				nodecov[i-1]/pnode->len() <1000 && nodecov[i]*DROP>nodecov[i-1])  { // adjacent to parent
			exclude=true;
		}
		else {
			pathpat[inode->parent[p]]=1;
			int *pos=gpos[edge(inode->parent[p],i,gno)];
			if(pos) pathpat[*pos]=1;
			else GError("Found parent-child %d-%d not linked by edge\n",inode->parent[p],i);

			for(int j=0;j<pnode->trf.Count();j++) { // for all transfrags going through parent
				int t=pnode->trf[j];
				if(transfrag[t]->abundance<epsilon) { // this transfrag was used before -> needs to be deleted
					pnode->trf.Delete(j);
					j--;
				}
				else if(transfrag[t]->real && transfrag[t]->nodes[0]<=inode->parent[p] && transfrag[t]->nodes.Last()>=i &&   // transfrag goes from p to i
						(transfrag[t]->pattern[inode->parent[p]]||transfrag[t]->pattern[i]) &&    // transfrag is not incomplete through these nodes
						onpath(transfrag[t]->pattern,transfrag[t]->nodes,pathpat,inode->parent[p],path[0],no2gnode,gno,gpos)) { // transfrag is compatible with path
					// comment: I need to check for the parent to make sure I am not going through another node!!
					//fprintf(stderr,"parent=%d add transfr[%d]->abund=%f\n",inode->parent[p],t,transfrag[t]->abundance);
					parentcov+=transfrag[t]->abundance;
				}
			}

			if(parentcov>maxcov) {
				maxcov=parentcov;
				maxp=inode->parent[p];
			}

			pathpat[inode->parent[p]]=0;
			if(pos) pathpat[*pos]=0;
		}
	}
	if(maxp==-1) {
		if(exclude && nodecov[i-1]) {
			CGraphnode *pnode=no2gnode[i-1];
			float parentcov=0;
			pathpat[i-1]=1;
			int *pos=gpos[edge(i-1,i,gno)];
			if(pos) pathpat[*pos]=1;
			else GError("Found parent-child %d-%d not linked by edge\n",i-1,i);

			for(int j=0;j<pnode->trf.Count();j++) { // for all transfrags going through parent
				int t=pnode->trf[j];
				if(transfrag[t]->abundance<epsilon) { // this transfrag was used before -> needs to be deleted
					pnode->trf.Delete(j);
					j--;
				}
				else if(transfrag[t]->real && transfrag[t]->nodes[0]<=i-1 && transfrag[t]->nodes.Last()>=i &&   // transfrag goes from p to i
						(transfrag[t]->pattern[i-1]||transfrag[t]->pattern[i]) &&    // transfrag is not incomplete through these nodes
						onpath(transfrag[t]->pattern,transfrag[t]->nodes,pathpat,i-1,path[0],no2gnode,gno,gpos)) { // transfrag is compatible with path
					// comment: I need to check for the parent to make sure I am not going through another node!!
					//fprintf(stderr,"parent=%d add transfr[%d]->abund=%f\n",i-1,t,transfrag[t]->abundance);
					parentcov+=transfrag[t]->abundance;
				}
			}
			pathpat[i-1]=0;
			if(pos) pathpat[*pos]=0;

			if(parentcov) {
				maxp=i-1;
			}
			//else return false;
			else {
				//fprintf(stderr,"return false\n");
				return false;
			}
		}
		//else return false; //maxp=maxparent;
		else {
			//fprintf(stderr,"return false maxp=maxparent=%d\n",maxp);
			return false; //maxp=maxparent;
		}
	}

	/*
	{ // DEBUG ONLY
		printTime(stderr);
		fprintf(stderr,"maxp=%d\n",maxp);
	}
	*/

	if(maxp) { // add maxp to path only if not source
		path.Add(maxp);
	}

	pathpat[maxp]=1;                 // if maxp is source I added it in the pathpat
	int *pos=gpos[edge(maxp,i,gno)];
	if(pos) pathpat[*pos]=1;
	else GError("Found parent-child %d-%d not linked by edge\n",maxp,i);

	return back_to_source_fast(maxp,path,pathpat,transfrag,no2gnode,nodecov,gno,gpos);
}

bool back_to_source_path(int i,GVec<int>& path,GBitVec& pathpat,GVec<float>& pathincov, GVec<float>& pathoutcov,
		GBitVec& istranscript,GBitVec& removable,GPVec<CTransfrag>& transfrag,GHash<CComponent>& computed,
		GBitVec& compatible,GPVec<CGraphnode>& no2gnode,GVec<float>& nodecov,int gno,GIntHash<int>& gpos){

	// find all parents -> if parent is source then go back
	CGraphnode *inode=no2gnode[i];

	if(!inode->parent[0]) { // parent is source; what if it has other parent than source?
		pathpat[0]=1;
		int *pos=gpos[edge(0,i,gno)];
		if(pos) pathpat[*pos]=1;
		else GError("Found parent-child 0-%d not linked by edge\n",i);

		return true; // parent doesn't get added to the path, but it gets added to the pattern
	}

	GVec<CTrInfo> set;

	/*
	{ // DEBUG ONLY
		printTime(stderr);
		fprintf(stderr,"start back_to_source at %d\n",i);
	}
	*/

	// I should penalize somehow the fact that an adjacent parent has to big a drop in coverage;
	bool exclude=false;
	if(sensitivitylevel && i>1 && inode->start==1+no2gnode[i-1]->end && nodecov[i]*DROP>nodecov[i-1])  { // adjacent to parent
		//fprintf(stderr,"node %d to adjacent parent %d has a drop in cov from %f to %f\n",i,i-1,nodecov[i],nodecov[i-1]);
		exclude=true;
	}


	//collect all in transfrags that are compatible with path
	for(int j=0;j<inode->trf.Count();j++) {
		int t=inode->trf[j];
		if(transfrag[t]->abundance<epsilon) { // this transfrag was used before -> needs to be deleted
			inode->trf.Delete(j);
			j--;
		}
		else if(transfrag[t]->nodes[0]!=i && transfrag[t]->real) { // this is an in transfrag; I might change this to include the cases where I have a positive in for the transfrag, e.g. if(transfrag[inode->trf[j]]->nodes[0]!=i || transfrag[inode->trf[j]]->in)
			if(istranscript[t] || onpath(transfrag[t]->pattern,transfrag[t]->nodes,pathpat,i,path[0],no2gnode,gno,gpos)) {
				float abund=transfrag[t]->abundance;
				float penalty=0;
				if(!transfrag[t]->pattern[i]) abund=0; // this is an incomplete transfrag -> I might reconsider this case
				if(istranscript[t]) {
					//penalty=transfrag[t]->abundance;
					if(can_be_removed_back(transfrag[t]->abundance,transfrag[t]->pattern,penalty,transfrag[t]->nodes.Last(),path,pathincov,pathoutcov)) {
					//if(can_be_removed_back_all(transfrag[t]->abundance,transfrag[t]->nodes.Last(),path,pathincov,pathoutcov))
						removable[t]=1;
						if(exclude && transfrag[t]->pattern[i-1]) {
							update_transcript_to_path_back(-transfrag[t]->abundance,transfrag[t]->nodes,path,pathincov,pathoutcov);
							istranscript[t]=false;
						}
					}
					else removable[t]=0;
				}
				if(istranscript[t] || !exclude || !transfrag[t]->pattern[i-1]) {
					CTrInfo tr(t,abund,penalty);
					set.Add(tr);
				}
			}
		}
	}

	//  here I could collect all out transfrags -> this are all different transfrags since they are ending at different parents

	/*
	{ // DEBUG ONLY
		printTime(stderr);
		fprintf(stderr,"Compute max component out of %d transfrags\n",set.Count());
	}
	*/

	int maxp=-1; // maximum parent
	// compute the maximum component
	GVec<int> *maxset=NULL;
	if(set.Count()) {
		if(set.Count()<MAX_MAXCOMP) { // I can compute the max component in reasonable time
			// set.Sort(setCmp); no need to sort here -> only (maybe) if I called the out transfrag without binary insert
			float maxsize=MIN_VAL;
			maxset=max_compon_size_with_penalty(transfrag.Count(),maxsize,set,compatible,istranscript,removable,computed);
			if(!maxset) {
				if(i>1 && exclude && inode->parentpat[i-1]) maxp=i-1;
				else return(false);
			}
		}
		else { // I need to go fast
			if(back_to_source_fast(i,path,pathpat,transfrag,no2gnode,nodecov,gno,gpos)) {
				int n=path.Count();
				istranscript.reset();
				for(int z=0;z<n;z++) {
					if(path[z]<i) {
						pathincov.cAdd(0.0);
						pathoutcov.cAdd(0.0);
					}
					else {
						pathincov[z]=0;
						pathoutcov[z]=0;
					}
					int nt=no2gnode[path[z]]->trf.Count();
					for(int j=0;j<nt;j++) {
						int t=no2gnode[path[z]]->trf[j];
						if(transfrag[t]->real && transfrag[t]->abundance && onpath(transfrag[t]->pattern,transfrag[t]->nodes,pathpat,path.Last(),path[0],no2gnode,gno,gpos)) {

							istranscript[t]=1;
							if(transfrag[t]->nodes[0]<path[z]) { // transfrag enters this node
								pathincov[z]+=transfrag[t]->abundance;
							}
							if(transfrag[t]->nodes.Last()>path[z]) { // transfrag exits this node
								pathoutcov[z]+=transfrag[t]->abundance;
							}
						}
					}
				}
				// now deal with source node
				pathincov.cAdd(0.0);
				pathoutcov.cAdd(0.0);
				int nt=no2gnode[0]->trf.Count();
				for(int j=0;j<nt;j++) {
					int t=no2gnode[0]->trf[j];
					if(transfrag[t]->abundance && onpath(transfrag[t]->pattern,transfrag[t]->nodes,pathpat,0,path[0],no2gnode,gno,gpos)) {
						istranscript[t]=1;
						pathoutcov.Last()+=transfrag[t]->abundance;
					}
				}
				return(true);
			}
			else return(false);
		}
	}
	else if(i>1 && exclude && inode->parentpat[i-1]) maxp=i-1;
	else return(false);

	/*
	{ // DEBUG ONLY
		printTime(stderr);
		fprintf(stderr,"Max set[%d] bck: ",i);
		if(maxset) for(int t=0;t<maxset->Count();t++) fprintf(stderr," %d",maxset->Get(t));
		fprintf(stderr," from set: ");
		for(int t=0;t<set.Count();t++) fprintf(stderr," %d",set[t].trno);
		fprintf(stderr,"\n");
	}
	*/


	int possiblep=-1; // possible maximum parent
	if(maxset) {
		int j=0;
		int k=0;
		while(k<set.Count()) {
			while((j<maxset->Count() && set[k].trno<maxset->Get(j)) || (j==maxset->Count() && k<set.Count())) {
				if(istranscript[set[k].trno]) { // this transcript got removed
					update_transcript_to_path_back(-transfrag[set[k].trno]->abundance,transfrag[set[k].trno]->nodes,path,pathincov,pathoutcov);
					istranscript[set[k].trno]=0;
				}
				k++;
			}

			if(j<maxset->Count()) {
				int t=maxset->Get(j);
				if(transfrag[t]->nodes[0]) { // this transcript doesn't start from source
					if(!istranscript[t]) {
						istranscript[t]=1;
						update_transcript_to_path_back(transfrag[t]->abundance,transfrag[t]->nodes,path,pathincov,pathoutcov);
					}
					pathincov.Last()+=transfrag[t]->abundance; // I need to update the pathincov of node i as well

					if(maxp==-1 && transfrag[t]->nodes[0]<i) { // attempt to determine maxp
						bool found=false;
						for(int p=0;p<inode->parent.Count();p++) {
							int *pos=gpos[edge(inode->parent[p],i,gno)];
							if(pos && transfrag[t]->pattern[*pos]) {
								maxp=inode->parent[p];
								found=true;
								break;
							}
						}
						if(!found) { // find maximum possible parent
							int p=0;
							while(p<transfrag[t]->nodes.Count() && transfrag[t]->nodes[p]<i) p++;
							if(p && transfrag[t]->nodes[p-1]>possiblep) possiblep=transfrag[t]->nodes[p-1];
						}
					}
				}
				else maxp=0;
				j++;
			} // end j<maxset->Count()
			k++;
		} // end while(k<set.Count())
	}
	else if(set.Count()) {
		for(int k=0;k<set.Count();k++)
			if(istranscript[set[k].trno]) { // this transcript got removed
				update_transcript_to_path_back(-transfrag[set[k].trno]->abundance,transfrag[set[k].trno]->nodes,path,pathincov,pathoutcov);
				istranscript[set[k].trno]=0;
			}
	}

	if(maxset) delete maxset;

	if(maxp==-1) {
		if(possiblep==-1) { // all transcripts are starting transcripts -> I need to find the node with maximum coverage to continue path
			float maxcov=nodecov[inode->parent[0]];
			maxp=inode->parent[0];
			for(int p=1;p<inode->parent.Count();p++) {
				if(nodecov[inode->parent[p]]>maxcov) {
					maxp=inode->parent[p];
					maxcov=nodecov[maxp];
				}
			}
		}
		else {
			float maxcov=-1;
			for(int p=0;p<inode->parent.Count();p++) {
				if(inode->parent[p]==possiblep) {
					maxp=possiblep;
					break;
				}
				else if(no2gnode[inode->parent[p]]->parentpat[possiblep]) {
					if(nodecov[inode->parent[p]]>maxcov) {
						maxp=inode->parent[p];
						maxcov=nodecov[maxp];
					}
				}
			}
		}
	}

	/*
	{ // DEBUG ONLY
		printTime(stderr);
		fprintf(stderr,"maxp=%d\n",maxp);
	}
	*/

	if(maxp==-1) return false;

	// need to update pathincov of i!!!
	pathincov.cAdd(0.0);
	pathoutcov.cAdd(0.0);
	for(int j=0;j<no2gnode[maxp]->trf.Count();j++) {
		int t=no2gnode[maxp]->trf[j];
		if(istranscript[t]) pathoutcov.Last()+=transfrag[t]->abundance;
	}

	pathpat[maxp]=1;
	int *pos=gpos[edge(maxp,i,gno)];
	if(pos) pathpat[*pos]=1;
	else GError("Found parent-child %d-%d not linked by edge\n",maxp,i);

	if(!maxp) return(true);

	// add maxp to path
	path.Add(maxp);

	return back_to_source_path(maxp,path,pathpat,pathincov,pathoutcov,istranscript,removable,transfrag,computed,compatible,no2gnode,nodecov,gno,gpos);
}

bool fwd_to_sink_path(int i,GVec<int>& path,GBitVec& pathpat,GVec<float>& pathincov, GVec<float>& pathoutcov,
		GBitVec& istranscript,GBitVec& removable,GPVec<CTransfrag>& transfrag,GHash<CComponent>& computed,
		GBitVec& compatible,GPVec<CGraphnode>& no2gnode,GVec<float>& nodecov,int gno,GIntHash<int>& gpos){

	// find all parents -> if parent is source then go back
	CGraphnode *inode=no2gnode[i];

	if(!inode->child.Count()) return true; // parent is sink (no children)

	GVec<CTrInfo> set;

	// I should penalize somehow the fact that an adjacent child has to big a drop in coverage;
	bool exclude=false;
	if(sensitivitylevel && i<gno-2 && inode->end+1==no2gnode[i+1]->start && nodecov[i]*DROP>nodecov[i+1])  { // adjacent to parent
		//fprintf(stderr,"node %d to adjacent child %d has a drop in cov from %f to %f\n",i,i+1,nodecov[i],nodecov[i+1]);
		exclude=true;
	}

	//collect all in transfrags that are compatible with path
	for(int j=0;j<inode->trf.Count();j++) {
		int t=inode->trf[j];
		if(transfrag[t]->abundance<epsilon) { // this transfrag was used before -> needs to be deleted
			inode->trf.Delete(j);
			j--;
		}
		else if(transfrag[t]->nodes.Last()!=i && transfrag[t]->real) { // this is an out transfrag; I might change this to include the cases where I have a positive out for the transfrag, e.g. if(transfrag[inode->trf[j]]->nodes.Last()!=i || transfrag[inode->trf[j]]->out)
			if(istranscript[t] || onpath(transfrag[t]->pattern,transfrag[t]->nodes,pathpat,path[0],i,no2gnode,gno,gpos)) {
				float abund=transfrag[t]->abundance;
				float penalty=0;
				if(!transfrag[t]->pattern[i]) abund=0; // this is an incomplete transfrag -> I might reconsider this case
				if(istranscript[t]) {
					//penalty=transfrag[t]->abundance;
					if(can_be_removed_fwd(transfrag[t]->abundance,transfrag[t]->pattern,penalty,transfrag[t]->nodes[0],path,pathincov,pathoutcov)) {
					//if(can_be_removed_fwd_all(transfrag[t]->abundance,transfrag[t]->nodes[0],path,pathincov,pathoutcov))
						removable[t]=1;
						if(exclude && transfrag[t]->pattern[i+1]) {
							update_transcript_to_path_fwd(-transfrag[t]->abundance,transfrag[t]->nodes,path,pathincov,pathoutcov);
							istranscript[t]=false;
						}
					}
					else removable[t]=0;
				}
				/*
				if(t==90 || t==69) {
					fprintf(stderr,"t=%d abund=%f",t,transfrag[t]->abundance);
					if(istranscript[t]) fprintf(stderr,"istranscript ");
					if(removable[t]) fprintf(stderr,"removable ");
					if(exclude) fprintf(stderr,"exclude ");
					if(transfrag[t]->pattern[i+1]) fprintf(stderr,"has %d child ",i+1);
					fprintf(stderr,"\n");
				}*/
				if(istranscript[t] || !exclude || !transfrag[t]->pattern[i+1]) {
					CTrInfo tr(t,abund,penalty);
					set.Add(tr);
				}
			}
		}
	}

	//  here I could collect all out transfrags -> this are all different transfrags since they are ending at different parents
	int maxc=0; // maximum child
	// compute the maximum component
	GVec<int> *maxset=NULL;
	if(set.Count()) {
		if(set.Count()<MAX_MAXCOMP) {
			// set.Sort(setCmp); no need to sort here -> only (maybe) if I collected the out transfrag without binary insert
			float maxsize=MIN_VAL;
			maxset=max_compon_size_with_penalty(transfrag.Count(),maxsize,set,compatible,istranscript,removable,computed);
			if(!maxset) {
				if(i<gno-2 && exclude && inode->childpat[i+1]) maxc=i+1;
				else return(false);
			}
		}
		else {
			if(fwd_to_sink_fast(i,path,pathpat,transfrag,no2gnode,nodecov,gno,gpos)) {
				istranscript.reset();
				return(true);
			}
			else return(false);

		}
	}
	else if(i<gno-2 && exclude && inode->childpat[i+1]) maxc=i+1;
	else {
		return(false);
	}

	/*
	{ // DEBUG ONLY
		printTime(stderr);
		fprintf(stderr,"Max set[%d] fwd: ",i);
		if(maxset) for(int t=0;t<maxset->Count();t++) fprintf(stderr," %d",maxset->Get(t));
		fprintf(stderr," from set: ");
		for(int t=0;t<set.Count();t++) fprintf(stderr," %d",set[t].trno);
		fprintf(stderr,"\n");
	}
	*/


	int possiblec=gno; // possible maximum parent
	if(maxset) {
		int j=0;
		int k=0;
		while(k<set.Count()) {
			while((j<maxset->Count() && set[k].trno<maxset->Get(j)) || (j==maxset->Count() && k<set.Count())) {
				if(istranscript[set[k].trno]) { // this transcript got removed
					update_transcript_to_path_fwd(-transfrag[set[k].trno]->abundance,transfrag[set[k].trno]->nodes,path,pathincov,pathoutcov);
					istranscript[set[k].trno]=0;
				}
				k++;
			}

			if(j<maxset->Count()) {
				int t=maxset->Get(j);
				if(!istranscript[t]) {
					istranscript[t]=1;
					update_transcript_to_path_fwd(transfrag[t]->abundance,transfrag[t]->nodes,path,pathincov,pathoutcov);
				}
				pathoutcov.Last()+=transfrag[t]->abundance; // I need to update the pathoutcov of node i as well
				if(!maxc && transfrag[t]->nodes.Last()>i) { // attempt to determine maxc = maximum child
					bool found=false;
					for(int c=0;c<inode->child.Count();c++) {
						int *pos=gpos[edge(i,inode->child[c],gno)];
						if(pos && transfrag[t]->pattern[*pos]) {
							maxc=inode->child[c];
							found=true;
							break;
						}
					}
					if(!found) { // find mimimum possible child
						int c=0;
						while(c<transfrag[t]->nodes.Count() && transfrag[t]->nodes[c]<=i) c++;
						if(c<transfrag[t]->nodes.Count() && transfrag[t]->nodes[c]<possiblec) possiblec=transfrag[t]->nodes[c];
					}
				}
				j++;
			} // end j<maxset->Count()
			k++;
		} // end while(k<set.Count())
	}
	else if(set.Count()) {
		for(int k=0;k<set.Count();k++)
			if(istranscript[set[k].trno]) { // this transcript got removed
				update_transcript_to_path_fwd(-transfrag[set[k].trno]->abundance,transfrag[set[k].trno]->nodes,path,pathincov,pathoutcov);
				istranscript[set[k].trno]=0;
			}
	}

	if(maxset) delete maxset;

	if(!maxc) {
		if(possiblec==gno) {
			float maxcov=nodecov[inode->child[0]];
			maxc=inode->child[0];
			for(int c=1;c<inode->child.Count();c++) {
				if(nodecov[inode->child[c]]>maxcov) {
					maxc=inode->child[c];
					maxcov=nodecov[maxc];
				}
			}
		}
		else {
			float maxcov=-1;
			for(int c=0;c<inode->child.Count();c++) {
				if(inode->child[c]==possiblec) {
					maxc=possiblec;
					break;
				}
				else if(no2gnode[inode->child[c]]->childpat[possiblec]) {
					if(nodecov[inode->child[c]]>maxcov) {
						maxc=inode->child[c];
						maxcov=nodecov[maxc];
					}
				}
			}
			maxc=possiblec;
		}
	}

	/*
	{ // DEBUG ONLY
		printTime(stderr);
		fprintf(stderr,"maxc=%d\n",maxc);
	}
	*/

	if(!maxc) return false;

	// add maxc to path
	path.Add(maxc);
	pathincov.cAdd(0.0);
	pathoutcov.cAdd(0.0);
	for(int j=0;j<no2gnode[maxc]->trf.Count();j++) {
		int t=no2gnode[maxc]->trf[j];
		if(istranscript[t]) pathincov.Last()+=transfrag[t]->abundance;
	}
	pathpat[maxc]=1;
	int *pos=gpos[edge(i,maxc,gno)];
	if(pos) pathpat[*pos]=1;
	else GError("Found parent-child %d-%d not linked by edge\n",i,maxc);

	return fwd_to_sink_path(maxc,path,pathpat,pathincov,pathoutcov,istranscript,removable,transfrag,computed,compatible,no2gnode,nodecov,gno,gpos);
}

void eliminate_tr_left(int t,float val,GVec<float>& used,GVec<int>& usedpos,GVec<float>& capacity,	GVec<int>& path,int i,
		GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,GBitVec& istranscript) {

	usedpos[i]=0;

	if(!transfrag[t]->pattern[path[i]]) // node is not among transfrag's nodes
		eliminate_tr_left(t,val,used,usedpos,capacity,path,i-1,no2gnode,transfrag,istranscript);
	else {
		capacity[i]+=val;
		if(transfrag[t]->nodes[0]<path[i]) {
			eliminate_tr_left(t,val,used,usedpos,capacity,path,i-1,no2gnode,transfrag,istranscript);
		}
		else {
			CGraphnode *inode=no2gnode[path[i]];
			if(i && inode->rate){
				float newval=0;
				int nt=inode->trf.Count();
				int j=0;
				while(j<nt && newval<val) {
					int t=inode->trf[j];
					if(istranscript[t]&&transfrag[t]->nodes.Last()==path[i]) { // transcript ends at this point
						float maxval=used[t];
						float allowed=(val-newval)/inode->rate;
						if(allowed<epsilon) break;
						if(maxval>allowed) maxval=allowed;
						if(maxval<epsilon) maxval=0;
						else {
							used[t]-=maxval;
							if(used[t]<epsilon) used[t]=0;
							transfrag[t]->abundance-=maxval;
							if(transfrag[t]->abundance<epsilon) {
								transfrag[t]->abundance=0;
								inode->trf.Delete(j);
								j--;
								nt--;
							}
							eliminate_tr_left(t,maxval,used,usedpos,capacity,path,i-1,no2gnode,transfrag,istranscript);
							newval+=maxval*inode->rate;
						}
					}
					j++;
				}
			}
		}
	}
}

void eliminate_tr_right(int t,float val,GVec<float>& used,GVec<int>& usedpos,GVec<float>& capacity,GVec<int>& path,int i,
		GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,GBitVec& istranscript) {

	usedpos[i]=0;
	if(!transfrag[t]->pattern[path[i]]) // node is not among transfrag's nodes
		eliminate_tr_right(t,val,used,usedpos,capacity,path,i+1,no2gnode,transfrag,istranscript);
	else {
		if(transfrag[t]->nodes.Last()>path[i]) {
			capacity[i]+=val;
			eliminate_tr_right(t,val,used,usedpos,capacity,path,i+1,no2gnode,transfrag,istranscript);
		}
		else {
			CGraphnode *inode=no2gnode[path[i]];
			if(inode->rate) capacity[i]+=val*inode->rate;
			else capacity[i]+=val;

			if(i<path.Count()-1 && inode->rate){
				CGraphnode *inode=no2gnode[path[i]];
				float newval=0;
				int nt=inode->trf.Count();
				int j=0;
				while(j<nt && newval<val) {
					int t=inode->trf[j];
					if(istranscript[t]&&transfrag[t]->nodes[0]==path[i]) { // transcript starts at this point
						float maxval=used[t];
						float allowed=(val-newval)*inode->rate;
						if(allowed<epsilon) break;
						if(maxval>allowed) maxval=allowed;
						if(maxval<epsilon) maxval=0;
						else {
							used[t]-=maxval;
							if(used[t]<epsilon) used[t]=0;
							transfrag[t]->abundance-=maxval;
							if(transfrag[t]->abundance<epsilon) {
								transfrag[t]->abundance=0;
								inode->trf.Delete(j);
								j--;
								nt--;
							}
							eliminate_tr_right(t,maxval,used,usedpos,capacity,path,i+1,no2gnode,transfrag,istranscript);
							newval+=maxval/inode->rate;
						}
					}
					j++;
				}
			}
		}
	}
}

float extend_tr_left(int t,float val,GVec<float>& used,GVec<int>& usedpos,int& update,GVec<int>& path,int i,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,
		GVec<int>& node2path,GBitVec& istranscript) {

	if(!transfrag[t]->nodes[0]) return(val);
	float newval=0;

	if(transfrag[t]->nodes[0]<path[i])
		newval=extend_tr_left(t,val,used,usedpos,update,path,node2path[transfrag[t]->nodes[0]],no2gnode,transfrag,node2path,istranscript);
	else { // transfrag starts at i
		CGraphnode *inode=no2gnode[path[i]];
		if(!inode->rate || !i) return(val); // transcript ends here
		int nt=inode->trf.Count();
		int j=usedpos[i];
		while(j<nt) {
			int t=inode->trf[j];
			usedpos[i]=j;
			if(j && i<update) update=i;
			if(istranscript[t] && transfrag[t]->nodes.Last()==path[i]) { // transcript ends at this point
				float maxval=transfrag[t]->abundance-used[t];
				if(maxval>0) {
					float allowed=(val-newval)/inode->rate;
					if(allowed<epsilon) return(newval);
					if(maxval>allowed) maxval=allowed;
					maxval=extend_tr_left(t,maxval,used,usedpos,update,path,node2path[transfrag[t]->nodes[0]],no2gnode,transfrag,node2path,istranscript);
					newval+=maxval*inode->rate;
					used[t]+=maxval;
					if(newval+epsilon>=val) return(val);
				}
			}
			j++;
		}
	}
	return(newval);
}

float extend_tr_right(int t,float val,GVec<float>& used,GVec<int>& usedpos,int& update,GVec<int>& path,int i,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,
		GVec<int>& node2path,GBitVec& istranscript) {

	float newval=0;

	if(transfrag[t]->nodes.Last()>path[i])
		newval=extend_tr_right(t,val,used,usedpos,update,path,node2path[transfrag[t]->nodes.Last()],no2gnode,transfrag,node2path,istranscript);
	else { // transfrag starts at i
		CGraphnode *inode=no2gnode[path[i]];
		if(!inode->rate || i==path.Count()-1) return(val); // transcript ends here
		int nt=inode->trf.Count();
		int j=usedpos[i];
		while(j<nt) {
			int t=inode->trf[j];
			usedpos[i]=j;
			if(j && i>update) update=i;
			if(istranscript[t] && transfrag[t]->nodes[0]==path[i]) { // transcript starts at this point
				float maxval=transfrag[t]->abundance-used[t];
				if(maxval>0) {
					float allowed=(val-newval)*inode->rate;
					if(allowed<epsilon) return(newval);
					if(maxval>=allowed) maxval=allowed;
					maxval=extend_tr_right(t,maxval,used,usedpos,update,path,node2path[transfrag[t]->nodes.Last()],no2gnode,transfrag,node2path,istranscript);
					newval+=maxval/inode->rate;
					used[t]+=maxval;
					if(newval+epsilon>=val) return(val);
				}
			}
			j++;
		}
	}
	return(newval);
}

void update_capacity(int start,CTransfrag *t,float val,GVec<float>& capacity,GVec<int>& node2path) {
	t->abundance-=val;
	if(t->abundance<epsilon) t->abundance=0;
	for(int j=start;j<t->nodes.Count()-1;j++) {
		int i=node2path[t->nodes[j]];
		capacity[i]+=val;
	}
}

void compute_capacity(int lastn, CTransfrag *t,float val,GVec<float>& capacity,GVec<int>& node2path,GPVec<CGraphnode>& no2gnode) {
	int j=0;
	while(j<t->nodes.Count()-1 && t->nodes[j]<lastn) {
		int i=node2path[t->nodes[j]];
		capacity[i]+=val;
		j++;
	}
	if(j<t->nodes.Count() && t->nodes[j]==lastn) {
		if(j<t->nodes.Count()-1)
			capacity[node2path[t->nodes[j]]]+=val;
		else capacity[node2path[t->nodes[j]]]+=val*no2gnode[t->nodes[j]]->rate;
	}
}


void compute_capacity_back(int firstn, CTransfrag *t,float val,GVec<float>& capacity,GVec<int>& node2path) {
	int j=t->nodes.Count()-2;
	while(j>=0 && t->nodes[j]>=firstn) {
		int i=node2path[t->nodes[j]];
		capacity[i]+=val;
		j--;
	}
}


bool bfs(int n,GVec<float> *capacity,GVec<float> *flow,GVec<int> *link,GVec<int>& pred) {
	GVec<int> color;
	color.Resize(n+2,0);
	int head=0;
	int tail=0;
	GVec<int> q;

	// enque 0 (source)
	q.cAdd(0);
	tail++;
	color[0]=1;
	pred[0]=-1;
	while(head!=tail) {
		// deque
		int u=q[head];
		head++;
		color[u]=2;
		//for(int v=0;v<link[u].Count();v++)
		for(int v=link[u].Count()-1;v>=0;v--) // do this so longer transcripts would have higher priority (be considered first)
			if(!color[link[u][v]] && capacity[u][link[u][v]]-flow[u][link[u][v]]>epsilon) {
				// enque v
				q.Add(link[u][v]);
				tail++;
				color[link[u][v]]=1;
				pred[link[u][v]]=u;
			}
	}

	return(color[n-1]==2);
}


bool weight_bfs(int n,GVec<float> *capacity,GVec<float> *flow,GVec<int> *link,GVec<int>& pred) {
	GVec<int> color;
	color.Resize(n,0);
	int head=0;
	int tail=0;
	GVec<int> q;

	// enque 0 (source)
	q.cAdd(0);
	tail++;
	color[0]=1;
	pred[0]=-1;
	while(head!=tail) {
		// deque
		int u=q[head];
		head++;
		color[u]=2;
		for(int v=0;v<link[u].Count();v++)
			if(!color[link[u][v]] && capacity[u][link[u][v]]-flow[u][link[u][v]]>epsilon && (link[u][v]<u || capacity[u][u]-flow[u][u]>epsilon)) {
				// enque v
				q.Add(link[u][v]);
				tail++;
				color[link[u][v]]=1;
				pred[link[u][v]]=u;
			}
	}

	return(color[n-1]==2);
}

/*
// I don't use this one : it doesn't work!
bool weight_bfs_EM(int n,GVec<float> *capacity,GVec<float> *flow,GVec<int> *link,GVec<int>& pred) {
	GVec<int> color;
	color.Resize(n+2,0);
	int head=0;
	int tail=0;
	GVec<int> q;

	// enque 0 (source)
	q.cAdd(0);
	tail++;
	color[0]=1;
	pred[0]=-1;
	while(head!=tail) {
		// deque
		int u=q[head];
		head++;
		color[u]=2;
		for(int v=0;v<link[u].Count();v++)
			if(!color[link[u][v]] && capacity[u][link[u][v]]-flow[u][link[u][v]]>epsilon && (link[u][v]<u || capacity[u][u]-flow[u][u]>epsilon || link[u][v]==n || (u==n-1 && link[u][v]==n+1))) {
				// enque v
				q.Add(link[u][v]);
				tail++;
				color[link[u][v]]=1;
				pred[link[u][v]]=u;
			}
	}

	return(color[n-1]==2);
}
*/

float find_capacity(int i,GVec<int>& path,GVec<float>& nodeflux,GPVec<CGraphnode>& no2gnode) {
	int j=i-1;
	while(j>=0 && !no2gnode[path[j]]->capacity) j--;
	if(j>=0) return(nodeflux[j]/no2gnode[path[j]]->capacity);
	j=i+1;
	while(j<path.Count() && !no2gnode[path[j]]->capacity) j++;
	if(j<path.Count()) return(nodeflux[j]/no2gnode[path[j]]->capacity);
	return(0.0);
}

float find_capacity_back(int i,GVec<int>& path,GVec<float>& nodeflux,GPVec<CGraphnode>& no2gnode) {
	int n=path.Count();
	int j=i-1;
	while(j>=0 && !no2gnode[path[n-1-j]]->capacity) j--;
	if(j>=0) return(nodeflux[j]/no2gnode[path[j]]->capacity);
	j=i+1;
	while(j<n && !no2gnode[path[n-1-j]]->capacity) j++;
	if(j<n) return(nodeflux[j]/no2gnode[path[n-1-j]]->capacity);
	return(0.0);
}

/*
// I don't use this one
float max_flow_partial(int lastn,int gno,GVec<int>& path,GBitVec& istranscript,GPVec<CTransfrag>& transfrag,GPVec<CGraphnode>& no2gnode,
		GVec<int>& node2path,GVec<float>& nodecapacity,GVec<int>& nodes,bool& all) {

	float flux=0;
	int n=path.Count();
	GVec<float> *capacity=new GVec<float>[n]; // capacity of edges in network
	GVec<float> *flow=new GVec<float>[n]; // flow in network
	GVec<int> *link=new GVec<int>[n]; // for each node remembers it's neighbours
	GVec<int> pred; // this stores the augmenting path
	pred.Resize(n,-1);
	nodecapacity.Clear();

	nodes.Clear();
	all=true;

	for(int i=0;i<n;i++) {
		node2path[path[i]]=i;
		nodecapacity.cAdd(0.0);
		capacity[i].Resize(n,0);
		flow[i].Resize(n,0);

		nodes.cAdd(0.0);

	}


	// establish capacities in the network
	for(int i=0;i<n;i++) {
		int nt=no2gnode[path[i]]->trf.Count();
		for(int j=0;j<nt;j++) {
			int t=no2gnode[path[i]]->trf[j];
			if(istranscript[t]) {
				if(transfrag[t]->nodes[0]==path[i]) { // transfrag starts at this node
					int n1=i;
					int n2=node2path[transfrag[t]->nodes.Last()];
					if(transfrag[t]->nodes.Last()>lastn) n2=n-1;
					if(!no2gnode[path[i]]->rate) n1=0;
					if(!no2gnode[transfrag[t]->nodes.Last()]->rate) n2=n-1;
					if(!capacity[n1][n2]) {
						link[n1].Add(n2);
						link[n2].Add(n1);
					}
					capacity[n1][n2]+=transfrag[t]->abundance;
				}
			}
		}
	}

	GVec<float> rate;
	rate.Resize(n,1);

	while(bfs(n,capacity,flow,link,pred)) {
		int r=0;
		float increment=FLT_MAX;
		rate[r++]=1;
		for(int u=n-1;pred[u]>=0;u=pred[u]) {
			float adjflux=(capacity[pred[u]][u]-flow[pred[u]][u])*rate[r-1];
			increment = increment < adjflux ? increment : adjflux;
			if(pred[pred[u]]>=0) {
				if(pred[u]<u) {
					if(pred[pred[u]]<pred[u]) rate[r]=rate[r-1]*no2gnode[path[pred[u]]]->rate;
					else rate[r]=rate[r-1];
				}
				else {
					if(pred[pred[u]]<pred[u]) rate[r]=rate[r-1];
					else rate[r]=rate[r-1]/no2gnode[path[pred[u]]]->rate;
				}
				r++;
			}
		}
		r=0;
		for(int u=n-1;pred[u]>=0;u=pred[u]) {
			flow[pred[u]][u]+=increment/rate[r];
			flow[u][pred[u]]-=increment/rate[r];
			r++;
		}
		flux+=increment;
	}

	//fprintf(stderr,"flux=%f ",flux);

	for(int i=0;i<n;i++) {
		int nt=no2gnode[path[i]]->trf.Count();
		for(int j=0;j<nt;j++) {
			int t=no2gnode[path[i]]->trf[j];
			if(istranscript[t] && transfrag[t]->abundance) {
				if(transfrag[t]->nodes[0]==path[i]) { // transfrag starts at this node
					int n1=i;
					int n2=node2path[transfrag[t]->nodes.Last()];
					if(transfrag[t]->nodes.Last()>lastn) n2=n-1;
					if(!no2gnode[path[i]]->rate) n1=0;
					if(!no2gnode[transfrag[t]->nodes.Last()]->rate) n2=n-1;
					if(flow[n1][n2]>0) {

						if(abs(flow[n1][n2]-capacity[n1][n2])<epsilon) { nodes[n1]=1;if(n1) all=false;}

						if(flow[n1][n2]<transfrag[t]->abundance) {
							compute_capacity(lastn,transfrag[t],flow[n1][n2],nodecapacity,node2path,no2gnode);
							flow[n1][n2]=0;
						}
						else {
							flow[n1][n2]-=transfrag[t]->abundance;
							compute_capacity(lastn,transfrag[t],transfrag[t]->abundance,nodecapacity,node2path,no2gnode);
						}
					}
				}
			}
		}
	}

	float cov=0;
	int len=0;
	CGraphnode *prevnode=NULL;
	// compute coverage
	for(int i=1;i<path.Count();i++) {

		if(path[i]<gno-1) {
			CGraphnode *node=no2gnode[path[i]];
			len+=node->end-node->start+1;
			float usedcov=node->cov;
			if(node->cov) {
				if(node->capacity) usedcov*=nodecapacity[i]/node->capacity;
				else if(prevnode && prevnode->capacity) usedcov*=nodecapacity[i-1]/prevnode->capacity;
				else {
					float rate=find_capacity(i,path,nodecapacity,no2gnode);
					usedcov*=rate;
				}
			}

			//if(node->end>=node->start) fprintf(stderr,"node[%d]=%f ",path[i],usedcov/(node->end-node->start+1));

			cov+=usedcov;
			prevnode=node;
		}
	}
	if(len) cov/=len;

	// clean up
	delete [] capacity;
	delete [] flow;
	delete [] link;

	return(cov);
}
*/

void get_rate(int n1, int n2,GVec<CNetEdge>& edg,GVec<float> *capacity,GVec<float> *rate,float noderate) {
	int k=0;
	int n=edg.Count();
	float abundance=capacity[n1][n2];
	capacity[n1][n2]=0;
	rate[n1][n2]=0;
	while(abundance && k<n) {
		int n3=edg[k].link;
		if(rate[n3][n1]) {
			float available=capacity[n3][n1]*noderate/rate[n3][n1];
			if(available<abundance) {
				capacity[n1][n2]+=capacity[n3][n1];
				abundance-=available;
				rate[n1][n2]+=available;
			}
			else {
				capacity[n1][n2]+=abundance*rate[n3][n1]/noderate;
				rate[n1][n2]+=abundance;
				abundance=0;
			}
		}
		else break;
		k++;
	}
	if(rate[n1][n2]) rate[n1][n2]=capacity[n1][n2]/rate[n1][n2];
}

/*
// I don't use this one
float weight_max_flow_partial(int lastn,int gno,GVec<int>& path,GBitVec& istranscript,GPVec<CTransfrag>& transfrag,
		GPVec<CGraphnode>& no2gnode,GVec<int>& node2path,GVec<float>& nodecapacity) {

	float flux=0;
	int n=path.Count();
	GVec<float> *capacity=new GVec<float>[n]; // capacity of edges in network
	GVec<float> *flow=new GVec<float>[n]; // flow in network
	GVec<float> *rate=new GVec<float>[n]; // edge rates
	GVec<int> *link=new GVec<int>[n]; // for each node remembers it's neighbours
	GVec<int> pred; // this stores the augmenting path
	pred.Resize(n,-1);
	nodecapacity.Clear();

	for(int i=0;i<n;i++) {
		node2path[path[i]]=i;
		nodecapacity.cAdd(0.0);
		capacity[i].Resize(n,0);
		flow[i].Resize(n,0);
		rate[i].Resize(n,1);
	}


	// establish capacities in the network
	for(int i=0;i<n;i++) {
		int nt=no2gnode[path[i]]->trf.Count();
		for(int j=0;j<nt;j++) {
			int t=no2gnode[path[i]]->trf[j];
			if(istranscript[t]) {
				if(transfrag[t]->nodes[0]==path[i]) { // transfrag starts at this node
					int n1=i;
					int n2=node2path[transfrag[t]->nodes.Last()];
					if(transfrag[t]->nodes.Last()>lastn) n2=n-1;
					else if(!no2gnode[transfrag[t]->nodes.Last()]->rate) n2=n-1;
					if(!no2gnode[path[i]]->rate) n1=0;
					if(!capacity[n1][n2]) {
						link[n1].Add(n2);
						link[n2].Add(n1);
					}
					capacity[n1][n2]+=transfrag[t]->abundance;
					capacity[n1][n1]+=transfrag[t]->abundance;
				}
			}
		}
	}

	// Now compute the rates and capacities
	for(int n1=1;n1<n;n1++) {
		GVec<CNetEdge> sortedg;
		for(int n2=0;n2<link[n1].Count();n2++) if(capacity[link[n1][n2]][n1]) { // incoming edge
			CNetEdge e(link[n1][n2],rate[link[n1][n2]][n1]);
			sortedg.Add(e);
		}
		sortedg.Sort(edgeCmp);
		get_rate(n1,n1,sortedg,capacity,rate,no2gnode[path[n1]]->rate);
		for(int n2=0;n2<link[n1].Count();n2++) if(capacity[n1][link[n1][n2]]) // outgoing edge
			get_rate(n1,link[n1][n2],sortedg,capacity,rate,no2gnode[path[n1]]->rate);
	}

	while(weight_bfs(n,capacity,flow,link,pred)) {
		float increment=FLT_MAX;
		for(int u=n-1;pred[u]>=0;u=pred[u]) {
			float adjflux=(capacity[pred[u]][u]-flow[pred[u]][u]);
			increment = increment < adjflux ? increment : adjflux;
			if(pred[u]<u) {
				adjflux=capacity[pred[u]][pred[u]]-flow[pred[u]][pred[u]];
				increment = increment < adjflux ? increment : adjflux; // don't allow to go over the node capacity
			}
		}

		for(int u=n-1;pred[u]>=0;u=pred[u]) {
			flow[pred[u]][u]+=increment;
			flow[u][pred[u]]-=increment;
			if(pred[u]<u) flow[pred[u]][pred[u]]+=increment;
			else flow[u][u]-=increment;
		}
		flux+=increment;
	}

	//fprintf(stderr,"flux=%f ",flux);

	for(int i=0;i<n;i++) {
		int nt=no2gnode[path[i]]->trf.Count();
		for(int j=0;j<nt;j++) {
			int t=no2gnode[path[i]]->trf[j];
			if(istranscript[t] && transfrag[t]->abundance) {
				if(transfrag[t]->nodes[0]==path[i]) { // transfrag starts at this node
					int n1=i;
					int n2=node2path[transfrag[t]->nodes.Last()];
					if(transfrag[t]->nodes.Last()>lastn) n2=n-1;
					else if(!no2gnode[transfrag[t]->nodes.Last()]->rate) n2=n-1;
					if(!no2gnode[path[i]]->rate) n1=0;
					if(flow[n1][n2]>0) {
						float flown1n2=flow[n1][n2]/rate[n1][n2];
						if(flown1n2<transfrag[t]->abundance) {
							compute_capacity(lastn,transfrag[t],flown1n2,nodecapacity,node2path,no2gnode);
							flow[n1][n2]=0;
						}
						else {
							flow[n1][n2]-=transfrag[t]->abundance*rate[n1][n2];
							compute_capacity(lastn,transfrag[t],transfrag[t]->abundance,nodecapacity,node2path,no2gnode);
						}
					}
				}
			}
		}
	}

	float cov=0;
	int len=0;
	CGraphnode *prevnode=NULL;
	// compute coverage
	for(int i=1;i<path.Count();i++) {

		if(path[i]<gno-1) {
			CGraphnode *node=no2gnode[path[i]];
			len+=node->end-node->start+1;
			float usedcov=node->cov;
			if(node->cov) {
				if(node->capacity) usedcov*=nodecapacity[i]/node->capacity;
				else if(prevnode && prevnode->capacity) usedcov*=nodecapacity[i-1]/prevnode->capacity;
				else {
					float rate=find_capacity(i,path,nodecapacity,no2gnode);
					usedcov*=rate;
				}
			}

			//if(node->end>=node->start) fprintf(stderr,"node[%d]=%f ",path[i],usedcov/(node->end-node->start+1));

			cov+=usedcov;
			prevnode=node;
		}
	}
	if(len) cov/=len;

	// clean up
	delete [] capacity;
	delete [] flow;
	delete [] link;
	delete [] rate;

	return(cov);
}
*/

/*
// I don't use this one
float max_flow_partial_back(int firstn,GVec<int>& path,GBitVec& istranscript,GPVec<CTransfrag>& transfrag,GPVec<CGraphnode>& no2gnode,
		GVec<int>& node2path,GVec<float>& nodecapacity) {

	float flux=0;
	int n=path.Count();
	GVec<float> *capacity=new GVec<float>[n]; // capacity of edges in network
	GVec<float> *flow=new GVec<float>[n]; // flow in network
	GVec<int> *link=new GVec<int>[n]; // for each node remembers it's neighbours
	GVec<int> pred; // this stores the augmenting path
	pred.Resize(n,-1);
	nodecapacity.Clear();

	for(int i=0;i<n;i++) {
		node2path[path[n-1-i]]=i;
		nodecapacity.cAdd(0.0);
		capacity[i].Resize(n,0);
		flow[i].Resize(n,0);
	}

	//fprintf(stderr,"Capacities from transcripts: ");
	// establish capacities in the network
	for(int i=0;i<n;i++) {
		int nt=no2gnode[path[n-1-i]]->trf.Count();
		for(int j=0;j<nt;j++) {
			int t=no2gnode[path[n-1-i]]->trf[j];
			if(istranscript[t]) {
				if((!i && transfrag[t]->nodes[0]<path[n-1]) || transfrag[t]->nodes[0]==path[n-1-i]) { // transfrag starts at this node

					int n1=node2path[transfrag[t]->nodes[0]];
					if(transfrag[t]->nodes[0]<path[n-1]) n1=0;
					int n2=node2path[transfrag[t]->nodes.Last()];

					//fprintf(stderr,"%d(%d-%d) ",t,n1,n2);

					if(!no2gnode[path[n-1-i]]->rate) n1=0;
					if(!no2gnode[transfrag[t]->nodes.Last()]->rate) n2=n-1;
					if(!capacity[n1][n2]) {
						link[n1].Add(n2);
						link[n2].Add(n1);
					}
					capacity[n1][n2]+=transfrag[t]->abundance;
				}
			}
		}
	}

	GVec<float> rate;
	rate.Resize(n,1);

	while(bfs(n,capacity,flow,link,pred)) {
		int r=0;
		float increment=FLT_MAX;
		rate[r++]=1;
		for(int u=n-1;pred[u]>=0;u=pred[u]) {
			float adjflux=(capacity[pred[u]][u]-flow[pred[u]][u])*rate[r-1];
			increment = increment < adjflux ? increment : adjflux;
			if(pred[pred[u]]>=0) {
				if(pred[u]<u) {
					if(pred[pred[u]]<pred[u]) rate[r]=rate[r-1]*no2gnode[path[n-1-pred[u]]]->rate;
					else rate[r]=rate[r-1];
				}
				else {
					if(pred[pred[u]]<pred[u]) rate[r]=rate[r-1];
					else rate[r]=rate[r-1]/no2gnode[path[n-1-pred[u]]]->rate;
				}
				r++;
			}
		}
		r=0;
		for(int u=n-1;pred[u]>=0;u=pred[u]) {
			flow[pred[u]][u]+=increment/rate[r];
			flow[u][pred[u]]-=increment/rate[r];
			r++;
		}
		flux+=increment;
	}

	//fprintf(stderr,"flux=%f ",flux);

	for(int i=0;i<n;i++) {
		int nt=no2gnode[path[n-1-i]]->trf.Count();
		for(int j=0;j<nt;j++) {
			int t=no2gnode[path[n-1-i]]->trf[j];
			if(istranscript[t] && transfrag[t]->abundance) {
				if((!i && transfrag[t]->nodes[0]<path[n-1]) || transfrag[t]->nodes[0]==path[n-1-i]) { // transfrag starts at this node
					int n1=node2path[transfrag[t]->nodes[0]];
					if(transfrag[t]->nodes[0]<path[n-1]) n1=0;
					int n2=node2path[transfrag[t]->nodes.Last()];
					if(!no2gnode[path[n-1-i]]->rate) n1=0;
					if(!no2gnode[transfrag[t]->nodes.Last()]->rate) n2=n-1;
					if(flow[n1][n2]>0) {
						if(flow[n1][n2]<transfrag[t]->abundance) {
							compute_capacity_back(firstn,transfrag[t],flow[n1][n2],nodecapacity,node2path);
							flow[n1][n2]=0;
						}
						else {
							flow[n1][n2]-=transfrag[t]->abundance;
							compute_capacity_back(firstn,transfrag[t],transfrag[t]->abundance,nodecapacity,node2path);
						}
					}
				}
			}
		}
	}

	float cov=0;
	int len=0;
	CGraphnode *prevnode=NULL;
	// compute coverage
	for(int i=0;i<n-1;i++) {

		if(path[n-1-i]) {
			CGraphnode *node=no2gnode[path[n-1-i]];
			len+=node->end-node->start+1;
			float usedcov=node->cov;
			if(node->cov) {
				if(node->capacity) usedcov*=nodecapacity[i]/node->capacity;
				else if(prevnode && prevnode->capacity) usedcov*=nodecapacity[i-1]/prevnode->capacity;
				else {
					float rate=find_capacity_back(i,path,nodecapacity,no2gnode);
					usedcov*=rate;
				}
			}

			cov+=usedcov;
			prevnode=node;
		}
	}
	if(len) cov/=len;

	// clean up
	delete [] capacity;
	delete [] flow;
	delete [] link;

	return(cov);
}
*/


float max_flow(int gno,GVec<int>& path,GBitVec& istranscript,GPVec<CTransfrag>& transfrag,GPVec<CGraphnode>& no2gnode,
		GVec<float>& nodecapacity,GBitVec& pathpat) { //,float& fragno) {

	float flux=0;
	int n=path.Count();
	GVec<float> *capacity=new GVec<float>[n]; // capacity of edges in network
	GVec<float> *flow=new GVec<float>[n]; // flow in network
	GVec<int> *link=new GVec<int>[n]; // for each node remembers it's neighbours
	GVec<int> pred; // this stores the augmenting path
	pred.Resize(n,-1);
	GVec<int> node2path;
	node2path.Resize(gno,-1);

	/*
	{ // DEBUG ONLY
		printTime(stderr);
		fprintf(stderr,"Start max flow algorithm for path ");
		printBitVec(pathpat);
		fprintf(stderr," :");
		for(int i=0;i<n;i++) fprintf(stderr," %d:%d",i,path[i]);
		fprintf(stderr,"\n");
		fprintf(stderr,"Used transcripts:");
		for(int i=0;i<transfrag.Count();i++) if(istranscript[i]) fprintf(stderr," %d",i);
		fprintf(stderr,"\n");
	}
	*/

	for(int i=0;i<n;i++) {
		node2path[path[i]]=i;
		nodecapacity.cAdd(0.0);
		capacity[i].Resize(n,0);
		flow[i].Resize(n,0);
	}

	// establish capacities in the network
	for(int i=0;i<n;i++) {
		int nt=no2gnode[path[i]]->trf.Count();
		for(int j=0;j<nt;j++) {
			int t=no2gnode[path[i]]->trf[j];
			if(transfrag[t]->abundance && (istranscript[t] || ((pathpat & transfrag[t]->pattern)==transfrag[t]->pattern))) {
				istranscript[t]=1;
				if(transfrag[t]->nodes[0]==path[i]) { // transfrag starts at this node
					int n1=i;
					int n2=node2path[transfrag[t]->nodes.Last()];
					if(!no2gnode[path[i]]->rate) n1=0;
					if(!no2gnode[transfrag[t]->nodes.Last()]->rate) n2=n-1;
					//fprintf(stderr,"t=%d n1=%d n2=%d ",t,n1,n2);
					if(!capacity[n1][n2]) { // haven't seen this link before
						link[n1].Add(n2);
						link[n2].Add(n1);
					}
					capacity[n1][n2]+=transfrag[t]->abundance;
					//fprintf(stderr,"capacity[%d][%d]=%f after adding transfrag[%d]->abundance=%f\n",n1,n2,capacity[n1][n2],t,transfrag[t]->abundance);
				}
			}
		}
	}

	for(int i=0;i<n;i++) link[i].Sort();

	GVec<float> rate;
	rate.Resize(n,1);

	while(bfs(n,capacity,flow,link,pred)) {
		int r=0;
		float increment=FLT_MAX;
		rate[r++]=1;
		for(int u=n-1;pred[u]>=0;u=pred[u]) {
			float adjflux=(capacity[pred[u]][u]-flow[pred[u]][u])*rate[r-1];
			increment = increment < adjflux ? increment : adjflux; // minimum flux increment on the path
			if(pred[pred[u]]>=0) {
				if(pred[u]<u) {
					if(pred[pred[u]]<pred[u]) rate[r]=rate[r-1]*no2gnode[path[pred[u]]]->rate;
					else rate[r]=rate[r-1];
				}
				else {
					if(pred[pred[u]]<pred[u]) rate[r]=rate[r-1];
					else rate[r]=rate[r-1]/no2gnode[path[pred[u]]]->rate;
				}
				r++;
			}
		}
		r=0;
		for(int u=n-1;pred[u]>=0;u=pred[u]) {
			flow[pred[u]][u]+=increment/rate[r];
			flow[u][pred[u]]-=increment/rate[r];
			r++;
		}
		flux+=increment;
	}

	/*
	{ // DEBUG ONLY
		fprintf(stderr,"Flow:");
		for(int n1=0;n1<n;n1++)
			for(int n2=n1+1;n2<n;n2++) if(flow[n1][n2]) fprintf(stderr," [%d][%d]=%f",n1,n2,flow[n1][n2]);
		fprintf(stderr,"\n");
	}
	*/

	// adjust transfrag abundances
	for(int i=0;i<n;i++) {
		int nt=no2gnode[path[i]]->trf.Count();
		float sumout=0;
		int pos=-1;
		for(int j=0;j<nt;j++) {
			int t=no2gnode[path[i]]->trf[j];
			if(istranscript[t] && transfrag[t]->abundance) {
				if(transfrag[t]->nodes[0]==path[i]) { // transfrag starts at this node
					int n1=i;
					int n2=node2path[transfrag[t]->nodes.Last()];
					if(!no2gnode[path[i]]->rate) n1=0;
					if(!no2gnode[transfrag[t]->nodes.Last()]->rate) n2=n-1;
					if(flow[n1][n2]>0) {
						if(flow[n1][n2]<transfrag[t]->abundance) {
							if(!i) sumout+=flow[n1][n2];
							update_capacity(0,transfrag[t],flow[n1][n2],nodecapacity,node2path);
							//if(path[i] && transfrag[t]->nodes.Last()!=gno-1) fragno+=flow[n1][n2];
							flow[n1][n2]=0;
						}
						else {
							if(!i) sumout+=transfrag[t]->abundance;
							flow[n1][n2]-=transfrag[t]->abundance;
							//if(path[i] && transfrag[t]->nodes.Last()!=gno-1) fragno+=transfrag[t]->abundance;
							update_capacity(0,transfrag[t],transfrag[t]->abundance,nodecapacity,node2path);
						}
					}
				}
				else if(!i && transfrag[t]->nodes.Last()==path[i]) pos=j;
			}
		}
		if(!i && pos>-1) { // this is first node -> adjust entering transfrag
			int t=no2gnode[path[i]]->trf[pos];
			float val=sumout/no2gnode[path[i]]->rate;
			transfrag[t]->abundance-=val;
			if(transfrag[t]->abundance<epsilon) transfrag[t]->abundance=0;
		}
	}

	// clean up
	delete [] capacity;
	delete [] flow;
	delete [] link;

	return(flux);
}

float guide_max_flow(bool adjust,GVec<int>& path,GBitVec& istranscript,GPVec<CTransfrag>& transfrag,GPVec<CGraphnode>& no2gnode,
		GVec<float>& nodecapacity,GBitVec& pathpat,GVec<float> *capacity,GVec<float> *flow,GVec<int> *link,GVec<int>& node2path) {
		//float& fragno) {

	float flux=0;
	int n=path.Count();

	/*
	{ // DEBUG ONLY
		fprintf(stderr,"Flow & Capacity:");
		for(int n1=0;n1<n;n1++)
			for(int n2=0;n2<link[n1].Count();n2++)
				if(link[n1][n2]>n1) fprintf(stderr," flow[%d][%d]=%g capacity[%d][%d]=%g",n1,link[n1][n2],flow[n1][link[n1][n2]],n1,link[n1][n2],capacity[n1][link[n1][n2]]);
		fprintf(stderr,"\n");
	}
	*/

	for(int i=0;i<n;i++) {
		nodecapacity.cAdd(0.0);
		if(adjust) for(int j=0;j<link[i].Count();j++) flow[i][link[i][j]]=0;
	}


	if(adjust) { // recompute the flow
		GVec<int> pred; // this stores the augmenting path
		pred.Resize(n,-1);

		GVec<float> rate;
		rate.Resize(n,1);

		while(bfs(n,capacity,flow,link,pred)) {
			int r=0;
			float increment=FLT_MAX;
			rate[r++]=1;
			for(int u=n-1;pred[u]>=0;u=pred[u]) {
				float adjflux=(capacity[pred[u]][u]-flow[pred[u]][u])*rate[r-1];
				increment = increment < adjflux ? increment : adjflux; // minimum flux increment on the path
				if(pred[pred[u]]>=0) {
					if(pred[u]<u) {
						if(pred[pred[u]]<pred[u]) rate[r]=rate[r-1]*no2gnode[path[pred[u]]]->rate;
						else rate[r]=rate[r-1];
					}
					else {
						if(pred[pred[u]]<pred[u]) rate[r]=rate[r-1];
						else rate[r]=rate[r-1]/no2gnode[path[pred[u]]]->rate;
					}
					r++;
				}
			}
			r=0;
			for(int u=n-1;pred[u]>=0;u=pred[u]) {
				flow[pred[u]][u]+=increment/rate[r];
				flow[u][pred[u]]-=increment/rate[r];
				r++;
			}
			flux+=increment;
		}

		/*
		{ // DEBUG ONLY
			fprintf(stderr,"Flow:");
			for(int n1=0;n1<n;n1++)
				for(int n2=0;n2<link[n1].Count();n2++)
					if(link[n1][n2]>n1 && flow[n1][link[n1][n2]]>0) fprintf(stderr," flow[%d][%d]=%f",n1,link[n1][n2],flow[n1][link[n1][n2]]);
			fprintf(stderr,"\n");
		}
		*/

	}


	// store flow in capacities so that I don't have to modify flow
	for(int n1=0;n1<n;n1++)
		for(int n2=0;n2<link[n1].Count();n2++)
			if(link[n1][n2]>n1) capacity[n1][link[n1][n2]]=flow[n1][link[n1][n2]];


	// adjust transfrag abundances
	for(int i=0;i<n;i++) {
		int nt=no2gnode[path[i]]->trf.Count();
		float sumout=0;
		int pos=-1;
		for(int j=0;j<nt;j++) {
			int t=no2gnode[path[i]]->trf[j];
			if(transfrag[t]->abundance && (istranscript[t] || ((pathpat & transfrag[t]->pattern)==transfrag[t]->pattern))) {
				istranscript[t]=1;
				if(transfrag[t]->nodes[0]==path[i]) { // transfrag starts at this node
					int n1=i;
					int n2=node2path[transfrag[t]->nodes.Last()];
					if(!no2gnode[path[i]]->rate) n1=0;
					if(!no2gnode[transfrag[t]->nodes.Last()]->rate) n2=n-1;
					if(capacity[n1][n2]>0) {
						if(capacity[n1][n2]<transfrag[t]->abundance) {
							if(!i) sumout+=capacity[n1][n2];
							update_capacity(0,transfrag[t],capacity[n1][n2],nodecapacity,node2path);
							//if(path[i] && transfrag[t]->nodes.Last()!=gno-1) fragno+=capacity[n1][n2];
							capacity[n1][n2]=0;
						}
						else {
							if(!i) sumout+=transfrag[t]->abundance;
							capacity[n1][n2]-=transfrag[t]->abundance;
							//if(path[i] && transfrag[t]->nodes.Last()!=gno-1) fragno+=transfrag[t]->abundance;
							update_capacity(0,transfrag[t],transfrag[t]->abundance,nodecapacity,node2path);
						}
					}
				}
				else if(!i && transfrag[t]->nodes.Last()==path[i]) pos=j;
			}
		}
		if(!i && pos>-1) { // this is first node -> adjust entering transfrag
			int t=no2gnode[path[i]]->trf[pos];
			float val=sumout/no2gnode[path[i]]->rate;
			transfrag[t]->abundance-=val;
			if(transfrag[t]->abundance<epsilon) transfrag[t]->abundance=0;
		}
	}


	return(flux);

}


float guideflow(int gno,GVec<int>& path,GBitVec& istranscript,GPVec<CTransfrag>& transfrag,GPVec<CGraphnode>& no2gnode,
		GBitVec& pathpat,GVec<float> *capacity,GVec<float> *flow,GVec<int> *link,GVec<int>& node2path) {

	float flux=0;
	int n=path.Count();
	GVec<int> pred; // this stores the augmenting path
	pred.Resize(n,-1);
	node2path.Resize(gno,-1);

	//fprintf(stderr,"path: ");

	for(int i=0;i<n;i++) {
		//fprintf(stderr,"%d ",path[i]);
		node2path[path[i]]=i;
		capacity[i].Resize(n,0);
		flow[i].Resize(n,0);
	}

	//fprintf(stderr,"n=%d ",n);

	// establish capacities in the network
	for(int i=0;i<n;i++) {
		int nt=no2gnode[path[i]]->trf.Count();
		for(int j=0;j<nt;j++) {
			int t=no2gnode[path[i]]->trf[j];
			if(transfrag[t]->abundance && (istranscript[t] || ((pathpat & transfrag[t]->pattern)==transfrag[t]->pattern))) {
				istranscript[t]=1;
				if(transfrag[t]->nodes[0]==path[i]) { // transfrag starts at this node

					//fprintf(stderr,"t=%d n1=%d n2=%d ",t,i,transfrag[t]->nodes.Last());

					int n1=i;
					int n2=node2path[transfrag[t]->nodes.Last()];
					if(!no2gnode[path[i]]->rate) n1=0;
					if(!no2gnode[transfrag[t]->nodes.Last()]->rate) n2=n-1;
					if(!capacity[n1][n2]) { // haven't seen this link before
						link[n1].Add(n2);
						link[n2].Add(n1);
					}
					capacity[n1][n2]+=transfrag[t]->abundance;
					//fprintf(stderr,"capacity[%d][%d]=%f after adding transfrag[%d]->abundance=%f\n",n1,n2,capacity[n1][n2],t,transfrag[t]->abundance);
				}
			}
		}
	}

	for(int i=0;i<n;i++) link[i].Sort();


	GVec<float> rate;
	rate.Resize(n,1);

	while(bfs(n,capacity,flow,link,pred)) {
		int r=0;
		float increment=FLT_MAX;
		rate[r++]=1;
		for(int u=n-1;pred[u]>=0;u=pred[u]) {
			float adjflux=(capacity[pred[u]][u]-flow[pred[u]][u])*rate[r-1];
			increment = increment < adjflux ? increment : adjflux;
			if(pred[pred[u]]>=0) {
				if(pred[u]<u) {
					if(pred[pred[u]]<pred[u]) rate[r]=rate[r-1]*no2gnode[path[pred[u]]]->rate;
					else rate[r]=rate[r-1];
				}
				else {
					if(pred[pred[u]]<pred[u]) rate[r]=rate[r-1];
					else rate[r]=rate[r-1]/no2gnode[path[pred[u]]]->rate;
				}
				r++;
			}
		}
		r=0;
		for(int u=n-1;pred[u]>=0;u=pred[u]) {
			flow[pred[u]][u]+=increment/rate[r];
			flow[u][pred[u]]-=increment/rate[r];
			r++;
		}
		flux+=increment;
	}

	/*
	{ // DEBUG ONLY
		fprintf(stderr,"Flow:");
		for(int n1=0;n1<n;n1++)
			for(int n2=n1+1;n2<n;n2++) if(flow[n1][n2]) fprintf(stderr," [%d][%d]=%f",n1,n2,flow[n1][n2]);
		fprintf(stderr,"\n");
	}
	*/

	return(flux);
}


float max_flow_EM(int gno,GVec<int>& path,GBitVec& istranscript,GPVec<CTransfrag>& transfrag,GPVec<CGraphnode>& no2gnode,
		GVec<float>& nodecapacity,GBitVec& pathpat) {//,float &fragno) {



	float flux=0;
	int n=path.Count();
	int m=n+2;
	GVec<float> *capacity=new GVec<float>[m]; // capacity of edges in network
	GVec<float> *flow=new GVec<float>[m]; // flow in network
	GVec<int> *link=new GVec<int>[m]; // for each node remembers it's neighbours
	GVec<int> pred; // this stores the augmenting path
	pred.Resize(m,-1);
	GVec<int> node2path;
	node2path.Resize(gno,-1);

	/*
	{ // DEBUG ONLY
		printTime(stderr);
		fprintf(stderr,"Start max flow algorithm for path ");
		//printBitVec(pathpat);
		fprintf(stderr," :");
		for(int i=0;i<n;i++) fprintf(stderr," %d:%d",i,path[i]);
		fprintf(stderr,"\n");
		fprintf(stderr,"Used transcripts:");
		for(int i=0;i<transfrag.Count();i++) if(istranscript[i]) fprintf(stderr," %d",i);
		fprintf(stderr,"\n");
	}
	*/

	GVec<float> through; // these are the capacity of the "trough" transfrags through each node in the path
	through.Resize(n,0);

	for(int i=0;i<m;i++) {
		if(i<n) node2path[path[i]]=i;
		if(i<n) nodecapacity.cAdd(0.0);
		capacity[i].Resize(m,0);
		flow[i].Resize(m,0);
	}

	// establish capacities in the network
	for(int i=0;i<n;i++) {
		int nt=no2gnode[path[i]]->trf.Count();
		for(int j=0;j<nt;j++) {
			int t=no2gnode[path[i]]->trf[j];
			if(transfrag[t]->abundance && (istranscript[t] || ((pathpat & transfrag[t]->pattern)==transfrag[t]->pattern))) {
				istranscript[t]=1;
				if(transfrag[t]->nodes[0]==path[i]) { // transfrag starts at this node
					int n1=i;
					int n2=node2path[transfrag[t]->nodes.Last()];
					if(!capacity[n1][n2]) { // haven't seen this link before
						link[n1].Add(n2);
						link[n2].Add(n1);
					}
					capacity[n1][n2]+=transfrag[t]->abundance;
				}
				else if(transfrag[t]->nodes[0]<path[i] && transfrag[t]->nodes.Last()>path[i] && transfrag[t]->pattern[path[i]]) { // through transfrag
					through[i]+=transfrag[t]->abundance;
				}
			}
		}

		if(i && i<n-1 && through[i]) { // not source or sink and I have transfrags going through the node
			// 0 -> n : source links to fake node n
			link[n].Add(i);
			link[i].Add(n);

			// n+1 -> sink : fake node n+1 links to sink
			int n1=n+1;
			link[n1].Add(i);
			link[i].Add(n1);

			capacity[n][i]+=through[i];
			capacity[i][n1]+=through[i];

			int sink=n-1;
			if(!capacity[n1][sink]) {
				link[n1].Add(sink);
				link[sink].Add(n1);
			}
			capacity[n1][sink]+=through[i];

			if(!capacity[0][n]) {
				link[0].Add(n);
				link[n].cAdd(0);
			}
			capacity[0][n]+=through[i];
			capacity[0][0]+=through[i];
		}

	}

	for(int i=0;i<n;i++) link[i].Sort();

	bool doEM=true;
	int iterations=0;
	GHash<float> tabund;

	GVec<float> rate;
	rate.Resize(m,1);

	while(doEM && iterations<10) {
		flux=0;
		while(bfs(n,capacity,flow,link,pred)) {
			int r=0;
			float increment=FLT_MAX;
			rate[r++]=1;
			for(int u=n-1;pred[u]>=0;u=pred[u]) {
				float adjflux=(capacity[pred[u]][u]-flow[pred[u]][u])*rate[r-1];
				increment = increment < adjflux ? increment : adjflux;
				if(pred[pred[u]]>=0) {
					if(pred[u]>=n) rate[r]=rate[r-1];
					else
						if(pred[u]<u) {
							if(pred[pred[u]]<pred[u]) rate[r]=rate[r-1]*no2gnode[path[pred[u]]]->rate;
							else rate[r]=rate[r-1];
						}
						else {
							if(pred[pred[u]]<pred[u]) rate[r]=rate[r-1];
							else rate[r]=rate[r-1]/no2gnode[path[pred[u]]]->rate;
						}
					r++;
				}
			}
			r=0;
			for(int u=n-1;pred[u]>=0;u=pred[u]) {
				flow[pred[u]][u]+=increment/rate[r];
				flow[u][pred[u]]-=increment/rate[r];
				r++;
			}
			flux+=increment;
		}

		/*
		{ // DEBUG ONLY
			printTime(stderr);
			fprintf(stderr,"Flow:");
			for(int n1=0;n1<n;n1++)
				for(int n2=n1+1;n2<n;n2++) if(flow[n1][n2]) fprintf(stderr," [%d][%d]=%f",n1,n2,flow[n1][n2]);
			fprintf(stderr,"\n");
		}
		*/

		tabund.Clear();
		doEM=false;

		for(int i=0;i<n;i++) {
			through[i]=0;
			int nt=no2gnode[path[i]]->trf.Count();
			for(int j=0;j<nt;j++) {
				int t=no2gnode[path[i]]->trf[j];
				if(istranscript[t] && transfrag[t]->abundance) {
					if(transfrag[t]->nodes[0]==path[i]) { // transfrag starts at this node
						int n1=i;
						int n2=node2path[transfrag[t]->nodes.Last()];
						if(flow[n1][n2]>0) {
							if(flow[n1][n2]<transfrag[t]->abundance) {
								GStr tid(t);
								tabund.Add(tid.chars(),new float(flow[n1][n2]));
								flow[n1][n2]=0;
							}
							else {
								flow[n1][n2]-=transfrag[t]->abundance;
								GStr tid(t);
								tabund.Add(tid.chars(),new float(transfrag[t]->abundance));
							}
						}
					}
					else if(transfrag[t]->nodes[0]<path[i] && transfrag[t]->nodes.Last()>path[i] && transfrag[t]->pattern[path[i]]) { // through transfrag
						GStr tid(t);
						const float *abund=tabund[tid.chars()];
						if(abund) through[i]+= *abund;
					}
				}
			}
			// now check if I should continue the EM algorithm
			if(flow[n][i]-through[i]>epsilon) {
				doEM=true;
				capacity[n][i]=through[i];
			}
			if(flow[i][n+1]-through[i]>epsilon) {
				doEM=true;
				capacity[i][n+1]=through[i];
			}

		}


		if(doEM)  // reset flow to 0
			for(int i=0;i<m;i++) flow[i].Resize(m,0);

		iterations++;

	}

	// adjust transfrag abundances
	for(int i=0;i<n;i++) {
		int nt=no2gnode[path[i]]->trf.Count();
		for(int j=0;j<nt;j++) {
			int t=no2gnode[path[i]]->trf[j];
			if(istranscript[t] && transfrag[t]->abundance && transfrag[t]->nodes[0]==path[i]) {
				GStr tid(t);
				const float *abund=tabund[tid.chars()];
				if(abund) {
					update_capacity(0,transfrag[t],*abund,nodecapacity,node2path);
					//if(path[i] && transfrag[t]->nodes.Last()!=gno-1) fragno+=*abund;
				}
			}
		}
	}

	// clean up
	delete [] capacity;
	delete [] flow;
	delete [] link;

	return(flux);
}




float weight_max_flow(int gno,GVec<int>& path,GBitVec& istranscript,GPVec<CTransfrag>& transfrag,GPVec<CGraphnode>& no2gnode,
		GVec<float>& nodecapacity,GBitVec& pathpat) {//,float& fragno) {

	int n=path.Count();

	/*
	{ // DEBUG ONLY
		printTime(stderr);
		fprintf(stderr,"Start max flow algorithm for path ");
		printBitVec(pathpat);
		fprintf(stderr," :");
		for(int i=0;i<n;i++) fprintf(stderr," %d:%d",i,path[i]);
		fprintf(stderr,"\n");
	}
	*/

	float flux=0;
	GVec<float> *capacity=new GVec<float>[n]; // capacity of edges in network
	GVec<float> *flow=new GVec<float>[n]; // flow in network
	GVec<float> *rate=new GVec<float>[n]; // edge rates
	GVec<int> *link=new GVec<int>[n]; // for each node remembers it's neighbours
	GVec<int> pred; // this stores the augmenting path
	pred.Resize(n,-1);
	GVec<int> node2path;
	node2path.Resize(gno,-1);

	for(int i=0;i<n;i++) {
		node2path[path[i]]=i;
		nodecapacity.cAdd(0.0);
		capacity[i].Resize(n,0);
		flow[i].Resize(n,0);
		rate[i].Resize(n,1);
	}

	// establish capacities in the network
	for(int i=0;i<n;i++) {
		int nt=no2gnode[path[i]]->trf.Count();
		for(int j=0;j<nt;j++) {
			int t=no2gnode[path[i]]->trf[j];
			if(transfrag[t]->abundance && (istranscript[t] || ((pathpat & transfrag[t]->pattern)==transfrag[t]->pattern))) {
				istranscript[t]=1;
				if(transfrag[t]->nodes[0]==path[i]) { // transfrag starts at this node
					int n1=i;
					int n2=node2path[transfrag[t]->nodes.Last()];
					if(!no2gnode[path[i]]->rate) n1=0;
					if(!no2gnode[transfrag[t]->nodes.Last()]->rate) n2=n-1;
					if(!capacity[n1][n2]) { // haven't seen this link before
						link[n1].Add(n2);
						link[n2].Add(n1);
					}
					capacity[n1][n2]+=transfrag[t]->abundance;
					capacity[n1][n1]+=transfrag[t]->abundance;
				}
			}
		}
	}

	/*
	{ //DEBUG ONLY
		printTime(stderr);
		fprintf(stderr,"Abundances:");
		for(int n1=0;n1<n;n1++) {
			for(int n2=n1+1;n2<n;n2++) if(capacity[n1][n2]) {
				fprintf(stderr," [%d][%d]=%f",n1,n2,capacity[n1][n2]);
			}
			fprintf(stderr," tr=");
			for(int j=0;j<no2gnode[path[n1]]->trf.Count();j++) {
				int t=no2gnode[path[n1]]->trf[j];
				if(transfrag[t]->nodes[0]==path[n1] && istranscript[t]) fprintf(stderr," %d(->%d)",t,transfrag[t]->nodes.Last());
			}
		}
		fprintf(stderr,"\n");
	}
	*/

	// Now compute the rates and capacities
	for(int n1=1;n1<n;n1++) {
		GVec<CNetEdge> sortedg;
		for(int n2=0;n2<link[n1].Count();n2++) if(capacity[link[n1][n2]][n1]) { // incoming edge
				CNetEdge e(link[n1][n2],rate[link[n1][n2]][n1]);
				sortedg.Add(e);
		}
		sortedg.Sort(edgeCmp); // largest rate comes first
		get_rate(n1,n1,sortedg,capacity,rate,no2gnode[path[n1]]->rate);
		for(int n2=0;n2<link[n1].Count();n2++) if(capacity[n1][link[n1][n2]]) // outgoing edge
			get_rate(n1,link[n1][n2],sortedg,capacity,rate,no2gnode[path[n1]]->rate);
	}

	/*
	{ //DEBUG ONLY
		printTime(stderr);
		fprintf(stderr,"Capacities:");
		for(int n1=0;n1<n;n1++) {
			fprintf(stderr," rate[%d]=%f",path[n1],no2gnode[path[n1]]->rate);
			for(int n2=n1;n2<n;n2++) if(capacity[n1][n2]) {
				fprintf(stderr," [%d][%d]=%f(%f)",n1,n2,capacity[n1][n2],rate[n1][n2]);
			}
		}
		fprintf(stderr,"\n");
	}
	*/

	while(weight_bfs(n,capacity,flow,link,pred)) {
		float increment=FLT_MAX;
		for(int u=n-1;pred[u]>=0;u=pred[u]) {
			float adjflux=capacity[pred[u]][u]-flow[pred[u]][u];
			increment = increment < adjflux ? increment : adjflux;
			if(pred[u]<u) {
				adjflux=capacity[pred[u]][pred[u]]-flow[pred[u]][pred[u]];
				increment = increment < adjflux ? increment : adjflux; // don't allow to go over the node capacity
			}
		}
		for(int u=n-1;pred[u]>=0;u=pred[u]) {
			flow[pred[u]][u]+=increment;
			flow[u][pred[u]]-=increment;
			if(pred[u]<u) flow[pred[u]][pred[u]]+=increment;
			else flow[u][u]-=increment;
		}
		flux+=increment;
	}

	/*
	{ // DEBUG ONLY
		printTime(stderr);
		fprintf(stderr,"Flow:");
		for(int n1=0;n1<n;n1++)
			for(int n2=n1+1;n2<n;n2++) if(flow[n1][n2]) fprintf(stderr," [%d][%d]=%f(%f)",n1,n2,flow[n1][n2],flow[n1][n2]/rate[n1][n2]);
		fprintf(stderr,"\n");
	}
	*/

	// adjust transfrag abundances
	for(int i=0;i<n;i++) {
		int nt=no2gnode[path[i]]->trf.Count();
		float sumout=0;
		int pos=-1;
		for(int j=0;j<nt;j++) {
			int t=no2gnode[path[i]]->trf[j];
			if(istranscript[t] && transfrag[t]->abundance) {
				if(transfrag[t]->nodes[0]==path[i]) { // transfrag starts at this node
					int n1=i;
					int n2=node2path[transfrag[t]->nodes.Last()];
					if(!no2gnode[path[i]]->rate) n1=0;
					if(!no2gnode[transfrag[t]->nodes.Last()]->rate) n2=n-1;
					if(flow[n1][n2]>0) {
						float flown1n2=flow[n1][n2]/rate[n1][n2];
						if(flown1n2<transfrag[t]->abundance) {
							if(!i) sumout+=flown1n2;
							update_capacity(0,transfrag[t],flown1n2,nodecapacity,node2path);
							//if(path[i] && transfrag[t]->nodes.Last()!=gno-1) fragno+=flown1n2;
							flow[n1][n2]=0;
						}
						else {
							if(!i) sumout+=transfrag[t]->abundance;
							flow[n1][n2]-=transfrag[t]->abundance*rate[n1][n2];
							//if(path[i] && transfrag[t]->nodes.Last()!=gno-1) fragno+=transfrag[t]->abundance;
							update_capacity(0,transfrag[t],transfrag[t]->abundance,nodecapacity,node2path);
						}
					}
				}
				else if(!i && transfrag[t]->nodes.Last()==path[i]) pos=j; // NOTE: this will never work if the pathpat doesn't include the link to source, because the transcript linking back to source is not on the path
			}
		}
		if(!i && pos>-1) { // this is first node -> adjust entering transfrag
			int t=no2gnode[path[i]]->trf[pos];
			float val=sumout/no2gnode[path[i]]->rate;
			transfrag[t]->abundance-=val;
			if(transfrag[t]->abundance<epsilon) transfrag[t]->abundance=0;
		}
	}

	// clean up
	delete [] capacity;
	delete [] flow;
	delete [] link;
	delete [] rate;

	return(flux);
}

/*
// I don't use this one: doesn't work
float weight_max_flow_EM(int gno,GVec<int>& path,GBitVec& istranscript,GPVec<CTransfrag>& transfrag,GPVec<CGraphnode>& no2gnode,
		GVec<float>& nodecapacity,GBitVec& pathpat) {


	{ // DEBUG ONLY
		printTime(stderr);
		fprintf(stderr,"Start max flow EM algorithm.\n");
	}


	float flux=0;
	int n=path.Count();
	int m=n+2; // this includes the two fake nodes I need to add
	GVec<float> *capacity=new GVec<float>[m]; // capacity of edges in network
	GVec<float> *flow=new GVec<float>[m]; // flow in network
	GVec<float> *rate=new GVec<float>[m]; // edge rates
	GVec<int> *link=new GVec<int>[m]; // for each node remembers it's neighbours
	GVec<int> pred; // this stores the augmenting path
	pred.Resize(m,-1);
	GVec<int> node2path;
	node2path.Resize(gno,-1);

	GVec<float> through; // these are the capacity of the "trough" transfrags through each node in the path
	through.Resize(n,0);

	for(int i=0;i<m;i++) {
		if(i<n) node2path[path[i]]=i;
		if(i<n) nodecapacity.cAdd(0.0);
		capacity[i].Resize(m,0);
		flow[i].Resize(m,0);
		rate[i].Resize(m,1);
	}

	// establish capacities in the network
	for(int i=0;i<n;i++) {
		int nt=no2gnode[path[i]]->trf.Count();
		for(int j=0;j<nt;j++) {
			int t=no2gnode[path[i]]->trf[j];
			if(transfrag[t]->abundance && (istranscript[t] || ((pathpat & transfrag[t]->pattern)==transfrag[t]->pattern))) {
				istranscript[t]=1;
				if(transfrag[t]->nodes[0]==path[i]) { // transfrag starts at this node
					int n1=i;
					int n2=node2path[transfrag[t]->nodes.Last()];
					if(!capacity[n1][n2]) { // haven't seen this link before
						link[n1].Add(n2);
						link[n2].Add(n1);
					}
					capacity[n1][n2]+=transfrag[t]->abundance;
					capacity[n1][n1]+=transfrag[t]->abundance; // why do I need this?
				}
				else if(transfrag[t]->nodes[0]<path[i] && transfrag[t]->nodes.Last()>path[i] && transfrag[t]->pattern[path[i]]) { // through transfrag
					through[i]+=transfrag[t]->abundance;
				}
			}
		}

		if(i && i<n-1 && through[i]) { // not source or sink and I have transfrags going through the node
			// 0 -> n : source links to fake node n  ; actually this will work only if I modify weight_bfs too to establish the pred[u]<u relationship like this: n<i for all i>0, and n+1<n
			link[n].Add(i);
			link[i].Add(n);

			// n+1 -> sink : fake node n+1 links to sink
			int n1=n+1;
			link[n1].Add(i);
			link[i].Add(n1);

			capacity[n][i]+=through[i];
			capacity[n][n]+=through[i];
			capacity[i][n1]+=through[i];
			capacity[i][i]+=through[i];

			int sink=n-1;
			if(!capacity[n1][sink]) {
				link[n1].Add(sink);
				link[sink].Add(n1);
			}
			capacity[n1][sink]+=through[i];
			capacity[n1][n1]+=through[i];

			if(!capacity[0][n]) {
				link[0].Add(n);
				link[n].cAdd(0);
			}
			capacity[0][n]+=through[i];
			capacity[0][0]+=through[i];
		}
	}


	{ //DEBUG ONLY
		printTime(stderr);
		fprintf(stderr,"Abundances:");
		for(int n1=0;n1<n;n1++) {
			for(int n2=n1+1;n2<n;n2++) if(capacity[n1][n2]) {
				fprintf(stderr," [%d][%d]=%f",n1,n2,capacity[n1][n2]);
			}
		}
		fprintf(stderr,"\n");
	}

	bool doEM=true;
	int iterations=0;
	GHash<float> tabund;

	while(doEM && iterations<10) {

		flux=0;
		// Compute the rates and capacities
		for(int n1=1;n1<m;n1++) {
			GVec<CNetEdge> sortedg;
			for(int n2=0;n2<link[n1].Count();n2++) if(capacity[link[n1][n2]][n1]) { // incoming edge
				CNetEdge e(link[n1][n2],rate[link[n1][n2]][n1]);
				sortedg.Add(e);
			}
			sortedg.Sort(edgeCmpEM); // sorts the edges with the larger rate first
			float noderate=1;
			if(n1<n) noderate=no2gnode[path[n1]]->rate;
			get_rate(n1,n1,sortedg,capacity,rate,noderate);

			for(int n2=0;n2<link[n1].Count();n2++) if(capacity[n1][link[n1][n2]]) // outgoing edge
				get_rate(n1,link[n1][n2],sortedg,capacity,rate,noderate);

		}



		{ //DEBUG ONLY
			printTime(stderr);
			fprintf(stderr,"Capacities:");
			for(int n1=0;n1<n;n1++) {
				fprintf(stderr," rate[%d]=%f",path[n1],no2gnode[path[n1]]->rate);
				for(int n2=n1;n2<n;n2++) if(capacity[n1][n2]) {
					fprintf(stderr," [%d][%d]=%f(%f)",n1,n2,capacity[n1][n2],rate[n1][n2]);
				}
			}
			fprintf(stderr,"\n");
		}



		while(weight_bfs_EM(n,capacity,flow,link,pred)) {
			float increment=FLT_MAX;
			for(int u=n-1;pred[u]>=0;u=pred[u]) {
				float adjflux=capacity[pred[u]][u]-flow[pred[u]][u];
				increment = increment < adjflux ? increment : adjflux;
				if(pred[u]<u) {
					adjflux=capacity[pred[u]][pred[u]]-flow[pred[u]][pred[u]];
					increment = increment < adjflux ? increment : adjflux; // don't allow to go over the node capacity
				}
			}
			for(int u=n-1;pred[u]>=0;u=pred[u]) {
				flow[pred[u]][u]+=increment;
				flow[u][pred[u]]-=increment;
				if(pred[u]<u) flow[pred[u]][pred[u]]+=increment;
				else flow[u][u]-=increment;
			}
			flux+=increment;
		}


		{ // DEBUG ONLY
			printTime(stderr);
			fprintf(stderr,"Flow:");
			for(int n1=0;n1<n;n1++)
				for(int n2=n1+1;n2<n;n2++) if(flow[n1][n2]) fprintf(stderr," [%d][%d]=%f(%f)",n1,n2,flow[n1][n2],flow[n1][n2]/rate[n1][n2]);
			fprintf(stderr,"\n");
		}


		tabund.Clear();
		doEM=false;

		for(int i=0;i<n;i++) {
			through[i]=0;
			int nt=no2gnode[path[i]]->trf.Count();
			for(int j=0;j<nt;j++) {
				int t=no2gnode[path[i]]->trf[j];
				if(istranscript[t] && transfrag[t]->abundance) {
					if(transfrag[t]->nodes[0]==path[i]) { // transfrag starts at this node
						int n1=i;
						int n2=node2path[transfrag[t]->nodes.Last()];
						if(flow[n1][n2]>0) {
							float flown1n2=flow[n1][n2]/rate[n1][n2];
							if(flown1n2<transfrag[t]->abundance) {
								GStr tid(t);
								tabund.Add(tid.chars(),new float(flown1n2));
								flow[n1][n2]=0;
							}
							else {
								flow[n1][n2]-=transfrag[t]->abundance*rate[n1][n2];
								GStr tid(t);
								tabund.Add(tid.chars(),new float(transfrag[t]->abundance));
							}
						}
					}
					else if(transfrag[t]->nodes[0]<path[i] && transfrag[t]->nodes.Last()>path[i] && transfrag[t]->pattern[path[i]]) { // through transfrag
						GStr tid(t);
						const float *abund=tabund[tid.chars()];
						if(abund) through[i]+= *abund;
					}
				}
			}
			// now check if I should continue the EM algorithm
			if(through[i]<flow[n][i]) {
				doEM=true;
				capacity[n][i]=through[i];
			}
			if(through[i]<flow[i][n+1]) {
				doEM=true;
				capacity[i][n+1]=through[i];
			}
		}

		if(doEM)  // reset flow to 0
			for(int i=0;i<m;i++) flow[i].Resize(m,0);

		iterations++;
	}

	// adjust transfrag abundances
	for(int i=0;i<n;i++) {
		int nt=no2gnode[path[i]]->trf.Count();
		for(int j=0;j<nt;j++) {
			int t=no2gnode[path[i]]->trf[j];
			if(istranscript[t] && transfrag[t]->abundance && transfrag[t]->nodes[0]==path[i]) {
				GStr tid(t);
				const float *abund=tabund[tid.chars()];
				if(abund) update_capacity(0,transfrag[t],*abund,nodecapacity,node2path);
			}
		}
	}


	// clean up
	delete [] capacity;
	delete [] flow;
	delete [] link;
	delete [] rate;

	return(flux);
}
*/

/*
// I don't use this one
float update_flux_fast(GVec<int>& path,GBitVec& istranscript,GPVec<CTransfrag>& transfrag,GPVec<CGraphnode>& no2gnode,
		GVec<float>& capacity,GBitVec& pathpat) {

	float flux=0;

	GVec<float> used;
	used.Resize(transfrag.Count(),0);

	int n=path.Count();


	{ // DEBUG ONLY
		fprintf(stderr,"Pathpat=");
		printBitVec(pathpat);
		fprintf(stderr,"\nTransfrags:");
	}


	for(int i=0;i<n;i++) {
		capacity.cAdd(0.0);
		int nt=no2gnode[path[i]]->trf.Count();
		GVec<int> out;
		float sumout=0;
		float sumin=0;
		float sumthrough=0;
		for(int t=0;t<nt;t++) {
			if(istranscript[no2gnode[path[i]]->trf[t]] ||
					((pathpat & transfrag[no2gnode[path[i]]->trf[t]]->pattern)==transfrag[no2gnode[path[i]]->trf[t]]->pattern)) {
				istranscript[no2gnode[path[i]]->trf[t]]=1;
				if(i<n-1) {
					if(transfrag[no2gnode[path[i]]->trf[t]]->nodes.Last()==path[i]) { // transfrag enters node
						sumin+=used[no2gnode[path[i]]->trf[t]];
					}
					else if(transfrag[no2gnode[path[i]]->trf[t]]->nodes[0]==path[i]) { // transfrag exits node
						sumout+=transfrag[no2gnode[path[i]]->trf[t]]->abundance;
						out.Add(no2gnode[path[i]]->trf[t]);
					}
					else sumthrough+=transfrag[no2gnode[path[i]]->trf[t]]->abundance;
				}
			}
		}
		if(i<n-1) {
			float frac=1;
			if(i && no2gnode[path[i]]->rate) {
				sumin*=no2gnode[path[i]]->rate; // adjust sumin for rate
				if(sumin<sumout) {
					if(sumout) frac=sumin/sumout;
					sumout=sumin;
				}
			}
			for(int t=0;t<out.Count();t++) {
				used[out[t]]=frac*transfrag[out[t]]->abundance;
			}
		}
	}
	for(int i=n-2;i>=0;i--) { // skip sink
		int nt=no2gnode[path[i]]->trf.Count();
		GVec<int> in;
		float sumout=0;
		float sumin=0;
		for(int j=0;j<nt;j++) {
			int t=no2gnode[path[i]]->trf[j];
			if(istranscript[t]) {
				if(transfrag[t]->nodes.Last()==path[i] && used[t]) { // transfrag enters node
					sumin+=used[t];
					in.Add(t);
				}
				else if(transfrag[t]->nodes[0]==path[i]) { // transfrag exits node
					sumout+=used[t];
					flux+=used[t];
					transfrag[t]->abundance-=used[t];
					if(transfrag[t]->abundance<epsilon) {
						transfrag[t]->abundance=0;
						no2gnode[path[i]]->trf.Delete(j);
						j--;nt--;
					}
				}
				else  if(transfrag[t]->pattern[path[i]]) { // through transfrag and not incomplete
					capacity[i]+=used[t];
				}
			}
		}

		capacity[i]+=sumout;
		float frac=1;
		if(no2gnode[path[i]]->rate) {
			sumout/=no2gnode[path[i]]->rate;
			if(sumin>sumout) {
				if(sumin) frac=sumout/sumin;
				sumin=sumout;
			}
			for(int t=0;t<in.Count();t++) {
				used[in[t]]=frac*used[in[t]];
			}
		}
		else capacity[i]+=sumin;
	}

	return(flux);
}
*/

/*
// I don't use this one
float update_flux(int gno,GVec<int>& path,GBitVec& istranscript,GPVec<CTransfrag>& transfrag,GBitVec& added,GPVec<CGraphnode>& no2gnode,
		GVec<float>& capacity,GBitVec& pathpat) {

	float abundance=0;

	int n=path.Count();

	// add transcripts that might be compatible with path and compute in and out capacities
	GVec<int> trf;
	GVec<int> node2path;
	node2path.Resize(gno,-1);

	for(int i=0;i<n;i++) {
		node2path[path[i]]=i;
		capacity.cAdd(0.0);
		int nt=no2gnode[path[i]]->trf.Count();
		for(int t=0;t<nt;t++) {
			if(istranscript[no2gnode[path[i]]->trf[t]] ||
					((pathpat & transfrag[no2gnode[path[i]]->trf[t]]->pattern)==transfrag[no2gnode[path[i]]->trf[t]]->pattern)) {
				if(!added[no2gnode[path[i]]->trf[t]]) {
					trf.Add(no2gnode[path[i]]->trf[t]);
					added[no2gnode[path[i]]->trf[t]]=1;
				}
				istranscript[no2gnode[path[i]]->trf[t]]=1;
			}
		}
	}

	trf.Sort();


	{ //DEBUG ONLY
		fprintf(stderr,"Pathpat=");
		printBitVec(pathpat);
		fprintf(stderr,"\n");
		fprintf(stderr,"Transfrags:");
		for(int i=0;i<trf.Count();i++)
			fprintf(stderr," %d(%f)",trf[i],transfrag[trf[i]]->abundance);
		fprintf(stderr,"\n");
	}


	GVec<float> used;
	used.Resize(transfrag.Count(),0);
	GVec<int> usedpos;
	usedpos.Resize(n,0);

	int i=0;
	while(i<trf.Count()) {
		int t=trf[i];
		if(istranscript[t] && transfrag[t]->abundance) { // transcript might have been deleted?
			int updateleft=n;
			int updateright=0;
			float val1=extend_tr_left(t,transfrag[t]->abundance,used,usedpos,updateleft,path,node2path[transfrag[t]->nodes[0]],no2gnode,transfrag,node2path,istranscript);
			if(val1<epsilon) val1=0;
			float val2=0;
			if(val1) val2=extend_tr_right(t,val1,used,usedpos,updateright,path,node2path[transfrag[t]->nodes.Last()],no2gnode,transfrag,node2path,istranscript);
			if(val2) {
				eliminate_tr_left(t,val2,used,usedpos,capacity,path,node2path[transfrag[t]->nodes[0]],no2gnode,transfrag,istranscript);
				eliminate_tr_right(t,val2,used,usedpos,capacity,path,node2path[transfrag[t]->nodes.Last()],no2gnode,transfrag,istranscript);
				update_capacity(1,transfrag[t],val2,capacity,node2path);
			}
			istranscript[t]=0;
			if(val2!=val1 || !val1) {
				if(val1) for(int j=i+1;j<trf.Count();j++) used[trf[j]]=0;
				for(int j=updateleft;j<=node2path[transfrag[t]->nodes[0]];j++) usedpos[j]=0;
				for(int j=node2path[transfrag[t]->nodes.Last()];j<updateright;j++) usedpos[j]=0;
			}
			abundance+=val2;
		}
		i++;
	}

	return(abundance);
}
*/

//float store_transcript(GList<CPrediction>& pred,GVec<int>& path,GVec<float>& nodeflux,GVec<float>& nodecov,
//		GPVec<CGraphnode>& no2gnode,int& geneno,bool& first,int strand,int gno,bool& included,GBitVec& prevpath,float fragno,char* id=NULL) {
float store_transcript(GList<CPrediction>& pred,GVec<int>& path,GVec<float>& nodeflux,GVec<float>& nodecov,
		GPVec<CGraphnode>& no2gnode,int& geneno,bool& first,int strand,int gno,GIntHash<int>& gpos, bool& included,
		GBitVec& prevpath,//float fragno, char* id=NULL) {
		   GffObj* t=NULL) {
	float cov=0;
	int len=0;
	CGraphnode *prevnode=NULL;
	GVec<GSeg> exons;
	GVec<float> exoncov;
	float excov=0;

	int s=0;
	if(!path[0]) s=1;

	//if(t) fprintf(stderr,"store transcript with id=%s\n",t->getID());

	bool firstex=true;
	bool lastex=false;

	for(int i=s;i<path.Count()-1;i++) {
		if(!prevpath[path[i]]) {
			included=false;
			prevpath[path[i]]=1;
		}
		int *pos=gpos[edge(path[i-1],path[i],gno)];
		if(i && (!pos || !prevpath[*pos])) {
			included=false;
			if(pos) prevpath[*pos]=1;
		}

		CGraphnode *node=no2gnode[path[i]];

		float usedcov=node->cov;
		if(node->cov) {
			if(node->capacity) usedcov*=nodeflux[i]/node->capacity;
			else if(prevnode && prevnode->capacity) usedcov*=nodeflux[i-1]/prevnode->capacity;
			else {
				float rate=find_capacity(i,path,nodeflux,no2gnode);
				usedcov*=rate;
			}
		}
		nodecov[path[i]]-=usedcov/(node->end-node->start+1);

		if(t && (node->end<t->start || lastex)) {
			prevnode=node; continue;
		}

		uint nodestart=node->start;
		uint nodeend=node->end;
		if(t) {
			if(firstex) nodestart=t->start;
			if(t->end<=node->end || i==path.Count()-2) {
				lastex=true;
				nodeend=t->end;
			}
		}

		if(!prevnode || firstex || node->start>prevnode->end+1) { // this is a new exon
			if(prevnode && !firstex) { // compute exon coverage
				excov/=exons.Last().end-exons.Last().start+1;
				exoncov.Add(excov);
				excov=0;
			}
			//fprintf(stderr,"Add exon %d-%d\n",nodestart,nodeend);
			GSeg exon(nodestart,nodeend);
			exons.Add(exon);
			firstex=false;
		}
		else if(!firstex) exons.Last().end=nodeend;

		len+=nodeend-nodestart+1;

		cov+=usedcov;
		excov+=usedcov;
		//if(node->cov) fragno+=node->frag*usedcov/node->cov;

		prevnode=node;
	}

	/*
	{ // DEBUG ONLY
		fprintf(stderr,"Predicted transcript cov=%f usedcov=%f len=%d fragno=%g included=%d path.count=%d\n",cov/len, cov,len,fragno,included,path.Count());
		if(t) fprintf(stderr,"Ref_id=%s\n",t->getID());
	}
	*/

	// Add last exon coverage
	if(prevnode) { // compute exon coverage
		excov/=exons.Last().end-exons.Last().start+1;
		exoncov.Add(excov);
	}
	if(len) cov/=len;


	//if(id || ((!included || sensitivitylevel) && cov>=readthr && len>=mintranscriptlen)) { // store transcript here
	if(t || ((!included || sensitivitylevel) && cov>=readthr && len>=mintranscriptlen)) { // store transcript here
	// //if(id || ( cov>=readthr && len>=mintranscriptlen)) { // store transcript here
		char sign='-';
		if(strand) { sign='+';}
		if(first) { geneno++;}
		//CPrediction *p=new CPrediction(geneno-1, id, exons[0].start, exons.Last().end, cov, sign, fragno, len);
		//if(t) fprintf(stderr,"store prediction with start=%d and end=%d\n",exons[0].start, exons.Last().end);
		float gcov=cov;
		if(t && t->exons.Count()==1) { // if single exon
			RC_TData* tdata=(RC_TData*)(t->uptr);
			if(len) gcov=(tdata->t_exons[0])->movlcount/len;
			if(cov<gcov) gcov=cov;
		}

		/*
		{ // DEBUG ONLY
			if(t) {
				fprintf(stderr,"Coverage for transcript %s = %f\n",t->getID(),gcov);
			}
		}
	    */

		CPrediction *p=new CPrediction(geneno-1, t, exons[0].start, exons.Last().end, gcov, sign, len);
		p->exons=exons;
		if(t && t->exons.Count()==1) exoncov[0]=gcov;
		p->exoncov=exoncov;
		pred.Add(p);
		first=false;

		//fprintf(stderr,"Transcript stored\n");
	}

	return(cov);
}

/*
// I don't use this one
float store_transcript(GList<CPrediction>& pred,GVec<int>& path,GVec<float>& nodeflux,GVec<float>& nodecov,
		GPVec<CGraphnode>& no2gnode,int& geneno,bool& first,int strand, char* id=NULL) {

	float cov=0;
	int len=0;
	CGraphnode *prevnode=NULL;
	GVec<GSeg> exons;

	int s=0;
	if(!path[0]) s=1;

	for(int i=s;i<path.Count()-1;i++) {

		CGraphnode *node=no2gnode[path[i]];
		if(!prevnode || node->start>prevnode->end+1) { // this is a new exon
			GSeg exon(node->start,node->end);
			exons.Add(exon);
		}
		else exons.Last().end=node->end;
		len+=node->end-node->start+1;
		float usedcov=node->cov;
		if(node->cov) {
			if(node->capacity) usedcov*=nodeflux[i]/node->capacity;
			else if(prevnode && prevnode->capacity) usedcov*=nodeflux[i-1]/prevnode->capacity;
			else {
				float rate=find_capacity(i,path,nodeflux,no2gnode);
				usedcov*=rate;
			}
		}
		nodecov[path[i]]-=usedcov/(node->end-node->start+1);
		cov+=usedcov;
		prevnode=node;
	}

	if(len) cov/=len;

	//fprintf(stderr,"Predicted transcript cov=%f\n",cov);

	if(id || (cov>=readthr && len>=mintranscriptlen)) { // store transcript here
		char sign='-';
		if(strand) { sign='+';}
		if(first) { geneno++;}
		CPrediction *p=new CPrediction(geneno-1,id,exons[0].start,exons.Last().end,cov,sign);
		p->exons=exons;
		pred.Add(p);
		first=false;
	}

	return(cov);
}
*/

/*
// I don't use this one
// this is the max flow path that works the best:
float find_max_flow_path(int gno,GVec<int> *path,GBitVec *pathpat,GBitVec *istranscript,GPVec<CTransfrag>& transfrag,
		GPVec<CGraphnode>&no2gnode,GVec<float>& nodeflux) {

	float pathcov=0;

	for(int j=0;j<no2gnode[0]->trf.Count();j++) {
		istranscript[0][no2gnode[0]->trf[j]]=1;
	}

	GVec<int> node2path;
	node2path.Resize(gno,-1);

	GVec<float> nodecapacity;

	GVec<int> nodes;
	GVec<int> nodetmp;
	bool all=false;
	bool alltmp=false;

	pathpat[0][0]=1;
	path[0].cAdd(0);

	for(int i=1;i<gno;i++) {
		CGraphnode *inode=no2gnode[i];
		int maxp=-1;
		float maxcov=0;
		GVec<int> maxrestore;
		if(i<gno-1 && inode->parent.Count()==1) {
			maxp=0;
			CGraphnode *pnode=no2gnode[inode->parent[maxp]];
			pathpat[inode->parent[maxp]][i]=1;
			pathpat[inode->parent[maxp]][edge(inode->parent[maxp],i,gno)]=1;
			for(int j=0;j<pnode->trf.Count();j++) {
				if(istranscript[inode->parent[maxp]][pnode->trf[j]] && transfrag[pnode->trf[j]]->nodes.Last()>inode->parent[maxp] &&
						!onpath(transfrag[pnode->trf[j]]->pattern,transfrag[pnode->trf[j]]->nodes,pathpat[inode->parent[maxp]],0,i,no2gnode,gno)) { // transfrag is not compatible with parent i
					maxrestore.Add(pnode->trf[j]);
				}
			}
			pathpat[inode->parent[maxp]][i]=0;
			pathpat[inode->parent[maxp]][edge(inode->parent[maxp],i,gno)]=0;
			//fprintf(stderr,"%d-%d skip cov\n",inode->parent[maxp],i);
		}
		else for(int p=0;p<inode->parent.Count();p++) {
			if(path[inode->parent[p]].Count()) {
				GVec<int> restore;
				CGraphnode *pnode=no2gnode[inode->parent[p]];
				pathpat[inode->parent[p]][i]=1;
				pathpat[inode->parent[p]][edge(inode->parent[p],i,gno)]=1;
				for(int j=0;j<pnode->trf.Count();j++) {
					if(istranscript[inode->parent[p]][pnode->trf[j]] && transfrag[pnode->trf[j]]->nodes.Last()>inode->parent[p] &&
							!onpath(transfrag[pnode->trf[j]]->pattern,transfrag[pnode->trf[j]]->nodes,pathpat[inode->parent[p]],0,i,no2gnode,gno)) { // transfrag is not compatible with parent i
						istranscript[inode->parent[p]][pnode->trf[j]]=0;
						restore.Add(pnode->trf[j]);
					}
				}
				path[inode->parent[p]].Add(i);
				float cov=max_flow_partial(i,gno,path[inode->parent[p]],istranscript[inode->parent[p]],transfrag,no2gnode,node2path,nodecapacity,nodetmp,alltmp);

				//fprintf(stderr,"%d-%d cov=%f all=%d\n",inode->parent[p],i,cov,alltmp);

				if(cov>maxcov) {
					maxp=p;
					maxcov=cov;
					maxrestore=restore;
					if(i==gno-1) { nodeflux=nodecapacity; nodes=nodetmp;all=alltmp;}
				}
				for(int t=0;t<restore.Count();t++) istranscript[inode->parent[p]][restore[t]]=1;
				path[inode->parent[p]].Pop();
				pathpat[inode->parent[p]][i]=0;
				pathpat[inode->parent[p]][edge(inode->parent[p],i,gno)]=0;
			}
		}

		if(maxp>-1) {
			pathpat[i]=pathpat[inode->parent[maxp]];
			pathpat[i][i]=1;
			pathpat[i][edge(inode->parent[maxp],i,gno)]=1;
			istranscript[i]=istranscript[inode->parent[maxp]];
			for(int t=0;t<maxrestore.Count();t++) {
				istranscript[i][maxrestore[t]]=0;
			}
			path[i]=path[inode->parent[maxp]];
			path[i].Add(i);
			for(int j=0;j<inode->trf.Count();j++) {
				if(transfrag[inode->trf[j]]->nodes[0]==i) {
					istranscript[i][inode->trf[j]]=1;
				}
			}
			if(i==gno-1) {
				pathcov=maxcov;
				// update transfrag coverages
				nodecapacity=nodeflux;
				for(int j=0;j<path[i].Count()-1;j++) {
					if(all || nodes[j]) {
					CGraphnode *jnode=no2gnode[path[i][j]];
					int t=0;
					while(t<jnode->trf.Count() && nodecapacity[j]>epsilon) {
						if(istranscript[i][jnode->trf[t]] && transfrag[jnode->trf[t]]->nodes[0]==path[i][j]) {
							if(nodecapacity[j]>transfrag[jnode->trf[t]]->abundance) {
								nodecapacity[j]-=transfrag[jnode->trf[t]]->abundance;
								transfrag[jnode->trf[t]]->abundance=0;
							}
							else {
								transfrag[jnode->trf[t]]->abundance-=nodecapacity[j];
								if(transfrag[jnode->trf[t]]->abundance<epsilon) transfrag[jnode->trf[t]]->abundance=0;
								nodecapacity[j]=0;
							}
						}
						t++;
					}
				}
				}
			}
		}
	}

	return(pathcov);
}
*/

/*
// I don't use this one
float find_weight_max_flow_path(int gno,GVec<int> *path,GBitVec *pathpat,GBitVec *istranscript,GPVec<CTransfrag>& transfrag,
		GPVec<CGraphnode>&no2gnode,GVec<float>& nodeflux) {

	float pathcov=0;

	for(int j=0;j<no2gnode[0]->trf.Count();j++) {
		istranscript[0][no2gnode[0]->trf[j]]=1;
	}

	GVec<int> node2path;
	node2path.Resize(gno,-1);

	GVec<float> nodecapacity;

	pathpat[0][0]=1;
	path[0].cAdd(0);

	for(int i=1;i<gno;i++) {
		CGraphnode *inode=no2gnode[i];
		int maxp=-1;
		float maxcov=0;
		GVec<int> maxrestore;
		if(i<gno-1 && inode->parent.Count()==1) {
			maxp=0;
			CGraphnode *pnode=no2gnode[inode->parent[maxp]];
			pathpat[inode->parent[maxp]][i]=1;
			pathpat[inode->parent[maxp]][edge(inode->parent[maxp],i,gno)]=1;
			for(int j=0;j<pnode->trf.Count();j++) {
				if(istranscript[inode->parent[maxp]][pnode->trf[j]] && transfrag[pnode->trf[j]]->nodes.Last()>inode->parent[maxp] &&
						!onpath(transfrag[pnode->trf[j]]->pattern,transfrag[pnode->trf[j]]->nodes,pathpat[inode->parent[maxp]],0,i,no2gnode,gno)) { // transfrag is not compatible with parent i
					maxrestore.Add(pnode->trf[j]);
				}
			}
			pathpat[inode->parent[maxp]][i]=0;
			pathpat[inode->parent[maxp]][edge(inode->parent[maxp],i,gno)]=0;
			//fprintf(stderr,"%d-%d skip cov\n",inode->parent[maxp],i);
		}
		else for(int p=0;p<inode->parent.Count();p++) {
			if(path[inode->parent[p]].Count()) {
				GVec<int> restore;
				CGraphnode *pnode=no2gnode[inode->parent[p]];
				pathpat[inode->parent[p]][i]=1;
				pathpat[inode->parent[p]][edge(inode->parent[p],i,gno)]=1;
				for(int j=0;j<pnode->trf.Count();j++) {
					if(istranscript[inode->parent[p]][pnode->trf[j]] && transfrag[pnode->trf[j]]->nodes.Last()>inode->parent[p] &&
							!onpath(transfrag[pnode->trf[j]]->pattern,transfrag[pnode->trf[j]]->nodes,pathpat[inode->parent[p]],0,i,no2gnode,gno)) { // transfrag is not compatible with parent i
						istranscript[inode->parent[p]][pnode->trf[j]]=0;
						restore.Add(pnode->trf[j]);
					}
				}
				path[inode->parent[p]].Add(i);
				float cov=weight_max_flow_partial(i,gno,path[inode->parent[p]],istranscript[inode->parent[p]],transfrag,no2gnode,node2path,nodecapacity);

				//fprintf(stderr,"%d-%d cov=%f all=%d\n",inode->parent[p],i,cov);

				if(cov>maxcov) {
					maxp=p;
					maxcov=cov;
					maxrestore=restore;
					if(i==gno-1) { nodeflux=nodecapacity;}
				}
				for(int t=0;t<restore.Count();t++) istranscript[inode->parent[p]][restore[t]]=1;
				path[inode->parent[p]].Pop();
				pathpat[inode->parent[p]][i]=0;
				pathpat[inode->parent[p]][edge(inode->parent[p],i,gno)]=0;
			}
		}

		if(maxp>-1) {
			pathpat[i]=pathpat[inode->parent[maxp]];
			pathpat[i][i]=1;
			pathpat[i][edge(inode->parent[maxp],i,gno)]=1;
			istranscript[i]=istranscript[inode->parent[maxp]];
			for(int t=0;t<maxrestore.Count();t++) {
				istranscript[i][maxrestore[t]]=0;
			}
			path[i]=path[inode->parent[maxp]];
			path[i].Add(i);
			for(int j=0;j<inode->trf.Count();j++) {
				if(transfrag[inode->trf[j]]->nodes[0]==i) {
					istranscript[i][inode->trf[j]]=1;
				}
			}
			if(i==gno-1) {
				pathcov=maxcov;
				// update transfrag coverages
				nodecapacity=nodeflux;
				for(int j=0;j<path[i].Count()-1;j++) {
					CGraphnode *jnode=no2gnode[path[i][j]];
					int t=0;
					while(t<jnode->trf.Count() && nodecapacity[j]>epsilon) {
						if(istranscript[i][jnode->trf[t]] && transfrag[jnode->trf[t]]->nodes[0]==path[i][j]) {
							if(nodecapacity[j]>transfrag[jnode->trf[t]]->abundance) {
								nodecapacity[j]-=transfrag[jnode->trf[t]]->abundance;
								transfrag[jnode->trf[t]]->abundance=0;
							}
							else {
								transfrag[jnode->trf[t]]->abundance-=nodecapacity[j];
								if(transfrag[jnode->trf[t]]->abundance<epsilon) transfrag[jnode->trf[t]]->abundance=0;
								nodecapacity[j]=0;
							}
						}
						t++;
					}
				}
			}
		}
	}

	return(pathcov);
}
*/

/*
// I don't use this one
float find_max_flow_path_back(int gno,GVec<int> *path,GBitVec *pathpat,GBitVec *istranscript,GPVec<CTransfrag>& transfrag,
		GPVec<CGraphnode>&no2gnode,GVec<float>& nodeflux) {

	float pathcov=0;

	for(int j=0;j<no2gnode[gno-1]->trf.Count();j++) {
		istranscript[gno-1][no2gnode[gno-1]->trf[j]]=1;
	}

	GVec<int> node2path;
	node2path.Resize(gno,-1);

	GVec<float> nodecapacity;

	pathpat[gno-1][gno-1]=1;
	path[gno-1].cAdd(gno-1);

	for(int i=gno-2;i>=0;i--) {
		CGraphnode *inode=no2gnode[i];
		int maxc=-1;
		float maxcov=0;
		GVec<int> maxrestore;
		if(i && inode->child.Count()==1) {
			maxc=0;
			CGraphnode *cnode=no2gnode[inode->child[maxc]];
			pathpat[inode->child[maxc]][i]=1;
			pathpat[inode->child[maxc]][edge(i,inode->child[maxc],gno)]=1;
			for(int j=0;j<cnode->trf.Count();j++) {
				if(istranscript[inode->child[maxc]][cnode->trf[j]] && transfrag[cnode->trf[j]]->nodes[0]<inode->child[maxc] &&
						!onpath(transfrag[cnode->trf[j]]->pattern,transfrag[cnode->trf[j]]->nodes,pathpat[inode->child[maxc]],i,gno-1,no2gnode,gno)) { // transfrag is not compatible with child of i
					maxrestore.Add(cnode->trf[j]);
				}
			}
			pathpat[inode->child[maxc]][i]=0;
			pathpat[inode->child[maxc]][edge(i,inode->child[maxc],gno)]=0;
			//fprintf(stderr,"%d-%d skip cov\n",i,inode->child[maxc]);
		}
		else for(int c=0;c<inode->child.Count();c++) {
			if(path[inode->child[c]].Count()) {
				GVec<int> restore;
				CGraphnode *cnode=no2gnode[inode->child[c]];
				pathpat[inode->child[c]][i]=1;
				pathpat[inode->child[c]][edge(i,inode->child[c],gno)]=1;
				for(int j=0;j<cnode->trf.Count();j++) {
					if(istranscript[inode->child[c]][cnode->trf[j]] && transfrag[cnode->trf[j]]->nodes[0]<inode->child[c] &&
							!onpath(transfrag[cnode->trf[j]]->pattern,transfrag[cnode->trf[j]]->nodes,pathpat[inode->child[c]],i,gno-1,no2gnode,gno)) { // transfrag is not compatible with child of i
						istranscript[inode->child[c]][cnode->trf[j]]=0;
						restore.Add(cnode->trf[j]);
					}
				}
				path[inode->child[c]].Add(i);
				float cov=max_flow_partial_back(i,path[inode->child[c]],istranscript[inode->child[c]],transfrag,no2gnode,node2path,nodecapacity);

				//fprintf(stderr,"%d-%d cov=%f\n",i,inode->child[c],cov);

				if(cov>maxcov) {
					maxc=c;
					maxcov=cov;
					maxrestore=restore;
					if(i==0) nodeflux=nodecapacity;
				}
				for(int t=0;t<restore.Count();t++) istranscript[inode->child[c]][restore[t]]=1;
				path[inode->child[c]].Pop();
				pathpat[inode->child[c]][i]=0;
				pathpat[inode->child[c]][edge(i,inode->child[c],gno)]=0;
			}
		}

		if(maxc>-1) {
			pathpat[i]=pathpat[inode->child[maxc]];
			pathpat[i][i]=1;
			pathpat[i][edge(i,inode->child[maxc],gno)]=1;
			istranscript[i]=istranscript[inode->child[maxc]];
			for(int t=0;t<maxrestore.Count();t++) {
				istranscript[i][maxrestore[t]]=0;
			}
			path[i]=path[inode->child[maxc]];
			path[i].Add(i);
			for(int j=0;j<inode->trf.Count();j++) {
				if(transfrag[inode->trf[j]]->nodes.Last()==i) {
					istranscript[i][inode->trf[j]]=1;
				}
			}
			if(!i) {
				pathcov=maxcov;
				// update transfrag coverages
				nodecapacity=nodeflux;
				int n=path[i].Count()-1;
				for(int j=0;j<n;j++) {
					CGraphnode *jnode=no2gnode[path[i][n-j]];
					int t=0;
					while(t<jnode->trf.Count() && nodecapacity[j]>epsilon) {
						if(istranscript[i][jnode->trf[t]] && transfrag[jnode->trf[t]]->nodes[0]==path[i][n-j]) {
							if(nodecapacity[j]>transfrag[jnode->trf[t]]->abundance) {
								nodecapacity[j]-=transfrag[jnode->trf[t]]->abundance;
								transfrag[jnode->trf[t]]->abundance=0;
							}
							else {
								transfrag[jnode->trf[t]]->abundance-=nodecapacity[j];
								if(transfrag[jnode->trf[t]]->abundance<epsilon) transfrag[jnode->trf[t]]->abundance=0;
								nodecapacity[j]=0;
							}
						}
						t++;
					}
				}
			}
		}
	}

	return(pathcov);
}
*/

/*
// I don't use this one
// this one uses max_compon_size to compute compatible transcripts -> probably overkill
float find_max_flow_path(int gno,GVec<int> *path,GBitVec *pathpat,GBitVec *istranscript,GPVec<CTransfrag>& transfrag,
		GPVec<CGraphnode>&no2gnode,GVec<float>& nodeflux,GVec<bool>& compatible) {

	GHash<CComponent> computed;
	float pathcov=0;

	for(int j=0;j<no2gnode[0]->trf.Count();j++) {
		istranscript[0][no2gnode[0]->trf[j]]=1;
	}

	GVec<int> node2path;
	node2path.Resize(gno,-1);

	GVec<float> nodecapacity;

	pathpat[0][0]=1;
	path[0].cAdd(0);

	for(int i=1;i<gno;i++) {
		CGraphnode *inode=no2gnode[i];
		int maxp=-1;
		float maxcov=0;
		GVec<int> maxrestore;
		if(i<gno-1 && inode->parent.Count()==1) {
			maxp=0;
			CGraphnode *pnode=no2gnode[inode->parent[maxp]];
			pathpat[inode->parent[maxp]][i]=1;
			pathpat[inode->parent[maxp]][edge(inode->parent[maxp],i,gno)]=1;
			for(int j=0;j<pnode->trf.Count();j++) {
				if(istranscript[inode->parent[maxp]][pnode->trf[j]] && transfrag[pnode->trf[j]]->nodes.Last()>inode->parent[maxp] &&
						!onpath(transfrag[pnode->trf[j]]->pattern,transfrag[pnode->trf[j]]->nodes,pathpat[inode->parent[maxp]],0,i,no2gnode,gno)) { // transfrag is not compatible with parent i
					maxrestore.Add(pnode->trf[j]);
				}
			}
			pathpat[inode->parent[maxp]][i]=0;
			pathpat[inode->parent[maxp]][edge(inode->parent[maxp],i,gno)]=0;

			//fprintf(stderr,"%d-%d skip cov\n",inode->parent[maxp],i);
		}
		else for(int p=0;p<inode->parent.Count();p++) {
			if(path[inode->parent[p]].Count()) {
				GVec<CTrInfo> set;
				GVec<int> restore;
				CGraphnode *pnode=no2gnode[inode->parent[p]];
				pathpat[inode->parent[p]][i]=1;
				pathpat[inode->parent[p]][edge(inode->parent[p],i,gno)]=1;
				for(int j=0;j<pnode->trf.Count();j++) {
					if(istranscript[inode->parent[p]][pnode->trf[j]] && transfrag[pnode->trf[j]]->nodes.Last()>inode->parent[p]) {
						if(onpath(transfrag[pnode->trf[j]]->pattern,transfrag[pnode->trf[j]]->nodes,pathpat[inode->parent[p]],0,i,no2gnode,gno)) {
							CTrInfo tr(pnode->trf[j],transfrag[pnode->trf[j]]->abundance,0);
							set.Add(tr);
						}
						else { // transfrag is not compatible with child i
							istranscript[inode->parent[p]][pnode->trf[j]]=0;
							restore.Add(pnode->trf[j]);
						}
					}
				}
				path[inode->parent[p]].Add(i);

				float cov=0;
				// compute the maximum component
				GVec<int> *maxset=NULL;
				if(set.Count()) { // there are transcripts compatible with this path to child i from parent p
					// set.Sort(setCmp); no need to sort here -> only (maybe) if I collected the out transfrag without binary insert
					float maxsize=MIN_VAL;
					maxset=max_compon_size(transfrag.Count(),maxsize,set,compatible,computed);


					{ // DEBUG ONLY
						fprintf(stderr,"Max set[%d]: ",i);
						if(maxset) for(int t=0;t<maxset->Count();t++) fprintf(stderr," %d",maxset->Get(t));
						fprintf(stderr," from set: ");
						for(int t=0;t<set.Count();t++) fprintf(stderr," %d",set[t].trno);
						fprintf(stderr,"\n");
					}


					if(maxset) { // found compatible transcripts -> take out the uncompatible ones
						int j=0;
						int k=0;
						GVec<int> ignore;
						while(k<set.Count()) {
							while(k<set.Count() && ((j<maxset->Count() && set[k].trno<maxset->Get(j)) || (j==maxset->Count()))) {
								if(istranscript[inode->parent[p]][set[k].trno]) { // this transcript got removed
									istranscript[inode->parent[p]][set[k].trno]=0;
									ignore.Add(set[k].trno);
								}
								k++;
							}
							if(j<maxset->Count()) j++;
							k++;
						}
						//cov=max_flow_partial(i,gno,path[inode->parent[p]],istranscript[inode->parent[p]],transfrag,no2gnode,node2path,nodecapacity);
						for(j=0;j<ignore.Count();j++) { // all of these transcripts are compatible with p->i path so I will keep them for now
							istranscript[inode->parent[p]][ignore[j]]=1;
						}
					}
				}

				//fprintf(stderr,"%d-%d cov=%f\n",inode->parent[p],i,cov);

				if(cov>maxcov) {
					maxp=p;
					maxcov=cov;
					maxrestore=restore;
					if(i==gno-1) nodeflux=nodecapacity;
				}
				for(int t=0;t<restore.Count();t++) istranscript[inode->parent[p]][restore[t]]=1;
				path[inode->parent[p]].Pop();
				pathpat[inode->parent[p]][i]=0;
				pathpat[inode->parent[p]][edge(inode->parent[p],i,gno)]=0;
			}
		}

		if(maxp>-1) {
			pathpat[i]=pathpat[inode->parent[maxp]];
			pathpat[i][i]=1;
			pathpat[i][edge(inode->parent[maxp],i,gno)]=1;
			istranscript[i]=istranscript[inode->parent[maxp]];
			for(int t=0;t<maxrestore.Count();t++) {
				istranscript[i][maxrestore[t]]=0;
			}
			path[i]=path[inode->parent[maxp]];
			path[i].Add(i);
			for(int j=0;j<inode->trf.Count();j++) {
				if(transfrag[inode->trf[j]]->nodes[0]==i) {
					istranscript[i][inode->trf[j]]=1;
				}
			}
			if(i==gno-1) {
				pathcov=maxcov;
				// update transfrag coverages
				nodecapacity=nodeflux;
				for(int j=0;j<path[i].Count()-1;j++) {
					CGraphnode *jnode=no2gnode[path[i][j]];
					int t=0;
					while(t<jnode->trf.Count() && nodecapacity[j]>epsilon) {
						if(istranscript[i][jnode->trf[t]] && transfrag[jnode->trf[t]]->nodes[0]==path[i][j]) {
							if(nodecapacity[j]>transfrag[jnode->trf[t]]->abundance) {
								nodecapacity[j]-=transfrag[jnode->trf[t]]->abundance;
								transfrag[jnode->trf[t]]->abundance=0;
							}
							else {
								transfrag[jnode->trf[t]]->abundance-=nodecapacity[j];
								if(transfrag[jnode->trf[t]]->abundance<epsilon) transfrag[jnode->trf[t]]->abundance=0;
								nodecapacity[j]=0;
							}
						}
						t++;
					}
				}
			}
		}
	}

	return(pathcov);
}
*/

/*
// I don't use this one
float find_max_flow(int gno,GVec<int>& path, GBitVec& pathpat,GPVec<CTransfrag>& transfrag,GPVec<CGraphnode>& no2gnode,GVec<float>& nodeflux)
{
	GVec<int> fwdpath[gno];
	GBitVec fwdistranscript[gno];
	GBitVec fwdpathpat[gno];

	for(int i=0;i<gno;i++) {
		fwdistranscript[i].resize(transfrag.Count(),false);
		fwdpathpat[i].resize(1+gno*(gno+1)/2,false);
	}

	float maxcov=find_max_flow_path(gno,fwdpath,fwdpathpat,fwdistranscript,transfrag,no2gnode,nodeflux);
	int maxpath=gno-1;

	if(maxcov) {
		GVec<int> backpath[gno];
		GBitVec backistranscript[gno];
		GBitVec backpathpat[gno];

		GVec<float> nodecapacity;

		for(int i=0;i<gno;i++) {
			backistranscript[i].resize(transfrag.Count(),false);
			backpathpat[i].resize(1+gno*(gno+1)/2,false);
		}
		float cov=find_max_flow_path_back(gno,backpath,backpathpat,backistranscript,transfrag,no2gnode,nodecapacity);

		if(cov>maxcov) {
			maxcov=cov;
			maxpath=0;
			backpath[0].Reverse();
			fwdpath[0]=backpath[0];
			fwdpathpat[0]=backpathpat[0];
			nodeflux=nodecapacity;
		}

		GVec<int> node2path;
		node2path.Resize(gno,-1);

		for(int i=1;i<gno-1;i++) {
			fwdpath[i].Pop();
			backpath[i].Reverse();
			fwdpath[i].Add(backpath[i]);
			fwdpathpat[i] = fwdpathpat[i] | backpathpat[i];
			CGraphnode *inode=no2gnode[i];
			for(int j=0;j<inode->trf.Count();j++) {
				int t=inode->trf[j];
				if(!fwdistranscript[i][t] || !backistranscript[i][t] ){
					fwdistranscript[i][t]=0;
					backistranscript[i][t]=0;
				}
			}
			fwdistranscript[i]=fwdistranscript[i]|backistranscript[i];
			for(int j=0;j<fwdpath[i].Count();j++) node2path[fwdpath[i][j]]=j;
			//cov=max_flow_partial(gno-1,gno,fwdpath[i],fwdistranscript[i],transfrag,no2gnode,node2path,nodecapacity);
			if(cov>maxcov) {
				maxcov=cov;
				maxpath=i;
				nodeflux=nodecapacity;
			}
		}

	}

	if(maxcov) {
		path=fwdpath[maxpath];
		pathpat=fwdpathpat[maxpath];
	}

	return(maxcov);
}
*/

/*
// I don't use this one
void parse_trf_max_flow(int gno,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,
		int& geneno,int strand,GList<CPrediction>& pred,GVec<float>& nodecov,GBitVec& prevpath) {


	GVec<int> path[gno];
	GBitVec istranscript[gno];
	GBitVec pathpat[gno];

	for(int i=0;i<gno;i++) {
		istranscript[i].resize(transfrag.Count(),false);
		pathpat[i].resize(1+gno*(gno+1)/2,false);
	}

	bool first=true;
	GVec<float> nodeflux;

	float cov;
	do {
		cov=find_max_flow_path(gno,path,pathpat,istranscript,transfrag,no2gnode,nodeflux);


		// Some previous tries:
		// 1: cov=find_max_flow_path_back(gno,path,pathpat,istranscript,transfrag,no2gnode,nodeflux);
		// 2: cov=find_max_flow(gno,path,pathpat,transfrag,no2gnode,nodeflux);
		// 3: cov=find_max_flow_path(gno,path,pathpat,istranscript,transfrag,no2gnode,nodeflux,compatible);



		{ // DEBUG ONLY
			fprintf(stderr,"cov=%f\n",cov);
			fprintf(stderr,"path=");
			for(int j=0;j<path[gno-1].Count();j++) fprintf(stderr," %d",path[gno-1][j]);
			//for(int j=path[0].Count()-1;j>=0;j--) fprintf(stderr," %d",path[0][j]);
			fprintf(stderr,"\n");
		}


		if(cov>=readthr) {
			bool included=true;
			store_transcript(pred,path[gno-1],nodeflux,nodecov,no2gnode,geneno,first,strand,gno,included,prevpath);
			//path[0].Reverse();store_transcript(pred,path[0],nodeflux,nodecov,no2gnode,geneno,first,strand,gno,included,prevpath);
			//store_transcript(pred,path[gno-1],nodeflux,nodecov,no2gnode,geneno,first,strand);
			//store_transcript(pred,path,nodeflux,nodecov,no2gnode,geneno,first,strand,gno,included,prevpath);


			{ // DEBUG ONLY
				fprintf(stderr,"\nAfter update:\n");
				for(int i=1;i<gno;i++) {
					fprintf(stderr,"Node %d: %f ",i,nodecov[i]);
					fprintf(stderr,"trf=");
					for(int t=0;t<no2gnode[i]->trf.Count();t++) fprintf(stderr," %d(%f)",no2gnode[i]->trf[t],transfrag[no2gnode[i]->trf[t]]->abundance);
					fprintf(stderr,"\n");
				}
			}


			int maxi=0;
			float maxcov=0;
			for(int i=0;i<gno;i++) {
				path[i].Clear();
				istranscript[i].reset();
				pathpat[i].reset();
				if(nodecov[i]>maxcov) {
					maxcov=nodecov[i];
					maxi=i;
				}
				cov=nodecov[maxi];
			}
		}
		else if(cov) { // try another path by eliminating maximum abundance edge
			float maxabund=0;
			int maxt=-1;
			for(int t=0;t<transfrag.Count();t++)
				if(transfrag[t]->abundance>maxabund) {
					maxabund=transfrag[t]->abundance;
					maxt=t;
				}
			if(maxt>-1) {
				transfrag[maxt]->abundance=0;
				int maxi=0;
				float maxcov=0;
				for(int i=0;i<gno;i++) {
					path[i].Clear();
					istranscript[i].reset();
					pathpat[i].reset();
					if(nodecov[i]>maxcov) {
						maxcov=nodecov[i];
						maxi=i;
					}
					cov=nodecov[maxi];
				}
			}
		}

	} while (cov>=readthr);
}
*/

/*
// I don't use this one
void parse_trf_weight_max_flow(int gno,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,
		int& geneno,int strand,GList<CPrediction>& pred,GVec<float>& nodecov,GBitVec& prevpath) {


	GVec<int> path[gno];
	GBitVec istranscript[gno];
	GBitVec pathpat[gno];

	for(int i=0;i<gno;i++) {
		istranscript[i].resize(transfrag.Count(),false);
		pathpat[i].resize(1+gno*(gno+1)/2,false);
	}

	bool first=true;
	GVec<float> nodeflux;

	float cov;
	do {
		cov=find_weight_max_flow_path(gno,path,pathpat,istranscript,transfrag,no2gnode,nodeflux);


		{ // DEBUG ONLY
			fprintf(stderr,"cov=%f\n",cov);
			fprintf(stderr,"path=");
			for(int j=0;j<path[gno-1].Count();j++) fprintf(stderr," %d",path[gno-1][j]);
			//for(int j=path[0].Count()-1;j>=0;j--) fprintf(stderr," %d",path[0][j]);
			fprintf(stderr,"\n");
		}


		if(cov>=readthr) {
			bool included=true;
			store_transcript(pred,path[gno-1],nodeflux,nodecov,no2gnode,geneno,first,strand,gno,included,prevpath);


			{ // DEBUG ONLY
				fprintf(stderr,"\nAfter update:\n");
				for(int i=1;i<gno;i++) {
					fprintf(stderr,"Node %d: %f ",i,nodecov[i]);
					fprintf(stderr,"trf=");
					for(int t=0;t<no2gnode[i]->trf.Count();t++) fprintf(stderr," %d(%f)",no2gnode[i]->trf[t],transfrag[no2gnode[i]->trf[t]]->abundance);
					fprintf(stderr,"\n");
				}
			}


			int maxi=0;
			float maxcov=0;
			for(int i=0;i<gno;i++) {
				path[i].Clear();
				istranscript[i].reset();
				pathpat[i].reset();
				if(nodecov[i]>maxcov) {
					maxcov=nodecov[i];
					maxi=i;
				}
				cov=nodecov[maxi];
			}
		}
		else if(cov) { // try another path by eliminating maximum abundance edge
			float maxabund=0;
			int maxt=-1;
			for(int t=0;t<transfrag.Count();t++)
				if(transfrag[t]->abundance>maxabund) {
					maxabund=transfrag[t]->abundance;
					maxt=t;
				}
			if(maxt>-1) {
				transfrag[maxt]->abundance=0;
				int maxi=0;
				float maxcov=0;
				for(int i=0;i<gno;i++) {
					path[i].Clear();
					istranscript[i].reset();
					pathpat[i].reset();
					if(nodecov[i]>maxcov) {
						maxcov=nodecov[i];
						maxi=i;
					}
					cov=nodecov[maxi];
				}
			}
		}

	} while (cov>=readthr);
}
*/


void parse_trf(int maxi,int gno,int edgeno, GIntHash<int> &gpos,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,
		GBitVec& compatible,	int& geneno,bool first,int strand,GList<CPrediction>& pred,GVec<float>& nodecov,
		GBitVec& istranscript,GBitVec& removable,GBitVec& usednode,float maxcov,GBitVec& prevpath,bool fast) {

	 GVec<int> path;
	 GVec<float> pathincov;
	 GVec<float> pathoutcov;
	 path.Add(maxi);
	 pathincov.cAdd(0.0);
	 pathoutcov.cAdd(0.0);
	 GBitVec pathpat(gno+edgeno);
	 pathpat[maxi]=1;
	 istranscript.reset();
	 GHash<CComponent> computed;

	 float flux=0;
	 //float fragno=0;
	 GVec<float> nodeflux;

	 /*
	 { // DEBUG ONLY
	 	 fprintf(stderr,"start parse_trf with maxi=%d\n",maxi);
		 //fprintf(stderr,"Transcripts before path:");
		 //for(int i=0;i<transfrag.Count();i++) if(istranscript[i]) fprintf(stderr," %d",i);
		 //fprintf(stderr,"\n");

#ifdef GMEMTRACE
	 	 double vm,rsm;
	 	 get_mem_usage(vm, rsm);
	 	 GMessage("\t\tM(s):parse_trf memory usage: rsm=%6.1fMB vm=%6.1fMB\n",rsm/1024,vm/1024);
#endif
	 }
	 */

	 //if(back_to_source_path(maxi,path,pathpat,pathincov,pathoutcov,istranscript,removable,transfrag,computed,compatible,no2gnode,nodecov,gno)) {
	 if((fast && back_to_source_fast(maxi,path,pathpat,transfrag,no2gnode,nodecov,gno,gpos)) ||
			 (!fast && back_to_source_path(maxi,path,pathpat,pathincov,pathoutcov,istranscript,removable,transfrag,computed,compatible,no2gnode,nodecov,gno,gpos))) {
		 	 if(includesource) path.cAdd(0);
	 		 path.Reverse(); // back to source adds the nodes at the end to avoid pushing the list all the time
	 		 pathincov.Reverse();
	 		 pathoutcov.Reverse();

	 		 //if(fwd_to_sink_path(maxi,path,pathpat,pathincov,pathoutcov,istranscript,removable,transfrag,computed,compatible,no2gnode,nodecov,gno)) {
	 		if((fast && fwd_to_sink_fast(maxi,path,pathpat,transfrag,no2gnode,nodecov,gno,gpos)) ||
	 				(!fast && fwd_to_sink_path(maxi,path,pathpat,pathincov,pathoutcov,istranscript,removable,transfrag,computed,compatible,no2gnode,nodecov,gno,gpos))) {
	 			 pathincov.Clear();
	 			 pathoutcov.Clear();

	 			 //removable.reset();
	 			 //flux=update_flux(gno,path,istranscript,transfrag,removable,no2gnode,nodeflux,pathpat);


	 			 if(EM) flux=max_flow_EM(gno,path,istranscript,transfrag,no2gnode,nodeflux,pathpat);
	 			 else if(weight)
	 				 	 //flux=weight_max_flow_EM(gno,path,istranscript,transfrag,no2gnode,nodeflux,pathpat);
	 				 	 flux=weight_max_flow(gno,path,istranscript,transfrag,no2gnode,nodeflux,pathpat);
	 			 else flux=max_flow(gno,path,istranscript,transfrag,no2gnode,nodeflux,pathpat);

	 			 /*
	 			 { // DEBUG ONLY
	 				 printTime(stderr);
	 				 fprintf(stderr,"flux=%g Path:",flux);
	 				 for(int i=0;i<path.Count();i++) fprintf(stderr," %d",path[i]);
	 				 fprintf(stderr,"\n");
	 			 }
	 			 */
	 		}
	 		else {
	 			//pathpat.reset();
	 			//pathpat[maxi]=1;
	 			pathincov.Clear();
	 			pathoutcov.Clear();
	 		}
	 }
	 else {
		 //pathpat.reset();
		 //pathpat[maxi]=1;
		 pathincov.Clear();
		 pathoutcov.Clear();
	 }

	 bool cont=true;

	 if(flux>epsilon) {
		 bool included=true;
		 float cov=store_transcript(pred,path,nodeflux,nodecov,no2gnode,geneno,first,strand,gno,gpos,included,prevpath);

		 /*
		 { // DEBUG ONLY
			 //fprintf(stderr,"Prevpath=");
			 //printBitVec(prevpath);
			 //fprintf(stderr,"\n");
		 	 fprintf(stderr,"cov=%f maxcov=%f\n",cov,maxcov);
		 }
	     */

		 if(cov<isofrac*maxcov) {
			 if(sensitivitylevel) usednode[maxi]=1;
			 else usednode = usednode | prevpath;
			 maxi=0;
			 maxcov=0;
			 cont=false;
		 }
		 else if(cov>maxcov) maxcov=cov;
	 }
	 else {
		 if(sensitivitylevel) usednode[maxi]=1; // start at different locations in graph
		 else {
			 usednode = usednode | prevpath;
			 usednode = usednode | pathpat;
		 }

		 maxi=0;
		 maxcov=0;
		 cont=false;
	 }

	 // Node coverages:
	 for(int i=1;i<gno;i++)
		 if(!usednode[i] && nodecov[i]>nodecov[maxi]) maxi=i;

	 //fprintf(stderr," maxi=%d maxcov=%f\n",maxi,nodecov[maxi]);

	 if(nodecov[maxi]>=readthr && (!specific || cont)) { // if I still have nodes that are above coverage threshold

		 /*
		 { // DEBUG ONLY
			 printTime(stderr);
			 fprintf(stderr,"\nAfter update:\n");
			 for(int i=0;i<gno;i++) {
				 fprintf(stderr,"Node %d: %f ",i,nodecov[i]);
				 fprintf(stderr,"trf=");
				 for(int t=0;t<no2gnode[i]->trf.Count();t++) fprintf(stderr," %d(%f)",no2gnode[i]->trf[t],transfrag[no2gnode[i]->trf[t]]->abundance);
				 fprintf(stderr," maxi=%d maxcov=%f\n",maxi,nodecov[maxi]);
			 }
		 }
		 */

		 path.Clear();
		 nodeflux.Clear();
		 computed.Clear();
		 parse_trf(maxi,gno,edgeno,gpos,no2gnode,transfrag,compatible,geneno,first,strand,pred,nodecov,istranscript,removable,usednode,maxcov,prevpath,fast);
	 }

}


CTransfrag *find_guide_pat(GffObj *guide,GPVec<CGraphnode>& no2gnode,int gno,int edgeno,GIntHash<int> &gpos) {

	CTransfrag *trguide=NULL;

	int i=1;
	while(i<gno-1) {
		if(no2gnode[i]->overlap(guide->exons[0])) { // first exon of guide overlaps node i
			int j=i+1;
			while(j<gno-1 && no2gnode[j]->overlap(guide->exons[0]) && no2gnode[j]->start==no2gnode[j-1]->end+1
					&& no2gnode[j-1]->child.Count() && no2gnode[j-1]->child[0]==j) j++;
			if(j<gno-1 && no2gnode[j]->overlap(guide->exons[0])) i=j-1; // I am skipping nodes that are linked by introns to the last node that is overlapped by the first exon in guide
			else break;
		}
		i++;
	}

	if(i<gno-1) { // found start node
		GBitVec guidepat(gno+edgeno);
		GVec<int> nodes;
		guidepat[i]=1;
		nodes.Add(i);
		CGraphnode *inode=no2gnode[i];
		int j=0; // guide coordinate I am looking at
		int n=guide->exons.Count();
		while(j<n) {
			if(guide->exons[j]->end<inode->end) {
				if(j==n-1) trguide=new CTransfrag(nodes,guidepat); // only find guide abundance later if necessary
				return(trguide);
			}
			else if(guide->exons[j]->end==inode->end) {
				if(j==n-1) {
					trguide=new CTransfrag(nodes,guidepat);
					return(trguide);
				}
				j++;
				int k=0;
				int nc=inode->child.Count();
				while(k<nc && no2gnode[inode->child[k]]->start != guide->exons[j]->start) k++;
				if(k==nc) return(trguide);
				guidepat[inode->child[k]]=1;
				int *pos=gpos[edge(i,inode->child[k],gno)];
				if(pos) guidepat[*pos]=1;
				else GError("Found parent-child %d-%d not linked by edge\n",i,inode->child[k]);
				i=inode->child[k];
				nodes.Add(i);
				inode=no2gnode[i];
			}
			else {
				if(i<gno-1 && (no2gnode[i+1]->start==inode->end+1)) { // !!!  maybe here I need to consider bundledist
					if(inode->child.Count() && inode->child[0]==i+1) {
						guidepat[i+1]=1;
						int *pos=gpos[edge(i,i+1,gno)];
						if(pos) guidepat[*pos]=1;
						else GError("Found parent-child %d-%d not linked by edge\n",i,i+1);
						i++;
						nodes.Add(i);
						inode=no2gnode[i];
					}
					else return(trguide);
				}
				else { // node doesn't continue
					if(j==n-1) trguide=new CTransfrag(nodes,guidepat); // only find guide abundance later if necessary
					return(trguide);
				}
			}
		}
	}
	return(trguide);
}


/*
// I don't use this one
int process_guides(int gno,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,GVec<bool>& compatible,int& geneno,int s,
		GPVec<GffObj>& guides,GList<CPrediction>& pred,GVec<float>& nodecov,GBitVec& istranscript,GBitVec& removable,GBitVec& pathpat) {

	int maxi=1;
	bool cov=false; // tells me if max node coverage was determined

	char strand='-';
	if(s) strand='+';

	// find guides' patterns
	GVec<CTransfrag> guidetrf;
	for(int g=0;g<guides.Count();g++) {
		if((guides[g]->strand==strand) && (guides[g]->overlap(no2gnode[1]->start,no2gnode[gno-2]->end))) {
			CTransfrag *trguide=find_guide_pat(guides[g],no2gnode,gno);
			if(trguide) { // the guide can be found among the graph nodes
				guidetrf.Add(trguide);



				{ // DEBUG ONLY
					fprintf(stderr,"s=%d strand = %c trguide[%d]=",s,strand,g);
					printBitVec(trguide->pattern);
					fprintf(stderr,"\n");
				}

			}
		}
	}

	// compute guides' abundances
	for(int t=0;t<transfrag.Count();t++)
		for(int g=0;g<guidetrf.Count();g++) {
			if(((transfrag[t]->pattern) & guidetrf[g].pattern) == transfrag[t]->pattern) {
				guidetrf[g].abundance+=transfrag[t]->abundance;
			}
		}


	{ // DEBUG ONLY
		for(int g=0;g<guidetrf.Count();g++) {
			fprintf(stderr,"Abundance of guide[%d]=%f with nodes:",g,guidetrf[g].abundance);
			for(int i=0;i<guidetrf[g].nodes.Count();i++) fprintf(stderr," %d",guidetrf[g].nodes[i]);
			fprintf(stderr,"\n");
		}
	}


	guidetrf.Sort(guideabundCmp);

	int g=0;
	bool first=true;
	while(g<guidetrf.Count()) {
		// check if guide is included in previous guide paths
		int p=0;
		bool included=false;
		if(!complete) { // if guides are incomplete exclude the ones that are included into the more complete ones
			while(p<g) {
				//CTransfrag guideg=guidetrf[g];
				//CTransfrag guidep=guidetrf[p];
				if((guidetrf[g].pattern & guidetrf[p].pattern)==guidetrf[g].pattern) {
					included=true;
					guidetrf.Delete(g);
					break;
				}
				p++;
			}
		}

		if(!included) {
			// build guidepath
			istranscript.reset();
			GHash<CComponent> computed;
			GVec<int> path;
			GVec<float> pathincov;
			GVec<float> pathoutcov;
			int n=guidetrf[g].nodes.Count();
			for(int i=n-1;i>0;i--) {
				path.Add(guidetrf[g].nodes[i]);
				pathincov.cAdd(0.0);
				pathoutcov.cAdd(0.0);
				CGraphnode *inode=no2gnode[guidetrf[g].nodes[i]];
				//collect all in transfrags that are compatible with path
				for(int j=0;j<inode->trf.Count();j++) {
					int t=inode->trf[j];
					if(!transfrag[t]->abundance) { // this transfrag was used before -> needs to be deleted
						inode->trf.Delete(j);
						j--;
					}
					else if(transfrag[t]->nodes[0]!=guidetrf[g].nodes[i]) { // this is an in transfrag; I might change this to include the cases where I have a positive in for the transfrag, e.g. if(transfrag[inode->trf[j]]->nodes[0]!=i || transfrag[inode->trf[j]]->in)
						if(istranscript[t] || onpath(transfrag[t]->pattern,transfrag[t]->nodes,guidetrf[g].pattern,guidetrf[g].nodes[0],guidetrf[g].nodes[n-1],no2gnode,gno)) {
							istranscript[t]=1;
							pathincov[path.Count()-1]+=transfrag[t]->abundance;
						}
					}
					else if(transfrag[t]->nodes.Last()!=guidetrf[g].nodes[i]) { // this is an out transfrag; I might change this to include the cases where I have a positive out for the transfrag, e.g. if(transfrag[inode->trf[j]]->nodes.Last()!=i || transfrag[inode->trf[j]]->out)
						if(istranscript[t] || onpath(transfrag[t]->pattern,transfrag[t]->nodes,guidetrf[g].pattern,guidetrf[g].nodes[0],guidetrf[g].nodes[n-1],no2gnode,gno)) {
							istranscript[t]=1;
							pathoutcov[path.Count()-1]+=transfrag[t]->abundance;
						}
					}
				}
			}
			path.Add(guidetrf[g].nodes[0]);
			pathincov.cAdd(0.0);
			pathoutcov.cAdd(0.0);
			if(n>1) {
				CGraphnode *inode=no2gnode[guidetrf[g].nodes[0]];
				//collect all in transfrags that are compatible with path
				for(int j=0;j<inode->trf.Count();j++) {
					int t=inode->trf[j];
					if(!transfrag[t]->abundance) { // this transfrag was used before -> needs to be deleted
						inode->trf.Delete(j);
						j--;
					}
					else if(transfrag[t]->nodes.Last()!=guidetrf[g].nodes[0]) { // this is an out transfrag; I might change this to include the cases where I have a positive out for the transfrag, e.g. if(transfrag[inode->trf[j]]->nodes.Last()!=i || transfrag[inode->trf[j]]->out)
						if(istranscript[t] || onpath(transfrag[t]->pattern,transfrag[t]->nodes,guidetrf[g].pattern,guidetrf[g].nodes[0],guidetrf[g].nodes[n-1],no2gnode,gno)) {
							istranscript[t]=1;
							pathoutcov[path.Count()-1]+=transfrag[t]->abundance;
						}
					}
				}
			}

			if(back_to_source_path(guidetrf[g].nodes[0],path,guidetrf[g].pattern,pathincov,pathoutcov,istranscript,removable,transfrag,computed,compatible,no2gnode,nodecov,gno)) {
					 path.Reverse(); // back to source adds the nodes at the end to avoid pushing the list all the time
					 pathincov.Reverse();
					 pathoutcov.Reverse();

					 if(fwd_to_sink_path(guidetrf[g].nodes.Last(),path,guidetrf[g].pattern,pathincov,pathoutcov,istranscript,removable,transfrag,computed,compatible,no2gnode,nodecov,gno)) {
						 pathincov.Clear();
						 pathoutcov.Clear();

						 GVec<float> nodeflux;
						 removable.reset();
						 //float flux = update_flux(gno,path,istranscript,transfrag,removable,no2gnode,nodeflux,guidetrf[g].pattern);

						 //float flux = max_flow(gno,path,istranscript,transfrag,no2gnode,nodeflux,guidetrf[g].pattern);
						 float flux = weight_max_flow(gno,path,istranscript,transfrag,no2gnode,nodeflux,guidetrf[g].pattern);


			 			 { // DEBUG ONLY
			 				 fprintf(stderr,"flux[%d]=%g Path:",g,flux);
			 				 for(int i=0;i<path.Count();i++) fprintf(stderr," %d",path[i]);
			 				 fprintf(stderr,"\n");
			 			 }


						 if(flux>epsilon) {
							 bool include=true;
							 store_transcript(pred,path,nodeflux,nodecov,no2gnode,geneno,first,s,gno,include,pathpat);
							 cov=true;
							 // Node coverages:
							 for(int i=1;i<gno-1;i++)
								 if(nodecov[i]>nodecov[maxi]) maxi=i;

							 if(nodecov[maxi]<readthr) break; // no need to find other paths since they are under allowed read threshold


							 { // DEBUG ONLY
								 fprintf(stderr,"\nAfter update:\n");
								 for(int i=1;i<gno;i++) {
									 fprintf(stderr,"Node %d: %f ",i,nodecov[i]);
									 fprintf(stderr,"trf=");
									 for(int t=0;t<no2gnode[i]->trf.Count();t++) fprintf(stderr," %d(%f)",no2gnode[i]->trf[t],transfrag[no2gnode[i]->trf[t]]->abundance);
									 fprintf(stderr," maxi=%d maxcov=%f\n",maxi,nodecov[maxi]);
								 }
							 }


							 // clean up what's possible
							 path.Clear();
							 nodeflux.Clear();
							 computed.Clear();
						 }
					 } // end if(fwd_to_sink_path)
				 } // end if(back_to_source_path)
			g++;
		} // end if(!included)
	} // end while(g<guidetrf.Count()

	if(!cov) for(int i=2;i<gno-1;i++)
		if(nodecov[i]>nodecov[maxi]) maxi=i;

	return(maxi);
}
*/

void process_refguides(int gno,int edgeno,GIntHash<int>& gpos,int& lastgpos,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,int s,GPVec<GffObj>& guides,
		GVec<CGuide>& guidetrf) {

	char strand='-';
	if(s) strand='+';

	// find guides' patterns

	for(int g=0;g<guides.Count();g++) {
		//fprintf(stderr,"Consider guide[%d out of %d] %s\n",g,guides.Count(),guides[g]->getID());
		if((guides[g]->strand==strand) && (guides[g]->overlap(no2gnode[1]->start,no2gnode[gno-2]->end))) {
			CTransfrag *trguide=find_guide_pat(guides[g],no2gnode,gno,edgeno,gpos);
			if(trguide) { // the guide can be found among the graph nodes
				//CGuide newguide(trguide,guides[g]->getID());
				CGuide newguide(trguide,guides[g]);
				guidetrf.Add(newguide);

				//fprintf(stderr,"Added guidetrf with ID=%s overlapping transcript interval %d - %d\n",guides[g]->getID(),no2gnode[1]->start,no2gnode[gno-2]->end);

				/*
				{ // DEBUG ONLY
					fprintf(stderr,"s=%d strand = %c trguide[%d]=",s,strand,g);
					//printBitVec(trguide->pattern);
					fprintf(stderr,"\n");
				}
				*/
			}
		}
	}

	// compute guides' abundances
	for(int t=0;t<transfrag.Count();t++) {
		//if(transfrag[t]->nodes.Count()> 1)
		for(int g=0;g<guidetrf.Count();g++) {
			if(((transfrag[t]->pattern) & guidetrf[g].trf->pattern) == transfrag[t]->pattern) {
				guidetrf[g].trf->abundance+=transfrag[t]->abundance;
			}
		}
	}

	/*
	{ // DEBUG ONLY
		for(int g=0;g<guidetrf.Count();g++) {
			fprintf(stderr,"Abundance of guide[%d]=%f with nodes:",g,guidetrf[g].trf->abundance);
			for(int i=0;i<guidetrf[g].trf->nodes.Count();i++) fprintf(stderr," %d",guidetrf[g].trf->nodes[i]);
			fprintf(stderr,"\n");
		}
	}
	*/

	guidetrf.Sort(guidedabundCmp);

	int g=0;
	while(g<guidetrf.Count()) {
		// check if guide is included in previous guide paths
		int p=0;

		//fprintf(stderr,"Process guide=%s\n",guidetrf[g].t->getID());

		while(p<g) {
			//CTransfrag guideg=guidetrf[g];
			//CTransfrag guidep=guidetrf[p];
			if((guidetrf[g].trf->pattern & guidetrf[p].trf->pattern)==guidetrf[g].trf->pattern) {
				guidetrf[g].trf->real=false;  // this marks a guide that is included in another one
				if(!complete) { // if guides are incomplete exclude the ones that are included into the more complete ones
					GFREE(guidetrf[g].trf);
					guidetrf.Delete(g);
					break;
				}
			}
			p++;
		}

		// find if first node of guidetrf extends to source
		bool sourcestart=true;
		int lasti=guidetrf[g].trf->nodes[0]; // first node of guide
		CGraphnode *inode=no2gnode[lasti];
		GVec<int> extendpath;
		float maxabund=trthr; // maximum abundance of source to node transfrag
		int maxnode=lasti; // node where to add the source to
		while(inode->parent.Count()) { // node still has parents
			int i=-1;
			for(int p=0;p<inode->parent.Count();p++) {
				if(!inode->parent[p]) { i=0; break;} // found link to source -> break from path
				if(no2gnode[inode->parent[p]]->end==inode->start-1) { // no splice site between parent and inode
					i=inode->parent[p];
					float abundin=0;
					for(int j=0;j<inode->trf.Count();j++){
						int t=inode->trf[j];
						if(transfrag[t]->nodes[0]<i) { // transfrag doesn't start at this node (in or through transfrag)
							abundin+=transfrag[t]->abundance;
						}
					}

					float nodecov=inode->cov/inode->len();
					if(nodecov-abundin>maxabund) {
						maxabund=nodecov-abundin;
						maxnode=lasti;
					}
				}
			}
			if(i==-1) { sourcestart=false; break;} // did not find a continuing node to source
			extendpath.Add(i);
			inode=no2gnode[i];
			lasti=i;
		}
		int j=guidetrf[g].trf->nodes[0];
		if(j) { // if source is not already among nodes in guide pattern ( it shouldn't be!!!)
			if(sourcestart) { // I can extend path
				for(int i=0;i<extendpath.Count();i++) {
					guidetrf[g].trf->pattern[extendpath[i]]=1;
					int *pos=gpos[edge(extendpath[i],j,gno)];
					if(pos) guidetrf[g].trf->pattern[*pos]=1;
					else GError("Found parent-child %d-%d not linked by edge\n",extendpath[i],j);
					j=extendpath[i];
				}
				if(!includesource) extendpath.Pop(); // last added node is source -> I need to pop it out for weight flow
				extendpath.Reverse();
				extendpath.Add(guidetrf[g].trf->nodes);
				guidetrf[g].trf->nodes.Clear();
				guidetrf[g].trf->nodes.Add(extendpath);
			}
			else { // if no path to source was found I need to add one
				int source=0;
				int *pos=gpos[edge(0,maxnode,gno)];
				if(!pos) { // here is where I need lastgpos
					int key=edge(0,maxnode,gno);
					gpos.Add(key,lastgpos);
					lastgpos++;
					pos=gpos[key];
				}

				// extend guide to maxnode
				if(maxnode<j) {
					int i=0;
					GVec<int> tmpextend;
					while(i<extendpath.Count() && maxnode<=extendpath[i]) {
						tmpextend.Add(extendpath[i]);
						guidetrf[g].trf->pattern[extendpath[i]]=1;
						int *pos1=gpos[edge(extendpath[i],j,gno)];
						if(pos1) guidetrf[g].trf->pattern[*pos1]=1;
						else GError("Found parent-child %d-%d not linked by edge\n",extendpath[i],j);
						j=extendpath[i];
						i++;
					}
					if(includesource) tmpextend.cAdd(0); // I need to comment this if I need path to include the source
					tmpextend.Reverse();
					tmpextend.Add(guidetrf[g].trf->nodes);
					guidetrf[g].trf->nodes.Clear();
					guidetrf[g].trf->nodes.Add(tmpextend);
					guidetrf[g].trf->pattern[0]=1;
					guidetrf[g].trf->pattern[*pos]=1;
				}
				else {
					if(includesource) guidetrf[g].trf->nodes.Insert(0,source); // I need to comment this if I need path to include the source
					guidetrf[g].trf->pattern[0]=1;
					guidetrf[g].trf->pattern[*pos]=1;
				}

				// add source to maxnode transfrag
				GVec<int> nodes;
				nodes.cAdd(0);
				nodes.Add(maxnode);
				GBitVec trpat(gno+edgeno);
				trpat[0]=1;
				trpat[maxnode]=1;
				trpat[*pos]=1;
				CTransfrag *tr=new CTransfrag(nodes,trpat,maxabund);
				transfrag.Add(tr);

				// add source among maxnode parents
				inode=no2gnode[maxnode];
				inode->parent.Insert(0,source);
				inode->parentpat[*pos]=1; // source should be already among parents of maxnode but not the edge to source
			}
			extendpath.Clear();
		}
		// find if I can extend guidetrf to sink
		lasti=guidetrf[g].trf->nodes.Last();

		if(lasti != gno-1) { // sink is not already present in guide pattern
			bool sinkend=true;
			inode=no2gnode[lasti];
			maxabund=trthr;
			maxnode=lasti;
			while(inode->child.Count()) { // node still has children
				int i=-1;
				for(int c=0;c<inode->child.Count();c++) {
					if(inode->child[c]==gno-1) { i=gno-1; break;} // found link to source
					if(no2gnode[inode->child[c]]->start==inode->end+1) { // no splice site
						i=inode->child[c];
						float abundout=0;
						for(int j=0;j<inode->trf.Count();j++){
							int t=inode->trf[j];
							if(transfrag[t]->nodes.Last()>i) { // transfrag doesn't end at this node (out or through transfrag)
								abundout+=transfrag[t]->abundance;
							}
						}

						float nodecov=inode->cov/inode->len();
						if(nodecov-abundout>maxabund) {
							maxabund=nodecov-abundout;
							maxnode=lasti;
						}
					}
				}
				//fprintf(stderr," i=%d maxabund=%g maxnode=%d\n",i,maxabund,maxnode);
				if(i==-1) { sinkend=false; break;} // did not find a continuing node to source
				extendpath.Add(i);
				inode=no2gnode[i];
				lasti=i;
			}
			j=guidetrf[g].trf->nodes.Last();
			if(sinkend) { // extend path
				for(int i=0;i<extendpath.Count();i++) {
					guidetrf[g].trf->nodes.Add(extendpath[i]);
					guidetrf[g].trf->pattern[extendpath[i]]=1;
					int *pos=gpos[edge(j,extendpath[i],gno)];
					if(pos) guidetrf[g].trf->pattern[*pos]=1;
					else GError("Found parent-child %d-%d not linked by edge\n",j,extendpath[i]);
					j=extendpath[i];
				}
			}
			else { // if no path to sink was found I need to add one

				//fprintf(stderr,"no path to sink was found when maxnode=%d j=%d maxabund=%g\n",maxnode,j,maxabund);

				int sink=gno-1;

				// extend guide to maxnode
				if(maxnode>j) {
					int i=0;
					while(i<extendpath.Count() && maxnode>=extendpath[i]) {
						guidetrf[g].trf->nodes.Add(extendpath[i]);
						guidetrf[g].trf->pattern[extendpath[i]]=1;
						int *pos=gpos[edge(j,extendpath[i],gno)];
						if(pos) guidetrf[g].trf->pattern[*pos]=1;
						else GError("Found parent-child %d-%d not linked by edge\n",j,extendpath[i]);
						j=extendpath[i];
						i++;
					}
				}
				guidetrf[g].trf->nodes.Add(sink);
				guidetrf[g].trf->pattern[sink]=1;
				int *pos=gpos[edge(maxnode,sink,gno)];
				if(pos) guidetrf[g].trf->pattern[*pos]=1;
				else {
					int key=edge(maxnode,sink,gno);
					gpos.Add(key,lastgpos);
					lastgpos++;
					pos=gpos[key];
				}

				// add maxnode to sink transfrag
				GVec<int> nodes;
				nodes.Add(maxnode);
				nodes.Add(sink);
				GBitVec trpat(gno+edgeno);
				trpat[maxnode]=1;
				trpat[sink]=1;
				trpat[*pos]=1;
				//fprintf(stderr,"introduce node from %d to sink=%d wih abundance=%g\n",maxnode,sink,maxabund);
				CTransfrag *tr=new CTransfrag(nodes,trpat,maxabund);
				transfrag.Add(tr);

				// add sink among maxnode children
				inode=no2gnode[maxnode];
				inode->child.Add(sink);
				inode->childpat[*pos]=1; // source should be already among parents of maxnode but not the edge to source
			}
		}

		g++;
	}

}

int guides_flow(int gno,GIntHash<int>& gpos,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,GVec<CGuide>& guidetrf,int& geneno,
		int s,GList<CPrediction>& pred,GVec<float>& nodecov,GBitVec& istranscript,GBitVec& pathpat) {

	int maxi=1;
	bool cov=false; // tells me if max node coverage was determined

	GVec<int> lastg; // keeps guides that are included in other ones for last
	bool included=true;

	//***int g=0; // this is the option that starts from the most abundant guide to the least abundant
	int gi=guidetrf.Count()-1; // in this case I start with the least abundant so that I won't lose it -> seems counterintuitive
	bool first=true;
	//***while(g<guidetrf.Count()) { // this is the option that starts from the most abundant guide to the least abundant
	while(gi>=0) {
		int g=gi;
		if(included) {
			if(!guidetrf[g].trf->real) {
				lastg.Add(g);
				gi--;
				continue; // skip included guides first and process them later
			}
		}
		else { g=lastg[gi];}

		// weight the transcript
		GVec<float> nodeflux;
		float flux=0;
		//float fragno=0;

		//fprintf(stderr,"guide=%d ",g);

		if(EM) flux= max_flow_EM(gno,guidetrf[g].trf->nodes,istranscript,transfrag,no2gnode,nodeflux,guidetrf[g].trf->pattern);
		else if(weight) flux= weight_max_flow(gno,guidetrf[g].trf->nodes,istranscript,transfrag,no2gnode,nodeflux,guidetrf[g].trf->pattern);
		else flux= max_flow(gno,guidetrf[g].trf->nodes,istranscript,transfrag,no2gnode,nodeflux,guidetrf[g].trf->pattern);

		istranscript.reset();

		/*
		{ // DEBUG ONLY
			//fprintf(stderr,"flux[%d]=%g Path:",g,flux);for(int i=0;i<path.Count();i++) fprintf(stderr," %d",path[i]);
			fprintf(stderr,"flux[%d]=%g\n",g,flux);
			//fprintf(stderr,"\n");
		}
		*/

		if(flux>epsilon) {
			bool include=true;
			/*char *predid=NULL;
			int idlen=strlen(guidetrf[g].id);
			if(idlen) {
				predid=Gstrdup(guidetrf[g].id);
			}
			*/
			//store_transcript(pred,guidetrf[g].trf->nodes,nodeflux,nodecov,no2gnode,geneno,first,s,gno,include,pathpat,fragno,predid);
			store_transcript(pred,guidetrf[g].trf->nodes,nodeflux,nodecov,no2gnode,geneno,first,s,gno,gpos,include,pathpat,guidetrf[g].t);
			cov=true;
			// Node coverages:
			for(int i=1;i<gno-1;i++)
				if(nodecov[i]>nodecov[maxi]) maxi=i;

			//if(nodecov[maxi]<readthr) break; // no need to find other paths since they aren't any above allowed read threshold
			//if(nodecov[maxi]<1) break; // I shouldn't be restricting this at all?

			/*
			{ // DEBUG ONLY
				fprintf(stderr,"\nAfter update:\n");
				for(int i=0;i<gno;i++) {
					fprintf(stderr,"Node %d: %f ",i,nodecov[i]);
					fprintf(stderr,"trf=");
					for(int t=0;t<no2gnode[i]->trf.Count();t++) fprintf(stderr," %d(%f)",no2gnode[i]->trf[t],transfrag[no2gnode[i]->trf[t]]->abundance);
					fprintf(stderr," maxi=%d maxcov=%f\n",maxi,nodecov[maxi]);
				}
			}
			*/

			nodeflux.Clear();
		}
		//***g++; // this is the option that starts from the most abundant guide to the least abundant
		gi--;
		if(gi<0 && included) {
			gi=lastg.Count()-1;
			included=false;
		}
	} // end while(g<guidetrf.Count()

	if(!cov) for(int i=2;i<gno-1;i++)
		if(nodecov[i]>nodecov[maxi]) maxi=i;

	return(maxi);
}

bool is_reference_transcript(GVec<CGuide>& guidetrf,GBitVec& pattern) {
	int g=0;
	while(g<guidetrf.Count()){
		if((guidetrf[g].trf->pattern & pattern)==pattern) return true;
		g++;
	}
	return false;
}

bool sharedlink(int gno,GIntHash<int>& gpos,int g1,int g2,GVec<CGuide>& guidetrf,int i, int j) { // the two guides share the path between n1 and n2

	if(!guidetrf[g2].trf->pattern[guidetrf[g1].trf->nodes[i]]) return false;

	for(int k=i+1;k<=j;k++) {
		if(!guidetrf[g2].trf->pattern[guidetrf[g1].trf->nodes[k]]) return false;
		int *pos=gpos[edge(guidetrf[g1].trf->nodes[k-1],guidetrf[g1].trf->nodes[k],gno)];
		if(!pos || !guidetrf[g2].trf->pattern[*pos])
			return false;
	}

	return true;
}

int guides_maxflow(int gno,GIntHash<int>& gpos,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,GVec<CGuide>& guidetrf,int& geneno,
		int s,GList<CPrediction>& pred,GVec<float>& nodecov,GBitVec& istranscript,GBitVec& pathpat,bool &first) {


	int maxi=1;

	int ng=guidetrf.Count();

	if(ng==1) { // if only one guide I do not need to do the 2 pass
		GVec<float> nodeflux;
		//float fragno=0;
		float flux= max_flow(gno,guidetrf[0].trf->nodes,istranscript,transfrag,no2gnode,nodeflux,guidetrf[0].trf->pattern);
		istranscript.reset();

		/*
		{ // DEBUG ONLY
			fprintf(stderr,"guide=%s flux[0]=%g\n",guidetrf[0].t->getID(),flux);
		}
		*/

		if(flux>epsilon) {
			bool include=true;
			store_transcript(pred,guidetrf[0].trf->nodes,nodeflux,nodecov,no2gnode,geneno,first,s,gno,gpos,include,pathpat,guidetrf[0].t);

			nodeflux.Clear();
		}

		// Node coverages:
		for(int i=1;i<gno-1;i++)
			if(nodecov[i]>nodecov[maxi]) maxi=i;

		//if(nodecov[maxi]<readthr) break; // no need to find other paths since they aren't any above allowed read threshold
		//if(nodecov[maxi]<1) break; // I shouldn't be restricting this at all?


		/*
		{ // DEBUG ONLY
		  fprintf(stderr,"\nAfter update:\n");
		  for(int i=0;i<gno;i++) {
			  fprintf(stderr,"Node %d: %f ",i,nodecov[i]);
			  fprintf(stderr,"trf=");
			  for(int t=0;t<no2gnode[i]->trf.Count();t++) fprintf(stderr," %d(%f)",no2gnode[i]->trf[t],transfrag[no2gnode[i]->trf[t]]->abundance);
			  fprintf(stderr," maxi=%d maxcov=%f\n",maxi,nodecov[maxi]);
		  }
		}
		*/

		return(maxi);
	}

	bool cov=false; // tells me if max node coverage was determined

	GVec<float> flux;
	GVec<float> **capacity=new GVec<float>*[ng];
	GVec<float> **flow=new GVec<float>*[ng];
	GVec<int> **link=new GVec<int>*[ng];
	GVec<int> *node2path=new GVec<int>[ng];

	// calculate maximum flow for each guide
	for(int g=0;g<ng;g++) {
		int n=guidetrf[g].trf->nodes.Count();
		capacity[g]=new GVec<float>[n];
		flow[g]=new GVec<float>[n];
		link[g]=new GVec<int>[n];

		//fprintf(stderr,"guide=%s ",guidetrf[g].t->getID());

		float initflux=guideflow(gno,guidetrf[g].trf->nodes,istranscript,transfrag,no2gnode,guidetrf[g].trf->pattern,capacity[g],flow[g],link[g],node2path[g]);
		flux.Add(initflux);
		istranscript.reset();

		/*
		{ // DEBUG ONLY
			fprintf(stderr," flux[%d]=%g\n",g,flux[g]);
		}
		*/
	}


	// adjust guide flow
	for(int g=ng-1;g>=0;g--)
		if(flux[g]>epsilon) {
			//fprintf(stderr,"consider guide %d\n",g);
			bool adjust=false;
			for(int i=0;i<guidetrf[g].trf->nodes.Count();i++) { // for all nodes in guide g, recompute the capacities allowed
				int n1=guidetrf[g].trf->nodes[i];
				//fprintf(stderr,"n1=%d count=%d\n",n1,link[g][i].Count());
				for(int j=0;j<link[g][i].Count();j++) {
					int n2=guidetrf[g].trf->nodes[link[g][i][j]];
					//fprintf(stderr,"  path %d-%d capacity=%g\n",n1,n2,capacity[g][i][link[g][i][j]]);
					if(capacity[g][i][link[g][i][j]]) { // decrease capacity if there is any left already
						float usedflow=0;
						float varflow=0;
						for(int r=0;r<ng;r++) {
							if(r==g || sharedlink(gno,gpos,g,r,guidetrf,i,link[g][i][j])){ // the two guides share the path between n1 and n2

								//fprintf(stderr,"Guides %d and %d share the path between %d and %d with flows %g and %g\n",g,r,n1,n2,flow[g][i][link[g][i][j]],flow[r][node2path[r][n1]][node2path[r][n2]]);

								if(r>g) { // this is fixed flow
									usedflow+=flow[r][node2path[r][n1]][node2path[r][n2]];
								}
								else { // this is variable flow
									varflow+=flow[r][node2path[r][n1]][node2path[r][n2]];
								}
							}
						}
						if(varflow+usedflow>capacity[g][i][link[g][i][j]]) { // I need to adjust the flow
							if(varflow) capacity[g][i][link[g][i][j]]=flow[g][i][link[g][i][j]]*(capacity[g][i][link[g][i][j]]-usedflow)/varflow;
							else capacity[g][i][link[g][i][j]]=0; // this should never be the case
							adjust=true;
						}
						else {
							capacity[g][i][link[g][i][j]]=flow[g][i][link[g][i][j]]+capacity[g][i][link[g][i][j]]-usedflow-varflow; // this is how much I allow the flow to increase
						}
						//fprintf(stderr,"new capacity[%d][%d]=%g\n",i,link[g][i][j],capacity[g][i][link[g][i][j]]);
					}
				}
			}

			// recompute maxflow for the guide with adjustment for the new computed capacities only if adjust
			// is true, otherwise there is no need to, but I still need to update the abundances
			GVec<float> nodeflux;
			//float fragno=0;
			float newflux=guide_max_flow(adjust,guidetrf[g].trf->nodes,istranscript,transfrag,no2gnode,nodeflux,guidetrf[g].trf->pattern,capacity[g],flow[g],link[g],node2path[g]);
			if(!newflux) newflux=flux[g];
			istranscript.reset();

			//fprintf(stderr,"newflux=%g\n",newflux);

			if(newflux>epsilon) {
				bool include=true;
				store_transcript(pred,guidetrf[g].trf->nodes,nodeflux,nodecov,no2gnode,geneno,first,s,gno,gpos,include,pathpat,guidetrf[g].t);
				cov=true;
				// Node coverages:
				for(int i=1;i<gno-1;i++)
					if(nodecov[i]>nodecov[maxi]) maxi=i;

				//if(nodecov[maxi]<readthr) break; // no need to find other paths since they aren't any above allowed read threshold
				//if(nodecov[maxi]<1) break; // I shouldn't be restricting this at all?

				/*
				{ // DEBUG ONLY
				  fprintf(stderr,"\nAfter update:\n");
				  for(int i=0;i<gno;i++) {
					  fprintf(stderr,"Node %d: %f ",i,nodecov[i]);
					  fprintf(stderr,"trf=");
					  for(int t=0;t<no2gnode[i]->trf.Count();t++) fprintf(stderr," %d(%f)",no2gnode[i]->trf[t],transfrag[no2gnode[i]->trf[t]]->abundance);
					  fprintf(stderr," maxi=%d maxcov=%f\n",maxi,nodecov[maxi]);
				  }
				}
				*/

				nodeflux.Clear();
			}
		}

	// clean up memory

	for(int g=0;g<ng;g++) {
		delete [] capacity[g];
		delete [] flow[g];
		delete [] link[g];
	}
	delete [] capacity;
	delete [] flow;
	delete [] link;
	delete [] node2path;

	if(!cov) for(int i=2;i<gno-1;i++)
		if(nodecov[i]>nodecov[maxi]) maxi=i;

	return(maxi);
}

int find_transcripts(int gno,int edgeno, GIntHash<int> &gpos,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,GBitVec& compatible,
		int geneno,int strand,GVec<CGuide>& guidetrf,GList<CPrediction>& pred,bool fast) {

	// process in and out coverages for each node
	int maxi=0; // node with maximum coverage
	GVec<float> nodecov; // node coverages


	for(int i=0;i<gno;i++) {
		CGraphnode *inode=no2gnode[i]; // this is here only because of the DEBUG option below
		nodecov.cAdd(0.0);

		if(i) { // for all nodes but the source

		    if(i<gno-1 && inode->len()) nodecov[i]=inode->cov/inode->len(); // sink also has 0 coverage
		    if(nodecov[i]>nodecov[maxi]) maxi=i;
		    int nn=inode->trf.Count();
		    float abundin=0;
		    float abundout=0;
		    float abundthrough=0;
		    for(int j=0;j<nn;j++){
		    	int t=inode->trf[j];
		    	if(transfrag[t]->nodes.Last()==i) { // transfrag ends at this node (in transfrag)
		    		abundin+=transfrag[t]->abundance;
		    	}
		    	else if(transfrag[t]->nodes[0]==i) { // transfrag starts at this node (out transfrag)
		    		abundout+=transfrag[t]->abundance;
		    	}
		    	else if(transfrag[t]->pattern[i]) { // through transfrag (here I checked that the transfrag clearly goes through the node)
		    		abundthrough+=transfrag[t]->abundance;
		    	}
		    }

		    // EM case:
		    if(EM) {
		    	inode->rate=(abundout+abundthrough)/(abundin+abundthrough);
		    	inode->capacity=abundout+abundthrough;
		    }
		    else {
		    	if(abundin) inode->rate=abundout/abundin;
		    	if(abundout) inode->capacity=abundout+abundthrough; // node capacity tells me how much of that node coverage I can use give how many transfrags leave the node
		    	else inode->capacity=abundin+abundthrough;
		    }
		} // end if i

		/*
		{ // DEBUG ONLY
			printTime(stderr);
			fprintf(stderr,"Node %d: cov=%f capacity=%f rate=%f ",i,inode->cov/(inode->end-inode->start+1),inode->capacity,inode->rate);
			fprintf(stderr,"trf=");
			for(int t=0;t<inode->trf.Count();t++) fprintf(stderr," %d(%f)",inode->trf[t],transfrag[inode->trf[t]]->abundance);
			fprintf(stderr," maxi=%d maxcov=%f\n",maxi,nodecov[maxi]);
		}
		*/

	} // end for i

	GBitVec istranscript(transfrag.Count());
	GBitVec pathpat(gno+edgeno);

/*
#ifdef GMEMTRACE
	double vm,rsm;
	get_mem_usage(vm, rsm);
	GMessage("\t\tM(after istranscript and pathpat init):find_transcripts memory usage: rsm=%6.1fMB vm=%6.1fMB\n",rsm/1024,vm/1024);
#endif
*/

	// process guides first
	//fprintf(stderr,"guidetrf.count=%d\n",guidetrf.Count());
	//if(guidetrf.Count()) maxi=guides_flow(gno,no2gnode,transfrag,guidetrf,geneno,strand,pred,nodecov,istranscript,pathpat);

	bool first=true;

	//fprintf(stderr,"guide count=%d\n",guidetrf.Count());

	if(guidetrf.Count()) maxi=guides_maxflow(gno,gpos,no2gnode,transfrag,guidetrf,geneno,strand,pred,nodecov,istranscript,pathpat,first);


	if(nodecov[maxi]>=readthr) {

		// process rest of the transfrags
		if(nodecov[maxi]>=readthr) {
			GBitVec removable(transfrag.Count(),true);

			// 1:
			// parse_trf_weight_max_flow(gno,no2gnode,transfrag,geneno,strand,pred,nodecov,pathpat);
			// 2:
			GBitVec usednode(gno+edgeno);
			parse_trf(maxi,gno,edgeno,gpos,no2gnode,transfrag,compatible,geneno,first,strand,pred,nodecov,istranscript,removable,usednode,0,pathpat,fast);

		}
	}

	return(geneno);
}

/*
void get_trims(GVec<CTrimPoint>& trims,CBundlenode *currbnode,int refstart,GVec<float>& bpcov) {

	uint sourcestart;
	uint sinkend;
	float sinkabundance;
	float sourceabundance;

	while(currbnode!=NULL) {
		sourcestart=0;
		sinkend=0;
		sinkabundance=0;
		sourceabundance=0;
		find_trims(refstart,currbnode->start,currbnode->end,bpcov,sourcestart,sourceabundance,sinkend,sinkabundance);
		if(sourcestart<sinkend) {
			if(sourcestart) {
				CTrimPoint point(sourcestart,sourceabundance,true);
				trims.Add(point);
			}
			CTrimPoint point(sinkend,sinkabundance,false); // sinkend is always positive since it's > sourcestart
			trims.Add(point);
		}
		else {
			if(sinkend) {
				CTrimPoint point(sinkend,sinkabundance,false);
				trims.Add(point);
			}
			if(sourcestart) {
				CTrimPoint point(sourcestart,sourceabundance,true);
				trims.Add(point);
			}
		}
		currbnode=currbnode->nextnode;
	}
}
*/

void exon_covered(int ex,GffObj *guide,int &b,GPVec<CBundle>& bundle,GPVec<CBundlenode>& bnode,
		int& maxlen,int& leftlen,int& rightlen) {

	while(b<bundle.Count()){
		if(guide->exons[ex]->end<bnode[bundle[b]->startnode]->start) break;
		int ovlp=guide->exons[ex]->overlapLen(bnode[bundle[b]->startnode]->start,bnode[bundle[b]->lastnodeid]->end);
		if(ovlp>maxlen) maxlen=ovlp;
		if(guide->exons[ex]->start>=bnode[bundle[b]->startnode]->start) leftlen=ovlp;
		if(guide->exons[ex]->end<=bnode[bundle[b]->lastnodeid]->end) {
			rightlen=ovlp;
			break;
		}
		b++;
	}
}

void get_partial_covered(GffObj *guide,GPVec<CBundle>& bundle,GPVec<CBundlenode>& bnode,GList<CJunction>& junction) {

	int ntr=0;

	int s=-1;
	if(guide->strand=='+') s=1;

	int fex=0;

	int nj=0; // index of junctions
	int njunctions=junction.Count();
	int b=0;

	while(fex<guide->exons.Count()) {

		int maxlen=0;
		int leftlen=0;
		int rightlen=0;
		// first determine if exon fex is covered
		exon_covered(fex,guide,b,bundle,bnode,maxlen,leftlen,rightlen);

		if(rightlen) { // exon is covered up to the left
			int lex=fex;
			int trlen=rightlen;
			while(lex<guide->exons.Count()) {
				while(nj<njunctions && junction[nj]->start<guide->exons[lex]->end) nj++;
				if(nj==njunctions || junction[nj]->start>guide->exons[lex]->end) { // junction starting at lex is not covered
					if(trlen>mintranscriptlen || (lex==fex && maxlen>mintranscriptlen)) {
						bool partial=true;
						if(lex-fex+1==guide->exons.Count()) partial=false;
						for(int i=fex;i<=lex;i++)
							if(partial) fprintf(c_out,"%s\tpartial\texon\t%u\t%u\t.\t%c\t.\ttranscript_id \"%s_part%d\";\n",guide->getGSeqName(),
									guide->exons[i]->start,guide->exons[i]->end,guide->strand,guide->getID(),ntr);
							else fprintf(c_out,"%s\tcomplete\texon\t%u\t%u\t.\t%c\t.\ttranscript_id \"%s\";\n",guide->getGSeqName(),
									guide->exons[i]->start,guide->exons[i]->end,guide->strand,guide->getID());
						ntr++;
					}
					fex=lex+1;
					break;
				}
				else { // junction might be covered
					bool found=false;
					while(nj<njunctions && junction[nj]->start==guide->exons[lex]->end){
						if(lex+1<guide->exons.Count() && junction[nj]->strand==s && junction[nj]->end==guide->exons[lex+1]->start) {
							found=true;
							break;
						}
						nj++;
					}
					if(found) { // junction is covered
						lex++;
						maxlen=0;
						leftlen=0;
						rightlen=0;
						exon_covered(lex,guide,b,bundle,bnode,maxlen,leftlen,rightlen);
						trlen+=leftlen;
						if(!rightlen || leftlen!=rightlen) { // exon is not covered from left to right
							bool partial=true;
							if(lex-fex+1==guide->exons.Count()) partial=false;
							if(trlen>mintranscriptlen) {
								for(int i=fex;i<=lex;i++)
									if(partial) fprintf(c_out,"%s\tpartial\texon\t%u\t%u\t.\t%c\t.\ttranscript_id \"%s_part%d\";\n",guide->getGSeqName(),
											guide->exons[i]->start,guide->exons[i]->end,guide->strand,guide->getID(),ntr);
									else fprintf(c_out,"%s\tcomplete\texon\t%u\t%u\t.\t%c\t.\ttranscript_id \"%s\";\n",guide->getGSeqName(),
											guide->exons[i]->start,guide->exons[i]->end,guide->strand,guide->getID());
								ntr++;
							}
							if(!rightlen || !partial) fex=lex+1;
							else fex=lex;
							break;
						}
					}
					else { // couldn't find junction to the right
						if(trlen>mintranscriptlen || (lex==fex && maxlen>mintranscriptlen)) {
							bool partial=true;
							if(lex-fex+1==guide->exons.Count()) partial=false;
							for(int i=fex;i<=lex;i++)
								if(partial) fprintf(c_out,"%s\tpartial\texon\t%u\t%u\t.\t%c\t.\ttranscript_id \"%s_part%d\";\n",guide->getGSeqName(),
										guide->exons[i]->start,guide->exons[i]->end,guide->strand,guide->getID(),ntr);
								else fprintf(c_out,"%s\tcomplete\texon\t%u\t%u\t.\t%c\t.\ttranscript_id \"%s\";\n",guide->getGSeqName(),
										guide->exons[i]->start,guide->exons[i]->end,guide->strand,guide->getID());
							ntr++;
						}
						fex=lex+1;
						break;
					}
				}
			}
		}
		else { // current single exon is not covered up to the left -> maxlen is the exon length
			if(maxlen>mintranscriptlen) {
				bool partial=true;
				if(guide->exons.Count()==1) partial=false;
				if(partial) fprintf(c_out,"%s\tpartial\texon\t%u\t%u\t.\t%c\t.\ttranscript_id \"%s_part%d\";\n",guide->getGSeqName(),
						guide->exons[fex]->start,guide->exons[fex]->end,guide->strand,guide->getID(),ntr);
				else fprintf(c_out,"%s\tcomplete\texon\t%u\t%u\t.\t%c\t.\ttranscript_id \"%s\";\n",guide->getGSeqName(),
						guide->exons[fex]->start,guide->exons[fex]->end,guide->strand,guide->getID());
				ntr++;
			}
			fex++;
		}
	}

}

bool get_covered(GffObj *guide,GPVec<CBundle>& bundle,GPVec<CBundlenode>& bnode,GList<CJunction>& junction,
		GVec<int>* bnodeguides,int g) {

	bool covered=true;

	int s=-1;
	if(guide->strand=='+') s=1;

	// first check if all exons are covered by junctions
	int nj=0; // index of junctions
	int njunctions=junction.Count();
	for(int i=1;i<guide->exons.Count();i++) {
		while(nj<njunctions && junction[nj]->start<guide->exons[i-1]->end) nj++;
		if(nj==njunctions || junction[nj]->start>guide->exons[i-1]->end) {
			covered=false;
			break;
		}
		bool found=false;
		while(nj<njunctions && junction[nj]->start==guide->exons[i-1]->end){
			if(junction[nj]->strand==s && junction[nj]->end==guide->exons[i]->start) {
				found=true;
				break;
			}
			nj++;
		}
		if(!found) {
			covered=false;
			break;
		}
	}

	// now check if the exons are covered
	if(covered) {
		covered=false;
		int b=0;
		uint guidestart=0;
		uint maxguidelen=0;
		while(!covered && b<bundle.Count()) {
			if(guide->end<bnode[bundle[b]->startnode]->start) break;
			if(guide->overlap(bnode[bundle[b]->startnode]->start,bnode[bundle[b]->lastnodeid]->end)) {
				CBundlenode *currbnode=bnode[bundle[b]->startnode];
				covered=true;
				int i=0;
				while(currbnode!=NULL && i<guide->exons.Count()) {
					if(guide->exons[i]->end<currbnode->start) break;
					if(guide->exons.Count()==1) {
						uint ovlp=(uint)guide->overlapLen(currbnode->start,currbnode->end);
						if(maxguidelen<ovlp) {
								maxguidelen=ovlp;
								covered=true;
								i++;
								// break; // should I have a break in here?
						}
						if(ovlp && bnodeguides) {
							bnodeguides[currbnode->bid].Add(g);
						}
					}
					else {
						if(i<guide->exons.Count()-1) { // not last exon
							if(guide->exons[i]->end<=currbnode->end) { // exon end within currbnode
								if(i) { // not first exon
									if(currbnode->start<=guide->exons[i]->start && guide->exons[i]->start<=currbnode->end) { // exon included within currbnode
										i++;
										continue;
									}
									else break;
								}
								else {
									if(!guidestart) {
										if(currbnode->start<guide->exons[i]->start) guidestart=guide->exons[i]->start;
										else guidestart=currbnode->start;
									}
									i++;
									continue;
								}
							}
						}
						else { // last exon -> only start needs to be within currbnode
							if(currbnode->start<=guide->exons[i]->start && guide->exons[i]->start<=currbnode->end) {
								if(!maxguidelen) {
									if(guide->exons[i]->end<currbnode->end) maxguidelen=guide->exons[i]->end-guidestart+1;
									else maxguidelen=currbnode->end-guidestart+1;
								}
								i++; break;
							}
						}
					}
					currbnode=currbnode->nextnode;
				}
				if(i<guide->exons.Count()) covered=false;
			}
			b++;
		}

		//fprintf(stderr,"covered=%d guidelen=%d\n",covered,maxguidelen);

		if(covered && maxguidelen>=(uint)mintranscriptlen) { if(c_out) guide->printTranscriptGff(c_out);}
		else covered=false;
	}

	return(covered);
}

bool guide_exon_overlap(GPVec<GffObj>& guides,int sno,uint start,uint end) {

	// maybe this shouldn't link groups together that are clear borders of other nodes because then I don't have any read spanning the edges

	char strand='.';
	if(sno==2) strand='+';
	else if(sno==0) strand='-';

	for(int g=0;g<guides.Count();g++) {
		if((sno==1 || guides[g]->strand==strand) && guides[g]->overlap(start,end)) { // guide overlaps than look at the exons
			for(int i=0;i<guides[g]->exons.Count();i++) {
				if(end<guides[g]->exons[i]->start) break;  // there won't be any further overlap
				if(start<guides[g]->exons[i]->start) break;	// overlap is before start of exon
				if(end<=guides[g]->exons[i]->end) {

					//fprintf(stderr,"overlap btw %d-%d and exon %d-%d of guide %s\n",start,end,guides[g]->exons[i]->start,guides[g]->exons[i]->end,guides[g]->getID());

					return true; // start is biger than exon start
				}
			}
		}
	}

	return false;
}

bool good_junc(CJunction& jd,int refstart, GVec<float>* bpcov,GPVec<GffObj>& guides) {

	if(eonly && !jd.guide_match) { // this way I am using only reads that are compatible to annotated transcripts
		jd.strand=0;
		return false;
	}

	if(guides.Count() && jd.guide_match) return true; // this junction is covered by at least one read: the one that calls good_junc
													  // ^ KEEP THIS IN MIND IF WE CHANGE HOW WE USE GOOD_JUNC

	/*
	{ // DEBUG ONLY
		fprintf(stderr,"Junction(%d): %d-%d with support=%g with out coverages: %g-%g ",jd.strand,jd.start,jd.end,jd.nreads_good,
			bpcov[jd.start-refstart-1],bpcov[jd.end-refstart]);
		if(bpcov[jd.start-refstart-1] && bpcov[jd.end-refstart])
			fprintf(stderr,"(%g-%g) ",jd.nreads_good*100/bpcov[jd.start-refstart-1],jd.nreads_good*100/bpcov[jd.end-refstart]);
		else fprintf(stderr,"(x-x) ");
		fprintf(stderr," and in coverages: %g-%g ",bpcov[jd.start-refstart],bpcov[jd.end-refstart-1]);
		if(bpcov[jd.start-refstart] && bpcov[jd.end-refstart-1])
			fprintf(stderr,"(%g-%g) ",jd.nreads_good*100/bpcov[jd.start-refstart],jd.nreads_good*100/bpcov[jd.end-refstart-1]);
		else fprintf(stderr,"(x-x) ");
		fprintf(stderr," and drop in cov: ");
		if(bpcov[jd.start-refstart-1] && bpcov[jd.end-refstart])
			fprintf(stderr,"%g-%g ",bpcov[jd.start-refstart]*100/bpcov[jd.start-refstart-1],bpcov[jd.end-refstart-1]*100/bpcov[jd.end-refstart]);
		else fprintf(stderr,"x-x ");
		fprintf(stderr,"\n");
	}
	*/

	if (jd.nreads_good<junctionthr) {
		jd.strand=0;
		return false;
	}

	// if(!jd.strand) return false; the strand has to be non-zero when we call good_junc -> KEEP THIS IN MIND IF WE CHANGE HOW WE USE GOOD_JUNC

	// don't trust spliced reads that have a very low coverage:
	int sno=(int)jd.strand+1;
	float leftcov=bpcov[sno][jd.start-refstart-1];
	if(bpcov[1][jd.start-refstart-1])
		leftcov+=(bpcov[1][jd.start-refstart-1]-bpcov[0][jd.start-refstart-1]-bpcov[2][jd.start-refstart-1])/bpcov[1][jd.start-refstart-1];
	float rightcov=bpcov[sno][jd.start-refstart-1];
	if(bpcov[1][jd.start-refstart])
		rightcov+=(bpcov[1][jd.start-refstart]-bpcov[0][jd.start-refstart]-bpcov[2][jd.start-refstart])/bpcov[1][jd.start-refstart];
	if(leftcov && jd.nreads_good*100/leftcov<isofrac && rightcov/leftcov>1-isofrac) {
		jd.strand=0;
		return false;
	}

	leftcov=bpcov[sno][jd.end-refstart];
	if(bpcov[1][jd.end-refstart])
		leftcov+=(bpcov[1][jd.end-refstart]-bpcov[0][jd.end-refstart]-bpcov[2][jd.end-refstart])/bpcov[1][jd.end-refstart];
	rightcov=bpcov[sno][jd.end-refstart-1];
	if(bpcov[1][jd.end-refstart-1])
		rightcov+=(bpcov[1][jd.end-refstart-1]-bpcov[0][jd.end-refstart-1]-bpcov[2][jd.end-refstart-1])/bpcov[1][jd.end-refstart-1];
	if(leftcov && jd.nreads_good*100/leftcov<isofrac && rightcov/leftcov>1-isofrac) {
		jd.strand=0;
		return false;
	}

	/*
	if ((int)(jd.end-jd.start)>longintron) { // I might want to reconsider this
		int leftreach = jd.start-longintronanchor-refstart;
		if(leftreach<0) leftreach=0;
		int rightreach = jd.end+longintronanchor-refstart;
		if(rightreach>=bpcov.Count()) rightreach=bpcov.Count()-1;

		//fprintf(stderr,"bpcov[%d]=%g bpcov[%d]=%g bpcov[%d]=%g bpcov[%d]=%g\n",leftreach,bpcov[leftreach],jd.start-refstart-1,bpcov[jd.start-refstart-1],
		//		rightreach,bpcov[rightreach],jd.end-refstart,bpcov[jd.end-refstart]);

		if((bpcov[leftreach]<1 && bpcov[jd.start-refstart-1]<2) ||
				(bpcov[rightreach]<1 && bpcov[jd.end-refstart]<2)) {
			jd.strand=0; // the split read becomes false but some previous reads might have made it - why??
			return false;
		}
	}
	*/

	return true;
}

/*
// this function can disrupt the sorted order of reads in the readlist
bool continue_read(GList<CReadAln>& readlist,int n,int leftlen,int idx) {
	int rightlen=0;
	for(int i=idx+1;i<readlist[n]->segs.Count();i++) {
		rightlen+=readlist[n]->segs[i].len();
	}
	// keep longest part of read
	if(leftlen<rightlen) { // discard left part of read
		readlist[n]->start=readlist[n]->segs[idx+1].start;
		for(int i=idx;i>=0;i--) {
			readlist[n]->segs.Delete(i);
			readlist[n]->juncs.Delete(i);
		}
	}
	else { // discar right part of read
		readlist[n]->end=readlist[n]->segs[idx].end;
		for(int i=readlist[n]->segs.Count()-1;i>idx;i--) {
			readlist[n]->segs.Delete(i);
			readlist[n]->juncs.Delete(i-1);
		}
		return false;
	}
	return true;
}
*/


void continue_read(GList<CReadAln>& readlist,int n,int idx) {
	// keep longest part of read
	readlist[n]->end=readlist[n]->segs[idx].end;
	for(int i=readlist[n]->segs.Count()-1;i>idx;i--) {
		readlist[n]->segs.Delete(i);
		readlist[n]->juncs.Delete(i-1);
	}
}

/* this version will not work because the bpcov now is on three strands
void update_junction_counts(CReadAln & rd,GVec<float>& bpcov,int refstart) {
	int nex=rd.segs.Count();
	GVec<uint> leftsup;
	GVec<uint> rightsup;
	uint maxleftsupport=0;
	uint maxrightsupport=0;
	for(int i=0;i<nex;i++) {
		if(i) {
			if(rd.segs[i-1].len()>maxleftsupport) maxleftsupport=rd.segs[i-1].len();
			if(rd.segs[nex-i].len()>maxrightsupport) maxrightsupport=rd.segs[nex-i].len();
			leftsup.Add(maxleftsupport);
			rightsup.Add(maxrightsupport);
		}
		cov_add(bpcov,rd.segs[i].start-refstart,rd.segs[i].end-refstart,0-rd.read_count);
	}
	for(int i=1;i<nex;i++) {
		//uint anchor=junctionsupport;
		//if(anchor<longintron && rd.juncs[i-1]->len()>longintron) anchor=longintronanchor; // this if I want to use a different anchor for long introns
		//if(leftsup[i-1]>=anchor && rightsup[nex-i-1]>=anchor) {
		if(leftsup[i-1]>=junctionsupport && rightsup[nex-i-1]>=junctionsupport) { // this was a good read to support this junction
			rd.juncs[i-1]->nreads_good-=rd.read_count;
			if(!rd.juncs[i-1]->nreads_good) rd.juncs[i-1]->strand=0;
			//else if(!rd.juncs[i-1]->guide_match && rd.juncs[i-1]->nreads_good<junctionthr) rd.juncs[i-1]->strand=0; // I am deleteing this one because
		}
		else if(!rd.juncs[i-1]->nreads_good) rd.juncs[i-1]->strand=0;
	}
}
*/

void update_junction_counts(CReadAln & rd) {
	int nex=rd.segs.Count();
	GVec<uint> leftsup;
	GVec<uint> rightsup;
	uint maxleftsupport=0;
	uint maxrightsupport=0;
	for(int i=0;i<nex;i++) {
		if(i) {
			if(rd.segs[i-1].len()>maxleftsupport) maxleftsupport=rd.segs[i-1].len();
			if(rd.segs[nex-i].len()>maxrightsupport) maxrightsupport=rd.segs[nex-i].len();
			leftsup.Add(maxleftsupport);
			rightsup.Add(maxrightsupport);
		}
		//cov_add(bpcov,rd.segs[i].start-refstart,rd.segs[i].end-refstart,0-rd.read_count); removed this in order not to have inconsistencies between reads
	}
	for(int i=1;i<nex;i++) {

		rd.juncs[i-1]->nreads-=rd.read_count;
		if(!rd.juncs[i-1]->nreads) rd.juncs[i-1]->strand=0;

		if(leftsup[i-1]>=junctionsupport && rightsup[nex-i-1]>=junctionsupport) { // this was a good junction but for some reason the read was thrown away
			rd.juncs[i-1]->nreads_good-=rd.read_count;
			//if(!rd.juncs[i-1]->nreads_good) rd.juncs[i-1]->strand=0;
		}
		//else if(!rd.juncs[i-1]->nreads_good) rd.juncs[i-1]->strand=0;
	}
}


int build_graphs(BundleData* bdata, bool fast) {
	int refstart = bdata->start;
	GList<CReadAln>& readlist = bdata->readlist;
	GList<CJunction>& junction = bdata->junction;
	GPVec<GffObj>& guides = bdata->keepguides;
	GVec<float>* bpcov = bdata->bpcov; // I might want to use a different type of data for bpcov to save memory in the case of very long bundles
	GList<CPrediction>& pred = bdata->pred;
	// form groups on strands: all groups below are like this: 0 = negative strand; 1 = unknown strand; 2 = positive strand
	GPVec<CGroup> group;
	CGroup *currgroup[3]={NULL,NULL,NULL}; // current group of each type
	CGroup *startgroup[3]={NULL,NULL,NULL}; // start group of each type
	int color=0; // next color to assign
	GVec<int> merge; // remembers merged groups
	GVec<int> equalcolor; // remembers colors for the same bundle
	GVec<int> *readgroup=new GVec<int>[readlist.Count()]; // remebers groups for each read; don't forget to delete it when no longer needed

	//int **readgroup = new int*[readlist.Count()];

/*
#ifdef GMEMTRACE
	double vm,rsm;
	get_mem_usage(vm, rsm);
	GMessage("\t\tM(s):build_graphs memory usage: rsm=%6.1fMB vm=%6.1fMB\n",rsm/1024,vm/1024);
#endif
*/

	//float fraglen=0;
	//uint fragno=0;

	GHash<bool> boundaryleft;
	GHash<bool> boundaryright;

	for (int n=0;n<readlist.Count();n++) {
		CReadAln & rd=*(readlist[n]);

		//if(rd.juncs.Count()) fprintf(stderr,"read[%d] %d %d %d [", n, rd.start, rd.end, (int)rd.strand);

		bool keep=true;
		int i=0;
		//int leftlen=0;
		while(i<rd.juncs.Count()) {
			//leftlen=rd.segs[i].len();
			CJunction& jd=*(rd.juncs[i]);
			//fprintf(stderr, " %d-%d:%4.2f", jd.start, jd.end, jd.nreads_good);
			if(!jd.strand || !good_junc(jd,refstart,bpcov,guides)) { // found a bad junction
				//update_junction_counts(rd,bpcov,refstart);
				update_junction_counts(rd);
				keep=false;
				rd.nh=0;
				//fprintf(stderr," - %c -",jd.strand);
				break;
				/*
				// this is one version but is kind of ad-hoc --> don't like it very much: a better option might be to update junction counts
				if(i>0) { // the read has some junctions that were good before --> might want to keep them
					continue_read(readlist,n,i);
				}
				else {
					keep=false;
					rd.nh=0;
					// I should also update the junction counts here because otherwise I will have problems later
					//if(!continue_read(readlist,n,leftlen,i))
					break;
				}
				*/
			}
			else if(guides.Count()){ // need to remember boundary
				bool exist=true;
		    	GStr bs((int)jd.start);
		    	if(!boundaryleft[bs.chars()]) boundaryleft.Add(bs.chars(),new bool(exist));
		    	GStr be((int)jd.end);
		    	if(!boundaryright[be.chars()]) boundaryright.Add(be.chars(),new bool(exist));
			}
			i++;
		}
		//if(rd.juncs.Count()) fprintf(stderr,"] keep=%d\n",keep);
		if(keep) { // if it's a good read that needs to be kept
			//fprintf(stderr,"add read %d:%d-%d w/count=%g for color=%d with npairs=%d\n",n,readlist[n]->start,readlist[n]->end,readlist[n]->read_count,color,readlist[n]->pair_idx.Count());
			color=add_read_to_group(n,readlist,color,group,currgroup,startgroup,readgroup,equalcolor,merge);

			// count fragments
			bdata->frag_len+=rd.len*rd.read_count;
			double single_count=rd.read_count;
			for(int i=0;i<rd.pair_idx.Count();i++) {
				// I am not counting the fragment if I saw the pair before and it wasn't deleted
				if(rd.pair_idx[i]!=-1 && n>rd.pair_idx[i] && readlist[rd.pair_idx[i]]->nh) {// only if read is paired and it comes first in the pair I cound the fragments
					single_count-=rd.pair_count[i];
				}
			}
			if(single_count>epsilon) {
				bdata->num_fragments+=single_count;
			}


			//fprintf(stderr,"now color=%d\n",color);
			/* this part is now in add_read_to_group -> hopefully faster
			int np=readlist[n]->pair_idx; // pair read number
			if(np>n) { // I keep curent read
				for(int i=0;i<readlist[np]->segs.Count();i++) {
					fraglen+=readlist[np]->segs[i].len();
				}
				for(int i=0;i<readlist[n]->segs.Count();i++) {
					fraglen+=readlist[n]->segs[i].len();
				}
				fragno++;
			}
			else if(np==-1){
				for(int i=0;i<readlist[n]->segs.Count();i++) {
					fraglen+=readlist[n]->segs[i].len();
				}
				fragno++;
			}
			*/
		}
		//else clean_read_junctions(readlist[n]);
	}

	//fprintf(stderr,"fragno=%d fraglen=%g\n",fragno,fraglen);

	//if(fragno) fraglen/=fragno;

	// merge groups that are close together or groups that are within the same exon of a reference gene
	if(bundledist || guides.Count()) {
		for(int sno=0;sno<3;sno++) {
			CGroup *lastgroup=NULL;
			CGroup *procgroup=startgroup[sno];
			while(procgroup!=NULL) {

				if(lastgroup) {

					GStr bstart((int)lastgroup->end);
			    	GStr bend((int)procgroup->start);

			    	if(!boundaryleft[bstart.chars()] && !boundaryright[bend.chars()] && (procgroup->start-lastgroup->end<=bundledist ||
						(guides.Count()  && guide_exon_overlap(guides,sno,lastgroup->end,procgroup->start)))) {

			    		//fprintf(stderr,"sno=%d merge groups btw %d and %d dist=%d\n",sno,lastgroup->end,procgroup->start,procgroup->start-lastgroup->end);

						merge_fwd_groups(group,lastgroup,procgroup,merge,equalcolor);
						procgroup=lastgroup->next_gr;
						continue;
					}
				}
				lastgroup=procgroup;
				procgroup=procgroup->next_gr;
			}
		}
	}

	if(guides.Count()) {
		boundaryleft.Clear();
		boundaryright.Clear();
	}

	/*
	{ // DEBUG ONLY
		fprintf(stderr,"%d groups created!\n",group.Count());
	    for(int sno=0;sno<3;sno++) {
	    	fprintf(stderr, "Groups on strand %d:\n",sno);
	    	CGroup *procgroup=startgroup[sno];
	    	while(procgroup!=NULL) {
	    		fprintf(stderr, " gr %d(%d,%.6f): %d-%d",procgroup->grid,procgroup->color,procgroup->cov_sum,procgroup->start,procgroup->end);
	    		procgroup=procgroup->next_gr;
	    	}
	    	fprintf(stderr,"\n");
	    }
	}
	*/

/*
#ifdef GMEMTRACE
	//double vm,rsm;
	get_mem_usage(vm, rsm);
	GMessage("\t\tM(after groups created):build_graphs memory usage: rsm=%6.1fMB vm=%6.1fMB\n",rsm/1024,vm/1024);
#endif
*/

	// ### form bundles here

    // first do color assignment
	for (int i=0;i<3;i++) currgroup[i]=startgroup[i];
	CGroup *prevgroup[3]={NULL,NULL,NULL};

	GVec<int> eqposcol(true);
	GVec<int> eqnegcol(true);

	eqposcol.Resize(equalcolor.Count(),-1);
	eqnegcol.Resize(equalcolor.Count(),-1);

	/*
	GVec<int> bundlecol(true); // will need this later: associates a bundle number to a group color

	// init strand colors
	int grcol=-1;
	for(int i=0;i<equalcolor.Count();i++) {
		eqposcol.Add(grcol);
		eqnegcol.Add(grcol);
		bundlecol.Add(grcol);
	}
	*/

	// each unstranded group needs to remember what proportion of stranded group it overlaps so that it can distribute reads later on -> maybe I can do this in the following while?

	while(currgroup[0]!=NULL || currgroup[1]!=NULL || currgroup[2]!=NULL) { // there are still groups to process

		int nextgr=get_min_start(currgroup); // gets the index of currgroup with the left most begining

		int grcol = currgroup[nextgr]->color;    // set smallest color for currgroup[$nextgr]

		while(equalcolor[grcol]!=grcol) {
			grcol=equalcolor[grcol];
		}
		currgroup[nextgr]->color=grcol;

		//print STDERR "nextgr=$nextgr grcol=$grcol current group is at coords: ",$currgroup[$nextgr][0],"-",$currgroup[$nextgr][1],"\n";
		//fprintf(stderr,"group %d id=%d: %u-%u col=%d from col=%d\n",nextgr,currgroup[nextgr]->grid,currgroup[nextgr]->start,currgroup[nextgr]->end,grcol,prevcol);

		if(nextgr == 1) { // unknown strand group

			if(prevgroup[0]!=NULL && currgroup[nextgr]->start <= prevgroup[0]->end) { // overlaps previous negative group
				//fprintf(stderr,"\tovlp to neg group: %u-%u\n",prevgroup[0]->start,prevgroup[0]->end);
				set_strandcol(currgroup[nextgr],prevgroup[0],prevgroup[0]->color,eqnegcol,equalcolor);
				uint maxstart = currgroup[nextgr]->start > prevgroup[0]->start ? currgroup[nextgr]->start : prevgroup[0]->start;
				uint minend = currgroup[nextgr]->end < prevgroup[0]->end ? currgroup[nextgr]->end : prevgroup[0]->end;
				currgroup[nextgr]->neg_prop+=prevgroup[0]->cov_sum*(minend-maxstart+1)/prevgroup[0]->len();
			}

			while(currgroup[0]!=NULL && currgroup[nextgr]->start <= currgroup[0]->end && currgroup[0]->start <= currgroup[nextgr]->end) { // overlaps current negative strand group
				//fprintf(stderr,"\tovlp to neg group: %u-%u\n",currgroup[0]->start,currgroup[0]->end);

				int grcol = currgroup[0]->color;    // set smallest color for currgroup[$nextgr]
				while(equalcolor[grcol]!=grcol) {
					grcol=equalcolor[grcol];
				}
				currgroup[0]->color=grcol;
				currgroup[0]->neg_prop=1;

				set_strandcol(currgroup[nextgr],currgroup[0],currgroup[0]->color,eqnegcol,equalcolor);
				uint maxstart = currgroup[nextgr]->start > currgroup[0]->start ? currgroup[nextgr]->start : currgroup[0]->start;
				uint minend = currgroup[nextgr]->end < currgroup[0]->end ? currgroup[nextgr]->end : currgroup[0]->end;
				currgroup[nextgr]->neg_prop+=currgroup[0]->cov_sum*(minend-maxstart+1)/currgroup[0]->len();


				prevgroup[0]=currgroup[0];
				currgroup[0]=currgroup[0]->next_gr;
			}

			float pos_prop=0;
			if(prevgroup[2]!=NULL && currgroup[nextgr]->start <= prevgroup[2]->end) { // overlaps positive strand group
				//fprintf(stderr,"\tovlp to pos group: %u-%u\n",prevgroup[2]->start,prevgroup[2]->end);
				set_strandcol(currgroup[nextgr],prevgroup[2],prevgroup[2]->color,eqposcol,equalcolor);
				if(currgroup[nextgr]->neg_prop) {
					uint maxstart = currgroup[nextgr]->start > prevgroup[2]->start ? currgroup[nextgr]->start : prevgroup[2]->start;
					uint minend = currgroup[nextgr]->end < prevgroup[2]->end ? currgroup[nextgr]->end : prevgroup[2]->end;
					pos_prop+=prevgroup[2]->cov_sum*(minend-maxstart+1)/prevgroup[2]->len();
				}
			}

			while(currgroup[2]!=NULL && currgroup[nextgr]->start <= currgroup[2]->end && currgroup[2]->start <= currgroup[nextgr]->end) { // overlaps positive strand group
				//fprintf(stderr,"\tovlp to pos group: %u-%u\n",currgroup[2]->start,currgroup[2]->end);

				int grcol = currgroup[2]->color;    // set smallest color for currgroup[$nextgr]
				while(equalcolor[grcol]!=grcol) {
					grcol=equalcolor[grcol];
				}
				currgroup[2]->color=grcol;

				set_strandcol(currgroup[nextgr],currgroup[2],currgroup[2]->color,eqposcol,equalcolor);
				if(currgroup[nextgr]->neg_prop) {
					uint maxstart = currgroup[nextgr]->start > currgroup[2]->start ? currgroup[nextgr]->start : currgroup[2]->start;
					uint minend = currgroup[nextgr]->end < currgroup[2]->end ? currgroup[nextgr]->end : currgroup[2]->end;
					pos_prop+=currgroup[2]->cov_sum*(minend-maxstart+1)/currgroup[2]->len();
				}

				prevgroup[2]=currgroup[2];
				currgroup[2]=currgroup[2]->next_gr;
			}

			if(pos_prop) {
				currgroup[nextgr]->neg_prop/=(currgroup[nextgr]->neg_prop+pos_prop);
			}
			else if(currgroup[nextgr]->neg_prop) currgroup[nextgr]->neg_prop=1;
			//fprintf(stderr,"neg_prop=%g pos_prop=%g\n",currgroup[nextgr]->neg_prop,pos_prop);
		}
		else if(nextgr == 0) { // negative strand group
			currgroup[nextgr]->neg_prop=1;
		}


		/*
		if(nextgr == 0) { // negative strand group
			if(prevgroup[1]!=NULL && currgroup[nextgr]->start <= prevgroup[1]->end) { // overlaps an unknown strand group

				fprintf(stderr,"neg group: %u-%u overlaps neutral group: %u-%u\n",currgroup[nextgr]->start,currgroup[nextgr]->end,prevgroup[1]->start,prevgroup[1]->end);

				set_strandcol(prevgroup[1],currgroup[nextgr],grcol,eqnegcol,equalcolor);

				if(currgroup[nextgr]->grid==1789) fprintf(stderr,"new grcol[1789]=%d from intersect with neutral group %d\n",grcol,prevgroup[1]->grid);
				else if(grcol==2910492 && grcol!=currgroup[nextgr]->color) fprintf(stderr,"changed col 2910492 to %d at group %d intersect with neutral group %d\n",currgroup[nextgr]->color,currgroup[nextgr]->grid,prevgroup[1]->grid);

			}
		}
		else if(nextgr == 2) { // positive strand group
			if(prevgroup[1]!=NULL && currgroup[nextgr]->start <= prevgroup[1]->end) { // overlaps an unknown strand group

				fprintf(stderr,"pos group: %u-%u overlaps neutral group: %u-%u\n",currgroup[nextgr]->start,currgroup[nextgr]->end,prevgroup[1]->start,prevgroup[1]->end);

				set_strandcol(prevgroup[1],currgroup[nextgr],grcol,eqposcol,equalcolor);

			}
		}
		else { // unknown strand group
			if(prevgroup[0]!=NULL && currgroup[nextgr]->start <= prevgroup[0]->end) { // overlaps negative strand group

				fprintf(stderr,"neutral group: %u-%u overlaps neg group: %u-%u\n",currgroup[nextgr]->start,currgroup[nextgr]->end,prevgroup[0]->start,prevgroup[0]->end);
				int prevcol=prevgroup[0]->color;
				set_strandcol(currgroup[nextgr],prevgroup[0],prevgroup[0]->color,eqnegcol,equalcolor);

				if(prevgroup[0]->grid==1789) fprintf(stderr,"grcol[1789]=%d from prev intersect with neutral group %d\n",prevgroup[0]->color,currgroup[nextgr]->grid);
				else if(prevcol==2910492 && prevcol!=prevgroup[0]->color) fprintf(stderr,"changed col 2910492 to %d at group %d intersect with prev neutral group %d\n",prevgroup[0]->color,prevgroup[0]->grid,currgroup[nextgr]->grid);
			}
			if(prevgroup[2]!=NULL && currgroup[nextgr]->start <= prevgroup[2]->end) { // overlaps positive strand group

				fprintf(stderr,"neutral group: %u-%u overlaps pos group: %u-%u\n",currgroup[nextgr]->start,currgroup[nextgr]->end,prevgroup[2]->start,prevgroup[2]->end);

				set_strandcol(currgroup[nextgr],prevgroup[2],prevgroup[2]->color,eqposcol,equalcolor);
			}
		}
		*/

		prevgroup[nextgr]=currgroup[nextgr];
		currgroup[nextgr]=currgroup[nextgr]->next_gr;

    }

	/*
    { // DEBUG ONLY
    	fprintf(stderr,"Colors assigned!\n");
    	for(int sno=0;sno<3;sno++) {
    		fprintf(stderr, "Colors of groups on strand %d:\n",sno);
    		CGroup *procgroup=startgroup[sno];
    		while(procgroup!=NULL) {

    			int grcol = procgroup->color;

    			while(equalcolor[grcol]!=grcol) {
    				grcol=equalcolor[grcol];
    			}

    			int negcol=eqnegcol[grcol];
    			if(eqnegcol[grcol]!=-1){
    				while(equalcolor[negcol]!=negcol) {
    					negcol=equalcolor[negcol];
    				}
    			}
    			int poscol=eqposcol[grcol];
    			if(eqposcol[grcol]!=-1){
    				while(equalcolor[poscol]!=poscol) {
    					poscol=equalcolor[poscol];
    				}
    			}

    			fprintf(stderr, " gr %d(%d,%d,%d,%.6f): %d-%d\n",procgroup->grid,grcol,negcol,poscol,procgroup->cov_sum,procgroup->start,procgroup->end);
    			procgroup=procgroup->next_gr;
    		}
    		//fprintf(stderr,"\n");
    	}
    	exit(0);
    }
    */


	// create bundles : bundles collect connected groups (with same color)
	for (int i=0;i<3;i++) {
		currgroup[i]=startgroup[i];
		prevgroup[i]=NULL;
	}

	GPVec<CBundle> bundle[3]; // all bundles on all strands: 0,1,2
	GPVec<CBundlenode> bnode[3]; // last bnodes on all strands: 0,1,2 for each bundle : this might be the key for overalps

	GVec<int> group2bundle[3]; // to retrace reads from group no to bundle
	for(int sno=0;sno<3;sno++) {
		group2bundle[sno].Resize(group.Count(),-1);  // for a given group id we get a certain bnode id
		bnode[sno].setFreeItem(false);
	}

	GVec<int> bundlecol(true); // associates a bundle number to a group color
	bundlecol.Resize(equalcolor.Count(),-1);

	while(currgroup[0]!=NULL || currgroup[1]!=NULL || currgroup[2]!=NULL) { // there are still groups to process

		int nextgr=get_min_start(currgroup);  // next group based on starting position

		// get group color; I need to redo this to ensure I equalize all colors -> they could still be hanged by set_strandcol
		int grcol = currgroup[nextgr]->color;

		while(equalcolor[grcol]!=grcol) {
			grcol=equalcolor[grcol];
		}
		currgroup[nextgr]->color=grcol;

		if(nextgr == 0 || nextgr ==2 || (nextgr==1 &&(eqnegcol[grcol]==-1) && (eqposcol[grcol]==-1))) { // negative or positive strand bundle or unstranded bundle

			int bno=bundlecol[grcol];

			if(bno>-1) { // bundle for group has been created before
				add_group_to_bundle(currgroup[nextgr],bundle[nextgr][bno],bnode[nextgr]);
			}
			else { // create new bundle
				bno=create_bundle(bundle[nextgr],currgroup[nextgr],bnode[nextgr]);
				bundlecol[grcol]=bno;
			}

			group2bundle[nextgr][currgroup[nextgr]->grid]=bundle[nextgr][bno]->lastnodeid;

			/*
			GStr id(currgroup[nextgr]->grid);id+=':';id+=nextgr;
			group2bundle.Add(id.chars(),bnode[nextgr][bno]);
			*/
		}
		else { // unknown strand : here is where I should compute positive and negative proportions

			if(eqnegcol[grcol]!=-1){
				int negcol=eqnegcol[grcol];
				while(equalcolor[negcol]!=negcol) {
					negcol=equalcolor[negcol];
				}

				int bno=bundlecol[negcol];
				if(bno>-1) { // bundle for group has been created before
					// print STDERR "Add neutral group ",$currgroup[$nextgr][3]," to neg bundle $bno\n";
					add_group_to_bundle(currgroup[nextgr],bundle[0][bno],bnode[0]);
				}
				else { // create new bundle
					bno=create_bundle(bundle[0],currgroup[nextgr],bnode[0]);
					//print STDERR "Added neutral group ",$currgroup[$nextgr][3]," to new neg bundle $bno\n";
					bundlecol[negcol]=bno;
				}
				group2bundle[0][currgroup[nextgr]->grid]=bundle[0][bno]->lastnodeid;
				/*
				GStr id(currgroup[nextgr]->grid);id+=':';id+='0';
				group2bundle.Add(id.chars(),bnode[0][bno]);
				*/
			} // if(eqnegcol[grcol]!=-1)

			if(eqposcol[grcol]!=-1){
				int poscol=eqposcol[grcol];
				while(equalcolor[poscol]!=poscol) {
					poscol=equalcolor[poscol];
				}

				int bno=bundlecol[poscol];
				if(bno>-1) { // bundle for group has been created before
					//print STDERR "Add neutral group ",$currgroup[$nextgr][3]," to pos bundle $bno\n";
					add_group_to_bundle(currgroup[nextgr],bundle[2][bno],bnode[2]);
				}
				else { // create new bundle
					bno=create_bundle(bundle[2],currgroup[nextgr],bnode[2]);
					//print STDERR "Added neutral group ",$currgroup[$nextgr][3]," to new pos bundle $bno\n";
					bundlecol[poscol]=bno;
				}
				group2bundle[2][currgroup[nextgr]->grid]=bundle[2][bno]->lastnodeid;
				/*
				GStr id(currgroup[nextgr]->grid);id+=':';id+='2';
				group2bundle.Add(id.chars(),bnode[2][bno]);
				*/
			}
		}

		// print STDERR "Done with groups: ",$currgroup[0],",",$currgroup[1],",",$currgroup[2],"\n";
		currgroup[nextgr]=currgroup[nextgr]->next_gr;

	} // while(currgroup[0]!=NULL || currgroup[1]!=NULL || currgroup[2]!=NULL)

	//print STDERR "Bundles created!\n";


	// Clean up no longer needed variables
	// group.Clear(); maybe I still need this?
	equalcolor.Clear();
	eqposcol.Clear();
	eqnegcol.Clear();
	bundlecol.Clear();

/*
#ifdef GMEMTRACE
	//double vm,rsm;
	get_mem_usage(vm, rsm);
	GMessage("\t\tM(after bundles created):build_graphs memory usage: rsm=%6.1fMB vm=%6.1fMB\n",rsm/1024,vm/1024);
#endif
*/

	// next variables are in order to remember if I get guide coverage
	GVec<int>* bnodeguides=NULL;
	if(bundle[1].Count() && bnode[1].Count()) {
		bnodeguides = new GVec<int>[bnode[1].Count()];
	}

	//if(guides.Count()) fprintf(stderr,"No of guides=%d partialcov=%d\n",guides.Count(),partialcov);

	/* if(partialcov) {
		for(int g=0;g<guides.Count();g++) {
			int s=0;
			if(guides[g]->strand=='+') s=2;
			get_partial_covered(guides[g],bundle[s],bnode[s],junction);
		}
		return(0);
	}
    */
	if(c_out || (bundle[1].Count() && bnode[1].Count())) // coverage is needed
		for(int g=0;g<guides.Count();g++) {
			//fprintf(stderr,"consider guide %d\n",g);
			int s=0;
			if(guides[g]->strand=='+') s=2;
			if((c_out && !get_covered(guides[g],bundle[s],bnode[s],junction,NULL,0)) ||
					(bundle[1].Count() && bnode[1].Count() && guides[g]->exons.Count()==1))
				get_covered(guides[g],bundle[1],bnode[1],junction,bnodeguides,g);
		}

	/*
	{ // DEBUG ONLY
		printTime(stderr);
		fprintf(stderr, "There are %d unstranded bundles %d negative bundles and %d positive bundles\n",bundle[1].Count(),bundle[0].Count(),bundle[2].Count());
		for(int sno=0;sno<3;sno++) {
			fprintf(stderr, "Bundles on strand %d:\n",sno);
			for(int b=0;b<bundle[sno].Count();b++) {
				int elen=bnode[sno][bundle[sno][b]->lastnodeid]->end-bnode[sno][bundle[sno][b]->startnode]->start+1;
				CBundlenode *currbnode=bnode[sno][bundle[sno][b]->startnode];
				fprintf(stderr,"***Bundle %d with len=%d:",b,elen);
				while(currbnode!=NULL) {
					fprintf(stderr, " %d-%d cov=%.6f",currbnode->start,currbnode->end,currbnode->cov/(currbnode->end-currbnode->start+1));
					currbnode=currbnode->nextnode;
				}
				currbnode=bnode[sno][bundle[sno][b]->lastnodeid];
				fprintf(stderr," last node:%d-%d\n",currbnode->start,currbnode->end);
			}
		}
	}
	*/

	int geneno=0;


    // ### predict transcripts for unstranded bundles here
	//if(fraglen)
	for(int b=0;b<bundle[1].Count();b++) {

    	if(bundle[1][b]->nread && (bundle[1][b]->multi/bundle[1][b]->nread)<=mcov && (guides.Count() || bundle[1][b]->len > mintranscriptlen)) { // there might be small transfrags that are worth showing, but here I am ignoring them

    		// bundle might contain multiple fragments of a transcript but since we don't know the complete structure -> print only the pieces that are well represented
    		CBundlenode *currbnode=bnode[1][bundle[1][b]->startnode];
    		int t=1;
    		while(currbnode!=NULL) {
    			int len=currbnode->end-currbnode->start+1;
    			float cov=currbnode->cov/(currbnode->end-currbnode->start+1);

    			bool printguides=false;

    			for(int i=0;i<bnodeguides[currbnode->bid].Count();i++) {
    				int g=bnodeguides[currbnode->bid][i];
    				geneno++;
    				int glen=guides[g]->end-guides[g]->start+1;
    				if(glen && guides[g]->exons.Count()==1) {
    					RC_TData* tdata=(RC_TData*)(guides[g]->uptr);
    					float gcov=(tdata->t_exons[0])->movlcount/glen;
    					// if(cov<gcov) gcov=cov; WHY DO I DO THIS?? CHECK!!!
    					CPrediction *p=new CPrediction(geneno-1, guides[g], guides[g]->start, guides[g]->end, gcov, guides[g]->strand, glen);
    					GSeg exon(guides[g]->start, guides[g]->end);
    					p->exons.Add(exon);
    					p->exoncov.Add(gcov);
    					pred.Add(p);
    					printguides=true;
    				}
    			}

    			/* char *predid=NULL;
    			if(!id.is_empty()) {
    				predid=Gstrdup(id.chars());
    			}
    			if((!eonly && cov>=readthr && len>mintranscriptlen) || (!id.is_empty())) {
    			*/
    			if(!printguides && cov>=readthr && len>mintranscriptlen) {
    				if(t==1) { geneno++;}
    				char sign='.';
    				//CPrediction *p=new CPrediction(geneno-1,predid,currbnode->start,currbnode->end,cov,sign,cov,fraglen);
    				CPrediction *p=new CPrediction(geneno-1, NULL, currbnode->start, currbnode->end, cov, sign, cov);
    				GSeg exon(currbnode->start,currbnode->end);
    				p->exons.Add(exon);
    				p->exoncov.Add(cov);
    				pred.Add(p);
    				t++;
    			}
    			currbnode=currbnode->nextnode;
    		}
    	}
    }

    //fprintf(stderr,"Done with unstranded bundles\n");
    if (bnodeguides) delete[] bnodeguides;

	// ### build graphs for stranded bundles here
    if(startgroup[0]!=NULL || startgroup[2]!=NULL) { //# there are stranded groups to process

    	// I don't need the groups here anymore : I only use their numbers
    	// group.Clear(); // I will need the proportions later on

    	// sort junctions -> junctions are sorted already according with their start, but not their end
    	GList<CJunction> ejunction(junction);
    	ejunction.setFreeItem(false);
    	if(ejunction.Count()) ejunction.setSorted(juncCmpEnd);

    	GVec<CGraphinfo> *bundle2graph[2]; // should I keep the neutral strand for consistency ? -> remember not to delete it
    	GVec<int> graphno[2];  // how many nodes are in a certain graph g, on strand s: graphno[s][g]
    	GVec<int> edgeno[2];  // how many edges are in a certain graph g, on strand s: edgeno[s][g]
    	// GVec<int> trnumber[2]; // how many transfrags are on a strand s, in a graph g -> I can find out this from transfrag[s][g].Count()
    	// int ngraph[2]={0,0};   // how many graphs are in each strand: negative (0), or positive(1) -> keep one for each bundle
    	GPVec<CTransfrag> *transfrag[2]; // for each transfrag t on a strand s, in a graph g, transfrag[s][g][t] gives it's abundance and it's pattern
    	GPVec<CGraphnode> *no2gnode[2]; // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i
    	CTreePat **tr2no[2]; // for each graph g, on a strand s, tr2no[s][g] keeps the tree pattern structure for quick retrieval of the index t of a tansfrag
    	GIntHash<int> *gpos[2]; // for each graph g, on a strand s, gpos[s][g] keeps the hash between edges and positions in the bitvec associated to a pattern
    	GVec<int> lastgpos[2];

    	int bno[2]={0,0};


    	// build graph structure
    	for(int sno=0;sno<3;sno+=2) { // skip neutral bundles -> those shouldn't have junctions

        	// guides appear to be sorted by start --> CHECK THIS!!
        	int g=0;
        	int ng=guides.Count();

    		int s=sno/2; // adjusted strand due to ignoring neutral strand
    		char strnd='-';
    		if(s) strnd='+';

    		bundle2graph[s]=NULL;
    		if(bnode[sno].Count()) bundle2graph[s]=new GVec<CGraphinfo>[bnode[sno].Count()];
    		transfrag[s]=NULL;
    		no2gnode[s]=NULL;
    		tr2no[s]=NULL;
    		gpos[s]=NULL;

    		if(bundle[sno].Count()) {
    			transfrag[s]=new GPVec<CTransfrag>[bundle[sno].Count()]; // for each bundle I have a graph ? only if I don't ignore the short bundles
    			no2gnode[s]=new GPVec<CGraphnode>[bundle[sno].Count()];

    			gpos[s]=new GIntHash<int>[bundle[sno].Count()];


    			GCALLOC(tr2no[s],bundle[sno].Count()*sizeof(CTreePat *));
    			bno[s]=bundle[sno].Count();

    			for(int b=0;b<bundle[sno].Count();b++) {
    				graphno[s].cAdd(0);
    				edgeno[s].cAdd(0);
    				lastgpos[s].cAdd(0);
    				// I am overestmating the edgeno below, hopefully not by too much

    				while(g<ng && guides[g]->end<bnode[sno][bundle[sno][b]->startnode]->start) g++;
    				int cg=g;

    				while(cg<ng && guides[cg]->start<=bnode[sno][bundle[sno][b]->lastnodeid]->end) { // this are potential guides that might overlap the current bundle, and they might introduce extra edges
    					if(guides[cg]->strand==strnd) edgeno[s][b]+=2; // this is an overestimate: possibly I have both an extra source and an extra sink link
    					cg++;
    				}

    				/*
    				{ // DEBUG ONLY
    				if(bundle[sno][b]->nread) {
    					fprintf(stderr,"proc bundle[%d][%d] %f/%f is %f len=%d\n",sno,b,bundle[sno][b]->multi,bundle[sno][b]->nread,(float)bundle[sno][b]->multi/bundle[sno][b]->nread,bundle[sno][b]->len);
    				} }
    				*/

    				// here I can add something in stringtie to lower the mintranscript len if there are guides?

    				if(bundle[sno][b]->nread && (bundle[sno][b]->multi/bundle[sno][b]->nread)<=mcov && bundle[sno][b]->len >= mintranscriptlen ) { // bundle is worth processing: it might be that there are small transfrags from source to sink that are worth processing

    					/*
    					// first identify drops/spikes in coverage points
    					GVec<CTrimPoint> trims;
    					if(trim) {
    						get_trims(trims,bnode[sno][bundle[sno][b]->startnode],refstart,bpcov);
    						{ // DEBUG ONLY
    							fprintf(stderr,"Trims found:");
    							for(int z=0;z<trims.Count();z++) fprintf(stderr," %d(%d,%f)",trims[z].pos,trims[z].start,trims[z].abundance);
    							fprintf(stderr,"\n");
    						}
    					}
    					graphno[s][b]=create_graph(s,b,bundle[sno][b],bnode[sno],junction,ejunction,
    							bundle2graph,no2gnode,transfrag,trims); // also I need to remember graph coverages somewhere -> probably in the create_graph procedure
    					*
    					*/
    					// create graph then

    					graphno[s][b]=create_graph(refstart,s,b,bundle[sno][b],bnode[sno],junction,ejunction,
    							bundle2graph,no2gnode,transfrag,gpos,bpcov,edgeno[s][b],lastgpos[s][b]); // also I need to remember graph coverages somewhere -> probably in the create_graph procedure

    					if(graphno[s][b]) tr2no[s][b]=construct_treepat(graphno[s][b],gpos[s][b],transfrag[s][b]);
    					else tr2no[s][b]=NULL;
    				}
    				else tr2no[s][b]=NULL;
    			}
    		}
    	}
    	//fprintf(stderr,"Done creating graphs\n");

    	/*
    	{ // DEBUG ONLY
    		printTime(stderr);
    		for(int s=0;s<2;s++) {
    			fprintf(stderr, "There are %d stranded[%d] graphs\n",bno[s],int(2*s));
    			for(int b=0;b<bno[s];b++) {
    				if(graphno[s][b]) {
    					GStr pat;
    					fprintf(stderr,"Graph[%d][%d] with %d nodes and %d edges with lastgpos=%d:",int(2*s),b,graphno[s][b],edgeno[s][b],lastgpos[s][b]);
    					for(int nd=1;nd<graphno[s][b]-1;nd++)
    						fprintf(stderr," %d(%d-%d)",nd,no2gnode[s][b][nd]->start,no2gnode[s][b][nd]->end);
    					fprintf(stderr,"\n");
    					//print_pattern(tr2no[s][b],pat,graphno[s][b]);
    				}
    			}
    		}
    	}
    	*/

/*
#ifdef GMEMTRACE
    	double vm,rsm;
    	get_mem_usage(vm, rsm);
    	GMessage("\t\tM(after graphs created):build_graphs memory usage: rsm=%6.1fMB vm=%6.1fMB\n",rsm/1024,vm/1024);
#endif
*/

		// I can clean up some data here:
    	for(int sno=0;sno<3;sno++) {
    		int n=bnode[sno].Count();
    		for(int b=0;b<n;b++) delete bnode[sno][b];
    		bnode[sno].Clear();
    		bundle[sno].Clear();
    	}


    	//fprintf(stderr,"Analyzing memory for %d read clusters\n",readlist.Count());

    	// because of this going throu
    	// compute probabilities for stranded bundles
    	for (int n=0;n<readlist.Count();n++) {

    		float single_count=readlist[n]->read_count;
    		for(int j=0; j<readlist[n]->pair_idx.Count();j++) {
    			int np=readlist[n]->pair_idx[j];
    			if(np>-1) {
    				single_count-=readlist[n]->pair_count[j];
    				if(n<np) {
    					get_fragment_pattern(readlist,n,np,readlist[n]->pair_count[j],readgroup,merge,group2bundle,bundle2graph,graphno,edgeno,gpos,no2gnode,transfrag,tr2no,group);
    				}
    			}
    		}
    		if(single_count>epsilon) {
    			get_fragment_pattern(readlist,n,-1,single_count,readgroup,merge,group2bundle,bundle2graph,graphno,edgeno,gpos,no2gnode,transfrag,tr2no,group);
    		}
    	}


/*
#ifdef GMEMTRACE
    	//double vm,rsm;
    	get_mem_usage(vm, rsm);
    	GMessage("\t\tM(read patterns counted):build_graphs memory usage: rsm=%6.1fMB vm=%6.1fMB\n",rsm/1024,vm/1024);
#endif
*/
    	// shouldn't readlist be also cleared up here? maybe bpcov too?

    	// don't forget to clean up the allocated data here
    	delete [] readgroup;
    	group.Clear();



    	// parse graph
    	for(int s=0;s<2;s++) {

    		for(int b=0;b<bno[s];b++) {
    			if(graphno[s][b]) {

    				// include source to guide starts links
    				GVec<CGuide> guidetrf;

    				if(guides.Count()) process_refguides(graphno[s][b],edgeno[s][b],gpos[s][b],lastgpos[s][b],no2gnode[s][b],transfrag[s][b],s,guides,guidetrf);

    				//process transfrags to eliminate noise, and set compatibilities, and node memberships
    				GBitVec compatible((1+transfrag[s][b].Count())*transfrag[s][b].Count()/2); // I might want to change this to gbitvec
    				process_transfrags(graphno[s][b],no2gnode[s][b],transfrag[s][b],tr2no[s][b],compatible,edgeno[s][b],gpos[s][b]);

    				/*
    				{ //DEBUG ONLY
    					//printTime(stderr);
    					fprintf(stderr,"There are %d nodes for graph[%d][%d]:\n",graphno[s][b],s,b);
    					for(int i=0;i<graphno[s][b];i++) {
    						fprintf(stderr,"%d (%d-%d): %f len=%d cov=%f",i,no2gnode[s][b][i]->start,no2gnode[s][b][i]->end,no2gnode[s][b][i]->cov,no2gnode[s][b][i]->len(),no2gnode[s][b][i]->cov/no2gnode[s][b][i]->len());
    						fprintf(stderr," parents:");
    						for(int j=0;j<no2gnode[s][b][i]->parent.Count();j++) fprintf(stderr," %d",no2gnode[s][b][i]->parent[j]);
    						fprintf(stderr," trf=");
    						for(int j=0;j<no2gnode[s][b][i]->trf.Count();j++) fprintf(stderr," %d",no2gnode[s][b][i]->trf[j]);
    						fprintf(stderr,"\n");
    					}
    					fprintf(stderr,"There are %d transfrags[%d][%d]:\n",transfrag[s][b].Count(),s,b);
    					for(int t=0;t<transfrag[s][b].Count();t++) {
    						fprintf(stderr,"%d: ",t);
    						printBitVec(transfrag[s][b][t]->pattern);
    						fprintf(stderr," %f\n",transfrag[s][b][t]->abundance);
    					}

    				}
    				*/

/*
#ifdef GMEMTRACE
    				double vm,rsm;
    				get_mem_usage(vm, rsm);
    				GMessage("\t\tM(after process_transfrags):build_graphs memory usage: rsm=%6.1fMB vm=%6.1fMB\n",rsm/1024,vm/1024);
#endif
*/

    				// find transcripts now
    				geneno=find_transcripts(graphno[s][b],edgeno[s][b],gpos[s][b],no2gnode[s][b],transfrag[s][b],compatible,
    						geneno,s,guidetrf,pred,fast);

    				for(int g=0;g<guidetrf.Count();g++) {
    					//GFREE(guidetrf[g].trf);
    					delete guidetrf[g].trf;
    				}

    				/*
    				{ //DEBUG ONLY
    					printTime(stderr);
    					fprintf(stderr,"Processed transcripts for s=%d b=%d\n",s,b);
    				}
    				*/

/*
#ifdef GMEMTRACE
    				//double vm,rsm;
    				get_mem_usage(vm, rsm);
    				GMessage("\t\tM(after processed transcripts):build_graphs memory usage: rsm=%6.1fMB vm=%6.1fMB\n",rsm/1024,vm/1024);
#endif
*/


    			}
    			// clean up what can be cleaned
    			if(tr2no[s][b]) free_treepat(tr2no[s][b]);
    		}

    		// final clean up: no2gnode, no2tr, transfrag, bundle2graph
    		if(bundle2graph[s]) delete [] bundle2graph[s];
    		if(transfrag[s]) delete [] transfrag[s];
    		if(no2gnode[s]) delete [] no2gnode[s];
    		if(gpos[s]) delete [] gpos[s];
    		if(tr2no[s]) GFREE(tr2no[s]);
    	}

	} // end if(startgroup[0]!=NULL || startgroup[2]!=NULL)
    else {

    	delete [] readgroup;
    	// clean up readgroup, bundle
    	for(int sno=0;sno<3;sno++) {
    		int n=bnode[sno].Count();
    		for(int b=0;b<n;b++) delete bnode[sno][b];
    		bnode[sno].Clear();
    	}
    }

/*
#ifdef GMEMTRACE
    //double vm,rsm;
    get_mem_usage(vm, rsm);
	GMessage("\t\tM(e):build_graphs memory usage: rsm=%6.1fMB vm=%6.1fMB\n",rsm/1024,vm/1024);
#endif
*/

    // don't forget to clean up the allocated data here
    return(geneno);
}

/* initial version that does not clean long introns
void clean_junctions(GList<CJunction>& junction) {

	fprintf(stderr,"Clean junctions:\n");
	for(int i=0;i<junction.Count();i++) {
		CJunction& jd=*(junction[i]);
		if(jd.nreads_good<junctionthr) {
			jd.strand=0;
		}
		{ // DEBUG ONLY
			fprintf(stderr,"Junction %d: %d %d %d\n",i,jd.start,jd.end,jd.strand);
		}
	}
}
*/

void count_good_junctions(GList<CReadAln>& readlist, int refstart, GVec<float>* bpcov) {

	for(int n=0;n<readlist.Count();n++) {
		CReadAln & rd=*(readlist[n]);
		int nex=rd.segs.Count();
		GVec<uint> leftsup;
		GVec<uint> rightsup;
		uint maxleftsupport=0;
		uint maxrightsupport=0;
		int sno=(int)rd.strand+1;
		for(int i=0;i<nex;i++) {
			if(i) {
				if(rd.segs[i-1].len()>maxleftsupport) maxleftsupport=rd.segs[i-1].len();
				if(rd.segs[nex-i].len()>maxrightsupport) maxrightsupport=rd.segs[nex-i].len();
				leftsup.Add(maxleftsupport);
				rightsup.Add(maxrightsupport);
				rd.juncs[i-1]->nreads+=rd.read_count;
			}
			cov_add(bpcov,sno,rd.segs[i].start-refstart,rd.segs[i].end-refstart,rd.read_count);
		}
		for(int i=1;i<nex;i++) {
			//uint anchor=junctionsupport;
			//if(anchor<longintron && rd.juncs[i-1]->len()>longintron) anchor=longintronanchor; // this if I want to use a different anchor for long introns
			//if(leftsup[i-1]>=anchor && rightsup[nex-i-1]>=anchor) rd.juncs[i-1]->nreads_good+=rd.read_count;
			if(leftsup[i-1]>=junctionsupport && rightsup[nex-i-1]>=junctionsupport) rd.juncs[i-1]->nreads_good+=rd.read_count;
		}
	}

}

/*
void clean_junctions(GList<CJunction>& junction, int refstart, GVec<float>& bpcov,GPVec<GffObj>& guides) {

	//fprintf(stderr,"Clean junctions:\n");
	for(int i=0;i<junction.Count();i++) {
		CJunction& jd=*(junction[i]);
		//if(jd.nreads_good<junctionthr && (!guides.Count() || guideintrons.IndexOf(jd)==-1)) {
		if (jd.nreads_good<junctionthr && (!guides.Count() || !jd.guide_match)) {
			//fprintf(stderr,"deleted junction: %d-%d (%d)\n",jd.start,jd.end,jd.strand);
			jd.strand=0;
		}
		else //if((int)(jd.end-jd.start)>longintron && (!guides.Count() || guideintrons.IndexOf(jd)==-1)) {
			// very long intron -> hard to trust unless it's well covered
		 if ((int)(jd.end-jd.start)>longintron && (!guides.Count() || !jd.guide_match)) {
			int leftreach = jd.start-longintronanchor-refstart;
			if(leftreach<0) leftreach=0;
			int rightreach = jd.end+longintronanchor-refstart;
			if(rightreach>=bpcov.Count()) rightreach=bpcov.Count()-1;
			//fprintf(stderr,"cov[%d]=%f cov[%d]=%f cov[%d][=%f cov[%d]=%f\n",leftreach+refstart,bpcov[leftreach],jd.start,bpcov[jd.start-refstart-1],rightreach+refstart,bpcov[rightreach],jd.end,bpcov[jd.end-refstart]);
			if((bpcov[leftreach]<1 && bpcov[jd.start-refstart-1]<2) ||
					(bpcov[rightreach]<1 && bpcov[jd.end-refstart]<2)) {
				jd.strand=0;
			}
		}
		//{ // DEBUG ONLY
		//	fprintf(stderr,"Junction %d: %d %d %d %g %g\n",i,jd.start,jd.end,jd.strand,jd.nreads,jd.nreads_good);
		//}
	} //for each read junction
}
*/

//debug funcs
/*
void showReads(GStr& refname, GList<CReadAln>& readlist) {
	if (readlist.Count()==0) return;
	GStr fname(outfname);
	fname+=".reads";
	FILE* f=fopen(fname.chars(),"a");
	for (int i=0;i<readlist.Count();i++) {
		CReadAln & rd=*(readlist[i]);
		fprintf(f,"%s %s %d %d %d [", rd.name.chars(), refname.chars(), rd.start, rd.end, (int)rd.strand);
		for (int j=0;j<rd.juncs.Count();j++) {
			CJunction& jd=*(rd.juncs[j]);
			if (j) fprintf(f, ",");
			fprintf(f, "%d-%d:%4.2f", jd.start, jd.end, jd.nreads_good);
		}
		fprintf(f,"]\n");
	}
	fclose(f);
}
*/

//int infer_transcripts(int refstart, GList<CReadAln>& readlist,
		//GList<CJunction>& junction, GPVec<GffObj>& guides, GVec<float>& bpcov, GList<CPrediction>& pred, bool fast) {
int infer_transcripts(BundleData* bundle, bool fast) {
	int geneno=0;

	//DEBUG ONLY: 	showReads(refname, readlist);

/*
#ifdef GMEMTRACE
	double vm,rsm;
	get_mem_usage(vm, rsm);
	GMessage("\t\tM(s):infer_transcripts memory usage: rsm=%6.1fMB vm=%6.1fMB\n",rsm/1024,vm/1024);
#endif
*/

	if(bundle->keepguides.Count() || !eonly) {
		count_good_junctions(bundle->readlist, bundle->start, bundle->bpcov);
		geneno = build_graphs(bundle, fast);
	}

/*
#ifdef GMEMTRACE
	//double vm,rsm;
	get_mem_usage(vm, rsm);
	GMessage("\t\tM(e):infer_transcripts memory usage: rsm=%6.1fMB vm=%6.1fMB\n",rsm/1024,vm/1024);
#endif
*/

	return(geneno);
}

void printGff3Header(FILE* f, GArgs& args) {
  fprintf(f, "# ");
  args.printCmdLine(f);
  fprintf(f, "##gff-version 3\n");
}


int predCmp(const pointer p1, const pointer p2) {
	CPrediction *a=(CPrediction*)p1;
	CPrediction *b=(CPrediction*)p2;
	if(a->start < b->start) return -1;
	if(a->start > b->start) return 1;
	//if(sensitivitylevel!=1 || sensitivitylevel!=2) { // try to see if I correct this if it makes any difference (it should be && instead of || in the if)
		if(a->exons.Count() < b->exons.Count()) return- 1;
		if(a->exons.Count() > b->exons.Count()) return 1;
		int i=0;
		uint a1=0;
		uint b1=0;
		while(i<a->exons.Count()) {
			if(a->exons[i].start<b->exons[i].start) {
				a1=a->exons[i].start;
				b1=b->exons[i].start;
				break;
			}
			if(a->exons[i].end<b->exons[i].end) {
				a1=a->exons[i].end;
				b1=b->exons[i].end;
				break;
			}
			i++;
		}
		if(a1) {
			if(a1<b1) return -1;
			if(a1>b1) return 1;
		}
	//}
	return 0;
}

int predexCmp(const pointer p1, const pointer p2) {
	CPrediction *a=(CPrediction*)p1;
	CPrediction *b=(CPrediction*)p2;
	if(a->exons.Count() < b->exons.Count()) return -1;
	if(a->exons.Count() > b->exons.Count()) return 1;
	if(a->tlen < b->tlen) return -1;
	if(a->tlen > b->tlen) return 1;
	return 0;
}

int predcovCmp(const pointer p1, const pointer p2) {
	CPrediction *a=(CPrediction*)p1;
	CPrediction *b=(CPrediction*)p2;
	if(a->cov < b->cov) return 1;
	if(a->cov > b->cov) return -1;
	return 0;
}

bool equal_pred(GList<CPrediction>& pred,int n1,int n2){

	if(pred[n1]->strand!=pred[n2]->strand) return false;

	//if(pred[n1]->start!=pred[n2]->start) return(false); // this allows genes with different start/ends to be merged together
	//if(pred[n1]->end!=pred[n2]->end) return(false);     // but I need to check if they overlap in case of single exons

	if((pred[n1]->end < pred[n2]->start) || (pred[n2]->end<pred[n1]->start)) return(false); // genes don't overlap

	int nex=pred[n1]->exons.Count();
	if(nex!=pred[n2]->exons.Count()) return(false);
	for(int i=0;i<nex;i++) {
		//if(pred[n1]->exons[i].start!=pred[n2]->exons[i].start) { fprintf(stderr,"ret false start[%d]: %d vs %d\n",i,pred[n1]->exons[i].start,pred[n2]->exons[i].start); return(false);}
		//if(pred[n1]->exons[i].end!=pred[n2]->exons[i].end) { fprintf(stderr,"ret false end[%d]: %d vs %d\n",i,pred[n1]->exons[i].end,pred[n2]->exons[i].end); return(false);}
		if(i>0 && (pred[n1]->exons[i].start!=pred[n2]->exons[i].start)) return(false);
		if(i<nex-1 && (pred[n1]->exons[i].end!=pred[n2]->exons[i].end)) return(false);
	}

	return(true);
}


CInterval *add_pred_to_cov(CInterval *maxcov, CPrediction* pred, bool *abundant=NULL) { // maybe I can eliminate some genes here


	//fprintf(stderr,"add pred: %d-%d with cov=%f to maxpos\n",pred->start,pred->end,pred->cov);

	if(maxcov==NULL || pred->end<maxcov->pos) { // prediction before current intervals or no current intervals
		uint linkstart=pred->end+1;
		CInterval *link=maxcov;
		if(maxcov && linkstart==maxcov->pos) { // next interval starts immediately after this one ends
			maxcov = new CInterval(pred->start,pred->cov,link);
		}
		else {
			CInterval *interval=new CInterval(linkstart,0,maxcov);
			maxcov = new CInterval(pred->start,pred->cov,interval);
		}
	}
	else { // I need to place current prediction
		CInterval *lastinterv=NULL;
		CInterval *interv=maxcov;
		while(interv && pred->start>=interv->pos) {
			lastinterv=interv;
			interv=interv->next;
		}
		if(interv){ // pred->start < interv->pos
			if(lastinterv) { // pred->start >= lastinterv->pos

				if(pred->end<interv->pos) {
					if(lastinterv->val<pred->cov) { // only in this case I am interested to do something, otherwise I am done
						if(pred->end+1<interv->pos) {
							lastinterv->next=new CInterval(pred->end+1,lastinterv->val,interv);
						}
						if(lastinterv->pos==pred->start) lastinterv->val=pred->cov;
						else {
							interv=lastinterv->next;
							lastinterv->next=new CInterval(pred->start,pred->cov,interv);
						}
					}
					else if(abundant && pred->cov<isofrac*lastinterv->val) *abundant=false;
				}
				else { // pred->end>=interv->pos; but pred->start >= lastinterv->pos
					float lastcov=lastinterv->val; // lastcov might be 0 here
					if(pred->cov>lastinterv->val) { // create new interval unless the position is in the same place as lastinterv
						if(pred->start==lastinterv->pos) lastinterv->val=pred->cov;
						else {
							lastinterv->next=new CInterval(pred->start,pred->cov,interv);
							lastinterv=lastinterv->next;
							lastcov=lastinterv->val;
						}
					}
					else if(abundant && pred->cov<isofrac*lastinterv->val) *abundant=false;
					while(interv && pred->end>=interv->pos) {
						if(interv->val<pred->cov) {
							if(lastinterv->val==pred->cov) { // I have to skip current interval
								lastinterv->next=interv->next;
								lastcov=interv->val;
								//free(interv);
								delete interv;
								interv=lastinterv;
							}
							else { lastcov=interv->val; interv->val=pred->cov;lastinterv=interv;}
						}
						else {
							if(abundant && pred->cov<isofrac*interv->val) *abundant=false;
							lastinterv=interv;
							lastcov=lastinterv->val;
						}
						interv=interv->next;
					}
					if(interv) { // pred->end<interv->pos: I might need to create new interval
						if(lastinterv->val<pred->cov) lastinterv->val=pred->cov;
						else if(abundant && pred->cov<isofrac*lastinterv->val) *abundant=false;
						uint newstart=pred->end+1;
						if(newstart<interv->pos) {
							lastinterv->next=new CInterval(newstart,lastcov,interv);
						}
					}
					else { // pred->end >= lastinterv->pos
						if(lastinterv->val<pred->cov) lastinterv->val=pred->cov;
						else if(abundant && pred->cov<isofrac*lastinterv->val) *abundant=false;
						lastinterv->next=new CInterval(pred->end+1);
					}
				}
			}
			else { // lastinterv == NULL
				maxcov=new CInterval(pred->start,pred->cov,interv);
				if(pred->end<interv->pos) {
					if(pred->end+1<interv->pos) maxcov->next=new CInterval(pred->end+1,0,interv);
				}
				else { // pred->end >= interv->pos
					lastinterv=maxcov;
					while(interv && pred->end>=interv->pos) {
						if(interv->val<pred->cov) {
							if(lastinterv->val==pred->cov) { // I have to skip current interval
								lastinterv->next=interv->next;
								//free(interv);
								delete interv; interv=lastinterv;
							}
							else { interv->val=pred->cov;lastinterv=interv;}
						}
						else {
							if(abundant && pred->cov<isofrac*interv->val) *abundant=false;
							lastinterv=interv;
						}
						interv=interv->next;
					}
					if(interv) { // pred->end<interv->pos I might need to create new interval
						float lastcov=lastinterv->val;
						if(lastinterv->val<pred->cov) lastinterv->val=pred->cov;
						else if(abundant && pred->cov<isofrac*lastinterv->val) *abundant=false;
						uint newstart=pred->end+1;
						if(newstart<interv->pos) lastinterv->next=new CInterval(newstart,lastcov,interv);
					}
					else {
						if(lastinterv->val<pred->cov) lastinterv->val=pred->cov;
						else if(abundant && pred->cov<isofrac*lastinterv->val) *abundant=false;
						lastinterv->next=new CInterval(pred->end+1);
					}
				}
			}
		}
		else { // prediction is after the end of maxcov, or is equal to lastinterv->pos
			interv=new CInterval(pred->end+1);
			if(pred->start==lastinterv->pos) {
				if(pred->cov>=lastinterv->val) {
					lastinterv->val=pred->cov;
					lastinterv->next=interv;
				}
				else {
					if(abundant && pred->cov<isofrac*lastinterv->val) *abundant=false;
					lastinterv->next= new CInterval(pred->start+1,pred->cov,interv);
				}
			}
			else lastinterv->next = new CInterval(pred->start,pred->cov,interv);
		}
	}

	return(maxcov);
}

bool is_pred_above_frac(CInterval *maxcov,CPrediction* pred) {

	CInterval *lastinterv=NULL;
	while(maxcov && pred->start>=maxcov->pos) {
		lastinterv=maxcov;
		maxcov=maxcov->next;
	}
	if(lastinterv && ((pred->exons.Count()==1 && pred->cov<lastinterv->val) || pred->cov<isofrac*lastinterv->val)) return(false); // I need to deal with single exons too here
	while(maxcov && pred->end>=maxcov->pos) {
		if((pred->exons.Count()==1 && pred->cov<maxcov->val) || pred->cov<isofrac*maxcov->val) return(false);
		maxcov=maxcov->next;
	}
	return(true);
}

void delete_interval(CInterval *interv){
	if(interv) {
		if(interv->next) delete_interval(interv->next);
		delete interv;
	}
}


/*
// I don't use this
int print_cluster(GList<CPrediction>& pred,GVec<int>& genes,GVec<int>& transcripts,
		int nstart, int nend, int geneno,GStr& refname) {

  GVec<int> keep;

  CInterval *maxposcov=NULL; //remembers intervals of maximum positive coverage
  CInterval *maxnegcov=NULL; //remembers intervals of maximum positive coverage

  bool pos;
  bool neg;
  int lastadded=0;
  for(int n=nstart;n<=nend;n++) {
	  if(pred[n]->strand=='+') { pos=true;neg=false;}
	  else if(pred[n]->strand=='-') { pos=false;neg=true;}
	  else { pos=true;neg=true;}
	  if(n>nstart) {
		  if(equal_pred(pred,lastadded,n)) {
			  pred[lastadded]->cov+=pred[n]->cov;
			  if(pos) maxposcov=add_pred_to_cov(maxposcov,pred[lastadded]);
			  if(neg) maxnegcov=add_pred_to_cov(maxnegcov,pred[lastadded]);
			  continue;
		  }
	  }
	  if(pos) maxposcov=add_pred_to_cov(maxposcov,pred[n]);
	  if(neg) maxnegcov=add_pred_to_cov(maxnegcov,pred[n]);
	  keep.Add(n);
	  lastadded=n;
  }

  for(int i=0;i<keep.Count();i++) {
	  int n=keep[i];
	  if(pred[n]->strand=='+') { pos=true;neg=false;}
	  else if(pred[n]->strand=='-') { pos=false;neg=true;}
	  else { pos=false;neg=false;}
	  if((pos && is_pred_above_frac(maxposcov,pred[n])) || (neg && is_pred_above_frac(maxnegcov,pred[n])) ||
			  (!pos && !neg && is_pred_above_frac(maxnegcov,pred[n]))) {
		  if(genes[pred[n]->geneno]==-1) genes[pred[n]->geneno]=++geneno;
		  transcripts[pred[n]->geneno]++;
		  fprintf(f_out,"%s\tStringTie\ttranscript\t%d\t%d\t1000\t%c\t.\tgene_id \"%s.%d\"; transcript_id \"%s.%d.%d\"; cov \"%.6f\";\n",
				  refname.chars(),pred[n]->start,pred[n]->end,pred[n]->strand,label.chars(),genes[pred[n]->geneno],
				  label.chars(),genes[pred[n]->geneno],transcripts[pred[n]->geneno],pred[n]->cov);
		  for(int j=0;j<pred[n]->exons.Count();j++)
			  fprintf(f_out,"%s\tStringTie\texon\t%d\t%d\t1000\t%c\t.\tgene_id \"%s.%d\"; transcript_id \"%s.%d.%d\"; exon_number \"%d\"; cov \"%.6f\";\n",
			  		 refname.chars(),pred[n]->exons[j].start,pred[n]->exons[j].end,pred[n]->strand,label.chars(),genes[pred[n]->geneno],
			  		 label.chars(),genes[pred[n]->geneno],transcripts[pred[n]->geneno],j+1,pred[n]->cov); // maybe add exon coverage here
	  }
  }

  delete_interval(maxposcov);
  delete_interval(maxnegcov);

  return(geneno);
}
*/

/*
// I don't use this
int print_transcript_cluster(GList<CPrediction>& pred,GVec<int>& genes,GVec<int>& transcripts,
		int nstart, int nend, int geneno,GStr& refname) {

  GVec<int> keep;

  float maxcovpos=0;
  float maxcovneg=0;
  bool pos;
  bool neg;
  int lastadded=0;
  for(int n=nstart;n<=nend;n++) {
	  if(pred[n]->strand=='+') { pos=true;neg=false;}
	  else if(pred[n]->strand=='-') { pos=false;neg=true;}
	  else { pos=true;neg=true;}
	  if(n>nstart) {
		  if(equal_pred(pred,lastadded,n)) {
			  pred[lastadded]->cov+=pred[n]->cov;
			  if(pos && pred[lastadded]->cov > maxcovpos) {
				  maxcovpos=pred[lastadded]->cov;
			  }
			  if(neg && pred[lastadded]->cov > maxcovneg) {
				  maxcovneg=pred[lastadded]->cov;
			  }
			  continue;
		  }
	  }
	  if(pos && pred[n]->cov>maxcovpos) {
		  maxcovpos=pred[n]->cov;
	  }
	  if(neg && pred[n]->cov>maxcovneg) {
		  maxcovneg=pred[n]->cov;
	  }
	  keep.Add(n);
	  lastadded=n;
  }

  for(int i=0;i<keep.Count();i++) {
	  int n=keep[i];
	  if(pred[n]->strand=='+') { pos=true;neg=false;}
	  else if(pred[n]->strand=='-') { pos=false;neg=true;}
	  else { pos=false;neg=false;}
	  if((pos && ((pred[n]->exons.Count()==1 && pred[n]->cov>=maxcovpos) || (pred[n]->exons.Count()>1 && pred[n]->cov/maxcovpos>=isofrac))) ||
			  (neg && ((pred[n]->exons.Count()==1 && pred[n]->cov>=maxcovneg) || (pred[n]->exons.Count()>1 && pred[n]->cov/maxcovneg>=isofrac))) ||
			  (!pos && !neg && ((pred[n]->exons.Count()==1 && pred[n]->cov>=maxcovneg && pred[n]->cov>=maxcovpos) ||
					  (pred[n]->exons.Count()>1 && pred[n]->cov/maxcovpos>=isofrac && pred[n]->cov/maxcovneg>=isofrac)))) { // print this transcript
		  if(genes[pred[n]->geneno]==-1) genes[pred[n]->geneno]=++geneno;
		  transcripts[pred[n]->geneno]++;
		  fprintf(f_out,"%s\tStringTie\ttranscript\t%d\t%d\t1000\t%c\t.\tgene_id \"%s.%d\"; transcript_id \"%s.%d.%d\"; cov \"%.6f\";\n",
				  refname.chars(),pred[n]->start,pred[n]->end,pred[n]->strand,label.chars(),genes[pred[n]->geneno],
				  label.chars(),genes[pred[n]->geneno],transcripts[pred[n]->geneno],pred[n]->cov);
		  for(int j=0;j<pred[n]->exons.Count();j++)
			  fprintf(f_out,"%s\tStringTie\texon\t%d\t%d\t1000\t%c\t.\tgene_id \"%s.%d\"; transcript_id \"%s.%d.%d\"; exon_number \"%d\"; cov \"%.6f\";\n",
			  		 refname.chars(),pred[n]->exons[j].start,pred[n]->exons[j].end,pred[n]->strand,label.chars(),genes[pred[n]->geneno],
			  		 label.chars(),genes[pred[n]->geneno],transcripts[pred[n]->geneno],j+1,pred[n]->cov); // maybe add exon coverage here
	  }
  }

  return(geneno);
}
*/

int print_signcluster(char strand,GList<CPrediction>& pred,GVec<int>& genes,GVec<int>& transcripts,
		int nstart, int nend, int geneno,GStr& refname) {

  GVec<int> keep;

  CInterval *maxpos=NULL; //remembers intervals of maximum coverage

  int lastadded=0;
  for(int n=nstart;n<=nend;n++) if(pred[n]->strand==strand || pred[n]->strand=='.'){
	  if(n>nstart) {
		  if(equal_pred(pred,lastadded,n)) {
			  if(pred[n]->cov>pred[lastadded]->cov){
				  pred[lastadded]->flag=pred[n]->flag;
				  if(pred[lastadded]->exons[0].start != pred[n]->exons[0].start ||
						  pred[lastadded]->exons.Last().end!=pred[n]->exons.Last().end) { // new prediction has to replace old one
					  pred[lastadded]->tlen=pred[n]->tlen;
					  pred[lastadded]->exons[0].start=pred[n]->exons[0].start;
					  pred[lastadded]->exons.Last().end=pred[n]->exons.Last().end;
				  }
			  }
			  for(int j=0;j<pred[n]->exons.Count();j++) {
				  pred[lastadded]->exoncov[j]+=pred[n]->exoncov[j];
			  }
			  pred[lastadded]->cov+=pred[n]->cov;
			  //pred[lastadded]->frag+=pred[n]->frag;

			  maxpos=add_pred_to_cov(maxpos,pred[lastadded]);
			  if(pred[n]->t_eq && !pred[lastadded]->t_eq) { pred[lastadded]->t_eq=pred[n]->t_eq;}
			  //if(pred[n]->id && !pred[lastadded]->id) { pred[lastadded]->id=Gstrdup(pred[n]->id);}
			  continue;
		  }
	  }
	  bool abundant=true;
	  maxpos=add_pred_to_cov(maxpos,pred[n],&abundant);
	  ////maxpos=add_pred_to_cov(maxpos,pred[n]);
	  if(pred[n]->t_eq || abundant) {
	  //if(pred[n]->id || abundant) {
		  keep.Add(n);
		  lastadded=n;
	  }
  }

  for(int i=0;i<keep.Count();i++) {
	  int n=keep[i];
	  if(pred[n]->t_eq || is_pred_above_frac(maxpos,pred[n])) { // print this transcript
		  if(genes[pred[n]->geneno]==-1) genes[pred[n]->geneno]=++geneno;
		  transcripts[pred[n]->geneno]++;
		  if(pred[n]->flag) {

			  uint t_id=0;
			  if (pred[n]->t_eq && pred[n]->t_eq->uptr) {
				  t_id = ((RC_TData*)pred[n]->t_eq->uptr)->t_id;
			  }
			  //fprintf(f_out,"%d %d %d %.6f %.6f\n",pred[n]->exons.Count()+1,pred[n]->tlen, t_id, pred[n]->frag,pred[n]->cov);
			  fprintf(f_out,"%d %d %d %.6f\n",pred[n]->exons.Count()+1,pred[n]->tlen, t_id, pred[n]->cov);
			  fprintf(f_out,"%s\tStringTie\ttranscript\t%d\t%d\t1000\t%c\t.\tgene_id \"%s.%d\"; transcript_id \"%s.%d.%d\"; ",
					  refname.chars(),pred[n]->start,pred[n]->end,pred[n]->strand,label.chars(),genes[pred[n]->geneno],
					  label.chars(),genes[pred[n]->geneno],transcripts[pred[n]->geneno]);
			  if(pred[n]->t_eq) {
				  fprintf(f_out,"reference_id \"%s\"; ",pred[n]->t_eq->getID());
				  if (pred[n]->t_eq->getGeneID())
				    	 fprintf(f_out,"ref_gene_id \"%s\"; ",pred[n]->t_eq->getGeneID());
				  if (pred[n]->t_eq->getGeneName())
				    	 fprintf(f_out,"ref_gene_name \"%s\"; ",pred[n]->t_eq->getGeneName());
			  }
			  fprintf(f_out,"cov \"%.6f\";\n",pred[n]->cov);
			  for(int j=0;j<pred[n]->exons.Count();j++) {
				  fprintf(f_out,"%s\tStringTie\texon\t%d\t%d\t1000\t%c\t.\tgene_id \"%s.%d\"; transcript_id \"%s.%d.%d\"; exon_number \"%d\"; ",
						  refname.chars(),pred[n]->exons[j].start,pred[n]->exons[j].end,pred[n]->strand,label.chars(),genes[pred[n]->geneno],
						  label.chars(),genes[pred[n]->geneno],transcripts[pred[n]->geneno],j+1); // maybe add exon coverage here
				  if(pred[n]->t_eq) {
					  fprintf(f_out,"reference_id \"%s\"; ",pred[n]->t_eq->getID());
					  if (pred[n]->t_eq->getGeneID())
					    	 fprintf(f_out,"ref_gene_id \"%s\"; ",pred[n]->t_eq->getGeneID());
					  if (pred[n]->t_eq->getGeneName())
					    	 fprintf(f_out,"ref_gene_name \"%s\"; ",pred[n]->t_eq->getGeneName());
				  }
				  fprintf(f_out,"cov \"%.6f\";\n",pred[n]->exoncov[j]);
			  }
		  }
		  else pred[n]->flag=true;
	  }
	  else pred[n]->flag=false;
  }

  delete_interval(maxpos);

  return(geneno);
}

uint min(uint n1,uint n2) {
	if(n1<n2) return(n1);
	return(n2);
}

uint max(uint n1,uint n2) {
	if(n1<n2) return(n2);
	return n1;
}

bool included_pred(GPVec<CPrediction>& pred,int n1,int n2) { // check if the small prediction is included in the larger prediction

	if(pred[n1]->start > pred[n2]->end || pred[n2]->start>pred[n1]->end) return false;

	int big=n1;
	int small=n2;
	if(pred[n1]->exons.Count()<pred[n2]->exons.Count()) {
		big=n2;
		small=n1;
	}

	if(pred[small]->t_eq) { // this is a reference guide
		if(pred[big]->t_eq  && pred[small]->t_eq!=pred[big]->t_eq) return false;
		if(pred[small]->exons.Count()!=pred[big]->exons.Count()) return false;
	}

	/*
	if(pred[big]->t_eq && !pred[small]->t_eq && pred[small]->exons.Count()==1 &&
			pred[small]->start>=pred[big]->start && pred[small]->end<=pred[big]->end) { // bigger prediction is a reference transcript which includes a smaller single exon gene
		int bex=0;
		while(bex<pred[big]->exons.Count() && pred[small]->start>pred[big]->exons[bex].end) {
			bex++;
		}
		if(bex<pred[big]->exons.Count() && pred[small]->start>=pred[big]->exons[bex].start
				&& pred[small]->end<=pred[big]->exons[bex].end) return true;
		return false;
	}
	*/

	int bex=0;
	while(bex<pred[big]->exons.Count()) {
		if(pred[small]->exons[0].start>pred[big]->exons[bex].end) bex++;
		else { // now pred[small]->exons[0].start <= pred[big]->exons[bex].end
			if(pred[small]->exons[0].end<pred[big]->exons[bex].start) return false; // no overlap
			int sex=0;
			while(sex<pred[small]->exons.Count() && bex<pred[big]->exons.Count()) {
				if(sex==pred[small]->exons.Count()-1) { // I am at end of small pred
					if(bex==pred[big]->exons.Count()-1) return true; // same intron structure and there is overlap
					// here not at last exon in big prediction
					if(pred[small]->exons[sex].end>=pred[big]->exons[bex+1].start) return false;
					return true;
				}
				// sex is not last exon in small prediction but overlaps bex
				if(bex==pred[big]->exons.Count()-1) return false; // small pred extends past big pred
				if(pred[small]->exons[sex].end != pred[big]->exons[bex].end) return false;
				bex++;
				sex++;
				if(pred[small]->exons[sex].start != pred[big]->exons[bex].start) return false;
			}
			return false;
		}
	}

	return false;
}

void update_cov(GPVec<CPrediction>& pred,int big,int small,float frac=1) { // small gets included into big

	if(pred[big]->strand=='.') pred[big]->strand=pred[small]->strand;

	if(pred[small]->t_eq && !pred[big]->t_eq) {
		pred[big]->tlen=pred[small]->tlen;
		pred[big]->exons[0].start=pred[small]->exons[0].start;
		pred[big]->start=pred[small]->start;
		pred[big]->exons.Last().end=pred[small]->exons.Last().end;
		pred[big]->end=pred[small]->end;
	}

	int bex=0;
	while(pred[small]->exons[0].start>pred[big]->exons[bex].end) bex++;

	if(!pred[big]->t_eq && !bex && pred[small]->exons.Count()>1 && pred[big]->exons[0].start<pred[small]->exons[0].start) { // adjust start to account for trimming
		pred[big]->tlen-=pred[small]->exons[0].start-pred[big]->exons[0].start;
		pred[big]->exons[0].start=pred[small]->exons[0].start;
		pred[big]->start=pred[small]->start;
	}

	int sex=0;
	int overlap=0;
	while(sex<pred[small]->exons.Count()) {
		int exovlp=(pred[small]->exons[sex].end<pred[big]->exons[bex].end ? pred[small]->exons[sex].end : pred[big]->exons[bex].end)-
				(pred[small]->exons[sex].start>pred[big]->exons[bex].start ? pred[small]->exons[sex].start : pred[big]->exons[bex].start)+1;

		if(!pred[big]->t_eq && bex==pred[big]->exons.Count()-1 && sex>=1 && pred[big]->exons[bex].end>pred[small]->exons[sex].end) { // adjust end
			pred[big]->tlen-=pred[big]->exons[bex].end-pred[small]->exons[sex].end;
			pred[big]->exons[bex].end=pred[small]->exons[sex].end;
			pred[big]->end=pred[small]->end;
		}

		pred[big]->exoncov[bex]=(pred[big]->exoncov[bex]*pred[big]->exons[bex].len()+frac*pred[small]->exoncov[sex]*exovlp)/pred[big]->exons[bex].len();
		overlap+=exovlp;
		sex++;bex++;
	}

	pred[big]->cov=(pred[big]->tlen*pred[big]->cov+overlap*frac*pred[small]->cov)/pred[big]->tlen;

}

void merge_exons(CGene& gene,GList<GffExon>& exons) {
	int ig=0;
	int ie=0;
	while(ie<exons.Count()) {
		if(ig==gene.exons.Count() || exons[ie]->end<gene.exons[ig].start) {
			GSeg ex(exons[ie]->start,exons[ie]->end);
			gene.exons.Insert(ig,ex);
			ie++;
			ig++;
			continue;
		}
		while(ig<gene.exons.Count() && exons[ie]->start>gene.exons[ig].end) ig++;
		if(ig<gene.exons.Count()) { // here exons[ie]->start<=gene.exons[ig].end and exons[ie]->end>=gene.exons[ig].start
			if(exons[ie]->start<=gene.exons[ig].start) gene.exons[ig].start=exons[ie]->start;
			if(exons[ie]->end>=gene.exons[ig].end) {
				gene.exons[ig].end=exons[ie]->end;
				ig++;
				while(ig<gene.exons.Count() && exons[ie]->end>=gene.exons[ig].start) {
					if(gene.exons[ig].end>exons[ie]->end) gene.exons[ig-1].end=gene.exons[ig].end;
					gene.exons.Delete(ig);
				}
			}
			ie++;
		}
	}
}

void merge_exons(CGene& gene,GVec<GSeg>& exons) {
	int ig=0;
	int ie=0;
	while(ie<exons.Count()) {
		if(ig==gene.exons.Count() || exons[ie].end<gene.exons[ig].start) {
			gene.exons.Insert(ig,exons[ie]);
			ie++;
			ig++;
			continue;
		}
		while(ig<gene.exons.Count() && exons[ie].start>gene.exons[ig].end) ig++;
		if(ig<gene.exons.Count()) { // here exons[ie]->start<=gene.exons[ig].end and exons[ie]->end>=gene.exons[ig].start
			if(exons[ie].start<=gene.exons[ig].start) gene.exons[ig].start=exons[ie].start;
			if(exons[ie].end>=gene.exons[ig].end) {
				gene.exons[ig].end=exons[ie].end;
				ig++;
				while(ig<gene.exons.Count() && exons[ie].end>=gene.exons[ig].start) {
					if(gene.exons[ig].end>exons[ie].end) gene.exons[ig-1].end=gene.exons[ig].end;
					gene.exons.Delete(ig);
				}
			}
			ie++;
		}
	}
}


// default printing function for sensitivitylevel=1: all the others are deprecated for now
int print_cluster(GPVec<CPrediction>& pred,GVec<int>& genes,GVec<int>& transcripts, int geneno,GStr& refname,
		GVec<CGene>& refgene, GHash<int>& hashgene, GVec<CGene>& predgene, int startgno) {

	// sort predictions from the most abundant to the least:
	pred.Sort(predcovCmp);
	GVec<int> keep;

	CInterval *maxpos=NULL; //remembers intervals of maximum coverage

	for(int n=0;n<pred.Count();n++) {

		/*
		{ // DEBUG ONLY
			fprintf(stderr,"Consider prediction[%d] %c cov=%f:",n,pred[n]->strand,pred[n]->cov);
			for(int i=0;i<pred[n]->exons.Count();i++) fprintf(stderr," %d-%d",pred[n]->exons[i].start,pred[n]->exons[i].end);
			if(pred[n]->t_eq) fprintf(stderr," ref_id=%s",pred[n]->t_eq->getID());
			fprintf(stderr,"\n");
		}
		*/

		int k=0;
		bool included=false;
		while(!included && k<keep.Count() && keep[k]<n) {

			if(included_pred(pred,keep[k],n)) {

				//fprintf(stderr,"%d prediction is included in prediction %d\n",keep[k],n);

				bool checkall=false;

				if(pred[keep[k]]->exons.Count()<pred[n]->exons.Count()) {
					//if(pred[keep[k]]->cov>pred[n]->cov) break; // this is new and improves specificity but I loose some things -> TO CHECK WHAT IT ACT: also this should always happen because of the sort procedure
					//if(pred[keep[k]]->exons.Count()>2) break; // this is what I had before but I don't think it makes sense completely so I introduce the next one for those cases where some single exons are still included in here
					//*** if(pred[keep[k]]->exons.Count()>2) { k++; break;} // I need to check this how it affects performance in general
					if(pred[keep[k]]->exons.Count()>2) { k++; if(pred[n]->t_eq) continue; else break;} // I need to check this how it affects performance in general
					update_cov(pred,n,keep[k]);
					pred[keep[k]]->cov=pred[n]->cov;
					pred[keep[k]]->exons.Clear();
					pred[keep[k]]->exons.Add(pred[n]->exons);
					pred[keep[k]]->exoncov.Clear();
					pred[keep[k]]->exoncov.Add(pred[n]->exoncov);
					pred[keep[k]]->flag=true;
					if(!pred[keep[k]]->t_eq) pred[keep[k]]->start=pred[n]->start; // only adjust start if it's not already known
					if(!pred[keep[k]]->t_eq) pred[keep[k]]->end=pred[n]->end; // only adjust end if it's not already known
					pred[keep[k]]->tlen=pred[n]->tlen;
					if(pred[n]->t_eq) checkall=true;
				}
				else update_cov(pred,keep[k],n);

				if(pred[keep[k]]->strand=='.' && pred[n]->strand!='.') pred[keep[k]]->strand=pred[n]->strand;

				//if(pred[n]->id && !pred[keep[k]]->id) pred[keep[k]]->id=Gstrdup(pred[n]->id);
				if(pred[n]->t_eq && !pred[keep[k]]->t_eq) pred[keep[k]]->t_eq=pred[n]->t_eq;
				//pred[keep[k]]->frag+=pred[n]->frag;


				if(checkall) { // I need to test this too
					int j=k+1;
					while(j<keep.Count()){
						if(pred[keep[j]]->exons.Count()==1 && included_pred(pred,keep[k],keep[j])) { // if it's a single exon and it's included in keep[k], i might want to remove it because it is of higher value
							update_cov(pred,keep[k],keep[j]);
							if(pred[keep[j]]->t_eq && !pred[keep[k]]->t_eq) pred[keep[k]]->t_eq=pred[keep[j]]->t_eq;
							//pred[keep[k]]->frag+=pred[keep[j]]->frag;
							if(pred[keep[k]]->strand=='.' && pred[keep[j]]->strand!='.') pred[keep[k]]->strand=pred[keep[j]]->strand;
							keep.Delete(j);
						}
						else j++;
					}
				}


				//fprintf(stderr,"...included in prediction[%d] with cov=%f\n",keep[k],pred[keep[k]]->cov);
				included=true;
				break; // if it's included than I am done with the while loop because n got used
			}
			k++;
		}
		if(included) {
			maxpos=add_pred_to_cov(maxpos,pred[keep[k]]);

			/*
			{ // DEBUG ONLY
				fprintf(stderr,"Maxpos is:");
				CInterval *interval=maxpos;
				while(interval!=NULL) {
					fprintf(stderr," pos=%d val=%f",interval->pos,interval->val);
					interval=interval->next;
				}
				fprintf(stderr,"\n");
			}
			*/

			continue;
		}
		bool abundant=true;
		maxpos=add_pred_to_cov(maxpos,pred[n],&abundant);
		//maxpos=add_pred_to_cov(maxpos,pred[n]);

		/*
		{ // DEBUG ONLY
			fprintf(stderr,"Maxpos is:");
			CInterval *interval=maxpos;
			while(interval!=NULL) {
				fprintf(stderr," pos=%d val=%f",interval->pos,interval->val);
				interval=interval->next;
			}
			fprintf(stderr,"\n");
		}
		*/

		if(pred[n]->t_eq || abundant) {
		//if(pred[n]->id || abundant) {
			keep.Add(n);
			//fprintf(stderr,"...keep prediction %d\n",n);
		}
	}

  for(int i=0;i<keep.Count();i++) {
	  int n=keep[i];
	  if(pred[n]->t_eq || (!eonly && is_pred_above_frac(maxpos,pred[n]))) { // print this transcript

		  /*
		  { // DEBUG ONLY
			  fprintf(stderr,"print prediction %d",n);
			  if(pred[n]->flag) fprintf(stderr," with true flag");
			  fprintf(stderr,"\n");
		  }
		  */

		  if(pred[n]->flag) {
			  if(genes[pred[n]->geneno]==-1) genes[pred[n]->geneno]=++geneno;
			  transcripts[pred[n]->geneno]++;
			  uint t_id=0;
			  if (pred[n]->t_eq && pred[n]->t_eq->uptr) {
				  t_id = ((RC_TData*)pred[n]->t_eq->uptr)->t_id;
			  }
			  //fprintf(f_out,"%d %d %d %.6f %.6f\n",pred[n]->exons.Count()+1,pred[n]->tlen, t_id, pred[n]->frag,pred[n]->cov);
			  fprintf(f_out,"1 %d %d %d %.6f\n",pred[n]->exons.Count()+1,pred[n]->tlen, t_id,pred[n]->cov);
			  fprintf(f_out,"%s\tStringTie\ttranscript\t%d\t%d\t1000\t%c\t.\tgene_id \"%s.%d\"; transcript_id \"%s.%d.%d\"; ",
					  refname.chars(),pred[n]->start,pred[n]->end,pred[n]->strand,label.chars(),genes[pred[n]->geneno],
					  label.chars(),genes[pred[n]->geneno],transcripts[pred[n]->geneno]);
			  if(pred[n]->t_eq) {
				  fprintf(f_out,"reference_id \"%s\"; ",pred[n]->t_eq->getID());
				  if (pred[n]->t_eq->getGeneID())
				    	 fprintf(f_out,"ref_gene_id \"%s\"; ",pred[n]->t_eq->getGeneID());
				  if (pred[n]->t_eq->getGeneName())
				    	 fprintf(f_out,"ref_gene_name \"%s\"; ",pred[n]->t_eq->getGeneName());
			  }
			  fprintf(f_out,"cov \"%.6f\";\n",pred[n]->cov);
			  for(int j=0;j<pred[n]->exons.Count();j++) {
				  fprintf(f_out,"%s\tStringTie\texon\t%d\t%d\t1000\t%c\t.\tgene_id \"%s.%d\"; transcript_id \"%s.%d.%d\"; exon_number \"%d\"; ",
						  refname.chars(),pred[n]->exons[j].start,pred[n]->exons[j].end,pred[n]->strand,label.chars(),genes[pred[n]->geneno],
						  label.chars(),genes[pred[n]->geneno],transcripts[pred[n]->geneno],j+1); // maybe add exon coverage here
				  if(pred[n]->t_eq) {
					  fprintf(f_out,"reference_id \"%s\"; ",pred[n]->t_eq->getID());
					  if (pred[n]->t_eq->getGeneID())
					    	 fprintf(f_out,"ref_gene_id \"%s\"; ",pred[n]->t_eq->getGeneID());
					  if (pred[n]->t_eq->getGeneName())
					    	 fprintf(f_out,"ref_gene_name \"%s\"; ",pred[n]->t_eq->getGeneName());
				  }
				  fprintf(f_out,"cov \"%.6f\";\n",pred[n]->exoncov[j]);
			  }
			  pred[n]->flag=false;

			  // now deal with the genes
			  // predicted:
			  int gno=geneno-startgno;
			  if(!eonly) {
				  if(gno>=predgene.Count()) { // I did not see this gene before
					  CGene g(pred[n]->start,pred[n]->end,pred[n]->strand);
					  for(int i=0;i<pred[n]->exons.Count();i++) {
						  g.exons.Add(pred[n]->exons[i]);
					  }
					  g.cov+=pred[n]->cov*pred[n]->tlen;
					  g.covsum+=pred[n]->cov;
					  predgene.Add(g);
				  }
				  else { // I've seen this gene before
					  if(pred[n]->start<predgene[gno].start) predgene[gno].start=pred[n]->start;
					  if(pred[n]->end>predgene[gno].end) predgene[gno].end=pred[n]->end;
					  merge_exons(predgene[gno],pred[n]->exons);
					  predgene[gno].cov+=pred[n]->cov*pred[n]->tlen;
					  predgene[gno].covsum+=pred[n]->cov;
				  }
			  }
			  // annotated
			  if(pred[n]->t_eq && (geneabundance ||eonly)) {
				  GStr gid(pred[n]->t_eq->getGeneID());
				  const int *ng=hashgene[gid.chars()];
				  if(ng) { // this should always be true because we parsed all predictions in printResults
					  gno=*ng;
					  /* //I don't need to do this because I already did it in printResults
					  if(pred[n]->start<refgene[gno].start) refgene[gno].start=pred[n]->start;
					  if(pred[n]->end>predgene[gno].end) predgene[gno].end=pred[n]->end;
					  merge_exons(predgene[gno],pred[n]->exons);
					  */
					  refgene[gno].cov+=pred[n]->cov*pred[n]->tlen;
					  refgene[gno].covsum+=pred[n]->cov;
				  }
			  }
		  }
		  else pred[n]->flag=true;
	  }
	  else pred[n]->flag=false;
  }

  delete_interval(maxpos);

  return(geneno);
}


int print_cluster_inclusion(GPVec<CPrediction>& pred,GVec<int>& genes,GVec<int>& transcripts, int geneno,GStr& refname, int limit=3) {

	//fprintf(stderr,"start print cluster...\n");
	// sort predictions from the one with the most exons to the one with the least:
	pred.Sort(predexCmp);

	GVec< GVec<int> > included(pred.Count()); included.Resize(pred.Count());
	GVec<float> maxcov;
	GVec<float> totalcov;

	for(int n1=0;n1<pred.Count()-1;n1++) {
		float elem=0;
		maxcov.Add(elem);
		totalcov.Add(elem);
		bool equal=false;
		for(int n2=n1+1;n2<pred.Count();n2++)
			if(included_pred(pred,n1,n2)) {
				if(equal && pred[n1]->exons.Count()<pred[n2]->exons.Count()) break;
				if(pred[n1]->exons.Count()==pred[n2]->exons.Count()) {
					if(maxcov[n1]<pred[n1]->cov) maxcov[n1]=pred[n1]->cov;
					equal=true;
				}
				included[n1].Add(n2);
				if(pred[n2]->cov>maxcov[n1]) maxcov[n1]=pred[n2]->cov;
				totalcov[n1]+=pred[n2]->cov;
			}
	}

	GVec<int> keep;

	CInterval *maxpos=NULL; //remembers intervals of maximum coverage

	for(int n=0;n<pred.Count();n++) { // don't need this anymore since I already took care of it before: if(pred[n]->strand==strand || pred[n]->strand=='.'){

		/*
		{ // DEBUG ONLY
			fprintf(stderr,"Consider prediction[%d] %c cov=%f:",n,pred[n]->strand,pred[n]->cov);
			for(int i=0;i<pred[n]->exons.Count();i++) fprintf(stderr," %d-%d",pred[n]->exons[i].start,pred[n]->exons[i].end);
			fprintf(stderr," included in predictions:");
			for(int i=0;i<included[n].Count();i++) fprintf(stderr," %d",included[n][i]);
			if(n<pred.Count()-1) fprintf(stderr," maxcov=%f totalcov=%f",maxcov[n],totalcov[n]);
			fprintf(stderr,"\n");
		}
		*/

		if(included[n].Count() && (maxcov[n]>=pred[n]->cov || pred[n]->exons.Count()<limit)) { // this prediction is included in others, and is less abundant or has very few exons: make it <=2 if two exon genes also are to be ignored
			for(int k=0;k<included[n].Count();k++) {
				update_cov(pred,included[n][k],n,pred[included[n][k]]->cov/totalcov[n]);
			}
		}
		else {

			bool abundant=true;
			maxpos=add_pred_to_cov(maxpos,pred[n],&abundant);
			//maxpos=add_pred_to_cov(maxpos,pred[n]);

			/*
			{ // DEBUG ONLY
				fprintf(stderr,"Maxpos is:");
				CInterval *interval=maxpos;
				while(interval!=NULL) {
					fprintf(stderr," pos=%d val=%f",interval->pos,interval->val);
					interval=interval->next;
				}
				fprintf(stderr,"\n");
			}
			*/

			if(pred[n]->t_eq || abundant) {
			//if(pred[n]->id || abundant) {
				keep.Add(n);
				//fprintf(stderr,"...keep prediction\n");
			}
		}
	}

	for(int i=0;i<keep.Count();i++) {
	  int n=keep[i];
	  if(pred[n]->t_eq || (!eonly && is_pred_above_frac(maxpos,pred[n]))) { // print this transcript

		  //fprintf(stderr,"print prediction %d",n);
		  //if(pred[n]->flag) fprintf(stderr," with true flag");
		  //fprintf(stderr,"\n");

		  if(pred[n]->flag) {
			  if(genes[pred[n]->geneno]==-1) genes[pred[n]->geneno]=++geneno;
			  transcripts[pred[n]->geneno]++;
			  uint t_id=0;
			  if (pred[n]->t_eq && pred[n]->t_eq->uptr) {
				  t_id = ((RC_TData*)pred[n]->t_eq->uptr)->t_id;
			  }
			  //fprintf(f_out,"%d %d %d %.6f %.6f\n",pred[n]->exons.Count()+1,pred[n]->tlen, t_id, pred[n]->frag,pred[n]->cov);
			  fprintf(f_out,"%d %d %d %.6f\n",pred[n]->exons.Count()+1,pred[n]->tlen, t_id, pred[n]->cov);
			  fprintf(f_out,"%s\tStringTie\ttranscript\t%d\t%d\t1000\t%c\t.\tgene_id \"%s.%d\"; transcript_id \"%s.%d.%d\"; ",
					  refname.chars(),pred[n]->start,pred[n]->end,pred[n]->strand,label.chars(),genes[pred[n]->geneno],
					  label.chars(),genes[pred[n]->geneno],transcripts[pred[n]->geneno]);
			  if(pred[n]->t_eq) {
				  fprintf(f_out,"reference_id \"%s\"; ",pred[n]->t_eq->getID());
				  if (pred[n]->t_eq->getGeneID())
				    	 fprintf(f_out,"ref_gene_id \"%s\"; ",pred[n]->t_eq->getGeneID());
				  if (pred[n]->t_eq->getGeneName())
				    	 fprintf(f_out,"ref_gene_name \"%s\"; ",pred[n]->t_eq->getGeneName());
			  }
			  fprintf(f_out,"cov \"%.6f\";\n",pred[n]->cov);
			  for(int j=0;j<pred[n]->exons.Count();j++) {
				  fprintf(f_out,"%s\tStringTie\texon\t%d\t%d\t1000\t%c\t.\tgene_id \"%s.%d\"; transcript_id \"%s.%d.%d\"; exon_number \"%d\"; ",
						  refname.chars(),pred[n]->exons[j].start,pred[n]->exons[j].end,pred[n]->strand,label.chars(),genes[pred[n]->geneno],
						  label.chars(),genes[pred[n]->geneno],transcripts[pred[n]->geneno],j+1); // maybe add exon coverage here
				  if(pred[n]->t_eq) {
					  fprintf(f_out,"reference_id \"%s\"; ",pred[n]->t_eq->getID());
					  if (pred[n]->t_eq->getGeneID())
					    	 fprintf(f_out,"ref_gene_id \"%s\"; ",pred[n]->t_eq->getGeneID());
					  if (pred[n]->t_eq->getGeneName())
					    	 fprintf(f_out,"ref_gene_name \"%s\"; ",pred[n]->t_eq->getGeneName());
				  }
				  fprintf(f_out,"cov \"%.6f\";\n",pred[n]->exoncov[j]);
			  }
		  }
		  else pred[n]->flag=true;
	  }
	  else pred[n]->flag=false;
  }

  delete_interval(maxpos);

  return(geneno);
}


int print_transcript_signcluster(char strand,GList<CPrediction>& pred,GVec<int>& genes,GVec<int>& transcripts,
		int nstart, int nend, int geneno,GStr& refname) {

  GVec<int> keep;

  float maxcov=0;
  int lastadded=0;
  for(int n=nstart;n<=nend;n++) if(pred[n]->strand==strand || pred[n]->strand=='.'){
	  if(n>nstart) {
		  if(equal_pred(pred,lastadded,n)) {
			  /*
			  fprintf(stderr,"pred %d is equal to pred %d btw %d-%d\n",lastadded,n,pred[lastadded]->start,pred[lastadded]->end);
			  for(int j=0;j<pred[n]->exons.Count();j++) {
				  fprintf(stderr,"%d-%d vs %d-%d\n",pred[lastadded]->exons[j].start,pred[lastadded]->exons[j].end,pred[n]->exons[j].start,pred[n]->exons[j].end);
			  }
			  */
		    if(pred[n]->cov>pred[lastadded]->cov){
		      pred[lastadded]->flag=pred[n]->flag;
		      if(pred[lastadded]->exons[0].start != pred[n]->exons[0].start ||
		    		  pred[lastadded]->exons.Last().end!=pred[n]->exons.Last().end) { // new prediction has to replace old one
		    	  pred[lastadded]->tlen=pred[n]->tlen;
		    	  pred[lastadded]->exons[0].start=pred[n]->exons[0].start;
		    	  pred[lastadded]->exons.Last().end=pred[n]->exons.Last().end;

		      }
		    }
		    for(int j=0;j<pred[n]->exons.Count();j++) {
		      pred[lastadded]->exoncov[j]+=pred[n]->exoncov[j];
		    }
		    pred[lastadded]->cov+=pred[n]->cov;
		    //pred[lastadded]->frag+=pred[n]->frag;
		    if(pred[lastadded]->cov > maxcov) {
		      maxcov=pred[lastadded]->cov;
		    }
		    if(pred[n]->t_eq && !pred[lastadded]->t_eq) { pred[lastadded]->t_eq=pred[n]->t_eq;}
		    //if (pred[n]->id && !pred[lastadded]->id) { pred[lastadded]->id=Gstrdup(pred[n]->id);}
		    continue;
		  }
	  }
	  if(pred[n]->cov>maxcov) {
		  maxcov=pred[n]->cov;
	  }
	  keep.Add(n);
	  lastadded=n;
  }

  //fprintf(stderr,"keepcount=%d\n",keep.Count());

  for(int i=0;i<keep.Count();i++) {
	  int n=keep[i];
	  if(pred[n]->t_eq || (pred[n]->exons.Count()==1 && pred[n]->cov>=maxcov) || (pred[n]->exons.Count()>1 && pred[n]->cov/maxcov>=isofrac))	{ // print this transcript
		  if(genes[pred[n]->geneno]==-1) genes[pred[n]->geneno]=++geneno;
		  transcripts[pred[n]->geneno]++;
		  if(pred[n]->flag) {
			  uint t_id=0;
			  if (pred[n]->t_eq && pred[n]->t_eq->uptr) {
				  t_id = ((RC_TData*)pred[n]->t_eq->uptr)->t_id;
			  }
			  //fprintf(f_out,"%d %d %d %.6f %.6f\n",pred[n]->exons.Count()+1,pred[n]->tlen, t_id, pred[n]->frag,pred[n]->cov);
			  fprintf(f_out,"%d %d %d %.6f\n",pred[n]->exons.Count()+1,pred[n]->tlen, t_id,pred[n]->cov);
			  fprintf(f_out,"%s\tStringTie\ttranscript\t%d\t%d\t1000\t%c\t.\tgene_id \"%s.%d\"; transcript_id \"%s.%d.%d\"; ",
				  refname.chars(),pred[n]->start,pred[n]->end,pred[n]->strand,label.chars(),genes[pred[n]->geneno],
				  label.chars(),genes[pred[n]->geneno],transcripts[pred[n]->geneno]);
			  if(pred[n]->t_eq) {
				  fprintf(f_out,"reference_id \"%s\"; ",pred[n]->t_eq->getID());
				  if (pred[n]->t_eq->getGeneID())
				    	 fprintf(f_out,"ref_gene_id \"%s\"; ",pred[n]->t_eq->getGeneID());
				  if (pred[n]->t_eq->getGeneName())
				    	 fprintf(f_out,"ref_gene_name \"%s\"; ",pred[n]->t_eq->getGeneName());
			  }
			  fprintf(f_out,"cov \"%.6f\";\n",pred[n]->cov);
			  for(int j=0;j<pred[n]->exons.Count();j++) {
				  fprintf(f_out,"%s\tStringTie\texon\t%d\t%d\t1000\t%c\t.\tgene_id \"%s.%d\"; transcript_id \"%s.%d.%d\"; exon_number \"%d\"; ",
			  		 refname.chars(),pred[n]->exons[j].start,pred[n]->exons[j].end,pred[n]->strand,label.chars(),genes[pred[n]->geneno],
			  		 label.chars(),genes[pred[n]->geneno],transcripts[pred[n]->geneno],j+1); // maybe add exon coverage here
				  if(pred[n]->t_eq) {
					  fprintf(f_out,"reference_id \"%s\"; ",pred[n]->t_eq->getID());
					  if (pred[n]->t_eq->getGeneID())
					    	 fprintf(f_out,"ref_gene_id \"%s\"; ",pred[n]->t_eq->getGeneID());
					  if (pred[n]->t_eq->getGeneName())
					    	 fprintf(f_out,"ref_gene_name \"%s\"; ",pred[n]->t_eq->getGeneName());
				  }
				  fprintf(f_out,"cov \"%.6f\";\n",pred[n]->exoncov[j]);
			  }
		  }
		  else pred[n]->flag=true;
	  }
	  else pred[n]->flag=false;
  }

  return(geneno);
}


void add_pred(GList<CPrediction>& pred,int x,int y, float cov) { // add single exon prediction y to prediction x

	if(pred[y]->start<pred[x]->exons[0].end) { // add y to first exon in x
		int addlen=0;
		if(pred[y]->end<pred[x]->start)  // predictions do not overlap
			addlen=pred[x]->start-pred[y]->end-1;
		pred[x]->cov=(pred[x]->cov*pred[x]->tlen+cov*pred[y]->tlen)/(pred[x]->tlen+addlen+pred[y]->tlen);
		pred[x]->exoncov[0]= (pred[x]->exoncov[0]*(pred[x]->exons[0].end-pred[x]->exons[0].start+1)+
						cov*pred[y]->tlen)/(pred[x]->exons[0].end-pred[x]->exons[0].start+1+addlen+pred[y]->tlen);
		if(pred[y]->start<pred[x]->start) {
			pred[x]->start=pred[y]->start;
			pred[x]->exons[0].start=pred[y]->start;
		}
		pred[x]->tlen+=addlen;
	}
	else { // add y to last exon in x
		int addlen=0;
		if(pred[x]->end<pred[y]->start) // predictions do not overlap
			addlen=pred[y]->start-pred[x]->end-1;
		pred[x]->cov=(pred[x]->cov*pred[x]->tlen+cov*pred[y]->tlen)/(pred[x]->tlen+addlen+pred[y]->tlen);
		pred[x]->exoncov.Last()= (pred[x]->exoncov.Last()*(pred[x]->exons.Last().end-pred[x]->exons.Last().start+1)+
						cov*pred[y]->tlen)/(pred[x]->exons.Last().end-pred[x]->exons.Last().start+1+addlen+pred[y]->tlen);
		if(pred[x]->end<pred[y]->end) {
			pred[x]->end=pred[y]->end;
			pred[x]->exons.Last().end=pred[y]->end;
		}
		pred[x]->tlen+=addlen;
	}

	//pred[x]->frag+=pred[y]->frag;

}

int printResults(BundleData* bundleData, int ngenes, int geneno, GStr& refname) {


	// print transcripts including the necessary isoform fraction cleanings
	GList<CPrediction>& pred = bundleData->pred;
	int npred=pred.Count();
	pred.setSorted(predCmp);

	// this are needed for gene abundance estimations
	GVec<CGene> predgene;
	GVec<CGene> refgene;
	GHash<int> hashgene;
	int startgno=geneno+1;

	// process predictions that equal the same single exon guide and stich them together
	GPVec<GffObj>& guides = bundleData->keepguides;
	if(guides.Count()) {

		// first create reference genes
		int gno=0;
		for(int i=0;i<guides.Count();i++) {
			GStr gid(guides[i]->getGeneID());
			const int *n=hashgene[gid.chars()];
			if(n) { // I've seen the gene before
				if(guides[i]->start<refgene[*n].start) refgene[*n].start=guides[i]->start;
				if(guides[i]->end>refgene[*n].end) refgene[*n].end=guides[i]->start;
				merge_exons(refgene[*n],guides[i]->exons); // to write this
			}
			else { // create gene and hash
				hashgene.Add(gid.chars(),new int(gno));
				CGene g(guides[i]->start,guides[i]->end,guides[i]->strand,guides[i]->getGeneID(),guides[i]->getGeneName());
				// now add the exons
				for(int j=0;j<guides[i]->exons.Count();j++) {
					GSeg ex(guides[i]->exons[j]->start,guides[i]->exons[j]->end);
					g.exons.Add(ex);
				}
				refgene.Add(g);
				gno++;
			}
		}


		/*
		// this version only stiches together single exon predictions that overlap a reference guide
		GHash<int> seenguide;
		for(int n=0;n<pred.Count();n++) {
			if(pred[n]->t_eq && pred[n]->exons.Count()==1) { // this can only happen for single exon genes
				if(!seenguide[pred[n]->t_eq->getID()]) seenguide.Add(pred[n]->t_eq->getID(),new int(n+1));
				else {
					const int* p=seenguide[pred[n]->t_eq->getID()];
					*p--;
					uint start= pred[*p]->start < pred[n]->start ? pred[*p]->start : pred[n]->start;
					uint end= pred[*p]->end > pred[n]->end ? pred[*p]->end : pred[n]->end;
					pred[*p]->cov=(pred[*p]->cov*pred[*p]->tlen+pred[n]->cov*pred[n]->tlen)/(end-start+1);
					pred[*p]->start=start;
					pred[*p]->end=end;
					pred[*p]->exons[0].start=start;
					pred[*p]->exons[0].end=end;
					pred[*p]->exoncov[0]=pred[*p]->cov;
					pred[*p]->tlen=end-start+1;
					pred[*p]->frag+=pred[n]->frag;
					pred.Exchange(n,pred.Count()-1);
					pred.Delete(pred.Count()-1);
				}
			}
		}
		seenguide.Clear();
		*/

		// this version is more inclusive by stiching together single exons to reference guides that overlap them, but doesn't print them -> this is still done later on
		GVec< GVec<int> > reflink(npred); reflink.Resize(npred);
		for(int n=0;n<npred;n++) {

			if(pred[n] && pred[n]->t_eq) {

			  //fprintf(stderr,"pred[%d]: start=%d end=%d ID=%s strand=%c refstrand=%c refstart=%d\n",n,pred[n]->start,pred[n]->end,pred[n]->t_eq->getID(),pred[n]->strand,pred[n]->t_eq->strand,pred[n]->t_eq->start);


				// check if there are single exon on the left of the predicted transcript
				int i=n-1;
				while(i>=0 && pred[i]->start>=pred[n]->t_eq->start) {
				  //fprintf(stderr,"%d: excnot=%d end=%d start=%d str=%c\n",i,pred[i]->exons.Count(), pred[i]->end,pred[n]->start, pred[i]->strand);
					if(pred[i]->exons.Count()==1 && pred[i]->end<pred[n]->start && (pred[i]->strand==pred[n]->strand || pred[i]->strand=='.')) {
						if(!pred[i]->t_eq || !strcmp(pred[i]->t_eq->getID(),pred[n]->t_eq->getID())) reflink[i].Add(n);
					}
					i--;
				}
				// check if there are single exon on the right of the predicted transcript
				i=n+1;
				while(i<npred && pred[i]->start<pred[n]->t_eq->end) {
					if(pred[i]->exons.Count()==1 && pred[i]->start>pred[n]->end && (pred[i]->strand==pred[n]->strand || pred[i]->strand=='.')) {
						if(pred[i]->t_eq && !strcmp(pred[i]->t_eq->getID(),pred[n]->t_eq->getID())) { // same ID is pred[n]
							reflink[i].Add(n);
						}
						else if(pred[i]->end<=pred[n]->t_eq->end && !pred[i]->t_eq) reflink[i].Add(n);
					}
					i++;
				}
			}
		}

		for(int n=0;n<npred;n++) {

		  //fprintf(stderr,"pred[%d] reflink_cnt=%d\n",n,reflink[n].Count());
			if(pred[n] && reflink[n].Count()) {
				int mindist=100000;
				float sumcov=0;
				uint start=0;
				for(int i=0;i<reflink[n].Count();i++) {
					int ri=reflink[n][i];
					//fprintf(stderr,"pred[%d] is linked to pred[%d]\n",n,ri);
					if(pred[ri]) {
						int d;
						if(pred[ri]->start<pred[n]->start) d= pred[n]->start-pred[ri]->end;
						else d=pred[ri]->start-pred[n]->end;
						if(d<mindist) {
							mindist=d;
							sumcov=pred[ri]->cov;
							start=pred[ri]->start;
						}
						else if(d==mindist) sumcov+=pred[ri]->cov;
					}
				}
				if(sumcov) {
					for(int i=0;i<reflink[n].Count();i++) {
						int ri=reflink[n][i];
						if(pred[ri] && pred[ri]->start==start) add_pred(pred,ri,n,pred[n]->cov*pred[ri]->cov/sumcov);
					}
					// now replace pred[y] with null
					CPrediction *p=pred[n];
					pred.Forget(n);
					delete p;
				}
			}
		}
		pred.Pack();
		pred.setSorted(predCmp);
		npred=pred.Count();

	}

	//pred.setSorted(predCmp);
	//pred.setSorted(true);

	int currentstartpos=-1;
	uint currentendpos=0;
	int nstartpos=0;
	int nendpos=0;
	int currentstartneg=-1;
	uint currentendneg=0;
	int nstartneg=0;
	int nendneg=0;
	GVec<int> genes(true); // for each gene remembers it's geneno
	genes.Resize(ngenes,-1);
	GVec<int> transcripts(true); // for each gene remembers how many transcripts were printed
	transcripts.Resize(ngenes,0);

	GPVec<CPrediction> pospred(false);
	GPVec<CPrediction> negpred(false);

	for(int n=0;n<npred;n++) {

		if(pred[n]->strand=='.') pred[n]->flag=false; // only let it print if it passes threshold for printing

		if(pred[n]->strand=='+' || pred[n]->strand=='.') {
			if(pred[n]->start > currentendpos) { // begin new cluster
				// first print predictions I've seen so far
				if(currentstartpos>-1) { // I've seen a cluster before
					switch (sensitivitylevel) {
					case 0: geneno=print_transcript_signcluster('+',pred,genes,transcripts,nstartpos,nendpos,geneno,refname);break;
					case 1: geneno=print_cluster(pospred,genes,transcripts,geneno,refname,refgene,hashgene,predgene,startgno);break;
					case 2: geneno=print_cluster_inclusion(pospred,genes,transcripts,geneno,refname);break;
					case 3: geneno=print_signcluster('+',pred,genes,transcripts,nstartpos,nendpos,geneno,refname);break;
					}
					pospred.Clear();
				}

				pospred.Add(pred[n]);

				currentstartpos=pred[n]->start;
				currentendpos=pred[n]->end;
				nstartpos=n;
				nendpos=n;
			}
			else {
				if(pred[n]->end > currentendpos) currentendpos=pred[n]->end;
				nendpos=n;
				pospred.Add(pred[n]);

			}
		}

		if(pred[n]->strand=='-' || pred[n]->strand=='.') {
			if(pred[n]->start > currentendneg) { // begin new cluster

				// first print predictions I've seen so far
				if(currentstartneg>-1) { // I've seen a cluster before

					switch (sensitivitylevel) {
					case 0: geneno=print_transcript_signcluster('-',pred,genes,transcripts,nstartneg,nendneg,geneno,refname);break;
					case 1: geneno=print_cluster(negpred,genes,transcripts,geneno,refname,refgene,hashgene,predgene,startgno);break;
					case 2: geneno=print_cluster_inclusion(negpred,genes,transcripts,geneno,refname);break;
					case 3: geneno=print_signcluster('-',pred,genes,transcripts,nstartneg,nendneg,geneno,refname);break;
					}
					negpred.Clear();

				}

				negpred.Add(pred[n]);

				currentstartneg=pred[n]->start;
				currentendneg=pred[n]->end;
				nstartneg=n;
				nendneg=n;
			}
			else {
				negpred.Add(pred[n]);
				if(pred[n]->end > currentendneg) currentendneg=pred[n]->end;
				nendneg=n;
			}
		}
	}

	if(currentstartpos>-1) { // I've seen a cluster before

		switch (sensitivitylevel) {
		case 0: geneno=print_transcript_signcluster('+',pred,genes,transcripts,nstartpos,nendpos,geneno,refname);break;
		case 1: geneno=print_cluster(pospred,genes,transcripts,geneno,refname,refgene,hashgene,predgene,startgno);break;
		case 2: geneno=print_cluster_inclusion(pospred,genes,transcripts,geneno,refname);break;
		case 3: geneno=print_signcluster('+',pred,genes,transcripts,nstartpos,nendpos,geneno,refname);break;
		}
		pospred.Clear();

	}

	if(currentstartneg>-1) { // I've seen a cluster before

		switch (sensitivitylevel) {
		case 0: geneno=print_transcript_signcluster('-',pred,genes,transcripts,nstartneg,nendneg,geneno,refname);break;
		case 1: geneno=print_cluster(negpred,genes,transcripts,geneno,refname,refgene,hashgene,predgene,startgno);break;
		case 2: geneno=print_cluster_inclusion(negpred,genes,transcripts,geneno,refname);break;
		case 3: geneno=print_signcluster('-',pred,genes,transcripts,nstartneg,nendneg,geneno,refname);break;
		}
		negpred.Clear();

	}

	hashgene.Clear();
	// I am done printing all transcripts, now evaluate/print the gene abundances
	GVec<float>* bpcov = bundleData->bpcov;
	int refstart=bundleData->start;
	if(eonly  || geneabundance) { // I only need to evaluate the refgene coverages if I geneabundance is required, or these are the only gene coverages
		for(int i=0;i<refgene.Count();i++) {
			float cov=0;
			int s=1; // strand of gene
			if(refgene[i].strand=='+') s=2;
			else if(refgene[i].strand=='-') s=0;
			int glen=0;
			for(int j=0;j<refgene[i].exons.Count();j++) { // evaluate unused coverage
				glen+=refgene[i].exons[j].len();
				int start=(int)refgene[i].exons[j].start-refstart;
				int end=(int)refgene[i].exons[j].end-refstart+1;
				if(start<0) start=0;
				if(end>=bpcov[1].Count()) end=bpcov[1].Count()-1;
				//fprintf(stderr,"start=%d\nend=%d\ngene_start=%d gene_end=%d\nguides.count=%d refstart=%d",refgene[i].exons[j].start,refgene[i].exons[j].end,refgene[i].start,refgene[i].end,guides.Count(),refstart);
				for(int k=start;k<end;k++) {
					switch(s) {
					case 0: cov+=bpcov[1][k]-bpcov[2][k];break;
					case 1: cov+=bpcov[1][k]-bpcov[2][k]-bpcov[0][k];break;
					case 2: cov+=bpcov[1][k]-bpcov[0][k];break;
					}
				}
			}
			refgene[i].cov=cov-refgene[i].cov;
			if(refgene[i].cov>epsilon) refgene[i].covsum+=refgene[i].cov/glen;
			if(eonly) bundleData->sum_cov+=refgene[i].covsum;
			if(geneabundance) {
				refgene[i].cov=cov/glen; // only if I want to store the real gene coverage
				fprintf(f_out,"0 1 %d 0 %.6f\n",glen, refgene[i].covsum);
				fprintf(f_out,"%s\t",refgene[i].geneID);
				if(refgene[i].geneName) fprintf(f_out,"%s\t",refgene[i].geneName);
				else fprintf(f_out,"-\t");
				fprintf(f_out,"%c\t%d\t%d\t%d\t%.6f\n",refgene[i].strand,refgene[i].start,refgene[i].end,glen,refgene[i].cov);
			}
		}
	}
	if(!eonly) {
		for(int i=0;i<predgene.Count();i++) {
			float cov=0;
			int s=1; // strand of gene
			if(predgene[i].strand=='+') s=2;
			else if(predgene[i].strand=='-') s=0;
			int glen=0;
			for(int j=0;j<predgene[i].exons.Count();j++) { // evaluate unused coverage
				glen+=predgene[i].exons[j].len();
				int start=(int)predgene[i].exons[j].start-refstart;
				int end=(int)predgene[i].exons[j].end-refstart+1;
				// predgene start and end might have been also adjusted to reflect the annotation
				if(start<0) start=0;
				if(end>=bpcov[1].Count()) end=bpcov[1].Count()-1;
				for(int k=start;k<end;k++) {
					switch(s) {
					case 0: cov+=bpcov[1][k]-bpcov[2][k];break;
					case 1: cov+=bpcov[1][k]-bpcov[2][k]-bpcov[0][k];break;
					case 2: cov+=bpcov[1][k]-bpcov[0][k];break;
					}
				}
			}
			predgene[i].cov=cov-predgene[i].cov;
			if(predgene[i].cov>epsilon) predgene[i].covsum+=predgene[i].cov/glen;
			bundleData->sum_cov+=predgene[i].covsum;
			if(geneabundance) {
				predgene[i].cov=cov/glen; // only if I want to store the real gene coverage
				fprintf(f_out,"0 1 %d 0 %.6f\n",glen, predgene[i].covsum);
				fprintf(f_out,"%s.%d\t",label.chars(),startgno+i);
				fprintf(f_out,"-\t");
				fprintf(f_out,"%c\t%d\t%d\t%d\t%.6f\n",predgene[i].strand,predgene[i].start,predgene[i].end,glen,predgene[i].cov);
			}
		}
	}


	if (c_out) {
		for (int i=0;i<bundleData->covguides.Count();i++)
			bundleData->covguides[i]->print(c_out);
	}
	//rc_write_counts(refname.chars(), *bundleData);
	return(geneno);
}
