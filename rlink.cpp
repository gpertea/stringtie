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

//extern bool specific;
extern bool trim;
//extern bool partialcov;
//extern bool complete; // this is false if I use guides that are fragments, but I want them to be true
extern bool eonly;
extern bool nomulti;

extern bool fast;
//extern int maxReadCov;

extern float isofrac;
extern float mcov;
extern int mintranscriptlen; // minimum number for a transcript to be printed
//extern int sensitivitylevel;
extern uint junctionsupport; // anchor length for junction to be considered well supported <- consider shorter??
extern int junctionthr; // number of reads needed to support a particular junction
extern float readthr;     // read coverage per bundle bp to accept it; otherwise considered noise
extern uint bundledist;  // reads at what distance should be considered part of separate bundles
extern bool includesource;
//extern bool EM;
//extern bool weight;
extern bool geneabundance; // need to compute the gene abundance

extern float fpkm_thr;
extern float tpm_thr;
extern bool enableNames;
extern bool includecov;
extern bool retained_intron;

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

/*
void cov_add(GVec<float>& bpcov, int i, int j, float v) {
	if (j>=bpcov.Count()) bpcov.Resize(j+1, 0);
	for (int k=i;k<j;k++)
		bpcov[k]+=v;
}
*/

void cov_add(GVec<float>* bpcov, int sno, int i, int j, float v) {
	bool neutral=false;
	if(sno!=1) neutral=true; // why is neutral true here: because if the sno is -/+ than I want to add their counts to bpcov[1] too
	if (j>=bpcov[sno].Count())
		for(int s=0;s<3;s++) bpcov[s].Resize(j+1, 0);
	for (int k=i;k<=j;k++) {
		bpcov[sno][k]+=v;
		if(neutral) bpcov[1][k]+=v; // neutral (stranded) gets added twice here
	}
}

/*
float getBCov(GVec<float>& bpcov, int p) {
	//if (p<0) GMessage("Error: invalid bpcov index (%d)!\n", p);
	if (p>=bpcov.Count()) return 0;
	else return bpcov[p];
}


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

	bool match=false;  // true if current read matches a previous read
	int n=readlist.Count()-1;
	if(!mergeMode) while(n>-1 && readlist[n]->start==brec.start) {
		if(strand==readlist[n]->strand) match=exonmatch(readlist[n]->segs,brec.exons);
		if(match) break; // this way I make sure that I keep the n of the matching readlist
		n--;
	}
	else if((alndata.tinfo->cov>=0 && alndata.tinfo->cov<readthr) || (alndata.tinfo->fpkm>=0 && alndata.tinfo->fpkm<fpkm_thr) ||
			(alndata.tinfo->tpm>=0 && alndata.tinfo->tpm < tpm_thr)) return; // do not store 'read' if it doesn't meet minimum criteria
	else { //mergeMode but the read is above thresholds
		if(alndata.tinfo->cov>=0) bdata.covflags |= IS_COV_FLAG;
		if(alndata.tinfo->fpkm>=0) bdata.covflags |= IS_FPKM_FLAG;
		if(alndata.tinfo->tpm>=0) bdata.covflags |= IS_TPM_FLAG;
	}

	if (bdata.end<currentend) {// I am not sure why this is done here?
		bdata.start=currentstart;
		bdata.end=currentend;
	}
	bdata.numreads++;                // number of reads (alignments) actually considered
	//bdata.wnumreads+=float(1)/nh;

	if (!match) { // if this is a new read I am seeing I need to set it up
		if(mergeMode && mintranscriptlen) {
			int len=0;
			for (int i=0;i<brec.exons.Count();i++) len+=brec.exons[i].len();
			if(len<mintranscriptlen) return;
		}
		readaln=new CReadAln(strand, nh, brec.start, brec.end, alndata.tinfo);
		alndata.tinfo=NULL; //alndata.tinfo was passed to CReadAln
		for (int i=0;i<brec.exons.Count();i++) {
			readaln->len+=brec.exons[i].len();
			if(i) {
				CJunction* nj=add_junction(brec.exons[i-1].end, brec.exons[i].start, junction, strand);
				if (alndata.juncs.Count())
					nj->guide_match=alndata.juncs[i-1]->guide_match;
				if (nj) {
					readaln->juncs.Add(nj);
				}
			}
			readaln->segs.Add(brec.exons[i]);
		}
		n=readlist.Add(readaln);  // reset n for the case there is no match
	}
	else { //redundant read alignment matching a previous alignment
		// keep shortest nh so that I can see for each particular read the multi-hit proportion
		if(nh<readlist[n]->nh) readlist[n]->nh=nh;
		/*
		//for mergeMode, we have to free the transcript info:
		if (alndata.tinfo!=NULL) {
			 delete alndata.tinfo;
			 alndata.tinfo=NULL;
		}
		*/
	}

	if((int)brec.end>currentend) {
			currentend=brec.end;
	  	bdata.end=currentend;
	}

	float rdcount=1;
	float unitig_cov=brec.tag_float("YK");
	if(unitig_cov) rdcount=unitig_cov;

	if(!nomulti) rdcount/=nh;
	readlist[n]->read_count+=rdcount; // increase single count just in case I don't see the pair

	// store the mismatch count per junction so that I can eliminate it later
	int nm=brec.tag_int("NM"); // read mismatch
	if(!nm) {
		nm=brec.tag_int("nM"); // paired mismatch : big problem with STAR alignments
		if(brec.isPaired()) nm/=2;
	}
	// this counts reads over allowed fraction of mismatches and multi-mapped reads
	if(nh>1 || nm/readlist[n]->len>mismatchfrac) for(int i=0;i<readlist[n]->juncs.Count();i++) {
		readlist[n]->juncs[i]->nm+=rdcount;
	}

	// now set up the pairing
	if (brec.refId()==brec.mate_refId()) {  //only consider mate pairing data if mates are on the same chromosome/contig and are properly paired
	//if (brec.refId()==brec.mate_refId() && brec.isProperlyPaired()) {  //only consider mate pairing data if mates are on the same chromosome/contig and are properly paired
	//if (brec.isProperlyPaired()) {  //only consider mate pairing data if mates  are properly paired
		int pairstart=brec.mate_start();
		if (currentstart<=pairstart) { // if pairstart is in a previous bundle I don't care about it

			//GStr readname(brec.name());
			GStr id(brec.name(), 16); // init id with readname + append buffer
			if(pairstart<=readstart) { // if I've seen the pair already <- I might not have seen it yet because the pair starts at the same place
				id+='-';id+=pairstart;
				id+=".=";id+=hi; // (!) this useless suffix actually speeds up the hash by improving distribution!
				const int* np=hashread[id.chars()];
				if(np) { // the pair was stored --> why wouldn't it be? : only in the case that the pair starts at the same position
					if(readlist[*np]->nh>nh && !nomulti) rdcount=float(1)/readlist[*np]->nh;
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
				id+='-';id+=readstart; // this is the correct way
				id+=".=";id+=hi;
				hashread.fAdd(id.chars(), new int(n));
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

void add_group_to_bundle(CGroup *group, CBundle *bundle, GPVec<CBundlenode>& bnode, uint bundledist){

	CBundlenode *currlastnode=bnode[bundle->lastnodeid];
	int bid=bnode.Count();

	if(group->start > currlastnode->end + bundledist) { // group after last bnode
		CBundlenode *currbnode=new CBundlenode(group->start,group->end,group->cov_sum,bid);
		currlastnode->nextnode=currbnode;
		bnode.Add(currbnode);
		bundle->lastnodeid=bid;
		bundle->len+=group->end-group->start+1;
		bundle->cov+=group->cov_sum;
		bundle->nread+=group->nread;
		bundle->multi+=group->multi;
	}
	else { // group overlaps bnode within bundledist
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
	if(a->nodes.Last()-a->nodes[0]<b->nodes.Last()-b->nodes[0]) return 1; // transfrag with further reach comes first
	if(a->nodes.Last()-a->nodes[0]>b->nodes.Last()-b->nodes[0]) return -1;
	if(a->nodes.Count()<b->nodes.Count()) return 1;  // transfrag with more nodes comes first
	if(a->nodes.Count()>b->nodes.Count()) return -1;
	if(a->pattern.count()<b->pattern.count()) return 1;
	if(a->pattern.count()>b->pattern.count()) return -1;
	if(a->abundance<b->abundance) return 1;
	if(a->abundance>b->abundance) return -1;
	return 0;
}


int mgtrnodeCmp(const pointer p1, const pointer p2) {
	CMTransfrag *a=(CMTransfrag*)p1;
	CMTransfrag *b=(CMTransfrag*)p2;
	if(a->transfrag->pattern.count()<b->transfrag->pattern.count()) return 1;
	if(a->transfrag->pattern.count()>b->transfrag->pattern.count()) return -1;
	if(!a->transfrag->real && b->transfrag->real) return -1;
	if(a->transfrag->real && !b->transfrag->real) return 1;
	/*if(a->transfrag->abundance<b->transfrag->abundance) return 1;
	if(a->transfrag->abundance>b->transfrag->abundance) return -1;*/
	return 0;
}

int mgtrabundCmp(const pointer p1, const pointer p2) {
	CMTransfrag *a=(CMTransfrag*)p1;
	CMTransfrag *b=(CMTransfrag*)p2;
	if(!a->transfrag->real && b->transfrag->real) return 1;    // guides come first
	if(a->transfrag->real && !b->transfrag->real) return -1;
	if(a->transfrag->abundance<b->transfrag->abundance) return 1; // most abundant transcript comes first
	if(a->transfrag->abundance>b->transfrag->abundance) return -1;
	if(a->transfrag->pattern.count()<b->transfrag->pattern.count()) return 1;  // the one with largest pattern comes first
	if(a->transfrag->pattern.count()>b->transfrag->pattern.count()) return -1;

	// something to decide equal cases
	if(a->transfrag->nodes.Count()<b->transfrag->nodes.Count()) return 1;  // the one with most nodes comes first
	if(a->transfrag->nodes.Count()>b->transfrag->nodes.Count()) return -1;

	int i=0;
	while(i<a->transfrag->nodes.Count()) {
		if(a->transfrag->nodes[i]<b->transfrag->nodes[i]) return 1; // the transcript starting more to the left comes first
		if(a->transfrag->nodes[i]>b->transfrag->nodes[i]) return -1;
		i++;
	}

	return 0;
}

int partguideCmp(const pointer p1, const pointer p2) { // sorting partial guides in order of preference
	CPartGuide *a=(CPartGuide*)p1;
	CPartGuide *b=(CPartGuide*)p2;
	if(a->cov<b->cov) return 1;
	if(a->cov>b->cov) return -1;
	if(a->gcount<b->gcount) return 1;
	if(a->gcount>b->gcount) return -1;
	if(a->olen<b->olen) return 1;
	if(a->olen>b->olen) return -1;
	if(a->allolen<b->allolen) return 1;
	if(a->allolen>b->allolen) return -1;
	if(a->glen<b->glen) return 1;
	if(a->glen>b->glen) return -1;
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
	if(a->trf->real && !b->trf->real) return 1; // this ensures included guides are treated last
	if(!a->trf->real && b->trf->real) return -1;
	if(a->trf->abundance<b->trf->abundance) return 1; // most abundant guide takes precedence
	if(a->trf->abundance>b->trf->abundance) return -1;
	if(a->trf->pattern.count()<b->trf->pattern.count()) return 1; // longest guide takes precedence
	if(a->trf->pattern.count()>b->trf->pattern.count()) return -1;
	return 0;
}

int guideCmp(const pointer p1, const pointer p2) {
	CGuide *a=(CGuide*)p1;
	CGuide *b=(CGuide*)p2;
	if(a->trf->pattern.count()<b->trf->pattern.count()) return 1; // longest guide takes precedence
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

		    		//fprintf(stderr,"2 Split pairs: %d-%d and %d-%d on strand %d\n",readlist[np]->start,readlist[np]->end,readlist[n]->start,readlist[n]->end,readlist[n]->strand);
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

			//fprintf(stderr,"3 Split pairs: %d-%d and %d-%d on strand %d\n",readlist[np]->start,readlist[np]->end,readlist[n]->start,readlist[n]->end,readlist[n]->strand);
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

	if(!mergeMode) for(int p=0;p<readlist[n]->pair_idx.Count();p++) {

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
	if(mergeMode || single_count>epsilon) { // my way of controlling for rounding errors
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
		// I wasn't squaring before but it doesn't make sense not to because then I would add a bunch of negatives
		if(mul) chi+= (winleft[j]-mul)*(winleft[j]-mul)/mul;
		if(mur) chi+=(winright[j]-mur)*(winright[j]-mur)/mur;
	}
	return(chi);
}


float compute_chi2(GArray<float>& winleft, GArray<float>& winright, float sumleft, float sumright) {

	float chi=0;
	float basecost=0;
	float leftcost=0;
	float rightcost=0;
	float mubase=(sumleft+sumright)/(2*CHI_WIN);
	float muleft=sumleft/CHI_WIN;
	float muright=sumright/CHI_WIN;
	for(int j=0;j<CHI_WIN;j++) {
		basecost+=(winleft[j]-mubase)*(winleft[j]-mubase)+(winright[j]-mubase)*(winright[j]-mubase);
		leftcost+=(winleft[j]-muleft)*(winleft[j]-muleft);
		rightcost+=(winright[j]-muright)*(winright[j]-muright);
	}
	chi=basecost/(2*CHI_WIN)-(leftcost+rightcost)/CHI_WIN;
	return(chi);
}

float compute_cost(GVec<float>& wincov,float sumleft,float sumright,int leftstart,int rightstart,int rightend) {
	float cost=0;
	float basecost=0;
	float leftcost=0;
	float rightcost=0;
	int len=rightend-leftstart+1;
	int leftlen=rightstart-leftstart;
	int rightlen=rightend-rightstart+1;
	float mubase=(sumleft+sumright)/len;
	float muleft=sumleft/leftlen;
	float muright=sumright/rightlen;
	for(int j=leftstart;j<rightstart;j++) {
		//basecost+=(wincov[j]-mubase)*(wincov[j]-mubase);
		basecost+=fabs(wincov[j]-mubase);
		//leftcost+=(wincov[j]-muleft)*(wincov[j]-muleft);
		leftcost+=fabs(wincov[j]-muleft);
	}
	for(int j=rightstart;j<=rightend;j++) {
		//basecost+=(wincov[j]-mubase)*(wincov[j]-mubase);
		basecost+=fabs(wincov[j]-mubase);
		//rightcost+=(wincov[j]-muright)*(wincov[j]-muright);
		rightcost+=fabs(wincov[j]-muright);
	}

	cost=basecost/len-leftcost/leftlen-rightcost/rightlen;

	/*
	fprintf(stderr,"leftlen=%d leftsum=%f leftcost=%f muleft=%f rightlen=%d righttsum=%f rightcost=%f muright=%f baselen=%d basesum=%f basecost=%f mubase=%f\n",
			leftlen,sumleft,leftcost/leftlen,muleft,
			rightlen,sumright,rightcost/rightlen,muright,
			len,sumleft+sumright,basecost/len,mubase);
	*/
	return(cost);
}

float compute_cost(GArray<float>& wincov,int i,float sumleft,float sumright,int len) {
	float cost=0;
	float basecost=0;
	float leftcost=0;
	float rightcost=0;
	float mubase=(sumleft+sumright)/len;
	float muleft=sumleft/i;
	float muright=sumright/(len-i);
	for(int j=0;j<len;j++) {
		basecost+=(wincov[j]-mubase)*(wincov[j]-mubase);
		if(j<i) leftcost+=(wincov[j]-muleft)*(wincov[j]-muleft);
		else rightcost+=(wincov[j]-muright)*(wincov[j]-muright);
	}
	cost=basecost/len-leftcost/i-rightcost/(len-i);
	return(cost);
}


void find_trims(int refstart,int sno,uint start,uint end,GVec<float>* bpcov,uint& sourcestart,float& maxsourceabundance,uint& sinkend,
		float& maxsinkabundance, float& sourcecovleft,float& sourcecovright,float& sinkcovleft, float& sinkcovright){


	int len=end-start+1;

	//fprintf(stderr,"find_trimms for start=%d end=%d len=%d sno=%d\n",start,end,len,sno);

	if(len<2*CHI_WIN) return;

	int bw=3; // bandwidth to consider for changes in coverage; should always be less than chi_win

	GVec<float> wincov;
	GVec<float> upleft;
	GVec<float> downleft;
	GVec<int> upleftstart;
	GVec<int> downleftstart;
	GVec<float> currcovleft;

	float cov=0;
	int j=0;

	float prevcov=0;
	float chileft=0;
	float chiright=0;

	for(int i=0;i<bw;i++) { // store prevcov
		upleft.Add(prevcov);
		downleft.Add(prevcov);
		j=0;
		upleftstart.Add(j);
		downleftstart.Add(j);

		j=i+start-refstart;
		cov=bpcov[1][j];
		if(sno!=1) cov-=bpcov[2-sno][j];

		wincov.Add(cov);
		prevcov+=cov;
		currcovleft.Add(prevcov);

		if(i>len-1-CHI_WIN) chiright+=cov;
		else if(i>len-1-2*CHI_WIN) chileft+=cov;

	}

	float nextcov=0;

	for(int i=bw; i<2*bw; i++) { // compute nextcov
		j=i+start-refstart;
		cov=bpcov[1][j];
		if(sno!=1) cov-=bpcov[2-sno][j];

		nextcov+=cov;
		wincov.Add(cov);
		cov+=currcovleft.Last();
		currcovleft.Add(cov);

		if(i>len-1-CHI_WIN) chiright+=cov;
		else if(i>len-1-2*CHI_WIN) chileft+=cov;

	}

	for(int i=2*bw; i<len; i++) { // compute left trends and store wincov

		// check trends
		float add=0;
		if(prevcov<nextcov) { // downtrend
			upleft.Add(add);
			j=i-bw;
			upleftstart.Add(j);
			add=wincov[i-bw-1]+downleft.Last();
			downleft.Add(add);
			j=downleftstart.Last();
			downleftstart.Add(j);
		}
		else if(prevcov>nextcov) { // uptrend
			downleft.Add(add);
			j=i-bw;
			downleftstart.Add(j);
			add=wincov[i-bw-1]+upleft.Last();
			upleft.Add(add);
			j=upleftstart.Last();
			upleftstart.Add(j);
		}
		else { // equal
			add=wincov[i-bw-1]+downleft.Last();
			downleft.Add(add);
			j=downleftstart.Last();
			downleftstart.Add(j);
			add=wincov[i-bw-1]+upleft.Last();
			upleft.Add(add);
			j=upleftstart.Last();
			upleftstart.Add(j);
		}

		j=i+start-refstart;
		cov=bpcov[1][j];
		if(sno!=1) cov-=bpcov[2-sno][j];

		wincov.Add(cov);
		prevcov-=wincov[i-2*bw];
		prevcov+=wincov[i-bw];
		nextcov-=wincov[i-bw];
		nextcov+=cov;

		if(i>len-1-CHI_WIN) chiright+=cov;
		else if(i>len-1-2*CHI_WIN) chileft+=cov;

		cov+=currcovleft.Last();
		currcovleft.Add(cov);

	}

	float upright=0;
	float downright=0;
	int uprightend=len-1;
	int downrightend=len-1;

	nextcov=0;
	prevcov=0;
	for(int i=len-1;i>=len-bw;i--) {
		nextcov+=wincov[i];
	}
	for(int i=len-bw-1;i>len-2*bw-1;i--) {
		prevcov+=wincov[i];
	}

	float opt_source=0;
	float opt_sink=0;

	float currcovright=chiright;

	for(int i=len-2*CHI_WIN;i>=0;i--) {


		if(25<(chiright-chileft)/CHI_WIN) { // possible cut to source
			float sumleft=chileft+downleft[i];
			float sumright=chiright+upright;
			float muleft=sumleft/(i+CHI_WIN-downleftstart[i]);
			float muright=sumright/(uprightend-i-CHI_WIN+1);

			//fprintf(stderr,"Pos=%d chileft=%f chiright=%f muleft=%f muright=%f\n",i+CHI_WIN+start,chileft,chiright,muleft,muright);

			if(muright*isofrac>muleft) {
				//float opt_cost=muright-muleft; // b30
				float opt_cost=(muright+1)/(muleft+1); //b31: so that I don't get into under 1 roundings
				if(opt_cost>opt_source) {
					//fprintf(stderr,"sourcestart=%d sumleft=%f sumright=%f lenleft=%d lenright=%d muleft=%f muright=%f opt_cost=%f\n",
					//i+CHI_WIN+start,sumleft,sumright,CHI_WIN+i-downleftstart[i],uprightend-i-CHI_WIN+1,muleft,muright,opt_cost);
					//maxsourceabundance=opt_cost;// b21
				  //maxsourceabundance=2*opt_cost;// b30
					maxsourceabundance=2*(muright-muleft);// b31

				  //maxsourceabundance=(chiright-chileft)/CHI_WIN; // b22
					opt_source=opt_cost;
					sourcestart=i+CHI_WIN+start;
					sourcecovleft=currcovleft[i+CHI_WIN-1];
					sourcecovright=currcovright;
				}
			}
		}
		else if(25<(chileft-chiright)/CHI_WIN) { // possible cut to sink
			float sumleft=chileft+upleft[i];
			float sumright=chiright+downright;

			float muleft=sumleft/(i+CHI_WIN-upleftstart[i]);
			float muright=sumright/(downrightend-i-CHI_WIN+1);
			//fprintf(stderr,"sinkend=%d sumleft=%f sumright=%f lenleft=%d lenright=%d muleft=%f muright=%f upleftstart=%d\n",
			//		i+CHI_WIN-1+start,sumleft,sumright,i+CHI_WIN-upleftstart[i],downrightend-i-CHI_WIN+1,muleft,muright,upleftstart[i]);
			if(muleft*isofrac>muright) {
				// float opt_cost=muleft-muright; // b30
				float opt_cost=(muleft+1)/(muright+1); // b31: so that I don't get into under 1 roundings

				//fprintf(stderr,"sinkend=%d sumleft=%f sumright=%f lenleft=%d lenright=%d muleft=%f muright=%f opt_cost=%f\n",
				//		i+CHI_WIN-1+start,sumleft,sumright,i+CHI_WIN-upleftstart[i],downrightend-i-CHI_WIN+1,muleft,muright,opt_cost);
				if(opt_cost>opt_sink) {
					//fprintf(stderr,"optimal\n");
				  //maxsinkabundance=opt_cost;// b21
				  //maxsinkabundance=2*opt_cost;// b30
					maxsinkabundance=2*(muleft-muright);// b31

				  //maxsinkabundance=(chileft-chiright)/CHI_WIN;//b22
					opt_sink=opt_cost;
					sinkend=i+CHI_WIN-1+start;
					sinkcovleft=currcovleft[i+CHI_WIN-1];
					sinkcovright=currcovright;
				}
			}
		}

		j=i+2*CHI_WIN-1;

		//fprintf(stderr,"j=%d cov=%f prevcov=%f nextcov=%f downright=%f downrightend=%d\n",j+start,wincov[j],prevcov,nextcov,downright,downrightend+start);

		if(j>len-bw) { // no need to update prevcov and next cov
			upright+=wincov[j];
			downright+=wincov[j];
		}
		else { // I have all info necessary for next time
			if(prevcov<nextcov) { // uptrend
				downright=0;
				downrightend=j-1;
				upright+=wincov[j];
			}
			else if(prevcov>nextcov) { // downtrend
				upright=0;
				uprightend=j-1;
				downright+=wincov[j];
			}
			else { // equal
				upright+=wincov[j];
				downright+=wincov[j];
			}

			nextcov-=wincov[j+bw-1];
			nextcov+=wincov[j-1];
			prevcov-=wincov[j-1];
			prevcov+=wincov[j-bw-1];

			if(prevcov<0) prevcov=0;
			if(nextcov<0) nextcov=0;

		}

		if(i) {
			chiright-=wincov[j];
			chiright+=wincov[i+CHI_WIN-1];
			chileft-=wincov[i+CHI_WIN-1];
			chileft+=wincov[i-1];

			currcovright+=wincov[i+CHI_WIN-1];

			if(chileft<0) chileft=0;
			if(chiright<0) chiright=0;
		}

	}

}



/*
// seems reasonable
void find_trims(int refstart,int sno,uint start,uint end,GVec<float>* bpcov,uint& sourcestart,float& maxsourceabundance,uint& sinkend,
		float& maxsinkabundance){


	int len=end-start+1;

	fprintf(stderr,"find_trimms for start=%d end=%d len=%d\n",start,end,len);

	if(len<2*CHI_WIN) return;

	GVec<float> wincov;
	GVec<float> upleft;
	GVec<float> downleft;
	GVec<int> upleftstart;
	GVec<int> downleftstart;

	float cov=0;
	upleft.Add(cov);
	downleft.Add(cov);

	int j=0;
	upleftstart.Add(j);
	downleftstart.Add(j);

	float chileft=0;
	float chiright=0;

	for(uint i=start;i<=end;i++) {
		cov=bpcov[sno][i-refstart];
		if(bpcov[1][i-refstart]>cov) { // it means bpcov[1]>0 since bpcov can not be negative and sno!=1 because it can't be neutral group
			if(bpcov[2-sno][i-refstart]) // I have stranded coverage on different strands
				cov+=(bpcov[1][i-refstart]-bpcov[0][i-refstart]-bpcov[2][i-refstart])*cov/(cov+bpcov[2-sno][i-refstart]);
			else cov=bpcov[1][i-refstart];
		}

		if(i>start) {
			float prevcov=wincov.Last();
			float add=0;
			if(prevcov<cov) { // downtrend
				upleft.Add(add);
				j=i-start;
				upleftstart.Add(j);
				add=prevcov+downleft.Last();
				downleft.Add(add);
				j=downleftstart.Last();
				downleftstart.Add(j);
			}
			else if(prevcov>cov) { // uptrend
				downleft.Add(add);
				j=i-start;
				downleftstart.Add(j);
				add=prevcov+upleft.Last();
				upleft.Add(add);
				j=upleftstart.Last();
				upleftstart.Add(j);
			}
			else { // equal
				add=prevcov+downleft.Last();
				downleft.Add(add);
				j=downleftstart.Last();
				downleftstart.Add(j);
				add=prevcov+upleft.Last();
				upleft.Add(add);
				j=upleftstart.Last();
				upleftstart.Add(j);
			}
		}
		wincov.Add(cov);

		if(i>end-CHI_WIN) chiright+=cov;
		else if(i>end-2*CHI_WIN) chileft+=cov;

		//if(i<start+CHI_WIN) chileft+=cov;
		//else if(i<start+2*CHI_WIN) chiright+=cov;
	}

	float upright=0;
	float downright=0;
	int uprightend=len-1;
	int downrightend=len-1;

	//float opt_source=0;
	//float opt_sink=0;

	float nextcov=wincov[len-1];

	for(int i=len-2*CHI_WIN;i>=0;i--) {

		// *
		if(chileft*isofrac>chiright || chiright*isofrac>chileft) { // possible drop (sink cut) or increase (source cut)
			if(chileft<chiright) { // possible cut to source
				float sumleft=chileft+downleft[i];
				float sumright=chiright+upright;
				float opt_cost=compute_cost(wincov,sumleft,sumright,downleftstart[i],i+CHI_WIN,uprightend);
				fprintf(stderr,"source optcost=%f at position=%d\n",opt_cost,i+CHI_WIN+start);
				if(opt_cost>opt_source) {
					maxsourceabundance=sumright/(uprightend-i-CHI_WIN+1)-sumleft/(i-CHI_WIN-downleftstart[i]);
					opt_source=opt_cost;
					sourcestart=i+CHI_WIN+start;
				}
			}
			else if(chiright<chileft) { // possible cut to sink
				float sumleft=chileft+upleft[i];
				float sumright=chiright+downright;
				float opt_cost=compute_cost(wincov,sumleft,sumright,upleftstart[i],i+CHI_WIN,downrightend);
				fprintf(stderr,"sink optcost=%f at position=%d\n",opt_cost,i+CHI_WIN-1+start);
				if(opt_cost>opt_sink) {
					maxsinkabundance=sumleft/(i+CHI_WIN-upleftstart[i])-sumright/(downrightend-i-CHI_WIN+1);
					opt_sink=opt_cost;
					sinkend=i+CHI_WIN-1+start;
				}
			}
		}
		* //


		if(chileft<chiright) { // possible cut to source
			float sumleft=chileft+downleft[i];
			float sumright=chiright+upright;
			float muleft=sumleft/(i+CHI_WIN-downleftstart[i]);
			float muright=sumright/(uprightend-i-CHI_WIN+1);
			if(muright*isofrac>muleft) {
				float opt_cost=muright-muleft;
				fprintf(stderr,"sourcestart=%d sumleft=%f sumright=%f lenleft=%d lenright=%d muleft=%f muright=%f opt_cost=%f\n",
						i+CHI_WIN+start,sumleft,sumright,i+CHI_WIN-downleftstart[i],uprightend-i-CHI_WIN+1,muleft,muright,opt_cost);
				if(opt_cost>maxsourceabundance) {
					maxsourceabundance=opt_cost;
					sourcestart=i+CHI_WIN+start;
				}
			}
		}
		else if(chiright<chileft) { // possible cut to sink
			float sumleft=chileft+upleft[i];
			float sumright=chiright+downright;
			float muleft=sumleft/(i+CHI_WIN-upleftstart[i]);
			float muright=sumright/(downrightend-i-CHI_WIN+1);
			fprintf(stderr,"sinkend=%d sumleft=%f sumright=%f lenleft=%d lenright=%d muleft=%f muright=%f\n",
					i+CHI_WIN-1+start,sumleft,sumright,i+CHI_WIN-upleftstart[i],downrightend-i-CHI_WIN+1,muleft,muright);
			if(muleft*isofrac>muright) {
				float opt_cost=muleft-muright;
				fprintf(stderr,"sinkend=%d sumleft=%f sumright=%f lenleft=%d lenright=%d muleft=%f muright=%f opt_cost=%f\n",
						i+CHI_WIN-1+start,sumleft,sumright,i+CHI_WIN-upleftstart[i],downrightend-i-CHI_WIN+1,muleft,muright,opt_cost);
				if(opt_cost>maxsinkabundance) {
					maxsinkabundance=opt_cost;
					sinkend=i+CHI_WIN-1+start;
				}
			}
		}

		cov=wincov[i+2*CHI_WIN-2]; // this is the next ending coverage
		if(cov<nextcov) { // uptrend
			downright=0;
			downrightend=i+2*CHI_WIN-2;
			upright+=nextcov;
		}
		else if(cov>nextcov) { // downtrend
			upright=0;
			uprightend=i+2*CHI_WIN-2;
			downright+=nextcov;
		}
		else { // equal
			upright+=nextcov;
			downright+=nextcov;
		}

		if(i) {
			chiright-=nextcov;
			chiright+=wincov[i+CHI_WIN-1];
			chileft-=wincov[i+CHI_WIN-1];
			chileft+=wincov[i-1];
		}

		nextcov=cov;

	}

}
*/


/*
void find_trims(int refstart,int sno,uint start,uint end,GVec<float>* bpcov,uint& sourcestart,float& maxsourceabundance,uint& sinkend,
		float& maxsinkabundance){

	if(end-start<2*CHI_WIN-1) return;

	float sumleft=0;
	float sumright=0;

	float maxsinkchi=0;
	float maxsourcechi=0;

	//float sinkabundance=0;
	//float sourceabundance=0;

	GArray<float> winleft(CHI_WIN,false); // not auto-sort
	GArray<float> winright(CHI_WIN,false); // not auto-sort

	float cov;

	for(uint i=start;i<=end;i++) {

		cov=bpcov[sno][i-refstart];
		if(bpcov[1][i-refstart]>cov) { // it means bpcov[1]>0 since bpcov can not be negative and sno!=1 because it can't be neutral group
			//cov+=(bpcov[1][i-refstart]-bpcov[0][i-refstart]-bpcov[2][i-refstart])/bpcov[1][i-refstart];
			if(bpcov[2-sno][i-refstart]) // I have stranded coverage on different strand
				cov+=(bpcov[1][i-refstart]-bpcov[0][i-refstart]-bpcov[2][i-refstart])*cov/(cov+bpcov[2-sno][i-refstart]);
			else cov=bpcov[1][i-refstart];
		}

		if(i-start<2*CHI_WIN-1)  { // I have to compute the sumleft and sumright first
			if(i-start<CHI_WIN) {
				sumleft+=cov;
				winleft.Add(cov);
			}
			else {
				sumright+=cov;
				winright.Add(cov);
				if(i-start==2*CHI_WIN-2) {
					winleft.setSorted(true);
					winright.setSorted(true);
				}
			}
	    }
	    else { // I can do the actual sumleft, sumright comparision
	    	sumright+=cov;
	    	winright.Add(cov);

		if(fabs(sumleft-sumright)/CHI_WIN>50) { // what is the threshold at which the difference becomes significant?
			float chi=0;
			if(sumleft!=sumright) chi=compute_chi(winleft,winright,sumleft,sumright);

			if(chi>epsilon) { // there is a significant difference

				//fprintf(stderr,"i=%d chi=%f sumleft=%f sumright=%f\n",i,chi,sumleft,sumright);

				if(sumleft*isofrac>sumright) { // possible drop (sink cut)
					if(chi>maxsinkchi) {
					  maxsinkabundance=(sumleft-sumright)/CHI_WIN;
					  // * how much bigger
					  // if(sumright) maxsinkabundance=sumleft/sumright;
					  // else maxsinkabundance=100;
					  // *
					  maxsinkchi=chi;
					  sinkend=i-CHI_WIN;

					}
				}
				else if(sumright*isofrac>sumleft) { // increase (source cut)

					//fprintf(stderr,"i=%d chi=%f sumleft=%f sumright=%f sno=%d\n",i,chi,sumleft,sumright,sno);

					if(chi>maxsourcechi) {
					  maxsourceabundance=(sumright-sumleft)/CHI_WIN;
					  // * how much bigger
					  // if(sumleft) maxsourceabundance=sumright/sumleft;
					  // else maxsourceabundance=100;
					  // *
					  sourcestart=i-CHI_WIN+1;
					  maxsourcechi=chi;
					}
					}
				}
			}

			cov=bpcov[sno][i-refstart-2*CHI_WIN+1];
			if(bpcov[1][i-refstart-2*CHI_WIN+1]>cov) {
				//cov+=(bpcov[1][i-refstart-2*CHI_WIN+1]-bpcov[0][i-refstart-2*CHI_WIN+1]-bpcov[2][i-refstart-2*CHI_WIN+1])/bpcov[1][i-refstart-2*CHI_WIN+1];
				if(bpcov[2-sno][i-refstart-2*CHI_WIN+1]) // I have stranded coverage on different strand
					cov+=(bpcov[1][i-refstart-2*CHI_WIN+1]-bpcov[0][i-refstart-2*CHI_WIN+1]-bpcov[2][i-refstart-2*CHI_WIN+1])*cov/(cov+bpcov[2-sno][i-refstart-2*CHI_WIN+1]);
				else cov=bpcov[1][i-refstart-2*CHI_WIN+1];
			}

	    	sumleft-=cov;
	    	int idx=winleft.IndexOf(cov);
	    	if(idx>=0) winleft.Delete(idx);

			cov=bpcov[sno][i-refstart-CHI_WIN+1];
			if(bpcov[1][i-refstart-CHI_WIN+1]>cov) {
				//cov+=(bpcov[1][i-refstart-CHI_WIN+1]-bpcov[0][i-refstart-CHI_WIN+1]-bpcov[2][i-refstart-CHI_WIN+1])/bpcov[1][i-refstart-CHI_WIN+1];
				if(bpcov[2-sno][i-refstart-CHI_WIN+1]) // I have stranded coverage on different strand
					cov+=(bpcov[1][i-refstart-CHI_WIN+1]-bpcov[0][i-refstart-CHI_WIN+1]-bpcov[2][i-refstart-CHI_WIN+1])*cov/(cov+bpcov[2-sno][i-refstart-CHI_WIN+1]);
				else cov=bpcov[1][i-refstart-CHI_WIN+1];
			}

	    	sumleft+=cov;
	    	winleft.Add(cov);
	    	sumright-=cov;
	    	idx=winright.IndexOf(cov);
	    	if(idx>=0) winright.Delete(idx);
	    }
	}

}
*/

CGraphnode *add_trim_to_graph(int s, int g,uint lastpos,CTrimPoint& mytrim,CGraphnode *graphnode,CGraphnode *source,CGraphnode *sink,GVec<float>& futuretr,
		int& graphno,CBundlenode *bundlenode,GVec<CGraphinfo> **bundle2graph,GPVec<CGraphnode> **no2gnode) {

	if(mytrim.start) { // this is a source link
		float tmp=graphno-1;
		CGraphnode *prevnode=NULL;
		if(mytrim.pos>graphnode->start) { // there is place for another node
			uint prevend=graphnode->end;
			graphnode->end=mytrim.pos-1;
			prevnode=graphnode;
			//fprintf(stderr,"add trim 1\n");
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
			//fprintf(stderr,"add trim 2\n");
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

CGraphnode *source2guide(int s, int g, int refstart,uint newstart,uint newend, CGraphnode *graphnode,CGraphnode *source,
		GVec<float>* bpcov,GVec<float>& futuretr, int& graphno,CBundlenode *bundlenode,GVec<CGraphinfo> **bundle2graph,
		GPVec<CGraphnode> **no2gnode, int &edgeno) {

	//fprintf(stderr,"source2guide for newstart=%d and newend=%d\n",newstart,newend);

	// compute maxabund
	float leftcov=0;
	float rightcov=0;
	for(uint p=graphnode->start;p<newstart;p++) {
		int j=p-refstart;
		leftcov+=bpcov[1][j]-bpcov[2-2*s][j];
	}
	for(uint p=newstart;p<newend;p++) {
		int j=p-refstart;
		rightcov+=bpcov[1][j]-bpcov[2-2*s][j];
	}
	float maxabund=rightcov-leftcov;
	if(maxabund<trthr) maxabund=trthr;

	uint prevend=graphnode->end;
	graphnode->end=newstart-1;
	CGraphnode *prevnode=graphnode;
	graphnode=create_graphnode(s,g,newstart,prevend,graphno,bundlenode,bundle2graph,no2gnode);
	graphno++;
	source->child.Add(graphnode->nodeid);  // this node is the child of source
	graphnode->parent.Add(source->nodeid); // this node has source as parent
	prevnode->child.Add(graphnode->nodeid); // this node is the child of previous node
	graphnode->parent.Add(prevnode->nodeid); // this node has as parent the previous node
	float tmp=graphno-1;
	futuretr.cAdd(0.0);
	futuretr.Add(tmp);
	futuretr.Add(maxabund);
	tmp=prevnode->nodeid;futuretr.Add(tmp);
	tmp=graphnode->nodeid;futuretr.Add(tmp);
	tmp=trthr;futuretr.Add(tmp);
	// COUNT 1 EDGE HERE because the source to guide edge was already included in our count
	edgeno++;

	return(graphnode);

}

CGraphnode *guide2sink(int s, int g, int refstart,uint newstart,uint newend, CGraphnode *graphnode,CGraphnode *sink,
		GVec<float>* bpcov,GVec<float>& futuretr, int& graphno,CBundlenode *bundlenode,GVec<CGraphinfo> **bundle2graph,
		GPVec<CGraphnode> **no2gnode, int &edgeno) {

	//fprintf(stderr,"guide2sink for newstart=%d and newend=%d\n",newstart,newend);

	// compute maxabund
	float leftcov=0;
	float rightcov=0;
	for(uint p=graphnode->start;p<=newstart;p++) {
		int j=p-refstart;
		leftcov+=bpcov[1][j]-bpcov[2-2*s][j];
	}
	for(uint p=newstart+1;p<newend;p++) {
		int j=p-refstart;
		rightcov+=bpcov[1][j]-bpcov[2-2*s][j];
	}
	float maxabund=leftcov-rightcov;
	if(maxabund<trthr) maxabund=trthr;

	float tmp=graphno-1;
	uint prevend=graphnode->end;
	graphnode->end=newstart;
	CGraphnode *prevnode=graphnode;
	graphnode=create_graphnode(s,g,newstart+1,prevend,graphno,bundlenode,bundle2graph,no2gnode);
	graphno++;
	prevnode->child.Add(graphnode->nodeid); // this node is the child of previous node
	graphnode->parent.Add(prevnode->nodeid); // this node has as parent the previous node
	sink->parent.Add(prevnode->nodeid); // prevnode is the parent of sink
	futuretr.Add(tmp);
	futuretr.cAdd(-1.0);
	futuretr.Add(maxabund);
	tmp=prevnode->nodeid;futuretr.Add(tmp);
	tmp=graphnode->nodeid;futuretr.Add(tmp);
	tmp=trthr;futuretr.Add(tmp);
	// COUNT 1 EDGE HERE because the source to guide edge was already included in our count
	edgeno++;

	return(graphnode);

}

CGraphnode *trimnode(int s, int g, int refstart,uint newend, CGraphnode *graphnode,CGraphnode *source, CGraphnode *sink, GVec<float>* bpcov,
		GVec<float>& futuretr, int& graphno,CBundlenode *bundlenode,GVec<CGraphinfo> **bundle2graph,GPVec<CGraphnode> **no2gnode, int &edgeno) {

	uint sourcestart=0;
	uint sinkend=0;
	float sinkabundance=0;
	float sourceabundance=0;
	float sourcecovleft=0;
	float sourcecovright=0;
	float sinkcovleft=0;
	float sinkcovright=0;
	find_trims(refstart,2*s,graphnode->start,newend,bpcov,sourcestart,sourceabundance,sinkend,sinkabundance,sourcecovleft,sourcecovright,sinkcovleft,sinkcovright);

	if(sourcestart < sinkend) { // source trimming comes first

		if(sourcestart) { // there is evidence of graphnode trimming from source
			graphnode->end=sourcestart-1;
			CGraphnode *prevnode=graphnode;
			//fprintf(stderr,"trimnode 1\n");
			graphnode=create_graphnode(s,g,sourcestart,newend,graphno,bundlenode,bundle2graph,no2gnode);
			graphno++;
			source->child.Add(graphnode->nodeid);  // this node is the child of source
			graphnode->parent.Add(source->nodeid); // this node has source as parent
			//fprintf(stderr,"trim edge 0-%d with sourceabund=%f at position=%d\n",graphnode->nodeid,sourceabundance,graphnode->start);
			prevnode->child.Add(graphnode->nodeid); // this node is the child of previous node
			graphnode->parent.Add(prevnode->nodeid); // this node has as parent the previous node
			//fprintf(stderr,"trim edge %d-%d\n",prevnode->nodeid,graphnode->nodeid);
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
		//fprintf(stderr,"trimnode 2\n");
		graphnode=create_graphnode(s,g,sinkend+1,newend,graphno,bundlenode,bundle2graph,no2gnode);
		graphno++;
		prevnode->child.Add(graphnode->nodeid); // this node is the child of previous node
		graphnode->parent.Add(prevnode->nodeid); // this node has as parent the previous node
		//fprintf(stderr,"trim edge %d-%d with sinkabundance=%f at position=%d\n",prevnode->nodeid,graphnode->nodeid,sinkabundance,prevnode->end);
		sink->parent.Add(prevnode->nodeid); // prevnode is the parent of sink
		// remember to create transfrag as well -> I don't know the gno yet, so I can not create it here
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
			//fprintf(stderr,"trimnode 3\n");
			graphnode=create_graphnode(s,g,sinkend+1,newend,graphno,bundlenode,bundle2graph,no2gnode);
			graphno++;
			prevnode->child.Add(graphnode->nodeid); // this node is the child of previous node
			graphnode->parent.Add(prevnode->nodeid); // this node has as parent the previous node
			//fprintf(stderr,"trim edge %d-%d\n",prevnode->nodeid,graphnode->nodeid);
			sink->parent.Add(prevnode->nodeid); // prevnode is the parent of sink
			//fprintf(stderr,"trim edge %d-sink with sinkabundance=%f at position=%d\n",prevnode->nodeid,sinkabundance,prevnode->end);
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
		//fprintf(stderr,"trimnode 3\n");
		graphnode=create_graphnode(s,g,sourcestart,newend,graphno,bundlenode,bundle2graph,no2gnode);
		graphno++;
		source->child.Add(graphnode->nodeid);  // this node is the child of source
		graphnode->parent.Add(source->nodeid); // this node has source as parent
		//fprintf(stderr,"trim edge 0-%d with sourceabundance=%f at position=%d\n",graphnode->nodeid,sourceabundance,graphnode->start);
		prevnode->child.Add(graphnode->nodeid); // this node is the child of previous node
		graphnode->parent.Add(prevnode->nodeid); // this node has as parent the previous node
		//fprintf(stderr,"trim edge %d-%d\n",prevnode->nodeid,graphnode->nodeid);
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
		//fprintf(stderr,"Add %d children of node %d (%d-%d): ",n,node->nodeid,node->start,node->end);
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

			//fprintf(stderr,"min=%d max=%d key=%d lastgpos=%d\n",min,max,key,lastgpos);

			if(pos!=NULL) {
				childparents[*pos]=1; // add edge from node to child to the set of parents from child
				node->childpat[*pos]=1; // add edge from node to child to the set of node children
			}
			else {
				gpos[s][g].Add(key,lastgpos);
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

int create_graph(int refstart,int s,int g,CBundle *bundle,GPVec<CBundlenode>& bnode,
		GList<CJunction>& junction,GList<CJunction>& ejunction,GVec<CGraphinfo> **bundle2graph,
		GPVec<CGraphnode> **no2gnode,GPVec<CTransfrag> **transfrag,GIntHash<int> **gpos,GVec<float>* bpcov,int &edgeno,
		int &lastgpos,GArray<GEdge>& guideedge, int refend=0){

	CGraphnode* source=new CGraphnode(0,0,0);
	no2gnode[s][g].Add(source);
	CGraphnode* sink=new CGraphnode();

	int njunctions=junction.Count();

	//fprintf(stderr,"Start with %d edgeno and lastgpos=%d\n",edgeno,lastgpos);

	/*
	{ // DEBUG ONLY
		fprintf(stderr,"Junctions[%d][%d]: ",s,g);
		for(int i=0;i<njunctions;i++) fprintf(stderr," %d-%d:%d",junction[i]->start,junction[i]->end,junction[i]->strand);
		fprintf(stderr,"\n");
		fprintf(stderr,"eJunctions[%d][%d]: ",s,g);
		for(int i=0;i<njunctions;i++) fprintf(stderr," %d-%d:%d",ejunction[i]->start,ejunction[i]->end,ejunction[i]->strand);
		fprintf(stderr,"\n");
	}
	*/

	int nge=0;
	bool processguide=false;
	if(nge<guideedge.Count()) processguide=true;

	int njs=0; // index of sorted junction starts
	int nje=0; // index of sorted junction ends

	int graphno=1; // number of nodes in graph
	GHash<GVec<int> > ends; // keeps ids of all nodes ending at a certain position; OR ALL NODES THAT ARE LINKED BY JUNCTIONS TO A CERTAIN POSITION

	CBundlenode *bundlenode=bnode[bundle->startnode];

	GVec<float> futuretr;

	if(mergeMode) { // I have a bunch of junctions at the start for which I need to create ends

		while(njs<njunctions && !junction[njs]->start ) { // remember ends here for source node
			if((junction[njs]->strand+1) == 2*s) {
				GStr je((int)junction[njs]->end);
				GVec<int> *e=ends[je.chars()];
				if(!e) {
					e = new GVec<int>();
					ends.Add(je.chars(),e);
				}
				e->cAdd(0);
			}
			njs++;
		}

	}

	//int seenjunc=0;

	while(bundlenode!=NULL) {

		//fprintf(stderr,"process bundlenode %d-%d\n",bundlenode->start,bundlenode->end);

	    uint currentstart=bundlenode->start; // current start is bundlenode's start
	    uint endbundle=bundlenode->end; // initialize end with bundlenode's end for now

	    //fprintf(stderr,"create graph 1\n");
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
	    			//fprintf(stderr,"1 Edge %d-%d, edgeno=%d\n",node->nodeid,graphnode->nodeid,edgeno);
	    		}
	    	}
	    	else { // I haven't seen nodes before that finish here => link to source
		    	source->child.Add(graphnode->nodeid);  // this node is the child of source
		    	graphnode->parent.Add(source->nodeid); // this node has source as parent
		    	// COUNT EDGE HERE
		    	edgeno++;
		    	//fprintf(stderr,"2 Edge 0-%d, edgeno=%d\n",graphnode->nodeid,edgeno);
	    	}
	    }
	    else { // this node comes from source directly
	    	source->child.Add(graphnode->nodeid);  // this node is the child of source
	    	graphnode->parent.Add(source->nodeid); // this node has source as parent
	    	// COUNT EDGE HERE
			edgeno++;
			//fprintf(stderr,"3 Edge 0-%d, edgeno=%d\n",graphnode->nodeid,edgeno);
	    }


	    bool completed=false;

	    do {

	    	while(nje<njunctions && (((int)ejunction[nje]->strand+1) != 2*s)) nje++; // skip junctions that don't have the same strand
	    	while(njs<njunctions && ((((int)junction[njs]->strand+1)!= 2*s) || (junction[njs]->start<currentstart))) njs++; // junctions that start before the current graphnode and I haven't seen them before are part of a different bundle


	    	int minjunction = -1; // process next junction -> either a start or an ending whichever has the first position on the genome; if they have same position then process ending first
	    	if((nje<njunctions && (ejunction[nje]->end<=endbundle)) || (njs<njunctions && (junction[njs]->start<=endbundle))) {
	    		if(nje<njunctions) { // there are still junctions endings
	    			if(njs<njunctions) { // there are still junctions starting
	    				minjunction = junction[njs]->start >= ejunction[nje]->end ? 1 : 0; // one of them is clearly before the endbundle from the initial if
	    			}
	    			else minjunction = 1;
	    		}
	    		else minjunction = 0;
	    	}

	    	//fprintf(stderr,"minjunction=%d\n",minjunction);
	    	//if(nje<njunctions) fprintf(stderr,"Found junction:%d-%d(%d)\n",ejunction[nje]->start,ejunction[nje]->end,ejunction[nje]->strand);

	    	if(minjunction == 0 ) { // found a start junction here

	    		// add guide starts/ends first
	    		if(processguide) {
	    			while(nge<guideedge.Count() && guideedge[nge].val<=graphnode->start) nge++;
	    			if(nge<guideedge.Count()) {

	    				while(true) {

	    					while(nge<guideedge.Count() && guideedge[nge].strand!=s) nge++;

	    					if(nge>=guideedge.Count() || guideedge[nge].val>=junction[njs]->start) break;

	    					//fprintf(stderr,"guideedge val=%d endval=%d strand=%d currentstart=%d endbundle=%d\n",guideedge[nge].val,guideedge[nge].endval,guideedge[nge].strand,currentstart,endbundle);

	    					uint start=guideedge[nge].val;
	    					uint end=junction[njs]->start;
	    					bool sourceguide=false;
	    					if(guideedge[nge].val<guideedge[nge].endval) sourceguide=true;
	    					nge++;
	    					if(sourceguide) { if(guideedge[nge-1].endval>endbundle) continue;}
	    					else if(guideedge[nge-1].endval<currentstart) continue;

	    					while(nge<guideedge.Count() && guideedge[nge].strand!=s) nge++;
	    					if(nge<guideedge.Count() && guideedge[nge].val<junction[njs]->start) end=guideedge[nge].val;

	    					if(sourceguide)	graphnode=source2guide(s,g,refstart,start,end,graphnode,source,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);
	    					else graphnode=guide2sink(s,g,refstart,start,end,graphnode,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);

	    				}
	    			}
	    		}

	    		if(trim && !mergeMode) graphnode=trimnode(s,g,refstart,junction[njs]->start,graphnode,source,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);// do something to find intermediate nodes; alternatively, I could only do this for end nodes

	    		// if no trimming required just set the end of the node
	    		graphnode->end=junction[njs]->start; // set the end of current graphnode to here; introduce smaller nodes if trimming is activated
	    		uint pos=junction[njs]->start;
	    		while(njs<njunctions && junction[njs]->start==pos ) { // remember ends here
	    			if((junction[njs]->strand+1) == 2*s) {

	    				//seenjunc++;
	    				if(mergeMode && (int)junction[njs]->end==refend) { // this node goes straight to sink
	    					sink->parent.Add(graphnode->nodeid); // graphnode is the parent of sink: check to see if I have a conflict with this
	    					edgeno++; // count edge here
	    				}
	    				else {
	    					GStr je((int)junction[njs]->end);
	    					GVec<int> *e=ends[je.chars()];
	    					if(!e) {
	    						e = new GVec<int>();
	    						ends.Add(je.chars(),e);
	    					}
	    					e->Add(graphnode->nodeid);
	    				}
	    			}
	    			njs++;
	    		}

	    		if(pos<endbundle) { // there is still place for another node in this bundle (I might put a limit of length here for the graphnode -> because otherwise one can assume this is just a pre-mRNA fragment)
	    			//fprintf(stderr,"create graph 2\n");
	    			CGraphnode *nextnode = create_graphnode(s,g,pos+1,endbundle,graphno,bundlenode,bundle2graph,no2gnode);
	    			graphno++;
	    			graphnode->child.Add(nextnode->nodeid); // make nextnode a child of current graphnode
	    			nextnode->parent.Add(graphnode->nodeid);// make graphnode a parent of nextnode
	    			// COUNT EDGE HERE
	    			edgeno++;
	    			//fprintf(stderr,"4 Edge %d-%d, edgeno=%d nextnode: %u-%u pos=%d\n",graphnode->nodeid,nextnode->nodeid,edgeno,nextnode->start,nextnode->end,pos);
	    			graphnode=nextnode;
	    		}
	    		else completed=true;
	    	}
	    	else if(minjunction == 1) { // found a junction end here

	    		uint pos=ejunction[nje]->end;
	    		while(nje<njunctions && ejunction[nje]->end==pos) { // read all junction ends at the current start
	    			nje++;
	    		}

	    		if(graphnode->start<pos) { // last created node starts before the position of the new node I want to create

	    			// add guide starts/ends first
	    			if(processguide) {
	    				while(nge<guideedge.Count() && guideedge[nge].val<=graphnode->start) nge++;
	    				if(nge<guideedge.Count()) {

	    					while(true) {

	    						while(nge<guideedge.Count() && guideedge[nge].strand!=s) nge++;

	    						if(nge>=guideedge.Count() || guideedge[nge].val>=pos-1) break;

	    						//fprintf(stderr,"guideedge val=%d endval=%d strand=%d currentstart=%d endbundle=%d\n",guideedge[nge].val,guideedge[nge].endval,guideedge[nge].strand,currentstart,endbundle);

	    						uint start=guideedge[nge].val;
	    						uint end=pos-1;
	    						bool sourceguide=false;
	    						if(guideedge[nge].val<guideedge[nge].endval) sourceguide=true;
	    						nge++;
	    						if(sourceguide) { if(guideedge[nge-1].endval>endbundle) continue;}
	    						else if(guideedge[nge-1].endval<currentstart) continue;

	    						while(nge<guideedge.Count() && guideedge[nge].strand!=s) nge++;
	    						if(nge<guideedge.Count() && guideedge[nge].val<pos-1) end=guideedge[nge].val;

	    						if(sourceguide)	graphnode=source2guide(s,g,refstart,start,end,graphnode,source,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);
	    						else graphnode=guide2sink(s,g,refstart,start,end,graphnode,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);

	    					}
	    				}
	    			}

	    			if(trim && !mergeMode) graphnode=trimnode(s,g,refstart,pos-1,graphnode,source,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);// do something to find intermediate nodes; alternatively, I could only do this for end nodes

	    			graphnode->end=pos-1; // set end of current graphnode here
	    			//fprintf(stderr,"create graph 3\n");
	    			CGraphnode *nextnode = create_graphnode(s,g,pos,endbundle,graphno,bundlenode,bundle2graph,no2gnode);
	    			graphno++;
	    			graphnode->child.Add(nextnode->nodeid); // make nextnode a child of current graphnode
	    			nextnode->parent.Add(graphnode->nodeid);// make graphnode a parent of nextnode

	    			// COUNT EDGE HERE
	    			edgeno++;
	    			//fprintf(stderr,"5 Edge %d-%d, edgeno=%d\n",graphnode->nodeid,nextnode->nodeid,edgeno);

	    			graphnode=nextnode;
	    		}

	    		GStr spos((int)pos);
	    		GVec<int> *e=ends[spos.chars()]; // WHY DOESN'T THIS REPEAT THE SAME THING IN CASE THE START HASN'T BEEN ADJUSTED? because nje is bigger now than the ones that end at the currentstart
	    		if(e) for(int i=0;i<e->Count();i++) {
	    			CGraphnode *node=no2gnode[s][g][e->Get(i)];
	    			node->child.Add(graphnode->nodeid);  // this node is the child of previous node
	    			graphnode->parent.Add(node->nodeid); // this node has as parent the previous node
	    			// COUNT EDGE HERE
	    			edgeno++;
	    			//fprintf(stderr,"6 Edge %d-%d, edgeno=%d\n",node->nodeid,graphnode->nodeid,edgeno);
	    		}
	    	}

	    	//if(nje<njunctions) fprintf(stderr,"ejunc: %d-%d\n",ejunction[nje]->start,ejunction[nje]->end);
	    	//if(njs<njunctions) fprintf(stderr,"junc: %d-%d\n",junction[njs]->start,junction[njs]->end);

	    } while((nje<njunctions && (ejunction[nje]->end<=endbundle)) || (njs<njunctions && (junction[njs]->start<=endbundle)));


	    if(!completed) { // I did not finish node --> this will be an ending node

	    	// add guide starts/ends first
	    	if(processguide) {
	    		while(nge<guideedge.Count() && guideedge[nge].val<=graphnode->start) nge++;
	    		if(nge<guideedge.Count()) {

	    			while(true) {

	    				while(nge<guideedge.Count() && guideedge[nge].strand!=s) nge++;

	    				if(nge>=guideedge.Count() || guideedge[nge].val>=endbundle) break;

						//fprintf(stderr,"guideedge val=%d endval=%d strand=%d currentstart=%d endbundle=%d\n",guideedge[nge].val,guideedge[nge].endval,guideedge[nge].strand,currentstart,endbundle);

	    				uint start=guideedge[nge].val;
	    				uint end=endbundle;
	    				bool sourceguide=false;
	    				if(guideedge[nge].val<guideedge[nge].endval) sourceguide=true;
	    				nge++;
	    				if(sourceguide) { if(guideedge[nge-1].endval>endbundle) continue;}
	    				else if(guideedge[nge-1].endval<currentstart) continue;

	    				while(nge<guideedge.Count() && guideedge[nge].strand!=s) nge++;
	    				if(nge<guideedge.Count() && guideedge[nge].val<endbundle) end=guideedge[nge].val;

	    				if(sourceguide)	graphnode=source2guide(s,g,refstart,start,end,graphnode,source,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);
	    				else graphnode=guide2sink(s,g,refstart,start,end,graphnode,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);

	    			}
	    		}
	    	}

	    	if(trim && !mergeMode) // do something to find intermediate nodes; alternatively, I could only do this for end nodes
	    		graphnode=trimnode(s,g,refstart,endbundle,graphnode,source,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);

	    	graphnode->end=endbundle;
	    	// COUNT EDGE HERE (this is an edge to sink)
	    	edgeno++;
	    	//fprintf(stderr,"7 Edge to sink from %d, edgeno=%d\n",graphnode->nodeid,edgeno);
	    }

	    bundlenode=bundlenode->nextnode; // advance to next bundle
	} // end while(bundlenode!=NULL)

	sink->nodeid=graphno;
	no2gnode[s][g].Add(sink);
	graphno++;

	if(mergeMode) { // I might have a bunch of sink's parents that are not linked to sink
		for(int i=0;i<sink->parent.Count();i++) {
			CGraphnode *node=no2gnode[s][g][sink->parent[i]];
			node->child.Add(sink->nodeid);
		}
	}

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

	//fprintf(stderr,"traverse graph[%d][%d] now with %d nodes, %d edges and lastgpos=%d....\n",s,g,graphno,edgeno,lastgpos);//edgeno=0;
	traverse_dfs(s,g,source,sink,parents,graphno,visit,no2gnode,transfrag,edgeno,gpos,lastgpos);
	//fprintf(stderr,"done traversing with edgeno=%d lastgpos=%d\n",edgeno,lastgpos);

	// delete variables created here, like e.g. ends; do I need to delete the GVec<int> elements created too?
	ends.Clear();

	return(graphno);

}


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
	CTreePat *root=new CTreePat(0,gno-1); // if links from source to nodes are desired source==1 and all nodes are treated as +1

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

/* if I ever want to use the source
CTreePat *construct_treepat(int gno, GIntHash<int>& gpos,GPVec<CTransfrag>& transfrag) {

	int source=0; // if I ever need source set this to 1

	gno+=source;

	// create root CTreePat first
	CTreePat *root=new CTreePat(0,gno-1); // if links from source to nodes are desired source==1 and all nodes are treated as +1

	// now construct all child CTreePat's
	for(int t=0;t<transfrag.Count();t++)
		if(transfrag[t]->nodes[0]+source){ // don't include transfrags from source -> not needed
			CTreePat *tree=root;
			int m=0; // previous node in pattern that was set in pattern
			for(int n=1;n<gno;n++){
				if(transfrag[t]->pattern[n-source]) {
					CTreePat *child;
					if(m) { // there is a node m that was seen before
						int *pos=gpos[edge(m-source,n-source,gno-source)];
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
}*/

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

/* if I want to use the source
CTransfrag *findtrf_in_treepat(int gno,GIntHash<int>& gpos,GVec<int>& node,GBitVec& pattern,CTreePat *tr2no) { // doesn't work for patterns including source node

	int source=0; // if I ever need source, set this to 1

	gno+=source;

	if(!tr2no) return(NULL);

	CTreePat *tree=tr2no;
	for(int n=0;n<node.Count();n++) {
		if(n) { // not the first node in pattern
			int *pos=gpos[edge(node[n-1],node[n],gno-source)];
			if(pos && pattern[*pos]) // there is an edge between node[n-1] and node[n]
				tree=tree->nextpat[gno-1-node[n-1]-source+node[n]-node[n-1]-1];
			else tree=tree->nextpat[node[n]-node[n-1]-1];
		}
		else tree=tree->nextpat[node[n]+source-1];
		if(!tree) return(NULL);
	}

	return(tree->tr);
}*/


CTransfrag *update_abundance(int s,int g,int gno,GIntHash<int>&gpos,GBitVec& pattern,float abundance,GVec<int>& node,GPVec<CTransfrag> **transfrag,
		CTreePat ***tr2no){

	/*
	{ // DEBUG ONLY
		fprintf(stderr," (check nodes:");
		for(int i=0;i<node.Count();i++) fprintf(stderr," %d",node[i]);
		fprintf(stderr,")");
	}
	*/

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

	return(t);
}

void get_fragment_pattern(GList<CReadAln>& readlist,int n, int np,float readcov,GVec<int> *readgroup,GVec<int>& merge,
		GVec<int> *group2bundle,GVec<CGraphinfo> **bundle2graph,GVec<int> *graphno,GVec<int> *edgeno, GIntHash<int> **gpos,GPVec<CGraphnode> **no2gnode,
		GPVec<CTransfrag> **transfrag,CTreePat ***tr2no,GPVec<CGroup> &group) {

	//fprintf(stderr,"get fragment for read[%d]:%d-%d-%d with pair[%d]\n",n,readlist[n]->start,readlist[n]->end,int(readlist[n]->strand),np);

	GBitVec rpat[2];
	int rgno[2]={-1,-1};
	GVec<int> rnode[2];
	float rprop[2]={1,1};
	//bool goodfrag=false;

	// compute proportions of read associated to strands
	if(readlist[n]->nh && !readlist[n]->strand && np>-1 && readlist[np]->nh && !readlist[np]->strand) { // both reads are unstranded
		int gr1=readgroup[n][0]; // read is unstranded => it should belong to one group only
		while(merge[gr1]!=gr1) gr1=merge[gr1];
		int gr2=readgroup[np][0]; // read is unstranded => it should belong to one group only
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

void get_read_to_transfrag(GList<CReadAln>& readlist,int n,GVec<int> *readgroup,GVec<int>& merge,GVec<int> *group2bundle,
		GVec<CGraphinfo> **bundle2graph,GVec<int> *graphno,GVec<int> *edgeno,GIntHash<int> **gpos,GPVec<CGraphnode> **no2gnode,
		GPVec<CTransfrag> **transfrag,GPVec<CMTransfrag> **mgt,CTreePat ***tr2no,GPVec<CGroup> &group) {

	// compute proportions of read associated to strands
	float rprop[2]={1,1};
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


	GBitVec rpat[2];
	int rgno[2]={-1,-1};
	GVec<int> rnode[2];

	if(readlist[n]->nh) {
		get_read_pattern(readlist[n]->read_count,rpat[0],rpat[1],rgno,rprop,rnode,readlist,n,readgroup,merge,group2bundle,bundle2graph,graphno,edgeno,gpos,no2gnode);
	}

	for(int s=0;s<2;s++)
		if(rgno[s]>-1) { // read is valid (has pattern) on strand s
			// update transfrag
			CTransfrag *t=update_abundance(s,rgno[s],graphno[s][rgno[s]],gpos[s][rgno[s]],rpat[s],rprop[s]*readlist[n]->read_count,rnode[s],transfrag,tr2no);

			// if I want to inset source/sink maybe here would be the place

			// mark transfrag as guide
			if(readlist[n]->tinfo->g>-1) t->real=true;

			// update MTransfrag
			int i=0;
			while(i<mgt[s][rgno[s]].Count()) {
				if(mgt[s][rgno[s]][i]->transfrag==t) {
					for(int x=0;x<readlist[n]->pair_idx.Count();x++) {
						mgt[s][rgno[s]][i]->read.Add(readlist[n]->pair_idx[x]);
					}
					//mgt[s][rgno[s]][i]->read.Add(n);
					mgt[s][rgno[s]][i]->len=readlist[n]->len;
					break;
				}
				i++;
			}
			if(i==mgt[s][rgno[s]].Count()) { // MTransfrag not found
				CMTransfrag *mt=new CMTransfrag(t);
				for(int x=0;x<readlist[n]->pair_idx.Count();x++) {
					mt->read.Add(readlist[n]->pair_idx[x]);
				}
				//mt->read.Add(n);
				mt->len=readlist[n]->len;
				mgt[s][rgno[s]].Add(mt);
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
		for(int t=transfrag.Count()-1;t>=0;t--)
			if(transfrag[t]->abundance<threshold && transfrag[t]->nodes[0] && transfrag[t]->nodes.Last()<gno-1) { // need to delete transfrag that doesn't come from source
				settrf_in_treepat(NULL,gno,gpos,transfrag[t]->nodes,transfrag[t]->pattern,tr2no); // this should be eliminated if I want to store transcripts from 0 node
				transfrag.Exchange(t,transfrag.Count()-1);
				transfrag.Delete(transfrag.Count()-1);
			}
		threshold++;
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

// returns true if at list one node gets assigned the incomplete transfrag
bool binary_insert_trf_to_node(int t, GVec<int>& trf,int first,int last) {

	while(first<=last) {
		int mid=(first + last) / 2;  // compute mid point.
		if(t>trf[mid]) first=mid+1;
		else if(t<trf[mid]) last=mid-1;
		else return false; // transcript already inserted
	}

	int pos=first;

	trf.idxInsert(pos,t);
	return true;

}

bool assign_incomplete_trf_to_nodes(int t,int n1, int n2,GPVec<CGraphnode>& no2gnode){
	// a better version of this function could take into account also the fragments inner length

	bool assigned=false;

	for(int i=n1+1;i<n2;i++) {
		if(no2gnode[n1]->childpat[i] && no2gnode[i]->childpat[n2]) { // i is a child of n1, and n2 is a child of i
			assigned = binary_insert_trf_to_node(t,no2gnode[i]->trf,0,no2gnode[i]->trf.Count()-1) or assigned;
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

	return assigned;
}

inline int comptbl_pos(int t1,int t2,int n){
	return(t2+t1*(2*n-t1-1)/2);
}

// this function inspects if there are multiple ways between two nodes in the transfrag; finds the first node like this and distributes the abundance of the transfrag between the two nodes in proportion to the edges that leave the first node
bool trf_real(int t,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag, GIntHash<int> &gpos,int gno) {

	int i=0; // index of the current node (the one I am looking at) in transfrag
	int nt=transfrag[t]->nodes.Count()-1; // index of last node in transfrag
	int n=transfrag[t]->nodes[i]; // start at first node in transfrag
	float totalcov=0;
	while(i<nt) {
		//fprintf(stderr,"trf=%d n=%d nextnode=%d\n",t,n,transfrag[t]->nodes[i+1]);
		if(n+1==transfrag[t]->nodes[i+1]) { // if next node in transfrag is adjacent -> there is only one way to continue on transfrag path
			n=transfrag[t]->nodes[i+1];
			i++;
		}
		else { // there might be multiple ways to reach next node from node
			int *pos=gpos[edge(n,transfrag[t]->nodes[i+1],gno)];
			if(pos && transfrag[t]->pattern[*pos]) { // if there is an edge between node and nextnode in transfrag -> I don't have discontinuity there
				n=transfrag[t]->nodes[i+1];
				i++;
			}
			else { // found a potential discontinuity in transfrag pattern --> try to see if there are different branches going on from here
				CGraphnode *node=no2gnode[n];
				int nextnode=transfrag[t]->nodes[i+1];
				int ntrf=node->trf.Count(); // total number of transfrags in node
				int np=0;  // upper bound on number of paths I can find from node to nextnode
				totalcov=0; // total abundance of transfrags going from node to children of node

				//fprintf(stderr,"%d transcripts of node %d:",ntrf,n);
				//for(int f=0;f<ntrf;f++) fprintf(stderr," %d",node->trf[f]);
				//fprintf(stderr,"\n");

				// I check if there are different paths to reach nextnode from node
				for(int j=0;j<node->child.Count();j++) {
					if(node->child[j]>nextnode) break; // if child of node is after nextnode -> I don't care
					if(node->child[j]==nextnode || no2gnode[node->child[j]]->childpat[nextnode]) { // I found a node on the path to nextnode
						CPath p(n,node->child[j]);
						int *pos=gpos[edge(n,node->child[j],gno)];
						//fprintf(stderr,"Consider n=%d and child=%d with pos=%d and transcripts: ",n,node->child[j],*pos);
						if(pos) for(int f=0;f<ntrf;f++) {
							//fprintf(stderr," %d(%f)",node->trf[f],transfrag[node->trf[f]]->abundance);
							if(transfrag[node->trf[f]]->pattern[*pos]) { // this is a transfrag that goes through node and child; ATTENTION: it might not reach nextnode!! or be compatible with nextnode
								p.abundance+=transfrag[node->trf[f]]->abundance;
								totalcov+=transfrag[node->trf[f]]->abundance;
							}
						}
						//fprintf(stderr,"\n");
						if(p.abundance) {
							//fprintf(stderr,"Found abundance btw n=%d and child=%d\n",n,node->child[j]);
							transfrag[t]->path.Add(p);
							np++;
						}
					}
				}
				if(!totalcov || !np) // there is no way to continue
					return true;
				if(np==1) { // only one way to go forward
					n=transfrag[t]->path[0].contnode;
					if(n==nextnode) i++;
					transfrag[t]->path.Clear();
				}
				else break; // break at first discontinue point in transfrag -> might be many other ways to continue but this functions doesn't care
			}
		}
	}
	if(i==nt) return true;

	// this implementation only assumes one break in transfrag
	for(int j=0;j<transfrag[t]->path.Count();j++) transfrag[t]->path[j].abundance/=totalcov;

	return false;
}

void process_transfrags(int gno,int edgeno,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,CTreePat *tr2no,
		GBitVec& compatible,GIntHash<int> &gpos,GVec<CGuide>& guidetrf) {

	/*
	{ // DEBUG ONLY
		printTime(stderr);
		fprintf(stderr,"There are %d transfrags before clean up:\n",transfrag.Count());
		for(int i=0;i<transfrag.Count();i++) {
			fprintf(stderr,"transfrag[%d]:",i);
			for(int j=0;j<transfrag[i]->nodes.Count();j++) fprintf(stderr," %d",transfrag[i]->nodes[j]);
			fprintf(stderr,"\n");
		}
	}
	*/

	// eliminate transfrags below threshold (they represent noise) if they don't come from source
	eliminate_transfrags_under_thr(gno,gpos,transfrag,tr2no,trthr);

	// add all guide patterns to the set of transfrags so that I can have a "backbone" for each guide
	// I need this because there might be an uncompatible transfrag connecting the nodes in the guide
	for(int i=0;i<guidetrf.Count();i++) {
		CTransfrag *t=new CTransfrag(guidetrf[i].trf->nodes,guidetrf[i].trf->pattern,trthr);

		/*
		{ // DEBUG ONLY
			fprintf(stderr,"Add guidetrf with nodes:");
			for(int j=0;j<guidetrf[i].trf->nodes.Count();j++) fprintf(stderr," %d",guidetrf[i].trf->nodes[j]);
			//fprintf(stderr," and pattern: ");
			//printBitVec(guidetrf[i].trf->pattern);
			fprintf(stderr,"\n");
		}
		*/

		// I might not need to do this for the normal max_flow, but for the push I am counting on having only one transcript linking back to source from the first node
		t->pattern[0]=0;
		t->pattern[gno-1]=0;
		int *pos=gpos[edge(0,t->nodes[1],gno)];
		if(pos) t->pattern[*pos]=0;
		pos=gpos[edge(t->nodes[t->nodes.Count()-2],t->nodes.Last(),gno)];
		if(pos) t->pattern[*pos]=0;
		t->nodes.Pop();
		t->nodes.Shift();

		transfrag.Add(t);
	}


	// add edges between disconnected parent-child nodes
	GBitVec allpat(gno+edgeno);
	for(int t=0;t<transfrag.Count();t++) allpat=allpat | transfrag[t]->pattern;

	for(int i=1;i<gno-1;i++) { // for all nodes check if there is a connection to child
		CGraphnode *n=no2gnode[i];
		for(int c=0;c<n->child.Count();c++) {
			int *pos=gpos[edge(i,n->child[c],gno)];
			if(pos && !allpat[*pos]) {
				GVec<int> nodes;
				nodes.Add(i);
				nodes.Add(n->child[c]);
				GBitVec trpat(gno+edgeno);
				trpat[i]=1;
				trpat[n->child[c]]=1;
				trpat[*pos]=1;
				CTransfrag *t=new CTransfrag(nodes,trpat,trthr);
				transfrag.Add(t);
			}
		}
	}


	// sort transfrag with smallest being the one that has the most nodes, and ties are decided by the abundance (largest abundance first); last transfrags all have 1 node
	transfrag.Sort(trCmp);

	if(!fast) compatible.resize((1+transfrag.Count())*transfrag.Count()/2); // I might want to change this to gbitvec

	/*
	{ // DEBUG ONLY
		printTime(stderr);
		fprintf(stderr,"There are %d transfrags that remained\n",transfrag.Count());
	}
	*/

	GVec<int> incompletetrf; //remembers incomplete transfrags (the ones that don't have edges between two consecutive nodes

	// create compatibilities
	for(int t1=0;t1<transfrag.Count();t1++) { // transfrags are processed in increasing order -> important for the later considerations

		// update nodes
		int n1=transfrag[t1]->nodes.Count();

		if(n1>1) { // add transfrag to nodes' in and out; if a transfrag only has one node then it is node added to a node; I might want to change this for the computation of fpkm's
			bool incomplete = false;
			for(int n=0;n<n1;n++) { // for all nodes in transfrag

				if(n && n<transfrag[t1]->nodes.Count()-1) {// not first or last node
					// add t1 to in and out of node
					no2gnode[transfrag[t1]->nodes[n]]->trf.Add(t1);

					// check if transfrag t1 is incomplete between node[n-1] and node [n]
					int *pos=gpos[edge(transfrag[t1]->nodes[n-1],transfrag[t1]->nodes[n],gno)];
					if(!pos || !transfrag[t1]->pattern[*pos]) // there is no edge between node[n-1] and node[n]
						incomplete = assign_incomplete_trf_to_nodes(t1,transfrag[t1]->nodes[n-1],transfrag[t1]->nodes[n],no2gnode) or incomplete; 	// this ensures that I still have compatibilities when going through a certain path: !!! THIS IS NEW COMPARED TO PERL
				}
				else if(n) { // last but not first node
					// add t1 to in of node
					no2gnode[transfrag[t1]->nodes[n]]->trf.Add(t1);

					// check if transfrag t1 is incomplete between node[n-1] and node [n]
					int *pos=gpos[edge(transfrag[t1]->nodes[n-1],transfrag[t1]->nodes[n],gno)];
					if(!pos || !transfrag[t1]->pattern[*pos]) // there is no edge between node[n-1] and node[n]
						incomplete = assign_incomplete_trf_to_nodes(t1,transfrag[t1]->nodes[n-1],transfrag[t1]->nodes[n],no2gnode) or incomplete; 	// this ensures that I still have compatibilities when going through a certain path: !!! THIS IS NEW COMPARED TO PERL
				}
				else { // first node -> only add transfrag to out of node
					no2gnode[transfrag[t1]->nodes[n]]->trf.Add(t1);
				}
			}

			if(incomplete) incompletetrf.Add(t1);
			else transfrag[t1]->real=true;

		}
		/*
		else { // this transcript is included completely in node
			no2gnode[transfrag[t1]->nodes[0]]->frag+=transfrag[t1]->abundance;
		}
		*/

		if(!fast) {
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
		} // and if(!fast)
	} // end for(int t1=0;t1<transfrag.Count();t1++)

	// set source-to-child transfrag abundances: optional in order not to keep these abundances too low:
	// update the abundances of the transfrags coming in from source and going to a node that doesn't have other parents than source
	// * this part was removed to improve performance
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
	// */

	for(int t=0;t<incompletetrf.Count();t++)
		transfrag[incompletetrf[t]]->real=trf_real(incompletetrf[t],no2gnode,transfrag,gpos,gno);

}


void process_merge_transfrags(int gno,GPVec<CGraphnode>& no2gnode, GPVec<CMTransfrag>& mgt,GBitVec& compatible,GIntHash<int> &gpos) {

	// sort transfrag with smallest being the one that has the most nodes
	mgt.Sort(mgtrnodeCmp);

	for(int t1=1;t1<mgt.Count();t1++)
		for(int t2=0;t2<t1;t2++) {
			/*
			GBitVec p(mgt[t1]->transfrag->pattern);
			p[0]=0;
			p[gno-1]=0;
			int *pos=gpos[edge(0,mgt[t1]->transfrag->nodes[1],gno)];
			if(pos) p[pos]=0;
			if(mgt[t1]->transfrag->nodes.Count()>1) {
				*pos=gpos[edge(mgt[t1]->transfrag->nodes[mgt[t1]->transfrag->nodes.Count()-2],gno-1,gno)];
				if(pos) p[pos]=0;
			}
			*/

			if(!mgt[t1]->transfrag->real &&
					((mgt[t2]->transfrag->pattern & mgt[t1]->transfrag->pattern) == mgt[t1]->transfrag->pattern)) { // t1 is included in t2, and it's not a guide
				mgt[t2]->read.Add(mgt[t1]->read);
				mgt[t2]->transfrag->abundance+=mgt[t1]->transfrag->abundance;
				mgt.Exchange(t1,mgt.Count()-1);
				mgt.Delete(mgt.Count()-1);
			}
		}

	mgt.Sort(mgtrabundCmp); // sort transfrags from the most abundant to the least, with guides coming in first

	// create compatibilities
	for(int t1=0;t1<mgt.Count();t1++) { // transfrags are processed in increasing order -> important for the later considerations

		// update nodes
		int n1=mgt[t1]->transfrag->nodes.Count();

		for(int n=0;n<n1;n++) { // for all nodes in transfrag
			no2gnode[mgt[t1]->transfrag->nodes[n]]->trf.Add(t1);
		}

		// set compatibilities
		// add t1 to t1 compatibility
		bool comp=true;
		compatible[comptbl_pos(t1,t1,mgt.Count())]=1;
		for(int t2=t1+1;t2<mgt.Count();t2++) {
			// here check compatibility between t1 and t2;
			int n2=mgt[t2]->transfrag->nodes.Count();
			int i1=0;
			int i2=0;
			comp=true;
			while(i1<n1 && i2<n2) {
				if(mgt[t1]->transfrag->nodes[i1]==mgt[t2]->transfrag->nodes[i2]) {
					if(i1==n1-1 || i2==n2-1) { // one transcript finishes -> no need to check anymore
						i1=n1;i2=n2;
					}
					else { // advance the smallest one
						if(mgt[t1]->transfrag->nodes[i1+1]<mgt[t2]->transfrag->nodes[i2+1]) i1++;
						else i2++;
					}
				}
				else if(mgt[t1]->transfrag->nodes[i1]<mgt[t2]->transfrag->nodes[i2]) {
					i1++;
					if(conflict(i1,mgt[t2]->transfrag->nodes[i2],mgt[t1]->transfrag->nodes,n1,no2gnode,mgt[t1]->transfrag->pattern,gno,gpos)) {
						comp=false;
						break;
					}
				}
				else {
					i2++;
					if(conflict(i2,mgt[t1]->transfrag->nodes[i1],mgt[t2]->transfrag->nodes,n2,no2gnode,mgt[t2]->transfrag->pattern,gno,gpos)) {
						comp=false;
						break;
					}
				}
			}
			if(comp) compatible[comptbl_pos(t1,t2,mgt.Count())]=1;

		} // end for(int t2=t1+1;t2<transfrag.Count();t2++)
	} // end for(int t1=0;t1<transfrag.Count();t1++)

}



void extend_path_left(int n, int gno, GVec<int>& keept, GBitVec &unionpat,GPVec<CGraphnode>& no2gnode,
		GPVec<CMTransfrag>& mgt,GIntHash<int> &gpos,GVec<int>& leftpath) {

	CGraphnode *node=no2gnode[n];

	int maxp=-1;
	float maxabund=0;

	/*
	{ // DEBUG ONLY
		fprintf(stderr," continue left node %d with %d parents and %d compatible transcripts:",n,node->parent.Count(),keept.Count());
		//for(int i=0;i<keept.Count();i++) fprintf(stderr," %d",keept[i]);
		fprintf(stderr,"\n");
	}
	*/

	// this version takes the most abundant path only
	for(int i=0;i<node->parent.Count();i++) {
		int p=node->parent[i];
		if(p) { // parent not source
			if(unionpat[p]) { // there might be transcripts on the path: this will not be true for fuzzy continuations
				float sumabund=0;
				for(int t=0;t<keept.Count();t++) if(mgt[keept[t]]->transfrag->pattern[p]) { // parent is on transcript
					if(mgt[keept[t]]->transfrag->nodes[mgt[keept[t]]->nl]==p) {
						sumabund+=mgt[keept[t]]->transfrag->abundance;
					}
					else {
						int *pos=gpos[edge(p,n,gno)];
						if(pos && mgt[keept[t]]->transfrag->pattern[*pos]) {
							sumabund+=mgt[keept[t]]->transfrag->abundance;
						}
					}
				}
				//else if(mgt[keept[t]]->transfrag->nodes.Last()<p) sumabund+=mgt[keept[t]]->transfrag->abundance; this transcript shouldn't contribute to the path chosen
				if(sumabund>maxabund) {
					maxabund=sumabund;
					maxp=p;
				}
			}
		}
	}

	if(maxp>-1) { // I found a non-fuzzy way to continue
		GVec<int> ekeept; // extend keept;
		for(int t=0;t<keept.Count();t++) {
			if(mgt[keept[t]]->transfrag->pattern[maxp]) { // parent is on transcript
				if(mgt[keept[t]]->transfrag->nodes[mgt[keept[t]]->nl]==maxp) ekeept.Add(keept[t]);
				else {
					int *pos=gpos[edge(maxp,n,gno)];
					if(pos && mgt[keept[t]]->transfrag->pattern[*pos]) ekeept.Add(keept[t]);
				}
			}
			else if(mgt[keept[t]]->transfrag->nodes.Last()<maxp) ekeept.Add(keept[t]); // keep upstream transcript just in case it might turn useful
			else if(mgt[keept[t]]->transfrag->nodes[0]<=maxp && mgt[keept[t]]->transfrag->nodes[mgt[keept[t]]->nf]>maxp){ // the transcript is fuzzy
				ekeept.Add(keept[t]);
			}
		}
		if(ekeept.Count()) { // if I can extend path
			leftpath.Add(maxp);
			extend_path_left(maxp,gno,ekeept,unionpat,no2gnode,mgt,gpos,leftpath);
		}
	}
	else { // the path either stops here or it's fuzzy
		int min=n;
		GVec<int> path;
		for(int t=0;t<keept.Count();t++) { // we need to check if transcripts are fuzzy because they might end further upstream
			if(mgt[keept[t]]->transfrag->nodes[mgt[keept[t]]->nf] >= n && mgt[keept[t]]->transfrag->nodes[0]<min) { // transcript is fuzzy and starts before min
				for(int j=0;j<mgt[keept[t]]->transfrag->nodes.Count();j++) {
					if(mgt[keept[t]]->transfrag->nodes[j]==min) break;
					path.Add(mgt[keept[t]]->transfrag->nodes[j]);
				}
				min=mgt[keept[t]]->transfrag->nodes[0];
			}
		}
		for(int i=path.Count()-1;i>=0;i--) leftpath.Add(path[i]);
	}
}

/*   // needs more work
void extend_path_left_EM(int n, int gno, GVec<int>& keept, GBitVec &unionpat,GPVec<CGraphnode>& no2gnode,
		GPVec<CMTransfrag>& mgt,GIntHash<int> &gpos,GVec<int>& leftpath,int& start) {

	CGraphnode *node=no2gnode[n];

	int maxp=-1;
	float maxabund=0;

	// this version takes the most abundant path only
	for(int i=0;i<node->parent.Count();i++) {
		int p=node->parent[i];
		if(p) { // parent not source
			if(unionpat[p]) { // there might be transcripts on the path: this will not be true for fuzzy continuations
				float sumabund=0;
				for(int t=0;t<keept.Count();t++) if(mgt[keept[t]]->transfrag->pattern[p]) { // parent is on transcript
					if(mgt[keept[t]]->transfrag->nodes[mgt[keept[t]]->nl]==p) {
						sumabund+=mgt[keept[t]]->transfrag->abundance;
					}
					else {
						int *pos=gpos[edge(p,n,gno)];
						if(pos && mgt[keept[t]]->transfrag->pattern[*pos]) {
							sumabund+=mgt[keept[t]]->transfrag->abundance;
						}
					}
				}
				//else if(mgt[keept[t]]->transfrag->nodes.Last()<p) sumabund+=mgt[keept[t]]->transfrag->abundance; this transcript shouldn't contribute to the path chosen
				if(sumabund>maxabund) {
					maxabund=sumabund;
					maxp=p;
				}
			}
		}
	}


	if(maxp>-1) { // I found a non-fuzzy way to continue
		GVec<int> ekeept; // extend keept;
		for(int t=0;t<keept.Count();t++) {
			if(mgt[keept[t]]->transfrag->pattern[maxp]) { // parent is on transcript
				if(mgt[keept[t]]->transfrag->nodes[mgt[keept[t]]->nl]==maxp) ekeept.Add(keept[t]);
				else {
					int *pos=gpos[edge(maxp,n,gno)];
					if(pos && mgt[keept[t]]->transfrag->pattern[*pos]) ekeept.Add(keept[t]);
				}
			}
			else if(mgt[keept[t]]->transfrag->nodes.Last()<maxp) ekeept.Add(keept[t]); // keep upstream transcript just in case it might turn useful
			else if(mgt[keept[t]]->transfrag->nodes[0]<=maxp && mgt[keept[t]]->transfrag->nodes[mgt[keept[t]]->nf]>maxp){ // the transcript is fuzzy
				ekeept.Add(keept[t]);
			}
		}
		if(ekeept.Count()) { // if I can extend path
			leftpath.Add(maxp);
			extend_path_left(maxp,gno,ekeept,unionpat,no2gnode,mgt,gpos,leftpath,start);
		}
	}
	else { // the path either stops here or it's fuzzy
		start=leftpath.Count()-1;
		int min=n;
		GVec<int> path;
		for(int t=0;t<keept.Count();t++) { // we need to check if transcripts are fuzzy because they might end further upstream
			if(mgt[keept[t]]->transfrag->nodes[mgt[keept[t]]->nf] >= n && mgt[keept[t]]->transfrag->nodes[0]<min) { // transcript is fuzzy and starts before min
				for(int j=0;j<mgt[keept[t]]->transfrag->nodes.Count();j++) {
					if(mgt[keept[t]]->transfrag->nodes[j]==min) break;
					path.Add(mgt[keept[t]]->transfrag->nodes[j]);
				}
				min=mgt[keept[t]]->transfrag->nodes[0];
			}
		}
		for(int i=path.Count()-1;i>=0;i--) leftpath.Add(path[i]);
	}
}
*/

CPrediction* store_merge_prediction(GVec<int>& alltr,GPVec<CMTransfrag>& mgt,GVec<int>& path,
		GPVec<CGraphnode>& no2gnode,int strand,int& geneno,bool& first,GList<CReadAln>& readlist,GPVec<GffObj>& guides) {

	/*
	{ // DEBUG ONLY
		fprintf(stderr,"  store prediction:");
		for(int i=0;i<path.Count();i++) fprintf(stderr," %d",path[i]);
		fprintf(stderr,"\n");
	}
	*/

	GffObj *t=NULL;
	GStr name;

	// first compute coverage
	float cov=0;
	for(int i=0;i<alltr.Count();i++) {
		int a=alltr[i];
		cov+=mgt[alltr[i]]->transfrag->abundance;
		//fprintf(stderr,"Add transfrag %d(%d %d) with cov=%f\n",alltr[i],mgt[alltr[i]]->nf,mgt[alltr[i]]->nl-mgt[alltr[i]]->transfrag->nodes.Count()+1,mgt[alltr[i]]->transfrag->abundance);
		mgt[alltr[i]]->transfrag->abundance=0;
		if(enableNames || (alltr.Count()==1 && mgt[alltr[0]]->transfrag->real)) for(int j=0;j<mgt[a]->read.Count();j++) {
			int r=mgt[a]->read[j];
			if(readlist[r]->tinfo->g == -1) {
				int fidx=1+readlist[r]->tinfo->fileidx;
				GStr fid(fidx);
				if(!name.is_empty()) {
					name+=", ";
				}
				name+=fid;
				name+=':';
				name+=readlist[r]->tinfo->name;
			}
			else {
				t=guides[readlist[r]->tinfo->g];
				if(!enableNames) break;
			}
		}
	}

	int len=0;
	CGraphnode *prevnode=NULL;
	GVec<GSeg> exons;

	bool firstex=true;

	for(int i=0;i<path.Count();i++) {

		CGraphnode *node=no2gnode[path[i]];

		uint nodestart=node->start;
		uint nodeend=node->end;

		if(!prevnode || firstex || node->start>prevnode->end+1) { // this is a new exon
			GSeg exon(nodestart,nodeend);
			exons.Add(exon);
			firstex=false;
		}
		else if(!firstex) exons.Last().end=nodeend;

		len+=nodeend-nodestart+1;

		prevnode=node;
	}

	char sign='-';
	if(strand) { sign='+';}
	if(first) { geneno++;}

	CPrediction *p=new CPrediction(geneno-1, t, exons[0].start, exons.Last().end, cov, sign, len);
	p->exons=exons;
	if(enableNames) p->mergename=name;
	first=false;

	return(p);
}

// used in merge_transfrags
CPrediction* store_merge_prediction(float cov,GVec<int>& alltr,GPVec<CMTransfrag>& mgt,GVec<int>& path,
		GPVec<CGraphnode>& no2gnode,int strand,int& geneno,bool& first,GList<CReadAln>& readlist,GPVec<GffObj>& guides, int g) {

	/*
	{ // DEBUG ONLY
		fprintf(stderr,"  store prediction with cov=%f:",cov);
		for(int i=0;i<path.Count();i++) fprintf(stderr," %d",path[i]);
		fprintf(stderr,"\n");
	}
	*/


	GffObj *t=NULL;
	if(g>-1) t=guides[g];
	GStr name;

	if(enableNames) for(int i=0;i<alltr.Count();i++) {
		int a=alltr[i];
		for(int j=0;j<mgt[a]->read.Count();j++) {
			int r=mgt[a]->read[j];
			if(!name.is_empty()) {
				name+=", ";
			}
			if(readlist[r]->tinfo->g == -1) {
				int fidx=1+readlist[r]->tinfo->fileidx;
				GStr fid(fidx);
				name+=fid;
				name+=':';
				name+=readlist[r]->tinfo->name;
			}
			else {
				name+=t->getID();
			}
		}
	}

	int len=0;
	CGraphnode *prevnode=NULL;
	GVec<GSeg> exons;

	bool firstex=true;

	for(int i=0;i<path.Count();i++) {

		CGraphnode *node=no2gnode[path[i]];

		uint nodestart=node->start;
		uint nodeend=node->end;

		if(!prevnode || firstex || node->start>prevnode->end+1) { // this is a new exon
			GSeg exon(nodestart,nodeend);
			exons.Add(exon);
			firstex=false;
		}
		else if(!firstex) exons.Last().end=nodeend;

		len+=nodeend-nodestart+1;

		prevnode=node;
	}

	char sign='-';
	if(strand) { sign='+';}
	if(first) { geneno++;}

	CPrediction *p=new CPrediction(geneno-1, t, exons[0].start, exons.Last().end, cov/len, sign, len);
	p->exons=exons;
	if(enableNames) p->mergename=name;
	first=false;

	return(p);
}

/*   // needs more work
int store_merge_prediction_EM(GList<CPrediction>& pred,GVec<int>& alltr,GPVec<CMTransfrag>& mgt,GVec<int>& path,
		GPVec<CGraphnode>& no2gnode,int strand,int geneno,bool& first,GList<CReadAln>& readlist,GPVec<GffObj>& guides,float cov) {


	{ // DEBUG ONLY
		fprintf(stderr,"  geneno=%d first=",geneno);
		if(first) fprintf(stderr,"1\n");
		else fprintf(stderr,"0\n");
		fprintf(stderr,"  store prediction:");
		for(int i=0;i<path.Count();i++) fprintf(stderr," %d",path[i]);
		fprintf(stderr,"\n");
	}


	GffObj *t=NULL;
	GStr name;

	// first compute coverage
	//float cov=0;
	for(int i=0;i<alltr.Count();i++) {
		int a=alltr[i];
		//cov+=mgt[alltr[0]]->transfrag->abundance;
		if(enableNames || (alltr.Count()==1 && mgt[alltr[0]]->transfrag->real)) for(int j=0;j<mgt[a]->read.Count();j++) {
			int r=mgt[a]->read[j];
			if(readlist[r]->tinfo->g == -1) {
				int fidx=1+readlist[r]->tinfo->fileidx;
				GStr fid(fidx);
				if(!name.is_empty()) {
					name+=", ";
				}
				name+=fid;
				name+=':';
				name+=readlist[r]->tinfo->name;
			}
			else {
				t=guides[readlist[r]->tinfo->g];
				if(!enableNames) break;
			}
		}
	}

	int len=0;
	CGraphnode *prevnode=NULL;
	GVec<GSeg> exons;

	bool firstex=true;

	for(int i=0;i<path.Count();i++) {

		CGraphnode *node=no2gnode[path[i]];

		uint nodestart=node->start;
		uint nodeend=node->end;

		if(!prevnode || firstex || node->start>prevnode->end+1) { // this is a new exon
			GSeg exon(nodestart,nodeend);
			exons.Add(exon);
			firstex=false;
		}
		else if(!firstex) exons.Last().end=nodeend;

		len+=nodeend-nodestart+1;

		prevnode=node;
	}

	char sign='-';
	if(strand) { sign='+';}
	if(first) { geneno++;}

	CPrediction *p=new CPrediction(geneno-1, t, exons[0].start, exons.Last().end, cov, sign, len);
	p->exons=exons;
	if(enableNames) p->mergename=name;
	pred.Add(p);
	first=false;

	return(geneno);
}


// this version only extends with the most abundant path
void extend_path_right(int n, int gno, GVec<int>& keept,GVec<int>& path, GVec<int>& uniont, GBitVec &unionpat,
		GPVec<CGraphnode>& no2gnode,GPVec<CMTransfrag>& mgt,GIntHash<int> &gpos,GVec<int>& leftpath,
		int& geneno,int strand,GVec<CMPrediction>& pred,GList<CReadAln>& readlist,GPVec<GffObj>& guides,bool& first) {

	CGraphnode *node=no2gnode[n];


	{ // DEBUG ONLY
		fprintf(stderr," continue node %d with %d children and %d compatible transcripts:",n,node->child.Count(),keept.Count());
		//for(int i=0;i<keept.Count();i++) fprintf(stderr," %d",keept[i]);
		fprintf(stderr,"\n");
	}


	int maxc=-1;
	float maxabund=0;

	for(int i=0;i<node->child.Count();i++) {
		int c=node->child[i];
		if(c<gno-1) { // child is no sink
			if(unionpat[c]) { // there might be transcripts on the path
				float sumabund=0;
				for(int t=0;t<keept.Count();t++) if(mgt[keept[t]]->transfrag->pattern[c]) { // child is on transcript or transcript starts further downstream
					if(mgt[keept[t]]->transfrag->nodes[mgt[keept[t]]->nf]==c) {
						sumabund+=mgt[keept[t]]->transfrag->abundance;
					}
					else {
						int *pos=gpos[edge(n,c,gno)];
						if(pos && mgt[keept[t]]->transfrag->pattern[*pos]) {
							sumabund+=mgt[keept[t]]->transfrag->abundance;
						}
					}
				}
				//else if(mgt[keept[t]]->transfrag->nodes[0]>c) sumabund+=mgt[keept[t]]->transfrag->abundance;
				if(sumabund>maxabund) { // if I can extend path
					maxabund=sumabund;
					maxc=c;
				}
			}
		}
	}

	if(maxc>-1) { // I found a non-fuzzy way to continue
		GVec<int> ekeept; // extend keept;
		for(int t=0;t<keept.Count();t++) {
			if(mgt[keept[t]]->transfrag->pattern[maxc]) { // child is on transcript or
				if(mgt[keept[t]]->transfrag->nodes[mgt[keept[t]]->nl]==maxc)  ekeept.Add(keept[t]);
				else {
					int *pos=gpos[edge(n,maxc,gno)];
					if(pos && mgt[keept[t]]->transfrag->pattern[*pos])  ekeept.Add(keept[t]);
				}
			}
			else if(mgt[keept[t]]->transfrag->nodes[0]>maxc) ekeept.Add(keept[t]); // transcript starts further downstream
			else if(mgt[keept[t]]->transfrag->nodes.Last()>=maxc && mgt[keept[t]]->transfrag->nodes[mgt[keept[t]]->nl]<maxc)
				ekeept.Add(keept[t]);
		}
		if(ekeept.Count()) { // if I can extend path
			path.Add(maxc);
			extend_path_right(maxc,gno,ekeept,path,uniont,unionpat,no2gnode,mgt,gpos,leftpath,geneno,strand,pred,readlist,guides,first);
		}
	}
	else { // only possible child is sink or fuzzy continuation

		GBitVec pathpat(mgt[uniont[0]]->transfrag->pattern);
		int nnodes=leftpath.Count()+mgt[uniont[0]]->transfrag->nodes.Count()+path.Count();
		GVec<int> printpath(nnodes); // nnodes just sets up an initial capacity -> our best estimate so far
		for(int j=leftpath.Count()-1;j>=0;j--) {
			printpath.Add(leftpath[j]);
			pathpat[leftpath[j]]=1;
			int *pos=NULL;
			if(j) pos=gpos[edge(leftpath[j],leftpath[j-1],gno)];
			else pos=gpos[edge(leftpath[0],mgt[uniont[0]]->transfrag->nodes[mgt[uniont[0]]->nf],gno)];
			if(pos) pathpat[*pos]=1;
		}
		for(int j=mgt[uniont[0]]->nf;j<=mgt[uniont[0]]->nl;j++) printpath.Add(mgt[uniont[0]]->transfrag->nodes[j]);

		for(int j=0;j<path.Count();j++) {
			printpath.Add(path[j]);
			pathpat[path[j]]=1;
			int *pos=NULL;
			if(j) pos=gpos[edge(path[j-1],path[j],gno)];
			else pos=gpos[edge(mgt[uniont[0]]->transfrag->nodes[mgt[uniont[0]]->nl],path[0],gno)];
			if(pos) pathpat[*pos]=1;
		}

		// extend path with fuzzy nodes if present
		path.Clear();
		int max=n;
		for(int t=0;t<keept.Count();t++) {
			if(mgt[keept[t]]->transfrag->nodes[mgt[keept[t]]->nl] <= n && mgt[keept[t]]->transfrag->nodes.Last()>max) { // transcript is fuzzy and starts before min
				for(int j=mgt[keept[t]]->transfrag->nodes.Count()-1;j>=0;j--) {
					if(mgt[keept[t]]->transfrag->nodes[j]==max) break;
					path.Add(mgt[keept[t]]->transfrag->nodes[j]);
				}
				max=mgt[keept[t]]->transfrag->nodes.Last();
			}
		}

		for(int j=path.Count()-1;j>=0;j--) {
			int *pos=NULL;
			if(j==path.Count()-1) pos=gpos[edge(printpath.Last(),path[j],gno)];
			else pos=gpos[edge(path[j],path[j+1],gno)];
			if(pos) pathpat[*pos]=1;
			printpath.Add(path[j]);
		}

		// only print transfrag if it's not real (and not printed already) or if it's not entirely equal to the path
		if(!mgt[uniont[0]]->transfrag->real ||
				(mgt[uniont[0]]->transfrag->pattern) != pathpat) {

			GVec<int> alltr(1,uniont[0]);
			for(int j=1;j<uniont.Count();j++) if((pathpat & mgt[uniont[j]]->transfrag->pattern) == mgt[uniont[j]]->transfrag->pattern) {
				//mgt[uniont[j]]->transfrag->abundance=0;
				alltr.Add(uniont[j]);
			}

			CPrediction *p=store_merge_prediction(alltr,mgt,printpath,no2gnode,strand,geneno,first,readlist,guides);
			CMPrediction mp(p);
			mp.b=pathpat;
			for(int i=0;i<printpath.Count();i++) {
				mp.b[printpath[i]]=0;
				if(i && no2gnode[printpath[i]]->start==no2gnode[printpath[i-1]]->end+1) {
					int *pos=gpos[edge(printpath[i-1],printpath[i],gno)];
					if(pos) mp.b[*pos]=0;
				}
			}
			pred.Add(mp);
		}

	}

}
*/


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

/*   // needs more work
float merge_max_flow(int gno,GVec<int>& path,GBitVec& istranscript,GPVec<CMTransfrag>& mgt,GPVec<CGraphnode>& no2gnode,
		GVec<int>& alltr) {

	float flux=0;
	int n=path.Count();
	GVec<float> *capacity=new GVec<float>[n]; // capacity of edges in network
	GVec<float> *flow=new GVec<float>[n]; // flow in network
	GVec<int> *link=new GVec<int>[n]; // for each node remembers it's neighbours
	GVec<int> pred; // this stores the augmenting path
	pred.Resize(n,-1);
	GVec<int> node2path;
	node2path.Resize(gno,-1);


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


	for(int i=0;i<n;i++) {
		node2path[path[i]]=i;
		capacity[i].Resize(n,0);
		flow[i].Resize(n,0);
	}

	// establish capacities in the network
	for(int i=0;i<n;i++) {
		int nt=no2gnode[path[i]]->trf.Count();
		//fprintf(stderr,"Node %d has %d transcripts\n",path[i],nt);
		for(int j=0;j<nt;j++) {
			int t=no2gnode[path[i]]->trf[j];
			if(mgt[t]->transfrag->abundance && istranscript[t]) {
				if(mgt[t]->transfrag->nodes[mgt[t]->nf]==path[i]) { // transfrag fuzzy starts at this node
					int n1=i;
					int n2=node2path[mgt[t]->transfrag->nodes[mgt[t]->nl]];
					if(!no2gnode[path[i]]->rate) n1=0;
					if(!no2gnode[mgt[t]->transfrag->nodes[mgt[t]->nl]]->rate) n2=n-1;
					//fprintf(stderr,"path[%d]=%d t=%d n1=%d n2=%d ",i,path[i],t,n1,n2);
					if(!capacity[n1][n2]) { // haven't seen this link before
						link[n1].Add(n2);
						link[n2].Add(n1);
					}
					capacity[n1][n2]+=mgt[t]->transfrag->abundance;
					//fprintf(stderr,"capacity[%d][%d]=%f after adding transfrag[%d]->abundance=%f\n",n1,n2,capacity[n1][n2],t,mgt[t]->transfrag->abundance);
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


	{ // DEBUG ONLY
		fprintf(stderr,"Flow:");
		for(int n1=0;n1<n;n1++)
			for(int n2=n1+1;n2<n;n2++) if(flow[n1][n2]) fprintf(stderr," [%d][%d]=%f",n1,n2,flow[n1][n2]);
		fprintf(stderr,"\n");
	}


	float fluxcov=0;

	// adjust transfrag abundances
	for(int i=0;i<n;i++) {
		int nt=no2gnode[path[i]]->trf.Count();
		float sumout=0;
		int pos=-1;
		for(int j=0;j<nt;j++) {
			int t=no2gnode[path[i]]->trf[j];
			if(istranscript[t] && mgt[t]->transfrag->abundance) {
				if(mgt[t]->transfrag->nodes[mgt[t]->nf]==path[i]) { // transfrag starts at this node
					int n1=i;
					int n2=node2path[mgt[t]->transfrag->nodes[mgt[t]->nl]];
					if(!no2gnode[path[i]]->rate) n1=0;
					if(!no2gnode[mgt[t]->transfrag->nodes[mgt[t]->nl]]->rate) n2=n-1;
					if(flow[n1][n2]>0) {
						alltr.Add(t);
						if(flow[n1][n2]<mgt[t]->transfrag->abundance) {
							if(!i) sumout+=flow[n1][n2];
							fluxcov+=flow[n1][n2];
							mgt[t]->transfrag->abundance-=flow[n1][n2];
							if(mgt[t]->transfrag->abundance<epsilon) mgt[t]->transfrag->abundance=0;
							flow[n1][n2]=0;
						}
						else {
							if(!i) sumout+=mgt[t]->transfrag->abundance;
							fluxcov+=mgt[t]->transfrag->abundance;
							flow[n1][n2]-=mgt[t]->transfrag->abundance;
							mgt[t]->transfrag->abundance=0;
						}
					}
				}
				else if(!i && mgt[t]->transfrag->nodes[mgt[t]->nl]==path[i]) pos=j;
			}
		}

		if(!i && pos>-1) { // this is first node -> adjust entering transfrag
			int t=no2gnode[path[i]]->trf[pos];
			float val=sumout/no2gnode[path[i]]->rate;
			mgt[t]->transfrag->abundance-=val;
			if(mgt[t]->transfrag->abundance<epsilon) mgt[t]->transfrag->abundance=0;
		}

	}

	// clean up
	delete [] capacity;
	delete [] flow;
	delete [] link;

	return(flux);
}


bool extend_path_right_EM(int n, int gno, GVec<int>& keept,GVec<int>& path, GVec<int>& uniont, GBitVec &unionpat,
		GPVec<CGraphnode>& no2gnode,GPVec<CMTransfrag>& mgt,GIntHash<int> &gpos,GVec<int>& leftpath,int start,
		int& geneno,int strand,GList<CPrediction>& pred,GList<CReadAln>& readlist,GPVec<GffObj>& guides,bool& first) {

	bool result=false;

	CGraphnode *node=no2gnode[n];


	{ // DEBUG ONLY
		fprintf(stderr," continue node %d with %d children and %d compatible transcripts:",n,node->child.Count(),keept.Count());
		//for(int i=0;i<keept.Count();i++) fprintf(stderr," %d",keept[i]);
		fprintf(stderr,"\n");
	}


	int maxc=-1;
	float maxabund=0;

	for(int i=0;i<node->child.Count();i++) {
		int c=node->child[i];
		if(c<gno-1) { // child is no sink
			if(unionpat[c]) { // there might be transcripts on the path
				float sumabund=0;
				for(int t=0;t<keept.Count();t++) if(mgt[keept[t]]->transfrag->pattern[c]) { // child is on transcript or transcript starts further downstream
					if(mgt[keept[t]]->transfrag->nodes[mgt[keept[t]]->nf]==c) {
						sumabund+=mgt[keept[t]]->transfrag->abundance;
					}
					else {
						int *pos=gpos[edge(n,c,gno)];
						if(pos && mgt[keept[t]]->transfrag->pattern[*pos]) {
							sumabund+=mgt[keept[t]]->transfrag->abundance;
						}
					}
				}
				//else if(mgt[keept[t]]->transfrag->nodes[0]>c) sumabund+=mgt[keept[t]]->transfrag->abundance;
				if(sumabund>maxabund) { // if I can extend path
					maxabund=sumabund;
					maxc=c;
				}
			}
		}
	}

	if(maxc>-1) { // I found a non-fuzzy way to continue
		GVec<int> ekeept; // extend keept;
		for(int t=0;t<keept.Count();t++) {
			if(mgt[keept[t]]->transfrag->pattern[maxc]) { // child is on transcript or
				if(mgt[keept[t]]->transfrag->nodes[mgt[keept[t]]->nf]==maxc)  ekeept.Add(keept[t]);
				else {
					int *pos=gpos[edge(n,maxc,gno)];
					if(pos && mgt[keept[t]]->transfrag->pattern[*pos])  ekeept.Add(keept[t]);
				}
			}
			else if(mgt[keept[t]]->transfrag->nodes[0]>maxc) ekeept.Add(keept[t]); // transcript starts further downstream
			else if(mgt[keept[t]]->transfrag->nodes.Last()>=maxc && mgt[keept[t]]->transfrag->nodes[mgt[keept[t]]->nl]<maxc)
				ekeept.Add(keept[t]);
		}
		if(ekeept.Count()) { // if I can extend path
			path.Add(maxc);
			if(extend_path_right(maxc,gno,ekeept,path,uniont,unionpat,no2gnode,mgt,gpos,leftpath,start,geneno,strand,pred,readlist,guides,first))
				result=true;
		}
	}
	else { // only possible child is sink or fuzzy continuation

		GBitVec pathpat(mgt[uniont[0]]->transfrag->pattern);
		int nnodes=leftpath.Count()+mgt[uniont[0]]->transfrag->nodes.Count()+path.Count();
		GVec<int> printpath(nnodes); // nnodes just sets up an initial capacity -> our best estimate so far
		for(int j=start;j>=0;j--) {
			printpath.Add(leftpath[j]);
			pathpat[leftpath[j]]=1;
			int *pos=NULL;
			if(j) pos=gpos[edge(leftpath[j],leftpath[j-1],gno)];
			else pos=gpos[edge(leftpath[0],mgt[uniont[0]]->transfrag->nodes[mgt[uniont[0]]->nf],gno)];
			if(pos) pathpat[*pos]=1;
		}
		for(int j=mgt[uniont[0]]->nf;j<=mgt[uniont[0]]->nl;j++) printpath.Add(mgt[uniont[0]]->transfrag->nodes[j]);

		for(int j=0;j<path.Count();j++) {
			printpath.Add(path[j]);
			pathpat[path[j]]=1;
			int *pos=NULL;
			if(j) pos=gpos[edge(path[j-1],path[j],gno)];
			else pos=gpos[edge(mgt[uniont[0]]->transfrag->nodes.Last(),path[0],gno)];
			if(pos) pathpat[*pos]=1;
		}

		// compute flow
		GBitVec istranscript(mgt.Count());
		for(int j=0;j<uniont.Count();j++) if((pathpat & mgt[uniont[j]]->transfrag->pattern) == mgt[uniont[j]]->transfrag->pattern) {
			istranscript[uniont[j]]=1;
			//fprintf(stderr,"Add t=%d\n",uniont[j]);
		}
		GVec<int> alltr;
		float fluxcov=merge_max_flow(gno,printpath,istranscript,mgt,no2gnode,alltr);
		//fprintf(stderr,"flux=%f\n",fluxcov);
		if(fluxcov) {

			// extend path with fuzzy nodes if present
			path.Clear();
			for(int j=leftpath.Count()-1;j>start;j--) path.Add(leftpath[j]);
			leftpath.Clear();
			leftpath.Add(path);
			leftpath.Add(printpath);

			path.Clear();
			int max=n;
			for(int t=0;t<keept.Count();t++) {
				if(mgt[keept[t]]->transfrag->nodes[mgt[keept[t]]->nl] <= n && mgt[keept[t]]->transfrag->nodes.Last()>max) { // transcript is fuzzy and starts before min
					for(int j=mgt[keept[t]]->transfrag->nodes.Count()-1;j>=0;j--) {
						if(mgt[keept[t]]->transfrag->nodes[j]==max) break;
						path.Add(mgt[keept[t]]->transfrag->nodes[j]);
					}
					max=mgt[keept[t]]->transfrag->nodes.Last();
				}
			}

			for(int j=path.Count()-1;j>=0;j--) {
				leftpath.Add(path[j]);
			}

			geneno=store_merge_prediction(pred,alltr,mgt,leftpath,no2gnode,strand,geneno,first,readlist,guides,fluxcov);

			if(!result && (mgt[uniont[0]]->transfrag->pattern == pathpat)) result=true;
		}
	}

	return(result);

}


int add_transfrag_to_prediction(GVec<int>& uniont, GBitVec &unionpat,int gno,GPVec<CGraphnode>& no2gnode, GPVec<CMTransfrag>& mgt,GIntHash<int> &gpos,
		int geneno,int strand,GVec<CMPrediction>& pred,GList<CReadAln>& readlist,GPVec<GffObj>& guides,bool& first) {

	GVec<int> leftpath; // initially this has no elements
	GVec<int> path;

	// first print transfrag if it's real so that I can give it all credit
	if(mgt[uniont[0]]->transfrag->real) {
		GVec<int> alltr(1,uniont[0]); // all transfrags used in creating path
		CPrediction *p=store_merge_prediction(alltr,mgt,mgt[uniont[0]]->transfrag->nodes,no2gnode,strand,geneno,first,readlist,guides);
		CMPrediction mp(p);
		mp.b=mgt[uniont[0]]->transfrag->pattern;
		for(int i=0;i<mgt[uniont[0]]->transfrag->nodes.Count();i++) {
			mp.b[mgt[uniont[0]]->transfrag->nodes[i]]=0;
			if(i && no2gnode[mgt[uniont[0]]->transfrag->nodes[i]]->start==no2gnode[mgt[uniont[0]]->transfrag->nodes[i-1]]->end+1) {
				int *pos=gpos[edge(mgt[uniont[0]]->transfrag->nodes[i-1],mgt[uniont[0]]->transfrag->nodes[i],gno)];
				if(pos) mp.b[*pos]=0;
			}
		}
		pred.Add(mp);
	}

	//extend_path_left(mgt[uniont[0]]->transfrag->nodes[0],gno,uniont,path,unionpat,no2gnode,mgt,gpos,leftpath);
	extend_path_left(mgt[uniont[0]]->transfrag->nodes[mgt[uniont[0]]->nf],gno,uniont,unionpat,no2gnode,mgt,gpos,leftpath);

	path.Clear();
	extend_path_right(mgt[uniont[0]]->transfrag->nodes[mgt[uniont[0]]->nl],gno,uniont,path,uniont,unionpat,no2gnode,mgt,gpos,
			leftpath,geneno,strand,pred,readlist,guides,first);

	return(geneno);
}

// needs more work
int add_transfrag_to_prediction_EM(GVec<int>& uniont, GBitVec &unionpat,int gno,GPVec<CGraphnode>& no2gnode, GPVec<CMTransfrag>& mgt,GIntHash<int> &gpos,
		int geneno,int strand,GList<CPrediction>& pred,GList<CReadAln>& readlist,GPVec<GffObj>& guides,bool& first) {

	int t=uniont[0];
	float abund=mgt[t]->transfrag->abundance;

	GVec<int> leftpath; // initially this has no elements
	GVec<int> path;

	bool printtr=false;

	while(mgt[t]->transfrag->abundance>isofrac*abund) {

		leftpath.Clear();
		path.Clear();
		int start=-1;
		extend_path_left(mgt[uniont[0]]->transfrag->nodes[mgt[uniont[0]]->nf],gno,uniont,unionpat,no2gnode,mgt,gpos,leftpath,start);

		path.Clear();
		int npred=pred.Count();
		printtr=printtr or extend_path_right(mgt[uniont[0]]->transfrag->nodes[mgt[uniont[0]]->nl],gno,uniont,path,uniont,unionpat,no2gnode,mgt,gpos,
			leftpath,start,geneno,strand,pred,readlist,guides,first);
		if(npred==pred.Count()) break;
	}

	if(!printtr && mgt[uniont[0]]->transfrag->real) {
		GVec<int> alltr(1,uniont[0]); // all transfrags used in creating path
		geneno=store_merge_prediction(pred,alltr,mgt,mgt[uniont[0]]->transfrag->nodes,no2gnode,strand,geneno,first,readlist,guides);

	}

	return(geneno);
}
*/

/*
int nodecapCmp(const pointer p1, const pointer p2) {
	CNodeCapacity *a=(CNodeCapacity*)p1;
	CNodeCapacity *b=(CNodeCapacity*)p2;
	if(a->perc > b->perc) return 1;
	if(a->perc < b->perc) return -1;
	return 0;
}
*/

int mpredCmp(const pointer p1, const pointer p2) {
	CMPrediction *a=(CMPrediction*)p1;
	CMPrediction *b=(CMPrediction*)p2;

	if(a->p->exons.Count() < b->p->exons.Count()) return 1;     // more exons come first
	if(a->p->exons.Count() > b->p->exons.Count()) return -1;

	if(a->p->t_eq==NULL && b->p->t_eq!=NULL) return 1;     // known gene come first
	if(a->p->t_eq!=NULL && b->p->t_eq==NULL) return -1;

	if(a->p->cov < b->p->cov) return 1;      // more coverage come first
	if(a->p->cov > b->p->cov) return -1;

	// add something that guarantees same order
	if(a->p->start < b->p->start) return -1; // order based on start
	if(a->p->start > b->p->start) return 1;

	// same start

	// same number of exons
	int i=0;
	uint a1=0;
	uint b1=0;
	while(i<a->p->exons.Count()) {
		if(a->p->exons[i].start<b->p->exons[i].start) {
			a1=a->p->exons[i].start;
			b1=b->p->exons[i].start;
			break;
		}
		if(a->p->exons[i].end<b->p->exons[i].end) {
			a1=a->p->exons[i].end;
			b1=b->p->exons[i].end;
			break;
		}
		i++;
	}
	if(a1) {
		if(a1<b1) return -1; // the one with the first exon comes first
		if(a1>b1) return 1;
	}
	return 0;
}

bool merge_onpath(CTransfrag *t,int nf, int nl,GBitVec& pathpattern,int mini,int maxi,GPVec<CGraphnode>& no2gnode,int gno,
		GIntHash<int>& gpos) {

	if(t->nodes[nf]<mini) // mini can be reached through transcript
		 if(!no2gnode[mini]->parentpat[t->nodes[nf]])	return false;

	if(t->nodes[nl]>maxi) // from maxi I can reach end of transcript
	    if(!no2gnode[maxi]->childpat[t->nodes[nl]]) return false;

	int first=1;

	for(int i=nf;i<=nl;i++) {
		if(t->nodes[i]>=mini && t->nodes[i]<=maxi) {
	      if(!pathpattern[t->nodes[i]]) return false;
    	  int *pos=NULL;
    	  if(i>nf) pos=gpos[edge(t->nodes[i-1],t->nodes[i],gno)];
	      if(first) {
	    	  first=0;
	    	  if(pos && t->nodes[i]>mini && t->pattern[*pos])
	    		  return false;
	    	  if(i>nf && !no2gnode[mini]->parentpat[t->nodes[i-1]]) return false; // I can not reach mini from previous node in transfrag
	      }
	      else if(pos && t->pattern[*pos] && !pathpattern[*pos]) return false;
	    }
	    if(t->nodes[i]>maxi) {
	    	int *pos=NULL;
	    	if(i>nf) pos=gpos[edge(t->nodes[i-1],t->nodes[i],gno)];
	    	if(pos && t->nodes[i-1]<maxi && t->pattern[*pos])
	    		return false;
	    	if(!no2gnode[maxi]->childpat[t->nodes[i]]) return false; // I can not reach this node from maxi
	    	break;
	    }
	}

	return true;
}


bool merge_topath(CTransfrag *t,int nf, int nl,GBitVec& pathpattern,int mini,int maxi,GPVec<CGraphnode>& no2gnode,int gno,
		GIntHash<int>& gpos) {


	int mint=nf; // first node on path
	while(mint<=nl && !pathpattern[t->nodes[mint]]) mint++;
	if(mint>nl) return false;

	// weird possiblities:
	// 1. path = ooo---ooo----
	//    t    =     ooooo----
	//
	// 2. path =     ooooo---
	//    t    = ooo---ooo---

	bool firstone=false;
	bool firsttwo=false;
	if(mint>nf) { // first nodes in t are not on path
		if(mini<t->nodes[nf]) { // check case 1
			int i=mint;
			while(i>nf) {
				if(t->nodes[i-1]==t->nodes[i]-1 && no2gnode[t->nodes[i-1]]->end+1==no2gnode[t->nodes[i]]->start) i--;
				else return false;
			}
			firstone=true;
		}
		else if(t->nodes[nf]<mini) { // check case 2
			int i=t->nodes[mint];
			while(i>mini) {
				if(pathpattern[i-1] && no2gnode[i-1]->end+1==no2gnode[i]->start) i--;
				else return false;
			}
			firsttwo=true;
		}
	}

	int maxt=nl; // last node on path
	while(maxt>mint && !pathpattern[t->nodes[maxt]]) maxt--;

	// weird possiblities:
	// 1. path = ---ooo----ooo
	//    t    = ---ooooo
	//
	// 2. path = ---oooooo
	//    t    = ---oooo----ooo

	bool lastone=false;
	bool lasttwo=false;
	if(maxt<nl) { // last nodes in t are not on path
		if(maxi>t->nodes[nl]) { // check case 1
			int i=maxt;
			while(i<nl) {
				if(t->nodes[i]==t->nodes[i+1]-1 && no2gnode[t->nodes[i]]->end+1==no2gnode[t->nodes[i+1]]->start) i++;
				else return false;
			}
			lastone=true;
		}
		else if(t->nodes[nl]>maxi) { // check case 2
			int i=t->nodes[maxt];
			while(i<maxi) {
				if(pathpattern[i+1] && no2gnode[i]->end+1==no2gnode[i+1]->start) i++;
				else return false;
			}
			lasttwo=true;
		}
	}

	// check if all nodes and introns between mint and maxt are in pathpattern
	for(int i=mint;i<maxt;i++) {
		if(!pathpattern[t->nodes[i]]) return false;
		int *pos=gpos[edge(t->nodes[i],t->nodes[i+1],gno)];
		if(pos && t->pattern[*pos] && !pathpattern[*pos]) return false;
	}

	if(firstone) { // modify start of t
		for(int i=nf;i<mint;i++) {
			t->pattern[t->nodes[i]]=0;
			int *pos=gpos[edge(t->nodes[i],t->nodes[i+1],gno)];
			if(pos) t->pattern[*pos]=0;
		}
	}
	else if(firsttwo) { // modify start of pathpattern
		for(int i=mini;i<t->nodes[mint];i++) {
			pathpattern[i]=0;
			int *pos=gpos[edge(i,i+1,gno)];
			if(pos) pathpattern[*pos]=0;
		}
	}

	if(lastone) { // modify end of t
		for(int i=maxt;i<nl;i++) {
			t->pattern[t->nodes[i+1]]=0;
			int *pos=gpos[edge(t->nodes[i],t->nodes[i+1],gno)];
			if(pos) t->pattern[*pos]=0;
		}
	}
	else if(lasttwo) { // modify start of pathpattern
		for(int i=t->nodes[maxt];i<maxi;i++) {
			pathpattern[i+1]=0;
			int *pos=gpos[edge(i,i+1,gno)];
			if(pos) pathpattern[*pos]=0;
		}
	}

	return true;
}


void add_transfrag_to_path(int t,GBitVec& tforlater,int& ntforlater, GPVec<CMTransfrag>& mgt,GBitVec& pathpat,int& min,int &max,
		GPVec<CGraphnode>& no2gnode,int gno,GIntHash<int>& gpos,GVec<int>& alltr,float &cov){

	//if(merge_onpath(mgt[t]->transfrag,mgt[t]->nf,mgt[t]->nl,pathpat,min,max,no2gnode,gno,gpos)) { // I can add the transfrag if it's on path
	if(merge_topath(mgt[t]->transfrag,mgt[t]->nf,mgt[t]->nl,pathpat,min,max,no2gnode,gno,gpos)) { // I can add the transfrag if it's on path

		//fprintf(stderr,"  Add t=%d\n",t);

		pathpat = pathpat | mgt[t]->transfrag->pattern;
		cov+=mgt[t]->transfrag->abundance*mgt[t]->len;
		mgt[t]->transfrag->abundance=0;
		alltr.Add(t);
		int prevmin=min;
		int prevmax=max;
		if(mgt[t]->transfrag->nodes[mgt[t]->nf]<min) min=mgt[t]->transfrag->nodes[mgt[t]->nf];
		if(mgt[t]->transfrag->nodes[mgt[t]->nl]>max) max=mgt[t]->transfrag->nodes[mgt[t]->nl];
		if(min<prevmin) {
			int i=mgt[t]->nf;
			while(ntforlater && mgt[t]->transfrag->nodes[i]<prevmin) { // for all nodes in transfrags that come before prevmin
				CGraphnode *node=no2gnode[mgt[t]->transfrag->nodes[i]];
				int nn=node->trf.Count();
				for(int j=0;j<nn;j++){
					int tj=node->trf[j];
					if(tforlater[tj]) {
						tforlater[tj]=0;
						ntforlater--;
						add_transfrag_to_path(tj,tforlater,ntforlater,mgt,pathpat,min,max,no2gnode,gno,gpos,alltr,cov);
					}
				}
				i++;
			}
		}
		if(mgt[t]->transfrag->nodes[mgt[t]->nl]>prevmax) { // my max could have been updated -> I do not need to check everything!
			int i=mgt[t]->nl;
			while(ntforlater && mgt[t]->transfrag->nodes[i]>prevmax) { // for all nodes in transfrags that come after prevmax
				CGraphnode *node=no2gnode[mgt[t]->transfrag->nodes[i]];
				int nn=node->trf.Count();
				for(int j=0;j<nn;j++){
					int tj=node->trf[j];
					if(tforlater[tj]) {
						tforlater[tj]=0;
						ntforlater--;
						add_transfrag_to_path(tj,tforlater,ntforlater,mgt,pathpat,min,max,no2gnode,gno,gpos,alltr,cov);
						if(!ntforlater) break;
					}
				}
				i--;
			}
		}
	}

}

// compute overlap of a small prediction to a larger one with the same exon structure, except for intron retention
int compute_ovlp(CPrediction *small,CPrediction *big)
{
	int len=0;
	int i=0;
	int ovp=0;
	while(i<big->exons.Count() && !(ovp=small->exons[0].overlapLen(big->exons[i].start,big->exons[i].end))) i++;

	if(ovp) { // if ovp is zero than there is no way I will find anything overlapping
		len+=ovp;
		ovp=0;
		int j=big->exons.Count()-1;
		while(j>i && !(ovp=small->exons.Last().overlapLen(big->exons[j].start,big->exons[j].end))) j--;
		len+=ovp;
		for(int k=i+1; k<j;k++) len+=big->exons[k].len();
	}

	return len;
}

bool has_retained_intron(CMPrediction &p1,CMPrediction &p2,int gno,GIntHash<int> &gpos) { // this assumes that p2.b & p1.b == p2.b

	for(int i=1;i<p1.nodes.Count();i++) {
		if(p2.pat[p1.nodes[i-1]] && p2.pat[p1.nodes[i]]) { // both nodes from p1 are also present in p2
			int *pos=gpos[edge(p1.nodes[i-1],p1.nodes[i],gno)];
			if(pos && p1.b[*pos] && !p2.b[*pos]) return true; // there is an intron from n[i-1] to n[i] in p1 but not in p2
		}
	}

	return false;
}

bool overlaps_one_exon_only(CMPrediction &p1,CMPrediction &p2) { // this assumes that p2 has only one exon

	bool exonovlp=false;
	bool gap=false;
	for(int i=0;i<p2.nodes.Count();i++) {
		if(p1.pat[p2.nodes[i]]) { // node in p2 also in p1
			if(gap) return false;
			if(!exonovlp) exonovlp=true;
		}
		else if(exonovlp) gap=true;
	}
	return exonovlp;
}

int merge_transfrags(int gno,GPVec<CGraphnode>& no2gnode, GPVec<CMTransfrag>& mgt,GIntHash<int> &gpos,
		int geneno,int strand,GList<CPrediction>& pred,GList<CReadAln>& readlist,GPVec<GffObj>& guides) {

	mgt.Sort(mgtrabundCmp); // sort transfrags from the most abundant to the least, with guides coming in first
	GVec<CMPrediction> localpred;
	bool first=true;

	// determine "fuzzy" starts and ends
	int ng=0; // stores the number of guides
	for(int t=0;t<mgt.Count();t++) {
		mgt[t]->nf=0;
		mgt[t]->nl=mgt[t]->transfrag->nodes.Count()-1;
		if(mgt[t]->transfrag->real) { // this is a guide => it has no fuzzy endings
			ng++;
			GVec<int> alltr(1,t); // all transfrags used in creating path
			int g=-1;
			for(int j=0;j<mgt[t]->read.Count();j++) {
				int r=mgt[t]->read[j];
				if(readlist[r]->tinfo->g > -1) {
					g=readlist[r]->tinfo->g;
					break;
				}
			}
			CPrediction *p=store_merge_prediction(mgt[t]->transfrag->abundance,alltr,mgt,mgt[t]->transfrag->nodes,no2gnode,strand,geneno,first,readlist,guides,g);
			mgt[t]->transfrag->abundance=0; // by making it 0 we only consider extentions worth taking
			CMPrediction mp(p,mgt[t]->transfrag->nodes,mgt[t]->transfrag->pattern,mgt[t]->transfrag->pattern);
			for(int i=0;i<mgt[t]->transfrag->nodes.Count();i++) {
				mp.b[mgt[t]->transfrag->nodes[i]]=0;
				if(i && no2gnode[mgt[t]->transfrag->nodes[i]]->start==no2gnode[mgt[t]->transfrag->nodes[i-1]]->end+1) {
					int *pos=gpos[edge(mgt[t]->transfrag->nodes[i-1],mgt[t]->transfrag->nodes[i],gno)];
					if(pos) mp.b[*pos]=0;
				}
			}
			localpred.Add(mp);
		}
		/*else {
			int i=0;
			while(i<mgt[t]->transfrag->nodes.Count()) {
				if(no2gnode[mgt[t]->transfrag->nodes[i]]->parent.Count()>2) { mgt[t]->nf=i; break; }
				if(no2gnode[mgt[t]->transfrag->nodes[i]]->parent.Count()>1) { // two parents
					// see if I can still continue
					mgt[t]->nf=i; // I have to stop because both parents are not source
					if(no2gnode[mgt[t]->transfrag->nodes[i]]->parent[0] && no2gnode[mgt[t]->transfrag->nodes[i]]->parent[1]) {
						break;
					} // else: if it's source, just ignore the source and go to the other checks below
				}
				if(no2gnode[mgt[t]->transfrag->nodes[i]]->child.Count()>2) break;
				if(no2gnode[mgt[t]->transfrag->nodes[i]]->child.Count()>1) { // two children -> see if I can ignore one
					if(no2gnode[mgt[t]->transfrag->nodes[i]]->child[0]<gno-1 && no2gnode[mgt[t]->transfrag->nodes[i]]->child[1]<gno-1) break;
				}
				if(i<mgt[t]->transfrag->nodes.Count()-1 && no2gnode[mgt[t]->transfrag->nodes[i]]->end+1 != no2gnode[mgt[t]->transfrag->nodes[i+1]]->start) break;
				i++;
			}


			i=mgt[t]->transfrag->nodes.Count()-1;
			while(i>=0) {
				if(no2gnode[mgt[t]->transfrag->nodes[i]]->child.Count()>2) {mgt[t]->nl=i; break;}
				if(no2gnode[mgt[t]->transfrag->nodes[i]]->child.Count()>1) { // two children
					// see if I can still continue
					mgt[t]->nl=i; // I have to stop because both parents are not sink
					if(no2gnode[mgt[t]->transfrag->nodes[i]]->child[0]<gno-1 && no2gnode[mgt[t]->transfrag->nodes[i]]->child[1]<gno-1) {
						break;
					} // else: if on is sink, just ignore it and go to the other checks below
				}
				if(no2gnode[mgt[t]->transfrag->nodes[i]]->parent.Count()>2) break;
				if(no2gnode[mgt[t]->transfrag->nodes[i]]->parent.Count()>1) { // two parents -> see if I can ignore one
					if(no2gnode[mgt[t]->transfrag->nodes[i]]->parent[0] && no2gnode[mgt[t]->transfrag->nodes[i]]->parent[1]) break;
				}
				if(i && no2gnode[mgt[t]->transfrag->nodes[i]]->start != no2gnode[mgt[t]->transfrag->nodes[i-1]]->end+1) break;
				i--;
			}

			// erase fuzzy pattern
			for(i=0;i<mgt[t]->nf;i++) {
				mgt[t]->transfrag->pattern[mgt[t]->transfrag->nodes[i]]=0;
				int *pos=NULL;
				if(i<mgt[t]->nf-1) {
					pos=gpos[edge(mgt[t]->transfrag->nodes[i],mgt[t]->transfrag->nodes[i+1],gno)];
				}
				else {
					pos=gpos[edge(mgt[t]->transfrag->nodes[i],mgt[t]->transfrag->nodes[mgt[t]->nf],gno)];
				}
				if(pos) mgt[t]->transfrag->pattern[*pos]=0;
			}
			for(i=mgt[t]->nl+1;i<mgt[t]->transfrag->nodes.Count();i++) {
				mgt[t]->transfrag->pattern[mgt[t]->transfrag->nodes[i]]=0;
				int *pos=NULL;
				if(i>mgt[t]->nl+1) {
					pos=gpos[edge(mgt[t]->transfrag->nodes[i-1],mgt[t]->transfrag->nodes[i],gno)];
				}
				else {
					pos=gpos[edge(mgt[t]->transfrag->nodes[mgt[t]->nl],mgt[t]->transfrag->nodes[i],gno)];
				}
				if(pos) mgt[t]->transfrag->pattern[*pos]=0;
			}
		}*/
		// add transfrag to node
		//fprintf(stderr,"Add trf %d with abund=%f to nodes nf=%d(%d) and nl=%d(%d)\n",t,mgt[t]->transfrag->abundance, mgt[t]->transfrag->nodes[mgt[t]->nf],mgt[t]->nf,mgt[t]->transfrag->nodes[mgt[t]->nl],mgt[t]->nl);
		no2gnode[mgt[t]->transfrag->nodes[mgt[t]->nf]]->trf.Add(t);
		no2gnode[mgt[t]->transfrag->nodes[mgt[t]->nl]]->trf.Add(t);
	}



	// get compatibilities
	for(int t1=0;t1<mgt.Count();t1++) if(mgt[t1]->transfrag->nodes[0] && mgt[t1]->transfrag->nodes.Last()<gno-1 && (mgt[t1]->transfrag->abundance || mgt[t1]->transfrag->real)) { // transfrags are processed in increasing order -> important for the later considerations

		/*
		{ // DEBUG ONLY
			fprintf(stderr,"Start transcript %d(%d-%d) cov=%g ",t1,mgt[t1]->nf,mgt[t1]->nl,mgt[t1]->transfrag->abundance);
			//if(mgt[t1]->transfrag->real) fprintf(stderr,"1 out of %d :",mgt.Count());
			//else fprintf(stderr,"0 out of %d:",mgt.Count());
			//for(int i=0;i<mgt[t1]->transfrag->nodes.Count();i++) fprintf(stderr," %d",mgt[t1]->transfrag->nodes[i]);
			fprintf(stderr,"\n");
		}
		*/

		GVec<int> alltr(1,t1);
		GBitVec pathpat(mgt[t1]->transfrag->pattern);
		float cov=mgt[t1]->transfrag->abundance*mgt[t1]->len;
		mgt[t1]->transfrag->abundance=0;

		GBitVec tforlater(mgt.Count());
		int ntforlater=0;

		int min=mgt[t1]->transfrag->nodes[mgt[t1]->nf];
		int max=mgt[t1]->transfrag->nodes[mgt[t1]->nl];

		int start=ng;
		if(t1+1>ng) start=t1+1;

		for(int t2=start;t2<mgt.Count();t2++) if(mgt[t2]->transfrag->nodes[0] && mgt[t2]->transfrag->nodes.Last()<gno-1 && mgt[t2]->transfrag->abundance) { // only if the transfrag wasn't used already
			if(mgt[t2]->transfrag->nodes[mgt[t2]->nl]<min || mgt[t2]->transfrag->nodes[mgt[t2]->nf]>max) { tforlater[t2]=1; ntforlater++;}
			else if((mgt[t2]->transfrag->pattern & pathpat) == mgt[t2]->transfrag->pattern) { // t2 is included in path seen so far
				//fprintf(stderr,"t2=%d included\n",t2);
				cov+=mgt[t2]->transfrag->abundance*mgt[t2]->len;
				mgt[t2]->transfrag->abundance=0;
				alltr.Add(t2);
			}
			else {
				add_transfrag_to_path(t2,tforlater,ntforlater,mgt,pathpat,min,max,no2gnode,gno,gpos,alltr,cov);
			}
		} // end for(int t2=t1+1;t2<transfrag.Count();t2++)

		/*
		{ // DEBUG ONLY
			fprintf(stderr,"min=%d max=%d pathpat=",min,max);
			printBitVec(pathpat);
			fprintf(stderr,"\n");
		}
		*/


		if(!mgt[t1]->transfrag->real ||
				(min<mgt[t1]->transfrag->nodes[mgt[t1]->nf]) || (max>mgt[t1]->transfrag->nodes[mgt[t1]->nl])) { // if transfrag is not real or is extending beyond the initial transfrag
			GVec<int> printpath;
			if(alltr.Count()==1) printpath=mgt[t1]->transfrag->nodes; // only initial transfrag is in path
			else {
				for(int i=0;i<alltr.Count();i++) {
					int t=alltr[i];
					if(mgt[t]->transfrag->nodes[0]<min) { // only fuzzy continuations here
						int j=0;
						while(mgt[t]->transfrag->nodes[j]<min) { pathpat[j]=1;j++;}
						min=mgt[t]->transfrag->nodes[0];
					}
					if(mgt[t]->transfrag->nodes.Last()>max) { // only fuzzy continuations here
						int j=mgt[t]->transfrag->nodes.Count()-1;
						while(mgt[t]->transfrag->nodes[j]>max) { pathpat[j]=1;j--;}
						max=mgt[t]->transfrag->nodes.Last();
					}
				}

				for(int i=1;i<gno-1;i++) {
					if(pathpat[i]) {
						printpath.Add(i);
					}
				}

			}

			CPrediction *p=store_merge_prediction(cov,alltr,mgt,printpath,no2gnode,strand,geneno,first,readlist,guides,-1);
			CMPrediction mp(p,printpath,pathpat,pathpat);
			for(int i=0;i<printpath.Count();i++) {
				mp.b[printpath[i]]=0;
				if(i && no2gnode[printpath[i]]->start==no2gnode[printpath[i-1]]->end+1) {
					int *pos=gpos[edge(printpath[i-1],printpath[i],gno)];
					if(pos) mp.b[*pos]=0;
				}
			}
			localpred.Add(mp);
		}
		/*// this will increase the guide's coverage and I am not sure it is fair since I am then eliminating based on isofrac fraction
		else if(mgt[t1]->transfrag->real && cov) { // the transfrag is real but didn't extend past its limits
			// I already added the transfrag to the localpred
			localpred[t1].p->cov+=cov;
		}
		*/

	} // end for(int t1=0;t1<transfrag.Count();t1++)

	// process predictions
	localpred.Sort(mpredCmp);
	for(int i=localpred.Count()-1;i>=0;i--) {

		/*
		{ // DEBUG ONLY
			fprintf(stderr,"Get ");
			if(localpred[i].p->t_eq) fprintf(stderr,"guide ");
			fprintf(stderr,"prediction %d %d-%d cov=%g ",i,localpred[i].p->start,localpred[i].p->end,localpred[i].p->cov);
			//printBitVec(localpred[i].b);
			fprintf(stderr,"\n");
		}
		*/

		if(localpred[i].p->t_eq==NULL) { // if a guide -> just store it
			for(int j=0;j<i;j++) if(localpred[j].p!=NULL && localpred[j].p->overlap(localpred[i].p)) {

				bool remove=false;

				if(localpred[j].p->exons.Count()==1) remove=true;
				else {
					if(localpred[i].p->exons.Count()==1) {
						if(overlaps_one_exon_only(localpred[j],localpred[i])) remove=true;
						else if(!retained_intron && localpred[i].p->cov<localpred[j].p->cov) remove=true;
					}
					else if((localpred[i].b & localpred[j].b) == localpred[i].b) {
						if(has_retained_intron(localpred[j],localpred[i],gno,gpos)) {
							if(!retained_intron && localpred[i].p->cov<localpred[j].p->cov) remove=true;
						}
						else remove=true;
					}
				}

				/*
				if( (localpred[i].p->exons.Count()>1 &&
						//(localpred[i].p->exons.Count() <= localpred[j].p->exons.Count()) &&  // this should be true from the sorting procedure
						((localpred[i].b & localpred[j].b) == localpred[i].b)) ||  // the way this works gets read of the retained intron transcripts whih could be very long -> maybe I should revise it
						//(localpred[i].p->exons.Count()==1 && localpred[j].p->exons.Count()==1 && // if localpred[j] has one exon so does locapred[i] from the way we sorted
						(localpred[j].p->exons.Count()==1 &&
								localpred[j].p->overlap(localpred[i].p)) ) {

					if(localpred[j].p->exons.Count()==1 || retained_intron
							|| has_retained_intron(localpred[j],localpred[i],gno,gpos)
							|| localpred[i].p->cov<localpred[j].p->cov) {
			    */

				if(remove) {
					//fprintf(stderr,"Pred %d includes pred %d\n",j,i);
					localpred[j].p->cov*=localpred[j].p->tlen; // I need to adjust coverage later
					if(localpred[j].p->t_eq==NULL) {
						if(localpred[i].p->start<localpred[j].p->start) {
							localpred[j].p->tlen+=localpred[j].p->start-localpred[i].p->start;
							localpred[j].p->start=localpred[i].p->start;
							localpred[j].p->exons[0].start=localpred[i].p->exons[0].start;
						}
						if(localpred[i].p->end>localpred[j].p->end) {
							localpred[j].p->tlen+=localpred[i].p->end-localpred[j].p->end;
							localpred[j].p->end=localpred[i].p->end;
							localpred[j].p->exons.Last().end=localpred[i].p->exons.Last().end;
						}
					}

					// need to update weighted coverage here
					int ovplen=compute_ovlp(localpred[i].p,localpred[j].p);
					localpred[j].p->cov=(ovplen*localpred[i].p->cov+localpred[j].p->cov)/localpred[j].p->tlen;
					localpred[i].p->cov=0;

					break;
				}
			}
		}


		if(localpred[i].p->t_eq || localpred[i].p->cov) pred.Add(localpred[i].p);
		else { // I have a prediction that needs to be cleaned up
			//fprintf(stderr,"delete prediction\n");
			delete localpred[i].p;
			localpred[i].p=NULL;
		}
	}

	return(geneno);
}

/*
// this version just continues a transcript with the most abundant path
int merge_transfrags(int gno,GPVec<CGraphnode>& no2gnode, GPVec<CMTransfrag>& mgt,GIntHash<int> &gpos,
		int geneno,int strand,GList<CPrediction>& pred,GList<CReadAln>& readlist,GPVec<GffObj>& guides) {

	mgt.Sort(mgtrabundCmp); // sort transfrags from the most abundant to the least, with guides coming in first
	GVec<CMPrediction> localpred;

	// determine "fuzzy" starts and ends
	int ng=0; // stores the number of guides
	for(int t=0;t<mgt.Count();t++) {
		mgt[t]->nf=0;
		mgt[t]->nl=mgt[t]->transfrag->nodes.Count()-1;
		if(mgt[t]->transfrag->real) { // this is a guide => it has no fuzzy endings
			ng++;
		}
		else {

			// see if previous nodes are all adjacent and reach the source --> why do I care about this??
			// *  bool cont=true;
			int n=mgt[t]->transfrag->nodes[0];

			while(n) {
				if(no2gnode[n]->parent.Count()>2) { cont=false; break;} // if more than one parent -> I can not continue
				if(no2gnode[n]->parent.Count()>1) { // two parents -> check if one is source
					if(no2gnode[n]->parent[0] && no2gnode[n]->parent[1]) { cont=false;break;} // both parents are not source
					if(no2gnode[n]->parent[0] && no2gnode[no2gnode[n]->parent[0]]->end+1 != no2gnode[n]->start) {cont=false;break;} // the parent that continues is not adjacent
					if(no2gnode[n]->parent[1] && no2gnode[no2gnode[n]->parent[1]]->end+1 != no2gnode[n]->start) {cont=false;break;}
				}
				else { // only one parent here
					if(!no2gnode[n]->parent[0]) break; // if it's source then I am done
					if(no2gnode[no2gnode[n]->parent[0]]->end+1 != no2gnode[n]->start) {cont=false;break;} // it's not adjacent
				}
				n--;
			}

			int i=0;
			if(cont) while(i<mgt[t]->transfrag->nodes.Count()) {*   //
			int i=0;
			while(i<mgt[t]->transfrag->nodes.Count()) {
				if(no2gnode[mgt[t]->transfrag->nodes[i]]->parent.Count()>2) { mgt[t]->nf=i; break; }
				if(no2gnode[mgt[t]->transfrag->nodes[i]]->parent.Count()>1) { // two parents
					// see if I can still continue
					mgt[t]->nf=i; // I have to stop because both parents are not source
					if(no2gnode[mgt[t]->transfrag->nodes[i]]->parent[0] && no2gnode[mgt[t]->transfrag->nodes[i]]->parent[1]) {
						break;
					} // else: if it's source, just ignore the source and go to the other checks below
				}
				if(no2gnode[mgt[t]->transfrag->nodes[i]]->child.Count()>2) break;
				if(no2gnode[mgt[t]->transfrag->nodes[i]]->child.Count()>1) { // two children -> see if I can ignore one
					if(no2gnode[mgt[t]->transfrag->nodes[i]]->child[0]<gno-1 && no2gnode[mgt[t]->transfrag->nodes[i]]->child[1]<gno-1) break;
				}
				if(i<mgt[t]->transfrag->nodes.Count()-1 && no2gnode[mgt[t]->transfrag->nodes[i]]->end+1 != no2gnode[mgt[t]->transfrag->nodes[i+1]]->start) break;
				i++;
			}


			// see if next nodes are all adjacent and reach the sink -- why??
			//  * cont=true;
			n=mgt[t]->transfrag->nodes.Last();
			while(n<gno-1) {
				if(no2gnode[n]->child.Count()>2) { cont=false; break;} // if more than one child -> I can not continue
				if(no2gnode[n]->child.Count()>1) { // two children -> check if one is sink
					if(no2gnode[n]->child[0]<gno-1 && no2gnode[n]->child[1]<gno-1) { cont=false;break;} // both children are not sink
					if(no2gnode[n]->child[0]<gno-1 && no2gnode[no2gnode[n]->child[0]]->start != no2gnode[n]->end) {cont=false;break;} // the parent that continues is not adjacent
					if(no2gnode[n]->child[1]<gno-1 && no2gnode[no2gnode[n]->child[1]]->start != no2gnode[n]->end) {cont=false;break;} // the parent that continues is not adjacent
				}
				else { // only one child here
					if(no2gnode[n]->child[0]==gno-1) break; // if it's sink then I am done
					if(no2gnode[no2gnode[n]->child[0]]->start != no2gnode[n]->end+1) {cont=false;break;} // it's not adjacent
				}
				n++;
			}

			i=mgt[t]->transfrag->nodes.Count()-1;
			if(cont) while(i>=0) {*   //
			i=mgt[t]->transfrag->nodes.Count()-1;
			while(i>=0) {
				if(no2gnode[mgt[t]->transfrag->nodes[i]]->child.Count()>2) {mgt[t]->nl=i; break;}
				if(no2gnode[mgt[t]->transfrag->nodes[i]]->child.Count()>1) { // two children
					// see if I can still continue
					mgt[t]->nl=i; // I have to stop because both parents are not sink
					if(no2gnode[mgt[t]->transfrag->nodes[i]]->child[0]<gno-1 && no2gnode[mgt[t]->transfrag->nodes[i]]->child[1]<gno-1) {
						break;
					} // else: if on is sink, just ignore it and go to the other checks below
				}
				if(no2gnode[mgt[t]->transfrag->nodes[i]]->parent.Count()>2) break;
				if(no2gnode[mgt[t]->transfrag->nodes[i]]->parent.Count()>1) { // two parents -> see if I can ignore one
					if(no2gnode[mgt[t]->transfrag->nodes[i]]->parent[0] && no2gnode[mgt[t]->transfrag->nodes[i]]->parent[1]) break;
				}
				if(i && no2gnode[mgt[t]->transfrag->nodes[i]]->start != no2gnode[mgt[t]->transfrag->nodes[i-1]]->end+1) break;
				i--;
			}

			// erase fuzzy pattern
			for(i=0;i<mgt[t]->nf;i++) {
				mgt[t]->transfrag->pattern[mgt[t]->transfrag->nodes[i]]=0;
				int *pos=NULL;
				if(i<mgt[t]->nf-1) {
					pos=gpos[edge(mgt[t]->transfrag->nodes[i],mgt[t]->transfrag->nodes[i+1],gno)];
				}
				else {
					pos=gpos[edge(mgt[t]->transfrag->nodes[i],mgt[t]->transfrag->nodes[mgt[t]->nf],gno)];
				}
				if(pos) mgt[t]->transfrag->pattern[*pos]=0;
			}
			for(i=mgt[t]->nl+1;i<mgt[t]->transfrag->nodes.Count();i++) {
				mgt[t]->transfrag->pattern[mgt[t]->transfrag->nodes[i]]=0;
				int *pos=NULL;
				if(i>mgt[t]->nl+1) {
					pos=gpos[edge(mgt[t]->transfrag->nodes[i-1],mgt[t]->transfrag->nodes[i],gno)];
				}
				else {
					pos=gpos[edge(mgt[t]->transfrag->nodes[mgt[t]->nl],mgt[t]->transfrag->nodes[i],gno)];
				}
				if(pos) mgt[t]->transfrag->pattern[*pos]=0;
			}

		}
	}

	bool first=true;

	// get compatibilities
	for(int t1=0;t1<mgt.Count();t1++) if(mgt[t1]->transfrag->nodes[0] && mgt[t1]->transfrag->nodes.Last()<gno-1 && (mgt[t1]->transfrag->abundance || mgt[t1]->transfrag->real)) { // transfrags are processed in increasing order -> important for the later considerations

		int n1=mgt[t1]->nl+1;


		{ // DEBUG ONLY
			fprintf(stderr,"Start transcript %d(%d-%d) ",t1,mgt[t1]->nf,mgt[t1]->nl);
			if(mgt[t1]->transfrag->real) fprintf(stderr,"1 out of %d :",mgt.Count());
			else fprintf(stderr,"0 out of %d:",mgt.Count());
			for(int i=0;i<n1;i++) fprintf(stderr," %d",mgt[t1]->transfrag->nodes[i]);
			fprintf(stderr,"\n");
		}


		GVec<int> uniont(1,t1);
		GBitVec unionpat(mgt[t1]->transfrag->pattern);

		// set compatibilities with transfrags that are not guides
		bool comp;
		int max=ng;
		if(t1+1>ng) max=t1+1;

		for(int t2=max;t2<mgt.Count();t2++) if(mgt[t2]->transfrag->nodes[0] && mgt[t2]->transfrag->nodes.Last()<gno-1 && mgt[t2]->transfrag->abundance) { // only if the transfrag wasn't used already
			if((mgt[t2]->transfrag->pattern & mgt[t1]->transfrag->pattern) == mgt[t2]->transfrag->pattern) { // t2 is included in t1
				mgt[t1]->transfrag->abundance+=mgt[t2]->transfrag->abundance;
				mgt[t2]->transfrag->abundance=0;
			}
			else {
				// here check compatibility between t1 and t2;
				int n2=mgt[t2]->nl+1;
				int i1=mgt[t1]->nf;
				int i2=mgt[t2]->nf;
				comp=true;
				while(i1<n1 && i2<n2) {
					if(mgt[t1]->transfrag->nodes[i1]==mgt[t2]->transfrag->nodes[i2]) {
						if(i1==n1-1 || i2==n2-1) { // one transcript finishes -> no need to check anymore
							i1=n1;i2=n2;
						}
						else { // advance the smallest one
							if(mgt[t1]->transfrag->nodes[i1+1]<mgt[t2]->transfrag->nodes[i2+1]) i1++;
							else i2++;
						}
					}
					else if(mgt[t1]->transfrag->nodes[i1]<mgt[t2]->transfrag->nodes[i2]) {
						i1++;
						if(conflict(i1,mgt[t2]->transfrag->nodes[i2],mgt[t1]->transfrag->nodes,n1,no2gnode,mgt[t1]->transfrag->pattern,gno,gpos)) {
							comp=false;
							break;
						}
					}
					else {
						i2++;
						if(conflict(i2,mgt[t1]->transfrag->nodes[i1],mgt[t2]->transfrag->nodes,n2,no2gnode,mgt[t2]->transfrag->pattern,gno,gpos)) {
							comp=false;
							break;
						}
					}
				}
				if(comp) {
					uniont.Add(t2);
					unionpat = unionpat | mgt[t2]->transfrag->pattern;
					//fprintf(stderr," add transfrag %d\n",t2);
				}
			}
		} // end for(int t2=t1+1;t2<transfrag.Count();t2++)

		geneno=add_transfrag_to_prediction(uniont,unionpat,gno,no2gnode,mgt,gpos,geneno,strand,localpred,readlist,guides,first);

	} // end for(int t1=0;t1<transfrag.Count();t1++)

	// process predictions
	localpred.Sort(mpredCmp);
	for(int i=localpred.Count()-1;i>=0;i--) if(localpred[i].p->exons.Count()>1) { // only for two exon predictions


		fprintf(stderr,"Get ");
		if(localpred[i].p->t_eq) fprintf(stderr,"guide ");
		fprintf(stderr,"prediction %d-%d ",localpred[i].p->start,localpred[i].p->end);
		printBitVec(localpred[i].b);
		fprintf(stderr,"\n");


		if(localpred[i].p->t_eq==NULL) { // if a guide -> just store it
			for(int j=0;j<i;j++) {
				if((localpred[i].b & localpred[j].b) == localpred[i].b) {
					localpred[j].p->cov+=localpred[i].p->cov;
					localpred[i].p->cov=0;
					if(localpred[j].p->t_eq==NULL) { // see if I need to adjust start/end
						if(localpred[i].p->start<localpred[j].p->start) {
							localpred[j].p->start=localpred[i].p->start;
							localpred[j].p->exons[0].start=localpred[i].p->exons[0].start;
						}
						if(localpred[i].p->end>localpred[j].p->end) {
							localpred[j].p->end=localpred[i].p->end;
							localpred[j].p->exons.Last().end=localpred[i].p->exons.Last().end;
						}
					}
				}
			}

		}

		if(localpred[i].p->t_eq || localpred[i].p->cov) pred.Add(localpred[i].p);
	}

	return(geneno);
}
*/

/*   // needs more work
// this version uses the maximum flow algorithm to merge the transcripts: i.e. if 2 transcripts continue another one, the
// abundance will be share proportionally among them: it might introduce too many errors (false positives) while dropping others
int merge_transfrags_EM(int gno,GPVec<CGraphnode>& no2gnode, GPVec<CMTransfrag>& mgt,GIntHash<int> &gpos,
		int geneno,int strand,GList<CPrediction>& pred,GList<CReadAln>& readlist,GPVec<GffObj>& guides) {

	mgt.Sort(mgtrabundCmp); // sort transfrags from the most abundant to the least, with guides coming in first

	bool first=true;

	// determine "fuzzy" starts and ends
	int ng=0; // stores the number of guides
	for(int t=0;t<mgt.Count();t++) {
		mgt[t]->nf=0;
		mgt[t]->nl=mgt[t]->transfrag->nodes.Count()-1;
		if(mgt[t]->transfrag->real) { // this is a guide => it has no fuzzy endings
			ng++;
			GVec<int> alltr(1,t); // all transfrags used in creating path
			geneno=store_merge_prediction(pred,alltr,mgt,mgt[t]->transfrag->nodes,no2gnode,strand,geneno,first,readlist,guides);
		}
		else {

			// see if previous nodes are all adjacent and reach the source
			bool cont=true;
			int n=mgt[t]->transfrag->nodes[0];
			while(n) {
				if(no2gnode[n]->parent.Count()>2) { cont=false; break;} // if more than one parent -> I can not continue
				if(no2gnode[n]->parent.Count()>1) { // two parents -> check if one is source
					if(no2gnode[n]->parent[0] && no2gnode[n]->parent[1]) { cont=false;break;} // both parents are not source
					if(no2gnode[n]->parent[0] && no2gnode[no2gnode[n]->parent[0]]->end+1 != no2gnode[n]->start) {cont=false;break;} // the parent that continues is not adjacent
					if(no2gnode[n]->parent[1] && no2gnode[no2gnode[n]->parent[1]]->end+1 != no2gnode[n]->start) {cont=false;break;}
				}
				else { // only one parent here
					if(!no2gnode[n]->parent[0]) break; // if it's source then I am done
					if(no2gnode[no2gnode[n]->parent[0]]->end+1 != no2gnode[n]->start) {cont=false;break;} // it's not adjacent
				}
				n--;
			}

			int i=0;
			if(cont) while(i<mgt[t]->transfrag->nodes.Count()) {
				if(no2gnode[mgt[t]->transfrag->nodes[i]]->parent.Count()>2) { mgt[t]->nf=i; break; }
				if(no2gnode[mgt[t]->transfrag->nodes[i]]->parent.Count()>1) { // two parents
					// see if I can still continue
					mgt[t]->nf=i; // I have to stop because both parents are not source
					if(no2gnode[mgt[t]->transfrag->nodes[i]]->parent[0] && no2gnode[mgt[t]->transfrag->nodes[i]]->parent[1]) {
						break;
					} // else: if it's source, just ignore the source and go to the other checks below
				}
				if(no2gnode[mgt[t]->transfrag->nodes[i]]->child.Count()>2) break;
				if(no2gnode[mgt[t]->transfrag->nodes[i]]->child.Count()>1) { // two children -> see if I can ignore one
					if(no2gnode[mgt[t]->transfrag->nodes[i]]->child[0]<gno-1 && no2gnode[mgt[t]->transfrag->nodes[i]]->child[1]<gno-1) break;
				}
				if(i<mgt[t]->transfrag->nodes.Count()-1 && no2gnode[mgt[t]->transfrag->nodes[i]]->end+1 != no2gnode[mgt[t]->transfrag->nodes[i+1]]->start) break;
				i++;
			}


			// see if next nodes are all adjacent and reach the sink
			cont=true;
			n=mgt[t]->transfrag->nodes.Last();
			while(n<gno-1) {
				if(no2gnode[n]->child.Count()>2) { cont=false; break;} // if more than one child -> I can not continue
				if(no2gnode[n]->child.Count()>1) { // two children -> check if one is sink
					if(no2gnode[n]->child[0]<gno-1 && no2gnode[n]->child[1]<gno-1) { cont=false;break;} // both children are not sink
					if(no2gnode[n]->child[0]<gno-1 && no2gnode[no2gnode[n]->child[0]]->start != no2gnode[n]->end) {cont=false;break;} // the parent that continues is not adjacent
					if(no2gnode[n]->child[1]<gno-1 && no2gnode[no2gnode[n]->child[1]]->start != no2gnode[n]->end) {cont=false;break;} // the parent that continues is not adjacent
				}
				else { // only one child here
					if(no2gnode[n]->child[0]==gno-1) break; // if it's sink then I am done
					if(no2gnode[no2gnode[n]->child[0]]->start != no2gnode[n]->end+1) {cont=false;break;} // it's not adjacent
				}
				n++;
			}

			i=mgt[t]->transfrag->nodes.Count()-1;
			if(cont) while(i>=0) {
				if(no2gnode[mgt[t]->transfrag->nodes[i]]->child.Count()>2) {mgt[t]->nl=i; break;}
				if(no2gnode[mgt[t]->transfrag->nodes[i]]->child.Count()>1) { // two children
					// see if I can still continue
					mgt[t]->nl=i; // I have to stop because both parents are not sink
					if(no2gnode[mgt[t]->transfrag->nodes[i]]->child[0]<gno-1 && no2gnode[mgt[t]->transfrag->nodes[i]]->child[1]<gno-1) {
						break;
					} // else: if on is sink, just ignore it and go to the other checks below
				}
				if(no2gnode[mgt[t]->transfrag->nodes[i]]->parent.Count()>2) break;
				if(no2gnode[mgt[t]->transfrag->nodes[i]]->parent.Count()>1) { // two parents -> see if I can ignore one
					if(no2gnode[mgt[t]->transfrag->nodes[i]]->parent[0] && no2gnode[mgt[t]->transfrag->nodes[i]]->parent[1]) break;
				}
				if(i && no2gnode[mgt[t]->transfrag->nodes[i]]->start != no2gnode[mgt[t]->transfrag->nodes[i-1]]->end+1) break;
				i--;
			}

			// erase fuzzy pattern
			for(i=0;i<mgt[t]->nf;i++) {
				mgt[t]->transfrag->pattern[mgt[t]->transfrag->nodes[i]]=0;
				int *pos=NULL;
				if(i<mgt[t]->nf-1) {
					pos=gpos[edge(mgt[t]->transfrag->nodes[i],mgt[t]->transfrag->nodes[i+1],gno)];
				}
				else {
					pos=gpos[edge(mgt[t]->transfrag->nodes[i],mgt[t]->transfrag->nodes[mgt[t]->nf],gno)];
				}
				if(pos) mgt[t]->transfrag->pattern[*pos]=0;
			}
			for(i=mgt[t]->nl+1;i<mgt[t]->transfrag->nodes.Count();i++) {
				mgt[t]->transfrag->pattern[mgt[t]->transfrag->nodes[i]]=0;
				int *pos=NULL;
				if(i>mgt[t]->nl+1) {
					pos=gpos[edge(mgt[t]->transfrag->nodes[i-1],mgt[t]->transfrag->nodes[i],gno)];
				}
				else {
					pos=gpos[edge(mgt[t]->transfrag->nodes[mgt[t]->nl],mgt[t]->transfrag->nodes[i],gno)];
				}
				if(pos) mgt[t]->transfrag->pattern[*pos]=0;
			}

		}

		// add transfrag to node
		//fprintf(stderr,"Add trf %d with abund=%f to nodes nf=%d(%d) and nl=%d(%d)\n",t,mgt[t]->transfrag->abundance, mgt[t]->transfrag->nodes[mgt[t]->nf],mgt[t]->nf,mgt[t]->transfrag->nodes[mgt[t]->nl],mgt[t]->nl);
		no2gnode[mgt[t]->transfrag->nodes[mgt[t]->nf]]->trf.Add(t);
		no2gnode[mgt[t]->transfrag->nodes[mgt[t]->nl]]->trf.Add(t);

	}

	// compute node capacities
	for(int i=1;i<gno-1;i++) {
		CGraphnode *inode=no2gnode[i]; // this is here only because of the DEBUG option below
		int nn=inode->trf.Count();
		float abundin=0;
		float abundout=0;
		for(int j=0;j<nn;j++){
			int t=inode->trf[j];
			if(mgt[t]->transfrag->nodes[mgt[t]->nl]==i) { // transfrag ends at this node (in transfrag)
				abundin+=mgt[t]->transfrag->abundance;
			}
			else if(mgt[t]->transfrag->nodes[mgt[t]->nf]==i) { // transfrag starts at this node (out transfrag)
				abundout+=mgt[t]->transfrag->abundance;
			}
		}

		if(abundin) inode->rate=abundout/abundin;
		if(abundout) inode->capacity=abundout; // node capacity tells me how much of that node coverage I can use give how many transfrags leave the node
		else inode->capacity=abundin;
	} // end for i

	// get compatibilities
	for(int t1=0;t1<mgt.Count();t1++) if(mgt[t1]->transfrag->nodes[0] && mgt[t1]->transfrag->nodes.Last()<gno-1 && (mgt[t1]->transfrag->abundance || mgt[t1]->transfrag->real)) { // transfrags are processed in increasing order -> important for the later considerations

		int n1=mgt[t1]->nl+1;


		{ // DEBUG ONLY
			fprintf(stderr,"Start transcript %d ",t1);
			if(mgt[t1]->transfrag->real) fprintf(stderr,"1 out of %d :",mgt.Count());
			else fprintf(stderr,"0 out of %d:",mgt.Count());
			for(int i=0;i<n1;i++) fprintf(stderr," %d",mgt[t1]->transfrag->nodes[i]);
			fprintf(stderr,"\n");
		}


		GVec<int> uniont(1,t1);
		GBitVec unionpat(mgt[t1]->transfrag->pattern);

		// set compatibilities with transfrags that are not guides
		bool comp;
		int max=t1+1;
		//if(ng>max) max=ng;

		for(int t2=max;t2<mgt.Count();t2++) if(mgt[t2]->transfrag->nodes[0] && mgt[t2]->transfrag->nodes.Last()<gno-1 && mgt[t2]->transfrag->abundance) { // only if the transfrag wasn't used already
			// here check compatibility between t1 and t2;
			int n2=mgt[t2]->nl+1;
			int i1=mgt[t1]->nf;
			int i2=mgt[t2]->nf;
			comp=true;
			while(i1<n1 && i2<n2) {
				if(mgt[t1]->transfrag->nodes[i1]==mgt[t2]->transfrag->nodes[i2]) {
					if(i1==n1-1 || i2==n2-1) { // one transcript finishes -> no need to check anymore
						i1=n1;i2=n2;
					}
					else { // advance the smallest one
						if(mgt[t1]->transfrag->nodes[i1+1]<mgt[t2]->transfrag->nodes[i2+1]) i1++;
						else i2++;
					}
				}
				else if(mgt[t1]->transfrag->nodes[i1]<mgt[t2]->transfrag->nodes[i2]) {
					i1++;
					if(conflict(i1,mgt[t2]->transfrag->nodes[i2],mgt[t1]->transfrag->nodes,n1,no2gnode,mgt[t1]->transfrag->pattern,gno,gpos)) {
						comp=false;
						break;
					}
				}
				else {
					i2++;
					if(conflict(i2,mgt[t1]->transfrag->nodes[i1],mgt[t2]->transfrag->nodes,n2,no2gnode,mgt[t2]->transfrag->pattern,gno,gpos)) {
						comp=false;
						break;
					}
				}
			}
			if(comp) {
				uniont.Add(t2);
				unionpat = unionpat | mgt[t2]->transfrag->pattern;
				//fprintf(stderr," add transfrag %d\n",t2);
			}

		} // end for(int t2=t1+1;t2<transfrag.Count();t2++)

		geneno=add_transfrag_to_prediction(uniont,unionpat,gno,no2gnode,mgt,gpos,geneno,strand,pred,readlist,guides,first);

	} // end for(int t1=0;t1<transfrag.Count();t1++)

	return(geneno);
}
*/

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

bool onpath(GBitVec& trpattern,GVec<int>& trnode,GBitVec& pathpattern,int mini,int maxi,GPVec<CGraphnode>& no2gnode,int gno,
		GIntHash<int>& gpos) {

	if(trnode[0]<mini) // mini can be reached through transcript
	    if(!no2gnode[mini]->parentpat[trnode[0]])	return false;


	if(trnode.Last()>maxi) // from maxi I can reach end of transcript
	    if(!no2gnode[maxi]->childpat[trnode.Last()]) return false;

	int first=1;

	for(int i=0;i<trnode.Count();i++) {
		if(trnode[i]>=mini && trnode[i]<=maxi) {
	      if(!pathpattern[trnode[i]]) return false;
    	  int *pos=NULL;
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
		//if(sensitivitylevel && inode->child[c]==i+1 && i<gno-2 && inode->end+1==no2gnode[i+1]->start && cnode->len()
		if(inode->child[c]==i+1 && i<gno-2 && inode->end+1==no2gnode[i+1]->start && cnode->len()
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
				else if(transfrag[t]->nodes[0]<=i && transfrag[t]->nodes.Last()>=inode->child[c] &&   // transfrag goes from i to c
						//(transfrag[t]->pattern[inode->child[c]] || transfrag[t]->pattern[i]) &&    // transfrag is not incomplete through these nodes
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
				else if(transfrag[t]->nodes[0]<=i && transfrag[t]->nodes.Last()>=i+1 &&   // transfrag goes from i to c
						//(transfrag[t]->pattern[i+1] || transfrag[t]->pattern[i]) &&    // transfrag is not incomplete through these nodes
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
		//if(sensitivitylevel && inode->parent[p]==i-1 && i>1 && inode->start==no2gnode[i-1]->end+1 && pnode->len() &&
		if(inode->parent[p]==i-1 && i>1 && inode->start==no2gnode[i-1]->end+1 && pnode->len() &&
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
				else if(transfrag[t]->nodes[0]<=inode->parent[p] && transfrag[t]->nodes.Last()>=i &&   // transfrag goes from p to i
						//(transfrag[t]->pattern[inode->parent[p]]||transfrag[t]->pattern[i]) &&    // transfrag is not incomplete through these nodes
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
				else if(transfrag[t]->nodes[0]<=i-1 && transfrag[t]->nodes.Last()>=i &&   // transfrag goes from p to i
						//(transfrag[t]->pattern[i-1]||transfrag[t]->pattern[i]) &&    // transfrag is not incomplete through these nodes
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
	//if(sensitivitylevel && i>1 && inode->start==1+no2gnode[i-1]->end && nodecov[i]*DROP>nodecov[i-1])  { // adjacent to parent
	if(i>1 && inode->start==1+no2gnode[i-1]->end && nodecov[i]*DROP>nodecov[i-1])  { // adjacent to parent
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
		else if(transfrag[t]->nodes[0]!=i) { // this is an in transfrag; I might change this to include the cases where I have a positive in for the transfrag, e.g. if(transfrag[inode->trf[j]]->nodes[0]!=i || transfrag[inode->trf[j]]->in)
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
						if(transfrag[t]->abundance && onpath(transfrag[t]->pattern,transfrag[t]->nodes,pathpat,path.Last(),path[0],no2gnode,gno,gpos)) {

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
	//if(sensitivitylevel && i<gno-2 && inode->end+1==no2gnode[i+1]->start && nodecov[i]*DROP>nodecov[i+1])  { // adjacent to parent
	if(i<gno-2 && inode->end+1==no2gnode[i+1]->start && nodecov[i]*DROP>nodecov[i+1])  { // adjacent to parent
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
		else if(transfrag[t]->nodes.Last()!=i) { // this is an out transfrag; I might change this to include the cases where I have a positive out for the transfrag, e.g. if(transfrag[inode->trf[j]]->nodes.Last()!=i || transfrag[inode->trf[j]]->out)
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

/*
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
*/

void update_capacity(int start,CTransfrag *t,float val,GVec<float>& capacity,GVec<int>& node2path) {
	t->abundance-=val;
	if(t->abundance<epsilon) t->abundance=0;
	for(int j=start;j<t->nodes.Count()-1;j++) { // for all nodes but the last
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


// seems like this would work even if source is not first, and last is not sink
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

	/*
	fprintf(stderr,"Used transcripts:");
	for(int i=0;i<transfrag.Count();i++) if(istranscript[i]) fprintf(stderr," %d(%f)",i,transfrag[i]->abundance);
	fprintf(stderr,"\n");
	*/

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
							//fprintf(stderr,"Update capacity of transfrag[%d] with value %f to %f\n",t,transfrag[t]->abundance, transfrag[t]->abundance-flow[n1][n2]);
							update_capacity(0,transfrag[t],flow[n1][n2],nodecapacity,node2path);
							//if(path[i] && transfrag[t]->nodes.Last()!=gno-1) fragno+=flow[n1][n2];
							flow[n1][n2]=0;
						}
						else {
							if(!i) sumout+=transfrag[t]->abundance;
							flow[n1][n2]-=transfrag[t]->abundance;
							//if(path[i] && transfrag[t]->nodes.Last()!=gno-1) fragno+=transfrag[t]->abundance;
							//fprintf(stderr,"Update capacity of transfrag[%d] with value=%f to 0\n",t,transfrag[t]->abundance);
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

/*
// compute by BFS all coverages of paths that an incomplete transfrag could follow
void bfs_cov1(float &totalcov,int len, float cov, int n, int inode, int t,GPVec<CGraphnode>& no2gnode,
		GPVec<CTransfrag>& transfrag,GVec<float>& nodecov) {

	fprintf(stderr,"bfs[n=%d len=%d cov=%f inode=%d]\n",n,len,cov,inode);

	CGraphnode *node=no2gnode[n];
	int nextnode=transfrag[t]->nodes[inode];

	for(int i=0;i<node->child.Count();i++) {
		if(node->child[i]<nextnode) {
			bfs_cov1(totalcov,len+no2gnode[node->child[i]]->len(),cov+nodecov[node->child[i]],node->child[i],inode,t,no2gnode,transfrag,nodecov);
		}
		else if(node->child[i]==nextnode) { // reached the next fixed node in transfrag
			int newlen=len+no2gnode[nextnode]->len();
			float newcov=cov+nodecov[nextnode];
			if(nextnode==transfrag[t]->nodes.Last()) {
				totalcov+=newcov/newlen;
				fprintf(stderr,"cov[%d]=%f len[%d]=%d newcov=%f newlen=%d\n",nextnode,nodecov[nextnode],nextnode,no2gnode[nextnode]->len(),newcov,newlen);
			}
			else bfs_cov1(totalcov,newlen,newcov,nextnode,inode+1,t,no2gnode,transfrag,nodecov);
			break;
		}
		else break;
	}

}

// compute by BFS all coverages of paths that an incomplete transfrag could follow
void bfs_cov(float &totalcov,int len, float cov, int n, int inode, int t,GPVec<CGraphnode>& no2gnode,
		GPVec<CTransfrag>& transfrag,GVec<float>& nodecov, GIntHash<int> &gpos,int gno) {

	CGraphnode *node=no2gnode[n];
	cov+=nodecov[n];
	len+=node->len();

	if(n==transfrag[t]->nodes.Last()) {
		totalcov+=cov/len;
		return;
	}

	int nextnode=transfrag[t]->nodes[inode];
	if(n==nextnode) { // n is not the last node in transfrag => I can increase inode without stepping over the boundaries
		inode++;
		nextnode=transfrag[t]->nodes[inode];
		int *pos=gpos[edge(n,nextnode,gno)];
		if(!pos || !transfrag[t]->pattern[*pos]) { // there is no edge between n and nextnode -> I might still reach nextnode
			for(int i=0;i<node->child.Count();i++) {
				if(node->child[i]>nextnode) return;
				bfs_cov(totalcov,len,cov,node->child[i],inode,t,no2gnode,transfrag,nodecov,gpos,gno);
			}
		}
		else { // I have to advance to nextnode because the pattern forces me
			bfs_cov(totalcov,len,cov,nextnode,inode,t,no2gnode,transfrag,nodecov,gpos,gno);
		}
	}
	else { // n is before nextnode
		for(int i=0;i<node->child.Count();i++) {
			if(node->child[i]>nextnode) return;
			bfs_cov(totalcov,len,cov,node->child[i],inode,t,no2gnode,transfrag,nodecov,gpos,gno);
		}
	}
}


float push_max_flow(int gno,GVec<int>& path,GBitVec& istranscript,GPVec<CTransfrag>& transfrag,GPVec<CGraphnode>& no2gnode,
		GVec<float>& nodeflux,GBitVec& pathpat) { //,float& fragno) {

	int n=path.Count();
	GVec<int> node2path;
	node2path.Resize(gno,-1);
	for(int i=0;i<n;i++) {
		node2path[path[i]]=i;
		nodeflux.cAdd(0.0);
	}
	GVec<float> capacityleft;   // how many transcripts compatible to path enter node
	GVec<float> capacityright;  // how many transcripts compatible to path exit node
	capacityleft.Resize(n,0);
	capacityright.Resize(n,0);
	GVec<float> sumleft;        // how many transcripts enter node
	GVec<float> sumright;       // how many transcripts exit node
	sumleft.Resize(n,0);
	sumright.Resize(n,0);
	//GVec<CNodeCapacity> nodecap; // all node capacity percentages observed
	*
	 * struct CNodeCapacity {
	 *	int id;
	 *	bool left;
	 *	float perc;
	 *	CNodeCapacity(int nid=0,bool leftnode=false,float p=0): id(nid),left(leftnode),perc(p) {}
	 * };
	 *


	{ // DEBUG ONLY
		printTime(stderr);
		fprintf(stderr,"Start max flow algorithm for path ");
		printBitVec(pathpat);
		fprintf(stderr," :");
		for(int i=0;i<n;i++) fprintf(stderr," %d:%d",i,path[i]);
		fprintf(stderr,"\n");
		fprintf(stderr,"Used transcripts:");
		for(int i=0;i<transfrag.Count();i++) if(istranscript[i]) fprintf(stderr," %d(%f)",i,transfrag[i]->abundance);
		fprintf(stderr,"\n");
	}


	// compute capacities and sums for all nodes
	for(int i=1;i<n-1;i++) {
		int nt=no2gnode[path[i]]->trf.Count();
		for(int j=0;j<nt;j++) {
			int t=no2gnode[path[i]]->trf[j];
			if(transfrag[t]->abundance) {
				if(istranscript[t] || ((pathpat & transfrag[t]->pattern)==transfrag[t]->pattern)) { // transcript on path
					istranscript[t]=1;
					if(transfrag[t]->nodes[0]<path[i]) { // transfrag starts before this node
						sumleft[i]+=transfrag[t]->abundance;
						capacityleft[i]+=transfrag[t]->abundance;
					}
					if(transfrag[t]->nodes.Last()>path[i]) { // transfrag ends after this node
						sumright[i]+=transfrag[t]->abundance;
						capacityright[i]+=transfrag[t]->abundance;
					}
					/ when my incomplete transfrafs are not mapped to all possible intermediary nodes :
					if(transfrag[t]->nodes[0]==path[i]) { // transfrag starts at this node
						int n2=node2path[transfrag[t]->nodes.Last()];
						sumright[i]+=transfrag[t]->abundance;
						capacityright[i]+=transfrag[t]->abundance;
						for(int k=i+1;k<n2;k++) {
							sumleft[k]+=transfrag[t]->abundance;
							capacityleft[k]+=transfrag[t]->abundance;
							sumright[k]+=transfrag[t]->abundance;
							capacityright[k]+=transfrag[t]->abundance;
						}
					} else if(transfrag[t]->nodes.Last()==path[i]) {
						sumleft[i]+=transfrag[t]->abundance;
						capacityleft[i]+=transfrag[t]->abundance;
					}/

				}
				else { // transfrag not on path
					if(path[i]>transfrag[t]->nodes[0]) sumleft[i]+=transfrag[t]->abundance;
					if(path[i]<transfrag[t]->nodes.Last()) sumright[i]+=transfrag[t]->abundance;
				}
			}
		}
		if(!capacityleft[i]) return(0); /else {
			CNodeCapacity ncleft(i,true,capacityleft[i]/sumleft[i]);
			nodecap.Add(ncleft);
		}/
		if(!capacityright[i]) return(0); /else {
			CNodeCapacity ncright(i,false,capacityright[i]/sumright[i]);
			nodecap.Add(ncright);
		}/
	}

	//nodecap.Sort(nodecapCmp);


	{ // DEBUG ONLY
		for(int i=1;i<n-1;i++) {
			fprintf(stderr,"Node %d LEFT: capacity=%f total=%f ",path[i],capacityleft[i],sumleft[i]);
			if(sumleft[i]) fprintf(stderr,"perc=%f ",capacityleft[i]/sumleft[i]);
			else fprintf(stderr,"perc=n/a ");
			fprintf(stderr,"RIGHT: capacity=%f total=%f ",capacityright[i],sumright[i]);
			if(sumright[i]) fprintf(stderr,"perc=%f\n",capacityright[i]/sumright[i]);
			else fprintf(stderr,"perc=n/a\n");
		}
		fprintf(stderr,"Used transcripts:");
		for(int i=0;i<transfrag.Count();i++) if(istranscript[i]) fprintf(stderr," %d(%f)",i,transfrag[i]->abundance);
		fprintf(stderr,"\n");
	}


	// compute flow
	float prevflow=capacityleft[1];
	for(int i=1;i<n-1;i++) {
		float percleft=prevflow/sumleft[i];
		float percright=capacityright[i]/sumright[i];
		if(percright>percleft) { // more transfrags leave node
			percright=percleft;
		}
		prevflow=percright*sumright[i];
	}
	if(!prevflow) return(0);

	for(int i=n-2;i>0;i--) {
		nodeflux[i]=prevflow/sumright[i];
		capacityright[i]=prevflow;
		prevflow=nodeflux[i]*sumleft[i];
		capacityleft[i]=prevflow;
	}

	/ *for(int i=1;i<n-1;i++) { sumleft[i]=0;sumright[i]=0;} // just for printing purposes at stderr

	for(int z=0;z<nodecap.Count();z++) { // for all capacities starting from the smallest one
		int i=nodecap[z].id; // i is the path position of the node I am considering
		float capacity=capacityright[i];
		if(nodecap[z].left) capacity=capacityleft[i];

		fprintf(stderr,"Consider node[%d]=%d and capacity=%f\n",i,path[i],capacity);

		if(capacity) {
			int nt=no2gnode[path[i]]->trf.Count();
			for(int j=0;j<nt;j++) {
				int t=no2gnode[path[i]]->trf[j];
				if(istranscript[t] && transfrag[t]->abundance) {
					if((nodecap[z].left && transfrag[t]->nodes[0]<path[i]) ||    // transfrag doesn't start at node if I am considering left capapacity
							(!nodecap[z].left && transfrag[t]->nodes.Last()>path[i])) { // transfrag doesn't end at node if I am considering right capacity
						if(capacity>transfrag[t]->abundance) {
							capacity-=transfrag[t]->abundance;
							int tn=transfrag[t]->nodes.Count();
							for(int k=0;k<tn;k++) {
								int l=node2path[transfrag[t]->nodes[k]];
								if(k<tn-1) { capacityright[l]-=transfrag[t]->abundance; sumright[l]+=transfrag[t]->abundance;}
								if(k) { capacityleft[l]-=transfrag[t]->abundance; sumleft[l]+=transfrag[t]->abundance;}
							}
							transfrag[t]->abundance=0;
						}
						else {
							transfrag[t]->abundance-=capacity;
							if(transfrag[t]->abundance<epsilon) transfrag[t]->abundance=0;
							int tn=transfrag[t]->nodes.Count();
							for(int k=0;k<tn;k++) {
								int l=node2path[transfrag[t]->nodes[k]];
								if(k<tn-1) { capacityright[l]-=capacity;sumright[l]+=capacity;}
								if(k) { capacityleft[l]-=capacity; sumleft[l]+=capacity;}
							}
							capacity=0;
							break;
						}
					}
				}
			}
		}
	}

	for(int i=1;i<n-1;i++) {
		fprintf(stderr,"NODE[%d]= %d used LEFT: %f RIGHT: %f\n",i,path[i],sumleft[i],sumright[i]);
	}* /

	// * here I don't care what node I treat first
	for(int i=1;i<n-1;i++) if(capacityright[i]){
		int nt=no2gnode[path[i]]->trf.Count();
		for(int j=0;j<nt;j++) {
			int t=no2gnode[path[i]]->trf[j];
			if(istranscript[t] && transfrag[t]->abundance) {
				if(transfrag[t]->nodes[0]==path[i]) { // transfrag starts at this node
					if(capacityright[i]>transfrag[t]->abundance) {
						fprintf(stderr,"Update capacity of transfrag[%d] with value=%f to 0\n",t,transfrag[t]->abundance);
						capacityright[i]-=transfrag[t]->abundance;
						int n2=node2path[transfrag[t]->nodes.Last()];
						for(int k=i+1;k<n2;k++) {
							capacityright[k]-=transfrag[t]->abundance;
						}
						transfrag[t]->abundance=0;
					}
					else {
						fprintf(stderr,"Update capacity of transfrag[%d] with value=%f to %f\n",t,transfrag[t]->abundance,transfrag[t]->abundance-capacityright[i]);
						transfrag[t]->abundance-=capacityright[i];
						if(transfrag[t]->abundance<epsilon) transfrag[t]->abundance=0;
						int n2=node2path[transfrag[t]->nodes.Last()];
						for(int k=i+1;k<n2;k++) {
							capacityright[k]-=capacityright[i];
						}
						capacityright[i]=0;
						break;
					}
				}
			}
		}
	}

	// I only have to deal with source transfrag
	int nt=no2gnode[path[0]]->trf.Count();
	for(int j=0;j<nt;j++) {
		int t=no2gnode[path[0]]->trf[j];
		if(istranscript[t] && transfrag[t]->abundance) {
			fprintf(stderr,"Update capacity of transfrag[%d] with value=%f to %f\n",t,transfrag[t]->abundance,transfrag[t]->abundance-prevflow);
			transfrag[t]->abundance-=prevflow;
			if(transfrag[t]->abundance<epsilon) transfrag[t]->abundance=0;
		}
	} // /

	/ *
	{ // DEBUG ONLY
		fprintf(stderr,"Flow:");
		for(int i=0;i<n;i++)
			fprintf(stderr,"Used %f of node %d[%d]\n",nodeflux[i],i,path[i]);
		fprintf(stderr,"\nTranscript abundances");
		for(int i=0;i<transfrag.Count();i++) if(istranscript[i]) fprintf(stderr," %d(%f)",i,transfrag[i]->abundance);
		fprintf(stderr,"\n");
	}
	* /

	return(nodeflux[1]);

}


void bfs_trf_cov(float &totalcov,GBitVec& pattern,int len, float cov, int n, int inode, int t,GPVec<CGraphnode>& no2gnode,
		GPVec<CTransfrag>& transfrag, GIntHash<int> &gpos,int gno) {


	fprintf(stderr,"bfs[trf=%d n0=%d nl=%d n=%d cov=%f len=%d\n",t,transfrag[t]->nodes[0],transfrag[t]->nodes.Last(),n,cov,len);

	CGraphnode *node=no2gnode[n];

	if(n==transfrag[t]->nodes.Last()) {
		cov/=len;
		CPath p(pattern,cov);
		transfrag[t]->path.Add(p);
		totalcov+=cov;
		fprintf(stderr,"Done bfs cov=%f len=%d\n",cov,len);
		return;
	}

	int nextnode=transfrag[t]->nodes[inode];
	if(n==nextnode) { // n is not the last node in transfrag => I can increase inode without stepping over the boundaries
		inode++;
		nextnode=transfrag[t]->nodes[inode];

		while(nextnode==n+1) {
			inode++;
			if(inode==transfrag[t]->nodes.Count()) {
				cov/=len;
				CPath p(pattern,cov);
				transfrag[t]->path.Add(p);
				totalcov+=cov;
				fprintf(stderr,"Done bfs cov=%f len=%d\n",cov,len);
				return;
			}
			n=nextnode;
			nextnode=transfrag[t]->nodes[inode];
		}

		int *pos=gpos[edge(n,nextnode,gno)];
		if(!pos || !transfrag[t]->pattern[*pos]) { // there is no edge between n and nextnode -> I might still reach nextnode

			int nt=node->trf.Count();

			for(int i=0;i<node->child.Count();i++) {
				if(node->child[i]>nextnode) return;
				if(node->child[i]==nextnode || no2gnode[node->child[i]]->childpat[nextnode]) { // only if I can reach nextnode from this child I can go forward
					pattern[node->child[i]]=1;

					pos=gpos[edge(n,node->child[i],gno)];
					int trabund=0;
					if(pos) for(int f=0;f<nt;f++) {
						if(transfrag[node->trf[f]]->pattern[*pos]) trabund+=transfrag[node->trf[f]]->abundance;
					}

					bfs_trf_cov(totalcov,pattern,len+1,cov+trabund,node->child[i],inode,t,no2gnode,transfrag,gpos,gno);
					pattern[node->child[i]]=0;
				}
			}
		}
		else { // I have to advance to nextnode because the pattern forces me
			pattern[nextnode]=1;
			bfs_trf_cov(totalcov,pattern,len,cov,nextnode,inode,t,no2gnode,transfrag,gpos,gno);
			pattern[nextnode]=0;
		}
	}
	else { // n is before nextnode
		int nt=node->trf.Count();
		for(int i=0;i<node->child.Count();i++) {
			if(node->child[i]>nextnode) return;
			if(node->child[i]==nextnode || no2gnode[node->child[i]]->childpat[nextnode]) { // only if I can reach nextnode from this child I can go forward
				pattern[node->child[i]]=1;

				int *pos=gpos[edge(n,node->child[i],gno)];
				int trabund=0;
				if(pos) for(int f=0;f<nt;f++) {
					if(transfrag[node->trf[f]]->pattern[*pos]) trabund+=transfrag[node->trf[f]]->abundance;
				}

				bfs_trf_cov(totalcov,pattern,len+1,cov+trabund,node->child[i],inode,t,no2gnode,transfrag,gpos,gno);
				pattern[node->child[i]]=0;
			}
		}
	}
}
*/



// version of push_max_flow where I weight the incomplete transfrags
float push_max_flow(int gno,GVec<int>& path,GBitVec& istranscript,GPVec<CTransfrag>& transfrag,GPVec<CGraphnode>& no2gnode,
		GVec<float>& nodeflux,GBitVec& pathpat, GIntHash<int> &gpos) {

	int n=path.Count();
	GVec<int> node2path;
	node2path.Resize(gno,-1);
	for(int i=0;i<n;i++) {
		node2path[path[i]]=i;
		nodeflux.cAdd(0.0);
	}
	GVec<float> capacityleft;   // how many transcripts compatible to path enter node
	GVec<float> capacityright;  // how many transcripts compatible to path exit node
	capacityleft.Resize(n,0);
	capacityright.Resize(n,0);
	GVec<float> sumleft;        // how many transcripts enter node
	GVec<float> sumright;       // how many transcripts exit node
	sumleft.Resize(n,0);
	sumright.Resize(n,0);

	/*
	{ // DEBUG ONLY
		printTime(stderr);
		fprintf(stderr,"Start push max flow algorithm for path ");
		//printBitVec(pathpat);
		fprintf(stderr," :");
		for(int i=0;i<n;i++) fprintf(stderr," %d:%d",i,path[i]);
		fprintf(stderr,"\n");
		//fprintf(stderr,"Used transcripts:");
		//for(int i=0;i<transfrag.Count();i++) if(istranscript[i]) fprintf(stderr," %d(%f)",i,transfrag[i]->abundance);
		//fprintf(stderr,"\n");
	}
	*/

	// compute capacities and sums for all nodes
	for(int i=1;i<n-1;i++) {
		int nt=no2gnode[path[i]]->trf.Count();
		for(int j=0;j<nt;j++) {
			int t=no2gnode[path[i]]->trf[j];
			if(transfrag[t]->abundance) {
				if(istranscript[t] || ((pathpat & transfrag[t]->pattern)==transfrag[t]->pattern)) { // transcript on path
					istranscript[t]=1;

					//fprintf(stderr,"istranscript[%d] with abund=%f and path[%d]=%d and nodes[0]=%d and nodes[last]=%d\n",t,transfrag[t]->abundance,i,path[i],transfrag[t]->nodes[0],transfrag[t]->nodes.Last());

					if(i==1 || transfrag[t]->nodes[0]==path[i]) { // first time I encounter transfrag I have to set what abundance to use
						if(!transfrag[t]->real) { // if I still didn't solve transfrag
							transfrag[t]->usepath=-1;
							for(int p=0;p<transfrag[t]->path.Count();p++) {
								int *pos=gpos[edge(transfrag[t]->path[p].node,transfrag[t]->path[p].contnode,gno)];
								if(pos && pathpat[*pos]) {
									transfrag[t]->usepath=p; // this is path dependent
									break;
								}
							}
						}
					}


					if(transfrag[t]->nodes[0]<path[i]) { // transfrag starts before this node
						sumleft[i]+=transfrag[t]->abundance;
						if(transfrag[t]->real) capacityleft[i]+=transfrag[t]->abundance;
						else if(transfrag[t]->usepath>-1)
							capacityleft[i]+=transfrag[t]->abundance*transfrag[t]->path[int(transfrag[t]->usepath)].abundance;

						//fprintf(stderr,"add transfrag t=%d i=%d sumleft=%f capacityleft=%f\n",t,i,sumleft[i],capacityleft[i]);
					}
					if(transfrag[t]->nodes.Last()>path[i]) { // transfrag ends after this node
						sumright[i]+=transfrag[t]->abundance;
						if(transfrag[t]->real) capacityright[i]+=transfrag[t]->abundance;
						else if(transfrag[t]->usepath>-1)
							capacityright[i]+=transfrag[t]->abundance*transfrag[t]->path[int(transfrag[t]->usepath)].abundance;
						//fprintf(stderr,"add transfrag t=%d i=%d sumright=%f capacityright=%f\n",t,i,sumright[i],capacityright[i]);
					}
				}
				else { // transfrag not on path
					if(path[i]>transfrag[t]->nodes[0]) sumleft[i]+=transfrag[t]->abundance;
					if(path[i]<transfrag[t]->nodes.Last()) sumright[i]+=transfrag[t]->abundance;
				}
			}
		}

		/*
		{ // DEBUG ONLY
			fprintf(stderr,"Node %d LEFT: capacity=%f total=%f ",path[i],capacityleft[i],sumleft[i]);
			if(sumleft[i]) fprintf(stderr,"perc=%f ",capacityleft[i]/sumleft[i]);
			else fprintf(stderr,"perc=n/a ");
			fprintf(stderr,"RIGHT: capacity=%f total=%f ",capacityright[i],sumright[i]);
			if(sumright[i]) fprintf(stderr,"perc=%f\n",capacityright[i]/sumright[i]);
			else fprintf(stderr,"perc=n/a\n");
		}
		*/

		if(!capacityleft[i]) return(0);
		if(!capacityright[i]) return(0);

	}

	/*
	{ // DEBUG ONLY
		for(int i=1;i<n-1;i++) {
			fprintf(stderr,"Node %d LEFT: capacity=%f total=%f ",path[i],capacityleft[i],sumleft[i]);
			if(sumleft[i]) fprintf(stderr,"perc=%f ",capacityleft[i]/sumleft[i]);
			else fprintf(stderr,"perc=n/a ");
			fprintf(stderr,"RIGHT: capacity=%f total=%f ",capacityright[i],sumright[i]);
			if(sumright[i]) fprintf(stderr,"perc=%f\n",capacityright[i]/sumright[i]);
			else fprintf(stderr,"perc=n/a\n");
		}
		fprintf(stderr,"Used transcripts:");
		for(int i=0;i<transfrag.Count();i++) if(istranscript[i]) fprintf(stderr," %d(%f)",i,transfrag[i]->abundance);
		fprintf(stderr,"\n");
	}
	*/

	// compute flow
	float prevflow=capacityleft[1];
	for(int i=1;i<n-1;i++) {
		float percleft=prevflow/sumleft[i];
		float percright=capacityright[i]/sumright[i];
		if(percright>percleft) { // more transfrags leave node
			percright=percleft;
		}
		prevflow=percright*sumright[i];
	}
	if(!prevflow) return(0);

	for(int i=n-2;i>0;i--) {
		//fprintf(stderr,"i=%d sumright=%f prevflow=%f\n",i,sumright[i],prevflow);
		nodeflux[i]=prevflow/sumright[i];
		//fprintf(stderr,"nodeflux=%f\n",nodeflux[i]);
		capacityright[i]=prevflow;
		prevflow=nodeflux[i]*sumleft[i];
		//fprintf(stderr,"i=%d sumright=%f sumleft=%f prevflow=%f capacityright=%f nodeflux=%f\n",i,sumright[i],sumleft[i],prevflow,capacityright[i],nodeflux[i]);
		//capacityleft[i]=prevflow; // I don't use this
	}

	// * here I don't care what node I treat first
	for(int i=1;i<n-1;i++) if(capacityright[i]){
		int nt=no2gnode[path[i]]->trf.Count();
		for(int j=0;j<nt;j++) {
			int t=no2gnode[path[i]]->trf[j];
			if(istranscript[t] && transfrag[t]->abundance) {

				float trabundance=transfrag[t]->abundance;
				if(!transfrag[t]->real) {
					if(transfrag[t]->usepath>-1)
						trabundance=transfrag[t]->abundance*transfrag[t]->path[int(transfrag[t]->usepath)].abundance;
					else trabundance=0;
				}

				if(trabundance && transfrag[t]->nodes[0]==path[i]) { // transfrag starts at this node
					if(capacityright[i]>trabundance) {
						//fprintf(stderr,"Update capacity of transfrag[%d] with value=%f to 0\n",t,transfrag[t]->abundance);
						capacityright[i]-=trabundance;
						int n2=node2path[transfrag[t]->nodes.Last()];
						for(int k=i+1;k<n2;k++) {
							capacityright[k]-=trabundance;
						}
						transfrag[t]->abundance-=trabundance;
						if(transfrag[t]->abundance<epsilon) transfrag[t]->abundance=0;
						else if(!transfrag[t]->real) {
							transfrag[t]->path[int(transfrag[t]->usepath)].abundance=0;
							if(transfrag[t]->path.Count()-1 < 2) transfrag[t]->real=true;
							else {
								int np=0;
								for(int p=0;p<transfrag[t]->path.Count();p++)
									if(transfrag[t]->path[int(transfrag[t]->usepath)].abundance) np++;
								if(np<2) transfrag[t]->real=true;
							}
						}
					}
					else {
						//fprintf(stderr,"Update capacity of transfrag[%d] with value=%f to %f\n",t,transfrag[t]->abundance,transfrag[t]->abundance-capacityright[i]);
						transfrag[t]->abundance-=capacityright[i];
						if(transfrag[t]->abundance<epsilon) {
							transfrag[t]->abundance=0;
						}
						else if(!transfrag[t]->real) {
							//transfrag[t]->path[int(transfrag[t]->usepath)].abundance-=capacityright[i]; // not needed anymore because this stores proportions not actual abundances
							//if(transfrag[t]->path[int(transfrag[t]->usepath)].abundance<epsilon) {
							if(transfrag[t]->path[int(transfrag[t]->usepath)].abundance*transfrag[t]->abundance-capacityright[i]<epsilon) {
								transfrag[t]->path[int(transfrag[t]->usepath)].abundance=0;
								if(transfrag[t]->path.Count()-1 < 2) transfrag[t]->real=true;
								else {
									int np=0;
									for(int p=0;p<transfrag[t]->path.Count();p++)
										if(transfrag[t]->path[int(transfrag[t]->usepath)].abundance) np++;
									if(np<2) transfrag[t]->real=true;
								}
							}
						}

						int n2=node2path[transfrag[t]->nodes.Last()];
						for(int k=i+1;k<n2;k++) {
							capacityright[k]-=capacityright[i];
						}
						capacityright[i]=0;
						break;
					}
				}


			}
		}
	}

	// I only have to deal with source transfrag
	int nt=no2gnode[path[0]]->trf.Count();
	for(int j=0;j<nt;j++) {
		int t=no2gnode[path[0]]->trf[j];
		if(istranscript[t] && transfrag[t]->abundance) {
			//fprintf(stderr,"Update capacity of transfrag[%d] with value=%f to %f\n",t,transfrag[t]->abundance,transfrag[t]->abundance-prevflow);
			transfrag[t]->abundance-=prevflow;
			if(transfrag[t]->abundance<epsilon) transfrag[t]->abundance=0;
			break; // there is no point in updating more than one transfrag from source
		}
	} //*/

	/*
	{ // DEBUG ONLY
		fprintf(stderr,"Flow:");
		for(int i=0;i<n;i++)
			fprintf(stderr,"Used %f of node %d[%d]\n",nodeflux[i],i,path[i]);
		fprintf(stderr,"\nTranscript abundances");
		for(int i=0;i<transfrag.Count();i++) if(istranscript[i]) fprintf(stderr," %d(%f)",i,transfrag[i]->abundance);
		fprintf(stderr,"\n");
	}
	*/

	return(nodeflux[1]);

}


float push_guide_maxflow(int gno,GVec<int>& path,GBitVec& istranscript,GPVec<CTransfrag>& transfrag,GPVec<CGraphnode>& no2gnode,GBitVec& pathpat) {

	float guideabundance=0;

	int n=path.Count();
	GVec<int> node2path;
	node2path.Resize(gno,-1);
	for(int i=0;i<n;i++) {
		node2path[path[i]]=i;
	}
	GVec<float> capacityleft;   // how many transcripts compatible to path enter node
	GVec<float> capacityright;  // how many transcripts compatible to path exit node
	capacityleft.Resize(n,0);
	capacityright.Resize(n,0);
	GVec<float> sumleft;        // how many transcripts enter node
	GVec<float> sumright;       // how many transcripts exit node
	sumleft.Resize(n,0);
	sumright.Resize(n,0);

	/*
	{ // DEBUG ONLY
		printTime(stderr);
		fprintf(stderr,"Start push guide max flow algorithm for path ");
		//printBitVec(pathpat);
		fprintf(stderr," :");
		//for(int i=0;i<n;i++) fprintf(stderr," %d:%d",i,path[i]);
		for(int i=0;i<n;i++) fprintf(stderr," %d",path[i]);
		fprintf(stderr,"\n");
		//fprintf(stderr,"Used transcripts:");
		//for(int i=0;i<transfrag.Count();i++) if(istranscript[i]) fprintf(stderr," %d(%f)",i,transfrag[i]->abundance);
		//fprintf(stderr,"\n");
	}
	*/

	// compute capacities and sums for all nodes
	for(int i=1;i<n-1;i++) {
		int nt=no2gnode[path[i]]->trf.Count();
		for(int j=0;j<nt;j++) {
			int t=no2gnode[path[i]]->trf[j];
			if(transfrag[t]->abundance) {
				//fprintf(stderr,"Consider transcript %d with abundance %f for node %d\n",t,transfrag[t]->abundance,path[i]);
				if(istranscript[t] || ((pathpat & transfrag[t]->pattern)==transfrag[t]->pattern)) { // transcript on path
					istranscript[t]=1;
					//fprintf(stderr,"...on path\n");

					if(transfrag[t]->nodes[0]<path[i]) { // transfrag starts before this node
						sumleft[i]+=transfrag[t]->abundance;
						capacityleft[i]+=transfrag[t]->abundance;
					}
					if(transfrag[t]->nodes.Last()>path[i]) { // transfrag ends after this node
						sumright[i]+=transfrag[t]->abundance;
						capacityright[i]+=transfrag[t]->abundance;
					}
				}
				else { // transfrag not on path
					if(path[i]>transfrag[t]->nodes[0]) sumleft[i]+=transfrag[t]->abundance;
					if(path[i]<transfrag[t]->nodes.Last()) sumright[i]+=transfrag[t]->abundance;
				}
			}
		}

		/*
		{ // DEBUG ONLY
			fprintf(stderr,"Node %d LEFT: capacity=%f total=%f ",path[i],capacityleft[i],sumleft[i]);
			if(sumleft[i]) fprintf(stderr,"perc=%f ",capacityleft[i]/sumleft[i]);
			else fprintf(stderr,"perc=n/a ");
			fprintf(stderr,"RIGHT: capacity=%f total=%f ",capacityright[i],sumright[i]);
			if(sumright[i]) fprintf(stderr,"perc=%f\n",capacityright[i]/sumright[i]);
			else fprintf(stderr,"perc=n/a\n");
		}
		*/

		if(!capacityleft[i]) return(0);
		if(!capacityright[i]) return(0);


	}

	/*
	{ // DEBUG ONLY
		for(int i=1;i<n-1;i++) {
			fprintf(stderr,"Node %d LEFT: capacity=%f total=%f ",path[i],capacityleft[i],sumleft[i]);
			if(sumleft[i]) fprintf(stderr,"perc=%f ",capacityleft[i]/sumleft[i]);
			else fprintf(stderr,"perc=n/a ");
			fprintf(stderr,"RIGHT: capacity=%f total=%f ",capacityright[i],sumright[i]);
			if(sumright[i]) fprintf(stderr,"perc=%f\n",capacityright[i]/sumright[i]);
			else fprintf(stderr,"perc=n/a\n");
		}
		fprintf(stderr,"Used transcripts:");
		for(int i=0;i<transfrag.Count();i++) if(istranscript[i]) fprintf(stderr," %d(%f)",i,transfrag[i]->abundance);
		fprintf(stderr,"\n");
	}
	*/

	// compute flow
	float prevflow=capacityleft[1];
	for(int i=1;i<n-1;i++) {
		float percleft=prevflow/sumleft[i];
		float percright=capacityright[i]/sumright[i];
		if(percright>percleft) { // more transfrags leave node
			percright=percleft;
		}
		prevflow=percright*sumright[i];
	}
	if(!prevflow) return(0);

	for(int i=n-2;i>0;i--) {
		guideabundance+=no2gnode[path[i]]->cov*prevflow/sumright[i];
		prevflow=prevflow*sumleft[i]/sumright[i];
	}

	return(guideabundance);

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



float guidepushflow(int g,GVec<CGuide>& guidetrf,int gno,GBitVec& istranscript,GPVec<CTransfrag>& transfrag,
		GPVec<CGraphnode>& no2gnode,GVec<float>& nodeflux) {

	int n=guidetrf[g].trf->nodes.Count();
	GVec<int> node2path;
	node2path.Resize(gno,-1);
	//fprintf(stderr,"path: ");

	//fprintf(stderr,"Process guide %d with pattern: ",g);
	//printBitVec(guidetrf[g].trf->pattern);

	for(int i=0;i<n;i++) {
		//fprintf(stderr,"%d ",guidetrf[g].trf->nodes[i]);
		node2path[guidetrf[g].trf->nodes[i]]=i;
		nodeflux.cAdd(0.0);
	}

	GVec<float> capacityleft;	// how many transcripts compatible to path enter node
	GVec<float> capacityright;  // how many transcripts compatible to path exit node
	capacityleft.Resize(n,0);
	capacityright.Resize(n,0);
	GVec<float> sumleft;        // how many transcripts enter node
	GVec<float> sumright;       // how many transcripts exit node
	sumleft.Resize(n,0);
	sumright.Resize(n,0);

	// compute capacities and sums for all nodes
	for(int i=1;i<n-1;i++) {
		int pathi=guidetrf[g].trf->nodes[i];

		//fprintf(stderr,"Process node %d on path=%d\n",i,pathi);

		int nt=no2gnode[pathi]->trf.Count();
		for(int j=0;j<nt;j++) {
			int t=no2gnode[pathi]->trf[j];
			if(transfrag[t]->abundance) {
				if(istranscript[t] || ((guidetrf[g].trf->pattern & transfrag[t]->pattern)==transfrag[t]->pattern)) { // transcript on path
					istranscript[t]=1;
					//fprintf(stderr,"Transcript %d of abund=%f is part of guides: %d",t,transfrag[t]->abundance,g);
					// check if there are other guides sharing this transcript so that I can allocate proportionally to guide abundances
					float totalcov=guidetrf[g].trf->abundance;
					for(int r=g-1;r>=0;r--) if((guidetrf[r].trf->pattern & transfrag[t]->pattern)==transfrag[t]->pattern) {
						totalcov+=guidetrf[r].trf->abundance;
						//fprintf(stderr," %d",r);
					}
					float prop=1;
					if(totalcov>guidetrf[g].trf->abundance) prop=guidetrf[g].trf->abundance/totalcov;

					//fprintf(stderr," prop=%f\n",prop);

					transfrag[t]->usepath=prop*transfrag[t]->abundance;

					if(transfrag[t]->nodes[0]<pathi) { // transfrag starts before this node
						sumleft[i]+=transfrag[t]->usepath;
						capacityleft[i]+=transfrag[t]->usepath;
					}
					if(transfrag[t]->nodes.Last()>pathi) { // transfrag ends after this node
						sumright[i]+=transfrag[t]->usepath;
						capacityright[i]+=transfrag[t]->usepath;
					}
				}
				else { // transfrag not on path
					if(pathi>transfrag[t]->nodes[0]) sumleft[i]+=transfrag[t]->abundance;
					if(pathi<transfrag[t]->nodes.Last()) sumright[i]+=transfrag[t]->abundance;
				}
			}
		}

		//fprintf(stderr,"capacityleft[%d]=%f capacityright[%d]=%f\n",i,capacityleft[i],i,capacityright[i]);

		if(!capacityleft[i]) return(0);
		if(!capacityright[i]) return(0);
	}

	// compute flow -> deal with rounding errors here
	float prevflow=capacityleft[1];
	for(int i=1;i<n-1;i++) {
		float percleft=prevflow/sumleft[i];
		float percright=capacityright[i]/sumright[i];
		if(percright>percleft) { // more transfrags leave node
			percright=percleft;
		}
		prevflow=percright*sumright[i];
	}
	if(!prevflow) return(0);

	for(int i=n-2;i>0;i--) {
		nodeflux[i]=prevflow/sumright[i];
		if(nodeflux[i]>1) nodeflux[i]=1; // because of rounding errors this could become more than 1
		capacityright[i]=prevflow;
		prevflow=prevflow*sumleft[i]/sumright[i];
		//capacityleft[i]=prevflow;
	}

	for(int i=1;i<n-1;i++) {
		int pathi=guidetrf[g].trf->nodes[i];
		int nt=no2gnode[pathi]->trf.Count();
		for(int j=0;j<nt;j++) {
			int t=no2gnode[pathi]->trf[j];
			if(istranscript[t] && transfrag[t]->abundance) {
				float trabundance=transfrag[t]->usepath;

				if(transfrag[t]->nodes[0]==pathi) { // transfrag starts at this node

					if(capacityright[i]>trabundance) {
						//fprintf(stderr,"Update capacity of transfrag[%d] with value=%f to 0\n",t,transfrag[t]->abundance);
						//if(!transfrag[t]->real) fprintf(stderr,"Trf[%d] not real.\n",t);
						capacityright[i]-=trabundance;
						int n2=node2path[transfrag[t]->nodes.Last()];
						for(int k=i+1;k<n2;k++) {
							capacityright[k]-=trabundance;
						}
						transfrag[t]->abundance-=trabundance;
						if(transfrag[t]->abundance<epsilon) transfrag[t]->abundance=0;
						else transfrag[t]->usepath=0; // I only need to delete this if I still have abundance left in transfrag
					}
					else if(capacityright[i]) {
						//fprintf(stderr,"Update capacity of transfrag[%d] with value=%f to %f\n",t,transfrag[t]->abundance,transfrag[t]->abundance-capacityright[i]);
						//if(!transfrag[t]->real) fprintf(stderr,"Trf[%d] not real.\n",t);
						transfrag[t]->abundance-=capacityright[i];
						if(transfrag[t]->abundance<epsilon) transfrag[t]->abundance=0;
						else {
							transfrag[t]->usepath-=capacityright[i];
							if(transfrag[t]->usepath<epsilon)
								transfrag[t]->usepath=0;
						}
						int n2=node2path[transfrag[t]->nodes.Last()];
						for(int k=i+1;k<n2;k++) {
							capacityright[k]-=capacityright[i];
						}
						capacityright[i]=0;
						// break; I can not break because I might still have transfrag that need the cleanup
					}

				}
			}
		}
	}

	// I only have to deal with source transfrag
	int nt=no2gnode[guidetrf[g].trf->nodes[0]]->trf.Count();
	for(int j=0;j<nt;j++) {
		int t=no2gnode[guidetrf[g].trf->nodes[0]]->trf[j];
		if(istranscript[t] && transfrag[t]->abundance) {
			//fprintf(stderr,"Update capacity of transfrag[%d] with value=%f to %f\n",t,transfrag[t]->abundance,transfrag[t]->abundance-prevflow);
			transfrag[t]->abundance-=prevflow;
			if(transfrag[t]->abundance<epsilon) transfrag[t]->abundance=0;
			break;
		}
	} //*/

	return(nodeflux[1]);
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

float store_transcript(GList<CPrediction>& pred,GVec<int>& path,GVec<float>& nodeflux,GVec<float>& nodecov,
		GPVec<CGraphnode>& no2gnode,int& geneno,bool& first,int strand,int gno,GIntHash<int>& gpos, bool& included,
		GBitVec& prevpath, BundleData *bdata=NULL, //float fragno, char* id=NULL) {
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
		if(!prevpath[path[i]]) { // if I can find a node that was not included previously in any path then this is a new path
			included=false;
			prevpath[path[i]]=1;
		}
		int *pos=gpos[edge(path[i-1],path[i],gno)]; // if I can find an edge that was not included in any previous path then this is a new path
		if(i && (!pos || !prevpath[*pos])) {
			included=false;
			if(pos) prevpath[*pos]=1;
		}

		CGraphnode *node=no2gnode[path[i]];

		/*
		// normal:
		float usedcov=node->cov; // this is all the node coverage
		if(node->cov) {
			if(node->capacity) usedcov*=nodeflux[i]/node->capacity;
			else if(prevnode && prevnode->capacity) usedcov*=nodeflux[i-1]/prevnode->capacity;
			else {
				float rate=find_capacity(i,path,nodeflux,no2gnode);
				usedcov*=rate;
			}
		}
		nodecov[path[i]]-=usedcov/(node->end-node->start+1); // this is what's left from the node coverage
		*/

		// push
		float usedcov=nodecov[path[i]]*nodeflux[i]*(node->end-node->start+1);

		//fprintf(stderr,"usedcov=%f for nodecov[path[%d]]=%f nodeflux[%d]=%f node->end=%d node->start=%d\n",usedcov,i,nodecov[path[i]],i,nodeflux[i],node->end,node->start);


		nodecov[path[i]]*=(1-nodeflux[i]); // don't allow this to be less than 0

		if(t && (node->end<t->start || lastex)) { // I am skipping the nodes that do not overlap the guide so that I don't add them up to coverage
			prevnode=node; continue;
		}

		uint nodestart=node->start;
		uint nodeend=node->end;

		if(t) { // I am adjusting the start/end of the exon but shouldn't I also adjust the coverage? -> I added two ifs below
			float firstprop=0;
			float lastprop=0;
			if(firstex) {
				if(nodestart<t->start) {
					float leftcov=0;
					float rightcov=0;
					for(uint p=node->start;p<t->start;p++) {
						int j=p-bdata->start;
						leftcov+=bdata->bpcov[1][j]-bdata->bpcov[2-2*strand][j];
					}
					for(uint p=t->start;p<=node->end;p++) {
						int j=p-bdata->start;
						rightcov+=bdata->bpcov[1][j]-bdata->bpcov[2-2*strand][j];
					}
					if(leftcov) firstprop=leftcov/(leftcov+rightcov);
				}
				nodestart=t->start;
			}
			if(t->end<=nodeend || i==path.Count()-2) {
				lastex=true;
				if(t->end<node->end) {
					//usedcov*=(t->end-nodestart+1)/(nodeend-nodestart+1); // this way I am keeping coverage proportions right
					float leftcov=0;
					float rightcov=0;
					for(uint p=node->start;p<=t->end;p++) {
						int j=p-bdata->start;
						leftcov+=bdata->bpcov[1][j]-bdata->bpcov[2-2*strand][j];
					}
					for(uint p=t->end+1;p<=node->end;p++) {
						int j=p-bdata->start;
						rightcov+=bdata->bpcov[1][j]-bdata->bpcov[2-2*strand][j];
					}
					if(rightcov) lastprop=rightcov/(leftcov+rightcov);
				}
				nodeend=t->end;
			}
			if(firstprop || lastprop) {
				usedcov-=(firstprop +lastprop)*usedcov;
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

		//fprintf(stderr,"Add usedcov=%g to excov=%g = %g\n",usedcov,excov,excov+usedcov);
		excov+=usedcov;

		//if(node->cov) fragno+=node->frag*usedcov/node->cov;

		prevnode=node;
	}

	/*
	{ // DEBUG ONLY
		fprintf(stderr,"Predicted transcript cov=%f usedcov=%f len=%d path.count=%d\n",cov/len, cov,len,path.Count());
		if(t) fprintf(stderr,"Ref_id=%s\n",t->getID());
	}
	*/

	// Add last exon coverage
	if(prevnode) { // compute exon coverage
		excov/=exons.Last().end-exons.Last().start+1;
		exoncov.Add(excov);
	}
	if(len) cov/=len;


	//if(t || ((!included || sensitivitylevel) && cov>=readthr && len>=mintranscriptlen)) { // store transcript here
	if(t || (cov>=readthr && len>=mintranscriptlen)) { // store transcript here
	// //if(id || ( cov>=readthr && len>=mintranscriptlen)) { // store transcript here
		char sign='-';
		if(strand) { sign='+';}
		if(first) { geneno++;}
		//CPrediction *p=new CPrediction(geneno-1, id, exons[0].start, exons.Last().end, cov, sign, fragno, len);
		//if(t) fprintf(stderr,"store prediction with start=%d and end=%d\n",exons[0].start, exons.Last().end);
		float gcov=cov;

		if (t && t->uptr) {
			RC_TData &td = *(RC_TData*) (t->uptr);
			td.in_bundle=3;
			//fprintf(stderr,"st guide %s is stored\n",t->getID());
		}

		/*
		// this was up to version 1.2.1 -> I am not sure about keeping it
		if(t && t->exons.Count()==1) { // if single exon
			RC_TData* tdata=(RC_TData*)(t->uptr);
			if(len) gcov=(tdata->t_exons[0])->movlcount/len;
			if(cov<gcov) gcov=cov;
		}
		*/

		/*
		{ // DEBUG ONLY
			fprintf(stderr,"Store transcript with coverage %f \n",gcov);
			fprintf(stderr,"And exons cov:");
			for(int e=0;e<exons.Count();e++) fprintf(stderr," %g",exoncov[e]);
			fprintf(stderr,"\n");
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

// pred[np] = prediction to update; path = the nodes in the prediction or the nodes used from the prediction?
// nodeflux = in store transcript: quantity of transfrags that's used going out from each node;
// nodeflux = here: proportion of the node that is used
// nodecov = coverage of nodes
void update_guide_pred(GList<CPrediction>& pred,int np, GVec<int>& path,GVec<float>& nodeflux,GVec<float>& nodecov,
		GPVec<CGraphnode>& no2gnode,int gno,bool update) {

	/*
	{ // DEBUG ONLY
		fprintf(stderr,"Update guide ");
		if(pred[np]->t_eq) fprintf(stderr,"%s ",pred[np]->t_eq->getID());
		fprintf(stderr,"pred[%d]:",np);
		for(int j=0;j<pred[np]->exons.Count();j++) fprintf(stderr," %d-%d",pred[np]->exons[j].start,pred[np]->exons[j].end);
		fprintf(stderr,"\n");
	}
	*/

	int e=0;
	int nex=pred[np]->exons.Count();

	for(int i=0;i<path.Count();i++) {
		if(path[i] && path[i]<gno-1) { // not first or last node in graph
			CGraphnode *node=no2gnode[path[i]];
			float addcov=nodeflux[i]*nodecov[path[i]];
			if(update) nodecov[path[i]]-=addcov;
			// addcov/=node->len(); this is not necessary because it was already divided to node length
			while(e<nex && pred[np]->exons[e].end<node->start) e++; // if exon end before start of node skip it
			while(e<nex && pred[np]->exons[e].start<=node->end) { // there is overlap between exon and node
				int ovplen=node->overlapLen(pred[np]->exons[e].start, pred[np]->exons[e].end);
				float excov=addcov*ovplen;
				pred[np]->exoncov[e]+=excov/pred[np]->exons[e].len();
				pred[np]->cov+=excov/pred[np]->tlen;
				e++;
				//fprintf(stderr,"guide=%d exoncov[%d]=%g\n",np,e,pred[np]->exoncov[e]);
			}
		}
	}
}



int store_guide_transcript(GList<CPrediction>& pred,GVec<int>& path,GVec<float>& nodeflux,GVec<float>& nodecov,
		GPVec<CGraphnode>& no2gnode,int& geneno,bool& first,int gno, GffObj* t,bool update) {

	// first create the prediction based on the GffObj and then update it's coverage
	GVec<GSeg> exons;
	GVec<float> exoncov;
	int len=0;

	for(int i=0;i<t->exons.Count();i++) {
		exons.Add(t->exons[i]);
		exoncov.cAdd(0.0);
		len+=t->exons[i]->len();
	}

	if(first) { geneno++;first=false;}
	int np=pred.Count();
	CPrediction *p=new CPrediction(geneno-1, t, exons[0].start, exons.Last().end, 0, t->strand, len);
	p->exons=exons;
	p->exoncov=exoncov;
	pred.Add(p);

	if (t && t->uptr) {
		RC_TData &td = *(RC_TData*) (t->uptr);
		td.in_bundle=3;
		//fprintf(stderr,"sg guide %s is stored\n",t->getID());
	}

	update_guide_pred(pred,np,path,nodeflux,nodecov,no2gnode,gno,update);

	return(np);
}

void parse_trf(int maxi,int gno,int edgeno, GIntHash<int> &gpos,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,
		GBitVec& compatible,	int& geneno,bool first,int strand,GList<CPrediction>& pred,GVec<float>& nodecov,
		GBitVec& istranscript,GBitVec& removable,GBitVec& usednode,float maxcov,GBitVec& prevpath) {

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
	 	 fprintf(stderr,"\n\n***Start parse_trf with maxi=%d and cov=%f\n",maxi,nodecov[maxi]);
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


	 			 /*if(EM) flux=max_flow_EM(gno,path,istranscript,transfrag,no2gnode,nodeflux,pathpat);
	 			 else if(weight)
	 				 	 //flux=weight_max_flow_EM(gno,path,istranscript,transfrag,no2gnode,nodeflux,pathpat);
	 				 	 flux=weight_max_flow(gno,path,istranscript,transfrag,no2gnode,nodeflux,pathpat);
	 			 else */

	 			 // normal:
	 			  //flux=max_flow(gno,path,istranscript,transfrag,no2gnode,nodeflux,pathpat);
	 			 // push:
	 			 flux=push_max_flow(gno,path,istranscript,transfrag,no2gnode,nodeflux,pathpat,gpos);

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

	 //bool cont=true;

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

		 if(included || cov<isofrac*maxcov) {
			 /*
			 if(sensitivitylevel) usednode[maxi]=1;
			 else usednode = usednode | prevpath;
			 */
			 usednode[maxi]=1;
			 //fprintf(stderr,"used maxi=%d\n",maxi);
			 maxi=0;
			 maxcov=0;
			 //cont=false;

		 }
		 else if(cov>maxcov) maxcov=cov;
	 }
	 else {
		 /*
		 if(sensitivitylevel) usednode[maxi]=1; // start at different locations in graph
		 else {
			 usednode = usednode | prevpath;
			 usednode = usednode | pathpat;
		 }
		 */
		 usednode[maxi]=1;

		 //fprintf(stderr,"set usednode of %d\n",maxi);

		 maxi=0;
		 maxcov=0;
		 //cont=false;
	 }

	 // Node coverages:
	 for(int i=1;i<gno;i++)
		 if(!usednode[i] && nodecov[i]>nodecov[maxi]) maxi=i;

	 //fprintf(stderr," maxi=%d nodecov=%f\n",maxi,nodecov[maxi]);

	 //if(nodecov[maxi]>=readthr && (!specific || cont)) { // if I still have nodes that are above coverage threshold
	 if(nodecov[maxi]) { // if I still have nodes to use

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
		 parse_trf(maxi,gno,edgeno,gpos,no2gnode,transfrag,compatible,geneno,first,strand,pred,nodecov,istranscript,removable,usednode,maxcov,prevpath);
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

CTransfrag *find_guide_partial_pat(GffObj *guide,GPVec<CGraphnode>& no2gnode,int gno,int edgeno,GIntHash<int> &gpos,GVec<int>& olen,int &olensum) {

	CTransfrag *trguide=NULL;
	GBitVec guidepat(gno+edgeno);
	GVec<int> nodes;

	olensum=0;

	int nex=guide->exons.Count();
	int e=0;
	bool terminal=false;
	int lastnode=0; // lastnode intersected by guide
	for(int i=1;i<gno-1;i++) {
		CGraphnode *inode=no2gnode[i];
		int len=inode->overlapLen(guide->exons[e]->start,guide->exons[e]->end);
		//fprintf(stderr,"...node %d-%d and guide exon %d-%d overlap with len=%d\n",no2gnode[i]->start,no2gnode[i]->end,guide->exons[e]->start,guide->exons[e]->end,len);
		if(len) { // if node i intersects exon e of guide
			guidepat[i]=1;
			nodes.Add(i);
			olen.Add(len);
			olensum+=len;
			if(lastnode && inode->parentpat[lastnode]) { // the guide has intersected another node in the graph which was i's parent
				if(guide->exons[e]->start==inode->start) {
					if(terminal) { // lastnode ends the same as guide's exon
						int *pos=gpos[edge(lastnode,i,gno)];
						if(pos) guidepat[*pos]=1;
					}
				}
				else if(guide->exons[e]->start<inode->start) {
					if(no2gnode[lastnode]->end+1==inode->start) { // the previous lastnode comes right before this one
						int *pos=gpos[edge(lastnode,i,gno)];
						if(pos) guidepat[*pos]=1;
					}
				}
			}
			lastnode=i;
			terminal=false;
			if(guide->exons[e]->end<=inode->end) {
				if(guide->exons[e]->end==inode->end) terminal=true;
				e++;
				if(e==nex) break; // I've seen all guide's exons
			}
		}
		else if(inode->start>guide->exons[e]->end){ // stay with the same node if node comes after exon
			e++;
			i--;
			if(e==nex) break; // I've seen all guide's exons
		}
	}

	if(olensum) { // guide intersects at least one node in graph
		trguide=new CTransfrag(nodes,guidepat,0,false);
	}

	return(trguide);
}


void process_refguides(int gno,int edgeno,GIntHash<int>& gpos,int& lastgpos,GPVec<CGraphnode>& no2gnode,
		GPVec<CTransfrag>& transfrag,int s,GPVec<GffObj>& guides,GVec<CGuide>& guidetrf,BundleData *bdata) {

	char strand='-';
	if(s) strand='+';

	// find guides' patterns

	for(int g=0;g<guides.Count();g++) {
		//fprintf(stderr,"Consider guide[%d out of %d] %s\n",g,guides.Count(),guides[g]->getID());
		if((guides[g]->strand==strand || guides[g]->strand=='.') && ((RC_TData*)(guides[g]->uptr))->in_bundle>=2 && (guides[g]->overlap(no2gnode[1]->start,no2gnode[gno-2]->end))) {
			CTransfrag *trguide=find_guide_pat(guides[g],no2gnode,gno,edgeno,gpos);
			if(trguide) { // the guide can be found among the graph nodes
				//CGuide newguide(trguide,guides[g]);
				CGuide newguide(trguide,g);
				guidetrf.Add(newguide);

				/*
				{ // DEBUG ONLY
					fprintf(stderr,"Added guidetrf %d with ID=%s overlapping transcript interval %d - %d with nodes:",g,guides[g]->getID(),no2gnode[1]->start,no2gnode[gno-2]->end);
					for(int i=0;i<trguide->nodes.Count();i++) fprintf(stderr," %d",trguide->nodes[i]);
					fprintf(stderr,"\n");
					//fprintf(stderr,"s=%d strand = %c trguide[%d]=",s,strand,g);
					//printBitVec(trguide->pattern);
					//fprintf(stderr,"\n");
				}
				*/
			}
		}
	}

	// compute guides' abundances
	/* --> this is too simple; I should do it based on flow!
	for(int t=0;t<transfrag.Count();t++) {
		//if(transfrag[t]->nodes.Count()> 1)
		for(int g=0;g<guidetrf.Count();g++) {
			if(((transfrag[t]->pattern) & guidetrf[g].trf->pattern) == transfrag[t]->pattern) {
				guidetrf[g].trf->abundance+=transfrag[t]->abundance;
			}
		}
	}
	*/

	guidetrf.Sort(guideCmp);

	int g=0;
	while(g<guidetrf.Count()) {
		// check if guide is included in previous guide paths
		int p=0;

		//fprintf(stderr,"Process guide=%s\n",guides[guidetrf[g].g]->getID());

		while(p<g) { // here I assume that the guides are sorted by pattern
			//CTransfrag guideg=guidetrf[g];
			//CTransfrag guidep=guidetrf[p];
			if((guidetrf[g].trf->pattern & guidetrf[p].trf->pattern)==guidetrf[g].trf->pattern) {
				guidetrf[g].trf->real=false;  // this marks a guide that is included in another one -> I might want to swap them if they have the same abundance
				/*
				if(!complete) { // if guides are incomplete exclude the ones that are included into the more complete ones
					GFREE(guidetrf[g].trf);
					guidetrf.Delete(g);
					break;
				}
				*/
				break;
			}
			p++;
		}

		// find if first node of guidetrf extends to source
		int nodei=guidetrf[g].trf->nodes[0]; // first node of guide
		if(nodei) { // if source it isn't already in guide's pattern -> it shouldn't be

			int *pos=gpos[edge(0,nodei,gno)];
			if(!pos) { // here is where I need lastgpos
				//fprintf(stderr,"Add source link to position %d for guide=%d\n",no2gnode[nodei]->start,guidetrf[g].g);
				int key=edge(0,nodei,gno);
				gpos.Add(key,lastgpos);
				lastgpos++;
				pos=gpos[key];
			}

			if(includesource)
				guidetrf[g].trf->nodes.Insert(0,0); // I need to comment this if I need path not to include the source
			guidetrf[g].trf->pattern[0]=1;
			guidetrf[g].trf->pattern[*pos]=1;

			bool sourcestart=false; // assume there is no extension to source
			CGraphnode *inode=no2gnode[nodei];
			float leftcov=0;
			float rightcov=0;
			int previ=-1;
			for(int p=0;p<inode->parent.Count();p++) {
				if(!inode->parent[p]) { sourcestart=true; break;} // found link to source -> break from path
				if(no2gnode[inode->parent[p]]->end==inode->start-1) { // no splice site between parent and inode
					previ=inode->parent[p];
					if(inode->start>=guides[guidetrf[g].g]->start) {
						rightcov=inode->cov/inode->len();
						leftcov=no2gnode[previ]->cov/no2gnode[previ]->len();
					}
				}
			}
			if(!sourcestart) { // source is not the parent of first node

				// if guide starts inside node -> I need to compute abundances
				if(inode->start<guides[guidetrf[g].g]->start) {
					for(uint p=inode->start;p<guides[guidetrf[g].g]->start;p++) {
						int j=p-bdata->start;
						leftcov+=bdata->bpcov[1][j]-bdata->bpcov[2-2*s][j];
					}
					leftcov/=(guides[guidetrf[g].g]->start-inode->start);
					uint thisend=inode->end+1;
					if(guides[guidetrf[g].g]->end<inode->end) thisend=guides[guidetrf[g].g]->end+1;
					for(uint p=guides[guidetrf[g].g]->start;p<thisend;p++) {
						int j=p-bdata->start;
						rightcov+=bdata->bpcov[1][j]-bdata->bpcov[2-2*s][j];
					}
					rightcov/=(thisend-guides[guidetrf[g].g]->start);
				}

				float maxabund=trthr;
				if(rightcov>leftcov) {
					maxabund=rightcov-leftcov;
				}

				// add source to nodei transfrag
				GVec<int> nodes;
				nodes.cAdd(0);
				nodes.Add(nodei);
				GBitVec trpat(gno+edgeno);
				trpat[0]=1;
				trpat[nodei]=1;
				trpat[*pos]=1;
				CTransfrag *tr=new CTransfrag(nodes,trpat,maxabund);
				transfrag.Add(tr);

				// add source among maxnode parents
				inode->parent.Insert(0,0);
				inode->parentpat[*pos]=1; // source should be already among parents of maxnode but not the edge to source
			}
		}

		// find if I can extend guidetrf to sink
		nodei=guidetrf[g].trf->nodes.Last();
		if(nodei != gno-1) { // sink is not already present in guide pattern
			int sink=gno-1;
			guidetrf[g].trf->nodes.Add(sink);
			guidetrf[g].trf->pattern[sink]=1;
			int *pos=gpos[edge(nodei,sink,gno)];
			if(!pos) {
				//fprintf(stderr,"Add sink link from position %d for guide=%d\n",no2gnode[nodei]->end,guidetrf[g].g);
				int key=edge(nodei,sink,gno);
				gpos.Add(key,lastgpos);
				lastgpos++;
				pos=gpos[key];
			}
			guidetrf[g].trf->pattern[*pos]=1;

			bool sinkend=false;
			CGraphnode *inode=no2gnode[nodei];
			float leftcov=0;
			float rightcov=0;
			int nexti=-1;
			for(int c=0;c<inode->child.Count();c++) {
				if(inode->child[c]==gno-1) { sinkend=true; break;} // found link to source
				if(no2gnode[inode->child[c]]->start==inode->end+1) { // no splice site
					nexti=inode->child[c];
					if(inode->end<=guides[guidetrf[g].g]->end) {
						leftcov=inode->cov/inode->len();
						rightcov=no2gnode[nexti]->cov/no2gnode[nexti]->len();
					}
				}
			}

			if(!sinkend) { // if no path to sink was found I need to add one

				// node is longer than guide
				if(inode->end>guides[guidetrf[g].g]->end) {
					uint thisstart=inode->start;
					if(guides[guidetrf[g].g]->start>inode->start) thisstart=guides[guidetrf[g].g]->start;
					for(uint p=thisstart;p<=guides[guidetrf[g].g]->end;p++) {
						int j=p-bdata->start;
						leftcov+=bdata->bpcov[1][j]-bdata->bpcov[2-2*s][j];
					}
					leftcov/=(guides[guidetrf[g].g]->end-thisstart+1);
					for(uint p=guides[guidetrf[g].g]->end;p<=inode->end;p++) {
						int j=p-bdata->start;
						rightcov+=bdata->bpcov[1][j]-bdata->bpcov[2-2*s][j];
					}
					rightcov/=(inode->end-guides[guidetrf[g].g]->end+1);
				}

				float maxabund=trthr;
				if(rightcov<leftcov) {
					maxabund=leftcov-rightcov;
				}

				// add maxnode to sink transfrag
				GVec<int> nodes;
				nodes.Add(nodei);
				nodes.Add(sink);
				GBitVec trpat(gno+edgeno);
				trpat[nodei]=1;
				trpat[sink]=1;
				trpat[*pos]=1;
				//fprintf(stderr,"introduce node from %d to sink=%d wih abundance=%g where pos=%d\n",nodei,sink,maxabund,*pos);
				CTransfrag *tr=new CTransfrag(nodes,trpat,maxabund);
				transfrag.Add(tr);

				// add sink among maxnode children
				inode->child.Add(sink);
				inode->childpat[*pos]=1; // source should be already among parents of maxnode but not the edge to source
			}

		}

		/* previous source code --> I don't think it's that correct

		// find if first node of guidetrf extends to source ---> I probably shouldn't extend to source, I should just cut here -YES!!!
		int lasti=guidetrf[g].trf->nodes[0]; // first node of guide
		bool sourcestart=true;
		CGraphnode *inode=no2gnode[lasti];
		GVec<int> extendpath;
		float maxabund=trthr; // maximum abundance of source to node transfrag
		if(inode->start<guides[guidetrf[g].g]->start) {
			float leftcov=0;
			for(uint p=inode->start;p<guides[guidetrf[g].g]->start;p++) {
				int j=p-bdata->start;
				leftcov+=bdata->bpcov[1][j]-bdata->bpcov[2-2*s][j];
			}
			leftcov/=(guides[guidetrf[g].g]->start-inode->start);
			float rightcov=0;
			uint thisend=inode->end+1;
			if(guides[guidetrf[g].g]->end<inode->end) thisend=guides[guidetrf[g].g]->end+1;
			for(uint p=guides[guidetrf[g].g]->start;p<thisend;p++) {
				int j=p-bdata->start;
				rightcov+=bdata->bpcov[1][j]-bdata->bpcov[2-2*s][j];
			}
			rightcov/=(thisend-guides[guidetrf[g].g]->start);
			if(rightcov-leftcov>maxabund) maxabund=rightcov-leftcov;
		}
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
						if(transfrag[t]->nodes[0]<lasti) { // transfrag doesn't start at this node (in or through transfrag)
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
					if(includesource) tmpextend.cAdd(0); // I need to comment this if I need path not to include the source
					tmpextend.Reverse();
					tmpextend.Add(guidetrf[g].trf->nodes);
					guidetrf[g].trf->nodes.Clear();
					guidetrf[g].trf->nodes.Add(tmpextend);
					guidetrf[g].trf->pattern[0]=1;
					guidetrf[g].trf->pattern[*pos]=1;
				}
				else {
					if(includesource)
						guidetrf[g].trf->nodes.Insert(0,source); // I need to comment this if I need path not to include the source
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
			if(inode->end>guides[guidetrf[g].g]->end) {
				float leftcov=0;
				uint thisstart=inode->start;
				if(guides[guidetrf[g].g]->start>inode->start) thisstart=guides[guidetrf[g].g]->start;
				for(uint p=thisstart;p<=guides[guidetrf[g].g]->end;p++) {
					int j=p-bdata->start;
					leftcov+=bdata->bpcov[1][j]-bdata->bpcov[2-2*s][j];
				}
				leftcov/=(guides[guidetrf[g].g]->end-thisstart+1);
				float rightcov=0;
				for(uint p=guides[guidetrf[g].g]->end;p<=inode->end;p++) {
					int j=p-bdata->start;
					rightcov+=bdata->bpcov[1][j]-bdata->bpcov[2-2*s][j];
				}
				rightcov/=(inode->end-guides[guidetrf[g].g]->end+1);
				if(leftcov-rightcov>maxabund) maxabund=leftcov-rightcov;
			}
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
							if(transfrag[t]->nodes.Last()>lasti) { // transfrag doesn't end at this node (out or through transfrag)
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
					guidetrf[g].trf->pattern[*pos]=1;
				}

				// add maxnode to sink transfrag
				GVec<int> nodes;
				nodes.Add(maxnode);
				nodes.Add(sink);
				GBitVec trpat(gno+edgeno);
				trpat[maxnode]=1;
				trpat[sink]=1;
				trpat[*pos]=1;
				//fprintf(stderr,"introduce node from %d to sink=%d wih abundance=%g where pos=%d\n",maxnode,sink,maxabund,*pos);
				CTransfrag *tr=new CTransfrag(nodes,trpat,maxabund);
				transfrag.Add(tr);

				// add sink among maxnode children
				inode=no2gnode[maxnode];
				inode->child.Add(sink);
				inode->childpat[*pos]=1; // source should be already among parents of maxnode but not the edge to source
			}
		}
		*/

		g++;
	}




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

int find_cguidepat(GBitVec& pat,GVec<CTrGuidePat>& patvec) {

	for(int i=0;i<patvec.Count();i++) {
		if(pat==patvec[i].pat) return(i);
	}

	return(-1);
}

/*
int guides_maxflow(int gno,int edgeno,GIntHash<int>& gpos,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,GVec<CGuide>& guidetrf,int& geneno,
		int s,GList<CPrediction>& pred,GVec<float>& nodecov,GBitVec& istranscript,GBitVec& pathpat,bool &first,GPVec<GffObj>& guides,GVec<int> &guidepred, BundleData *bdata) {

	int maxi=1;
	int ng=guidetrf.Count();

	if(ng==1) { // if only one guide I do not need to do the 2 pass
		GVec<float> nodeflux;
		//float fragno=0;
		float flux= max_flow(gno,guidetrf[0].trf->nodes,istranscript,transfrag,no2gnode,nodeflux,guidetrf[0].trf->pattern);
		istranscript.reset();


		{ // DEBUG ONLY
			fprintf(stderr,"guide=%s flux[0]=%g\n",guides[guidetrf[0].g]->getID(),flux);
		}


		if(flux>epsilon) {
			bool include=true;
			store_transcript(pred,guidetrf[0].trf->nodes,nodeflux,nodecov,no2gnode,geneno,first,s,gno,gpos,include,pathpat,bdata,guides[guidetrf[0].g]);
			if(eonly) {
				guidepred[guidetrf[0].g]=pred.Count()-1; // NEED TO TEST: if this doesn't work for single genes I might want to recombine with the previous prediction in store_transcript
				//fprintf(stderr,"guidepred[%d]=%d\n",guidetrf[0].g,guidepred[guidetrf[0].g]);
			}

			nodeflux.Clear();
		}

		// Node coverages:
		for(int i=1;i<gno-1;i++)
			if(nodecov[i]>nodecov[maxi]) maxi=i;

		//if(nodecov[maxi]<readthr) break; // no need to find other paths since they aren't any above allowed read threshold
		//if(nodecov[maxi]<1) break; // I shouldn't be restricting this at all?



		{ // DEBUG ONLY
		  fprintf(stderr,"\nAfter update:\n");
		  for(int i=0;i<gno;i++) {
			  fprintf(stderr,"Node %d: %f ",i,nodecov[i]);
			  fprintf(stderr,"trf=");
			  for(int t=0;t<no2gnode[i]->trf.Count();t++) fprintf(stderr," %d(%f)",no2gnode[i]->trf[t],transfrag[no2gnode[i]->trf[t]]->abundance);
			  fprintf(stderr," maxi=%d maxcov=%f\n",maxi,nodecov[maxi]);
		  }
		}


		//return(maxi);
	}
	else if(ng) {
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

			//fprintf(stderr,"guide=%s ",guides[guidetrf[g].g]->getID());

			float initflux=guideflow(gno,guidetrf[g].trf->nodes,istranscript,transfrag,no2gnode,guidetrf[g].trf->pattern,capacity[g],flow[g],link[g],node2path[g]);
			flux.Add(initflux);
			istranscript.reset();


			{ // DEBUG ONLY
			   fprintf(stderr," flux[%d]=%g\n",g,flux[g]);
			}

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
					store_transcript(pred,guidetrf[g].trf->nodes,nodeflux,nodecov,no2gnode,geneno,first,s,gno,gpos,include,pathpat,bdata,guides[guidetrf[g].g]);
					if(eonly) {
						guidepred[guidetrf[g].g]=pred.Count()-1; // NEED TO TEST: if this doesn't work for single genes I might want to recombine with the previous prediction in store_transcript
						//fprintf(stderr,"2 guidepred[%d]=%d\n",guidetrf[g].g,guidepred[guidetrf[g].g]);
					}
					cov=true;
					// Node coverages:
					for(int i=1;i<gno-1;i++)
						if(nodecov[i]>nodecov[maxi]) maxi=i;

					//if(nodecov[maxi]<readthr) break; // no need to find other paths since they aren't any above allowed read threshold
					//if(nodecov[maxi]<1) break; // I shouldn't be restricting this at all?


					{ // DEBUG ONLY
				  	  fprintf(stderr,"\nAfter update:\n");
				  	  for(int i=0;i<gno;i++) {
					  	  fprintf(stderr,"Node %d: %f ",i,nodecov[i]);
					  	  fprintf(stderr,"trf=");
					  	  for(int t=0;t<no2gnode[i]->trf.Count();t++) fprintf(stderr," %d(%f)",no2gnode[i]->trf[t],transfrag[no2gnode[i]->trf[t]]->abundance);
					  	  fprintf(stderr," maxi=%d maxcov=%f\n",maxi,nodecov[maxi]);
				  	  }
					}


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
	}

	if(eonly && nodecov[maxi]>epsilon) { // this is the end for eonly so I should make sure I use all reads

		for(int g=0;g<guidetrf.Count();g++) delete guidetrf[g].trf;
		guidetrf.Clear();

		char strand='-';
		if(s) strand='+';

		GVec<CNodeGuide> nodeinfo(gno);
		nodeinfo.setCount(gno);

		int ng=0;

		GVec<int> olen;
		// find guides' partial patterns
		for(int g=0;g<guides.Count();g++) {
			//fprintf(stderr,"Consider guide[%d out of %d] %s\n",g,guides.Count(),guides[g]->getID());
			if((guides[g]->strand==strand) && (guides[g]->overlap(no2gnode[1]->start,no2gnode[gno-2]->end))) { // if guide on the same strand and overlaps at all the graph
				//fprintf(stderr,"...there is overlap\n");
				olen.Clear();
				int olensum;
				CTransfrag *trguide=find_guide_partial_pat(guides[g],no2gnode,gno,edgeno,gpos,olen,olensum);
				if(trguide) { // the guide can be found among the graph nodes
					//fprintf(stderr,"...partial pattern found!\n");
					CGuide newguide(trguide,g);
					guidetrf.Add(newguide);
					for(int i=0;i<trguide->nodes.Count();i++) {
						CPartGuide pg(ng,olen[i],olensum,guides[g]->covlen);
						nodeinfo[trguide->nodes[i]].guide.Add(pg);
						int idx=nodeinfo[trguide->nodes[i]].guide.Count()-1; // position of guide in nodeinfo
						bool terminal=true;
						if(i) { // not first node in guide covered
							int *pos=gpos[edge(trguide->nodes[i-1],trguide->nodes[i],gno)];
							if(pos && trguide->pattern[*pos]) terminal=false; // if there is edge from previous node and it's present in guide
						}
						nodeinfo[trguide->nodes[i]].guide[idx].terminal_in=terminal;
						terminal=true;
						if(i<trguide->nodes.Count()-1) { // not last node in guide covered
							int *pos=gpos[edge(trguide->nodes[i],trguide->nodes[i+1],gno)];
							if(pos && trguide->pattern[*pos]) terminal=false; // if there is edge from previous node and it's present in guide
						}
						nodeinfo[trguide->nodes[i]].guide[idx].terminal_out=terminal;
					}
					ng++;
				}
			}
		}

		// future: compute maxflow for each guide starting at the most covered node maybe and going all the way to where I can do it.
		// or maybe compute maxflow for each guide from all terminal nodes and then assign coverages based on winners or the one that gets the most coverage
		// hint: max_flow doesn't need source/sink in the path so I can go with that
		// or: do some type EM algorithm like I am trying to do below

		GVec<int> path(1);
		path.setCount(1);
		GVec<float> nodeflux(1);
		nodeflux.setCount(1);

		GVec<int> coverednode; // array of all nodes that still have coverage _and_ guides

		// determine gcounts and the strict assignments to guides
		for(int n=1;n<gno-1;n++) if(nodecov[n]>epsilon && nodeinfo[n].guide.Count()) { // for all nodes for which there is still coverage and have guides
			ng=nodeinfo[n].guide.Count(); // number of guides in node
			if(ng==1) { // there is only one guide -> it gets all the coverage
				int g=guidetrf[nodeinfo[n].guide[0].idx].g;
				path[0]=n;
				nodeflux[0]=1;
				//fprintf(stderr,"1 g=%d\n",g);
				if(guidepred[g]!=-1) update_guide_pred(pred,guidepred[g],path,nodeflux,nodecov,no2gnode,gno,true);
				else // new prediction for this guide
					guidepred[g]=store_guide_transcript(pred,path,nodeflux,nodecov,no2gnode,geneno,first,gno,guides[g],true);
			}
			else { // there are more than one guide overlapping the node -> count the strict/loose counts; might want to rethink doing this for single nodes in the graph

				float incount=0; // should be distrtibuted between incoming guides that continue and terminal ones -> not doing this now
				float outcount=0;
				//float throughcount=0;
				float terminal_incount=0;
				float terminal_outcount=0;
				float guide_incount=0;
				float guide_outcount=0;

				GVec<float> in_strictcount(ng,float(0)); // these are abundances that are unique to a single guide -> this will get assigned to guides
				GVec<float> in_loosecount(ng,float(0)); // these are guide abundances that are shared with other guides
				GVec<float> in_looseterminalcount(ng,float(0)); // these are guide abundances that are shared with other guides
				GVec<float> out_strictcount(ng,float(0));
				GVec<float> out_loosecount(ng,float(0));
				GVec<float> out_looseterminalcount(ng,float(0));

				GVec<CTrGuidePat> out_trcount; // I can use the node trcount for the in_counts but I need for the outcount because it needs to be adjusted by the rate

				int ncovguides=0;
				bool allstrictzero=true;
				bool allloosezero=true;

				CGraphnode *inode=no2gnode[n];
				int nn=inode->trf.Count();
				for(int j=0;j<nn;j++){
					int t=inode->trf[j];
					if(transfrag[t]->abundance>epsilon) { // only if transfrag still has abundance it's worth considering
						GVec<int> compguide; // keeps all guides that are compatible with this transfrag
						GBitVec guidepat(guidetrf.Count()); // pattern of guides that are present in transfrag
						bool terminal_in=true;  // true only if all guides compatible with transfrag are terminal
						bool terminal_out=true;

						for(int i=0;i<ng;i++) {
							int g=nodeinfo[n].guide[i].idx;
							if(((transfrag[t]->pattern) & guidetrf[g].trf->pattern) == transfrag[t]->pattern) { // transfrag is compatible to guide
								compguide.Add(i); // make sure that later you change it to guidetrf indexes
								guidepat[g]=1;
								if(terminal_in && !nodeinfo[n].guide[i].terminal_in) terminal_in=false;
								if(terminal_out && !nodeinfo[n].guide[i].terminal_out) terminal_out=false;
								if(!nodeinfo[n].guide[i].gcount) {
									nodeinfo[n].guide[i].gcount=1;
									ncovguides++;
								}
							}
						}
						if(!compguide.Count()) {
							terminal_in=false;
							terminal_out=false;
						}

						if(transfrag[t]->nodes.Last()==n) { // transfrag ends at this node (in transfrag)
							if(terminal_in)  {
								if(compguide.Count()>1) { // all guides compatible with transfrag are terminal
									terminal_incount+=transfrag[t]->abundance;
									for(int i=0;i<compguide.Count();i++) {
										in_looseterminalcount[compguide[i]]+=transfrag[t]->abundance;
										allloosezero=false;
									}
								}
							}
							else incount+=transfrag[t]->abundance; // this biases in favor of continuing guides -> might want to rethink this

							if(compguide.Count()==1) { // only one guide compatible with transfrag
								in_strictcount[compguide[0]]+=transfrag[t]->abundance;
								allstrictzero=false;
							}
							else if(compguide.Count()) { // there are multiple guides in transfrag
								int p=find_cguidepat(guidepat,nodeinfo[n].trcount);
								if(p==-1) { // didn't find the guide pattern
									CTrGuidePat pat(guidepat,transfrag[t]->abundance,terminal_in);
									nodeinfo[n].trcount.Add(pat);
									for(int i=0;i<compguide.Count();i++) nodeinfo[n].trcount[nodeinfo[n].trcount.Count()-1].g.Add(nodeinfo[n].guide[compguide[i]].idx);
								}
								else nodeinfo[n].trcount[p].abund+=transfrag[t]->abundance;

								if(!terminal_in) {
									guide_incount+=transfrag[t]->abundance; // terminals can only be guide transfrags
									for(int i=0;i<compguide.Count();i++) {
										in_loosecount[compguide[i]]+=transfrag[t]->abundance;
										allloosezero=false;
									}
								}
							}
						}
						else if(transfrag[t]->nodes[0]==n) { // transfrag starts at this node (out transfrag)
							if(terminal_out) {
								if(compguide.Count()>1) {
									terminal_outcount+=transfrag[t]->abundance;
									for(int i=0;i<compguide.Count();i++) {
										out_looseterminalcount[compguide[i]]+=transfrag[t]->abundance;
										allloosezero=false;
									}
								}
							}
							else outcount+=transfrag[t]->abundance;
							if(compguide.Count()==1) {
								out_strictcount[compguide[0]]+=transfrag[t]->abundance; // only one guide compatible with transfrag
								allstrictzero=false;
							}
							else if(compguide.Count()){
								int p=find_cguidepat(guidepat,out_trcount);
								if(p==-1) { // didn't find the guide pattern
									CTrGuidePat pat(guidepat,transfrag[t]->abundance,terminal_out);
									out_trcount.Add(pat);
									for(int i=0;i<compguide.Count();i++) out_trcount[out_trcount.Count()-1].g.Add(nodeinfo[n].guide[compguide[i]].idx);
								}
								else out_trcount[p].abund+=transfrag[t]->abundance;

								if(!terminal_out) {
									guide_outcount+=transfrag[t]->abundance;

									for(int i=0;i<compguide.Count();i++) {
										out_loosecount[compguide[i]]+=transfrag[t]->abundance;
										allloosezero=false;
									}
								}
							}
						}
						else if(transfrag[t]->pattern[n]) { // through transfrag (here I checked that the transfrag clearly goes through the node)
							//throughcount+=transfrag[t]->abundance;
							incount+=transfrag[t]->abundance; // I do this instead of the above because the throughout could be very high and the in/outcounts shouldn't dominate in that case
							outcount+=transfrag[t]->abundance;

							if(compguide.Count()==1) {
								in_strictcount[compguide[0]]+=transfrag[t]->abundance; // only one guide compatible with transfrag
								out_strictcount[compguide[0]]+=transfrag[t]->abundance; // only one guide compatible with transfrag
								allstrictzero=false;
							}
							else if(compguide.Count()) {
								guide_incount+=transfrag[t]->abundance;
								guide_outcount+=transfrag[t]->abundance;
								int p=find_cguidepat(guidepat,nodeinfo[n].trcount);
								if(p==-1) { // didn't find the guide pattern
									CTrGuidePat pat(guidepat,transfrag[t]->abundance,terminal_in);
									nodeinfo[n].trcount.Add(pat);
									for(int i=0;i<compguide.Count();i++) nodeinfo[n].trcount[nodeinfo[n].trcount.Count()-1].g.Add(nodeinfo[n].guide[compguide[i]].idx);
								}
								else {
									nodeinfo[n].trcount[p].abund+=transfrag[t]->abundance;
								}

								p=find_cguidepat(guidepat,out_trcount);
								if(p==-1) { // didn't find the guide pattern
									CTrGuidePat pat(guidepat,transfrag[t]->abundance,terminal_out);
									out_trcount.Add(pat);
									for(int i=0;i<compguide.Count();i++) out_trcount[out_trcount.Count()-1].g.Add(nodeinfo[n].guide[compguide[i]].idx);
								}
								else {
									out_trcount[p].abund+=transfrag[t]->abundance;
								}

								for(int i=0;i<compguide.Count();i++) {
									in_loosecount[compguide[i]]+=transfrag[t]->abundance;
									out_loosecount[compguide[i]]+=transfrag[t]->abundance;
									allloosezero=false;
								}
							}
						}
					}
				} // done computing counts for node

				if(allstrictzero && allloosezero) { // only single guides or no transfrag compatible to guides
												// -> I need to choose some preferences: bias toward longest overlapping guide
					for(int i=0;i<ng;i++) { // ng>1
						nodeinfo[n].guide[i].gcount=1; // give them a uniform probability at node level: coverage probabilities take priority
					}
					nodeinfo[n].sumtrcount=1; // one guide is sufficient to explain all abundances;
					// Note: pattern of guides in trcount should all be with only one 1 for the actual guide but we don't need to use it since using the guide explains everything
				}
				else { // some guides are compatible to the node transfrags

					if(ncovguides==1) { // only one guide is compatible with the node transfrags -> takes all node coverage (this biases strongly against single exon guides, or incomplete guides -> might want to rethink)
						int g=guidetrf[nodeinfo[n].guide[0].idx].g;
						path[0]=n;
						nodeflux[0]=1;
						//fprintf(stderr,"2 g=%d\n",g);
						if(guidepred[g]!=-1) update_guide_pred(pred,guidepred[g],path,nodeflux,nodecov,no2gnode,gno,true);
						else // new prediction for this guide
							guidepred[g]=store_guide_transcript(pred,path,nodeflux,nodecov,no2gnode,geneno,first,gno,guides[g],true);
					}
					else { // more than one guide got counts

						// first get the rate
						float rate=1;
						nodeinfo[n].sumtrcount=0;
						if(incount && outcount) {
							rate=incount/outcount;
							nodeinfo[n].sumtrcount+=rate*terminal_outcount;
							nodeinfo[n].sumtrcount+=terminal_incount;
						}
						nodeinfo[n].sumtrcount+=guide_incount+rate*guide_outcount+terminal_incount+rate*terminal_outcount;

						float sumstrict=nodeinfo[n].sumtrcount;
						if(!allstrictzero) { // there are strict counts that need to be assigned to
							for(int i=0;i<ng;i++) {
								sumstrict+=in_strictcount[i]+rate*out_strictcount[i];
								if(nodeinfo[n].guide[i].terminal_in) sumstrict+=in_strictcount[i]; // if node is terminal the counts should be adjusted
								if(nodeinfo[n].guide[i].terminal_out) sumstrict+=rate*out_strictcount[i];
							}
						}

						// update guide coverages with unique reads and gcounts for later
						for(int i=0;i<ng;i++) {

							if(!allstrictzero) {
								float strictcount=in_strictcount[i]+rate*out_strictcount[i];
								if(incount && outcount) {
									if(nodeinfo[n].guide[i].terminal_in) strictcount+=in_strictcount[i];
									if(nodeinfo[n].guide[i].terminal_out) strictcount+=rate*out_strictcount[i];
								}
								if(strictcount) { // if there are reads unique to the guide
									int g=guidetrf[nodeinfo[n].guide[i].idx].g;
									path[0]=n;
									nodeflux[0]=strictcount/sumstrict;
									//fprintf(stderr,"3 g=%d\n",g);
									if(guidepred[g]!=-1) update_guide_pred(pred,guidepred[g],path,nodeflux,nodecov,no2gnode,gno,false);
									else // new prediction for this guide
										guidepred[g]=store_guide_transcript(pred,path,nodeflux,nodecov,no2gnode,geneno,first,gno,guides[g],false);
								}
							}

							// now set gcounts
							nodeinfo[n].guide[i].gcount=in_loosecount[i]+rate*out_loosecount[i]+
									in_looseterminalcount[i]+rate*out_looseterminalcount[i];
							if(incount && outcount) {
								if(nodeinfo[n].guide[i].terminal_in) nodeinfo[n].guide[i].gcount+=in_looseterminalcount[i];
								if(nodeinfo[n].guide[i].terminal_out) nodeinfo[n].guide[i].gcount+=rate*out_looseterminalcount[i];
							}
						}

						if(!allstrictzero) { // update node coverage here
							if(allloosezero) nodecov[n]=0;
							else nodecov[n]*=nodeinfo[n].sumtrcount/sumstrict;
						}

						if(!allloosezero && nodecov[n]) { // get the pattern counts
							if(incount && outcount) for(int i=0;i<nodeinfo[n].trcount.Count();i++)
								if(nodeinfo[n].trcount[i].terminal) nodeinfo[n].trcount[i].abund*=2;
							for(int i=0;i<out_trcount.Count();i++) {
								int t=find_cguidepat(out_trcount[i].pat,nodeinfo[n].trcount);
								out_trcount[i].abund*=rate;
								if(out_trcount[i].terminal && incount && outcount) out_trcount[i].abund*=2;
								//else out_trcount[i].terminal=false; // I can not have a trcount with the same guides that is terminal both for in and out
								if(t==-1) { // didn't find the guide pattern
									nodeinfo[n].trcount.Add(out_trcount[i]);
								}
								else nodeinfo[n].trcount[t].abund+=out_trcount[i].abund;
							}
						}

					}
				}
			}
			if(nodecov[n]>epsilon) { // there is still coverage left for the node
				coverednode.Add(n);
			}
		}

		if(coverednode.Count()) { // I still have covered nodes: asign reads to guides; I need to do some sort of EM algorithm here
			// my prior bias in assigning reads should be based on coverages allready found

			int nEM=10; // number of times to repeat the EM algorithm: need to see if this is a good count, or if I should use and epsilon change to stop
			ng=guidetrf.Count();
			GVec<float> initgcov(ng,float(0));

			// first set the initial coverages in covered nodes
			for(int i=0;i<ng;i++) {
				int np=guidepred[guidetrf[i].g];
				if(np!=-1) {
					initgcov[i]=pred[np]->cov*pred[np]->tlen;
				}
			}

			GVec<float> prevgcov(ng,float(0)); // the new gcov that will be computed as part of the EM algorithm

			// now set the cov's in CPartGuide -> for sort purposes and for the EM algorithm
			for(int i=0;i<coverednode.Count();i++) {
				int n=coverednode[i];
				for(int j=0;j<nodeinfo[n].guide.Count();j++) {
					int g=nodeinfo[n].guide[j].idx;
					nodeinfo[n].guide[j].cov=initgcov[g];
					prevgcov[g]=initgcov[g];
				}
			}

			GVec<float> gcov(ng,float(0)); // the new gcov that will be computed as part of the EM algorithm

			int m=0;
			while(m<nEM) { // do the EM

				for(int i=0;i<ng;i++) gcov[i]=initgcov[i];

				for(int i=0;i<coverednode.Count();i++) { // for each node recompute probabilities
					int n=coverednode[i];
					nodeinfo[n].guide.Sort(partguideCmp); // I need to sort at each step because the cov's get updated
					float totalcount=0;
					bool covnotreached=true;
					float sumgcov=0;
					for(int j=0;j<nodeinfo[n].guide.Count();j++) { // for each guide in the node, starting from the most abundant to the least
						nodeinfo[n].guide[j].ncov=0; // I need to set this for each guide -> this is why I not break from the loop below
						int g=nodeinfo[n].guide[j].idx; // index in guidetrf
						sumgcov+=prevgcov[g];
						if(covnotreached) for(int t=0;t<nodeinfo[n].trcount.Count();t++) { // for each trguidepat in the node
							if(nodeinfo[n].trcount[t].pat[g]) { // current guide is in this pattern
								float sum=0; // sum of all coverages of guides in the transcript
								for(int k=0;k<nodeinfo[n].trcount[t].g.Count();k++) {
									int l=nodeinfo[n].trcount[t].g[k];
									sum+=prevgcov[l];
								}
								float abund=nodeinfo[n].sumtrcount-totalcount; // this is the maximum abundance allowed for this guide
								if(sum) { // I have some guides that have coverage
									float newabund=nodeinfo[n].trcount[t].abund*prevgcov[g]/sum;
									if(newabund<abund) abund=newabund;
								}
								else { // no guide has assigned coverages yet -> winner takes all
									if(nodeinfo[n].trcount[t].abund<abund) abund=nodeinfo[n].trcount[t].abund;
								}
								nodeinfo[n].guide[j].ncov+=abund;
								gcov[g]+=abund*nodeinfo[n].guide[j].olen/no2gnode[n]->len();
								totalcount+=abund;
								if(nodeinfo[n].sumtrcount-totalcount<epsilon) {
									covnotreached=false;
									break;
								}
							}
						}
					}
					if(covnotreached) { // there were not enough transfrags to give coverages -> redistribute reads based on coverages
						float abund=nodeinfo[n].sumtrcount-totalcount; // this is how much is left
						if(sumgcov) for(int j=0;j<nodeinfo[n].guide.Count();j++) {
							if(!nodeinfo[n].guide[j].cov) break;
							float newabund=abund*nodeinfo[n].guide[j].cov/sumgcov;
							nodeinfo[n].guide[j].ncov+=newabund;
							int g=nodeinfo[n].guide[j].idx;
							gcov[g]+=newabund*nodeinfo[n].guide[j].olen/no2gnode[n]->len();
						}
						else { // all guides have cov zero -> winner takes all
							nodeinfo[n].guide[0].ncov+=abund;
							int g=nodeinfo[n].guide[0].idx;
							gcov[g]+=abund*nodeinfo[n].guide[0].olen/no2gnode[n]->len();
						}
					}
				} // end node

				bool nochange=true;
				// check coverages
				for(int i=0;i<coverednode.Count();i++) {
					int n=coverednode[i];
					for(int j=0;j<nodeinfo[n].guide.Count();j++) {
						int g=nodeinfo[n].guide[j].idx;
						float diff=nodeinfo[n].guide[j].cov-gcov[g];
						if(diff<0) diff=-diff;
						if(diff) {
							prevgcov[g]=gcov[g];
							nodeinfo[n].guide[j].cov=gcov[g];
							if(diff>1) nochange=false;
						}
					}
				}

				// see if there is a change in probability
				if(nochange) break;

				m++;
			} // end EM

			// new coverages are estimated so now I can assign the reads to the predictions
			for(int g=0;g<guidetrf.Count();g++) if(gcov[g]>initgcov[g]) { // for each guide that has more coverage to add
				path.Clear();
				nodeflux.Clear();
				for(int i=0;i<guidetrf[g].trf->nodes.Count();i++) if(nodecov[guidetrf[g].trf->nodes[i]]) { // for each covered node
					int n=guidetrf[g].trf->nodes[i];
					for(int j=0;j<nodeinfo[n].guide.Count();j++) if(g==nodeinfo[n].guide[j].idx){ // found my guide
						if(nodeinfo[n].guide[j].ncov) { // if there is coverage assigned to node -> i should add it to the path
							path.Add(n);
							float prop=nodeinfo[n].guide[j].ncov/nodeinfo[n].sumtrcount;
							nodeflux.Add(prop);
						}
						break;
					}
				}
				//fprintf(stderr,"4 g=%d\n",guidetrf[g].g);
				if(guidepred[guidetrf[g].g]!=-1) update_guide_pred(pred,guidepred[guidetrf[g].g],path,nodeflux,nodecov,no2gnode,gno,false);
				else // new prediction for this guide
					guidepred[guidetrf[g].g]=store_guide_transcript(pred,path,nodeflux,nodecov,no2gnode,geneno,first,gno,guides[guidetrf[g].g],false);
			}
		}

		maxi=0;
	}

	return(maxi);
}
*/

int guides_pushmaxflow(int gno,int edgeno,GIntHash<int>& gpos,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,GVec<CGuide>& guidetrf,int& geneno,
		int s,GList<CPrediction>& pred,GVec<float>& nodecov,GBitVec& istranscript,GBitVec& pathpat,bool &first,GPVec<GffObj>& guides,GVec<int> &guidepred, BundleData *bdata) {

	int maxi=1;
	int ng=guidetrf.Count();

	if(ng==1) { // if only one guide I do not need to do the 2 pass
		GVec<float> nodeflux;
		//float fragno=0;
		float flux= push_max_flow(gno,guidetrf[0].trf->nodes,istranscript,transfrag,no2gnode,nodeflux,guidetrf[0].trf->pattern,gpos);
		istranscript.reset();

		/*
		{ // DEBUG ONLY
			fprintf(stderr,"guide=%s flux[0]=%g\n",guides[guidetrf[0].g]->getID(),flux);
		}
		*/

		if(flux>epsilon) {
			bool include=true;
			if(guidepred[guidetrf[0].g]==-1) {

				store_transcript(pred,guidetrf[0].trf->nodes,nodeflux,nodecov,no2gnode,geneno,first,s,gno,gpos,include,pathpat,bdata,guides[guidetrf[0].g]);
				//if(eonly) { // this is not correct because it might have been assigned before
					guidepred[guidetrf[0].g]=pred.Count()-1; // NEED TO TEST: if this doesn't work for single genes I might want to recombine with the previous prediction in store_transcript
					//fprintf(stderr,"guidepred[%d]=%d\n",guidetrf[0].g,guidepred[guidetrf[0].g]);
				//}
			}
			else {
				update_guide_pred(pred,guidepred[guidetrf[0].g],guidetrf[0].trf->nodes,nodeflux,nodecov,no2gnode,gno,true);
			}

		}

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

		//return(maxi);
	}
	else if(ng) {


		// this is an abundance for guides based on maximum flow for each guide (how much can they each carry
		for(int g=0;g<guidetrf.Count();g++) {
			guidetrf[g].trf->abundance=push_guide_maxflow(gno,guidetrf[g].trf->nodes,istranscript,transfrag,no2gnode,guidetrf[g].trf->pattern);
			istranscript.reset();
		}
		guidetrf.Sort(guidedabundCmp);

		/*
		{ // DEBUG ONLY
			for(int g=0;g<guidetrf.Count();g++) {
				fprintf(stderr,"Abundance of guide[%d]=%f with nodes:",g,guidetrf[g].trf->abundance);
				for(int i=0;i<guidetrf[g].trf->nodes.Count();i++) fprintf(stderr," %d",guidetrf[g].trf->nodes[i]);
				fprintf(stderr,"\n");
			}
		}
		*/

		GVec<float> nodeflux;
		for(int g=ng-1;g>=0;g--) if(guidetrf[g].trf->abundance){ // calculate maximum push flow for each guide starting from the less covered one

			float flux=guidepushflow(g,guidetrf,gno,istranscript,transfrag,no2gnode,nodeflux);

			istranscript.reset();

			/*
			{ // DEBUG ONLY
				fprintf(stderr,"guide=%s flux[%d]=%f\n",guides[g]->getID(),g,flux);
			}
			*/

			bool include=true;
			if(flux>epsilon) {
				if(guidepred[guidetrf[g].g]==-1) {
					store_transcript(pred,guidetrf[g].trf->nodes,nodeflux,nodecov,no2gnode,geneno,first,s,gno,gpos,include,pathpat,bdata,guides[guidetrf[g].g]);
					//if(eonly) {
						guidepred[guidetrf[g].g]=pred.Count()-1; // NEED TO TEST: if this doesn't work for single genes I might want to recombine with the previous prediction in store_transcript
						//fprintf(stderr,"2 guidepred[%d]=%d\n",guidetrf[g].g,guidepred[guidetrf[g].g]);
						//}

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
				}
				else {
					update_guide_pred(pred,guidepred[guidetrf[g].g],guidetrf[g].trf->nodes,nodeflux,nodecov,no2gnode,gno,true);
				}

				nodeflux.Clear();

			}
			else { // it's possible that this is a single exon gene included in a much larger interval and this is why it didn't get predicted
				if(guides[guidetrf[g].g]->exons.Count()==1) { // single exon gene included in another prediction
					// check if it overlaps other single exon genes

					bool overlap=false;
					for(int r=ng-1;r>g;r--) if(guidepred[guidetrf[r].g]>-1){
						if((guidetrf[g].trf->pattern & guidetrf[r].trf->pattern)==guidetrf[g].trf->pattern) {
							if(guides[guidetrf[r].g]->exons.Count()>1) {
								overlap=true;
								break;
							}
							else { // overlaps single exon gene -> check if it's a true overlap
								if(guides[guidetrf[r].g]->exons[0]->overlap(guides[guidetrf[g].g]->exons[0])) {
									overlap=true;
									break;
								}
							}
						}
					}
					if(!overlap) { // no overlap detected
						for(int z=0;z<guidetrf[g].trf->nodes.Count();z++) nodeflux.cAdd(1.0);
						if(guidepred[guidetrf[g].g]==-1) {
							store_transcript(pred,guidetrf[g].trf->nodes,nodeflux,nodecov,no2gnode,geneno,first,s,gno,gpos,include,pathpat,bdata,guides[guidetrf[g].g]);
							//if(eonly) {
								guidepred[guidetrf[g].g]=pred.Count()-1; // NEED TO TEST: if this doesn't work for single genes I might want to recombine with the previous prediction in store_transcript
								//fprintf(stderr,"2 guidepred[%d]=%d\n",guidetrf[g].g,guidepred[guidetrf[g].g]);
							//}
						}
						else {
							update_guide_pred(pred,guidepred[guidetrf[g].g],guidetrf[0].trf->nodes,nodeflux,nodecov,no2gnode,gno,true);
						}
						nodeflux.Clear();
					}
				}
			}
		}

	}

	// Node coverages:
	for(int i=1;i<gno-1;i++)
		if(nodecov[i]>nodecov[maxi]) maxi=i;


	if(eonly && nodecov[maxi]>epsilon) { // this is the end for eonly so I should make sure I use all reads

		for(int g=0;g<guidetrf.Count();g++) delete guidetrf[g].trf;
		guidetrf.Clear();

		char strand='-';
		if(s) strand='+';

		GVec<CNodeGuide> nodeinfo(gno);
		nodeinfo.setCount(gno);

		int ng=0;

		GVec<int> olen;
		// find guides' partial patterns
		for(int g=0;g<guides.Count();g++) {
			//fprintf(stderr,"Consider guide[%d out of %d] %s\n",g,guides.Count(),guides[g]->getID());
			if((guides[g]->strand==strand) && (guides[g]->overlap(no2gnode[1]->start,no2gnode[gno-2]->end))) { // if guide on the same strand and overlaps at all the graph
				//fprintf(stderr,"...there is overlap\n");
				olen.Clear();
				int olensum;
				CTransfrag *trguide=find_guide_partial_pat(guides[g],no2gnode,gno,edgeno,gpos,olen,olensum);
				if(trguide) { // the guide can be found among the graph nodes
					//fprintf(stderr,"...partial pattern found!\n");
					CGuide newguide(trguide,g);
					guidetrf.Add(newguide);
					for(int i=0;i<trguide->nodes.Count();i++) {
						CPartGuide pg(ng,olen[i],olensum,guides[g]->covlen);
						nodeinfo[trguide->nodes[i]].guide.Add(pg);
						int idx=nodeinfo[trguide->nodes[i]].guide.Count()-1; // position of guide in nodeinfo
						bool terminal=true;
						if(i) { // not first node in guide covered
							int *pos=gpos[edge(trguide->nodes[i-1],trguide->nodes[i],gno)];
							if(pos && trguide->pattern[*pos]) terminal=false; // if there is edge from previous node and it's present in guide
						}
						nodeinfo[trguide->nodes[i]].guide[idx].terminal_in=terminal;
						terminal=true;
						if(i<trguide->nodes.Count()-1) { // not last node in guide covered
							int *pos=gpos[edge(trguide->nodes[i],trguide->nodes[i+1],gno)];
							if(pos && trguide->pattern[*pos]) terminal=false; // if there is edge from previous node and it's present in guide
						}
						nodeinfo[trguide->nodes[i]].guide[idx].terminal_out=terminal;
					}
					ng++;
				}
			}
		}

		// future: compute maxflow for each guide starting at the most covered node maybe and going all the way to where I can do it.
		// or maybe compute maxflow for each guide from all terminal nodes and then assign coverages based on winners or the one that gets the most coverage
		// hint: max_flow doesn't need source/sink in the path so I can go with that
		// or: do some type EM algorithm like I am trying to do below

		GVec<int> path(1);
		path.setCount(1);
		GVec<float> nodeflux(1);
		nodeflux.setCount(1);

		GVec<int> coverednode; // array of all nodes that still have coverage _and_ guides

		// determine gcounts and the strict assignments to guides
		for(int n=1;n<gno-1;n++) if(nodecov[n]>epsilon && nodeinfo[n].guide.Count()) { // for all nodes for which there is still coverage and have guides
			ng=nodeinfo[n].guide.Count(); // number of guides in node
			if(ng==1) { // there is only one guide -> it gets all the coverage
				int g=guidetrf[nodeinfo[n].guide[0].idx].g;
				path[0]=n;
				nodeflux[0]=1;
				//fprintf(stderr,"1 g=%d\n",g);
				if(guidepred[g]!=-1) update_guide_pred(pred,guidepred[g],path,nodeflux,nodecov,no2gnode,gno,true);
				else // new prediction for this guide
					guidepred[g]=store_guide_transcript(pred,path,nodeflux,nodecov,no2gnode,geneno,first,gno,guides[g],true);
			}
			else { // there are more than one guide overlapping the node -> count the strict/loose counts; might want to rethink doing this for single nodes in the graph

				float incount=0; // should be distrtibuted between incoming guides that continue and terminal ones -> not doing this now
				float outcount=0;
				//float throughcount=0;
				float terminal_incount=0;
				float terminal_outcount=0;
				float guide_incount=0;
				float guide_outcount=0;

				GVec<float> in_strictcount(ng,float(0)); // these are abundances that are unique to a single guide -> this will get assigned to guides
				GVec<float> in_loosecount(ng,float(0)); // these are guide abundances that are shared with other guides
				GVec<float> in_looseterminalcount(ng,float(0)); // these are guide abundances that are shared with other guides
				GVec<float> out_strictcount(ng,float(0));
				GVec<float> out_loosecount(ng,float(0));
				GVec<float> out_looseterminalcount(ng,float(0));

				GVec<CTrGuidePat> out_trcount; // I can use the node trcount for the in_counts but I need for the outcount because it needs to be adjusted by the rate

				int ncovguides=0;
				bool allstrictzero=true;
				bool allloosezero=true;

				CGraphnode *inode=no2gnode[n];
				int nn=inode->trf.Count();
				for(int j=0;j<nn;j++){
					int t=inode->trf[j];
					if(transfrag[t]->abundance>epsilon) { // only if transfrag still has abundance it's worth considering
						GVec<int> compguide; // keeps all guides that are compatible with this transfrag
						GBitVec guidepat(guidetrf.Count()); // pattern of guides that are present in transfrag
						bool terminal_in=true;  // true only if all guides compatible with transfrag are terminal
						bool terminal_out=true;

						for(int i=0;i<ng;i++) {
							int g=nodeinfo[n].guide[i].idx;
							if(((transfrag[t]->pattern) & guidetrf[g].trf->pattern) == transfrag[t]->pattern) { // transfrag is compatible to guide
								compguide.Add(i); // make sure that later you change it to guidetrf indexes
								guidepat[g]=1;
								if(terminal_in && !nodeinfo[n].guide[i].terminal_in) terminal_in=false;
								if(terminal_out && !nodeinfo[n].guide[i].terminal_out) terminal_out=false;
								if(!nodeinfo[n].guide[i].gcount) {
									nodeinfo[n].guide[i].gcount=1;
									ncovguides++;
								}
							}
						}
						if(!compguide.Count()) {
							terminal_in=false;
							terminal_out=false;
						}

						if(transfrag[t]->nodes.Last()==n) { // transfrag ends at this node (in transfrag)
							if(terminal_in)  {
								if(compguide.Count()>1) { // all guides compatible with transfrag are terminal
									terminal_incount+=transfrag[t]->abundance;
									for(int i=0;i<compguide.Count();i++) {
										in_looseterminalcount[compguide[i]]+=transfrag[t]->abundance;
										allloosezero=false;
									}
								}
							}
							else incount+=transfrag[t]->abundance; // this biases in favor of continuing guides -> might want to rethink this

							if(compguide.Count()==1) { // only one guide compatible with transfrag
								in_strictcount[compguide[0]]+=transfrag[t]->abundance;
								allstrictzero=false;
							}
							else if(compguide.Count()) { // there are multiple guides in transfrag
								int p=find_cguidepat(guidepat,nodeinfo[n].trcount);
								if(p==-1) { // didn't find the guide pattern
									CTrGuidePat pat(guidepat,transfrag[t]->abundance,terminal_in);
									nodeinfo[n].trcount.Add(pat);
									for(int i=0;i<compguide.Count();i++) nodeinfo[n].trcount[nodeinfo[n].trcount.Count()-1].g.Add(nodeinfo[n].guide[compguide[i]].idx);
								}
								else nodeinfo[n].trcount[p].abund+=transfrag[t]->abundance;

								if(!terminal_in) {
									guide_incount+=transfrag[t]->abundance; // terminals can only be guide transfrags
									for(int i=0;i<compguide.Count();i++) {
										in_loosecount[compguide[i]]+=transfrag[t]->abundance;
										allloosezero=false;
									}
								}
							}
						}
						else if(transfrag[t]->nodes[0]==n) { // transfrag starts at this node (out transfrag)
							if(terminal_out) {
								if(compguide.Count()>1) {
									terminal_outcount+=transfrag[t]->abundance;
									for(int i=0;i<compguide.Count();i++) {
										out_looseterminalcount[compguide[i]]+=transfrag[t]->abundance;
										allloosezero=false;
									}
								}
							}
							else outcount+=transfrag[t]->abundance;
							if(compguide.Count()==1) {
								out_strictcount[compguide[0]]+=transfrag[t]->abundance; // only one guide compatible with transfrag
								allstrictzero=false;
							}
							else if(compguide.Count()){
								int p=find_cguidepat(guidepat,out_trcount);
								if(p==-1) { // didn't find the guide pattern
									CTrGuidePat pat(guidepat,transfrag[t]->abundance,terminal_out);
									out_trcount.Add(pat);
									for(int i=0;i<compguide.Count();i++) out_trcount[out_trcount.Count()-1].g.Add(nodeinfo[n].guide[compguide[i]].idx);
								}
								else out_trcount[p].abund+=transfrag[t]->abundance;

								if(!terminal_out) {
									guide_outcount+=transfrag[t]->abundance;

									for(int i=0;i<compguide.Count();i++) {
										out_loosecount[compguide[i]]+=transfrag[t]->abundance;
										allloosezero=false;
									}
								}
							}
						}
						else if(transfrag[t]->pattern[n]) { // through transfrag (here I checked that the transfrag clearly goes through the node)
							//throughcount+=transfrag[t]->abundance;
							incount+=transfrag[t]->abundance; // I do this instead of the above because the throughout could be very high and the in/outcounts shouldn't dominate in that case
							outcount+=transfrag[t]->abundance;

							if(compguide.Count()==1) {
								in_strictcount[compguide[0]]+=transfrag[t]->abundance; // only one guide compatible with transfrag
								out_strictcount[compguide[0]]+=transfrag[t]->abundance; // only one guide compatible with transfrag
								allstrictzero=false;
							}
							else if(compguide.Count()) {
								guide_incount+=transfrag[t]->abundance;
								guide_outcount+=transfrag[t]->abundance;
								int p=find_cguidepat(guidepat,nodeinfo[n].trcount);
								if(p==-1) { // didn't find the guide pattern
									CTrGuidePat pat(guidepat,transfrag[t]->abundance,terminal_in);
									nodeinfo[n].trcount.Add(pat);
									for(int i=0;i<compguide.Count();i++) nodeinfo[n].trcount[nodeinfo[n].trcount.Count()-1].g.Add(nodeinfo[n].guide[compguide[i]].idx);
								}
								else {
									nodeinfo[n].trcount[p].abund+=transfrag[t]->abundance;
								}

								p=find_cguidepat(guidepat,out_trcount);
								if(p==-1) { // didn't find the guide pattern
									CTrGuidePat pat(guidepat,transfrag[t]->abundance,terminal_out);
									out_trcount.Add(pat);
									for(int i=0;i<compguide.Count();i++) out_trcount[out_trcount.Count()-1].g.Add(nodeinfo[n].guide[compguide[i]].idx);
								}
								else {
									out_trcount[p].abund+=transfrag[t]->abundance;
								}

								for(int i=0;i<compguide.Count();i++) {
									in_loosecount[compguide[i]]+=transfrag[t]->abundance;
									out_loosecount[compguide[i]]+=transfrag[t]->abundance;
									allloosezero=false;
								}
							}
						}
					}
				} // done computing counts for node

				if(allstrictzero && allloosezero) { // only single guides or no transfrag compatible to guides
												// -> I need to choose some preferences: bias toward longest overlapping guide
					for(int i=0;i<ng;i++) { // ng>1
						nodeinfo[n].guide[i].gcount=1; // give them a uniform probability at node level: coverage probabilities take priority
					}
					nodeinfo[n].sumtrcount=1; // one guide is sufficient to explain all abundances;
					// Note: pattern of guides in trcount should all be with only one 1 for the actual guide but we don't need to use it since using the guide explains everything
				}
				else { // some guides are compatible to the node transfrags

					if(ncovguides==1) { // only one guide is compatible with the node transfrags -> takes all node coverage (this biases strongly against single exon guides, or incomplete guides -> might want to rethink)
						int g=guidetrf[nodeinfo[n].guide[0].idx].g;
						path[0]=n;
						nodeflux[0]=1;
						//fprintf(stderr,"2 g=%d\n",g);
						if(guidepred[g]!=-1) update_guide_pred(pred,guidepred[g],path,nodeflux,nodecov,no2gnode,gno,true);
						else // new prediction for this guide
							guidepred[g]=store_guide_transcript(pred,path,nodeflux,nodecov,no2gnode,geneno,first,gno,guides[g],true);
					}
					else { // more than one guide got counts

						// first get the rate
						float rate=1;
						nodeinfo[n].sumtrcount=0;
						if(incount && outcount) {
							rate=incount/outcount;
							nodeinfo[n].sumtrcount+=rate*terminal_outcount;
							nodeinfo[n].sumtrcount+=terminal_incount;
						}
						nodeinfo[n].sumtrcount+=guide_incount+rate*guide_outcount+terminal_incount+rate*terminal_outcount;

						float sumstrict=nodeinfo[n].sumtrcount;
						if(!allstrictzero) { // there are strict counts that need to be assigned to
							for(int i=0;i<ng;i++) {
								sumstrict+=in_strictcount[i]+rate*out_strictcount[i];
								if(nodeinfo[n].guide[i].terminal_in) sumstrict+=in_strictcount[i]; // if node is terminal the counts should be adjusted
								if(nodeinfo[n].guide[i].terminal_out) sumstrict+=rate*out_strictcount[i];
							}
						}

						// update guide coverages with unique reads and gcounts for later
						for(int i=0;i<ng;i++) {

							if(!allstrictzero) {
								float strictcount=in_strictcount[i]+rate*out_strictcount[i];
								if(incount && outcount) {
									if(nodeinfo[n].guide[i].terminal_in) strictcount+=in_strictcount[i];
									if(nodeinfo[n].guide[i].terminal_out) strictcount+=rate*out_strictcount[i];
								}
								if(strictcount) { // if there are reads unique to the guide
									int g=guidetrf[nodeinfo[n].guide[i].idx].g;
									path[0]=n;
									nodeflux[0]=strictcount/sumstrict;
									//fprintf(stderr,"3 g=%d\n",g);
									if(guidepred[g]!=-1) update_guide_pred(pred,guidepred[g],path,nodeflux,nodecov,no2gnode,gno,false);
									else // new prediction for this guide
										guidepred[g]=store_guide_transcript(pred,path,nodeflux,nodecov,no2gnode,geneno,first,gno,guides[g],false);
								}
							}

							// now set gcounts
							nodeinfo[n].guide[i].gcount=in_loosecount[i]+rate*out_loosecount[i]+
									in_looseterminalcount[i]+rate*out_looseterminalcount[i];
							if(incount && outcount) {
								if(nodeinfo[n].guide[i].terminal_in) nodeinfo[n].guide[i].gcount+=in_looseterminalcount[i];
								if(nodeinfo[n].guide[i].terminal_out) nodeinfo[n].guide[i].gcount+=rate*out_looseterminalcount[i];
							}
						}

						if(!allstrictzero) { // update node coverage here
							if(allloosezero) nodecov[n]=0;
							else nodecov[n]*=nodeinfo[n].sumtrcount/sumstrict;
						}

						if(!allloosezero && nodecov[n]) { // get the pattern counts
							if(incount && outcount) for(int i=0;i<nodeinfo[n].trcount.Count();i++)
								if(nodeinfo[n].trcount[i].terminal) nodeinfo[n].trcount[i].abund*=2;
							for(int i=0;i<out_trcount.Count();i++) {
								int t=find_cguidepat(out_trcount[i].pat,nodeinfo[n].trcount);
								out_trcount[i].abund*=rate;
								if(out_trcount[i].terminal && incount && outcount) out_trcount[i].abund*=2;
								//else out_trcount[i].terminal=false; // I can not have a trcount with the same guides that is terminal both for in and out
								if(t==-1) { // didn't find the guide pattern
									nodeinfo[n].trcount.Add(out_trcount[i]);
								}
								else nodeinfo[n].trcount[t].abund+=out_trcount[i].abund;
							}
						}

					}
				}
			}
			if(nodecov[n]>epsilon) { // there is still coverage left for the node
				coverednode.Add(n);
			}
		}

		if(coverednode.Count()) { // I still have covered nodes: asign reads to guides; I need to do some sort of EM algorithm here
			// my prior bias in assigning reads should be based on coverages allready found

			int nEM=10; // number of times to repeat the EM algorithm: need to see if this is a good count, or if I should use and epsilon change to stop
			ng=guidetrf.Count();
			GVec<float> initgcov(ng,float(0));

			// first set the initial coverages in covered nodes
			for(int i=0;i<ng;i++) {
				int np=guidepred[guidetrf[i].g];
				if(np!=-1) {
					initgcov[i]=pred[np]->cov*pred[np]->tlen;
				}
			}

			GVec<float> prevgcov(ng,float(0)); // the new gcov that will be computed as part of the EM algorithm

			// now set the cov's in CPartGuide -> for sort purposes and for the EM algorithm
			for(int i=0;i<coverednode.Count();i++) {
				int n=coverednode[i];
				for(int j=0;j<nodeinfo[n].guide.Count();j++) {
					int g=nodeinfo[n].guide[j].idx;
					nodeinfo[n].guide[j].cov=initgcov[g];
					prevgcov[g]=initgcov[g];
				}
			}

			GVec<float> gcov(ng,float(0)); // the new gcov that will be computed as part of the EM algorithm

			int m=0;
			while(m<nEM) { // do the EM

				for(int i=0;i<ng;i++) gcov[i]=initgcov[i];

				for(int i=0;i<coverednode.Count();i++) { // for each node recompute probabilities
					int n=coverednode[i];
					nodeinfo[n].guide.Sort(partguideCmp); // I need to sort at each step because the cov's get updated
					float totalcount=0;
					bool covnotreached=true;
					float sumgcov=0;
					for(int j=0;j<nodeinfo[n].guide.Count();j++) { // for each guide in the node, starting from the most abundant to the least
						nodeinfo[n].guide[j].ncov=0; // I need to set this for each guide -> this is why I not break from the loop below
						int g=nodeinfo[n].guide[j].idx; // index in guidetrf
						sumgcov+=prevgcov[g];
						if(covnotreached) for(int t=0;t<nodeinfo[n].trcount.Count();t++) { // for each trguidepat in the node
							if(nodeinfo[n].trcount[t].pat[g]) { // current guide is in this pattern
								float sum=0; // sum of all coverages of guides in the transcript
								for(int k=0;k<nodeinfo[n].trcount[t].g.Count();k++) {
									int l=nodeinfo[n].trcount[t].g[k];
									sum+=prevgcov[l];
								}
								float abund=nodeinfo[n].sumtrcount-totalcount; // this is the maximum abundance allowed for this guide
								if(sum) { // I have some guides that have coverage
									float newabund=nodeinfo[n].trcount[t].abund*prevgcov[g]/sum;
									if(newabund<abund) abund=newabund;
								}
								else { // no guide has assigned coverages yet -> winner takes all
									if(nodeinfo[n].trcount[t].abund<abund) abund=nodeinfo[n].trcount[t].abund;
								}
								nodeinfo[n].guide[j].ncov+=abund;
								gcov[g]+=abund*nodeinfo[n].guide[j].olen/no2gnode[n]->len();
								totalcount+=abund;
								if(nodeinfo[n].sumtrcount-totalcount<epsilon) {
									covnotreached=false;
									break;
								}
							}
						}
					}
					if(covnotreached) { // there were not enough transfrags to give coverages -> redistribute reads based on coverages
						float abund=nodeinfo[n].sumtrcount-totalcount; // this is how much is left
						if(sumgcov) for(int j=0;j<nodeinfo[n].guide.Count();j++) {
							if(!nodeinfo[n].guide[j].cov) break;
							float newabund=abund*nodeinfo[n].guide[j].cov/sumgcov;
							nodeinfo[n].guide[j].ncov+=newabund;
							int g=nodeinfo[n].guide[j].idx;
							gcov[g]+=newabund*nodeinfo[n].guide[j].olen/no2gnode[n]->len();
						}
						else { // all guides have cov zero -> winner takes all
							nodeinfo[n].guide[0].ncov+=abund;
							int g=nodeinfo[n].guide[0].idx;
							gcov[g]+=abund*nodeinfo[n].guide[0].olen/no2gnode[n]->len();
						}
					}
				} // end node

				bool nochange=true;
				// check coverages
				for(int i=0;i<coverednode.Count();i++) {
					int n=coverednode[i];
					for(int j=0;j<nodeinfo[n].guide.Count();j++) {
						int g=nodeinfo[n].guide[j].idx;
						float diff=nodeinfo[n].guide[j].cov-gcov[g];
						if(diff<0) diff=-diff;
						if(diff) {
							prevgcov[g]=gcov[g];
							nodeinfo[n].guide[j].cov=gcov[g];
							if(diff>1) nochange=false;
						}
					}
				}

				// see if there is a change in probability
				if(nochange) break;

				m++;
			} // end EM

			// new coverages are estimated so now I can assign the reads to the predictions
			for(int g=0;g<guidetrf.Count();g++) if(gcov[g]>initgcov[g]) { // for each guide that has more coverage to add
				path.Clear();
				nodeflux.Clear();
				for(int i=0;i<guidetrf[g].trf->nodes.Count();i++) if(nodecov[guidetrf[g].trf->nodes[i]]) { // for each covered node
					int n=guidetrf[g].trf->nodes[i];
					for(int j=0;j<nodeinfo[n].guide.Count();j++) if(g==nodeinfo[n].guide[j].idx){ // found my guide
						if(nodeinfo[n].guide[j].ncov) { // if there is coverage assigned to node -> i should add it to the path
							path.Add(n);
							float prop=nodeinfo[n].guide[j].ncov/nodeinfo[n].sumtrcount;
							nodeflux.Add(prop);
						}
						break;
					}
				}
				//fprintf(stderr,"4 g=%d\n",guidetrf[g].g);
				if(guidepred[guidetrf[g].g]!=-1) update_guide_pred(pred,guidepred[guidetrf[g].g],path,nodeflux,nodecov,no2gnode,gno,false);
				else // new prediction for this guide
					guidepred[guidetrf[g].g]=store_guide_transcript(pred,path,nodeflux,nodecov,no2gnode,geneno,first,gno,guides[guidetrf[g].g],false);
			}
		}

		maxi=0;
	}

	return(maxi);
}

bool transcript_cont_path(GVec<int>& path,CTransfrag *transfrag, int &ni,GVec<int>& alltr,int t) {

	int np=path.Count();
	int n=0;
	while(path[n]<transfrag->nodes[0]) n++;

	for(int i=0;i<transfrag->nodes.Count();i++) {
		if(n==np) { // I reached the end of search and didn't find any mismatch
			ni=i;
			return true;
		}
		if(path[n]!=transfrag->nodes[i]) return false;
		n++;
	}

	if(n==np) { // transfrag is completely included in path
		transfrag->abundance=0;
		alltr.Add(t);
		return false;
	}

	return true;
}

void collect_path(GPVec<CMTransfrag>& mgt,GVec<int>& alltr,GVec<int>& path,GPVec<CGraphnode>& no2gnode,int n, int gno) {

	CGraphnode *node=no2gnode[path[n]];

	if(node->child.Count()==1 && node->child[0]==gno-1) return; // only child of node is sink

	float maxabund=0;
	int maxt=-1;
	bool found=false;
	int nexti=-1;
	for(int j=0;j<node->trf.Count();j++) {
		int t=node->trf[j];
		if(mgt[t]->transfrag->abundance && mgt[t]->transfrag->nodes[0] &&
				(!mgt[t]->transfrag->real || mgt[t]->transfrag->nodes[0]==path[n] )) { // if t is a guide it should start at this node to be included in transcript
			int ni=-1;
			if(!found) {
				if(mgt[t]->transfrag->real && transcript_cont_path(path,mgt[t]->transfrag,ni,alltr,t)) {
					maxt=t;
					nexti=ni;
					found=true; // found my continuation transcript -> no need to look further
				}
				else if(!mgt[t]->transfrag->real && mgt[t]->transfrag->abundance>maxabund && transcript_cont_path(path,mgt[t]->transfrag,ni,alltr,t)) {
					maxt=t;
					nexti=ni;
					maxabund=mgt[t]->transfrag->abundance;
				}
			}
			else transcript_cont_path(path,mgt[t]->transfrag,ni,alltr,t);
		}
	}

	if(maxt==-1) return; // path ends here
	// found maxt -> extend path
	alltr.Add(maxt);
	if(!mgt[maxt]->transfrag->real) mgt[maxt]->transfrag->abundance=0;
	n++;
	for(int i=nexti;i<mgt[maxt]->transfrag->nodes.Count();i++) {
		path.Add(mgt[maxt]->transfrag->nodes[i]);
	}

}

/*
int merge_transcripts(int gno, GPVec<CGraphnode>& no2gnode,GPVec<CMTransfrag>& mgt,
		int geneno,int strand,GList<CPrediction>& pred,GList<CReadAln>& readlist,GPVec<GffObj>& guides) {

	bool first=true;

	// process all children's source
	CGraphnode *source=no2gnode[0];
	for(int c=0;c<source->child.Count();c++) { // for all children
		CGraphnode *node=no2gnode[source->child[c]];

		// process all transfrags leaving from source if they haven't been used already
		for(int j=0;j<node->trf.Count();j++) {
			int t=node->trf[j];


			{ // DEBUG ONLY
				fprintf(stderr,"Source child %d with transcript(%f):",source->child[c],mgt[t]->transfrag->abundance);
				for(int i=0;i<mgt[t]->transfrag->nodes.Count();i++) fprintf(stderr," %d",mgt[t]->transfrag->nodes[i]);
				fprintf(stderr,"\n");
			}


			if(mgt[t]->transfrag->abundance && mgt[t]->transfrag->nodes[0]) { // if I didn't use this transfrag
				GVec<int> alltr(1,t); // all transfrags used in creating path
				if(mgt[t]->transfrag->real) { // this is a guide transfrag -> always print
					geneno=store_merge_prediction(pred,alltr,mgt,mgt[t]->transfrag->nodes,no2gnode,strand,geneno,first,readlist,guides);
				}
				// see if I can build another path for it
				GVec<int> path(mgt[t]->transfrag->nodes);
				if(!mgt[t]->transfrag->real) mgt[t]->transfrag->abundance=0;
				if(path.Count()>1) collect_path(mgt,alltr,path,no2gnode,1,gno);
				// else this is a single exon transcript that is not included in another one
				if(!mgt[t]->transfrag->real || path.Count()>mgt[t]->transfrag->nodes.Count()) // found a different path
					geneno=store_merge_prediction(pred,alltr,mgt,path,no2gnode,strand,geneno,first,readlist,guides);
			}
		}
	}
	return(geneno);
}
*/

int find_transcripts(int gno,int edgeno, GIntHash<int> &gpos,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,GBitVec& compatible,
		int geneno,int strand,GVec<CGuide>& guidetrf,GList<CPrediction>& pred,GPVec<GffObj>& guides,GVec<int>& guidepred,BundleData* bdata) {

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

		    if(abundin) inode->rate=abundout/abundin;
		    if(abundout) inode->capacity=abundout+abundthrough; // node capacity tells me how much of that node coverage I can use given how many transfrags leave the node
		    else inode->capacity=abundin+abundthrough;
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

	if(eonly || guidetrf.Count()) maxi=guides_pushmaxflow(gno,edgeno,gpos,no2gnode,transfrag,guidetrf,geneno,strand,pred,nodecov,istranscript,pathpat,first,guides,guidepred,bdata);


	if(nodecov[maxi]>=readthr) {

		// process rest of the transfrags
		if(!eonly && nodecov[maxi]>=readthr) {
			GBitVec removable(transfrag.Count(),true);

			// 1:
			// parse_trf_weight_max_flow(gno,no2gnode,transfrag,geneno,strand,pred,nodecov,pathpat);
			// 2:
			GBitVec usednode(gno+edgeno);
			parse_trf(maxi,gno,edgeno,gpos,no2gnode,transfrag,compatible,geneno,first,strand,pred,nodecov,istranscript,removable,usednode,0,pathpat);

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
					if(trlen>=mintranscriptlen || (lex==fex && maxlen>mintranscriptlen)) {
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
							if(trlen>=mintranscriptlen) {
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
						if(trlen>=mintranscriptlen || (lex==fex && maxlen>=mintranscriptlen)) {
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
			if(maxlen>=mintranscriptlen) {
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

bool good_junc(CJunction& jd,int refstart, GVec<float>* bpcov) {

	if(eonly && !jd.guide_match) { // this way I am using only reads that are compatible to annotated transcripts
		jd.strand=0;
		return false;
	}

	if(jd.guide_match) return true; // this junction is covered by at least one read: the one that calls good_junc
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
	//if (jd.leftsupport<junctionthr || jd.rightsupport<junctionthr || (jd.end-jd.start>longintron && jd.nreads_good<junctionthr)) {
		jd.strand=0;
		//fprintf(stderr,"nosupport: left=%f right=%f good=%f\n",jd.leftsupport,jd.rightsupport,jd.nreads_good);
		return false;
	}

	// if(!jd.strand) return false; the strand has to be non-zero when we call good_junc -> KEEP THIS IN MIND IF WE CHANGE HOW WE USE GOOD_JUNC

	// don't trust spliced reads that have a very low coverage:
	int sno=(int)jd.strand+1;


	float leftcov=bpcov[sno][jd.start-refstart-1];
	if(bpcov[1][jd.start-refstart-1]>leftcov) {
		if(bpcov[2-sno][jd.start-refstart-1])
			leftcov+=(bpcov[1][jd.start-refstart-1]-bpcov[0][jd.start-refstart-1]-bpcov[2][jd.start-refstart-1])*
				leftcov/(leftcov+bpcov[2-sno][jd.start-refstart-1]);
		else leftcov=bpcov[1][jd.start-refstart-1];
	}

	float rightcov=bpcov[sno][jd.start-refstart];
	if(bpcov[1][jd.start-refstart]>rightcov) {
		if(bpcov[2-sno][jd.start-refstart])
			rightcov+=(bpcov[1][jd.start-refstart]-bpcov[0][jd.start-refstart]-bpcov[2][jd.start-refstart])*
				rightcov/(rightcov+bpcov[2-sno][jd.start-refstart]);
		else rightcov=bpcov[1][jd.start-refstart];
	}

	//if(leftcov && jd.nreads_good*100/leftcov<isofrac && rightcov/leftcov>1-isofrac) {
	if(leftcov && jd.leftsupport*100/leftcov<isofrac && rightcov/leftcov>1-isofrac) {
		jd.strand=0;
		//fprintf(stderr,"no leftsupport\n");
		return false;
	}

	leftcov=bpcov[sno][jd.end-refstart];
	if(bpcov[1][jd.end-refstart]>leftcov) {
		if(bpcov[2-sno][jd.end-refstart]>leftcov)
			leftcov+=(bpcov[1][jd.end-refstart]-bpcov[0][jd.end-refstart]-bpcov[2][jd.end-refstart])*
				leftcov/(leftcov+bpcov[2-sno][jd.end-refstart]);
		else leftcov=bpcov[1][jd.end-refstart];
	}

	rightcov=bpcov[sno][jd.end-refstart-1];
	if(bpcov[1][jd.end-refstart-1]>rightcov) {
		if(bpcov[2-sno][jd.end-refstart-1])
			rightcov+=(bpcov[1][jd.end-refstart-1]-bpcov[0][jd.end-refstart-1]-bpcov[2][jd.end-refstart-1])*
				rightcov/(rightcov+bpcov[2-sno][jd.end-refstart-1]);
		else rightcov=bpcov[1][jd.end-refstart-1];
	}

	//if(leftcov && jd.nreads_good*100/leftcov<isofrac && rightcov/leftcov>1-isofrac) {
	if(leftcov && jd.rightsupport*100/leftcov<isofrac && rightcov/leftcov>1-isofrac) {
		jd.strand=0;
		//fprintf(stderr,"no rightsupport\n");
		return false;
	}

	if (jd.nm && round(jd.nm)==round(jd.nreads) && jd.end-jd.start>longintron) { // don't believe long intron if all junctions are from bad reads
		jd.strand=0;
		jd.nm=0;
		//fprintf(stderr,"bad longintron\n");
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


bool good_merge_junc(CJunction& jd,GList<CJunction>& junction) {

	if(jd.guide_match) return true; // this junction is covered by at least one read: the one that calls good_junc

	// if(!jd.strand) return false; the strand has to be non-zero when we call good_junc -> KEEP THIS IN MIND IF WE CHANGE HOW WE USE GOOD_JUNC

	// don't trust junctions that have a very low coverage: (note: junctions are stored by start)
	int nj=junction.IndexOf(&jd);
	if(nj<0 || !jd.nreads_good) {
	//if(nj<0 || !jd.leftsupport || !jd.rightsupport) {
		jd.strand=0;
		return false;
	}

	int j=nj-1;
	float sumjunc=jd.nreads_good;
	while(j>=0 && junction[j]->start==jd.start) {
		sumjunc+=junction[j]->nreads_good; j--;
		/*if(junction[j]->end==jd.end && junction[j]->strand && junction[j]->strand!=jd.strand) {
			fprintf(stderr,"Found junction on strand %d matching to %d-%d\n",junction[j]->strand,junction[j]->start,junction[j]->end);
			if(jd.nreads_good<junction[j]->nreads_good || junction[j]->guide_match) {
				jd.strand=0;
				return false;
			}
			else junction[j]->strand=0;
		}*/
	}
	j=nj+1;
	while(j<junction.Count() && junction[j]->start==jd.start) {
		sumjunc+=junction[j]->nreads_good; j++;
		/*if(junction[j]->end==jd.end && junction[j]->strand && junction[j]->strand!=jd.strand) {
			fprintf(stderr,"Found junction on strand %d matching to %d-%d\n",junction[j]->strand,junction[j]->start,junction[j]->end);
			if(jd.nreads_good<junction[j]->nreads_good || junction[j]->guide_match) {
				jd.strand=0;
				return false;
			}
			else junction[j]->strand=0;
		}*/
	}

	// this sorts some of the really bad junctions out
	if(jd.nreads_good*100/sumjunc<isofrac) {
		jd.strand=0;
		return false;
	}


	return true;
}


void continue_read(GList<CReadAln>& readlist,int n,int idx) {
	// keep longest part of read
	readlist[n]->end=readlist[n]->segs[idx].end;
	for(int i=readlist[n]->segs.Count()-1;i>idx;i--) {
		readlist[n]->segs.Delete(i);
		readlist[n]->juncs.Delete(i-1);
	}
}


void adjust_read(CReadAln &rd,int rstart,int rend, uint rlen) {
	rd.start=rd.segs[rstart].start;
	rd.end=rd.segs[rend].end;
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

	rd.len=rlen;

	// clean up read after rend
	for(int i=nex-1;i>rend;i--) {

		if(rd.juncs[i-1]->nm && long(rd.juncs[i-1]->nm)>=long(rd.juncs[i-1]->nreads) && (rd.juncs[i-1]->end - rd.juncs[i-1]->start > longintron)) { // only the first time I am interested in checking if a junction is bad
			rd.juncs[i-1]->strand=0;
		}
		rd.juncs[i-1]->nm=0;

		rd.juncs[i-1]->nreads-=rd.read_count;
		if(!rd.juncs[i-1]->nreads) rd.juncs[i-1]->strand=0;

		uint anchor=junctionsupport;
		if(anchor<longintronanchor && rd.juncs[i-1]->len()>longintron) anchor=longintronanchor; // I want to use a longer anchor for long introns to believe them
		if(leftsup[i-1]>=anchor) {
			rd.juncs[i-1]->leftsupport-=rd.read_count;
			if(rightsup[nex-i-1]>=anchor) { // this was a good junction but for some reason the read was thrown away
				rd.juncs[i-1]->nreads_good-=rd.read_count;
				rd.juncs[i-1]->rightsupport-=rd.read_count;
			}
		}
		else if(rightsup[nex-i-1]>=anchor) { // this was a good junction but for some reason the read was thrown away
			rd.juncs[i-1]->rightsupport-=rd.read_count;
		}

		rd.juncs.Pop();
		rd.segs.Delete(i);

	}

	// move down in arrays elements after rstart
	if(rstart) { // the portion I keep from the read doesn't start at first exon
		for(int i=0;i<rend-rstart+1;i++) {
			if(i<rstart) { // this is a part that has a junction that needs to be removed
				if(rd.juncs[i]->nm && long(rd.juncs[i]->nm)>=long(rd.juncs[i]->nreads) && (rd.juncs[i]->end - rd.juncs[i]->start > longintron)) { // only the first time I am interested in checking if a junction is bad
					rd.juncs[i]->strand=0;
				}
				rd.juncs[i]->nm=0;

				rd.juncs[i]->nreads-=rd.read_count;
				if(!rd.juncs[i]->nreads) rd.juncs[i]->strand=0;

				uint anchor=junctionsupport;
				if(anchor<longintronanchor && rd.juncs[i]->len()>longintron) anchor=longintronanchor; // I want to use a longer anchor for long introns to believe them
				if(leftsup[i]>=anchor) {
					rd.juncs[i]->leftsupport-=rd.read_count;
					if(rightsup[nex-i-2]>=anchor) { // this was a good junction but for some reason the read was thrown away
						rd.juncs[i]->nreads_good-=rd.read_count;
						rd.juncs[i]->rightsupport-=rd.read_count;
					}
				}
				else if(rightsup[nex-i-2]>=anchor) { // this was a good junction but for some reason the read was thrown away
					rd.juncs[i]->rightsupport-=rd.read_count;
				}
			}
			// copy the info from position i+rstart here
			rd.segs[i]=rd.segs[i+rstart];
			if(i+rstart<rend) rd.juncs.Put(i,rd.juncs[i+rstart]);
		}

		// clean up elements after rend-rstart
		for(int i=rend;i>rend-rstart;i--) {
			rd.juncs.Pop();
			rd.segs.Delete(i);
		}
	}

}

void update_junction_counts(CReadAln & rd) {
	int nex=rd.segs.Count();
	/* this part was removed to improve performance
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
	*/
	for(int i=1;i<nex;i++) {

		if(rd.juncs[i-1]->nm && long(rd.juncs[i-1]->nm)>=long(rd.juncs[i-1]->nreads) && (rd.juncs[i-1]->end - rd.juncs[i-1]->start > longintron)) { // only the first time I am interested in checking if a junction is bad
			rd.juncs[i-1]->strand=0;
		}
		rd.juncs[i-1]->nm=0;

		rd.juncs[i-1]->nreads-=rd.read_count;
		if(!rd.juncs[i-1]->nreads) rd.juncs[i-1]->strand=0;

		/* this part was removed to improve performance
		uint anchor=junctionsupport;
		if(anchor<longintronanchor && rd.juncs[i-1]->len()>longintron) anchor=longintronanchor; // I want to use a longer anchor for long introns to believe them
		if(leftsup[i-1]>=anchor) {
			rd.juncs[i-1]->leftsupport-=rd.read_count;
			if(rightsup[nex-i-1]>=anchor) { // this was a good junction but for some reason the read was thrown away
				rd.juncs[i-1]->nreads_good-=rd.read_count;
				rd.juncs[i-1]->rightsupport-=rd.read_count;
			}
		}
		else if(rightsup[nex-i-1]>=anchor) {
			rd.juncs[i-1]->rightsupport-=rd.read_count;
		}
		*/
	}
}

/*
void assign_strand(CReadAln &rd,GList<CReadAln>& readlist) {
	int newstrand=0;
}
*/

int build_graphs(BundleData* bdata) {
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
	GVec<int> guidepred; // for eonly keeps the prediction number associated with a guide
	GArray<GEdge> guideedge; // 0: negative starts; 1 positive starts

	if(guides .Count()) {

		guideedge.setSorted(true);
		guideedge.setUnique(true);

		//if(eonly)
		for(int g=0;g<guides.Count();g++) {
			guidepred.cAdd(-1);
			bool covered=true;
			RC_TData* tdata=(RC_TData*)(guides[g]->uptr);
			for(int i=0;i<tdata->t_introns.Count();i++) {
				if(!tdata->t_introns[i]->rcount) {
					covered=false;
					break;
				}
			}
			if(covered) {
				tdata->in_bundle=2;
				int s=-1; // unknown strand
				if(guides[g]->strand=='+') s=1; // guide on positive strand
				else if(guides[g]->strand=='-') s=0; // guide on negative strand

				//fprintf(stderr,"Look to add guide g=%d start %d-%d and end %d-%d on strand %d\n",g,guides[g]->start,guides[g]->exons[0]->end,guides[g]->end,guides[g]->exons.Last()->start,s);

				int uses=s;
				if(s<0) uses=0;
				GEdge ge(guides[g]->start,guides[g]->exons[0]->end,uses);
				int idx=guideedge.IndexOf(ge);

				//fprintf(stderr,"look for ge(%d,%d,%d) => start idx=%d\n",ge.val,ge.endval,ge.strand,idx);

				if(idx<0) guideedge.Add(ge);
				else if(guideedge[idx].endval>guides[g]->exons[0]->end) guideedge[idx].endval=guides[g]->exons[0]->end;
				if(s<0) {
					ge.strand=1;
					idx=guideedge.IndexOf(ge);
					if(idx<0) guideedge.Add(ge);
					else if(guideedge[idx].endval>guides[g]->exons[0]->end) guideedge[idx].endval=guides[g]->exons[0]->end;
				}

				ge.val=guides[g]->end;
				ge.endval=guides[g]->exons.Last()->start;
				idx=guideedge.IndexOf(ge);

				//fprintf(stderr,"look for ge(%d,%d,%d) => end idx=%d\n",ge.val,ge.endval,ge.strand,idx);

				if(idx<0) guideedge.Add(ge);
				else if(guideedge[idx].endval<guides[g]->exons.Last()->start) guideedge[idx].endval=guides[g]->exons.Last()->start;
				if(s<0) {
					ge.strand=0; // ge.strand was set as 1 before
					idx=guideedge.IndexOf(ge);
					if(idx<0) guideedge.Add(ge);
					else if(guideedge[idx].endval<guides[g]->exons.Last()->start) guideedge[idx].endval=guides[g]->exons.Last()->start;
				}
			}
		}
	}

	// this part is for adjusting leftsupport and rightsupport when considering all junctions that start at a given point
	// sort junctions -> junctions are sorted already according with their start, but not their end
	GList<CJunction> ejunction(junction);
	ejunction.setFreeItem(false);
	if(ejunction.Count()) ejunction.setSorted(juncCmpEnd);

	uint start=0;
	uint end=0;
	double leftsupport=0;
	double rightsupport=0;
	for(int i=0;i<junction.Count();i++) {

		if(junction[i]->start!=start) {
			int j=i-1;
			while(j>=0 && junction[j]->start==start) {
				junction[j]->leftsupport=leftsupport;
				j--;
			}
			leftsupport=junction[i]->leftsupport;
			start=junction[i]->start;
		}
		else leftsupport+=junction[i]->leftsupport;

		if(junction[i]->end!=end) {
			int j=i-1;
			while(j>=0 && junction[j]->end==end) {
				junction[j]->rightsupport=rightsupport;
				j--;
			}
			rightsupport=junction[i]->rightsupport;
			end=junction[i]->end;
		}
		else rightsupport+=junction[i]->rightsupport;

	}
	// end adjusting leftsupport and rightsupport

	/*
	{ // DEBUG ONLY
		for(int i=0;i<junction.Count();i++) {
			fprintf(stderr,"Junction %d-%d has leftsupport=%f and rightsupport=%f\n",junction[i]->start,junction[i]->end,junction[i]->leftsupport,junction[i]->rightsupport);
		}
	}
	*/

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

		//if(rd.strand==0 && rd.pair_idx.Count()) assign_strand(rd,readlist);

		//if(rd.juncs.Count())
		//fprintf(stderr,"process read[%d] %d-%d:%d and pairs:",n,rd.start,rd.end,int(rd.strand));
		//for(int i=0;i<rd.pair_idx.Count();i++) fprintf(stderr," %d",rd.pair_idx[i]);
		//fprintf(stderr,"\n");

		/*
		// this version tried to adjust the read by eliminating bad junctions but makes the add_read_to_group crash because reads are no longer sorted by start coordinate
		int i=0;
		uint maxlen=rd.segs[0].len();
		uint seenlen=maxlen;
		int maxstart=0;
		int maxend=0;
		start=0;
		end=0;
		GVec<int> badjunc;
		while(i<rd.juncs.Count()) {
			CJunction& jd=*(rd.juncs[i]);

			if(!jd.strand || !good_junc(jd,refstart,bpcov)) { // found bad junction

				fprintf(stderr,"bad junc:%d-%d\n",jd.start,jd.end);

				badjunc.Add(i);
				if(seenlen>maxlen) {
					maxstart=start;
					maxend=end;
					maxlen=seenlen;
				}
				seenlen=rd.segs[i+1].len();
				start=i+1;
				end=i+1;
			}
			else { // good junction
				seenlen+=rd.segs[i+1].len();
				end++;
				if(guides.Count()){ // need to remember boundary
					bool exist=true;
					GStr bs((int)jd.start);
					if(!boundaryleft[bs.chars()]) boundaryleft.Add(bs.chars(),new bool(exist));
					GStr be((int)jd.end);
					if(!boundaryright[be.chars()]) boundaryright.Add(be.chars(),new bool(exist));
				}
			}
			i++;
		}

		if(seenlen>maxlen) {
			maxstart=start;
			maxend=end;
			maxlen=seenlen;
		}

		if(badjunc.Count()) { // if there are any bad junctions, I need to adjust the read
			adjust_read(rd,maxstart,maxend,maxlen);
		}
		*/

		// version that does not adjust read but discards it completely
		bool keep=true;
		int i=0;
		//int leftlen=0;
		while(i<rd.juncs.Count()) {
			//leftlen=rd.segs[i].len();
			CJunction& jd=*(rd.juncs[i]);
			//fprintf(stderr, " %d-%d:%4.2f nreads=%f nm=%f\n", jd.start, jd.end, jd.nreads_good,jd.nreads,jd.nm);
			if(!jd.strand || !good_junc(jd,refstart,bpcov)) { // found a bad junction
				//update_junction_counts(rd,bpcov,refstart);
				update_junction_counts(rd);
				keep=false;
				rd.nh=0; // I could make it negative to mark somehow that the read was not added to a group -> define as unreliable
				//fprintf(stderr," - %c -",jd.strand);
				break;
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
		}
		//else { fprintf(stderr,"read[%d] is not kept\n",n);}
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

					//fprintf(stderr,"sno=%d lastgroup->end=%d procgroup->start=%d procgroup->end=%d\n",sno,lastgroup->end,procgroup->start,procgroup->end);

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

			if(prevgroup[0]!=NULL && currgroup[nextgr]->start <= prevgroup[0]->end+bundledist) { // overlaps previous negative group ; this needs bundledist
				//fprintf(stderr,"\tovlp to neg group: %u-%u\n",prevgroup[0]->start,prevgroup[0]->end);
				set_strandcol(currgroup[nextgr],prevgroup[0],prevgroup[0]->color,eqnegcol,equalcolor);
				uint maxstart = currgroup[nextgr]->start > prevgroup[0]->start ? currgroup[nextgr]->start : prevgroup[0]->start;
				uint minend = currgroup[nextgr]->end < prevgroup[0]->end ? currgroup[nextgr]->end : prevgroup[0]->end;
				if(minend<maxstart) minend=maxstart; // this can only happen if bundledist >0
				currgroup[nextgr]->neg_prop+=prevgroup[0]->cov_sum*(minend-maxstart+1)/prevgroup[0]->len();
			}

			while(currgroup[0]!=NULL && currgroup[nextgr]->start <= currgroup[0]->end+bundledist && currgroup[0]->start <= currgroup[nextgr]->end +bundledist) { // overlaps current negative strand group
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
				if(minend<maxstart) minend=maxstart;
				currgroup[nextgr]->neg_prop+=currgroup[0]->cov_sum*(minend-maxstart+1)/currgroup[0]->len();


				prevgroup[0]=currgroup[0];
				currgroup[0]=currgroup[0]->next_gr;
			}

			float pos_prop=0;
			if(prevgroup[2]!=NULL && currgroup[nextgr]->start <= prevgroup[2]->end + bundledist) { // overlaps positive strand group
				//fprintf(stderr,"\tovlp to pos group: %u-%u\n",prevgroup[2]->start,prevgroup[2]->end);
				set_strandcol(currgroup[nextgr],prevgroup[2],prevgroup[2]->color,eqposcol,equalcolor);
				if(currgroup[nextgr]->neg_prop) {
					uint maxstart = currgroup[nextgr]->start > prevgroup[2]->start ? currgroup[nextgr]->start : prevgroup[2]->start;
					uint minend = currgroup[nextgr]->end < prevgroup[2]->end ? currgroup[nextgr]->end : prevgroup[2]->end;
					if(minend<maxstart) minend=maxstart; // this can only happen if bundledist >0
					pos_prop+=prevgroup[2]->cov_sum*(minend-maxstart+1)/prevgroup[2]->len();
				}
			}

			while(currgroup[2]!=NULL && currgroup[nextgr]->start <= currgroup[2]->end +bundledist && currgroup[2]->start <= currgroup[nextgr]->end + bundledist) { // overlaps positive strand group
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
					if(minend<maxstart) minend=maxstart; // this can only happen if bundledist >0
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
				add_group_to_bundle(currgroup[nextgr],bundle[nextgr][bno],bnode[nextgr],bundledist);
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
					add_group_to_bundle(currgroup[nextgr],bundle[0][bno],bnode[0],bundledist); // this needs bundledist
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
					add_group_to_bundle(currgroup[nextgr],bundle[2][bno],bnode[2],bundledist);
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

		//print STDERR "Done with groups: ",$currgroup[0],",",$currgroup[1],",",$currgroup[2],"\n";
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

	//if(guides.Count()) fprintf(stderr,"No of guides=%d\n",guides.Count());

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
				int nodes=0;
				while(currbnode!=NULL) {
					fprintf(stderr, " %d-%d cov=%.6f",currbnode->start,currbnode->end,currbnode->cov/(currbnode->end-currbnode->start+1));
					currbnode=currbnode->nextnode;
					nodes++;
				}
				currbnode=bnode[sno][bundle[sno][b]->lastnodeid];
				nodes++;
				fprintf(stderr," last node:%d-%d total nodes=%d\n",currbnode->start,currbnode->end,nodes);
			}
		}
	}
	*/

	int geneno=0;


    // ### predict transcripts for unstranded bundles here
	//if(fraglen)
	for(int b=0;b<bundle[1].Count();b++) if (bundle[1][b]->nread) {
		double multi_perc=bundle[1][b]->multi/bundle[1][b]->nread;
    	if (multi_perc>mcov) {
    		if (verbose) GMessage("  bundle %s:%d-%d rejected due to multi-mapped reads %.3f% > threshold %.3f%\n",
    				 bdata->refseq.chars(), bnode[1][bundle[1][b]->startnode]->start,
					 bnode[1][bundle[1][b]->lastnodeid]->end, multi_perc, mcov);
    	} else if (guides.Count() || bundle[1][b]->len >= mintranscriptlen) {
		   // there might be small transfrags that are worth showing, but here I am ignoring them
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
    					tdata->in_bundle=3;
    					float gcov=(tdata->t_exons[0])->movlcount/glen;
    					// if(cov<gcov) gcov=cov; WHY DO I DO THIS?? CHECK!!!
    					CPrediction *p=new CPrediction(geneno-1, guides[g], guides[g]->start, guides[g]->end, gcov, guides[g]->strand, glen);
    					GSeg exon(guides[g]->start, guides[g]->end);
    					p->exons.Add(exon);
    					p->exoncov.Add(gcov);
    					pred.Add(p);
    					printguides=true;
    					guidepred[g]=pred.Count()-1;
    				}
    			}

    			/* char *predid=NULL;
    			if(!id.is_empty()) {
    				predid=Gstrdup(id.chars());
    			}
    			if((!eonly && cov>=readthr && len>mintranscriptlen) || (!id.is_empty())) {
    			*/
    			if(!printguides && cov>=readthr && len>=mintranscriptlen) {
    				if(t==1) { geneno++;}
    				char sign='.';

    				/* if I want to use trimming for the neutral single exons --> doesn't seem to work that well since single exons
    				 * tend to have sparse coverings

    				float sprop=0;
    				uint nodestart=0;
    				uint nodeend=0;
    				float sourceabundance=0;
    				float sinkabundance=0;
    				float sourcecovleft=0;
    				float sourcecovright=0;
    				float sinkcovleft=0;
    				float sinkcovright=0;
    				find_trims(refstart,1,currbnode->start,currbnode->end,bpcov,nodestart,sourceabundance,nodeend,sinkabundance,sourcecovleft,sourcecovright,sinkcovleft,sinkcovright);
    				if(nodestart) {
    					if(sourcecovleft) sprop=sourcecovleft/(sourcecovleft+sourcecovright);
    				}
    				else nodestart=currbnode->start;

    				float eprop=0;
    				if(nodeend) {
    					if(sinkcovright) eprop=sinkcovright/(sinkcovleft+sinkcovright);
    				}
    				else nodeend=currbnode->end;

    				if(nodeend>nodestart && (int)(nodeend-nodestart)>mintranscriptlen && (eprop || sprop)) {
    					cov-=(sprop+eprop)*cov;
    				}
    				else {
    					nodestart=currbnode->start;
    					nodeend=currbnode->end;
    				}

    				//CPrediction *p=new CPrediction(geneno-1,predid,currbnode->start,currbnode->end,cov,sign,cov,fraglen);
    				CPrediction *p=new CPrediction(geneno-1, NULL, nodestart, nodeend, cov, sign, cov);
    				GSeg exon(nodestart,nodeend);*/

    				CPrediction *p=new CPrediction(geneno-1, NULL, currbnode->start, currbnode->end, cov, sign, len);
    				GSeg exon(currbnode->start,currbnode->end);

    				p->exons.Add(exon);
    				p->exoncov.Add(cov);
    				pred.Add(p);
    				t++;
    			}
    			currbnode=currbnode->nextnode;
    		}
    	} //processing bundle
    } //non-empty bundle

    //fprintf(stderr,"Done with unstranded bundles\n");
    if (bnodeguides) delete[] bnodeguides;

	// ### build graphs for stranded bundles here
    if(startgroup[0]!=NULL || startgroup[2]!=NULL) { //# there are stranded groups to process

    	// I don't need the groups here anymore : I only use their numbers
    	// group.Clear(); // I will need the proportions later on

    	/* this is when I use nreads_good
    	// sort junctions -> junctions are sorted already according with their start, but not their end
    	GList<CJunction> ejunction(junction);
    	ejunction.setFreeItem(false);
    	if(ejunction.Count()) ejunction.setSorted(juncCmpEnd);
    	*/

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

    				//fprintf(stderr,"Bundle is: %d - %d start at g=%d sno=%d b=%d\n",bnode[sno][bundle[sno][b]->startnode]->start,bnode[sno][bundle[sno][b]->lastnodeid]->end,g,sno,b);

    				while(g<ng && guides[g]->end<bnode[sno][bundle[sno][b]->startnode]->start) g++;

    				int cg=g;
    				int nolap=0;
    				while(cg<ng && guides[cg]->start<=bnode[sno][bundle[sno][b]->lastnodeid]->end) { // this are potential guides that might overlap the current bundle, and they might introduce extra edges

    					//fprintf(stderr,"...consider guide cg=%d with strand=%c and in_bundle=%d\n",cg,guides[cg]->strand,((RC_TData*)(guides[cg]->uptr))->in_bundle);
    					if((guides[cg]->strand==strnd || guides[cg]->strand=='.') && ((RC_TData*)(guides[cg]->uptr))->in_bundle>=2) {
    						//fprintf(stderr,"Add guide g=%d with start=%d end=%d\n",cg,guides[cg]->start,guides[cg]->end);
    						edgeno[s][b]+=2; // this is an overestimate: possibly I have both an extra source and an extra sink link
    						nolap++;
    					}
    					cg++;
    				}

    				/*
    				{ // DEBUG ONLY
        				fprintf(stderr,"edgeno[%d][%d]=%d\n",s,b,edgeno[s][b]);
    					if(bundle[sno][b]->nread) {
    						fprintf(stderr,"proc bundle[%d][%d] %f/%f is %f len=%d and %d guides\n",sno,b,
    								bundle[sno][b]->multi,bundle[sno][b]->nread,(float)bundle[sno][b]->multi/bundle[sno][b]->nread,
    								bundle[sno][b]->len,nolap);
    					}
    				}
    				*/

    				// here I can add something in stringtie to lower the mintranscript len if there are guides?

    				if(bundle[sno][b]->nread &&
    						(((bundle[sno][b]->multi/bundle[sno][b]->nread)<=mcov && bundle[sno][b]->len >= mintranscriptlen)
    								||nolap)) { // bundle is worth processing: it might be that there are small transfrags from source to sink that are worth processing

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
    							bundle2graph,no2gnode,transfrag,gpos,bpcov,edgeno[s][b],lastgpos[s][b],guideedge); // also I need to remember graph coverages somewhere -> probably in the create_graph procedure

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

    				/*
    				{ // DEBUG ONLY
    					fprintf(stderr,"process refguides for s=%d b=%d edgeno=%d gno=%d lastgpos=%d\n",s,b,edgeno[s][b],graphno[s][b],lastgpos[s][b]);
    					fprintf(stderr,"There are %d nodes for graph[%d][%d]:\n",graphno[s][b],s,b);
    					for(int i=0;i<graphno[s][b];i++) {
    						fprintf(stderr,"%d (%d-%d): %f len=%d cov=%f",i,no2gnode[s][b][i]->start,no2gnode[s][b][i]->end,no2gnode[s][b][i]->cov,no2gnode[s][b][i]->len(),no2gnode[s][b][i]->cov/no2gnode[s][b][i]->len());
    						fprintf(stderr," parents:");
    						for(int j=0;j<no2gnode[s][b][i]->parent.Count();j++) fprintf(stderr," %d",no2gnode[s][b][i]->parent[j]);
    						fprintf(stderr," trf=");
    						for(int j=0;j<no2gnode[s][b][i]->trf.Count();j++) fprintf(stderr," %d",no2gnode[s][b][i]->trf[j]);
    						fprintf(stderr,"\n");
    					}
    				}
    				*/

    				if(guides.Count()) process_refguides(graphno[s][b],edgeno[s][b],gpos[s][b],lastgpos[s][b],no2gnode[s][b],transfrag[s][b],s,guides,guidetrf,bdata);

    				//process transfrags to eliminate noise, and set compatibilities, and node memberships
    				GBitVec compatible;
    				process_transfrags(graphno[s][b],edgeno[s][b],no2gnode[s][b],transfrag[s][b],tr2no[s][b],compatible,gpos[s][b],guidetrf);

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
    						//printBitVec(transfrag[s][b][t]->pattern);
    						fprintf(stderr," %f nodes=%d",transfrag[s][b][t]->abundance, transfrag[s][b][t]->nodes.Count());
    						if(!transfrag[s][b][t]->abundance) fprintf(stderr," *");
    						fprintf(stderr,"\n");
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

    				//fprintf(stderr,"guidetrf no=%d\n",guidetrf.Count());

    				// find transcripts now
    				geneno=find_transcripts(graphno[s][b],edgeno[s][b],gpos[s][b],no2gnode[s][b],transfrag[s][b],compatible,
    						geneno,s,guidetrf,pred,guides,guidepred,bdata);

    				for(int g=0;g<guidetrf.Count();g++) delete guidetrf[g].trf;


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

    /*
    { // DEBUG ONLY
    	for(int i=0;i<pred.Count();i++) {
    		if(pred[i]->t_eq) fprintf(stderr,"%s ",pred[i]->t_eq->getID());
    		fprintf(stderr,"pred[%d]:",i);
    		for(int j=0;j<pred[i]->exons.Count();j++) fprintf(stderr," %d-%d",pred[i]->exons[j].start,pred[i]->exons[j].end);
    		fprintf(stderr,"\n");
    	}
    }
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

CReadAln *guide_to_read(GffObj *t, int g, GList<CJunction>& junction, int refend) {

	char s=0;
	if(t->strand=='-') s=-1;
	else if(t->strand=='+') s=1;
	else s=0;

	TAlnInfo *tif=new TAlnInfo();
	tif->g=g;

	CReadAln *r=new CReadAln(s,1,t->start,t->end,tif);
	r->read_count=0;

	// add junction from start here --> needs to be specially handled in build_merge
	CJunction* njsource=add_junction(0, t->start, junction, r->strand);
	if (njsource) {
		r->juncs.Add(njsource);
		njsource->guide_match=1;
	}

	for(int i=0;i<t->exons.Count();i++) {
		r->segs.Add(t->exons[i]);
		r->len+=t->exons[i]->len();
		if(i) {
			CJunction* nj=add_junction(t->exons[i-1]->end, t->exons[i]->start, junction, r->strand);
			if (nj) {
				r->juncs.Add(nj);
				nj->guide_match=1;
			}
		}
	}

	// add junction to sink here as well
	CJunction* njsink=add_junction(t->end, refend, junction, r->strand);
	if (njsink) {
		r->juncs.Add(njsink);
		njsink->guide_match=1;
	}


	return(r);
}

int build_merge(BundleData* bdata) { // here a "read" is in fact a transcript
	int refstart = bdata->start;
	int refend = bdata->end+1;
	GList<CReadAln>& readlist = bdata->readlist;
	GList<CJunction>& junction = bdata->junction;
	GPVec<GffObj>& guides = bdata->keepguides;   // I need to insert the guides with the reads and process them separately
	GList<CPrediction>& pred = bdata->pred;
	// form groups on strands: all groups below are like this: 0 = negative strand; 1 = unknown strand; 2 = positive strand
	GPVec<CGroup> group;
	CGroup *currgroup[3]={NULL,NULL,NULL}; // current group of each type
	CGroup *startgroup[3]={NULL,NULL,NULL}; // start group of each type
	int color=0; // next color to assign
	GVec<int> merge; // remembers merged groups
	GVec<int> equalcolor; // remembers colors for the same bundle

	if(guides.Count()) {
		for(int i=0;i<guides.Count();i++) { // more guides left
			CReadAln *rg=guide_to_read(guides[i],i,junction,refend);
			readlist.Add(rg);
		}
		readlist.setSorted(true); // I need to do this before adding color otherwise I am loosing the color, readgroup ordering etc.
	}

	int nreads=readlist.Count(); // initial number of transcripts before I add the guides
	GVec<int> *readgroup=new GVec<int>[nreads]; // remebers groups for each read; don't forget to delete it when no longer needed

	for (int n=0;n<nreads;n++) if(readlist[n]->nh){

		//fprintf(stderr,"Add read %d\n",n);

		CReadAln & rd=*(readlist[n]);
		/* 	// if I want to eliminate transcripts that have very bad junctions this is the code to do it

		bool keep=true;
		if(rd.tinfo->g==-1) { // if this is not a guide
			int i=0;
			while(i<rd.juncs.Count()) {
				CJunction& jd=*(rd.juncs[i]);
				if(!jd.strand || !good_merge_junc(jd,junction)) { // found a bad junction: no need of this or count_junctions procedure probably -> just colapse reads instead
					keep=false;
					rd.nh=0;
					break;
				}
				i++;
			}
		}

		if(keep) {*/


		rd.pair_idx.Add(n);

		if(rd.tinfo->g==-1) { // if not guide : only collapse non guides
			// I need to adjust read abundance
			if((bdata->covflags & IS_TPM_FLAG)!=0) rd.read_count=rd.tinfo->tpm;
			else if((bdata->covflags & IS_FPKM_FLAG)!=0) rd.read_count=rd.tinfo->fpkm;
			else if((bdata->covflags & IS_COV_FLAG)!=0) rd.read_count=rd.tinfo->cov;

			//if(rd.juncs.Count()) {
			int m=n+1;
			while(m<nreads && readlist[m]->start<=rd.segs[0].end) {
				if(readlist[m]->nh && readlist[m]->strand==rd.strand && readlist[m]->tinfo->g==-1 && readlist[m]->juncs.Count()==rd.juncs.Count()) {
					bool samejuncs=true;
					for(int i=0;i<rd.juncs.Count();i++)
						if(rd.juncs[i]!=readlist[m]->juncs[i]) {
							samejuncs=false;
							break;
						}
					if(samejuncs) {
						rd.pair_idx.Add(m);
						readlist[m]->nh=0;
						if((bdata->covflags & IS_TPM_FLAG)!=0) readlist[m]->read_count=readlist[m]->tinfo->tpm;
						else if((bdata->covflags & IS_FPKM_FLAG)!=0) readlist[m]->read_count=readlist[m]->tinfo->fpkm;
						else if((bdata->covflags & IS_COV_FLAG)!=0) readlist[m]->read_count=readlist[m]->tinfo->cov;
						rd.read_count=rd.read_count*rd.len+readlist[m]->read_count*readlist[m]->len;
						if(readlist[m]->end>rd.end) {
							rd.len+=readlist[m]->end-rd.end;
							rd.end=readlist[m]->end;
							rd.segs.Last().end=rd.end;
						}
						rd.read_count/=rd.len;
					}
				}
				m++;
			}
		}

		// I need to add the source to start, and end to sink junctions to read as well
		char strand=-1;
		if(rd.strand) strand=rd.strand;
		CJunction* nj=add_junction(0, rd.start, junction, strand);
		if (nj) rd.juncs.Insert(0,nj);
		nj=NULL;
		nj=add_junction(rd.end, refend, junction, strand);
		if (nj) rd.juncs.Add(nj);
		if(!rd.strand) { // add junction to the other strand as well
			strand=1;
			nj=add_junction(0, rd.start, junction, strand);
			//if (nj) rd.juncs.Insert(0,nj); I don't do this anymore because I don't know the strand
			nj=add_junction(rd.end, refend, junction, strand);
			//if (nj) rd.juncs.Add(nj);
		}

		/*
			if(!rd.read_count) {
				rd.read_count=0.000001; // if fpkm is zero I don't want to ignore the transcript completely
				GError("rd[%d].read_count=%f name=%s\n",n,rd.read_count,rd.tinfo->name.chars());
			}
		*/

		color=add_read_to_group(n,readlist,color,group,currgroup,startgroup,readgroup,equalcolor,merge);
		//} //end if(keep)

		//else fprintf(stderr,"keep false for read %d\n",n);

	}

	/* // doesn't quite work yet from the way it's designed
	// merge groups that are close together or groups that are within the same exon of a reference gene
	if(bundledist) {
		for(int sno=0;sno<3;sno++) {
			CGroup *lastgroup=NULL;
			CGroup *procgroup=startgroup[sno];
			while(procgroup!=NULL) {

				if(lastgroup) {

					GStr bstart((int)lastgroup->end);
			    	GStr bend((int)procgroup->start);
			    	if(procgroup->start-lastgroup->end<=bundledist) {

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
	*/

	// ### form bundles here

    // first do color assignment
	for (int i=0;i<3;i++) currgroup[i]=startgroup[i];
	CGroup *prevgroup[3]={NULL,NULL,NULL};

	GVec<int> eqposcol(true);
	GVec<int> eqnegcol(true);

	eqposcol.Resize(equalcolor.Count(),-1);
	eqnegcol.Resize(equalcolor.Count(),-1);


	// each unstranded group needs to remember what proportion of stranded group it overlaps so that it can distribute reads later on -> maybe I can do this in the following while?
	while(currgroup[0]!=NULL || currgroup[1]!=NULL || currgroup[2]!=NULL) { // there are still groups to process

		int nextgr=get_min_start(currgroup); // gets the index of currgroup with the left most begining

		int grcol = currgroup[nextgr]->color;    // set smallest color for currgroup[$nextgr]

		while(equalcolor[grcol]!=grcol) {
			grcol=equalcolor[grcol];
		}
		currgroup[nextgr]->color=grcol;


		if(nextgr == 1) { // unknown strand group

			//if(prevgroup[0]!=NULL && currgroup[nextgr]->start <= prevgroup[0]->end+bundledist) { // overlaps previous negative group
			if(prevgroup[0]!=NULL && currgroup[nextgr]->start <= prevgroup[0]->end) { // overlaps previous negative group

				set_strandcol(currgroup[nextgr],prevgroup[0],prevgroup[0]->color,eqnegcol,equalcolor);
				uint maxstart = currgroup[nextgr]->start > prevgroup[0]->start ? currgroup[nextgr]->start : prevgroup[0]->start;
				uint minend = currgroup[nextgr]->end < prevgroup[0]->end ? currgroup[nextgr]->end : prevgroup[0]->end;
				currgroup[nextgr]->neg_prop+=prevgroup[0]->cov_sum*(minend-maxstart+1)/prevgroup[0]->len();
			}

			//while(currgroup[0]!=NULL && currgroup[nextgr]->start <= currgroup[0]->end+bundledist && currgroup[0]->start <= currgroup[nextgr]->end+bundledist) { // overlaps current negative strand group
			while(currgroup[0]!=NULL && currgroup[nextgr]->start <= currgroup[0]->end && currgroup[0]->start <= currgroup[nextgr]->end) { // overlaps current negative strand group

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
			//if(prevgroup[2]!=NULL && currgroup[nextgr]->start <= prevgroup[2]->end + bundledist) { // overlaps positive strand group
			if(prevgroup[2]!=NULL && currgroup[nextgr]->start <= prevgroup[2]->end) { // overlaps positive strand group
				set_strandcol(currgroup[nextgr],prevgroup[2],prevgroup[2]->color,eqposcol,equalcolor);
				if(currgroup[nextgr]->neg_prop) {
					uint maxstart = currgroup[nextgr]->start > prevgroup[2]->start ? currgroup[nextgr]->start : prevgroup[2]->start;
					uint minend = currgroup[nextgr]->end < prevgroup[2]->end ? currgroup[nextgr]->end : prevgroup[2]->end;
					pos_prop+=prevgroup[2]->cov_sum*(minend-maxstart+1)/prevgroup[2]->len();
				}
			}

			//while(currgroup[2]!=NULL && currgroup[nextgr]->start <= currgroup[2]->end + bundledist && currgroup[2]->start <= currgroup[nextgr]->end + bundledist) { // overlaps positive strand group
			while(currgroup[2]!=NULL && currgroup[nextgr]->start <= currgroup[2]->end && currgroup[2]->start <= currgroup[nextgr]->end) { // overlaps positive strand group

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
				add_group_to_bundle(currgroup[nextgr],bundle[nextgr][bno],bnode[nextgr],0); // last parameter was bundledist before
			}
			else { // create new bundle
				bno=create_bundle(bundle[nextgr],currgroup[nextgr],bnode[nextgr]);
				bundlecol[grcol]=bno;
			}

			group2bundle[nextgr][currgroup[nextgr]->grid]=bundle[nextgr][bno]->lastnodeid;

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
					add_group_to_bundle(currgroup[nextgr],bundle[0][bno],bnode[0],0); // last parameter was bundledist before
				}
				else { // create new bundle
					bno=create_bundle(bundle[0],currgroup[nextgr],bnode[0]);
					//print STDERR "Added neutral group ",$currgroup[$nextgr][3]," to new neg bundle $bno\n";
					bundlecol[negcol]=bno;
				}
				group2bundle[0][currgroup[nextgr]->grid]=bundle[0][bno]->lastnodeid;
			} // if(eqnegcol[grcol]!=-1)

			if(eqposcol[grcol]!=-1){
				int poscol=eqposcol[grcol];
				while(equalcolor[poscol]!=poscol) {
					poscol=equalcolor[poscol];
				}

				int bno=bundlecol[poscol];
				if(bno>-1) { // bundle for group has been created before
					//print STDERR "Add neutral group ",$currgroup[$nextgr][3]," to pos bundle $bno\n";
					add_group_to_bundle(currgroup[nextgr],bundle[2][bno],bnode[2],0); // last parameter was bundledist before
				}
				else { // create new bundle
					bno=create_bundle(bundle[2],currgroup[nextgr],bnode[2]);
					//print STDERR "Added neutral group ",$currgroup[$nextgr][3]," to new pos bundle $bno\n";
					bundlecol[poscol]=bno;
				}
				group2bundle[2][currgroup[nextgr]->grid]=bundle[2][bno]->lastnodeid;
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


    // ### store transcripts in unstranded bundles here; I can not have guides in here because those are signed
	for(int b=0;b<bundle[1].Count();b++) {

    	if(bundle[1][b]->nread) {
			char sign='.';

    		// bundle shouldn't contain multiple transcripts since there is no pairing inside merging
    		CBundlenode *currbnode=bnode[1][bundle[1][b]->startnode];
    		int t=1;
    		while(currbnode!=NULL) {
    			if(t==1) { geneno++;}

    			//fprintf(stderr,"Store unstranded prediction: geneno=%d start=%d end=%d cov=%f\n",geneno,currbnode->start,currbnode->end,currbnode->cov);
    			int len=currbnode->end-currbnode->start+1;
    			CPrediction *p=new CPrediction(geneno-1, NULL, currbnode->start, currbnode->end, currbnode->cov/len, sign, len);
    			GSeg exon(currbnode->start,currbnode->end);
    			p->exons.Add(exon);
    			/* not clear how to store retrieve info without wasting memory and time
    			if(enableNames) {
    				p->mergename
    			}
    			*/
    			pred.Add(p);
    			t++;
    			currbnode=currbnode->nextnode;
    		}
    	}
    }

    //fprintf(stderr,"Done with unstranded bundles geneno=%d\n",geneno);

	// ### build graphs for stranded bundles here
    if(startgroup[0]!=NULL || startgroup[2]!=NULL) { //# there are stranded groups to process

    	// sort junctions -> junctions are sorted already according with their start, but not their end
    	GList<CJunction> ejunction(junction);
    	ejunction.setFreeItem(false);
    	if(ejunction.Count()) ejunction.setSorted(juncCmpEnd);

    	GVec<CGraphinfo> *bundle2graph[2]; // should I keep the neutral strand for consistency ? -> remember not to delete it
    	GVec<int> graphno[2];  // how many nodes are in a certain graph g, on strand s: graphno[s][g]
    	GVec<int> edgeno[2];  // how many edges are in a certain graph g, on strand s: edgeno[s][g]
    	GPVec<CTransfrag> *transfrag[2]; // for each transfrag t on a strand s, in a graph g, transfrag[s][g][t] gives it's abundance and it's pattern
    	GPVec<CMTransfrag> *mgt[2]; // merged super-transfrags
    	GPVec<CGraphnode> *no2gnode[2]; // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i
    	CTreePat **tr2no[2]; // for each graph g, on a strand s, tr2no[s][g] keeps the tree pattern structure for quick retrieval of the index t of a tansfrag
    	GIntHash<int> *gpos[2]; // for each graph g, on a strand s, gpos[s][g] keeps the hash between edges and positions in the bitvec associated to a pattern
    	GVec<int> lastgpos[2];

    	int bno[2]={0,0};


    	// build graph structure
    	for(int sno=0;sno<3;sno+=2) { // skip neutral bundles -> those shouldn't have junctions

    		int s=sno/2; // adjusted strand due to ignoring neutral strand
    		//char strnd='-';
    		//if(s) strnd='+';

    		bundle2graph[s]=NULL;
    		if(bnode[sno].Count()) bundle2graph[s]=new GVec<CGraphinfo>[bnode[sno].Count()];
    		transfrag[s]=NULL;
    		mgt[s]=NULL;
    		no2gnode[s]=NULL;
    		tr2no[s]=NULL;
    		gpos[s]=NULL;

    		if(bundle[sno].Count()) {
    			transfrag[s]=new GPVec<CTransfrag>[bundle[sno].Count()]; // for each bundle I have a graph ? only if I don't ignore the short bundles
    			mgt[s]=new GPVec<CMTransfrag>[bundle[sno].Count()];
    			no2gnode[s]=new GPVec<CGraphnode>[bundle[sno].Count()];

    			gpos[s]=new GIntHash<int>[bundle[sno].Count()];


    			GCALLOC(tr2no[s],bundle[sno].Count()*sizeof(CTreePat *));
    			bno[s]=bundle[sno].Count();

    			for(int b=0;b<bundle[sno].Count();b++) {
    				graphno[s].cAdd(0);
    				edgeno[s].cAdd(0);
    				lastgpos[s].cAdd(0);

    				/*
    				{ // DEBUG ONLY
    				if(bundle[sno][b]->nread) {
    					fprintf(stderr,"proc bundle[%d][%d] %f/%f is %f len=%d\n",sno,b,bundle[sno][b]->multi,bundle[sno][b]->nread,(float)bundle[sno][b]->multi/bundle[sno][b]->nread,bundle[sno][b]->len);
    				} }
    				*/

    				// create graph

    				GArray<GEdge> unused;
    				graphno[s][b]=create_graph(refstart,s,b,bundle[sno][b],bnode[sno],junction,ejunction,
    						bundle2graph,no2gnode,transfrag,gpos,NULL,edgeno[s][b],lastgpos[s][b],unused,refend);

    				if(graphno[s][b]) tr2no[s][b]=construct_treepat(graphno[s][b],gpos[s][b],transfrag[s][b]);
    				else tr2no[s][b]=NULL;

    				for(int i=0; i<transfrag[s][b].Count();i++) {

    					CMTransfrag *tr=new CMTransfrag(transfrag[s][b][i]); // transfrags that are created in the graph process can not be associated with any read
        				mgt[s][b].Add(tr);

    				}

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


		// I can clean up some data here:
    	for(int sno=0;sno<3;sno++) {
    		int n=bnode[sno].Count();
    		for(int b=0;b<n;b++) delete bnode[sno][b];
    		bnode[sno].Clear();
    		bundle[sno].Clear();
    	}


    	for (int n=0;n<readlist.Count();n++) { // check if read is active: nh>0

    		/*
    		{ // DEBUG ONLY
    			fprintf(stderr,"Read[%d]:",n);
    			for(int i=0;i<readlist[n]->segs.Count();i++) {
    				fprintf(stderr," %d-%d",readlist[n]->segs[i].start,readlist[n]->segs[i].end);
    			}
    			fprintf(stderr," gets tr: ");
    		}
    		*/

    		if(readlist[n]->nh) get_read_to_transfrag(readlist,n,readgroup,merge,group2bundle,bundle2graph,graphno,edgeno,gpos,no2gnode,transfrag,mgt,tr2no,group);

    		/*
    		{ // DEBUG ONLY
    			fprintf(stderr,"\n");
    		}
    		*/

    	}


    	// don't forget to clean up the allocated data here
    	delete [] readgroup;
    	group.Clear();

    	// parse graph
    	for(int s=0;s<2;s++) {

    		for(int b=0;b<bno[s];b++) {
    			if(graphno[s][b]) {

    				//process transfrags to eliminate noise, and set compatibilities, and node memberships
    				/*
    				GBitVec compatible((1+transfrag[s][b].Count())*transfrag[s][b].Count()/2); // I might want to change this to gbitvec
    				process_merge_transfrags(graphno[s][b],no2gnode[s][b],mgt[s][b],compatible,gpos[s][b]);
    				geneno=merge_transcripts(graphno[s][b],no2gnode[s][b],mgt[s][b],
    						geneno,s,pred,readlist,guides);
    				*/

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
    					mgt[s][b].Sort(mgtrabundCmp); // sort transfrags from the most abundant to the least, with guides coming in first
    					fprintf(stderr,"There are %d transfrags[%d][%d]:\n",mgt[s][b].Count(),s,b);
    					for(int t=0;t<mgt[s][b].Count();t++) {
    						fprintf(stderr,"%d (%d-%d): ",t,no2gnode[s][b][mgt[s][b][t]->transfrag->nodes[0]]->start,no2gnode[s][b][mgt[s][b][t]->transfrag->nodes.Last()]->end);
    						//printBitVec(transfrag[s][b][t]->pattern);
    						fprintf(stderr," %f:",mgt[s][b][t]->transfrag->abundance);
    						for(int i=0;i<mgt[s][b][t]->transfrag->nodes.Count();i++) fprintf(stderr," %d",mgt[s][b][t]->transfrag->nodes[i]);
    						fprintf(stderr,"\n");
    					}
    				}
    				*/

    				/*
    				// reads already have length computed so I could already get from the get_read_to_transfrag instead
    				// for each mtransfrag compute its length
    				for(int t=0;t<mgt[s][b].Count();t++) {
    					for(int i=0;i<mgt[s][b][t]->transfrag->nodes.Count();i++) {
    						CGraphnode *inode=no2gnode[s][b][mgt[s][b][t]->transfrag->nodes[i]];
    						mgt[s][b][t]->len+=inode->end-inode->start+1;
    					}
    				}
					*/

    				geneno=merge_transfrags(graphno[s][b],no2gnode[s][b], mgt[s][b],gpos[s][b],geneno,s,pred,readlist,guides);
    				//geneno=merge_transfrags_EM(graphno[s][b],no2gnode[s][b], mgt[s][b],gpos[s][b],geneno,s,pred,readlist,guides);

    			}
    			// clean up what can be cleaned
    			if(tr2no[s][b]) free_treepat(tr2no[s][b]);
    		}

    		// final clean up: no2gnode, no2tr, transfrag, bundle2graph
    		if(bundle2graph[s]) delete [] bundle2graph[s];
    		if(transfrag[s]) delete [] transfrag[s];
    		if(mgt[s]) delete [] mgt[s];
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


void count_good_junctions(GList<CReadAln>& readlist, int refstart, GVec<float>* bpcov) {

	for(int n=0;n<readlist.Count();n++) {
		CReadAln & rd=*(readlist[n]);
		int nex=rd.segs.Count();
		GVec<uint> leftsup;
		GVec<uint> rightsup;
		uint maxleftsupport=0;
		uint maxrightsupport=0;
		int sno=(int)rd.strand+1; // 0(-),1(.),2(+)
		//if(nex>1) fprintf(stderr,"Process spliced read[%d]: ",n);
		for(int i=0;i<nex;i++) {
			//if(nex>1) fprintf(stderr," %d-%d",rd.segs[i].start,rd.segs[i].end);
			if(i) {
				if(rd.segs[i-1].len()>maxleftsupport) maxleftsupport=rd.segs[i-1].len();
				if(rd.segs[nex-i].len()>maxrightsupport) maxrightsupport=rd.segs[nex-i].len();
				leftsup.Add(maxleftsupport);
				rightsup.Add(maxrightsupport);
				rd.juncs[i-1]->nreads+=rd.read_count;
			}
			cov_add(bpcov,sno,rd.segs[i].start-refstart,rd.segs[i].end-refstart,rd.read_count);
		}
		//if(nex>1) fprintf(stderr," With anchors: ");
		for(int i=1;i<nex;i++) {
			uint anchor=junctionsupport;
			if(anchor<longintronanchor && rd.juncs[i-1]->len()>longintron) anchor=longintronanchor; // I want to use a longer anchor for long introns to believe them
			//if(leftsup[i-1]>=anchor && rightsup[nex-i-1]>=anchor) rd.juncs[i-1]->nreads_good+=rd.read_count;
			if(leftsup[i-1]>=anchor) {
				rd.juncs[i-1]->leftsupport+=rd.read_count;
				if(rightsup[nex-i-1]>=anchor) {
					rd.juncs[i-1]->rightsupport+=rd.read_count;
					rd.juncs[i-1]->nreads_good+=rd.read_count;
				}
			}
			else if(rightsup[nex-i-1]>=anchor) {
				rd.juncs[i-1]->rightsupport+=rd.read_count;
			}
			//fprintf(stderr," %d(%d-%d)",anchor,leftsup[i-1],rightsup[nex-i-1]);
		}
		//if(nex>1) fprintf(stderr,"\n");
	}

}

void count_merge_junctions(GList<CReadAln>& readlist,char covflags) { // why do I need to impose any thresholds here??

	for(int n=0;n<readlist.Count();n++) {
		CReadAln & rd=*(readlist[n]);
		int nex=rd.segs.Count();
		for(int i=0;i<nex;i++) {
			if(i) {
				rd.juncs[i-1]->nreads+=rd.read_count;

				if((covflags & IS_TPM_FLAG)!=0) rd.juncs[i-1]->nreads_good+=rd.tinfo->tpm;
				else if((covflags & IS_FPKM_FLAG)!=0) rd.juncs[i-1]->nreads_good+=rd.tinfo->fpkm;
				else if((covflags & IS_COV_FLAG)!=0) rd.juncs[i-1]->nreads_good+=rd.tinfo->cov;
				else rd.juncs[i-1]->nreads_good+=rd.read_count;
			}

		}
	}

}


//int infer_transcripts(int refstart, GList<CReadAln>& readlist,
		//GList<CJunction>& junction, GPVec<GffObj>& guides, GVec<float>& bpcov, GList<CPrediction>& pred, bool fast) {
int infer_transcripts(BundleData* bundle) {
	int geneno=0;

	//DEBUG ONLY: 	showReads(refname, readlist);

/*
#ifdef GMEMTRACE
	double vm,rsm;
	get_mem_usage(vm, rsm);
	GMessage("\t\tM(s):infer_transcripts memory usage: rsm=%6.1fMB vm=%6.1fMB\n",rsm/1024,vm/1024);
#endif
*/
	if(mergeMode) {
		//count_merge_junctions(bundle->readlist,bundle->covflags); // this make sense only if I want to count junctions
		geneno = build_merge(bundle);
	}
	else if(bundle->keepguides.Count() || !eonly) {
		//fprintf(stderr,"Process %d reads from %lu.\n",bundle->readlist.Count(),bundle->numreads);

		count_good_junctions(bundle->readlist, bundle->start, bundle->bpcov);

		geneno = build_graphs(bundle);
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

	if(a->start < b->start) return -1; // order based on start
	if(a->start > b->start) return 1;

	if(a->cov > b->cov) return -1; // the one with higher coverage comes first
	if(a->cov < b->cov) return 1;

	// same start

	//if(sensitivitylevel!=1 || sensitivitylevel!=2) { // try to see if I correct this if it makes any difference (it should be && instead of || in the if)
	if(a->exons.Count() < b->exons.Count()) return -1; // order based on number of exons
	if(a->exons.Count() > b->exons.Count()) return 1;

	// same number of exons
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
		if(a1<b1) return -1; // the one with the first exon comes first
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

int predcovCmp(const pointer p1, const pointer p2) { // sort from highest to lowest coverage
	CPrediction *a=(CPrediction*)p1;
	CPrediction *b=(CPrediction*)p2;
	if(a->cov < b->cov) return 1;
	if(a->cov > b->cov) return -1;
	return 0;
}

int predlenCmp(const pointer p1, const pointer p2) { // sort from highest to lowest coverage
	CPrediction *a=(CPrediction*)p1;
	CPrediction *b=(CPrediction*)p2;
	if(a->exons.Count() < b->exons.Count()) return 1;
	if(a->exons.Count() > b->exons.Count()) return -1;
	if(a->t_eq==NULL && b->t_eq!=NULL) return 1;
	if(a->t_eq!=NULL && b->t_eq==NULL) return -1;
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
*/

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

	if(pred[big]->cov<pred[small]->cov && !pred[big]->t_eq && !bex && pred[small]->exons.Count()>1 && pred[big]->exons[0].start<pred[small]->exons[0].start) { // adjust start to account for trimming
		pred[big]->tlen-=pred[small]->exons[0].start-pred[big]->exons[0].start;
		pred[big]->exons[0].start=pred[small]->exons[0].start;
		pred[big]->start=pred[small]->start;
	}

	int sex=0;
	int overlap=0;
	while(sex<pred[small]->exons.Count()) {
		int exovlp=(pred[small]->exons[sex].end<pred[big]->exons[bex].end ? pred[small]->exons[sex].end : pred[big]->exons[bex].end)-
				(pred[small]->exons[sex].start>pred[big]->exons[bex].start ? pred[small]->exons[sex].start : pred[big]->exons[bex].start)+1;

		if(pred[big]->cov<pred[small]->cov && !pred[big]->t_eq && bex==pred[big]->exons.Count()-1 && sex>=1 && pred[big]->exons[bex].end>pred[small]->exons[sex].end) { // adjust end
			pred[big]->tlen-=pred[big]->exons[bex].end-pred[small]->exons[sex].end;
			pred[big]->exons[bex].end=pred[small]->exons[sex].end;
			pred[big]->end=pred[small]->end;
		}

		pred[big]->exoncov[bex]=(pred[big]->exoncov[bex]*pred[big]->exons[bex].len()+frac*pred[small]->exoncov[sex]*exovlp)/pred[big]->exons[bex].len();
		overlap+=exovlp;
		sex++;bex++;
	}

	pred[big]->cov=(pred[big]->tlen*pred[big]->cov+overlap*frac*pred[small]->cov)/pred[big]->tlen;

	//fprintf(stderr,"New prediction[%d] cov=%f\n",big,pred[big]->cov);

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

/* *** this is the printing function used in ver5sion 1.2.3
// default printing function for sensitivitylevel=1: all the others are deprecated for now
int print_cluster_prev(GPVec<CPrediction>& pred,GVec<int>& genes,GVec<int>& transcripts, int geneno,GStr& refname,
		GVec<CGene>& refgene, GHash<int>& hashgene, GVec<CGene>& predgene, int startgno) {

	// sort predictions from the most abundant to the least:
	pred.Sort(predcovCmp);
	GVec<int> keep;

	CInterval *maxpos=NULL; //remembers intervals of maximum coverage

	for(int n=0;n<pred.Count();n++) {


		{ // DEBUG ONLY
			fprintf(stderr,"Consider prediction[%d] %c cov=%f:",n,pred[n]->strand,pred[n]->cov);
			for(int i=0;i<pred[n]->exons.Count();i++) fprintf(stderr," %d-%d",pred[n]->exons[i].start,pred[n]->exons[i].end);
			if(pred[n]->t_eq) fprintf(stderr," ref_id=%s",pred[n]->t_eq->getID());
			fprintf(stderr,"\n");
		}


		int k=0;
		bool included=false;
		while(!included && k<keep.Count() && keep[k]<n) {

			if(included_pred(pred,keep[k],n)) {

				fprintf(stderr,"%d prediction is included in prediction %d\n",keep[k],n);

				bool checkall=false;

				if(pred[keep[k]]->exons.Count()<pred[n]->exons.Count()) {
					//if(pred[keep[k]]->cov>pred[n]->cov) break; // this is new and improves specificity but I loose some things -> TO CHECK WHAT IT ACT: also this should always happen because of the sort procedure
					//if(pred[keep[k]]->exons.Count()>2) break; // this is what I had before but I don't think it makes sense completely so I introduce the next one for those cases where some single exons are still included in here
					// *** if(pred[keep[k]]->exons.Count()>2) { k++; break;} // I need to check this how it affects performance in general
					if(pred[keep[k]]->exons.Count()>2) { k++; if(pred[n]->t_eq) continue; else break;} // I need to check this how it affects performance in general
					update_cov(pred,n,keep[k]);
					// I might also update maxpos here when I engulf a smaller prediction otherwise I could lose them both

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


				fprintf(stderr,"...included in prediction[%d] with cov=%f\n",keep[k],pred[keep[k]]->cov);
				included=true;
				break; // if it's included than I am done with the while loop because n got used
			}
			k++;
		}
		if(included) continue; / *  I am doing maxpos after instead of here {
			maxpos=add_pred_to_cov(maxpos,pred[keep[k]]);


			{ // DEBUG ONLY
				fprintf(stderr,"Maxpos is:");
				CInterval *interval=maxpos;
				while(interval!=NULL) {
					fprintf(stderr," pos=%d val=%f",interval->pos,interval->val);
					interval=interval->next;
				}
				fprintf(stderr,"\n");
			}


			continue;
		}

		bool abundant=true;
		maxpos=add_pred_to_cov(maxpos,pred[n],&abundant);
		//maxpos=add_pred_to_cov(maxpos,pred[n]);


		{ // DEBUG ONLY
			fprintf(stderr,"Maxpos is:");
			CInterval *interval=maxpos;
			while(interval!=NULL) {
				fprintf(stderr," pos=%d val=%f",interval->pos,interval->val);
				interval=interval->next;
			}
			fprintf(stderr," abundant=%d\n",abundant);
		}


		if(pred[n]->t_eq || abundant) {
		//if(pred[n]->id || abundant) {
			keep.Add(n);
			fprintf(stderr,"...keep prediction %d\n",n);
		}* /
		keep.Add(n);

	}
	// the prediction might not be sorted by coverage by this point because of the joining but keep refers to he initial ordering
	for(int i=0;i<keep.Count();i++) {
		maxpos=add_pred_to_cov(maxpos,pred[keep[i]]);
		fprintf(stderr,"...keep prediction %d with cov=%f\n",keep[i],pred[keep[i]]->cov);
	}

	for(int i=0;i<keep.Count();i++) {
	  int n=keep[i];
	  if(pred[n]->t_eq || (!eonly && is_pred_above_frac(maxpos,pred[n]))) { // print this transcript only if it's matching known gene or if not eonly and above threshold

		  / *
		  { // DEBUG ONLY
			  fprintf(stderr,"print prediction %d",n);
			  if(pred[n]->flag) fprintf(stderr," with true flag");
			  fprintf(stderr,"\n");
		  }
		  * /

		  if(pred[n]->flag) {
			  if(genes[pred[n]->geneno]==-1) genes[pred[n]->geneno]=++geneno;
			  transcripts[pred[n]->geneno]++;
			  uint t_id=0;
			  if (pred[n]->t_eq && pred[n]->t_eq->uptr) {
				  t_id = ((RC_TData*)pred[n]->t_eq->uptr)->t_id;
			  }
			  //fprintf(f_out,"%d %d %d %.6f %.6f\n",pred[n]->exons.Count()+1,pred[n]->tlen, t_id, pred[n]->frag,pred[n]->cov);
			  fprintf(f_out,"1 %d %d %d %.6f\n",pred[n]->exons.Count()+1,pred[n]->tlen, t_id,pred[n]->cov);

			  GStr geneid(label);geneid+='.';geneid+=genes[pred[n]->geneno];
			  GStr trid(geneid); trid+='.';trid+=transcripts[pred[n]->geneno];
			  if(eonly && pred[n]->t_eq && pred[n]->t_eq->getGeneID()) geneid=pred[n]->t_eq->getGeneID();
			  if(eonly && pred[n]->t_eq) trid=pred[n]->t_eq->getID();

			  fprintf(f_out,"%s\tStringTie\ttranscript\t%d\t%d\t1000\t%c\t.\tgene_id \"%s\"; transcript_id \"%s\"; ",
						  refname.chars(),pred[n]->start,pred[n]->end,pred[n]->strand,geneid.chars(),trid.chars());

			  if(pred[n]->t_eq) {
				  if(!eonly) fprintf(f_out,"reference_id \"%s\"; ",pred[n]->t_eq->getID());
				  if (!eonly && pred[n]->t_eq->getGeneID())
				    	 fprintf(f_out,"ref_gene_id \"%s\"; ",pred[n]->t_eq->getGeneID());
				  if (pred[n]->t_eq->getGeneName())
				    	 fprintf(f_out,"ref_gene_name \"%s\"; ",pred[n]->t_eq->getGeneName());
			  }
			  fprintf(f_out,"cov \"%.6f\";\n",pred[n]->cov);
			  for(int j=0;j<pred[n]->exons.Count();j++) {
				  fprintf(f_out,"%s\tStringTie\texon\t%d\t%d\t1000\t%c\t.\tgene_id \"%s\"; transcript_id \"%s\"; exon_number \"%d\"; ",
						  refname.chars(),pred[n]->exons[j].start,pred[n]->exons[j].end,pred[n]->strand,geneid.chars(),
						  trid.chars(),j+1); // maybe add exon coverage here
				  if(pred[n]->t_eq) {
					  if(!eonly) fprintf(f_out,"reference_id \"%s\"; ",pred[n]->t_eq->getID());
					  if (!eonly && pred[n]->t_eq->getGeneID())
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
				  if(gid.is_empty()) gid=pred[n]->t_eq->getGeneName();
				  if(!gid.is_empty()) {
					  gid+=pred[n]->strand;
					  const int *ng=hashgene[gid.chars()];
					  if(ng) { // this should always be true because we parsed all predictions in printResults
						  gno=*ng;
						  / * //I don't need to do this because I already did it in printResults
					  	  if(pred[n]->start<refgene[gno].start) refgene[gno].start=pred[n]->start;
					  	  if(pred[n]->end>predgene[gno].end) predgene[gno].end=pred[n]->end;
					  	  merge_exons(predgene[gno],pred[n]->exons);
						   * /
						  refgene[gno].cov+=pred[n]->cov*pred[n]->tlen;
						  refgene[gno].covsum+=pred[n]->cov;

					  }
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
*/


// need to work on this
int print_cluster(GPVec<CPrediction>& pred,GVec<int>& genes,GVec<int>& transcripts, int geneno,GStr& refname,
		GVec<CGene>& refgene, GHash<int>& hashgene, GVec<CGene>& predgene, int startgno,double& bigsumcov) {

	// sort predictions from the most abundant to the least:
	pred.Sort(predlenCmp);
	GVec<int> keep;   // only keep predictions that are not included in other ones

	for(int n=0;n<pred.Count();n++) {

		/*
		{ // DEBUG ONLY
			fprintf(stderr,"Consider prediction[%d] %c cov=%f len=%d:",n,pred[n]->strand,pred[n]->cov,pred[n]->tlen);
			for(int i=0;i<pred[n]->exons.Count();i++) fprintf(stderr," %d-%d",pred[n]->exons[i].start,pred[n]->exons[i].end);
			if(pred[n]->t_eq) fprintf(stderr," ref_id=%s",pred[n]->t_eq->getID());
			fprintf(stderr,"\n");
		}
		*/

		int k=0;

		if(pred[n]->t_eq==NULL) { // if this is a known transcript I have to keep it
			GVec<int> inc;
			float sumcov=0;
			bool included=false;
			while(k<keep.Count() && keep[k]<n) {

				if(included_pred(pred,keep[k],n)) { // prediction n has less exons, or is less abundant

					//fprintf(stderr,"Prediction %d is included in prediction %d\n",n,keep[k]);
					inc.Add(keep[k]);
					sumcov+=pred[keep[k]]->cov;
					if(pred[n]->exons.Count()<=2 || pred[n]->cov<pred[keep[k]]->cov
							|| pred[keep[k]]->t_eq || pred[n]->exons.Count()==pred[keep[k]]->exons.Count())
							included=true;
				}

				k++;
			}

			if(included && sumcov) { // prediction n gets included in a bunch of other predictions
				for(int i=0;i<inc.Count();i++) {
					float frac=pred[inc[i]]->cov/sumcov;
					update_cov(pred,inc[i],n,frac);
				}
				continue;
			}
		}

		keep.Add(n);
	}

	CInterval *maxpos=NULL; //remembers intervals of maximum coverage
	// the prediction might not be sorted by coverage by this point because of the joining but keep refers to he initial ordering
	for(int i=0;i<keep.Count();i++) {
		maxpos=add_pred_to_cov(maxpos,pred[keep[i]]);
		//fprintf(stderr,"...keep prediction %d with cov=%f\n",keep[i],pred[keep[i]]->cov);
	}

	for(int i=0;i<keep.Count();i++) {
	  int n=keep[i];
	  if(pred[n]->t_eq || (!eonly && is_pred_above_frac(maxpos,pred[n]))) { // print this transcript only if it's matching known gene or if not eonly and above threshold

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

			  GStr geneid(label);geneid+='.';geneid+=genes[pred[n]->geneno];
			  GStr trid(geneid); trid+='.';trid+=transcripts[pred[n]->geneno];
			  if(eonly && pred[n]->t_eq && pred[n]->t_eq->getGeneID()) geneid=pred[n]->t_eq->getGeneID();
			  if(eonly && pred[n]->t_eq) trid=pred[n]->t_eq->getID();

			  fprintf(f_out,"%s\tStringTie\ttranscript\t%d\t%d\t1000\t%c\t.\tgene_id \"%s\"; transcript_id \"%s\"; ",
						  refname.chars(),pred[n]->start,pred[n]->end,pred[n]->strand,geneid.chars(),trid.chars());

			  if(pred[n]->t_eq) {
				  if(!eonly) fprintf(f_out,"reference_id \"%s\"; ",pred[n]->t_eq->getID());
				  if (!eonly && pred[n]->t_eq->getGeneID())
				    	 fprintf(f_out,"ref_gene_id \"%s\"; ",pred[n]->t_eq->getGeneID());
				  if (pred[n]->t_eq->getGeneName())
				    	 fprintf(f_out,"ref_gene_name \"%s\"; ",pred[n]->t_eq->getGeneName());
			  }
			  fprintf(f_out,"cov \"%.6f\";\n",pred[n]->cov);
			  for(int j=0;j<pred[n]->exons.Count();j++) {
				  fprintf(f_out,"%s\tStringTie\texon\t%d\t%d\t1000\t%c\t.\tgene_id \"%s\"; transcript_id \"%s\"; exon_number \"%d\"; ",
						  refname.chars(),pred[n]->exons[j].start,pred[n]->exons[j].end,pred[n]->strand,geneid.chars(),
						  trid.chars(),j+1); // maybe add exon coverage here
				  if(pred[n]->t_eq) {
					  if(!eonly) fprintf(f_out,"reference_id \"%s\"; ",pred[n]->t_eq->getID());
					  if (!eonly && pred[n]->t_eq->getGeneID())
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
					  //fprintf(stderr,"3:add pred[%d]'s coverage=%g to gene cov\n",n,pred[n]->cov);
					  predgene.Add(g);
				  }
				  else { // I've seen this gene before
					  if(pred[n]->start<predgene[gno].start) predgene[gno].start=pred[n]->start;
					  if(pred[n]->end>predgene[gno].end) predgene[gno].end=pred[n]->end;
					  merge_exons(predgene[gno],pred[n]->exons);
					  predgene[gno].cov+=pred[n]->cov*pred[n]->tlen;
					  predgene[gno].covsum+=pred[n]->cov;
					  //fprintf(stderr,"4:add pred[%d]'s coverage=%g to gene cov\n",n,pred[n]->cov);
				  }
			  }
			  // annotated
			  if(pred[n]->t_eq && (geneabundance ||eonly)) {
				  GStr gid(pred[n]->t_eq->getGeneID());
				  if(gid.is_empty()) gid=pred[n]->t_eq->getGeneName();
				  if(!gid.is_empty()) {
					  gid+=pred[n]->strand;
					  const int *ng=hashgene[gid.chars()];
					  if(ng) { // this should always be true because we parsed all predictions in printResults
						  gno=*ng;
						  refgene[gno].cov+=pred[n]->cov*pred[n]->tlen;
						  refgene[gno].covsum+=pred[n]->cov;

					  }
					  else {

						  //fprintf(stderr,"1:add pred[%d]'s coverage=%g\n",n,pred[n]->cov);

						  bigsumcov+=pred[n]->cov; // I am adding to cov_sum the reference transcript's coverage --> isn't this double counting?
					  }
				  }
				  else {
					  //fprintf(stderr,"2:add pred[%d]'s coverage=%g\n",n,pred[n]->cov);
					  bigsumcov+=pred[n]->cov;
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



/*
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


		{ // DEBUG ONLY
			fprintf(stderr,"Consider prediction[%d] %c cov=%f:",n,pred[n]->strand,pred[n]->cov);
			for(int i=0;i<pred[n]->exons.Count();i++) fprintf(stderr," %d-%d",pred[n]->exons[i].start,pred[n]->exons[i].end);
			fprintf(stderr," included in predictions:");
			for(int i=0;i<included[n].Count();i++) fprintf(stderr," %d",included[n][i]);
			if(n<pred.Count()-1) fprintf(stderr," maxcov=%f totalcov=%f",maxcov[n],totalcov[n]);
			fprintf(stderr,"\n");
		}


		if(included[n].Count() && (maxcov[n]>=pred[n]->cov || pred[n]->exons.Count()<limit)) { // this prediction is included in others, and is less abundant or has very few exons: make it <=2 if two exon genes also are to be ignored
			for(int k=0;k<included[n].Count();k++) {
				update_cov(pred,included[n][k],n,pred[included[n][k]]->cov/totalcov[n]);
			}
		}
		else {

			bool abundant=true;
			maxpos=add_pred_to_cov(maxpos,pred[n],&abundant);
			//maxpos=add_pred_to_cov(maxpos,pred[n]);


			{ // DEBUG ONLY
				fprintf(stderr,"Maxpos is:");
				CInterval *interval=maxpos;
				while(interval!=NULL) {
					fprintf(stderr," pos=%d val=%f",interval->pos,interval->val);
					interval=interval->next;
				}
				fprintf(stderr,"\n");
			}


			if(pred[n]->t_eq || abundant) {
			//if(pred[n]->id || abundant) {
				keep.Add(n);
				//fprintf(stderr,"...keep prediction\n");
			}
		}
	}

	for(int i=0;i<keep.Count();i++) {
	  int n=keep[i];
	  if(pred[n]->t_eq || (!eonly && is_pred_above_frac(maxpos,pred[n]))) { // print this transcript only if it's matching known gene or if not eonly and above threshold

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
*/

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

bool equal_strand(CPrediction *p1,CPrediction *p2) {

	if(p1->strand==p2->strand || p2->strand=='.') return true;
	if(p1->strand=='.') { p1->strand=p2->strand; return true;}

	return false;
}

int printMergeResults(BundleData* bundleData, int geneno, GStr& refname) {
	GList<CPrediction>& pred = bundleData->pred;
	int npred=pred.Count();

	// this version sorts the predictions according to their start and only eliminates predictions that are overlaped by others that are more abundant within the threshold
	GVec<bool> keep(npred,bool(true));
	pred.setSorted(predCmp);
	for(int n=0;n<npred;n++) if(keep[n]){
		//fprintf(stderr,"Evaluate prediction[%d] of cov %g %c:",n,pred[n]->cov,pred[n]->strand);
		//for(int e=0;e<pred[n]->exons.Count();e++) fprintf(stderr," %d-%d",pred[n]->exons[e].start,pred[n]->exons[e].end);
		//fprintf(stderr,"\n");
		int m=n+1;
		while(m<npred && pred[m]->start<=pred[n]->end+bundledist) { // prediction m overlaps prediction n or is very close by
			if(keep[m] && equal_strand(pred[n],pred[m])) {
				if(pred[m]->start<=pred[n]->end) { // predicions actually overlap
					if(!pred[n]->t_eq && ((pred[m]->exons.Count()>1 && pred[n]->cov<isofrac*pred[m]->cov) ||
							(pred[n]->exons.Count()==1 && pred[n]->cov<pred[m]->cov))) { // only eliminate prediction if it's not guide
						keep[n]=false;
						//fprintf(stderr,"  pred[%d] eliminated by pred[%d]\n",n,m);
						break;
					}
					if(!pred[m]->t_eq && ((pred[n]->exons.Count()>1 && pred[m]->cov<isofrac*pred[n]->cov) ||
							(pred[m]->exons.Count()==1 && pred[m]->cov<pred[n]->cov))) // don't believe an intronic single exon prediction
							keep[m]=false;

				}
				else // pred[m] is within bundledist of pred[n] -> see if we can merge them
					if(pred[n]->exons.Count()==1 || pred[m]->exons.Count()==1) {

						//fprintf(stderr," pred[%d] of len=%d eliminated pred[%d] of len=%d and cov %g %c:",n,pred[n]->tlen,m,pred[m]->tlen,pred[m]->cov,pred[m]->strand);
						//for(int e=0;e<pred[m]->exons.Count();e++) fprintf(stderr," %d-%d",pred[m]->exons[e].start,pred[m]->exons[e].end);
						//fprintf(stderr,"\n");

						if(pred[n]->exons.Count()==1 && pred[m]->exons.Count()==1 && !pred[n]->t_eq && !pred[m]->t_eq) { // merge predictions if both single exons (code also works if only one is single exon
							keep[m]=false;

							pred[n]->cov=pred[n]->cov*pred[n]->tlen+pred[m]->cov*pred[m]->tlen;
							pred[n]->tlen+=pred[m]->tlen+pred[m]->start-pred[n]->end-1;
							pred[n]->cov/=pred[n]->tlen;

							//fprintf(stderr," and got new len=%d and new cov=%g\n",pred[n]->tlen,pred[n]->cov);

							pred[n]->end=pred[m]->end;
							if(pred[m]->exons.Count()==1) {
								pred[n]->exons.Last().end=pred[n]->end;
							}
							else { // pred[m] has more than one exon => pred[n] has one exon
								pred[n]->exons[0].end=pred[m]->exons[0].end;
								for(int e=1;e<pred[m]->exons.Count();e++) pred[n]->exons.Add(pred[m]->exons[e]);
							}
						}
						else {
							if(pred[n]->exons.Count()==1 && !pred[n]->t_eq && pred[n]->cov<pred[m]->cov) { // polymerase run-on ?
								keep[n]=false;
								break;
							}
							while(m<npred && pred[m]->exons.Count()==1 && !pred[m]->t_eq && pred[m]->cov<pred[n]->cov && equal_strand(pred[n],pred[m])) { // polymerase run-on ?
								keep[m]=false;
								m++;
							}
						}

					}

			}
			m++;
		}
	}

	int n=0;
	while(n<npred) {
		if(keep[n]){
			uint currentend=pred[n]->end;
			int t=0;
			int m=n;
			while(m<npred && pred[m]->start<=currentend) {
				if(keep[m] && pred[m]->geneno ==pred[n]->geneno) { // this is a prediction worth printing
					if(!t) geneno++;
					t++;
					if(pred[m]->end>currentend) currentend=pred[m]->end;
					GStr trid;
					if(pred[m]->t_eq) trid=pred[m]->t_eq->getID();
					else {
						trid=label+'.';
						trid+=geneno;
						trid+='.';
						trid+=t;
					}

					fprintf(f_out,"%s\tStringTie\ttranscript\t%d\t%d\t1000\t%c\t.\tgene_id \"%s.%d\"; transcript_id \"%s\"; ",
						refname.chars(),pred[m]->start,pred[m]->end,pred[m]->strand,label.chars(),geneno,
						trid.chars());

					if(pred[m]->t_eq) {
						//fprintf(f_out,"reference_id \"%s\"; ",pred[m]->t_eq->getID());
						if (pred[m]->t_eq->getGeneName())
							fprintf(f_out,"gene_name \"%s\"; ",pred[m]->t_eq->getGeneName());
						if (pred[m]->t_eq->getGeneID())
							fprintf(f_out,"ref_gene_id \"%s\"; ",pred[m]->t_eq->getGeneID());
					}

					if(includecov) fprintf(f_out,"cov \"%g\"; ",pred[m]->cov);

					if(enableNames) fprintf(f_out,"input_transcripts \"%s\";\n",pred[m]->mergename.chars());
					else fprintf(f_out,"\n");

					for(int j=0;j<pred[m]->exons.Count();j++) {
						fprintf(f_out,"%s\tStringTie\texon\t%d\t%d\t1000\t%c\t.\tgene_id \"%s.%d\"; transcript_id \"%s\"; exon_number \"%d\"; ",
							refname.chars(),pred[m]->exons[j].start,pred[m]->exons[j].end,pred[m]->strand,label.chars(),geneno,
							trid.chars(),j+1);
						if(pred[m]->t_eq) {
							//fprintf(f_out,"reference_id \"%s\"; ",pred[m]->t_eq->getID());
							if (pred[m]->t_eq->getGeneName())
								fprintf(f_out,"gene_name \"%s\"; ",pred[m]->t_eq->getGeneName());
							if (pred[m]->t_eq->getGeneID())
								fprintf(f_out,"ref_gene_id \"%s\"; ",pred[m]->t_eq->getGeneID());
						}
						fprintf(f_out,"\n");
					}
					keep[m]=false; // no need to print it again
				} // end if keep[m]
				m++;
			} // end while(m<npred ...)
		} // end if keep[n]
		n++;
	}
	return(geneno);
}

/* this version will not always update the geneno correctly; plus it eliminates all transcripts with coverage under thresholds within the gene cluster:
 * e.g. if t1 overlaps t2, and t2 overlaps t3, but t3 does not overlap t1 but is larger than the isofrac than both t1, and t2 it will eliminate both t1 and t2
int printMergeResults(BundleData* bundleData, int ngenes, int geneno, GStr& refname) {

	GList<CPrediction>& pred = bundleData->pred;
	int npred=pred.Count();

	//
	GVec<float> maxabund(ngenes,float(0)); // compute maxabundance for each gene
	for(int n=0;n<npred;n++) {
		if(pred[n]->cov>maxabund[pred[n]->geneno]) {
			maxabund[pred[n]->geneno]=pred[n]->cov;
		}
	}

	// I should get the predictions in order of processing
	int t; // transcript number
	int prevgene=-1;
	for(int n=0;n<npred;n++) {
		if(pred[n]->geneno!=prevgene) {
			t=0;
			prevgene=pred[n]->geneno;
		}
		//fprintf(stderr,"maxabund=%f isofrac=%f pred[%d]->cov=%f\n",maxabund[pred[n]->geneno],isofrac,n,pred[n]->cov);
		if(pred[n]->t_eq || !maxabund[pred[n]->geneno] || pred[n]->cov/maxabund[pred[n]->geneno]>isofrac){
			if(!t) geneno++;
			t++;
			GStr trid;
			if(pred[n]->t_eq) trid=pred[n]->t_eq->getID();
			else {
				trid=label+'.';
				trid+=geneno;
				trid+='.';
				trid+=t;
			}

			fprintf(f_out,"%s\tStringTie\ttranscript\t%d\t%d\t1000\t%c\t.\tgene_id \"%s.%d\"; transcript_id \"%s\"; ",
					refname.chars(),pred[n]->start,pred[n]->end,pred[n]->strand,label.chars(),geneno,
					trid.chars());

			if(pred[n]->t_eq) {
				//fprintf(f_out,"reference_id \"%s\"; ",pred[n]->t_eq->getID());
				if (pred[n]->t_eq->getGeneName())
					fprintf(f_out,"gene_name \"%s\"; ",pred[n]->t_eq->getGeneName());
				if (pred[n]->t_eq->getGeneID())
					fprintf(f_out,"ref_gene_id \"%s\"; ",pred[n]->t_eq->getGeneID());
			}
			if(enableNames) fprintf(f_out,"input_transcripts \"%s\";\n",pred[n]->mergename.chars());
			else fprintf(f_out,"\n");
			for(int j=0;j<pred[n]->exons.Count();j++) {
				fprintf(f_out,"%s\tStringTie\texon\t%d\t%d\t1000\t%c\t.\tgene_id \"%s.%d\"; transcript_id \"%s\"; exon_number \"%d\"; ",
						refname.chars(),pred[n]->exons[j].start,pred[n]->exons[j].end,pred[n]->strand,label.chars(),geneno,
						trid.chars(),j+1);
				if(pred[n]->t_eq) {
					//fprintf(f_out,"reference_id \"%s\"; ",pred[n]->t_eq->getID());
					if (pred[n]->t_eq->getGeneName())
						fprintf(f_out,"gene_name \"%s\"; ",pred[n]->t_eq->getGeneName());
					if (pred[n]->t_eq->getGeneID())
						fprintf(f_out,"ref_gene_id \"%s\"; ",pred[n]->t_eq->getGeneID());
				}
				fprintf(f_out,"\n");
			}
		}
	}

	return(geneno);
}
*/

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

			if(eonly) { // if eonly I need to print all guides that were not printed yet
				if (guides[i]->uptr && ((RC_TData*)guides[i]->uptr)->in_bundle<3) {
					fprintf(f_out,"1 %d %d %d 0.0\n",guides[i]->exons.Count()+1,guides[i]->covlen, ((RC_TData*)guides[i]->uptr)->t_id);
					fprintf(f_out, "%s\t%s\ttranscript\t%d\t%d\t.\t%c\t.\t",refname.chars(),
							guides[i]->getTrackName(), guides[i]->start,guides[i]->end, guides[i]->strand);
					if (guides[i]->getGeneID())
						fprintf(f_out, "gene_id \"%s\"; ", guides[i]->getGeneID()); // this shouldn't come up empty
					fprintf(f_out, "transcript_id \"%s\";",guides[i]->getID());
					if (guides[i]->getGeneName())
						fprintf(f_out, " ref_gene_name \"%s\";", guides[i]->getGeneName());
					fprintf(f_out, " cov \"0.0\";\n");
					for (int e=0;e<guides[i]->exons.Count();++e) {
						fprintf(f_out,"%s\t%s\texon\t%d\t%d\t.\t%c\t.\t",refname.chars(),
								guides[i]->getTrackName(), guides[i]->exons[e]->start, guides[i]->exons[e]->end, guides[i]->strand);
						if (guides[i]->getGeneID())
							fprintf(f_out, "gene_id \"%s\"; ",  guides[i]->getGeneID());
						fprintf(f_out,"transcript_id \"%s\"; exon_number \"%d\";",guides[i]->getID(), e+1);
						if (guides[i]->getGeneName())
							fprintf(f_out, " ref_gene_name \"%s\";", guides[i]->getGeneName());
						fprintf(f_out, " cov \"0.0\";\n");
					}
					//((RC_TData*)guides[i]->uptr)->in_bundle=1;
				}
			}


			GStr gid(guides[i]->getGeneID()); // geneid and strand should determine gene cluster; what if no geneid is present?

			if(gid.is_empty()) {
				gid=guides[i]->getGeneName();
			}


			if(!gid.is_empty()) {
				gid+=guides[i]->strand;
				const int *n=hashgene[gid.chars()];
				if(n) { // I've seen the gene before
					if(guides[i]->start<refgene[*n].start) refgene[*n].start=guides[i]->start;
					if(guides[i]->end>refgene[*n].end) refgene[*n].end=guides[i]->end;
					merge_exons(refgene[*n],guides[i]->exons);
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
	int currentstartneg=-1;
	uint currentendneg=0;
	//int nstartpos=0;
	//int nendpos=0;
	//int nstartneg=0;
	//int nendneg=0;
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
					/*
					switch (sensitivitylevel) {
					case 0: geneno=print_transcript_signcluster('+',pred,genes,transcripts,nstartpos,nendpos,geneno,refname);break;
					case 1: geneno=print_cluster(pospred,genes,transcripts,geneno,refname,refgene,hashgene,predgene,startgno);break;
					case 2: geneno=print_cluster_inclusion(pospred,genes,transcripts,geneno,refname);break;
					case 3: geneno=print_signcluster('+',pred,genes,transcripts,nstartpos,nendpos,geneno,refname);break;
					}
					*/
					geneno=print_cluster(pospred,genes,transcripts,geneno,refname,refgene,hashgene,predgene,startgno,bundleData->sum_cov);
					pospred.Clear();
				}

				pospred.Add(pred[n]);

				currentstartpos=pred[n]->start;
				currentendpos=pred[n]->end;
				//nstartpos=n;
				//nendpos=n;
			}
			else {
				if(pred[n]->end > currentendpos) currentendpos=pred[n]->end;
				//nendpos=n;
				pospred.Add(pred[n]);

			}
		}

		if(pred[n]->strand=='-' || pred[n]->strand=='.') {
			if(pred[n]->start > currentendneg) { // begin new cluster

				// first print predictions I've seen so far
				if(currentstartneg>-1) { // I've seen a cluster before
					/*
					switch (sensitivitylevel) {
					case 0: geneno=print_transcript_signcluster('-',pred,genes,transcripts,nstartneg,nendneg,geneno,refname);break;
					case 1: geneno=print_cluster(negpred,genes,transcripts,geneno,refname,refgene,hashgene,predgene,startgno);break;
					case 2: geneno=print_cluster_inclusion(negpred,genes,transcripts,geneno,refname);break;
					case 3: geneno=print_signcluster('-',pred,genes,transcripts,nstartneg,nendneg,geneno,refname);break;
					}
					*/
					geneno=print_cluster(negpred,genes,transcripts,geneno,refname,refgene,hashgene,predgene,startgno,bundleData->sum_cov);
					negpred.Clear();

				}

				negpred.Add(pred[n]);

				currentstartneg=pred[n]->start;
				currentendneg=pred[n]->end;
				//nstartneg=n;
				//nendneg=n;
			}
			else {
				negpred.Add(pred[n]);
				if(pred[n]->end > currentendneg) currentendneg=pred[n]->end;
				//nendneg=n;
			}
		}
	}

	if(currentstartpos>-1) { // I've seen a cluster before
		/*
		switch (sensitivitylevel) {
		case 0: geneno=print_transcript_signcluster('+',pred,genes,transcripts,nstartpos,nendpos,geneno,refname);break;
		case 1: geneno=print_cluster(pospred,genes,transcripts,geneno,refname,refgene,hashgene,predgene,startgno);break;
		case 2: geneno=print_cluster_inclusion(pospred,genes,transcripts,geneno,refname);break;
		case 3: geneno=print_signcluster('+',pred,genes,transcripts,nstartpos,nendpos,geneno,refname);break;
		}
		*/
		geneno=print_cluster(pospred,genes,transcripts,geneno,refname,refgene,hashgene,predgene,startgno,bundleData->sum_cov);
		pospred.Clear();

	}

	if(currentstartneg>-1) { // I've seen a cluster before
		/*
		switch (sensitivitylevel) {
		case 0: geneno=print_transcript_signcluster('-',pred,genes,transcripts,nstartneg,nendneg,geneno,refname);break;
		case 1: geneno=print_cluster(negpred,genes,transcripts,geneno,refname,refgene,hashgene,predgene,startgno);break;
		case 2: geneno=print_cluster_inclusion(negpred,genes,transcripts,geneno,refname);break;
		case 3: geneno=print_signcluster('-',pred,genes,transcripts,nstartneg,nendneg,geneno,refname);break;
		}
		*/
		geneno=print_cluster(negpred,genes,transcripts,geneno,refname,refgene,hashgene,predgene,startgno,bundleData->sum_cov);
		negpred.Clear();

	}

	hashgene.Clear();
	// I am done printing all transcripts, now evaluate/print the gene abundances
	GVec<float>* bpcov = bundleData->bpcov;
	int refstart=bundleData->start;
	if(eonly  || geneabundance) { // I only need to evaluate the refgene coverages if geneabundance is required, or these are the only gene coverages
		for(int i=0;i<refgene.Count();i++) {
			float cov=0;
			int s=1; // strand of gene
			if(refgene[i].strand=='+') s=2;
			else if(refgene[i].strand=='-') s=0;
			int glen=0;
			for(int j=0;j<refgene[i].exons.Count();j++) { // evaluate unused coverage
				glen+=refgene[i].exons[j].len();
				if(!eonly) {
					int start=(int)refgene[i].exons[j].start-refstart;
					int end=(int)refgene[i].exons[j].end-refstart+1;
					if(start<0) start=0;
					if(end>=bpcov[1].Count()) end=bpcov[1].Count()-1;
					//fprintf(stderr,"geneid=%s ex[%d]: start=%d end=%d gene_start=%d gene_end=%d guides.count=%d refstart=%d\n",
						//refgene[i].geneID,j,refgene[i].exons[j].start,refgene[i].exons[j].end,refgene[i].start,refgene[i].end,guides.Count(),refstart);
					for(int k=start;k<end;k++) {
						switch(s) {
						case 0:
							//cov+=bpcov[1][k]-bpcov[2][k];
							if(bpcov[2][k]) cov+=bpcov[0][k]+(bpcov[1][k]-bpcov[0][k]-bpcov[2][k])*bpcov[0][k]/(bpcov[0][k]+bpcov[2][k]);
							else cov+=bpcov[1][k];
							break;
						case 1: cov+=bpcov[1][k]-bpcov[2][k]-bpcov[0][k];break;
						case 2:
							//cov+=bpcov[1][k]-bpcov[0][k];
							if(bpcov[0][k]) cov+=bpcov[2][k]+(bpcov[1][k]-bpcov[0][k]-bpcov[2][k])*bpcov[2][k]/(bpcov[0][k]+bpcov[2][k]);
							else cov+=bpcov[1][k];
							break;
						}
					}
				}
			}

			if(eonly) bundleData->sum_cov+=refgene[i].covsum;
			else { // for eonly I only use the annotated transcripts, nothing more
				refgene[i].cov=cov-refgene[i].cov;
				if(refgene[i].cov>epsilon) refgene[i].covsum+=refgene[i].cov/glen;
			}


			if(geneabundance) {
				if(!eonly) refgene[i].cov=cov/glen; // only if I want to store the real gene coverage
				else refgene[i].cov/=glen;
				fprintf(f_out,"0 1 %d 0 %.6f\n",glen, refgene[i].covsum);
				const char* geneID=refgene[i].geneID;
				if (geneID==NULL) geneID=".";
				fprintf(f_out,"%s\t",geneID);
				if(refgene[i].geneName) fprintf(f_out,"%s\t",refgene[i].geneName);
				else fprintf(f_out,"-\t");
				//fprintf(f_out,"%s\t%c\t%d\t%d\t%d\t%.6f\n",refname.chars(),refgene[i].strand,refgene[i].start,refgene[i].end,glen,refgene[i].cov);
				fprintf(f_out,"%s\t%c\t%d\t%d\t%.6f\n",refname.chars(),refgene[i].strand,refgene[i].start,refgene[i].end,refgene[i].cov);
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
					case 0:
						//cov+=bpcov[1][k]-bpcov[2][k];
						if(bpcov[2][k]) cov+=bpcov[0][k]+(bpcov[1][k]-bpcov[0][k]-bpcov[2][k])*bpcov[0][k]/(bpcov[0][k]+bpcov[2][k]);
						else cov+=bpcov[1][k];
						break;
					case 1: cov+=bpcov[1][k]-bpcov[2][k]-bpcov[0][k];break;
					case 2:
						//cov+=bpcov[1][k]-bpcov[0][k];
						if(bpcov[0][k]) cov+=bpcov[2][k]+(bpcov[1][k]-bpcov[0][k]-bpcov[2][k])*bpcov[2][k]/(bpcov[0][k]+bpcov[2][k]);
						else cov+=bpcov[1][k];
						break;
					}
				}
			}
			predgene[i].cov=cov-predgene[i].cov; // THIS is the read coverage that is left for the genes after all the predictions were taken into account
			if(predgene[i].cov>epsilon) predgene[i].covsum+=predgene[i].cov/glen;
			bundleData->sum_cov+=predgene[i].covsum;
			if(geneabundance) {
				predgene[i].cov=cov/glen; // only if I want to store the real gene coverage
				fprintf(f_out,"0 1 %d 0 %.6f\n",glen, predgene[i].covsum);
				fprintf(f_out,"%s.%d\t",label.chars(),startgno+i);
				fprintf(f_out,"-\t");
				//fprintf(f_out,"%s\t%c\t%d\t%d\t%d\t%.6f\n",refname.chars(),predgene[i].strand,predgene[i].start,predgene[i].end,glen,predgene[i].cov);
				fprintf(f_out,"%s\t%c\t%d\t%d\t%.6f\n",refname.chars(),predgene[i].strand,predgene[i].start,predgene[i].end,predgene[i].cov);
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
