#include "rlink.h"
#include "GBitVec.h"
#include <float.h>

//import globals from main program:

//extern GffNames* gseqNames;
extern FILE *c_out;         // file handle for the input transcripts that are fully covered by reads

extern bool specific;
extern bool trim;
extern bool partialcov;
extern bool complete;
//extern bool debugMode;
//extern bool verbose;
extern bool eonly;

extern int maxReadCov;

extern float isofrac;
extern float mcov;
extern int mintranscriptlen; // minimum number for a transcript to be printed
extern int sensitivitylevel;
extern int junctionsupport; // anchor length for junction to be considered well supported <- consider shorter??
extern int junctionthr; // number of reads needed to support a particular junction
extern float readthr;     // read coverage per bundle bp to accept it; otherwise considered noise
extern uint bundledist;  // reads at what distance should be considered part of separate bundles
                        // <- this is not addressed everywhere, e.g. in infer_transcripts -> look into this

extern bool includesource;
extern bool EM;
extern bool weight;

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
	nj->nreads+=float(1)/nh;
	if (leftsupport >= junctionsupport && rightsupport >=junctionsupport) {
		nj->nreads_good+=float(1)/nh;
	}
	return nj;
}

void cov_add(GVec<float>& bpcov, int i, int j, float v) {
	if (j>=bpcov.Count()) bpcov.Resize(j+1, 0);
	for (int k=i;k<j;k++)
		bpcov[k]+=v;
}


float getBCov(GVec<float>& bpcov, int p) {
	//if (p<0) GMessage("Error: invalid bpcov index (%d)!\n", p);
	if (p>=bpcov.Count()) return 0;
	else return bpcov[p];
}

bool maxCovReached(int currentstart, GBamRecord& brec, BundleData& bdata) {
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


int processRead(int currentstart, int currentend, BundleData& bdata,
		 GHash<int>& hashread, GBamRecord& brec, char strand, int nh, int hi) {
	GList<CReadAln>& readlist = bdata.readlist;
	GList<CJunction>& junction = bdata.junction;
	GVec<float>& bpcov = bdata.bpcov;
	static uint32_t BAM_R2SINGLE = BAM_FREAD2 | BAM_FMUNMAP ;

	int readstart=brec.start;
	CReadAln* readaln=NULL;
	bool covSaturated=false;
	if (bdata.end<currentend) {
		bdata.start=currentstart;
		bdata.end=currentend;
	}
	//bool first_mapping = (nh==1 || hi==0);
	if (hi==0) {
		for (int i=0;i<brec.exons.Count();i++) {
				bdata.frag_len+=brec.exons[i].len();
		}
		if (!brec.isPaired() || ((brec.flags()&BAM_FREAD1)!=0) ||

				((brec.flags()&BAM_R2SINGLE)==BAM_R2SINGLE ) ) {
			bdata.num_fragments++;
		}
	}
	/*
	if (hi==0 && (!brec.isPaired() || ((brec.flags()&BAM_FREAD1)!=0) ||
				((brec.flags()&BAM_R2SINGLE)==BAM_R2SINGLE ) ) ) {
			bdata.num_fragments++;
	}*/
	if (maxReadCov>0) {
		covSaturated=maxCovReached(currentstart, brec, bdata);
	//if (bdata.firstPass==2 && !covSaturated)
	//single pass only -> unconditional read creation
	  if (!covSaturated)
		  readaln=new CReadAln(strand, nh, brec.start, brec.end);
	}
	//if (bdata.firstPass) { //first pass or singlePass use
		bdata.numreads++;
		//process cigar string
		int leftsupport=0;
		int rightsupport=brec.exons.Last().len();
		int support=0;
		int maxleftsupport=0;
		int njunc=0;
		GVec<int> leftsup;
		GVec<int> rightsup;
		GPVec<CJunction> newjunction(false);
		for (int i=0;i<brec.exons.Count();i++) {
			CJunction* nj=NULL;
			if (i) {
				//deal with introns
				GSeg seg(brec.exons[i-1].end, brec.exons[i].start);
				leftsupport=brec.exons[i-1].len();
				if (leftsupport>maxleftsupport) {
					maxleftsupport=leftsupport;
				}
				leftsup.Add(maxleftsupport);//on the left, always add current max support (?)
				support=brec.exons[i].len();
				rightsup.Add(support);     //..but on the right, real support value is added

				/* previously:
				if (i==brec.exons.Count()-1)
					nj=add_junction(seg.start, seg.end, leftsupport, support, junction, strand, nh);
				else
					nj=add_junction(seg.start, seg.end, maxleftsupport, support, junction, strand, nh);
				*/
				int maxrightsupport = support > rightsupport ? support : rightsupport;
				nj=add_junction(seg.start, seg.end, maxleftsupport, maxrightsupport, junction, strand, nh);
				newjunction.Add(nj);
			}
			cov_add(bpcov, brec.exons[i].start-currentstart,
					brec.exons[i].end-currentstart, float(1)/nh);
			if (readaln) {
				if (nj) readaln->juncs.Add(nj);
				readaln->segs.Add(brec.exons[i]);
			}
		}
		njunc=newjunction.Count();
		//--
		if (njunc>1) { //2 or more introns
			int maxrightsupport=rightsup[njunc-1];
			for (int j=njunc-2;j>=0;j--) {
				if (rightsup[j]>maxrightsupport) { maxrightsupport=rightsup[j];}
				else { //adjust support on the right so that I don't exclude very short internal exons
					if (rightsup[j]<junctionsupport && leftsup[j]>=junctionsupport && maxrightsupport>=junctionsupport){
						newjunction[j]->nreads_good+=float(1)/nh;
					}
				}
			}
		}
	  if((int)brec.end>currentend) {
	  	//return currentend;
	  	currentend=brec.end;
	  	bdata.end=currentend;
	  }
		//if (bdata.firstPass == 1 || covSaturated) {
	  if (covSaturated) {
		  return(currentend); //for 1st-pass only, exit here
		}
	//} //<-- for single-pass or first pass
	//--
	//--> 2nd pass or single-pass run (firstPass==0 or firstPass==2)
		GStr readname(brec.name());
		if (readaln==NULL) {
			//2nd pass
			//TODO: we should try to see if we can collapse this read first
			readaln=new CReadAln(strand, nh, brec.start, brec.end);
		}
		int n=readlist.Add(readaln);
		if (brec.refId()==brec.mate_refId()) {
			//only consider mate pairing data if mates are on the same chromosome/contig
			//TODO: this should be reconsidered if we decide to care about FUSION transcripts!
			GStr readname(brec.name());
			int pairstart=brec.mate_start();
			GStr id(readname); // init id with readname
			id+=':';id+=readstart;id+=':';id+=hi;
			if (readstart<=pairstart) {
				//only add the first mate in a pair
			  if (!hashread[id.chars()])
			  	hashread.Add(id.chars(), new int(n));
			}
			if (readstart>pairstart) {
				//must have seen its pair earlier
			GStr pairid(readname);
			pairid+=':';pairid+=pairstart;pairid+=':';pairid+=hi;
			const int* np=hashread[pairid.chars()];
			if (np) {
				readlist[n]->pair_idx=*np;
				readlist[*np]->pair_idx=n;
						//we can now discard this pair info
						hashread.Remove(pairid.chars());
						hashread.Remove(id.chars()); //just in case it exists
					}
			}
		} //<-- if mate is mapped on the same chromosome
		//int end=brec.end;
	/*
	if (bdata.firstPass==0) {//2nd pass only
		for (int i=0;i<brec.exons.Count();i++) {
			readaln->segs.Add(brec.exons[i]);
			if (i) {
				int jidx=-1;
				CJunction jn(brec.exons[i-1].end, brec.exons[i].start, strand);
				if (junction.Found(&jn, jidx)) {
					readaln->juncs.Add(junction.Get(jidx));
				} else GError("Error: junction %d-%d not found!\n", brec.exons[i-1].end, brec.exons[i].start);
			}
		}
	}
	*/
	return currentend; //currentend already set earlier
	//}
}
/*
int process_read(int currentstart, int currentend, GList<CReadAln>& readlist, GHash<int>& hashread,
		GList<CJunction>& junction, GBamRecord& brec, char strand, int nh, int hi, GVec<float>& bpcov) {
	//uint32_t flag=brec.flags();
	GStr readname(brec.name());
	int readstart=brec.start;

	CReadAln* readaln=new CReadAln(strand, nh, brec.start, brec.end);
	int n=readlist.Add(readaln);
	if (brec.refId()==brec.mate_refId()) {
		//only consider mate pairing data if mates are on the same chromosome/contig
		//TODO: this should be reconsidered if we decide to care about FUSION transcripts!
		GStr readname(brec.name());
	int pairstart=brec.mate_start();
	GStr id(readname); // init id with readname
	id+=':';id+=readstart;id+=':';id+=hi;
		if (readstart<=pairstart) {
			//only add the first mate in a pair
		  if (!hashread[id.chars()])
		hashread.Add(id.chars(), new int(n));
		}
		if (readstart>pairstart) {
			//must have seen its pair earlier
	GStr pairid(readname);
	pairid+=':';pairid+=pairstart;pairid+=':';pairid+=hi;
	const int* np=hashread[pairid.chars()];
	if (np) {
		readlist[n]->pair_idx=*np;
		readlist[*np]->pair_idx=n;
				//we can now discard this pair info
				hashread.Remove(pairid.chars());
				hashread.Remove(id.chars()); //just in case it exists
			}
		}
	}

	bam1_t* b=brec.get_b();
	//process cigar string
	int start=readstart;
	int end=start-1;
	int leftsupport=0;
	int support=0;
	int rightstart=0;
	int leftend=0;
	int maxleftsupport=0;
	GPVec<CJunction> newjunction(false);
	GVec<int> leftsup;
	GVec<int> rightsup;
	int njunc=0;

	for (int i = 0; i < b->core.n_cigar; ++i) {
		int oplen = (bam1_cigar(b)[i]>>BAM_CIGAR_SHIFT);
		char op = "MIDNSHP=X"[bam1_cigar(b)[i]&BAM_CIGAR_MASK];
		switch (op) {
			case 'M':
			case '=':
				support+=oplen; //deliberate case fall through
			case 'X':
			case 'D': //deliberate case fall through
				end+=oplen;
				break;
			case 'N': //intron starts here
				if (start<=end) {
					GSeg seg(start, end);
					readaln->segs.Add(seg);
					cov_add(bpcov, start-currentstart, end-currentstart, float(1)/nh);
				}
				start=end+oplen+1;
				//deal with junction specific staff
				if(leftsupport) {
					//new complete junction
					if (leftsupport>maxleftsupport) {
						//read is split into more than one exon
						maxleftsupport=leftsupport;
					}
					CJunction* nj=add_junction(leftend, rightstart, maxleftsupport, support, junction, strand, nh);
					newjunction.Add(nj);
					rightsup.Add(support);
					leftsup.Add(maxleftsupport);
					njunc++;
					readaln->juncs.Add(nj);
				}
				rightstart=start;
				leftend=end;
				leftsupport=support;
				support=0;
				end=start-1;
				break;
		}//switch(cigar_op)
	} //for each cigar op

	//take care of the last intron:
	if (start<=end) {
		GSeg seg(start, end);
		readaln->segs.Add(seg);
		cov_add(bpcov, start-currentstart, end-currentstart, float(1)/nh);
	}
	if(leftsupport) { // update junction
		if (leftsupport>maxleftsupport) {
			maxleftsupport=leftsupport;
		}
		CJunction* nj=add_junction(leftend, rightstart, leftsupport, support, junction, strand, nh);
		newjunction.Add(nj);
		rightsup.Add(support);
		leftsup.Add(maxleftsupport);
		njunc++;
		readaln->juncs.Add(nj);
	}
	if (njunc>1) { //2 or more introns
		int maxrightsupport=rightsup[njunc-1];
		for (int j=njunc-2;j>=0;j--) {
			if (rightsup[j]>maxrightsupport) { maxrightsupport=rightsup[j];}
			else { //adjust support on the right so that I don't exclude very short internal exons
				if (rightsup[j]<junctionsupport && leftsup[j]>=junctionsupport && maxrightsupport>=junctionsupport){
					newjunction[j]->nreads_good+=float(1)/nh;
				}
			}
		}
	}
	if(end<currentend) { return currentend; }

	return(end);

}
*/

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

int add_read_to_group(int n,GList<CReadAln>& readlist,int color,GPVec<CGroup>& group,
		CGroup **allcurrgroup,CGroup **startgroup,GVec<int> *readgroup,GVec<int>& eqcol,GVec<int>& merge) {

	int sno=readlist[n]->strand+1; // 0: negative strand; 1: zero strand; 2: positive strand
	int readcol=color;

	// check if I've seen read's pair and if yes get its readcol; at the least get read's pair strand if available
	int np=readlist[n]->pair_idx; // pair read number

	if(np>-1 && readlist[np]->nh) { // read pair exists and it wasn't deleted
		if(np<n) { // there is a pair and it came before the current read in sorted order of position
		    // first group of pair read is: $$readgroup[$np][0]
		    int grouppair=readgroup[np][0];
		    while( merge[grouppair]!=grouppair) {
		    	grouppair=merge[grouppair];
		    }
		    readgroup[np][0]=grouppair;

		    readcol=group[readgroup[np][0]]->color;
		    while(eqcol[readcol]!=readcol) { // get smallest color
		    	readcol=eqcol[readcol];
		    }
		    //print STDERR "Adjust color of group ",$$readgroup[$np][0]," to $readcol\n";
		    group[readgroup[np][0]]->color=readcol;
		}
		else { // it's the first time I see the read in the fragment

		    // see if I have the correct read strand
		    char snop=readlist[np]->strand+1;
		    if(sno!=snop) { // different strands
		    	if(sno==1 && snop!=1) { // read n is on zero strand
		    		readlist[n]->strand=readlist[np]->strand;
		    		sno=snop;
		    	}
		    	else if(snop==1 && sno!=1) { // pair read np is on zero strand
		    		readlist[np]->strand=readlist[n]->strand;
		    	}
		    	else { // conflicting strands -> un-pair reads in the hope that one is right
		    		readlist[n]->pair_idx=-1;
		    		readlist[np]->pair_idx=-1;
		    	}
		    }
		    eqcol.Add(color);
		    color++;
		}
	} // if(np>-1 && readlist[np]->nh) : read pair exists and it wasn't deleted
	else {
		eqcol.Add(color);
		color++;
	}

	CGroup *currgroup=allcurrgroup[sno];

	if(currgroup != NULL) { // this type of group - negative, unknown, or positive - was created before

		//set currgroup first
		CGroup *lastgroup=NULL;
		while(currgroup!=NULL && readlist[n]->start > currgroup->end) {
		    lastgroup=currgroup;
		    currgroup=currgroup->next_gr;
		}

		if(currgroup==NULL || readlist[n]->segs[0].end < currgroup->start) // currgroup is null only if we reached end of currgroup list
			currgroup=lastgroup;

		// now process each group of coordinates individually
		CGroup *thisgroup=currgroup;
		int ncoord=readlist[n]->segs.Count();
		int lastpushedgroup=-1;
		bool added=false;

		for(int i=0;i<ncoord;i++) {

		    // skip groups that are left behind
		    while(thisgroup!=NULL && readlist[n]->segs[i].start > thisgroup->end) {
		    	lastgroup=thisgroup;
		    	thisgroup=thisgroup->next_gr;
		    }

		    if(thisgroup && readlist[n]->segs[i].end >= thisgroup->start) { // read overlaps group

		    	if(!added) {
		    		thisgroup->nread+=(float)1/readlist[n]->nh;
		    		if(readlist[n]->nh>1) thisgroup->multi+=(float)1/readlist[n]->nh;
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
		    		readgroup[n].Add(thisgroup->grid);
		    		lastpushedgroup=thisgroup->grid;
		    	}
		    	thisgroup->cov_sum+=(float)(readlist[n]->segs[i].end-readlist[n]->segs[i].start+1)/readlist[n]->nh;
		    } // end if(thisgroup && readlist[n]->segs[i].end >= thisgroup->start)
		    else { // read is at the end of groups, or read is not overlapping other groups -> $lastgroup should be not null here

		    	int ngroup=group.Count();
		    	float nread=0;
		    	float multi=0;
		    	if(!added) {
		    		nread=(float)1/readlist[n]->nh;
		    		if(readlist[n]->nh>1) multi=(float)1/readlist[n]->nh;
		    		added=true;
		    	}
		    	CGroup *newgroup=new CGroup(readlist[n]->segs[i].start,readlist[n]->segs[i].end,readcol,ngroup,(float)(readlist[n]->segs[i].end-readlist[n]->segs[i].start+1)/readlist[n]->nh,nread,multi);
		    	group.Add(newgroup);
		    	merge.Add(ngroup);
		    	lastgroup->next_gr=newgroup;
		    	newgroup->next_gr=thisgroup;
		    	readgroup[n].Add(ngroup);
		    	lastpushedgroup=ngroup;
				thisgroup=lastgroup;
		    }
		} // for(int i=0;i<ncoord;i++)
	} // if(currgroup != NULL)
	else { // create new group of this type

		int ncoord=readlist[n]->segs.Count();
		CGroup *lastgroup=NULL;
		float nread=(float)1/readlist[n]->nh;
		float multi=0;
		if(readlist[n]->nh>1) multi=(float)1/readlist[n]->nh;
		for(int i=0;i<ncoord;i++) {
			int ngroup=group.Count();
			CGroup *newgroup=new CGroup(readlist[n]->segs[i].start,readlist[n]->segs[i].end,readcol,ngroup,(float)(readlist[n]->segs[i].end-readlist[n]->segs[i].start+1)/readlist[n]->nh,nread,multi);
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

void find_trims(int refstart,uint start,uint end,GVec<float>& bpcov,uint& sourcestart,float& sourceabundance,uint& sinkend,
		float& sinkabundance){

	if(end-start<2*CHI_WIN-1) return;

	float sumleft=0;
	float sumright=0;

	//float maxsinkchi=0;
	//float maxsourcechi=0;

	float maxsinkabundance=0;
	float maxsourceabundance=0;

	GArray<float> winleft(CHI_WIN,false); // not auto-sort
	GArray<float> winright(CHI_WIN,false); // not auto-sort

	for(uint i=start;i<=end;i++) {

		if(i-start<2*CHI_WIN-1)  { // I have to compute the sumleft and sumright first
			if(i-start<CHI_WIN) {
				sumleft+=bpcov[i-refstart];
				winleft.Add(bpcov[i-refstart]);
			}
			else {
				sumright+=bpcov[i-refstart];
				winright.Add(bpcov[i-refstart]);
				if(i-start==2*CHI_WIN-2) {
					winleft.setSorted(true);
					winright.setSorted(true);
				}
			}
	    }
	    else { // I can do the actual sumleft, sumright comparision

	    	sumright+=bpcov[i-refstart];
			winright.Add(bpcov[i-refstart]);

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

	    	sumleft-=bpcov[i-refstart-2*CHI_WIN+1];
	    	int idx=winleft.IndexOf(bpcov[i-refstart-2*CHI_WIN+1]);
	    	winleft.Delete(idx);
	    	sumleft+=bpcov[i-refstart-CHI_WIN+1];
	    	winleft.Add(bpcov[i-refstart-CHI_WIN+1]);
	    	sumright-=bpcov[i-refstart-CHI_WIN+1];
	    	idx=winright.IndexOf(bpcov[i-refstart-CHI_WIN+1]);
	    	winright.Delete(idx);
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

CGraphnode *trimnode(int s, int g, int refstart,uint newend, CGraphnode *graphnode,CGraphnode *source, CGraphnode *sink, GVec<float>& bpcov,
		GVec<float>& futuretr, int& graphno,CBundlenode *bundlenode,GVec<CGraphinfo> **bundle2graph,GPVec<CGraphnode> **no2gnode) {

	uint sourcestart=0;
	uint sinkend=0;
	float sinkabundance=0;
	float sourceabundance=0;
	find_trims(refstart,graphnode->start,newend,bpcov,sourcestart,sourceabundance,sinkend,sinkabundance);

	if(sourcestart < sinkend) { // source trimming comes first

		if(sourcestart) { // there is evidence of graphnode trimming from source
			graphnode->end=sourcestart-1;
			CGraphnode *prevnode=graphnode;
			graphnode=create_graphnode(s,g,sourcestart,newend,graphno,bundlenode,bundle2graph,no2gnode);
			graphno++;
			source->child.Add(graphnode->nodeid);  // this node is the child of source
			graphnode->parent.Add(source->nodeid); // this node has source as parent
			prevnode->child.Add(graphnode->nodeid); // this node is the child of previous node
			graphnode->parent.Add(prevnode->nodeid); // this node has as parent the previous node
			float tmp=graphno-1;
			futuretr.Add(tmp);
			futuretr.cAdd(0.0);
			futuretr.Add(sourceabundance);
		}

		// sinkend is always positive since it's bigger than sourcestart
		float tmp=graphno-1;
		graphnode->end=sinkend;
		CGraphnode *prevnode=graphnode;
		graphnode=create_graphnode(s,g,sinkend+1,newend,graphno,bundlenode,bundle2graph,no2gnode);
		graphno++;
		prevnode->child.Add(graphnode->nodeid); // this node is the child of previous node
		graphnode->parent.Add(prevnode->nodeid); // this node has as parent the previous node
		sink->parent.Add(prevnode->nodeid); // prevnode is the parent of sink
		// remember to create transfrag as well -> I don't know the gno yet, so I can not create it here
		futuretr.Add(tmp);
		futuretr.cAdd(1.0);
		futuretr.Add(sinkabundance);
	}
	else if(sourcestart > sinkend) { // sink trimming comes first
		if(sinkend) {
			graphnode->end=sinkend;
			CGraphnode *prevnode=graphnode;
			graphnode=create_graphnode(s,g,sinkend+1,newend,graphno,bundlenode,bundle2graph,no2gnode);
			graphno++;
			prevnode->child.Add(graphnode->nodeid); // this node is the child of previous node
			graphnode->parent.Add(prevnode->nodeid); // this node has as parent the previous node
			sink->parent.Add(prevnode->nodeid); // prevnode is the parent of sink
			// remember to create transfrag as well -> I don't know the gno yet, so I can not create it here
			float tmp=graphno-2;
			futuretr.Add(tmp);
			futuretr.cAdd(1.0);
			futuretr.Add(sinkabundance);
		}

		// sourcestart is positive since it's bigger than sinkend
		graphnode->end=sourcestart-1;
		CGraphnode *prevnode=graphnode;
		graphnode=create_graphnode(s,g,sourcestart,newend,graphno,bundlenode,bundle2graph,no2gnode);
		graphno++;
		source->child.Add(graphnode->nodeid);  // this node is the child of source
		graphnode->parent.Add(source->nodeid); // this node has source as parent
		prevnode->child.Add(graphnode->nodeid); // this node is the child of previous node
		graphnode->parent.Add(prevnode->nodeid); // this node has as parent the previous node
		float tmp=graphno-1;
		futuretr.Add(tmp);
		futuretr.cAdd(0.0);
		futuretr.Add(sourceabundance);
	}
	// else both source and sink trimming are not present

	return(graphnode);
}

inline int edge(int min, int max, int gno) {
	//return((gno-1)*min-min*(min-1)/2+max-min); // this should be changed if source to node edges are also stored
	return((gno-1)*(min+1)-min*(min-1)/2+max-min); // this includes source to node edges
}

GBitVec traverse_dfs(int s,int g,CGraphnode *node,CGraphnode *sink,GBitVec parents,int gno,GVec<bool>& visit,
		GPVec<CGraphnode> **no2gnode,GPVec<CTransfrag> **transfrag){

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
		node->childpat.resize(1+gno*(gno+1)/2);
		node->parentpat.resize(1+gno*(gno+1)/2);
		node->parentpat = node->parentpat | parents;
		visit[node->nodeid]=true;
		parents[node->nodeid]=1; // add the node to the parents

		if(node->parent.Count()==1 && !node->parent[0]) { // node has source only as parent -> add transfrag from source to node
			GBitVec trpat(1+gno*(gno+1)/2);
			trpat[0]=1;
			trpat[node->nodeid]=1;
			trpat[edge(0,node->nodeid,gno)]=1;
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
			GBitVec trpat(1+gno*(gno+1)/2);
			trpat[node->nodeid]=1;
			trpat[gno-1]=1;
			trpat[edge(node->nodeid,gno-1,gno)]=1;
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

	    for(int i=0; i< n; i++) { // for all children
	    	GBitVec childparents=parents;
	    	int min=node->nodeid; // nodeid is always smaller than child node ?
	    	int max=node->child[i];
	    	if(min>max) {
	    		max=node->nodeid; // nodeid is always smaller than child node ?
	    		min=node->child[i];
	    	}
	    	int edg=edge(min,max,gno);
	    	childparents[edg]=1; // add edge from node to child to the set of parents from child
	    	node->childpat[edg]=1; // add edge from node to child to the set of node children
	    	node->childpat = node->childpat | traverse_dfs(s,g,no2gnode[s][g][node->child[i]],sink,childparents,gno,visit,no2gnode,transfrag);
	    }
	} // end else from if(visit[node->nodeid])

	GBitVec children = node->childpat;
	children[node->nodeid]=1;

	return(children);
}

int create_graph(int refstart,int s,int g,CBundle *bundle,GPVec<CBundlenode>& bnode, GList<CJunction>& junction,GList<CJunction>& ejunction,GVec<CGraphinfo> **bundle2graph,
		GPVec<CGraphnode> **no2gnode,GPVec<CTransfrag> **transfrag,GVec<float>& bpcov){

	CGraphnode* source=new CGraphnode(0,0,0);
	no2gnode[s][g].Add(source);
	CGraphnode* sink=new CGraphnode();

	int njunctions=junction.Count();

	/*
	{ // DEBUG ONLY
		fprintf(stderr,"Junctions[%d][%d]: ",s,g);
		for(int i=0;i<njunctions;i++) fprintf(stderr," %d-%d",junction[i]->start,junction[i]->end);
	}
	*/

	int njs=0; // index of sorted junction starts
	int nje=0; // index of sorted junction ends

	int graphno=1; // number of nodes in graph
	GHash<GVec<int> > ends;

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

	    	int minjunction = -1; // process next junction -> either a start or an ending whichever has the first position on the genome; if they have same position then process ending first
	    	if((nje<njunctions && (ejunction[nje]->end<endbundle)) || (njs<njunctions && (junction[njs]->start<=endbundle))) {
	    		if(nje<njunctions) { // there are still junctions endings
	    			if(njs<njunctions) { // there are still junctions starting
	    				minjunction = junction[njs]->start >= ejunction[nje]->end ? 1 : 0;
	    			}
	    			else minjunction = 1;
	    		}
	    		else minjunction = 0;
	    	}

	    	if(minjunction == 0 ) { // found a start junction here

	    		if(trim) graphnode=trimnode(s,g,refstart,junction[njs]->start,graphnode,source,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode);// do something to find intermediate nodes; alternatively, I could only do this for end nodes

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

	    			if(trim) graphnode=trimnode(s,g,refstart,pos-1,graphnode,source,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode);// do something to find intermediate nodes; alternatively, I could only do this for end nodes

	    			graphnode->end=pos-1; // set end of current graphnode here
	    			CGraphnode *nextnode = create_graphnode(s,g,pos,endbundle,graphno,bundlenode,bundle2graph,no2gnode);
	    			graphno++;
	    			graphnode->child.Add(nextnode->nodeid); // make nextnode a child of current graphnode
	    			nextnode->parent.Add(graphnode->nodeid);// make graphnode a parent of nextnode
	    			graphnode=nextnode;
	    		}

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

	      if(trim) // do something to find intermediate nodes; alternatively, I could only do this for end nodes
	    	  graphnode=trimnode(s,g,refstart,endbundle,graphnode,source,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode);

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
	GBitVec parents(1+(graphno+1)*graphno/2);
	traverse_dfs(s,g,source,sink,parents,graphno,visit,no2gnode,transfrag);

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

void get_read_pattern(GBitVec& pattern0,GBitVec& pattern1,int *rgno, GVec<int> *rnode,GList<CReadAln>& readlist,int n,
		GVec<int> *readgroup,GVec<int>& merge,GVec<int> *group2bundle,GVec<CGraphinfo> **bundle2graph,GVec<int> *graphno,GPVec<CGraphnode> **no2gnode) {

	int lastgnode[2]={-1,-1}; // lastgnode[0] is for - strand; [1] is for + strand -> I need these in order to add the edges to the read pattern; check this: if it's not correct than storage was wrong!
	int ncoord=readlist[n]->segs.Count();

	int k[2]={0,0}; // need to keep track of coordinates already added to coverages of graphnodes
    bool valid[2]={true,true};

    for(int i=0;i<readgroup[n].Count();i++)
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
    								node->cov+=float(bp)/readlist[n]->nh;
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
    								if(s) pattern1[edge(min,max,graphno[s][ngraph])]=1; // added edge from previous gnode to current one for the read
    								else pattern0[edge(min,max,graphno[s][ngraph])]=1; // added edge from previous gnode to current one for the read
    								rnode[s].Add(gnode);
    							}
    							else { // first time considering read
    								rgno[s]=ngraph;
    								rnode[s].Add(gnode);
    							}
    							lastgnode[s]=gnode;
    							if(s) {
    								pattern1.resize(1+(graphno[s][ngraph]+1)*graphno[s][ngraph]/2);
    								pattern1[gnode]=1; // here I could remember kids as well to speed things up
    							}
    							else {
    								pattern0.resize(1+(graphno[s][ngraph]+1)*graphno[s][ngraph]/2);
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

CTreePat *construct_treepat(int gno, GPVec<CTransfrag>& transfrag) {

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
						if(transfrag[t]->pattern[edge(m,n,gno)]) // there is an edge between m and n
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

CTransfrag *findtrf_in_treepat(int gno,GVec<int>& node,GBitVec& pattern,CTreePat *tr2no) { // doesn't work for patterns including source node

	if(!tr2no) return(NULL);

	CTreePat *tree=tr2no;
	for(int n=0;n<node.Count();n++) {
		if(n) { // not the first node in pattern
			if(pattern[edge(node[n-1],node[n],gno)]) // there is an edge between node[n-1] and node[n]
				tree=tree->nextpat[gno-1-node[n-1]+node[n]-node[n-1]-1];
			else tree=tree->nextpat[node[n]-node[n-1]-1];
		}
		else tree=tree->nextpat[node[n]-1];
		if(!tree) return(NULL);
	}

	return(tree->tr);
}

bool istrf_in_treepat(int gno,GVec<int>& node,GBitVec& pattern,CTreePat *tr2no) {

	if(!tr2no) return(false);

	CTreePat *tree=tr2no;
	for(int n=0;n<node.Count();n++) {
		if(n) { // not the first node in pattern
			if(pattern[edge(node[n-1],node[n],gno)]) // there is an edge between node[n-1] and node[n]
				tree=tree->nextpat[gno-1-node[n-1]+node[n]-node[n-1]-1];
			else tree=tree->nextpat[node[n]-node[n-1]-1];
		}
		else tree=tree->nextpat[node[n]-1];
		if(!tree) return(false);
	}

	return(true);
}

void update_abundance(int s,int g,int gno,GBitVec& pattern,float abundance,GVec<int>& node,GPVec<CTransfrag> **transfrag,
		CTreePat ***tr2no){

	CTransfrag *t=findtrf_in_treepat(gno,node,pattern,tr2no[s][g]);
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
				if(pattern[edge(node[n-1],node[n],gno)]) // there is an edge between node[n-1] and node[n]
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

void get_fragment_pattern(GList<CReadAln>& readlist,int n, int np,GVec<int> *readgroup,GVec<int>& merge,
		GVec<int> *group2bundle,GVec<CGraphinfo> **bundle2graph,GVec<int> *graphno,GPVec<CGraphnode> **no2gnode,
		GPVec<CTransfrag> **transfrag,CTreePat ***tr2no) {

	GBitVec rpat[2];
	int rgno[2]={-1,-1};
	GVec<int> rnode[2];
	if(readlist[n]->nh) get_read_pattern(rpat[0],rpat[1],rgno,rnode,readlist,n,readgroup,merge,group2bundle,bundle2graph,graphno,no2gnode);

	GBitVec ppat[2];
	int pgno[2]={-1,-1};
	GVec<int> pnode[2];
	// get pair pattern if pair exists and it hasn't been deleted
	if(np>-1 && readlist[np]->nh) {
		get_read_pattern(ppat[0],ppat[1],pgno,pnode,readlist,np,readgroup,merge,group2bundle,bundle2graph,graphno,no2gnode);
	}

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
						update_abundance(s,rgno[s],graphno[s][rgno[s]],rpat[s],float(1)/readlist[n]->nh,rnode[s],transfrag,tr2no);
					}
				}
				if(conflict) { // update both patterns separately
					update_abundance(s,rgno[s],graphno[s][rgno[s]],rpat[s],float(1)/readlist[n]->nh,rnode[s],transfrag,tr2no);
					update_abundance(s,pgno[s],graphno[s][pgno[s]],ppat[s],float(1)/readlist[np]->nh,pnode[s],transfrag,tr2no);
				}
			}
			else { // pair has no valid pattern
				update_abundance(s,rgno[s],graphno[s][rgno[s]],rpat[s],float(1)/readlist[n]->nh,rnode[s],transfrag,tr2no);
			}
		}
		else // read has no valid pattern but pair might
			if(pgno[s]>-1) {
				update_abundance(s,pgno[s],graphno[s][pgno[s]],ppat[s],float(1)/readlist[np]->nh,pnode[s],transfrag,tr2no);
			}
	}

}

void settrf_in_treepat(CTransfrag *t,int gno,GVec<int>& node,GBitVec& pattern,CTreePat *tr2no) {

	if(!tr2no) return;

	CTreePat *tree=tr2no;
	for(int n=0;n<node.Count();n++) {
		if(n) { // not the first node in pattern
			if(pattern[edge(node[n-1],node[n],gno)]) // there is an edge between node[n-1] and node[n]
				tree=tree->nextpat[gno-1-node[n-1]+node[n]-node[n-1]-1];
			else tree=tree->nextpat[node[n]-node[n-1]-1];
		}
		else tree=tree->nextpat[node[n]-1];
		if(!tree) return;
	}
	tree->tr=t;
}

void eliminate_transfrags_under_thr(int gno,GPVec<CTransfrag>& transfrag, CTreePat *tr2no,float threshold) {

	for(int t=transfrag.Count()-1;t>=0;t--)
		if(transfrag[t]->abundance<threshold && transfrag[t]->nodes[0]) { // need to delete transfrag that doesn't come from source
			settrf_in_treepat(NULL,gno,transfrag[t]->nodes,transfrag[t]->pattern,tr2no); // this should be eliminated if I want to store transcripts from 0 node
			transfrag.Exchange(t,transfrag.Count()-1);
			transfrag.Delete(transfrag.Count()-1);
		}
}

bool conflict(int &i,int node,GVec<int>& trnode,int n,GPVec<CGraphnode>& no2gnode,GBitVec& trpat,int gno) {

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
	  if(trpat[edge(node1,node2,gno)])
		  return(true);
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

void process_transfrags(int gno,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,CTreePat *tr2no,
		GVec<bool>& compatible) {

	/*
	{ // DEBUG ONLY
		printTime(stderr);
		fprintf(stderr,"There are %d transfrags before clean up\n",transfrag.Count());
	}
	*/

	// eliminate transfrags below threshold (they represent noise) if they don't come from source
	eliminate_transfrags_under_thr(gno,transfrag,tr2no,trthr);

	// introduce "fake" transcripts as holders to adjust the abundances, if the EM algorithm is not used -> I don't need to to this for source and sink

	//if(!EM) {
	if(0) {
		GVec<float> n[gno]; // abundances of all transfrags passing through node and leaving node
		GVec<float> f[gno]; // abundances of all transfrags leaving node
		//GVec<int> u[gno]; // transfrags that are uncommited to a specific child; this is a more complicated issue -> ignore for now
		for(int i=1;i<gno-1;i++) {
			n[i].Resize(gno,0);
			f[i].Resize(gno,0);
		}

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
					GBitVec trpat(1+(gno+1)*gno/2);
					trpat[i]=1;
					trpat[c]=1;
					trpat[edge(i,c,gno)]=1;
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
							GBitVec trpat(1+(gno+1)*gno/2);
							trpat[i]=1;
							trpat[no2gnode[i]->child[k]]=1;
							trpat[edge(i,no2gnode[i]->child[k],gno)]=1;
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
					if(!transfrag[t1]->pattern[edge(transfrag[t1]->nodes[n-1],transfrag[t1]->nodes[n],gno)]) // there is no edge between node[n-1] and node[n]
						assign_incomplete_trf_to_nodes(t1,transfrag[t1]->nodes[n-1],transfrag[t1]->nodes[n],no2gnode); 	// this ensures that I still have compatibilities when going through a certain path: !!! THIS IS NEW COMPARED TO PERL
				}
				else if(n) { // last but not first node
					// add t1 to in of node
					no2gnode[transfrag[t1]->nodes[n]]->trf.Add(t1);

					// check if transfrag t1 is incomplete between node[n-1] and node [n]
					if(!transfrag[t1]->pattern[edge(transfrag[t1]->nodes[n-1],transfrag[t1]->nodes[n],gno)]) // there is no edge between node[n-1] and node[n]
						assign_incomplete_trf_to_nodes(t1,transfrag[t1]->nodes[n-1],transfrag[t1]->nodes[n],no2gnode); 	// this ensures that I still have compatibilities when going through a certain path: !!! THIS IS NEW COMPARED TO PERL
				}
				else { // first node -> only add transfrag to out of node
					no2gnode[transfrag[t1]->nodes[n]]->trf.Add(t1);
				}
			}
		}
		else { // this transcript is included completely in node
			no2gnode[transfrag[t1]->nodes[0]]->frag+=transfrag[t1]->abundance;
		}

		// add t1 to t1 compatibility
		bool comp=true;
		compatible.Add(comp);
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
					if(conflict(i1,transfrag[t2]->nodes[i2],transfrag[t1]->nodes,n1,no2gnode,transfrag[t1]->pattern,gno)) {
						comp=false;
						break;
					}
				}
				else {
					i2++;
					if(conflict(i2,transfrag[t1]->nodes[i1],transfrag[t2]->nodes,n2,no2gnode,transfrag[t2]->pattern,gno)) {
						comp=false;
						break;
					}
				}
			}
			compatible.Add(comp);
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

inline int comptbl_pos(int t1,int t2,int n){
	return(t2+t1*(2*n-t1-1)/2);
}

bool is_compatible(int t1,int t2, int n,GVec<bool>& compatible) {
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

GVec<int> *max_compon_size_with_penalty(int trnumber,float &maxsize,GVec<CTrInfo>& set,GVec<bool>& compatible,
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

bool onpath(GBitVec& trpattern,GVec<int>& trnode,GBitVec& pathpattern,int mini,int maxi,GPVec<CGraphnode>& no2gnode,int gno) {

	if(trnode[0]<mini) // mini can be reached through transcript
	    if(!no2gnode[mini]->parentpat[trnode[0]])	return false;


	if(trnode.Last()>maxi) // from maxi I can reach end of transcript
	    if(!no2gnode[maxi]->childpat[trnode.Last()]) return false;

	int first=1;

	for(int i=0;i<trnode.Count();i++) {
		if(trnode[i]>=mini && trnode[i]<=maxi) {
	      if(!pathpattern[trnode[i]]) return false;
	      if(first) {
	    	  first=0;
	    	  if(i && trnode[i]>mini && trpattern[edge(trnode[i-1],trnode[i],gno)])
	    		  return false;
	    	  if(i && !no2gnode[mini]->parentpat[trnode[i-1]]) return false; // I can not reach mini from previous node in transfrag
	      }
	      else if(i && trpattern[edge(trnode[i-1],trnode[i],gno)] && !pathpattern[edge(trnode[i-1],trnode[i],gno)]) return false;
	    }
	    if(trnode[i]>maxi) {
	    	if(i && trnode[i-1]<maxi && trpattern[edge(trnode[i-1],trnode[i],gno)])
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
		GVec<float>& nodecov,int gno){

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
		if(sensitivitylevel && inode->child[c]==i+1 && i<gno-2 && inode->end+1==no2gnode[i+1]->start && nodecov[i]*DROP>nodecov[i+1])  { // adjacent to parent
			exclude=true;
		}
		else {
			pathpat[inode->child[c]]=1;
			pathpat[edge(i,inode->child[c],gno)]=1;
			for(int j=0;j<cnode->trf.Count();j++) { // for all transfrags going through child
				int t=cnode->trf[j];
				if(transfrag[t]->abundance<epsilon) { // this transfrag was used before -> needs to be deleted
					cnode->trf.Delete(j);
					j--;
				}
				else if(transfrag[t]->real && transfrag[t]->nodes[0]<=i && transfrag[t]->nodes.Last()>=inode->child[c] &&   // transfrag goes from i to c
						(transfrag[t]->pattern[inode->child[c]] || transfrag[t]->pattern[i]) &&    // transfrag is not incomplete through these nodes
						onpath(transfrag[t]->pattern,transfrag[t]->nodes,pathpat,path[0],inode->child[c],no2gnode,gno)) { // transfrag is compatible with path
					//fprintf(stderr,"add transfr[%d]->abund=%f\n",t,transfrag[t]->abundance);
					childcov+=transfrag[t]->abundance;
				}
			}

			if(childcov>maxcov) {
				maxcov=childcov;
				maxc=inode->child[c];
			}

			pathpat[inode->child[c]]=0;
			pathpat[edge(i,inode->child[c],gno)]=0;
		}
	}
	if(maxc==-1) {
		if(exclude && nodecov[i+1]) {
			CGraphnode *cnode=no2gnode[i+1];
			float childcov=0;
			pathpat[i+1]=1;
			pathpat[edge(i,i+1,gno)]=1;
			for(int j=0;j<cnode->trf.Count();j++) { // for all transfrags going through child
				int t=cnode->trf[j];
				if(transfrag[t]->abundance<epsilon) { // this transfrag was used before -> needs to be deleted
					cnode->trf.Delete(j);
					j--;
				}
				else if(transfrag[t]->real && transfrag[t]->nodes[0]<=i && transfrag[t]->nodes.Last()>=i+1 &&   // transfrag goes from i to c
						(transfrag[t]->pattern[i+1] || transfrag[t]->pattern[i]) &&    // transfrag is not incomplete through these nodes
						onpath(transfrag[t]->pattern,transfrag[t]->nodes,pathpat,path[0],i+1,no2gnode,gno)) { // transfrag is compatible with path
					childcov+=transfrag[t]->abundance;
				}
			}
			pathpat[i+1]=1;
			pathpat[edge(i,i+1,gno)]=1;

			if(childcov) {
				maxc=i+1;
			}
			else return false;
		}
		else return false; //maxc=maxchild;
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
	pathpat[edge(i,maxc,gno)]=1;

	return fwd_to_sink_fast(maxc,path,pathpat,transfrag,no2gnode,nodecov,gno);
}

bool back_to_source_fast(int i,GVec<int>& path,GBitVec& pathpat,GPVec<CTransfrag>& transfrag,GPVec<CGraphnode>& no2gnode,
		GVec<float>& nodecov,int gno){

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
		if(sensitivitylevel && inode->parent[p]==i-1 && i>1 && inode->start==no2gnode[i-1]->end+1 && nodecov[i]*DROP>nodecov[i-1])  { // adjacent to parent
			exclude=true;
		}
		else {
			pathpat[inode->parent[p]]=1;
			pathpat[edge(inode->parent[p],i,gno)]=1;
			for(int j=0;j<pnode->trf.Count();j++) { // for all transfrags going through parent
				int t=pnode->trf[j];
				if(transfrag[t]->abundance<epsilon) { // this transfrag was used before -> needs to be deleted
					pnode->trf.Delete(j);
					j--;
				}
				else if(transfrag[t]->real && transfrag[t]->nodes[0]<=inode->parent[p] && transfrag[t]->nodes.Last()>=i &&   // transfrag goes from p to i
						(transfrag[t]->pattern[inode->parent[p]]||transfrag[t]->pattern[i]) &&    // transfrag is not incomplete through these nodes
						onpath(transfrag[t]->pattern,transfrag[t]->nodes,pathpat,inode->parent[p],path[0],no2gnode,gno)) { // transfrag is compatible with path
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
			pathpat[edge(inode->parent[p],i,gno)]=0;
		}
	}
	if(maxp==-1) {
		if(exclude && nodecov[i-1]) {
			CGraphnode *pnode=no2gnode[i-1];
			float parentcov=0;
			pathpat[i-1]=1;
			pathpat[edge(i-1,i,gno)]=1;
			for(int j=0;j<pnode->trf.Count();j++) { // for all transfrags going through parent
				int t=pnode->trf[j];
				if(transfrag[t]->abundance<epsilon) { // this transfrag was used before -> needs to be deleted
					pnode->trf.Delete(j);
					j--;
				}
				else if(transfrag[t]->real && transfrag[t]->nodes[0]<=i-1 && transfrag[t]->nodes.Last()>=i &&   // transfrag goes from p to i
						(transfrag[t]->pattern[i-1]||transfrag[t]->pattern[i]) &&    // transfrag is not incomplete through these nodes
						onpath(transfrag[t]->pattern,transfrag[t]->nodes,pathpat,i-1,path[0],no2gnode,gno)) { // transfrag is compatible with path
					// comment: I need to check for the parent to make sure I am not going through another node!!
					//fprintf(stderr,"parent=%d add transfr[%d]->abund=%f\n",i-1,t,transfrag[t]->abundance);
					parentcov+=transfrag[t]->abundance;
				}
			}
			pathpat[i-1]=0;
			pathpat[edge(i-1,i,gno)]=0;

			if(parentcov) {
				maxp=i-1;
			}
			else return false;
		}
		else return false; //maxp=maxparent;
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
	pathpat[edge(maxp,i,gno)]=1;


	return back_to_source_fast(maxp,path,pathpat,transfrag,no2gnode,nodecov,gno);
}

bool back_to_source_path(int i,GVec<int>& path,GBitVec& pathpat,GVec<float>& pathincov, GVec<float>& pathoutcov,
		GBitVec& istranscript,GBitVec& removable,GPVec<CTransfrag>& transfrag,GHash<CComponent>& computed,
		GVec<bool>& compatible,GPVec<CGraphnode>& no2gnode,GVec<float>& nodecov,int gno){

	// find all parents -> if parent is source then go back
	CGraphnode *inode=no2gnode[i];

	if(!inode->parent[0]) { // parent is source; what if it has other parent than source?
		pathpat[0]=1;
		pathpat[edge(0,i,gno)]=1;
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
			if(istranscript[t] || onpath(transfrag[t]->pattern,transfrag[t]->nodes,pathpat,i,path[0],no2gnode,gno)) {
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
			if(back_to_source_fast(i,path,pathpat,transfrag,no2gnode,nodecov,gno)) {
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
						if(transfrag[t]->real && transfrag[t]->abundance && onpath(transfrag[t]->pattern,transfrag[t]->nodes,pathpat,path.Last(),path[0],no2gnode,gno)) {

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
					if(transfrag[t]->abundance && onpath(transfrag[t]->pattern,transfrag[t]->nodes,pathpat,0,path[0],no2gnode,gno)) {
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
							if(transfrag[t]->pattern[edge(inode->parent[p],i,gno)]) {
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
	pathpat[edge(maxp,i,gno)]=1;

	if(!maxp) return(true);

	// add maxp to path
	path.Add(maxp);

	return back_to_source_path(maxp,path,pathpat,pathincov,pathoutcov,istranscript,removable,transfrag,computed,compatible,no2gnode,nodecov,gno);
}

bool fwd_to_sink_path(int i,GVec<int>& path,GBitVec& pathpat,GVec<float>& pathincov, GVec<float>& pathoutcov,
		GBitVec& istranscript,GBitVec& removable,GPVec<CTransfrag>& transfrag,GHash<CComponent>& computed,
		GVec<bool>& compatible,GPVec<CGraphnode>& no2gnode,GVec<float>& nodecov,int gno){

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
			if(istranscript[t] || onpath(transfrag[t]->pattern,transfrag[t]->nodes,pathpat,path[0],i,no2gnode,gno)) {
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
			if(fwd_to_sink_fast(i,path,pathpat,transfrag,no2gnode,nodecov,gno)) {
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
						if(transfrag[t]->pattern[edge(i,inode->child[c],gno)]) {
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
	pathpat[edge(i,maxc,gno)]=1;

	return fwd_to_sink_path(maxc,path,pathpat,pathincov,pathoutcov,istranscript,removable,transfrag,computed,compatible,no2gnode,nodecov,gno);
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
		GVec<float>& nodecapacity,GBitVec& pathpat,float& fragno) {

	float flux=0;
	int n=path.Count();
	GVec<float> *capacity=new GVec<float>[n]; // capacity of edges in network
	GVec<float> *flow=new GVec<float>[n]; // flow in network
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
							if(path[i] && transfrag[t]->nodes.Last()!=gno-1) fragno+=flow[n1][n2];
							flow[n1][n2]=0;
						}
						else {
							if(!i) sumout+=transfrag[t]->abundance;
							flow[n1][n2]-=transfrag[t]->abundance;
							if(path[i] && transfrag[t]->nodes.Last()!=gno-1) fragno+=transfrag[t]->abundance;
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


float max_flow_EM(int gno,GVec<int>& path,GBitVec& istranscript,GPVec<CTransfrag>& transfrag,GPVec<CGraphnode>& no2gnode,
		GVec<float>& nodecapacity,GBitVec& pathpat,float &fragno) {



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
					if(path[i] && transfrag[t]->nodes.Last()!=gno-1) fragno+=*abund;
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
		GVec<float>& nodecapacity,GBitVec& pathpat,float& fragno) {

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
							if(path[i] && transfrag[t]->nodes.Last()!=gno-1)
								fragno+=flown1n2;
							flow[n1][n2]=0;
						}
						else {
							if(!i) sumout+=transfrag[t]->abundance;
							flow[n1][n2]-=transfrag[t]->abundance*rate[n1][n2];
							if(path[i] && transfrag[t]->nodes.Last()!=gno-1)
								fragno+=transfrag[t]->abundance;
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
		GPVec<CGraphnode>& no2gnode,int& geneno,bool& first,int strand,int gno,bool& included,GBitVec& prevpath,float fragno,char* id=NULL) {

	float cov=0;
	int len=0;
	CGraphnode *prevnode=NULL;
	GVec<GSeg> exons;
	GVec<float> exoncov;
	float excov=0;

	int s=0;
	if(!path[0]) s=1;

	for(int i=s;i<path.Count()-1;i++) {
		if(!prevpath[path[i]]) {
			included=false;
			prevpath[path[i]]=1;
		}
		if(i && !prevpath[edge(path[i-1],path[i],gno)]) {
			included=false;
			prevpath[edge(path[i-1],path[i],gno)]=1;
		}

		CGraphnode *node=no2gnode[path[i]];
		if(!prevnode || node->start>prevnode->end+1) { // this is a new exon
			if(prevnode) { // compute exon coverage
				excov/=exons.Last().end-exons.Last().start+1;
				exoncov.Add(excov);
				excov=0;
			}
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
		excov+=usedcov;
		if(node->cov) fragno+=node->frag*usedcov/node->cov;
		prevnode=node;
	}

	//fprintf(stderr,"Predicted transcript cov=%f usedcov=%f len=%d included=%d path.count=%d\n",cov/len, cov,len,included,path.Count());

	// Add last exon coverage
	if(prevnode) { // compute exon coverage
		excov/=exons.Last().end-exons.Last().start+1;
		exoncov.Add(excov);
	}
	if(len) cov/=len;

	if(id || ((!included || sensitivitylevel) && cov>=readthr && len>=mintranscriptlen)) { // store transcript here
	//if(id || ( cov>=readthr && len>=mintranscriptlen)) { // store transcript here
		char sign='-';
		if(strand) { sign='+';}
		if(first) { geneno++;}
		CPrediction *p=new CPrediction(geneno-1,id,exons[0].start,exons.Last().end,cov,sign,fragno,len);
		p->exons=exons;
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


void parse_trf(int maxi,int gno,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,
		GVec<bool>& compatible,	int& geneno,bool first,int strand,GList<CPrediction>& pred,GVec<float>& nodecov,
		GBitVec& istranscript,GBitVec& removable,GBitVec& usednode,float maxcov,GBitVec& prevpath,bool fast) {

	 GVec<int> path;
	 GVec<float> pathincov;
	 GVec<float> pathoutcov;
	 path.Add(maxi);
	 pathincov.cAdd(0.0);
	 pathoutcov.cAdd(0.0);
	 GBitVec pathpat(1+gno*(gno+1)/2);
	 pathpat[maxi]=1;
	 istranscript.reset();
	 GHash<CComponent> computed;

	 float flux=0;
	 float fragno=0;
	 GVec<float> nodeflux;

	 /*
	 { // DEBUG ONLY
	 	 fprintf(stderr,"start parse_trf with maxi=%d\n",maxi);
		 //fprintf(stderr,"Transcripts before path:");
		 //for(int i=0;i<transfrag.Count();i++) if(istranscript[i]) fprintf(stderr," %d",i);
		 //fprintf(stderr,"\n");
	 }
	 */

	 //if(back_to_source_path(maxi,path,pathpat,pathincov,pathoutcov,istranscript,removable,transfrag,computed,compatible,no2gnode,nodecov,gno)) {
	 if((fast && back_to_source_fast(maxi,path,pathpat,transfrag,no2gnode,nodecov,gno)) ||
			 (!fast && back_to_source_path(maxi,path,pathpat,pathincov,pathoutcov,istranscript,removable,transfrag,computed,compatible,no2gnode,nodecov,gno))) {
		 	 if(includesource) path.cAdd(0);
	 		 path.Reverse(); // back to source adds the nodes at the end to avoid pushing the list all the time
	 		 pathincov.Reverse();
	 		 pathoutcov.Reverse();

	 		 //if(fwd_to_sink_path(maxi,path,pathpat,pathincov,pathoutcov,istranscript,removable,transfrag,computed,compatible,no2gnode,nodecov,gno)) {
	 		if((fast && fwd_to_sink_fast(maxi,path,pathpat,transfrag,no2gnode,nodecov,gno)) ||
	 				(!fast && fwd_to_sink_path(maxi,path,pathpat,pathincov,pathoutcov,istranscript,removable,transfrag,computed,compatible,no2gnode,nodecov,gno))) {
	 			 pathincov.Clear();
	 			 pathoutcov.Clear();

	 			 //removable.reset();
	 			 //flux=update_flux(gno,path,istranscript,transfrag,removable,no2gnode,nodeflux,pathpat);


	 			 if(EM) flux=max_flow_EM(gno,path,istranscript,transfrag,no2gnode,nodeflux,pathpat,fragno);
	 			 else if(weight)
	 				 	 //flux=weight_max_flow_EM(gno,path,istranscript,transfrag,no2gnode,nodeflux,pathpat);
	 				 	 flux=weight_max_flow(gno,path,istranscript,transfrag,no2gnode,nodeflux,pathpat,fragno);
	 			 else flux=max_flow(gno,path,istranscript,transfrag,no2gnode,nodeflux,pathpat,fragno);

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
		 float cov=store_transcript(pred,path,nodeflux,nodecov,no2gnode,geneno,first,strand,gno,included,prevpath,fragno);

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
		 parse_trf(maxi,gno,no2gnode,transfrag,compatible,geneno,first,strand,pred,nodecov,istranscript,removable,usednode,maxcov,prevpath,fast);
	 }

}


CTransfrag *find_guide_pat(GffObj *guide,GPVec<CGraphnode>& no2gnode,int gno) {

	CTransfrag *trguide=NULL;

	int i=1;
	while(i<gno-1) {
		if(no2gnode[i]->overlap(guide)) break;
		i++;
	}

	if(i<gno-1) { // found start node
		GBitVec guidepat(1+gno*(gno+1)/2);
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
				guidepat[edge(i,inode->child[k],gno)]=1;
				i=inode->child[k];
				nodes.Add(i);
				inode=no2gnode[i];
			}
			else {
				if(i<gno-1 && (no2gnode[i+1]->start==inode->end+1)) { // !!!  maybe here I need to consider bundledist
					if(inode->child.Count() && inode->child[0]==i+1) {
						guidepat[i+1]=1;
						guidepat[edge(i,i+1,gno)]=1;
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

void process_refguides(int gno,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,int s,GPVec<GffObj>& guides,
		GVec<CGuide>& guidetrf) {

	char strand='-';
	if(s) strand='+';

	// find guides' patterns

	for(int g=0;g<guides.Count();g++) {
		if((guides[g]->strand==strand) && (guides[g]->overlap(no2gnode[1]->start,no2gnode[gno-2]->end))) {
			CTransfrag *trguide=find_guide_pat(guides[g],no2gnode,gno);
			if(trguide) { // the guide can be found among the graph nodes
				CGuide newguide(trguide,guides[g]->getID());
				guidetrf.Add(newguide);

				//fprintf(stderr,"Added guidetrf with ID=%s overlapping transcript interval %d - %d\n",guides[g]->getID())

				/*
				{ // DEBUG ONLY
					fprintf(stderr,"s=%d strand = %c trguide[%d]=",s,strand,g);
					printBitVec(trguide->pattern);
					fprintf(stderr,"\n");
				}
				*/
			}
		}
	}

	// compute guides' abundances
	for(int t=0;t<transfrag.Count();t++)
		for(int g=0;g<guidetrf.Count();g++) {
			if(((transfrag[t]->pattern) & guidetrf[g].trf->pattern) == transfrag[t]->pattern) {
				guidetrf[g].trf->abundance+=transfrag[t]->abundance;
			}
		}

	/*
	{ // DEBUG ONLY
		for(int g=0;g<guidetrf.Count();g++) {
			fprintf(stderr,"Abundance of guide[%d]=%f with nodes:",g,guidetrf[g].abundance);
			for(int i=0;i<guidetrf[g].nodes.Count();i++) fprintf(stderr," %d",guidetrf[g].nodes[i]);
			fprintf(stderr,"\n");
		}
	}
	*/

	guidetrf.Sort(guidedabundCmp);

	int g=0;
	while(g<guidetrf.Count()) {
		// check if guide is included in previous guide paths
		int p=0;
		bool included=false;
		if(!complete) { // if guides are incomplete exclude the ones that are included into the more complete ones
			while(p<g) {
				//CTransfrag guideg=guidetrf[g];
				//CTransfrag guidep=guidetrf[p];
				if((guidetrf[g].trf->pattern & guidetrf[p].trf->pattern)==guidetrf[g].trf->pattern) {
					included=true;
					GFREE(guidetrf[g].trf);
					guidetrf.Delete(g);
					break;
				}
				p++;
			}
		}

		if(!included) {
			// find if first node of guidetrf extends to source
			bool sourcestart=true;
			int lasti=guidetrf[g].trf->nodes[0];
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
						guidetrf[g].trf->pattern[edge(extendpath[i],j,gno)]=1;
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

					// extend guide to maxnode
					if(maxnode<j) {
						int i=0;
						GVec<int> tmpextend;
						while(i<extendpath.Count() && maxnode<=extendpath[i]) {
							tmpextend.Add(extendpath[i]);
							guidetrf[g].trf->pattern[extendpath[i]]=1;
							guidetrf[g].trf->pattern[edge(extendpath[i],j,gno)]=1;
							j=extendpath[i];
							i++;
						}
						if(includesource) tmpextend.cAdd(0); // I need to comment this if I need path to include the source
						tmpextend.Reverse();
						tmpextend.Add(guidetrf[g].trf->nodes);
						guidetrf[g].trf->nodes.Clear();
						guidetrf[g].trf->nodes.Add(tmpextend);
						guidetrf[g].trf->pattern[0]=1;
						guidetrf[g].trf->pattern[edge(0,maxnode,gno)]=1;
					}
					else {
						if(includesource) guidetrf[g].trf->nodes.Insert(0,source); // I need to comment this if I need path to include the source
						guidetrf[g].trf->pattern[0]=1;
						guidetrf[g].trf->pattern[edge(0,maxnode,gno)]=1;
					}

					// add source to maxnode transfrag
					GVec<int> nodes;
					nodes.cAdd(0);
					nodes.Add(maxnode);
					GBitVec trpat(1+(gno+1)*gno/2);
					trpat[0]=1;
					trpat[maxnode]=1;
					trpat[edge(0,maxnode,gno)]=1;
					CTransfrag *tr=new CTransfrag(nodes,trpat,maxabund);
					transfrag.Add(tr);

					// add source among maxnode parents
					inode=no2gnode[maxnode];
					inode->parent.Insert(0,source);
					inode->parentpat[edge(0,maxnode,gno)]=1; // source should be already among parents of maxnode but not the edge to source
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
						guidetrf[g].trf->pattern[edge(j,extendpath[i],gno)]=1;
						j=extendpath[i];
					}
				}
				else { // if no path to sink was found I need to add one

					int sink=gno-1;

					// extend guide to maxnode
					if(maxnode>j) {
						int i=0;
						while(i<extendpath.Count() && maxnode>=extendpath[i]) {
							guidetrf[g].trf->nodes.Add(extendpath[i]);
							guidetrf[g].trf->pattern[extendpath[i]]=1;
							guidetrf[g].trf->pattern[edge(j,extendpath[i],gno)]=1;
							j=extendpath[i];
							i++;
						}
					}
					guidetrf[g].trf->nodes.Add(sink);

					// add maxnode to sink transfrag
					GVec<int> nodes;
					nodes.Add(maxnode);
					nodes.Add(sink);
					GBitVec trpat(1+(gno+1)*gno/2);
					trpat[maxnode]=1;
					trpat[sink]=1;
					trpat[edge(maxnode,sink,gno)]=1;
					CTransfrag *tr=new CTransfrag(nodes,trpat,maxabund);
					transfrag.Add(tr);

					// add sink among maxnode children
					inode=no2gnode[maxnode];
					inode->child.Add(sink);
					inode->childpat[edge(maxnode,sink,gno)]=1; // source should be already among parents of maxnode but not the edge to source
				}
			}
		}
		g++;
	}
}

int guides_flow(int gno,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,GVec<CGuide>& guidetrf,int& geneno,
		int s,GList<CPrediction>& pred,GVec<float>& nodecov,GBitVec& istranscript,GBitVec& pathpat) {

	int maxi=1;
	bool cov=false; // tells me if max node coverage was determined

	int g=0;
	bool first=true;
	while(g<guidetrf.Count()) {

		// weight the transcript
		GVec<float> nodeflux;
		float flux=0;
		float fragno=0;

		if(EM) flux= max_flow_EM(gno,guidetrf[g].trf->nodes,istranscript,transfrag,no2gnode,nodeflux,guidetrf[g].trf->pattern,fragno);
		else if(weight) flux= weight_max_flow(gno,guidetrf[g].trf->nodes,istranscript,transfrag,no2gnode,nodeflux,guidetrf[g].trf->pattern,fragno);
		else flux= max_flow(gno,guidetrf[g].trf->nodes,istranscript,transfrag,no2gnode,nodeflux,guidetrf[g].trf->pattern,fragno);

		istranscript.reset();
		/*
		{ // DEBUG ONLY
			fprintf(stderr,"flux[%d]=%g Path:",g,flux);for(int i=0;i<path.Count();i++) fprintf(stderr," %d",path[i]);
			fprintf(stderr,"\n");
		}
	    */

		if(flux>epsilon) {
			bool include=true;
			char *predid=NULL;
			int idlen=strlen(guidetrf[g].id);
			if(idlen) {
				predid=Gstrdup(guidetrf[g].id);
			}
			store_transcript(pred,guidetrf[g].trf->nodes,nodeflux,nodecov,no2gnode,geneno,first,s,gno,include,pathpat,fragno,predid);
			cov=true;
			// Node coverages:
			for(int i=1;i<gno-1;i++)
				if(nodecov[i]>nodecov[maxi]) maxi=i;

			//if(nodecov[maxi]<readthr) break; // no need to find other paths since they are under allowed read threshold
			if(nodecov[maxi]<1) break; // no need to find other paths since I don't even have one read/bp coverage; it's more strict here

			/*
			{ // DEBUG ONLY
				fprintf(stderr,"\nAfter update:\n");
				for(int i=1;i<gno;i++) {
					fprintf(stderr,"Node %d: %f ",i,nodecov[i]);
					fprintf(stderr,"trf=");
					for(int t=0;t<no2gnode[i]->trf.Count();t++) fprintf(stderr," %d(%f)",no2gnode[i]->trf[t],transfrag[no2gnode[i]->trf[t]]->abundance);
					fprintf(stderr," maxi=%d maxcov=%f\n",maxi,nodecov[maxi]);
				}
			}
			*/

			nodeflux.Clear();
		}
		g++;
	} // end while(g<guidetrf.Count()

	if(!cov) for(int i=2;i<gno-1;i++)
		if(nodecov[i]>nodecov[maxi]) maxi=i;

	return(maxi);
}

int find_transcripts(int gno,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,GVec<bool>& compatible,
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

	if(nodecov[maxi]>=readthr) {
		GBitVec istranscript(transfrag.Count());
		GBitVec pathpat(1+gno*(gno+1)/2);
		// process guides first
		if(guidetrf.Count()) maxi=guides_flow(gno,no2gnode,transfrag,guidetrf,geneno,strand,pred,nodecov,istranscript,pathpat);

		// process rest of the transfrags
		if(!eonly && nodecov[maxi]>=readthr) {
			GBitVec removable(transfrag.Count(),true);

			// 1:
			// parse_trf_weight_max_flow(gno,no2gnode,transfrag,geneno,strand,pred,nodecov,pathpat);
			// 2:
			GBitVec usednode(1+gno*(gno+1)/2);
			parse_trf(maxi,gno,no2gnode,transfrag,compatible,geneno,true,strand,pred,nodecov,istranscript,removable,usednode,0,pathpat,fast);

		}
	}

	return(geneno);
}

void get_trims(GVec<CTrimPoint>& trims,CBundlenode *currbnode,int refstart,GVec<float>& bpcov) {

	uint sourcestart;
	uint sinkend;
	float sinkabundance;
	float sourceabundance;

	while(currbnode!=NULL) {
		sourcestart=0;
		sinkend=0;
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

		if(covered && maxguidelen>=(uint)mintranscriptlen) { if(c_out) guide->printTranscriptGff(c_out);}
		else covered=false;
	}

	return(covered);
}

void clean_read_junctions(CReadAln *read) {
	for(int i=0;i<read->juncs.Count();i++) {
		CJunction& jd=*(read->juncs[i]);
		if(jd.strand) { // junction wasn't deleted -> it might need to be deleted now
			int j=0;
			while(j<read->segs.Count() && read->segs[j].end<jd.start) j++;
			if(j<read->segs.Count()-1 && read->segs[j].end==jd.start && read->segs[j+1].start==jd.end) {
				if ((int)read->segs[j].len() >= junctionsupport && (int)read->segs[j+1].len() >=junctionsupport) {
					jd.nreads_good-=float(1)/read->nh;
					if(jd.nreads_good<junctionthr) jd.strand=0;
				}
			}
		}
	}
}

int build_graphs(int refstart, GList<CReadAln>& readlist,
		GList<CJunction>& junction, GPVec<GffObj>& guides, GVec<float>& bpcov, GList<CPrediction>& pred,bool fast) {

	// form groups on strands: all groups below are like this: 0 = negative strand; 1 = unknown strand; 2 = positive strand
	GPVec<CGroup> group;
	CGroup *currgroup[3]={NULL,NULL,NULL}; // current group of each type
	CGroup *startgroup[3]={NULL,NULL,NULL}; // start group of each type
	int color=0; // next color to assign
	GVec<int> merge; // remembers merged groups
	GVec<int> equalcolor; // remembers colors for the same bundle
	GVec<int> *readgroup=new GVec<int>[readlist.Count()]; // remebers groups for each read; don't forget to delete it when no longer needed

	//int **readgroup = new int*[readlist.Count()];

	float fraglen=0;
	uint fragno=0;

	for (int n=0;n<readlist.Count();n++) {
		CReadAln & rd=*(readlist[n]);

		//if(rd.juncs.Count()) fprintf(stderr,"read[%d] %d %d %d [", n, rd.start, rd.end, (int)rd.strand);

		bool keep=true;
		int i=0;
		while(i<rd.juncs.Count()) {
			CJunction& jd=*(rd.juncs[i]);
			//fprintf(stderr, " %d-%d:%4.2f", jd.start, jd.end, jd.nreads_good);
			if(!jd.strand) {
				keep=false;
				rd.nh=0;
				if(rd.pair_idx>-1)
					readlist[rd.pair_idx]->pair_idx=-1;
				break;
			}
			i++;
		}
		//if(rd.juncs.Count()) fprintf(stderr,"] keep=%d\n",keep);
		if(keep) {
			color=add_read_to_group(n,readlist,color,group,currgroup,startgroup,readgroup,equalcolor,merge);

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
		}
		//else clean_read_junctions(readlist[n]);
	}

	if(fragno) fraglen/=fragno;

	// merge groups that are close together
	if(bundledist) {
		for(int sno=0;sno<3;sno++) {
			CGroup *lastgroup=NULL;
			CGroup *procgroup=startgroup[sno];
			while(procgroup!=NULL) {
				if(lastgroup && procgroup->start-lastgroup->end<=bundledist) {
					merge_fwd_groups(group,lastgroup,procgroup,merge,equalcolor);
					procgroup=lastgroup->next_gr;
					continue;
				}
				lastgroup=procgroup;
				procgroup=procgroup->next_gr;
			}
		}
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

	while(currgroup[0]!=NULL || currgroup[1]!=NULL || currgroup[2]!=NULL) { // there are still groups to process

		int nextgr=get_min_start(currgroup); // gets the index of currgroup with the left most begining

		int grcol = currgroup[nextgr]->color;    // set smallest color for currgroup[$nextgr]

		while(equalcolor[grcol]!=grcol) {
			grcol=equalcolor[grcol];
		}
		currgroup[nextgr]->color=grcol;

		//print STDERR "nextgr=$nextgr grcol=$grcol current group is at coords: ",$currgroup[$nextgr][0],"-",$currgroup[$nextgr][1],"\n";

		if(nextgr == 0) { // negative strand group
			if(prevgroup[1]!=NULL && currgroup[nextgr]->start < prevgroup[1]->end) { // overlaps an unknown strand group
				set_strandcol(prevgroup[1],currgroup[nextgr],grcol,eqnegcol,equalcolor);
			}
		}
		else if(nextgr == 2) { // positive strand group
			if(prevgroup[1]!=NULL && currgroup[nextgr]->start < prevgroup[1]->end) { // overlaps an unknown strand group
				set_strandcol(prevgroup[1],currgroup[nextgr],grcol,eqposcol,equalcolor);
			}
		}
		else { // unknown strand group
			if(prevgroup[0]!=NULL && currgroup[nextgr]->start < prevgroup[0]->end) { // overlaps negative strand group
				set_strandcol(currgroup[nextgr],prevgroup[0],prevgroup[0]->color,eqnegcol,equalcolor);
			}
			if(prevgroup[2]!=NULL && currgroup[nextgr]->start < prevgroup[2]->end) { // overlaps positive strand group
				set_strandcol(currgroup[nextgr],prevgroup[2],prevgroup[2]->color,eqposcol,equalcolor);
			}
		}

		prevgroup[nextgr]=currgroup[nextgr];
		currgroup[nextgr]=currgroup[nextgr]->next_gr;

    }

    //print STDERR "Colors assigned!\n";

	// create bundles
	for (int i=0;i<3;i++) {
		currgroup[i]=startgroup[i];
		prevgroup[i]=NULL;
	}

	GPVec<CBundle> bundle[3]; // all bundles on all strands: 0,1,2
	GPVec<CBundlenode> bnode[3]; // last bnodes on all strands: 0,1,2 for each bundle

	GVec<int> group2bundle[3]; // to retrace reads from group no to bundle
	for(int sno=0;sno<3;sno++) {
		group2bundle[sno].Resize(group.Count(),-1);
		bnode[sno].setFreeItem(false);
	}

	GVec<int> bundlecol(true); // associates a bundle number to a group color
	bundlecol.Resize(equalcolor.Count(),-1);

	while(currgroup[0]!=NULL || currgroup[1]!=NULL || currgroup[2]!=NULL) { // there are still groups to process

		int nextgr=get_min_start(currgroup);

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
		else { // unknown strand

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

	// next variables are in order to remember if I get guide coverage
	GVec<int>* bnodeguides=NULL;
	if(bundle[1].Count() && bnode[1].Count()) {
		bnodeguides = new GVec<int>[bnode[1].Count()];
	}

	//if(guides.Count()) fprintf(stderr,"No of guides=%d\n",guides.Count());

	if(partialcov) {
		for(int g=0;g<guides.Count();g++) {
			int s=0;
			if(guides[g]->strand=='+') s=2;
			get_partial_covered(guides[g],bundle[s],bnode[s],junction);
		}
		return(0);
	}

	if(c_out || (bundle[1].Count() && bnode[1].Count())) // coverage is needed
		for(int g=0;g<guides.Count();g++) {
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
	if(fraglen) for(int b=0;b<bundle[1].Count();b++) {

    	if(bundle[1][b]->nread && (bundle[1][b]->multi/bundle[1][b]->nread)<=mcov && bundle[1][b]->len > mintranscriptlen) { // there might be small transfrags that are worth showing, but here I am ignoring them

    		// bundle might contain multiple fragments of a transcript but since we don't know the complete structure -> print only the pieces that are well represented
    		CBundlenode *currbnode=bnode[1][bundle[1][b]->startnode];
    		int t=1;
    		while(currbnode!=NULL) {
    			int len=currbnode->end-currbnode->start+1;
    			float cov=currbnode->cov/(currbnode->end-currbnode->start+1);
    			GStr id;
    			char sign=0;

    			for(int i=0;i<bnodeguides[currbnode->bid].Count();i++) {
    				int g=bnodeguides[currbnode->bid][i];
    				if(!id.is_empty()) id+=';';
    				id+=guides[g]->getID();
    				if(!sign) sign=guides[g]->strand;
    				else if(sign!='.' && guides[g]->strand!=sign) sign='.';
    			}

    			char *predid=NULL;
    			if(!id.is_empty()) {
    				predid=Gstrdup(id.chars());
    				//predid=new char[id.length()+1];
    				//strcpy(predid,id.chars());
    			}
    			if((!eonly && cov>=readthr && len>mintranscriptlen) || (!id.is_empty())) {
    				if(t==1) { geneno++;}
    				if(!sign) sign='.';
    				CPrediction *p=new CPrediction(geneno-1,predid,currbnode->start,currbnode->end,cov,sign,cov,fraglen);
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

    //print STDERR "Done with unstranded bundles\n";
    if (bnodeguides) delete[] bnodeguides;

	// ### build graphs for stranded bundles here
    if(startgroup[0]!=NULL || startgroup[2]!=NULL) { //# there are stranded groups to process

    	// I don't need the groups here anymore : I only use their numbers
    	group.Clear();

    	// sort junctions -> junctions are sorted already according with their start, but not their end
    	GList<CJunction> ejunction(junction);
    	ejunction.setFreeItem(false);
    	if(ejunction.Count()) ejunction.setSorted(juncCmpEnd);

    	GVec<CGraphinfo> *bundle2graph[2]; // should I keep the neutral strand for consistency ? -> remember not to delete it
    	GVec<int> graphno[2];  // how many nodes are in a certain graph g, on strand s: graphno[s][g]
    	// GVec<int> trnumber[2]; // how many transfrags are on a strand s, in a graph g -> I can find out this from transfrag[s][g].Count()
    	// int ngraph[2]={0,0};   // how many graphs are in each strand: negative (0), or positive(1) -> keep one for each bundle
    	GPVec<CTransfrag> *transfrag[2]; // for each transfrag t on a strand s, in a graph g, transfrag[s][g][t] gives it's abundance and it's pattern
    	GPVec<CGraphnode> *no2gnode[2]; // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i
    	CTreePat **tr2no[2]; // for each graph g, on a strand s, tr2no[s][g] keeps the tree pattern structure for quick retrieval of the index t of a tansfrag

    	int bno[2]={0,0};

    	// build graph structure
    	for(int sno=0;sno<3;sno+=2) { // skip neutral bundles -> those shouldn't have junctions

    		int s=sno/2; // adjusted strand due to ignoring neutral strand
    		bundle2graph[s]=NULL;
    		if(bnode[sno].Count()) bundle2graph[s]=new GVec<CGraphinfo>[bnode[sno].Count()];
    		transfrag[s]=NULL;
    		no2gnode[s]=NULL;
    		tr2no[s]=NULL;

    		if(bundle[sno].Count()) {
    			transfrag[s]=new GPVec<CTransfrag>[bundle[sno].Count()]; // for each bundle I have a graph ? only if I don't ignore the short bundles
    			no2gnode[s]=new GPVec<CGraphnode>[bundle[sno].Count()];
    			GCALLOC(tr2no[s],bundle[sno].Count()*sizeof(CTreePat *));
    			bno[s]=bundle[sno].Count();

    			for(int b=0;b<bundle[sno].Count();b++) {
    				graphno[s].cAdd(0);

    				//if(bundle[sno][b]->nread) fprintf(stderr,"proc %f/%f is %f\n",bundle[sno][b]->multi,bundle[sno][b]->nread,(float)bundle[sno][b]->multi/bundle[sno][b]->nread);

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
    							bundle2graph,no2gnode,transfrag,bpcov); // also I need to remember graph coverages somewhere -> probably in the create_graph procedure

    					if(graphno[s][b]) tr2no[s][b]=construct_treepat(graphno[s][b],transfrag[s][b]);
    					else tr2no[s][b]=NULL;
    				}
    				else tr2no[s][b]=NULL;
    			}
    		}
    	}
    	//fprintf(stderr,"Done creating graphs\n");

		// I can clean up some data here:
    	for(int sno=0;sno<3;sno++) {
    		int n=bnode[sno].Count();
    		for(int b=0;b<n;b++) delete bnode[sno][b];
    		bnode[sno].Clear();
    		bundle[sno].Clear();
    	}


    	// compute probabilities for stranded bundles
    	for (int n=0;n<readlist.Count();n++) {
    		int np=readlist[n]->pair_idx;
    		if(np==-1 || n<np) // read might have been deleted -> check this in get_fragment_pattern
    			get_fragment_pattern(readlist,n,np,readgroup,merge,group2bundle,bundle2graph,graphno,no2gnode,transfrag,tr2no);
    	}
    	//fprintf(stderr, "Done read parsing and updating frequencies\n");

    	/*
    	{ // DEBUG ONLY
    		printTime(stderr);
    		for(int s=0;s<2;s++) {
    			fprintf(stderr, "There are %d stranded[%d]\n",bno[s],int(2*s));
    			for(int b=0;b<bno[s];b++) {
    				if(graphno[s][b]) {
    					GStr pat;
    					fprintf(stderr,"Graph[%d][%d] with %d nodes:\n",int(2*s),b,graphno[s][b]);
    					//print_pattern(tr2no[s][b],pat,graphno[s][b]);
    				}
    			}
    		}
    	}
    	*/

    	// don't forget to clean up the allocated data here
    	delete [] readgroup;

    	// parse graph
    	for(int s=0;s<2;s++) {

    		for(int b=0;b<bno[s];b++) {
    			if(graphno[s][b]) {

    				// include source to guide starts links
    				GVec<CGuide> guidetrf;
    				if(guides.Count()) process_refguides(graphno[s][b],no2gnode[s][b],transfrag[s][b],s,guides,guidetrf);

    				//process transfrags to eliminate noise, and set compatibilities, and node memberships
    				GVec<bool> compatible; // I might want to change this to gbitvec
    				process_transfrags(graphno[s][b],no2gnode[s][b],transfrag[s][b],tr2no[s][b],compatible);

    				/*
    				{ //DEBUG ONLY
    					printTime(stderr);
    					fprintf(stderr,"There are %d nodes:\n",graphno[s][b]);
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

    				// find transcripts now
    				geneno=find_transcripts(graphno[s][b],no2gnode[s][b],transfrag[s][b],compatible,
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


    			}
    			// clean up what can be cleaned
    			if(tr2no[s][b]) free_treepat(tr2no[s][b]);
    		}

    		// final clean up: no2gnode, no2tr, transfrag, bundle2graph
    		if(bundle2graph[s]) delete [] bundle2graph[s];
    		if(transfrag[s]) delete [] transfrag[s];
    		if(no2gnode[s]) delete [] no2gnode[s];
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

    // don't forget to clean up the allocated data here
    return(geneno);
}


void clean_junctions(GList<CJunction>& junction) {

	for(int i=0;i<junction.Count();i++) {
		CJunction& jd=*(junction[i]);
		if(jd.nreads_good<junctionthr) {
			jd.strand=0;
		}
	}
}


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

int infer_transcripts(int refstart, GList<CReadAln>& readlist,
		GList<CJunction>& junction, GPVec<GffObj>& guides, GVec<float>& bpcov, GList<CPrediction>& pred, bool fast) {

	//DEBUG ONLY: 	showReads(refname, readlist);

	clean_junctions(junction);
	int geneno=build_graphs(refstart, readlist, junction, guides, bpcov, pred, fast);

	return(geneno);
}


void printGff3Header(FILE* f, GArgs& args) {
  fprintf(f, "# ");
  args.printCmdLine(f);
  fprintf(f, "##gff-version 3\n");
}


