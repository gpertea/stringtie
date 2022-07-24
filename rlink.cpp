#include "rlink.h"
#include "GBitVec.h"
#include <float.h>
#include "GThreads.h"

//#define GMEMTRACE 1  //debugging memory allocation
#ifdef GMEMTRACE
#include "proc_mem.h"
#endif

#define BSIZE 10000 // bundle size

//import globals from main program:

//extern GffNames* gseqNames;
extern FILE *c_out;         // file handle for the input transcripts that are fully covered by reads
extern GFastMutex printCovMutex;

extern bool trim;
extern bool eonly;
extern bool nomulti;
extern bool viral;
extern bool mixedMode;
extern bool guided;

extern int allowed_nodes;
extern float isofrac;
extern bool isunitig;
extern bool longreads;
extern bool rawreads;
extern float mcov;
extern int mintranscriptlen; // minimum number for a transcript to be printed
extern uint junctionsupport; // anchor length for junction to be considered well supported <- consider shorter??
extern uint sserror;
extern int junctionthr; // number of reads needed to support a particular junction
extern float readthr;     // read coverage per bundle bp to accept it; otherwise considered noise
extern float singlethr;     // read coverage per bundle bp to accept it; otherwise considered noise
extern uint bundledist;  // reads at what distance should be considered part of separate bundles
                        // <- this is not addressed everywhere, e.g. in infer_transcripts -> look into this

extern bool includesource;
extern bool geneabundance; // need to compute the gene abundance

extern float fpkm_thr;
extern float tpm_thr;
extern bool enableNames;
extern bool includecov;
extern bool retained_intron;

extern FILE* f_out;
extern GStr label;

static GStr _id("", 256); //to prevent repeated reallocation for each parsed read
//not thread safe -- to only be used in processRead() as long as that's the unique producer

CJunction* add_junction(int start, int end, GList<CJunction>& junction, char strand) {

	int oidx=-1;
	CJunction *nj=NULL;

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



void cov_edge_add(GVec<float> *bpcov, int sno, int start, int end, float v) {
	bool neutral=false;
	if(sno!=1) neutral=true; // why is neutral true here: because if the sno is -/+ than I want to add their counts to bpcov[1] too
	bpcov[sno][start+1]+=v; // if sno==1 then I add v to it here
	bpcov[sno][end+1]-=v;
	if(neutral) { // if neutral (i.e. stranded) gets added to bpcov[1] here too => bpcov[1]=all coverage
		bpcov[1][start+1]+=v;
		bpcov[1][end+1]-=v;
	}
}

void add_read_to_cov(GList<CReadAln>& rd,int n,GVec<float> *bpcov,int refstart) {

	int sno=(int)rd[n]->strand+1; // 0(-),1(.),2(+)

	float single_count=rd[n]->read_count;
	for(int i=0;i<rd[n]->pair_idx.Count();i++) {
		single_count-=rd[n]->pair_count[i];
		if(rd[n]->pair_idx[i]<=n) { // pair comes before
			int np=rd[n]->pair_idx[i];
			int pcount=rd[np]->segs.Count();
			int rcount=rd[n]->segs.Count();
			int snop=(int)rd[np]->strand+1; //v6
			int p=0;
			int r=0;
			int nsegs=pcount+rcount;
			while(nsegs) {
				int start;
				int end;
				if(r<rcount) {
					if(p==pcount || rd[np]->segs[p].start>rd[n]->segs[r].end) { // no more pair segments or pair segment comes after read segment
						start=rd[n]->segs[r].start;
						end=rd[n]->segs[r].end;
						r++;
						nsegs--;
					}
					else if(rd[np]->segs[p].end<rd[n]->segs[r].start) {
						start=rd[np]->segs[p].start;
						end=rd[np]->segs[p].end;
						p++;
						nsegs--;
					}
					else { // segments overlap
						start=rd[np]->segs[p].start;
						end=rd[np]->segs[p].end;
						if((int)rd[n]->segs[r].start<start) start=rd[n]->segs[r].start;
						p++;
						nsegs--;
						bool cont=true;
						while(cont) {
							cont=false;
							if(r<rcount && (int)(rd[n]->segs[r].start)<=end) {
								if((int)rd[n]->segs[r].end>end) end=rd[n]->segs[r].end;
								r++;
								nsegs--;
								cont=true;
							}
							if(p<pcount && (int)(rd[np]->segs[p].start)<=end) {
								if((int)rd[np]->segs[p].end>end) end=rd[np]->segs[p].end;
								p++;
								nsegs--;
								cont=true;
							}
						}
					}
				}
				else {
					start=rd[np]->segs[p].start;
					end=rd[np]->segs[p].end;
					p++;
					nsegs--;
				}
				int strand=sno;
				if(sno!=snop) {
					if(sno==1) strand=snop;
					else if(snop!=1) strand=1;
				}
				cov_edge_add(bpcov,strand,start-refstart,end-refstart+1,rd[n]->pair_count[i]);
			}
		}
	}
	if(single_count>epsilon) {
		for(int i=0;i<rd[n]->segs.Count();i++) {
			cov_edge_add(bpcov,sno,rd[n]->segs[i].start-refstart,rd[n]->segs[i].end-refstart+1,single_count);
		}
	}
}

void countFragment(BundleData& bdata, GSamRecord& brec, int nh) {
	static uint32_t BAM_R2SINGLE = BAM_FREAD2 | BAM_FMUNMAP ;


	for (int i=0;i<brec.exons.Count();i++) {
		bdata.frag_len+=float(1)*brec.exons[i].len()/nh;
	}
	if (!brec.isPaired() || ((brec.flags()&BAM_FREAD1)!=0) ||
					((brec.flags()&BAM_R2SINGLE)==BAM_R2SINGLE ) ) {
		bdata.num_fragments+=float(1)/nh;
	}

}

bool deljuncmatch(CReadAln *rd,GVec<GSeg> &jd) {

	if(rd->juncs.Count()!=jd.Count()) return(false);

	for(int i=0;i<rd->juncs.Count();i++) {
		uint startj=rd->segs[i].end-rd->juncs[i]->start;
		uint endj=rd->juncs[i]->end-rd->segs[i+1].start;
		if(startj!=jd[i].start || endj!=jd[i].end) return(false);
	}

	return(true);
}

bool exonmatch(GVec<GSeg> &prevexons, GVec<GSeg> &exons) {
	if(prevexons.Count() != exons.Count()) return false;
	for(int i=0;i<exons.Count();i++) {
		if(prevexons[i].end!=exons[i].end || prevexons[i].start!=exons[i].start) return false;
	}
	return true;
}

bool mismatch_anchor(CReadAln *rd,char *mdstr,int refstart, bam1_t *b) {
	if(mdstr==NULL) return false;
	//--make a copy of the string, in case the original is a const string
	// (because the parseUInt() function modifies the string temporarily
	char* mdstring=Gstrdup(mdstr);
	char *p=mdstring;

	int i=0;
	int parsedlen=refstart;
	int rdlen=0;

	while (*p!='\0') {
		unsigned int num_matches=0;
		if (*p>='0' && *p<='9') {
			parseUInt(p,num_matches);
			//if (num_matches>0) GMessage("%d matching bases\n", num_matches);
			parsedlen+=num_matches;
			continue;
		}
		if (*p=='^') { //deletion --> found a problem with deletion
			//GDynArray<char> deletion; //deletion string accumulates here (if needed)
			int del_length=0;//tracking deletion length
			char delbase=*(++p);
			while (delbase>='A' && delbase<='Z') {
				//deletion.Add(delbase);
				del_length++;
				delbase=*(++p);
			}
			while(i<rd->segs.Count() && rdlen+(int)rd->segs[i].len()<parsedlen) {
				rdlen+=rd->segs[i].len();
				i++;
			}
			if(i==rd->segs.Count()) break;

			if((i && parsedlen-rd->segs[i].start<junctionsupport) || (i<rd->segs.Count()-1 && rd->segs[i].end+1-parsedlen-del_length<junctionsupport)) {
				GFREE(mdstring);
				return true;
			}
			parsedlen+=del_length;

			/*GMessage("%d base(s) deletion [", del_length);
			for (uint i=0;i<deletion.Count();++i) GMessage("%c",deletion[i]);
			GMessage("]\n");*/
			continue;
		}
		if (*p>='A' && *p<='Z') {
			//GMessage("base mismatch [%c]\n",*p);
			while(i<rd->segs.Count() && rdlen+(int)rd->segs[i].len()<parsedlen) {
				rdlen+=rd->segs[i].len();
				i++;
			}
			if(i==rd->segs.Count()) break;
			if((i && parsedlen-rd->segs[i].start<junctionsupport) || (i<rd->segs.Count()-1 && rd->segs[i].end-parsedlen<junctionsupport)) {
				GFREE(mdstring);
				return true;
			}
			parsedlen++;
		}
		p++;
	}

	GFREE(mdstring);

	uint32_t *cigar = bam_get_cigar(b);
	rdlen=0;
	parsedlen=0;
	i=0;

	for (uint j = 0; j < b->core.n_cigar; ++j) {
		int op = bam_cigar_op(cigar[j]);
		if (op == BAM_CMATCH || op==BAM_CEQUAL ||
				op == BAM_CDIFF || op == BAM_CDEL) {
			parsedlen += bam_cigar_oplen(cigar[j]);
		}
		else if(op == BAM_CINS) {
			while(i<rd->segs.Count() && rdlen+(int)rd->segs[i].len()<parsedlen) {
				rdlen+=rd->segs[i].len();
				i++;
			}
			if(i==rd->segs.Count()) break;
			if((i && parsedlen-rd->segs[i].start<junctionsupport) || (i<rd->segs.Count()-1 && rd->segs[i].end-parsedlen<junctionsupport)) {
				return true;
			}
		}
	}

	return false;
}

void processRead(int currentstart, int currentend, BundleData& bdata,
		 GHash<int>& hashread,  GReadAlnData& alndata) { // some false positives should be eliminated here in order to break the bundle
	GSamRecord& brec=*(alndata.brec);			   // bam record
	if((longreads || (mixedMode && brec.uval)) && (brec.flags() & BAM_FSECONDARY)) return;
	GList<CReadAln>& readlist = bdata.readlist;    // list of reads gathered so far
	GList<CJunction>& junction = bdata.junction;   // junctions added so far
    char strand=alndata.strand;
    int nh=alndata.nh;
    int hi=alndata.hi;
	int readstart=brec.start;
	CReadAln* readaln=NULL;                        // readaln is initialized with NULL
	//bool covSaturated=false;                       // coverage is set to not saturated

	/*
	{ // DEBUG ONLY
		fprintf(stderr,"Process read %s with strand=%d and exons:",brec.name(),strand);
		for (int i=0;i<brec.exons.Count();i++) {
			fprintf(stderr," %d-%d", brec.exons[i].start, brec.exons[i].end);
		}
		fprintf(stderr,"\n");
	}
	*/

	double nm=(double)brec.tag_int("NM"); // read mismatch
	float unitig_cov=0;
	unitig_cov=brec.tag_float("YK");

	bool longr=false;
	if(longreads|| brec.uval) longr=true; // second alignment is always from mixed reads

	bool match=false;  // true if current read matches a previous read
	int n=readlist.Count()-1;

	if(!mergeMode) while(n>-1 && readlist[n]->start==brec.start) {
		if(strand==readlist[n]->strand && (readlist[n]->longread==longr) && (!isunitig || (unitig_cov>0) == readlist[n]->unitig)) {
			match=exonmatch(readlist[n]->segs,brec.exons);
			if(match && (longreads || mixedMode)) match=deljuncmatch(readlist[n],brec.juncsdel); //DEL AWARE
		}
		//if(strand==readlist[n]->strand) match=exonmatch(readlist[n]->segs,brec.exons);
		if(match) break; // this way I make sure that I keep the n of the matching readlist
		n--;
	}
	else if((alndata.tinfo->cov>=0 && alndata.tinfo->cov<readthr) || (alndata.tinfo->fpkm>=0 && alndata.tinfo->fpkm<fpkm_thr) ||
			(alndata.tinfo->tpm>=0 && alndata.tinfo->tpm < tpm_thr)) return; // do not store 'read' if it doesn't meet minimum criteria; this is mergeMode
	else { //mergeMode but the read is above thresholds
		if(alndata.tinfo->cov>=0) bdata.covflags |= IS_COV_FLAG;
		if(alndata.tinfo->fpkm>=0) bdata.covflags |= IS_FPKM_FLAG;
		if(alndata.tinfo->tpm>=0) bdata.covflags |= IS_TPM_FLAG;
	}

	if (bdata.end<currentend) {
		bdata.start=currentstart;
		bdata.end=currentend;
	}
	bdata.numreads++;                         // number of reads gets increased no matter what
	//bdata.wnumreads+=float(1)/nh;

	if (!match) { // if this is a new read I am seeing I need to set it up
		if(mergeMode && mintranscriptlen) {
			int len=0;
			for (int i=0;i<brec.exons.Count();i++) len+=brec.exons[i].len();
			if(len<mintranscriptlen) return;
		}
		readaln=new CReadAln(strand, nh, brec.start, brec.end, alndata.tinfo);
		readaln->longread=longr;
		alndata.tinfo=NULL; //alndata.tinfo was passed to CReadAln
		for (int i=0;i<brec.exons.Count();i++) {
			readaln->len+=brec.exons[i].len();
			if(i) {

				if(!junction.Count()) { // always add null junction first
					CJunction *nullj=new CJunction(0, 0, 0);
					junction.Add(nullj);
				}
				int jstrand=strand;
				uint jstart=brec.exons[i-1].end;
				uint jend=brec.exons[i].start;

				//fprintf(stderr,"exon count=%d junctiondel count=%d exonend=%d exonstart=%d\n",brec.exons.Count(),brec.juncsdel.Count(),jstart,jend);

				if((longreads || mixedMode) && (brec.juncsdel[i-1].start || brec.juncsdel[i-1].end)) { // deletion at the junction start/end //DEL AWARE
					if(!alndata.juncs.Count() || !alndata.juncs[i-1]->guide_match) { // if this junction matches a guide, I do not need to do anything
						jstrand=0;
						jstart-=brec.juncsdel[i-1].start;
						jend+=brec.juncsdel[i-1].end;
					}
				}
				CJunction* nj=junction.AddIfNew(new CJunction(jstart, jend, jstrand), true);
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

	/*
	{ // DEBUG ONLY
		fprintf(stderr,"Add read %s with strand=%d and exons:",brec.name(),strand);
		for (int i=0;i<brec.exons.Count();i++) {
			fprintf(stderr," %d-%d", brec.exons[i].start, brec.exons[i].end);
		}
		fprintf(stderr,"\n");
		//fprintf(stderr,"Read %s is at n=%d with unitig_cov=%f and strand=%d\n",brec.name(),n,unitig_cov,strand);
	}
	*/

	if((int)brec.end>currentend) {
			currentend=brec.end;
	  	bdata.end=currentend;
	}

	float rdcount=(float)brec.tag_int("YC"); // alignment count
	if(!rdcount) rdcount=1;
	if(unitig_cov) {
		rdcount=unitig_cov;
		if(isunitig) readlist[n]->unitig=true; // treat unitig reads differently when -U is set
	}

	if(!nomulti) rdcount/=nh;
	readlist[n]->read_count+=rdcount; // increase single count just in case I don't see the pair

	// store the mismatch count per junction so that I can eliminate it later
	if(!nm) {
		nm=(double)brec.tag_int("nM"); // paired mismatch : big problem with STAR alignments
		if(brec.isPaired()) nm/=2;
	}
	if(brec.clipL) nm++;
	if(brec.clipR) nm++;

	//nm+=brec.clipL; Note: this clippings were to aggressive
	//nm+=brec.clipR;

	if(readlist[n]->juncs.Count()) {
		bool mismatch=false;
		if(readlist[n]->longread) mismatch=true;
		else if(nm/readlist[n]->len>mismatchfrac) mismatch=true;
		else if(nm && readlist[n]->juncs.Count()) {
			if(brec.clipL && readlist[n]->segs[0].len()<junctionsupport+brec.clipL) mismatch=true; // penalize mismatch that's too close to ss
			else if(brec.clipR && readlist[n]->segs.Last().len()<junctionsupport+brec.clipR) mismatch=true;
			else if(mismatch_anchor(readlist[n],brec.tag_str("MD"),currentstart,brec.get_b())) mismatch=true; // this line was not initially present in vs1 or vs3 but I noticed it doesn't do any difference in real data, so far it only helped with the SR in simulation -> I might want to take it out
		}

		for(int i=0;i<readlist[n]->juncs.Count();i++) { // if read is PacBio I might want to increase the mismatch fraction, although the nm only gets used for longintrons
			if(mismatch || nh>2) readlist[n]->juncs[i]->nm+=rdcount;
			if(readlist[n]->segs[i].len()>longintronanchor && readlist[n]->segs[i+1].len()>longintronanchor)
				readlist[n]->juncs[i]->mm+=rdcount;
			//if(nh>2) readlist[n]->juncs[i]->mm+=rdcount; // vs3; vs2 only has nh>1
			readlist[n]->juncs[i]->nreads+=rdcount;
		}
	}


	// now set up the pairing
	if (brec.refId()==brec.mate_refId()) {  //only consider mate pairing data if mates are on the same chromosome/contig and are properly paired
	//if (brec.refId()==brec.mate_refId() && brec.isProperlyPaired()) {  //only consider mate pairing data if mates are on the same chromosome/contig and are properly paired
	//if (brec.isProperlyPaired()) {  //only consider mate pairing data if mates  are properly paired
		int pairstart=brec.mate_start();
		if (currentstart<=pairstart) { // if pairstart is in a previous bundle I don't care about it
			//GStr readname();
			//GStr id(brec.name(), 16); // init id with readname
			_id.assign(brec.name()); //assign can be forced to prevent shrinking of the string
			if(pairstart<=readstart) { // if I've seen the pair already <- I might not have seen it yet because the pair starts at the same place
				_id+='-';_id+=pairstart;
				_id+=".=";_id+=hi; // (!) this suffix actually speeds up the hash by improving distribution!
				const int* np=hashread[_id.chars()];
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
					hashread.Remove(_id.chars());
				}
			}
			else { // I might still see the pair in the future
				_id+='-';_id+=readstart; // this is the correct way
				_id+=".=";_id+=hi;
				hashread.Add(_id.chars(), n);
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

void add_group_to_bundle(CGroup *group, CBundle *bundle, GPVec<CBundlenode>& bnode, uint localdist){


	CBundlenode *currlastnode=bnode[bundle->lastnodeid];
	int bid=bnode.Count();

	if(group->start > currlastnode->end + localdist) { // group after last bnode
		CBundlenode *currbnode=new CBundlenode(group->start,group->end,group->cov_sum,bid);
		currlastnode->nextnode=currbnode;
		bnode.Add(currbnode);
		bundle->lastnodeid=bid;
		bundle->len+=group->end-group->start+1;
		bundle->cov+=group->cov_sum;
		bundle->multi+=group->multi;
	}
	else { // group overlaps bnode within bundledist
		if(currlastnode->end < group->end) {
		    bundle->len+= group->end - currlastnode->end;
		    currlastnode->end= group->end;
		}
		bundle->cov+=group->cov_sum;
		bundle->multi+=group->multi;
		currlastnode->cov+=group->cov_sum;
	}
}

int create_bundle(GPVec<CBundle>& bundle,CGroup *group,GPVec<CBundlenode>& bnode) {

	int bid=bnode.Count();
	int bno=bundle.Count();
	CBundlenode *startbnode=new CBundlenode(group->start,group->end,group->cov_sum,bid);

	CBundle *newbundle=new CBundle(group->end-group->start+1,group->cov_sum,group->multi,bid,bid);
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

int longtrCmp(const pointer p1, const pointer p2) {
	CTransfrag *a=(CTransfrag*)p1;
	CTransfrag *b=(CTransfrag*)p2;
	if(!a->guide && b->guide) return 1;
	if(a->guide && !b->guide) return -1;
	if(a->abundance<b->abundance) return 1;  // most abundant transcript comes first (I want to keep the one that clearly dominates)
	if(a->abundance>b->abundance) return -1;
	if(a->nodes.Count()<b->nodes.Count()) return 1;  // transfrag with more nodes comes first
	if(a->nodes.Count()>b->nodes.Count()) return -1;
	if(a->pattern.count()<b->pattern.count()) return 1; // more complete transfrag comes first
	if(a->pattern.count()>b->pattern.count()) return -1;
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

	//fprintf(stderr,"merge group=%d into group=%d\n",group2->grid,group1->grid);

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
	}
	else if(group1->color>group2->color) {
		eqcol[group1->color]=group2->color;
		group1->color=group2->color;
	}

	group1->cov_sum+=group2->cov_sum;
	group1->next_gr=group2->next_gr; // this is possible because group1->next_gr=group2

	merge[group2->grid]=group1->grid;

	group1->multi+=group2->multi;

	// delete group2
	group.freeItem(group2->grid);
}


int merge_read_to_group(int n,int np, int p, float readcov, int sno,int readcol,GList<CReadAln>& readlist,int color,GPVec<CGroup>& group,
		CGroup **allcurrgroup,CGroup **startgroup,GVec<int> *readgroup,GVec<int>& eqcol,GVec<int>& merge,int *usedcol) {

	//fprintf(stderr,"merge readcol=%d for read=%d:%d-%d with paird=%d and sno=%d\n",readcol,n,readlist[n]->start,readlist[n]->end,np,sno);
	// check if read was processed before as a fragment
	uint localdist=0;
	if(!longreads && !mixedMode) localdist=bundledist+longintronanchor;
	//if(!longreads) localdist=bundledist+longintronanchor;
	if(np>-1 && np<n && readlist[np]->end<readlist[n]->start && readlist[np]->end+localdist>readlist[n]->start &&
			(readlist[np]->segs.Last().len()>=junctionsupport || (readlist[np]->juncs.Count() && readlist[np]->juncs.Last()->strand))
			&& (readlist[n]->segs[0].len()>=junctionsupport || (readlist[n]->juncs.Count() && readlist[n]->juncs[0]->strand))) { // valid last/first exon and close
			return color;
	}

	bool first=true;

	CGroup *currgroup=allcurrgroup[sno]; // currgroup will be the group that contains readlist[n]->start; initially is the group where the previous read on this strand was added

	if(currgroup != NULL) { // this type of group - negative, unknown, or positive - was created before

		//set currgroup first
		CGroup *lastgroup=NULL;
		while(currgroup!=NULL && readlist[n]->start > currgroup->end) { // while read start after the current group's end advance group -> I might have more groups leaving from current group due to splicing
		    lastgroup=currgroup;
		    currgroup=currgroup->next_gr;
		}

		if(currgroup==NULL || readlist[n]->segs[0].end < currgroup->start) // currgroup is null only if we reached end of currgroup list
																		   // because currgroup is not NULL initially and it starts BEFORE read
			currgroup=lastgroup; // making sure read comes after currgroup

		// now process each group of coordinates individually
		CGroup *thisgroup=currgroup;
		int ncoord=readlist[n]->segs.Count(); // number of "exons" in read
		int lastpushedgroup=-1;
		// bool added=false;

		int i=0;
		while(i<ncoord) {

			bool keep=true;
			// determine if this exon is good enough to be kept (it is big enough and has good junctions)
			if(readlist[n]->segs[i].len()<junctionsupport || (readlist[n]->longread && readlist[n]->segs[i].len()<CHI_THR && readlist[n]->segs[i].len()<DROP*readlist[n]->len)) { // exon is too small to keep

				if(i<readlist[n]->juncs.Count()) {
					if(!readlist[n]->juncs[i]->strand) { // check first exon
						// see if first exon is big enough
						if(!i || !readlist[n]->juncs[i-1]->strand) {
							keep=false;
							if(!i && lastgroup) currgroup=lastgroup;
						}
					}
				}
				else if(readlist[n]->juncs.Count()){ // this is last exon
					if(!readlist[n]->juncs[i-1]->strand) keep=false;
				}

			}


			if(keep) {
				// skip groups that are left behind
				while(thisgroup!=NULL && readlist[n]->segs[i].start > thisgroup->end) { // find the group where "exon" fits
					lastgroup=thisgroup;
					thisgroup=thisgroup->next_gr;
				}

				if(thisgroup && readlist[n]->segs[i].end + localdist >= thisgroup->start) { // read overlaps group or it's close enough -> it should get group color

					// I need to split pairs here if color didn't reach this group: it means there is a gap between these groups and no reads joining them
					if(!i) {

						int gr=thisgroup->grid;
						while(merge[gr]!=gr) {
							gr=merge[gr];
						}
						thisgroup->grid=gr;

						for(int g=0;g<readgroup[n].Count();g++) {
							gr=readgroup[n][g];
							while(merge[gr]!=gr) {
								gr=merge[gr];
							}
							if(gr==thisgroup->grid) {
								first=false; // I did see this groups before
								break;
							}
						}

						if(np>-1 && readlist[np]->nh && np<n) { // I only consider first exon here because the rest of the groups need to get the same color if the read is still paired

							int thiscol=thisgroup->color;
							while(eqcol[thiscol]!=thiscol) { // get smallest color
								thiscol=eqcol[thiscol];
							}
							thisgroup->color=thiscol;

							if(thiscol!=readcol) { // pair color didn't reach this group

							  //fprintf(stderr,"Split pairs: %d-%d and %d-%d on strand %d thiscol=%d readcol=%d\n",readlist[np]->start,readlist[np]->end,readlist[n]->start,readlist[n]->end,readlist[n]->strand,thiscol,readcol);
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

					}

					if(readlist[n]->nh>1) {
						thisgroup->multi+=readcov*readlist[n]->segs[i].len(); // where do I add readcov to group coverage -> see below at the end of '}'
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

					if(readcol!=thisgroup->color)  { // read color is different from group color
						if(i && !readlist[n]->juncs[i-1]->strand) { // bad junction before this exon
							// before I equal readcol to thisgroup->color, I need to make sure that the read is not interrupted
							// readcol for next exons should still be the color of thisgroup
							readcol=thisgroup->color;
						}
						else  {
							if(readcol<thisgroup->color) { // set group color to current read color
								eqcol[thisgroup->color]=readcol;
								thisgroup->color=readcol;
							}
							else { // read color is bigger than group

								eqcol[readcol]=thisgroup->color;
								readcol=thisgroup->color;
							}
						}
					}

					if(thisgroup->grid != lastpushedgroup) {

						//fprintf(stderr,"Assign group %d:%d-%d with color=%d to read %d\n",thisgroup->grid,thisgroup->start,thisgroup->end,thisgroup->color,n);
						if(first) readgroup[n].Add(thisgroup->grid);   // readgroup for read n gets the id of group
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
		    			else if(i && !readlist[n]->juncs[i-1]->strand) { // if exons are interrupted use new color
		    				// I have to use new color for read and I have no idea if usedcol was already used or not => I have to ignore what was already assigned to usedcol
		    				usedcol[sno]=color;
		    				readcol=color;
		    				eqcol.Add(color);
		    				color++;
		    			}

		    			int ngroup=group.Count();
		    			float multi=0;
		    			if(readlist[n]->nh>1) multi=readcov*readlist[n]->segs[i].len();

		    			CGroup *newgroup=new CGroup(readlist[n]->segs[i].start,readlist[n]->segs[i].end,readcol,ngroup,readlist[n]->segs[i].len()*readcov,multi);
		    			group.Add(newgroup);
		    			merge.Add(ngroup);
		    			lastgroup->next_gr=newgroup; // can lastgroup be null here -> no from the way I got here
		    			newgroup->next_gr=thisgroup;

		    			//fprintf(stderr,"Assign/create group %d:%d-%d with color=%d to read %d\n",ngroup,readlist[n]->segs[i].start,readlist[n]->segs[i].end,readcol,n);

		    			readgroup[n].Add(ngroup); // since group was created it clearly is new
		    			lastpushedgroup=ngroup;
		    			thisgroup=newgroup;
				}

				/* version that uses paired reads in a fragment to do bundling*/
				if(i==ncoord-1) { // last segment in read
					if(np!=-1 && np>n && readlist[n]->end<readlist[np]->start && readlist[n]->end+localdist>readlist[np]->start &&
							(readlist[np]->segs[0].len()>=junctionsupport || (readlist[np]->juncs.Count() && readlist[np]->juncs[0]->strand))) { // valid last exon and close
						n=np; // next I will process read pair as it was part of the same fragment
						np=-1;
						ncoord=readlist[n]->segs.Count();
						first=true;
						if(readlist[n]->start>thisgroup->end) { // I need to extend group
							CGroup *nextgroup=thisgroup->next_gr;
							while(nextgroup!=NULL && readlist[n]->start >= nextgroup->start) {
								merge_fwd_groups(group,thisgroup,nextgroup,merge,eqcol);
								nextgroup=thisgroup->next_gr;
							}
							if(readlist[n]->start > thisgroup->end) {
								thisgroup->end=readlist[n]->start;
							}
						}
						lastpushedgroup=-1;
						i=-1;
					}
				}

			} // if(keep)
			i++;
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

			// even if I split pairs I do not need to update single_count because this readcov is taken care of
		}

		int ncoord=readlist[n]->segs.Count();
		CGroup *lastgroup=NULL;
		float multi=0;
		if(readlist[n]->nh>1) multi=readcov; // this will be adjusted by the length of each segment
		int i=0;
		while(i<ncoord) {

			bool keep=true;
			// determine if this exon is good enough to be kept (it is big enough and has good junctions)
			if(readlist[n]->segs[i].len()<junctionsupport) { // exon is too small to keep --> this exon will not be added to group if it doesn't pass some criteria

				if(i<readlist[n]->juncs.Count()) {
					if(!readlist[n]->juncs[i]->strand) { // check first exon
						// see if first exon is big enough
						if(!i || !readlist[n]->juncs[i-1]->strand) keep=false;
					}
				}
				else if(readlist[n]->juncs.Count()){ // this is last exon
					if(!readlist[n]->juncs[i-1]->strand) keep=false;
				}

			}


			if(!i || keep) {

				// assign new color if read is interrupted
				if(i && !readlist[n]->juncs[i-1]->strand) { // if exons are interrupted use new color
					usedcol[sno]=color;
					readcol=color;
					eqcol.Add(color);
					color++;
				}

				int ngroup=group.Count();
				CGroup *newgroup=new CGroup(readlist[n]->segs[i].start,readlist[n]->segs[i].end,readcol,ngroup,readlist[n]->segs[i].len()*readcov,readlist[n]->segs[i].len()*multi);
				group.Add(newgroup);
				merge.Add(ngroup);
				if(lastgroup!=NULL) {
					lastgroup->next_gr=newgroup;
				}
				else {
					currgroup=newgroup;
				}
				lastgroup=newgroup;
				if(keep) {
					//fprintf(stderr,"Assign/create2 group %d:%d-%d with color=%d to read %d\n",ngroup,readlist[n]->segs[i].start,readlist[n]->segs[i].end,readcol,n);

					readgroup[n].Add(ngroup); // add group to read so that I can remember what groups belong to a read

					/* version that uses paired reads in a fragment to do bundling*/
					if(i==ncoord-1) { // last segment in read
						if(np!=-1 && np>n && readlist[n]->end<readlist[np]->start && readlist[n]->end+localdist>readlist[np]->start &&
								(readlist[np]->segs[0].len()>=junctionsupport || (readlist[np]->juncs.Count() && readlist[np]->juncs[0]->strand))) { // valid last exon and close
							n=np; // next I will process read pair as it was part of the same fragment
							np=-1;
							ncoord=readlist[n]->segs.Count();
							if(readlist[n]->segs[0].end>lastgroup->end) { // I need to extend group
								//fprintf(stderr,"Group %d:%d-%d will get end at %d n=%d ncoord=%d\n",lastgroup->grid,lastgroup->start,lastgroup->end,readlist[n]->segs[0].end,n,ncoord);
								lastgroup->end=readlist[n]->segs[0].end;
								if(first) readgroup[n].Add(ngroup); // add group to read so that I can remember what groups belong to a read
								lastgroup->cov_sum+=readlist[n]->segs[0].len()*readcov; // coverage is different than number of reads
							}
							i=0;
						}
					}
				}
			}
			i++;
		}
	}

	allcurrgroup[sno]=currgroup;

	if(startgroup[sno]==NULL) startgroup[sno]=currgroup;

	return color;
}

int add_read_to_group(int n,GList<CReadAln>& readlist,int color,GPVec<CGroup>& group,CGroup **allcurrgroup,
		CGroup **startgroup,GVec<int> *readgroup,GVec<int>& eqcol,GVec<int>& merge) {

	int usedcol[3]={-1,-1,-1}; // remembers last usedcol for read so that I don't use different colors for the same read

	float single_count=readlist[n]->read_count; // need to compute read single count

	if(!mergeMode) for(int p=0;p<readlist[n]->pair_idx.Count();p++) { // for all pairs of read

		// for each pair I have to pretend is independent of the previous one since I grouped together several reads which might be independent
		int sno=readlist[n]->strand+1; // 0: negative strand; 1: zero strand; 2: positive strand (starting from -1,0,1) // I have to reset sno every time because it might change due to snop

		// check if I've seen read's pair and if yes get its readcol; at the least get read's pair strand if available
		int np=readlist[n]->pair_idx[p]; // pair read number

		if(np>-1 && readlist[np]->nh) { // read pair still exists and it wasn't deleted -> false keep elminates all these cases

			// see if I have the correct read strand
			char snop=readlist[np]->strand+1;  // snop is the strand of pair read
			if(sno!=snop) { // different strands for read and pair

				if(sno==1) { // read n is on zero (neutral) strand, but pair has strand
					// readlist[n]->strand=readlist[np]->strand; I can not update the read strand anymore -> REMEMBER this later
					sno=snop;  // assign strand of pair to read
				}
				else if(snop!=1) { // conflicting strands -> un-pair reads in the hope that one is right
					readlist[n]->pair_idx[p]=-1;
					for(int j=0;j<readlist[np]->pair_idx.Count();j++) // also unpair read np to n
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
					// first group of pair read is: readgroup[np][0]

					int grouppair=readgroup[np][0];
					while( merge[grouppair]!=grouppair) {
						grouppair=merge[grouppair];
					}
					readgroup[np][0]=grouppair;

					readcol=group[readgroup[np][0]]->color;    // readcol gets assigned the color of the pair's first group --> this might not
															   // work if I split  a read into multiple groups
					while(eqcol[readcol]!=readcol) { // get smallest color
						readcol=eqcol[readcol];
					}
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

				float readcov=0;
				if(!readlist[n]->unitig)
					readcov=readlist[n]->pair_count[p];
				color=merge_read_to_group(n,np,p,readcov,sno,readcol,readlist,color,group,
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
		float readcov=0;
		if(!readlist[n]->unitig)
			readcov=single_count;
		color=merge_read_to_group(n,-1,-1,readcov,readlist[n]->strand+1,readcol,readlist,color,group,
								allcurrgroup,startgroup,readgroup,eqcol,merge,usedcol);

	}

	return color;
}

CGraphnode *create_graphnode(int s, int g, uint start,uint end,int nodeno,CBundlenode *bundlenode,
		GVec<CGraphinfo> **bundle2graph,GPVec<CGraphnode> **no2gnode) {

	/*
	{ // DEBUG ONLY
		fprintf(stderr,"create_graphnode[%d][%d]:%d-%d nodeno=%d\n",s,g,start,end,nodeno);
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

float get_cov(int s,uint start,uint end,GVec<float>* bpcov) {
	int m=int((end+1)/BSIZE);
	int k=int(start/BSIZE);
	float cov=0;
	for(int i=k+1;i<=m;i++)
		cov+=bpcov[s][i*BSIZE-1];
	cov+=bpcov[s][end+1]-bpcov[s][start];
	if(cov<0) cov=0;
	return(cov);

}

float get_cov_sign(int s,uint start,uint end,GVec<float>* bpcov) {
	int m=int((end+1)/BSIZE);
	int k=int(start/BSIZE);
	int o=2-s;
	float cov=0;
	for(int i=k+1;i<=m;i++)
		cov+=bpcov[1][i*BSIZE-1]-bpcov[o][i*BSIZE-1];
	cov+=bpcov[1][end+1]-bpcov[1][start]-bpcov[o][end+1]+bpcov[o][start];
	if(cov<0) cov=0;
	return(cov);

}

// cummulative bpcov; simple find trims that only checks for significant sudden drops in coverage
//void find_trims(int refstart,int sno,uint start,uint end,GVec<float>* bpcov,uint& sourcestart,float& sourcecovleft, float& sourcecovright,uint& sinkend, float& sinkcovleft,float& sinkcovright){
//void find_trims(int refstart,int sno,uint start,uint end,GVec<float>* bpcov,uint& sourcestart,float& sourceabundance,uint& sinkend,float& sinkabundance) {
void find_trims(int refstart,uint start,uint end,GVec<float>* bpcov,uint& sourcestart,float& sourceabundance,uint& sinkend,float& sinkabundance) {


	int len=end-start+1; // length of region where I look for trims
	if(len<CHI_THR) return; // very short exon -> do not check
	//if(len<CHI_WIN) return; // very short exon -> do not check

	float sourcedrop=1;
	float sinkdrop=1;

	float localdrop=ERROR_PERC;
	// bigger trimming window
	if(len<2*(CHI_WIN+CHI_THR)+1) {
	  if(len<CHI_WIN) localdrop=ERROR_PERC/(10*DROP);
	  	  for(uint i=start+longintronanchor;i<end-longintronanchor;i++) {
			float covleft=get_cov(1,start-refstart,i-1-refstart,bpcov)/(i-start);
			float covright=get_cov(1,i-refstart,end-refstart,bpcov)/(end-i+1);
			if(covleft<covright) {
				float thisdrop=covleft/covright;
				if(thisdrop<localdrop && thisdrop<sourcedrop) {
					sourcedrop=thisdrop;
					sourcestart=i;
					sourceabundance=(covright-covleft)/DROP;
				}
			}
			else if(covright<covleft) { // possible sink trimming
				float thisdrop=covright/covleft;
				if(thisdrop<localdrop && thisdrop<sinkdrop) {
					sinkdrop=thisdrop;
					sinkend=i-1;
					sinkabundance=(covleft-covright)/DROP;
				}
			}
		}
		return;
	}

	int winlen=CHI_WIN+CHI_THR;

	localdrop=ERROR_PERC/DROP;

	len+=start-winlen;

	for(int i=start+winlen-1;i<len;i++) {
		float covleft=get_cov(1,i-winlen+1-refstart,i-refstart,bpcov); // I add 1 bp coverage to all so I can avoid 0 coverages
		float covright=get_cov(1,i+1-refstart,i+winlen-refstart,bpcov);
		if(covleft<covright) { // possible source trimming
			float thisdrop=covleft/covright;
			if(thisdrop<localdrop && thisdrop<sourcedrop) {
				sourcedrop=thisdrop;
				sourcestart=i+1;
				//maxsourceabundance=covright-covleft;
				//sourcecovright=covright/(DROP*winlen);
				//sourcecovleft=covleft/(DROP*winlen);
				//sourcecovright=(covright-covleft)/winlen;
				//sourcecovright=(covright-covleft)/(DROP*winlen);
				sourceabundance=(covright-covleft)/(DROP*winlen);
			}
		}
		else if(covright<covleft) { // possible sink trimming
			float thisdrop=covright/covleft;
			if(thisdrop<localdrop && thisdrop<sinkdrop) {
				sinkdrop=thisdrop;
				sinkend=i;
				//maxsinkabundance=covleft-covright;
				//sinkcovright=covright/(DROP*winlen);
				//sinkcovleft=covleft/(DROP*winlen);
				//sinkcovleft=(covleft-covright)/winlen;
				//sinkcovleft=(covleft-covright)/(DROP*winlen);
				sinkabundance=(covleft-covright)/(DROP*winlen);
			}
		}
	}

	// make sure the drops are consistent across the whole exon
	if(sourcestart) {
		if(sourcestart<sinkend) {
			float covleft=get_cov(1,start-refstart,sourcestart-1-refstart,bpcov)/(sourcestart-start);
			float covright=localdrop*get_cov(1,sourcestart-refstart,sinkend-refstart,bpcov)/(sinkend-sourcestart+1);
			if(covleft>covright) sourcestart=0;
			covleft=get_cov(1,sinkend+1-refstart,end-refstart,bpcov)/(end-sinkend);
			if(covleft>covright) sinkend=0;
		}
		else if(sinkend) {
			float covleft=get_cov(1,start-refstart,sinkend-refstart,bpcov)/(sinkend-start+1);
			float covright=get_cov(1,sinkend+1-refstart,sourcestart-1-refstart,bpcov)/(sourcestart-sinkend-1);
			if(covright>covleft*localdrop) sinkend=0;
			covleft=get_cov(1,sourcestart-refstart,end-refstart,bpcov)/(end-sourcestart+1);
			if(covright>covleft*localdrop) sourcestart=0;
		}
		else {
			float covleft=get_cov(1,start-refstart,sourcestart-1-refstart,bpcov)/(sourcestart-start);
			float covright=localdrop*get_cov(1,sourcestart-refstart,end-refstart,bpcov)/(end-sourcestart+1);
			if(covleft>covright) sourcestart=0;
		}
	}
	else if(sinkend) {
		float covright=localdrop*get_cov(1,start-refstart,sinkend-refstart,bpcov)/(sinkend-start+1);
		float covleft=get_cov(1,sinkend+1-refstart,end-refstart,bpcov)/(end-sinkend);
		if(covleft>covright) sinkend=0;
	}
}

void find_trims_wsign(int refstart,int sno,uint start,uint end,GVec<float>* bpcov,uint& sourcestart,float& sourceabundance,uint& sinkend,float& sinkabundance) {


  //fprintf(stderr,"find trims in region:%d-%d:%d\n",start,end,sno);

	int len=end-start+1; // length of region where I look for trims
	if(len<CHI_THR) return; // very short exon -> do not check
	//if(len<CHI_WIN) return; // very short exon -> do not check

	float sourcedrop=1;
	float sinkdrop=1;

	//int os=2-sno; // other strand

	float localdrop=ERROR_PERC;
	// bigger trimming window
	if(len<2*(CHI_WIN+CHI_THR)+1) {
	  if(len<CHI_WIN) localdrop=ERROR_PERC/(10*DROP);
	  for(uint i=start+longintronanchor;i<end-longintronanchor;i++) {
		  //float covleft=(get_cov(1,start-refstart,i-1-refstart,bpcov)-get_cov(os,start-refstart,i-1-refstart,bpcov))/(i-start);
		  //float covright=(get_cov(1,i-refstart,end-refstart,bpcov)-get_cov(os,i-refstart,end-refstart,bpcov))/(end-i+1);
		  float covleft=get_cov_sign(sno,start-refstart,i-1-refstart,bpcov)/(i-start);
		  float covright=get_cov_sign(sno,i-refstart,end-refstart,bpcov)/(end-i+1);
		  if(covleft<covright && covright>0) {
			  float thisdrop=covleft/covright;
			  if(thisdrop<localdrop && thisdrop<sourcedrop) {
				  sourcedrop=thisdrop;
				  sourcestart=i;
				  sourceabundance=(covright-covleft)/DROP;
			  }
		  }
		  else if(covright<covleft && covleft>0) { // possible sink trimming
			  float thisdrop=covright/covleft;
			  if(thisdrop<localdrop && thisdrop<sinkdrop) {
				  sinkdrop=thisdrop;
				  sinkend=i-1;
				  sinkabundance=(covleft-covright)/DROP;
			  }
		  }
	  }
	  return;
	}

	int winlen=CHI_WIN+CHI_THR;

	/*localdrop=ERROR_PERC*DROP; // sharp drop
	// look for short drops first
	for(uint i=start+CHI_THR-1;i<start+winlen-1;i++) {
		float covleft=get_cov_sign(sno,i-CHI_THR+1-refstart,i-refstart,bpcov); // I add 1 bp coverage to all so I can avoid 0 coverages
		float covright=get_cov_sign(sno,i+1-refstart,i+CHI_THR-refstart,bpcov);
		if(covleft<covright && covright>0) {
			float thisdrop=covleft/covright;
			if(thisdrop<localdrop && thisdrop<sourcedrop) {
				sourcedrop=thisdrop;
				sourcestart=i;
				sourceabundance=(covright-covleft)/DROP;
			}
		}
		else if(covright<covleft && covleft>0) { // possible sink trimming
			float thisdrop=covright/covleft;
			if(thisdrop<localdrop && thisdrop<sinkdrop) {
				sinkdrop=thisdrop;
				sinkend=i-1;
				sinkabundance=(covleft-covright)/DROP;
			}
		}
	}*/


	localdrop=ERROR_PERC/DROP;

	len+=start-winlen;

	for(int i=start+winlen-1;i<len;i++) {
		//float covleft=get_cov(1,i-winlen+1-refstart,i-refstart,bpcov)-get_cov(os,i-winlen+1-refstart,i-refstart,bpcov); // I add 1 bp coverage to all so I can avoid 0 coverages
		//float covright=get_cov(1,i+1-refstart,i+winlen-refstart,bpcov)-get_cov(os,i+1-refstart,i+winlen-refstart,bpcov);
		float covleft=get_cov_sign(sno,i-winlen+1-refstart,i-refstart,bpcov); // I add 1 bp coverage to all so I can avoid 0 coverages
		float covright=get_cov_sign(sno,i+1-refstart,i+winlen-refstart,bpcov);
		if(covleft<covright && covright>0) { // possible source trimming
			float thisdrop=covleft/covright;
			if(thisdrop<localdrop && thisdrop<sourcedrop) {
				sourcedrop=thisdrop;
				sourcestart=i+1;
				//maxsourceabundance=covright-covleft;
				//sourcecovright=covright/(DROP*winlen);
				//sourcecovleft=covleft/(DROP*winlen);
				//sourcecovright=(covright-covleft)/winlen;
				//sourcecovright=(covright-covleft)/(DROP*winlen);
				sourceabundance=(covright-covleft)/(DROP*winlen);
			}
		}
		else if(covright<covleft && covleft>0) { // possible sink trimming
			float thisdrop=covright/covleft;
			if(thisdrop<localdrop && thisdrop<sinkdrop) {
				sinkdrop=thisdrop;
				sinkend=i;
				//maxsinkabundance=covleft-covright;
				//sinkcovright=covright/(DROP*winlen);
				//sinkcovleft=covleft/(DROP*winlen);
				//sinkcovleft=(covleft-covright)/winlen;
				//sinkcovleft=(covleft-covright)/(DROP*winlen);
				sinkabundance=(covleft-covright)/(DROP*winlen);
			}
		}
	}

	// make sure the drops are consistent across the whole exon --> this might be important
	if(sourcestart) {
		if(sourcestart<sinkend) {
			//float covleft=(get_cov(1,start-refstart,sourcestart-1-refstart,bpcov)-get_cov(os,start-refstart,sourcestart-1-refstart,bpcov))/(sourcestart-start);
			//float covright=localdrop*(get_cov(1,sourcestart-refstart,sinkend-refstart,bpcov)-get_cov(os,sourcestart-refstart,sinkend-refstart,bpcov))/(sinkend-sourcestart+1);
			float covleft=get_cov_sign(sno,start-refstart,sourcestart-1-refstart,bpcov)/(sourcestart-start);
			float covright=localdrop*get_cov_sign(sno,sourcestart-refstart,sinkend-refstart,bpcov)/(sinkend-sourcestart+1);
			if(covleft>covright) {
				sourcestart=0;
			}
			//covleft=(get_cov(1,sinkend+1-refstart,end-refstart,bpcov)-get_cov(os,sinkend+1-refstart,end-refstart,bpcov))/(end-sinkend);
			covleft=get_cov_sign(sno,sinkend+1-refstart,end-refstart,bpcov)/(end-sinkend);
			if(covleft>covright) {
				sinkend=0;
			}
		}
		else if(sinkend) {
			//float covleft=(get_cov(1,start-refstart,sinkend-refstart,bpcov)-get_cov(os,start-refstart,sinkend-refstart,bpcov))/(sinkend-start+1);
			//float covright=(get_cov(1,sinkend+1-refstart,sourcestart-1-refstart,bpcov)-get_cov(os,sinkend+1-refstart,sourcestart-1-refstart,bpcov))/(sourcestart-sinkend-1);
			float covleft=get_cov_sign(sno,start-refstart,sinkend-refstart,bpcov)/(sinkend-start+1);
			float covright=get_cov_sign(sno,sinkend+1-refstart,sourcestart-1-refstart,bpcov)/(sourcestart-sinkend-1);
			if(covright>covleft*localdrop) {
				sinkend=0;
			}
			//covleft=(get_cov(1,sourcestart-refstart,end-refstart,bpcov)-get_cov(os,sourcestart-refstart,end-refstart,bpcov))/(end-sourcestart+1);
			covleft=get_cov_sign(sno,sourcestart-refstart,end-refstart,bpcov)/(end-sourcestart+1);
			if(covright>covleft*localdrop) {
				sourcestart=0;
			}
		}
		else {
			//float covleft=(get_cov(1,start-refstart,sourcestart-1-refstart,bpcov)-get_cov(os,start-refstart,sourcestart-1-refstart,bpcov))/(sourcestart-start);
			//float covright=localdrop*(get_cov(1,sourcestart-refstart,end-refstart,bpcov)-get_cov(os,sourcestart-refstart,end-refstart,bpcov))/(end-sourcestart+1);
			float covleft=get_cov_sign(sno,start-refstart,sourcestart-1-refstart,bpcov)/(sourcestart-start);
			float covright=localdrop*get_cov_sign(sno,sourcestart-refstart,end-refstart,bpcov)/(end-sourcestart+1);
			if(covleft>covright) {
				sourcestart=0;
			}
		}
	}
	else if(sinkend) {
		//float covright=localdrop*(get_cov(1,start-refstart,sinkend-refstart,bpcov)-get_cov(os,start-refstart,sinkend-refstart,bpcov))/(sinkend-start+1);
		//float covleft=(get_cov(1,sinkend+1-refstart,end-refstart,bpcov)-get_cov(os,sinkend+1-refstart,end-refstart,bpcov))/(end-sinkend);
		float covright=localdrop*get_cov_sign(sno,start-refstart,sinkend-refstart,bpcov)/(sinkend-start+1);
		float covleft=get_cov_sign(sno,sinkend+1-refstart,end-refstart,bpcov)/(end-sinkend);
		if(covleft>covright) {
			sinkend=0;
		}
	}
}

void find_all_trims(int refstart,int sno,uint start,uint end,GVec<float>* bpcov,GVec<CTrimPoint> &trimpoint) {

// this procedure returns all trims it finds and stores them in trimpoint

	//fprintf(stderr,"find all trims in region:%d-%d:%d\n",start,end,sno);

	int len=end-start+1; // length of region where I look for trims
	if(len<CHI_THR) return; // very short exon -> do not check

	float localdrop=ERROR_PERC; // very sharp drop
	if(len<2*(CHI_WIN+CHI_THR)+1) {
		if(len<CHI_WIN) { // check very tightly, and only accept one trim
			localdrop=ERROR_PERC/(10*DROP); // very sharp drop
		}
		float lastdrop=localdrop;
		for(uint i=start+longintronanchor;i<end-longintronanchor;i++) {
			float covleft=get_cov_sign(sno,start-refstart,i-1-refstart,bpcov)/(i-start);
			float covright=get_cov_sign(sno,i-refstart,end-refstart,bpcov)/(end-i+1);
			if(covleft<covright) { // possible source trimming
				//float thisdrop=(covleft+1)/(covright+1); // make sure I add one read to mitigate gaps in coverage
				float thisdrop=covleft/covright; // make sure I add one read to mitigate gaps in coverage
				//fprintf(stderr,"found source drop=%f covleft=%f covright=%f at i=%d\n",thisdrop,covleft,covright,i+1);
				if(thisdrop<localdrop) {
					if(!trimpoint.Count() || (!trimpoint.Last().start && i+1-(int)trimpoint.Last().pos>CHI_THR)) { // add this point to trimpoints
						CTrimPoint t(i+1,(covright-covleft)/DROP,true);
						trimpoint.Add(t);
						lastdrop=thisdrop;
					}
					else if(thisdrop<lastdrop){ // smaller drop than before --> replace drop
						lastdrop=thisdrop;
						trimpoint.Last().pos=i+1;
						trimpoint.Last().abundance=(covright-covleft)/DROP;
						trimpoint.Last().start=true;
					}
				}
			}
			else if(covright!=covleft) { // possible sink trimming: covright<covleft here
				float thisdrop=(covright+1)/(covleft+1);
				//fprintf(stderr,"found sink drop=%f covleft=%f covright=%f at i=%d\n",thisdrop,covleft,covright,i);
				if(thisdrop<localdrop) {
					if(!trimpoint.Count() || (trimpoint.Last().start && i-(int)trimpoint.Last().pos>CHI_THR)) { // add this point to trimpoints
						CTrimPoint t(i,(covleft-covright)/DROP,false);
						trimpoint.Add(t);
						lastdrop=thisdrop;
					}
					else if(thisdrop<lastdrop){ // smaller drop than before --> replace drop
						lastdrop=thisdrop;
						trimpoint.Last().pos=i;
						trimpoint.Last().abundance=(covleft-covright)/DROP;
						trimpoint.Last().start=false;
					}
				}
			}
		}
		return;
	}

	// now len >= CHI_WIN
	int winlen=CHI_WIN+CHI_THR;
	float lastdrop=localdrop;

	for(uint i=start+CHI_THR-1;i<start+winlen-1;i++) {
		float covleft=get_cov_sign(sno,start-refstart,i-refstart,bpcov)/(i-start);
		float covright=get_cov_sign(sno,i+1-refstart,i+winlen-refstart,bpcov)/winlen;
		if(covleft<covright) { // possible source trimming
			float thisdrop=covleft/covright;
			//fprintf(stderr,"found source drop=%f covleft=%f covright=%f at i=%d lastdrop=%f\n",thisdrop,covleft,covright,i+1,lastdrop);
			if(thisdrop<localdrop) { // this drop passes threshold
				if(!trimpoint.Count() || (!trimpoint.Last().start && i+1-(int)trimpoint.Last().pos>CHI_THR)) { // add this point to trimpoints
					CTrimPoint t(i+1,(covright-covleft)/DROP,true);
					trimpoint.Add(t);
					lastdrop=thisdrop;
				}
				else if(thisdrop<lastdrop){ // smaller drop than before --> replace drop; it also has to be smallest within a 50bp prevwindow
					lastdrop=thisdrop;
					trimpoint.Last().pos=i+1;
					trimpoint.Last().abundance=(covright-covleft)/DROP;
					trimpoint.Last().start=true;
				}
			}
		}
		else if(covleft!=covright) { // possible sink trimming: covright<covleft here
			float thisdrop=covright/covleft;
			//fprintf(stderr,"found sink drop=%f covleft=%f covright=%f at i=%d lastdrop=%f\n",thisdrop,covleft,covright,i,lastdrop);
			if(thisdrop<localdrop) { // this drop passes threshold
				if(!trimpoint.Count() || (trimpoint.Last().start && i-(int)trimpoint.Last().pos>CHI_THR)) { // add this point to trimpoints
					CTrimPoint t(i,(covleft-covright)/DROP,false);
					trimpoint.Add(t);
					lastdrop=thisdrop;
				}
				else if(thisdrop<lastdrop){ // smaller drop than before --> replace drop
					lastdrop=thisdrop;
					trimpoint.Last().pos=i;
					trimpoint.Last().abundance=(covleft-covright)/DROP;
					trimpoint.Last().start=false;
				}
			}
		}
	}


	localdrop=ERROR_PERC/DROP; // I might want to use different parameters here
	if(mixedMode) localdrop=DROP*DROP; // tolerate more for mixedMode
	if(!trimpoint.Count()) lastdrop=localdrop;
	len+=start-winlen;
	for(int i=start+winlen-1;i<len;i++) {
		float covleft=get_cov_sign(sno,i-winlen+1-refstart,i-refstart,bpcov); // I add 1 bp coverage to all so I can avoid 0 coverages
		float covright=get_cov_sign(sno,i+1-refstart,i+winlen-refstart,bpcov);
		if(covleft<covright) { // possible source trimming
			//float thisdrop=(covleft+winlen)/(covright+winlen); // one read is added in order to avoid gaps in coverage
			float thisdrop=covleft/covright;
			if(thisdrop<localdrop) { // this drop passes threshold
				//fprintf(stderr,"found source drop=%f covleft=%f covright=%f at i=%d window:%d-%d lastdrop=%f\n",thisdrop,covleft,covright,i+1,i-winlen+1,i+winlen,lastdrop);
				if(!trimpoint.Count() || (!trimpoint.Last().start && i+1-(int)trimpoint.Last().pos>CHI_THR)) { // add this point to trimpoints
					CTrimPoint t(i+1,(covright-covleft)/(DROP*winlen),true);
					trimpoint.Add(t);
					lastdrop=thisdrop;
				}
				else if(thisdrop<lastdrop){ // smaller drop than before --> replace drop; it also has to be smallest within a 50bp prevwindow
					lastdrop=thisdrop;
					trimpoint.Last().pos=i+1;
					trimpoint.Last().abundance=(covright-covleft)/(DROP*winlen);
					trimpoint.Last().start=true;
				}
			}
		}
		else if(covleft!=covright) { // possible sink trimming: covright<covleft here
			//float thisdrop=(covright+winlen)/(covleft+winlen); // one read is added in order to avoid gaps in coverage
			float thisdrop=covright/covleft;
			if(thisdrop<localdrop) { // this drop passes threshold
				//fprintf(stderr,"found sink drop=%f covleft=%f covright=%f at i=%d window:%d-%d lastdrop=%f\n",thisdrop,covleft,covright,i,i-winlen+1,i+winlen,lastdrop);
				if(!trimpoint.Count() || (trimpoint.Last().start && i-(int)trimpoint.Last().pos>CHI_THR)) { // add this point to trimpoints
					CTrimPoint t(i,(covleft-covright)/(DROP*winlen),false);
					trimpoint.Add(t);
					lastdrop=thisdrop;
				}
				else if(thisdrop<lastdrop){ // smaller drop than before --> replace drop
					lastdrop=thisdrop;
					trimpoint.Last().pos=i;
					trimpoint.Last().abundance=(covleft-covright)/(DROP*winlen);
					trimpoint.Last().start=false;
				}
			}
		}
	}

	localdrop=ERROR_PERC;

	for(uint i=len;i<end-CHI_THR-1;i++) {
		float covleft=get_cov_sign(sno,i-winlen+1-refstart,i-refstart,bpcov)/winlen;
		float covright=get_cov_sign(sno,i+1-refstart,end-refstart,bpcov)/(end-i);
		if(covleft<covright) { // possible source trimming
			float thisdrop=covleft/covright;
			//fprintf(stderr,"found source drop=%f covleft=%f covright=%f at i=%d lastdrop=%f\n",thisdrop,covleft,covright,i+1,lastdrop);
			if(thisdrop<localdrop) { // this drop passes threshold
				if(!trimpoint.Count() || (!trimpoint.Last().start && i+1-(int)trimpoint.Last().pos>CHI_THR)) { // add this point to trimpoints
					CTrimPoint t(i+1,(covright-covleft)/DROP,true);
					trimpoint.Add(t);
					lastdrop=thisdrop;
				}
				else if(thisdrop<lastdrop){ // smaller drop than before --> replace drop; it also has to be smallest within a 50bp prevwindow
					lastdrop=thisdrop;
					trimpoint.Last().pos=i+1;
					trimpoint.Last().abundance=(covright-covleft)/DROP;
					trimpoint.Last().start=true;
				}
			}
		}
		else if(covleft!=covright) { // possible sink trimming: covright<covleft here
			float thisdrop=covright/covleft;
			//fprintf(stderr,"found sink drop=%f covleft=%f covright=%f at i=%d lastdrop=%f\n",thisdrop,covleft,covright,i,lastdrop);
			if(thisdrop<localdrop) { // this drop passes threshold
				if(!trimpoint.Count() || (trimpoint.Last().start && i-(int)trimpoint.Last().pos>CHI_THR)) { // add this point to trimpoints
					CTrimPoint t(i,(covleft-covright)/DROP,false);
					trimpoint.Add(t);
					lastdrop=thisdrop;
				}
				else if(thisdrop<lastdrop){ // smaller drop than before --> replace drop
					lastdrop=thisdrop;
					trimpoint.Last().pos=i;
					trimpoint.Last().abundance=(covleft-covright)/DROP;
					trimpoint.Last().start=false;
				}
			}
		}
	}

	localdrop=ERROR_PERC/DROP; // I might want to use different parameters here
	// I might adjust localdrop again here
	uint laststart=start;
	for(int i=0;i<trimpoint.Count();i++) {
		uint midpos=trimpoint[i].pos;
		if(trimpoint[i].start) midpos--;
		float covleft=get_cov_sign(sno,laststart-refstart,midpos-refstart,bpcov)/(midpos-laststart+1); // I add 1 bp coverage to all so I can avoid 0 coverages
		uint endpos=end;
		if(i<trimpoint.Count()-1) endpos=trimpoint[i+1].pos;
		float covright=get_cov_sign(sno,midpos+1-refstart,endpos-refstart,bpcov)/(endpos-midpos);

		if((trimpoint[i].start && covleft<localdrop*covright)||(covright<localdrop*covleft)) { // source|sink above threshold
			laststart=trimpoint[i].pos+1;
		}
		else { // trimpoint is not above threshold --> do not change laststart, but I need to reevaluate previous trimpoints
		  //fprintf(stderr,"pos=%d failed global test btw %d-%d midpos=%d covlef=%.1f covright=%.1f\n",trimpoint[i].pos,laststart,endpos,midpos,covleft,covright);
			trimpoint[i].pos=0;
			int k=i-1;
			while(k>=0) {
				while(k>=0 && !trimpoint[k].pos) k--;
				if(k>=0) { // found a point above threshold
					midpos=trimpoint[k].pos;
					if(trimpoint[k].start) midpos--;
					uint leftstart=start;
					int j=k-1;
					while(j>=0) {
						if(trimpoint[j].pos) { leftstart=trimpoint[j].pos;break;}
						j--;
					}
					covleft=get_cov_sign(sno,leftstart-refstart,midpos-refstart,bpcov)/(midpos-leftstart+1); // I add 1 bp coverage to all so I can avoid 0 coverages
					covright=get_cov_sign(sno,midpos+1-refstart,endpos-refstart,bpcov)/(endpos-midpos);
					if((trimpoint[k].start && covleft<localdrop*covright)||(covright<localdrop*covleft)) { // source|sink above threshold
						break; // finished looking for k
					}
					else {
					  //fprintf(stderr,"pos=%d failed global test btw %d-%d midpos=%d covlef=%.1f covright=%.1f\n",trimpoint[k].pos,leftstart,endpos,midpos,covleft,covright);
						trimpoint[k].pos=0;
					}
					k--;
				}
			}
		}
	}

	/*
	{ // DEBUG ONLY
		fprintf(stderr,"All trim points:");
		for(int i=0;i<trimpoint.Count();i++) {
			if(trimpoint[i].start) fprintf(stderr," s:%d",trimpoint[i].pos);
			else fprintf(stderr," e:%d",trimpoint[i].pos);
		}
		fprintf(stderr,"\n");
	}
	*/
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

// cummulative bpcov
CGraphnode *source2guide(int s, int g, int refstart,uint newstart,uint newend, CGraphnode *graphnode,CGraphnode *source,
		GVec<float>* bpcov,GVec<float>& futuretr, int& graphno,CBundlenode *bundlenode,GVec<CGraphinfo> **bundle2graph,
		GPVec<CGraphnode> **no2gnode, int &edgeno) {

	if(graphnode->start+longintronanchor>newstart) { // newstart is very close to graphnode start
		for(int p=0;p<graphnode->parent.Count();p++) if(!graphnode->parent[p]) return(graphnode);
	}

	// compute maxabund
	float leftcov=0;
	float rightcov=0;

	//int os=2-2*s; // other strand

	if(!mergeMode) {
		if(newstart>graphnode->start) {
			uint gstart=graphnode->start;
			if(newstart-gstart > CHI_WIN) {
				gstart=newstart-CHI_WIN;
			}
			//leftcov=get_cov(1,gstart-refstart,newstart-1-refstart,bpcov)- get_cov(os,gstart-refstart,newstart-1-refstart,bpcov);
			leftcov=get_cov_sign(2*s,gstart-refstart,newstart-1-refstart,bpcov);
			leftcov/=newstart-gstart;
		}
		if(newstart<newend) {
			uint gend=newend;
			if(newend-newstart>=CHI_WIN) {
				gend=newstart+CHI_WIN-1;
			}
			//rightcov=get_cov(1,newstart-refstart,gend-refstart,bpcov)-get_cov(os,newstart-refstart,gend-refstart,bpcov);
			rightcov=get_cov_sign(2*s,newstart-refstart,gend-refstart,bpcov);
			rightcov/=gend-newstart+1;
		}
	}

	float maxabund=rightcov-leftcov;
	if(maxabund<trthr) maxabund=trthr;

	if(graphnode->start<=newstart-1) {
		uint prevend=graphnode->end;
		graphnode->end=newstart-1;
		CGraphnode *prevnode=graphnode;
		graphnode=create_graphnode(s,g,newstart,prevend,graphno,bundlenode,bundle2graph,no2gnode);
		graphno++;
		float tmp=prevnode->nodeid;futuretr.Add(tmp);
		tmp=graphnode->nodeid;futuretr.Add(tmp);
		tmp=trthr;futuretr.Add(tmp);
		prevnode->child.Add(graphnode->nodeid); // this node is the child of previous node
		graphnode->parent.Add(prevnode->nodeid); // this node has as parent the previous node
	}
	source->child.Add(graphnode->nodeid);  // this node is the child of source
	graphnode->parent.Add(source->nodeid); // this node has source as parent
	float tmp=graphno-1;
	futuretr.cAdd(0.0);
	futuretr.Add(tmp);
	futuretr.Add(maxabund);
	// COUNT 1 EDGE HERE because the source to guide edge was already included in our count
	edgeno++;

	return(graphnode);

}

// cummulative version
CGraphnode *guide2sink(int s, int g, int refstart,uint newstart,uint newend, CGraphnode *graphnode,CGraphnode *sink,
		GVec<float>* bpcov,GVec<float>& futuretr, int& graphno,CBundlenode *bundlenode,GVec<CGraphinfo> **bundle2graph,
		GPVec<CGraphnode> **no2gnode, int &edgeno) {

	// compute maxabund
	float leftcov=0;
	float rightcov=0;

	//int os=2-2*s;

	if(!mergeMode) {
		if(newstart>=graphnode->start) {
			uint gstart=graphnode->start;
			if(newstart-gstart >= CHI_WIN) {
				gstart=newstart-CHI_WIN+1;
			}
			//leftcov=get_cov(1,gstart-refstart,newstart-refstart,bpcov)-get_cov(os,gstart-refstart,newstart-refstart,bpcov);
			leftcov=get_cov_sign(2*s,gstart-refstart,newstart-refstart,bpcov);
			leftcov/=newstart-gstart+1;
		}
		if(newstart+1<newend) {
			uint gend=newend;
			if(newend-newstart>CHI_WIN) {
				gend=newstart+CHI_WIN;
			}
			//rightcov=get_cov(1,newstart+1-refstart,gend-refstart,bpcov)-get_cov(os,newstart+1-refstart,gend-refstart,bpcov);
			rightcov=get_cov_sign(2*s,newstart+1-refstart,gend-refstart,bpcov);
			rightcov/=gend-newstart;
		}
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

CGraphnode *longtrim(int s, int g, int refstart,int nodeend, int &nls, int &nle, bool &startcov, bool endcov, GVec<CPred> &lstart, GVec<CPred> &lend,
		CGraphnode *graphnode,CGraphnode *source, CGraphnode *sink, GVec<float>& futuretr, int& graphno, GVec<float>* bpcov,
		CBundlenode *bundlenode,GVec<CGraphinfo> **bundle2graph,GPVec<CGraphnode> **no2gnode, int &edgeno) {

	while(nls<lstart.Count() && lstart[nls].predno<(int)graphnode->start) nls++;
	while(nle<lend.Count() && lend[nle].predno<(int)graphnode->start) nle++;
	while((nls<lstart.Count() && lstart[nls].predno<nodeend) || (nle<lend.Count() && lend[nle].predno<nodeend)){
		if(nle>=lend.Count() || (nls<lstart.Count() && lstart[nls].predno<=lend[nle].predno)) { // start comes first
			float tmpcov=0;
			if((startcov || lstart[nls].predno>(int)(graphnode->start+longintronanchor)) &&(endcov || lstart[nls].predno<nodeend+(int)longintronanchor)) { // start and ends can not be too close to a junction
				int startpos=lstart[nls].predno-refstart;
				int winstart=startpos-CHI_THR;
				if(winstart<0) winstart=0;
				int winend=startpos+CHI_THR-1;
				if(winend>=bpcov->Count()) winend=bpcov->Count()-1;
				/*tmpcov=(get_cov(1,startpos,winend,bpcov)-get_cov(2-2*s,startpos,winend,bpcov)-
						get_cov(1,winstart,startpos-1,bpcov)+get_cov(2-2*s,winstart,startpos-1,bpcov))/(DROP*CHI_THR);*/
				tmpcov=(get_cov_sign(2*s,startpos,winend,bpcov)-get_cov_sign(2*s,winstart,startpos-1,bpcov))/(DROP*CHI_THR);
			}
			if(tmpcov<=0 && lstart[nls].cov<0) tmpcov=ERROR_PERC; // to re-estimate later in process_transfrags
			if(tmpcov>0) {
				tmpcov+=trthr;
				uint prevend=graphnode->end;
				graphnode->end=lstart[nls].predno-1;
				CGraphnode *prevnode=graphnode;
				graphnode=create_graphnode(s,g,lstart[nls].predno,prevend,graphno,bundlenode,bundle2graph,no2gnode);
				graphnode->hardstart=true;
				graphno++;
				source->child.Add(graphnode->nodeid);  // this node is the child of source
				graphnode->parent.Add(source->nodeid); // this node has source as parent
				prevnode->child.Add(graphnode->nodeid); // this node is the child of previous node
					graphnode->parent.Add(prevnode->nodeid); // this node has as parent the previous node
				float tmp=graphno-1;
				futuretr.cAdd(0.0);
				futuretr.Add(tmp);
				//tmp=lstart[nls].cov+trthr;
				futuretr.Add(tmpcov);
				//		graphnode->start,lstart[nls].cov,get_cov(2*s,startpos,startpos+CHI_THR-1,bpcov),get_cov(1,startpos,startpos+CHI_THR-1,bpcov),get_cov(2-2*s,startpos,startpos+CHI_THR-1,bpcov),
				//		get_cov(2*s,startpos-CHI_THR,startpos-1,bpcov),get_cov(1,startpos-CHI_THR,startpos-1,bpcov),get_cov(2-2*s,startpos-CHI_THR,startpos-1,bpcov));
				tmp=prevnode->nodeid;futuretr.Add(tmp);
				tmp=graphnode->nodeid;futuretr.Add(tmp);
				tmp=trthr;futuretr.Add(tmp);
				// COUNT 2 EDGES HERE
				edgeno+=2;
				startcov=false;
			}
			nls++;
		}
		else if(nls>=lstart.Count() || (nle<lend.Count() && lend[nle].predno<lstart[nls].predno)) { // end comes first
			float tmpcov=0;
			if((!startcov || lend[nle].predno>(int)(graphnode->start+longintronanchor)) &&(!endcov || lend[nle].predno<nodeend+(int)longintronanchor)) {
				int endpos=lend[nle].predno-refstart;
				int winstart=endpos-CHI_THR+1;
				if(winstart<0) winstart=0;
				int winend=endpos+CHI_THR;
				if(winend>=bpcov->Count()) winend=bpcov->Count()-1;
				/*tmpcov=(get_cov(1,winstart,endpos,bpcov)-get_cov(2-2*s,winstart,endpos,bpcov)-
						get_cov(1,endpos+1,winend,bpcov)+get_cov(2-2*s,endpos+1,winend,bpcov))/(DROP*CHI_THR);*/
				tmpcov=(get_cov_sign(2*s,winstart,endpos,bpcov)-get_cov_sign(2*s,endpos+1,winend,bpcov))/(DROP*CHI_THR);
			}
			if(tmpcov<=0 && lend[nle].cov<0) tmpcov=ERROR_PERC; // to re-estimate later in process_transfrags
			if(tmpcov>0) {
				tmpcov+=trthr;
				float tmp=graphno-1;
				uint prevend=graphnode->end;
				graphnode->end=lend[nle].predno;
				CGraphnode *prevnode=graphnode;
				graphnode->hardend=true;
				graphnode=create_graphnode(s,g,lend[nle].predno+1,prevend,graphno,bundlenode,bundle2graph,no2gnode);
				graphno++;
				prevnode->child.Add(graphnode->nodeid); // this node is the child of previous node
				graphnode->parent.Add(prevnode->nodeid); // this node has as parent the previous node
				sink->parent.Add(prevnode->nodeid); // prevnode is the parent of sink
				// remember to create transfrag as well -> I don't know the gno yet, so I can not create it here
				futuretr.Add(tmp);
				futuretr.cAdd(-1.0);
				//tmp=lend[nle].cov+trthr;
				futuretr.Add(tmpcov);
				tmp=prevnode->nodeid;futuretr.Add(tmp);
				tmp=graphnode->nodeid;futuretr.Add(tmp);
				tmp=trthr;futuretr.Add(tmp);
				// COUNT 2 EDGES HERE
				edgeno+=2;
				startcov=true;
			}
			nle++;
		}
	}
	return(graphnode);
}

CGraphnode *trimnode(int s, int g, int refstart,uint newend, CGraphnode *graphnode,CGraphnode *source, CGraphnode *sink, GVec<float>* bpcov,
		GVec<float>& futuretr, int& graphno,CBundlenode *bundlenode,GVec<CGraphinfo> **bundle2graph,GPVec<CGraphnode> **no2gnode, int &edgeno) {

	uint sourcestart=0;
	uint sinkend=0;
	float sinkabundance=0;
	float sourceabundance=0;
	//float sourcecovleft=0;
	//float sourcecovright=0;
	//float sinkcovleft=0;
	//float sinkcovright=0;
	//find_trims(refstart,2*s,graphnode->start,newend,bpcov,sourcestart,sourcecovleft,sourcecovright,sinkcovleft,sinkcovright);
	find_trims_wsign(refstart,2*s,graphnode->start,newend,bpcov,sourcestart,sourceabundance,sinkend,sinkabundance);
	//find_trims(refstart,graphnode->start,newend,bpcov,sourcestart,sourceabundance,sinkend,sinkabundance);

	if(sourcestart < sinkend) { // source trimming comes first

		if(sourcestart) { // there is evidence of graphnode trimming from source
			graphnode->end=sourcestart-1;
			CGraphnode *prevnode=graphnode;
			graphnode=create_graphnode(s,g,sourcestart,newend,graphno,bundlenode,bundle2graph,no2gnode);
			graphnode->hardstart=true;
			graphno++;
			source->child.Add(graphnode->nodeid);  // this node is the child of source
			graphnode->parent.Add(source->nodeid); // this node has source as parent
			prevnode->child.Add(graphnode->nodeid); // this node is the child of previous node
			graphnode->parent.Add(prevnode->nodeid); // this node has as parent the previous node
			float tmp=graphno-1;
			futuretr.cAdd(0.0);
			futuretr.Add(tmp);
			sourceabundance+=trthr;futuretr.Add(sourceabundance);
			//sourcecovright+=trthr;futuretr.Add(sourcecovright);
			// start TRY linking to end also here
			//tmp=prevnode->nodeid;futuretr.Add(tmp);futuretr.cAdd(-1.0);tmp=sourcecovleft+trthr;futuretr.Add(tmp);sink->parent.Add(prevnode->nodeid);edgeno++;
			// end TRY linking to end also here
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
		sink->parent.Add(prevnode->nodeid); // prevnode is the parent of sink
		prevnode->hardend=true;
		// remember to create transfrag as well -> I don't know the gno yet, so I can not create it here
		futuretr.Add(tmp);
		futuretr.cAdd(-1.0);
		sinkabundance+=trthr;futuretr.Add(sinkabundance);
		//sinkcovleft+=trthr;futuretr.Add(sinkcovleft);
		// start TRY linking to end also here
		//tmp=graphnode->nodeid;futuretr.cAdd(0.0);futuretr.Add(tmp);tmp=sinkcovright+trthr;futuretr.Add(tmp);source->child.Add(graphnode->nodeid);edgeno++;
		// end TRY linking to end also here
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
			sink->parent.Add(prevnode->nodeid); // prevnode is the parent of sink
			prevnode->hardend=true;
			// remember to create transfrag as well -> I don't know the gno yet, so I can not create it here
			float tmp=graphno-2;
			futuretr.Add(tmp);
			futuretr.cAdd(-1.0);
			sinkabundance+=trthr;futuretr.Add(sinkabundance);
			//sinkcovleft+=trthr;futuretr.Add(sinkcovleft);
			// start TRY linking to end also here
			//tmp=graphnode->nodeid;futuretr.cAdd(0.0);futuretr.Add(tmp);tmp=sinkcovright+trthr;futuretr.Add(tmp);source->child.Add(graphnode->nodeid);edgeno++;
			// end TRY linking to end also here
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
		graphnode->hardstart=true;
		graphnode->parent.Add(source->nodeid); // this node has source as parent
		prevnode->child.Add(graphnode->nodeid); // this node is the child of previous node
		graphnode->parent.Add(prevnode->nodeid); // this node has as parent the previous node
		float tmp=graphno-1;
		futuretr.cAdd(0.0);
		futuretr.Add(tmp);
		sourceabundance+=trthr;futuretr.Add(sourceabundance);
		//sourcecovright+=trthr;futuretr.Add(sourcecovright);
		// start TRY linking to end also here
		//tmp=prevnode->nodeid;futuretr.Add(tmp);futuretr.cAdd(-1.0);tmp=sourcecovleft+trthr;futuretr.Add(tmp);sink->parent.Add(prevnode->nodeid);edgeno++;
		// end TRY linking to end also here
		tmp=prevnode->nodeid;futuretr.Add(tmp);
		tmp=graphnode->nodeid;futuretr.Add(tmp);
		tmp=trthr;futuretr.Add(tmp);
		// COUNT 2 EDGES HERE
		edgeno+=2;
	}
	// else both source and sink trimming are not present

	return(graphnode);
}


CGraphnode *trimnode_all(int s, int g, int refstart,uint newend, CGraphnode *graphnode,CGraphnode *source, CGraphnode *sink, GVec<float>* bpcov,
		GVec<float>& futuretr, int& graphno,CBundlenode *bundlenode,GVec<CGraphinfo> **bundle2graph,GPVec<CGraphnode> **no2gnode, int &edgeno) {

	GVec<CTrimPoint> trimpoint;
	find_all_trims(refstart,2*s,graphnode->start,newend,bpcov,trimpoint);
	for(int i=0;i<trimpoint.Count();i++) if(trimpoint[i].pos){
		if(trimpoint[i].start) { // source trim
			graphnode->end=trimpoint[i].pos-1;
			//fprintf(stderr,"Create source trim:%d-%d and %d-%d\n",graphnode->start,graphnode->end,trimpoint[i].pos,newend);
			CGraphnode *prevnode=graphnode;
			graphnode=create_graphnode(s,g,trimpoint[i].pos,newend,graphno,bundlenode,bundle2graph,no2gnode);
			graphnode->hardstart=true;
			graphno++;
			source->child.Add(graphnode->nodeid);  // this node is the child of source
			graphnode->parent.Add(source->nodeid); // this node has source as parent
			prevnode->child.Add(graphnode->nodeid); // this node is the child of previous node
			graphnode->parent.Add(prevnode->nodeid); // this node has as parent the previous node
			float tmp=graphno-1;
			futuretr.cAdd(0.0);
			futuretr.Add(tmp);
			float sourceabundance=trimpoint[i].abundance+trthr;futuretr.Add(sourceabundance);
			tmp=prevnode->nodeid;futuretr.Add(tmp);
			tmp=graphnode->nodeid;futuretr.Add(tmp);
			tmp=trthr;futuretr.Add(tmp);
			// COUNT 2 EDGES HERE
			edgeno+=2;
		}
		else { // this is sink
			graphnode->end=trimpoint[i].pos;
			CGraphnode *prevnode=graphnode;
			graphnode=create_graphnode(s,g,trimpoint[i].pos+1,newend,graphno,bundlenode,bundle2graph,no2gnode);
			graphno++;
			prevnode->child.Add(graphnode->nodeid); // this node is the child of previous node
			graphnode->parent.Add(prevnode->nodeid); // this node has as parent the previous node
			sink->parent.Add(prevnode->nodeid); // prevnode is the parent of sink
			prevnode->hardend=true;
			// remember to create transfrag as well -> I don't know the gno yet, so I can not create it here
			float tmp=graphno-2;
			futuretr.Add(tmp);
			futuretr.cAdd(-1.0);
			float sinkabundance=trimpoint[i].abundance+trthr;futuretr.Add(sinkabundance);
			tmp=prevnode->nodeid;futuretr.Add(tmp);
			tmp=graphnode->nodeid;futuretr.Add(tmp);
			tmp=trthr;futuretr.Add(tmp);
			// COUNT 2 EDGES HERE
			edgeno+=2;
		}
	}

	return(graphnode);
}

inline int edge(int min, int max, int gno) {
	//return((gno-1)*min-min*(min-1)/2+max-min); // this should be changed if source to node edges are also stored
	return((gno-1)*(min+1)-min*(min-1)/2+max-min); // this includes source to node edges
}

GBitVec traverse_dfs(int s,int g,CGraphnode *node,CGraphnode *sink,GBitVec parents,int gno, GVec<bool>& visit,
		GPVec<CGraphnode> **no2gnode,GPVec<CTransfrag> **transfrag, int &edgeno,GIntHash<int> **gpos,int &lastgpos){

	//fprintf(stderr,"Traverse node %d\n",node->nodeid);

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
			/*if(mixedMode) {
				tr->abundance*=2;
			}*/

			transfrag[s][g].Add(tr);
			if(mixedMode) { // I need to add a long read as well
				CTransfrag *longtr=new CTransfrag(nodes,trpat,trthr);
				longtr->longread=true;
				transfrag[s][g].Add(longtr);
			} else
			if(longreads) // || mixedMode)
				transfrag[s][g].Last()->longread=true;

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

			/*if(mixedMode) {
				tr->abundance*=2;
			}*/

			transfrag[s][g].Add(tr);
			if(mixedMode) { // I need to add a long read as well
				CTransfrag *longtr=new CTransfrag(nodes,trpat,trthr);
				longtr->longread=true;
				transfrag[s][g].Add(longtr);
			}
			else
			if(longreads)// || mixedMode)
				transfrag[s][g].Last()->longread=true;
			n++;
	    }
		/*
		fprintf(stderr,"Add %d children of node %d (%d-%d): ",n,node->nodeid,node->start,node->end);
		for(int i=0;i<n;i++) fprintf(stderr," %d",node->child[i]);
		fprintf(stderr,"\n");
		*/

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
				childparents[lastgpos]=1;
				node->childpat[lastgpos]=1;
				lastgpos++;
			}

			//fprintf(stderr,"Call for child %d with id=%d\n",node->child[i],no2gnode[s][g][node->child[i]]->nodeid);
	    	node->childpat = node->childpat | traverse_dfs(s,g,no2gnode[s][g][node->child[i]],sink,childparents,gno,visit,no2gnode,transfrag,edgeno,gpos,lastgpos);
	    }
	} // end else from if(visit[node->nodeid])

	GBitVec children = node->childpat;
	children[node->nodeid]=1;

	return(children);
}


void delete_Ginterval(GInterval *interv){
	if(interv) {
		if(interv->next) delete_Ginterval(interv->next);
		delete interv;
	}
}

GInterval *add_interval(uint start,uint end,GInterval *guidecov) {

	if(guidecov==NULL) { // no intervals defined yet
		guidecov=new GInterval(start,end);
	}
	else {
		if(end<guidecov->start-1) { // interval before guidecov; if(end==guidecov->start-1) it's like the intervals overlap
			GInterval *interval=new GInterval(start,end,guidecov);
			guidecov=interval;
		}
		else { // end>=guidecov->start-1 I need to place the new interval
			GInterval *lastinterv=NULL; // last interval before overlap
			GInterval *interval=guidecov;
			while(interval && start>interval->end+1) { lastinterv=interval; interval = interval->next;}
			if(interval) { // start <= interval->end+1; interval could be guidecov here
				if(end<interval->start-1) { // (start,end) come before interval but does not intersect lastinterv
					GInterval *newinterval=new GInterval(start,end,interval);
					lastinterv->next=newinterval; // lastinterv can not be null
				}
				else { // (start,end) intersect interval
					if(start<interval->start) interval->start=start;
					lastinterv=interval; // start>=lastinterv->start, end>=lastinterv->start-1
					if(end>lastinterv->end) lastinterv->end=end;
					interval=interval->next;
					while(interval && end>=interval->start-1) { // (start,end) intersects interval
						if(lastinterv->end<interval->end) lastinterv->end=interval->end;
						lastinterv->next=interval->next;
						delete interval;
						interval=lastinterv->next;
					}
				}
			}
			else { // (start,end) at the end of intervals
				interval=new GInterval(start,end);
				lastinterv->next=interval;
			}
		}
	}
	return(guidecov);
}

bool overlap_interval(uint start,uint end,GInterval *guidecov) {

	if(guidecov==NULL || end<guidecov->start) return false; // interval before guide

	while(guidecov && start>guidecov->end) guidecov=guidecov->next;

	if(guidecov) { // start<=guidecov->end
		if(end<guidecov->start) return false;
		// end >= guidecov->start
		return true;
	}

	return false;
}



int count_kept_nodes(GVec<CGNode> &node,int gno,GPVec<CGraphnode>& no2gnode) {
	int new_graphno=1; // first include source and first node because this is always kept

	for(int n=1;n<gno;n++) {
		if(node[n].keep) {
			bool update_count=true; // by default I update the count
			if(!node[n].future && node[n-1].keep && no2gnode[n-1]->end+1==no2gnode[n]->start) { // only in the case that previous node is adjacent I can collapse nodes; never merge a future node

				if(no2gnode[n-1]->child.Count()==1 && no2gnode[n-1]->child[0]==n && no2gnode[n]->parent.Count()==1) { // the adjacent parent is the only parent, and this node is the only child of it
					update_count=false;
				}
			}

			if(update_count) {
				new_graphno++;
				node[n].merge=false;
			}
			else node[n].merge=true;
		}
	}
	return(new_graphno);
}


int gjuncCmp(const pointer p1, const pointer p2) {
	CGJunc *a=(CGJunc *)p1;
	CGJunc *b=(CGJunc *)p2;

	if(a->goodcov>b->goodcov) return 1; // less coverage come first
	if(a->goodcov<b->goodcov) return -1;

	if(a->cov>b->cov) return 1; // less coverage come first
	if(a->cov<b->cov) return -1;

	return 0;
}


void delete_connection(GVec<CGNode> &node,int n1, int n2,  GPVec<CGraphnode>& no2gnode,int &edgeno,int &estimate_graphno) {
	int i=0;
	while(i<no2gnode[n1]->child.Count()) {
		if(no2gnode[n1]->child[i]==n2) {
			no2gnode[n1]->child.Delete(i);
			break;
		}
		i++;
	}

	if(!no2gnode[n1]->child.Count() && !no2gnode[n1]->parent.Count()){
		node[n1].keep=false;
		estimate_graphno--;
	}
	else if(no2gnode[n1]->child.Count()==1 && no2gnode[n1]->child[0]==n1+1 && no2gnode[n1]->end==no2gnode[n1+1]->start-1) { // node n1 can be merged with next one because there is no other junction from n1
		estimate_graphno--;
	}

	i=0;
	while(i<no2gnode[n2]->parent.Count()) {
		if(no2gnode[n2]->parent[i]==n1) {
			no2gnode[n2]->parent.Delete(i);
			break;
		}
		i++;
	}

	if(!no2gnode[n2]->child.Count() && !no2gnode[n2]->parent.Count()) {
		node[n2].keep=false;
		estimate_graphno--;
	}
	else if(no2gnode[n2]->parent.Count()==1 && no2gnode[n2]->parent[0]==n2-1 && no2gnode[n2-1]->end==no2gnode[n2]->start-1) { // node n2 can be merged with previous one because there is no other junction to n2
		estimate_graphno--;
	}

	edgeno--;
}

int prune_graph_nodes(int graphno,int s,int g,GVec<CGraphinfo> **bundle2graph, int bnodecount,
		GPVec<CGraphnode> **no2gnode,GList<CJunction>& junction,int &edgeno,GVec<float> &futuretr,CGraphnode *sink){

	//fprintf(stderr,"start with edgeno=%d\n",edgeno);

	GVec<CGJunc> kjunc; // scan all junctions and keep the ones that are in the graph and I can potentially delete
	int n=1; // first node in graph
	for(int i=0;i<junction.Count();i++) {
		if(junction[i]->start>no2gnode[s][g][graphno-1]->start) break; // all future junctions are out of the graph
		if(!junction[i]->guide_match && junction[i]->strand+1==2*s && junction[i]->start>=no2gnode[s][g][1]->end && junction[i]->end<=no2gnode[s][g][graphno-1]->start) { // junction in graph
			while(n<graphno && no2gnode[s][g][n]->end<junction[i]->start) n++;
			if(n==graphno) break;
			if(no2gnode[s][g][n]->end==junction[i]->start) {
				int n2=0;
				for(int c=0;c<no2gnode[s][g][n]->child.Count();c++) {
					if(no2gnode[s][g][no2gnode[s][g][n]->child[c]]->start==junction[i]->end) {
						n2=no2gnode[s][g][n]->child[c];
						break;
					}
				}
				if(n2) {
					CGJunc jc(n,n2,junction[i]->nreads,junction[i]->nreads_good);
					kjunc.Add(jc);
				}
			}
		}
	}

	kjunc.Sort(gjuncCmp);

	GVec<CGNode> node;
	for(int i=0;i<graphno;i++) {
		bool sinklink=false;
		if(!no2gnode[s][g][i]->child.Count()) sinklink=true; // I could use this to parse nodes that do not link anywhere
		CGNode inode(i,sinklink);
		node.Add(inode); // node[0] has id 1
	}

	for(int i=0;i<futuretr.Count();i+=3) {
		int n2=int(futuretr[i+1]);
		if(n2>0) node[n2].future=true;
		else node[futuretr[i]].last=true; // when n2==-1 there is a link to sink from n1
	}

	/*
	for(int i=1;i<graphno;i++) {
		//fprintf(stderr,"Node[%d]:%d-%d with parents:",i,no2gnode[s][g][i]->start,no2gnode[s][g][i]->end);
		for(int p=0;p<no2gnode[s][g][i]->parent.Count();p++) fprintf(stderr," %d",no2gnode[s][g][i]->parent[p]);
		fprintf(stderr," and children:");
		for(int c=0;c<no2gnode[s][g][i]->child.Count();c++) fprintf(stderr," %d",no2gnode[s][g][i]->child[c]);
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"there are %d junctions\n",kjunc.Count());
	*/


	int i=0;
	/*** code to increase threshold requirements for junctions
	//int new_graphno=graphno;
	while(i<kjunc.Count() && kjunc[i].goodcov<junctionthr+1) { // require at least one more read for support
		delete_connection(node,kjunc[i].leftnode,kjunc[i].rightnode,no2gnode[s][g],edgeno);
		i++;
	}
	end code ***/

	/*** code that estimates new graphno ***/
	int estimate_graphno=graphno;
	while(estimate_graphno > allowed_nodes && i<kjunc.Count()) {
		delete_connection(node,kjunc[i].leftnode,kjunc[i].rightnode,no2gnode[s][g],edgeno,estimate_graphno);
		i++;
	}
	/*** end code ***/

	int new_graphno=count_kept_nodes(node,graphno,no2gnode[s][g]);

	/*
	while(new_graphno > allowed_nodes && i<kjunc.Count()) {
		// delete less covered junctions as long as they have the same coverage
		double jcov=kjunc[i].cov+kjunc[i].goodcov;
		int n=0;
		while(i<kjunc.Count() && kjunc[i].cov+kjunc[i].goodcov==jcov) {
			delete_connection(node,kjunc[i].leftnode,kjunc[i].rightnode,no2gnode[s][g]);
			i++;
			n++;
		}
		new_graphno=count_kept_nodes(node,graphno,no2gnode[s][g]);
	}
    */

	// there was a part here to re-create lost children likns, e.g. for 1->2->3 delete link from 1->2 and from 2->3. 1 might still want to link to 3

	// I might need to recreate links to source though

	// re-number nodes that are kept from 1 to new_graphno instead
	i=1;
	for(int n=1;n<graphno;n++) {
		if(node[n].keep) {
			CGraphnode *nnode=no2gnode[s][g][n];
			if(node[n].merge) { // this node needs to be merged
				int m=n-1;
				while(node[m].merge) m--;
				CGraphnode *mnode=no2gnode[s][g][m];
				mnode->end=nnode->end;
				mnode->child=nnode->child;
				node[n].keep=false;
				nnode=mnode;
			}
			else {
				node[n].id=i++;

				// also create links to source if needed --> the only one that does not recreate itself in traverse_dfs
				if(!nnode->parent.Count()) { // I need to link to source
					CGraphnode *source=no2gnode[s][g][0];
					source->child.Add(nnode->nodeid);
					nnode->parent.Add(0);
				}
			}
		}
	}


	// redo the structures in the graph to reflect the new node numbering
	// structures to deal with: graphnode->child, graphnode->parent (and edgeno), no2gnode, bundle2graph -> this also needs to be addressed in get_fragment/read_pattern

	for(int n=graphno-1;n>0;n--) { // I am doing it from last node to first in order for the deletions to be meaningful
		CGraphnode *inode=no2gnode[s][g][n];
		if(node[n].keep) {
			for(int c=inode->child.Count()-1;c>=0;c--) {
				if(node[inode->child[c]].keep) { // if I keep it I have to give it the new id
					inode->child[c]=node[inode->child[c]].id;
				}
				else { // I need to delete the child
					inode->child.Delete(c);
					edgeno--; // delete edge btwn node and child
				}
			}
			for(int p=inode->parent.Count()-1;p>=0;p--) {
				if(inode->parent[p]) { // if the parent is not the source
					if(node[inode->parent[p]].keep) { // if I keep it I have to give it the new id
						inode->parent[p]=node[inode->parent[p]].id;
					}
					else if(node[inode->parent[p]].merge) {
						int m=inode->parent[p]-1;
						while(node[m].merge) m--;
						inode->parent[p]=node[m].id;
					}
					else { // I need to delete the parent
						inode->parent.Delete(p);
						edgeno--;
					}
				}
			}
			inode->nodeid=node[n].id;
		}
		else { // we need to delete the node
			// need to update edge count first
			for(int c=0;c<inode->child.Count();c++) {
				if(!node[inode->child[c]].keep) edgeno--; // link to a deleted node needs to be removed
			}
			no2gnode[s][g].Delete(n); // this deletes the actual element too if gpvec was not created with free element =false
		}
	}

	// deal with the source
	CGraphnode *source=no2gnode[s][g][0];
	for(int c=source->child.Count()-1;c>=0;c--) {
		if(node[source->child[c]].keep) { // if I keep it I have to give it the new id
			source->child[c]=node[source->child[c]].id;
		}
		else { // I need to delete the child --> this should never happen as I should never delete a link to source
			source->child.Delete(c);
			edgeno--; // delete edge btwn node and child
		}
	}

	// deal with the sink --> do I have sink??
	for(int p=sink->parent.Count()-1;p>=0;p--) {
		if(node[sink->parent[p]].keep) { // if I keep it I have to give it the new id
			sink->parent[p]=node[sink->parent[p]].id;
		}
		else { // I need to delete the parent
			sink->parent.Delete(p);
			edgeno--; // delete edge btwn node and child
		}
	}

	for(int i=0;i<futuretr.Count();i+=3) {
		int n1=int(futuretr[i]);
		int n2=int(futuretr[i+1]);

		if(n1>0) {
			if(node[n1].keep || node[n1].merge) {
				while(node[n1].merge) n1--;
				if(n2>0) {
					if(node[n2].keep) {
						futuretr[i]=(float)node[n1].id;
						futuretr[i+1]=(float)node[n2].id;
					}
					else { futuretr[i]=-1; edgeno--;}
				}
				else futuretr[i]=(float)node[n1].id;
			}
			else { futuretr[i]=-1; edgeno--;} // no longer interested in this futuretr
		}
		else { // n1 is the source
			if(node[n2].keep) {
				futuretr[i+1]=(float)node[n2].id;
			}
			else { futuretr[i]=-1;edgeno--;}
		}

	}

	for(i=0;i<bnodecount;i++) {
		for(int b=bundle2graph[s][i].Count()-1;b>=0;b--) {
			if(bundle2graph[s][i][b].ngraph==g) {
				if(node[bundle2graph[s][i][b].nodeno].keep) {
					bundle2graph[s][i][b].nodeno=node[bundle2graph[s][i][b].nodeno].id;
				}
				else {
					bundle2graph[s][i].Delete(b);
				}
			}
		}
	}

	/*
	{ //DEBUG ONLY
		for(int i=0;i<new_graphno;i++) {
			fprintf(stderr,"Node %d with parents:",i);
			for(int p=0;p<no2gnode[s][g][i]->parent.Count();p++) fprintf(stderr," %d",no2gnode[s][g][i]->parent[p]);
			fprintf(stderr," and children:");
			for(int c=0;c<no2gnode[s][g][i]->child.Count();c++) fprintf(stderr," %d",no2gnode[s][g][i]->child[c]);
			fprintf(stderr,"\n");
		}
	}
	fprintf(stderr,"new edgeno=%d\n",edgeno);
	*/

	return(new_graphno);
}

int create_graph(int refstart,int s,int g,CBundle *bundle,GPVec<CBundlenode>& bnode,
		GList<CJunction>& junction,GList<CJunction>& ejunction,GVec<CGraphinfo> **bundle2graph,
		GPVec<CGraphnode> **no2gnode,GPVec<CTransfrag> **transfrag,GIntHash<int> **gpos,BundleData* bdata,
		int &edgeno,int &lastgpos,GArray<GEdge>& guideedge, int refend=0){

	GVec<float>* bpcov = bdata ? bdata->bpcov : NULL; // I might want to use a different type of data for bpcov to save memory in the case of very long bundles

	CGraphnode* source=new CGraphnode(0,0,0);
	no2gnode[s][g].Add(source);
	CGraphnode* sink=new CGraphnode();

	int njunctions=junction.Count();

	//fprintf(stderr,"Start graph[%d][%d] with %d edgeno and lastgpos=%d\n",s,g,edgeno,lastgpos);

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
	CBundlenode *bundlenode=bnode[bundle->startnode];
	while(nge<guideedge.Count() && bundlenode!=NULL) {
		uint start=guideedge[nge].val;
		uint end=guideedge[nge].endval;
		if(start>end) Gswap(start,end);
		if(bundlenode->end<start) bundlenode=bundlenode->nextnode;
		else if(guideedge[nge].strand==s && bundlenode->start<=end) { // bundlenode->end>=start
			nge=0;
			processguide=true;
			break;
		}
		else nge++;
	}
	bundlenode=bnode[bundle->startnode];

	int njs=0; // index of sorted junction starts
	int nje=0; // index of sorted junction ends

	int graphno=1; // number of nodes in graph
	//GHash<GVec<int>* > ends; // keeps ids of all nodes ending at a certain position; OR ALL NODES THAT ARE LINKED BY JUNCTIONS TO A CERTAIN POSITION
    GIntHash< GVec<int>* > ends;
	GVec<float> futuretr;

	if(mergeMode) { // I have a bunch of junctions at the start for which I need to create ends

		while(njs<njunctions && !junction[njs]->start ) { // remember ends here for source node
			if((junction[njs]->strand+1) == 2*s) {
				//GStr je((int)junction[njs]->end);
				//GVec<int> *e=ends[je.chars()];
				GVec<int> *e=ends[junction[njs]->end];
				if(!e) {
					e = new GVec<int>();
					//ends.Add(je.chars(),e);
					ends.Add(junction[njs]->end, e);
				}
				e->cAdd(0);
			}
			njs++;
		}

	}

	//int seenjunc=0;

	int f=0; // feature index

	uint bundle_end=bnode[bundle->lastnodeid]->end;
	while(bundlenode!=NULL) {

		//fprintf(stderr,"process bundlenode %d-%d:%d bpcov_count=%d refstart=%d\n",bundlenode->start,bundlenode->end,s,bpcov->Count(),refstart);

	    uint currentstart=bundlenode->start; // current start is bundlenode's start
	    uint endbundle=bundlenode->end; // initialize end with bundlenode's end for now

	    int end=0;
	    while(nje<njunctions && ejunction[nje]->end<=currentstart) { // read all junction ends at or before the current start -> assuming there are any (at this point, smaller junction ends should not be relevant to this bundle/currentstart
	      if(ejunction[nje]->end==currentstart && (ejunction[nje]->strand+1) == 2*s) { // junction ends at current start and is on the same strand and not deleted
	    	  end=1;
	      }
	      nje++;
	    }

    	GVec<CPred> lstart; // CPred: prediction point class
    	GVec<CPred> lend;
    	int fs=-1; // first start feature index in lstart
    	int fe=-1; // first end feature index in lend
	    if(longreads) { // this can not work the same in mixedMode!

	    	GPVec<GPtFeature>& feature = bdata->ptfs; // these are point features (confirmed starts/stops)

	    	if(feature.Count()) { // I already know the features
	    		while(f<feature.Count() && feature[f]->coord<currentstart) f++;
	    	}
	    	if(endbundle-currentstart<=CHI_WIN+CHI_THR) { // I only consider the given features for this case
	    		while(f<feature.Count() && feature[f]->coord<=endbundle) {
	    			if(feature[f]->strand==2*s-1) {
	    				if(feature[f]->ftype==GPFT_TSS) { // this is a start
	    					CPred p(feature[f]->coord,-1); // start here -> a -1 indicates this is a given point feature
	    					if(fs<0) fs=lstart.Count();
	    					lstart.Add(p);
	    				}
	    				else if(feature[f]->ftype==GPFT_CPAS) { // this is an end
	    					CPred p(feature[f]->coord,-1); // start here -> a -1 indicates this is a given point feature
	    					if(fe<0) fe=lend.Count();
	    					lend.Add(p);
	    				}
	    			}
	    			f++;
	    		}
	    	}
	    	else { // populate starts and ends for the bundlenode here if I have enough space to determine starts/ends in the bundlenode

	    		while(f<feature.Count() && feature[f]->coord<currentstart+CHI_WIN+CHI_THR+longintronanchor) { // see if there are any features before the start
	    			if(feature[f]->strand==2*s-1) {
	    				if(feature[f]->ftype==GPFT_TSS) { // this is a start
	    					CPred p(feature[f]->coord,-1); // start here -> a -1 indicates this is a given point feature
	    					if(fs<0) fs=lstart.Count();
	    					lstart.Add(p);
	    				}
	    				else if(feature[f]->ftype==GPFT_CPAS) { // this is an end
	    					CPred p(feature[f]->coord,-1); // start here -> a -1 indicates this is a given point feature
	    					if(fe<0) fe=lend.Count();
	    					lend.Add(p);
	    				}
	    			}
	    			f++;
	    		}

	    		GVec<float> diffval;
	    		float sumstartleft=0;
	    		float sumendleft=0;
	    		float sumstartright=0;
	    		float sumendright=0;
	    		//float lastcov=get_cov(2*s,longintronanchor-1,longintronanchor-1,bpcov); // this only computes signed coverage (should I include neutral coverage here too?)
	    		float lastcov=get_cov_sign(2*s,currentstart-refstart+longintronanchor-1,currentstart-refstart+longintronanchor-1,bpcov); // this includes neutral coverage
	    		for(int i=currentstart+longintronanchor;i<(int)currentstart+CHI_THR+(int)longintronanchor;i++) {
	    			//float icov=get_cov(2*s,i-refstart,i-refstart,bpcov);
	    			float icov=get_cov_sign(2*s,i-refstart,i-refstart,bpcov);
	    			float diff=icov-lastcov; // compute how many reads start at position i
	    			diffval.Add(diff);
	    			if(diff>0) sumstartleft+=diff;
	    			else sumendleft-=diff;
	    			lastcov=icov;
	    		}
	    		for(int i=currentstart+CHI_THR+longintronanchor;i<(int)currentstart+CHI_WIN+(int)longintronanchor;i++) {
	    			//float icov=get_cov(2*s,i-refstart,i-refstart,bpcov);
	    			float icov=get_cov_sign(2*s,i-refstart,i-refstart,bpcov);
	    			float diff=icov-lastcov;
	    			diffval.Add(diff);
	    			if(diff>0) sumstartright+=diff;
	    			else sumendright-=diff;
	    			lastcov=icov;
	    		}
	    		int covend=endbundle-refstart;
	    		if(bpcov->Count()<(int)endbundle) covend=bpcov->Count();
	    		covend-=longintronanchor;
	    		for(int i=currentstart+CHI_WIN+longintronanchor-refstart;i<covend;i++) {

	    			// add features at position i
	    			if((uint)i+CHI_THR+refstart<endbundle) while(f<feature.Count() && feature[f]->coord==(uint)i+CHI_THR+refstart) { // see if there are any features before the start
	    				if(feature[f]->strand==2*s-1) {
	    					if(feature[f]->ftype==GPFT_TSS) { // this is a start
	    						CPred p(feature[f]->coord,-1); // start here -> a -1 indicates this is a given point feature
	    						lstart.Add(p);
	    					}
	    					else if(feature[f]->ftype==GPFT_CPAS) { // this is an end
	    						CPred p(feature[f]->coord,-1); // start here -> a -1 indicates this is a given point feature
	    						lend.Add(p);
	    					}
	    				}
	    				f++;
	    			}

	    			int m=(i-CHI_THR-longintronanchor-currentstart+refstart)%CHI_WIN;
					if(diffval[m]>1/ERROR_PERC) { // starts here
	    				if(sumstartleft<sumstartright*ERROR_PERC) { // potential start here
	    					bool addstart=true;
	    					int istart=refstart+i-CHI_THR;
	    					if(lstart.Count()) {
	    						int j=lstart.Count()-1;
	    						if(abs(istart-lstart[j].predno)<CHI_WIN) {
	    							addstart=false;
	    							if(lstart[j].cov>0 && diffval[m]>lstart[j].cov) {
	    								lstart[j].predno=istart;
	    								lstart[j].cov=diffval[m];
	    							}
	    						}
	    					}
	    					if(addstart) {
	    						CPred p(istart,diffval[m]);
	    						lstart.Add(p);
	    					}
	    				}
	    			}
					//float icov=get_cov(2*s,i,i,bpcov);
					float icov=get_cov_sign(2*s,i,i,bpcov);
	    			int p=(m+CHI_THR)%CHI_WIN;
	    			float diff=icov-lastcov;
	    			if(diffval[p]>0) sumstartleft-=diffval[p];
	    			else sumendleft+=diffval[p];
	    			if(diffval[m]>0) {
	    				sumstartleft+=diffval[m];
	    				sumstartright-=diffval[m];
	    			}
	    			else {
	    				sumendleft-=diffval[m];
	    				sumendright+=diffval[m];
	    			}
	    			diffval[p]=diff;
	    			if(diff>0) sumstartright+=diff;
	    			else sumendright-=diff;
	    			lastcov=icov;
					if(diffval[m]<-1/ERROR_PERC) { // ends here
	    				if(sumendleft*ERROR_PERC>sumendright) {
	    					bool addend=true;
	    					int istart=refstart+i-CHI_THR-1;
	    					if(lend.Count()) {
	    						int j=lend.Count()-1;
	    						if(abs(istart-lend[j].predno)<CHI_WIN) {
	    							addend=false;
	    							if(lend[j].cov>0 && -diffval[m]>=lend[j].cov) {
	    								lend[j].predno=istart;
	    								lend[j].cov=-diffval[m];
	    							}
	    						}
	    					}
	    					if(addend) {
	    						CPred p(refstart+i-CHI_THR-1,-diffval[m]);
	    						lend.Add(p);
	    					}
	    				}
	    			}
	    		}
	    		//for(int i=0;i<lstart.Count();i++) fprintf(stderr,"s=%d start at coord %d with cov=%f\n",s,lstart[i].predno,lstart[i].cov);
	    		//for(int i=0;i<lend.Count();i++) fprintf(stderr,"s=%d end at coord %d with cov=%f\n",s,lend[i].predno,lend[i].cov);

	    	}
	    }

	    // see if I need to adjust the start to ignore little hanging pieces that make no sense
	    if(!end) {
	    	while(nje<njunctions && ejunction[nje]->strand+1!=2*s) nje++; // skip all junctions that are not on the same strand
	    	if(!mergeMode && (nje<njunctions && ejunction[nje]->end - currentstart < junctionsupport) &&
	    			(fs<0 || (uint)lstart[fs].predno>=ejunction[nje]->end) &&  // I do not want to miss any hard starts/ends
	    			(fe<0 || (uint)lend[fe].predno>=ejunction[nje]->end)) { // there is a junction ending soon here
	    		float covleft=get_cov(1,currentstart-refstart,ejunction[nje]->end-1-refstart,bpcov);
	    		float covright=get_cov(1,ejunction[nje]->end-refstart,2*ejunction[nje]->end - currentstart-1-refstart,bpcov);
	    		if(covleft<covright*(1-ERROR_PERC)) { // adjust start here if needed
	    			currentstart=ejunction[nje]->end;
	    			// I have to check ending junctions here again
	    			while(nje<njunctions && ejunction[nje]->end<=currentstart) { // read all junction ends at or before the current start -> assuming there are any (at this point, smaller junction ends should not be relevant to this bundle/currentstart
	    				if(ejunction[nje]->end==currentstart && (ejunction[nje]->strand+1) == 2*s) { // junction ends at current start and is on the same strand and not deleted
	    					end=1;
	    				}
	    				nje++;
	    			}
	    		}
	    	}
	    }

	    //fprintf(stderr,"create graph 1\n");
	    CGraphnode *graphnode=create_graphnode(s,g,currentstart,endbundle,graphno,bundlenode,bundle2graph,no2gnode); // creates a $graphno graphnode  with start at bundle start, and end at bundle end
	    graphno++;

	    if(end) { // I might have nodes finishing here; but I have a junction finishing here for sure
	    	//GStr cs((int)currentstart);
	    	//GVec<int> *e=ends[cs.chars()]; // HOW CAN I HAVE MORE THAN ONE NODE FINISHING HERE???; because this keeps all nodes that are linked by junctions here
	    	GVec<int> *e=ends[currentstart];
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
	    	else { // I haven't seen nodes before that finish here (maybe due to error correction?) => link to source
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

	    bool dropcov=false; // false(0) means start of bundle or junction end (raise in coverage); true(1) means junction start (drop in coverage)
	    int nls=0; // index in longstart
	    int nle=0; // index in longend

	    do {

	    	while(nje<njunctions && (((int)ejunction[nje]->strand+1) != 2*s)) nje++; // skip junctions that don't have the same strand
	    	while(njs<njunctions && ((((int)junction[njs]->strand+1)!= 2*s) || (junction[njs]->start<currentstart))) njs++; // junctions that start before the current graphnode and I haven't seen them before are part of a different bundle


	    	int minjunction = -1; // process next junction -> either a start or an ending whichever has the first position on the genome; if they have same position then process ending first
	    	if((nje<njunctions && (ejunction[nje]->end<=endbundle)) || (njs<njunctions && (junction[njs]->start<=endbundle))) {
	    		if(njs<njunctions && (junction[njs]->start<=endbundle) && junction[njs]->end>bundle_end) njs++;
	    		else {
	    			if(nje<njunctions) { // there are still junctions endings
	    				if(njs<njunctions) { // there are still junctions starting
	    					minjunction = junction[njs]->start >= ejunction[nje]->end ? 1 : 0; // one of them is clearly before the endbundle from the initial if
	    				}
	    				else minjunction = 1;
	    			}
	    			else minjunction = 0;
	    		}
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

	    					uint gstart=guideedge[nge].val;
	    					uint gend=junction[njs]->start;
	    					bool sourceguide=false;
	    					if(guideedge[nge].val<guideedge[nge].endval) sourceguide=true;
	    					nge++;
	    					if(sourceguide) { if(guideedge[nge-1].endval>endbundle) continue;}
	    					else if(guideedge[nge-1].endval<currentstart) continue;

	    					while(nge<guideedge.Count() && guideedge[nge].strand!=s) nge++;
	    					if(nge<guideedge.Count() && guideedge[nge].val<junction[njs]->start) gend=guideedge[nge].val;

	    					// I need to check there is no other trimming needed due to drops from longreads
	    					if(longreads && (lstart.Count() || lend.Count())) graphnode=longtrim(s,g,refstart,gstart,nls,nle,dropcov,!sourceguide,lstart,lend,
	    							graphnode,source,sink,futuretr,graphno,bpcov,bundlenode,bundle2graph,no2gnode,edgeno);

	    					if(sourceguide)	{
	    						graphnode=source2guide(s,g,refstart,gstart,gend,graphnode,source,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);
	    						dropcov=false;
	    					}
	    					else {
	    						graphnode=guide2sink(s,g,refstart,gstart,gend,graphnode,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);
	    						dropcov=true;
	    					}

	    				}
	    			}
	    		}
	    		//if(trim && !processguide && !mergeMode) graphnode=trimnode(s,g,refstart,junction[njs]->start,graphnode,source,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);// do something to find intermediate nodes; alternatively, I could only do this for end nodes
	    		else if(longreads && (lstart.Count() || lend.Count())) graphnode=longtrim(s,g,refstart,junction[njs]->start,nls,nle,dropcov,true,lstart,lend,
							graphnode,source,sink,futuretr,graphno,bpcov,bundlenode,bundle2graph,no2gnode,edgeno);
	    		if(trim && !longreads && !mergeMode) graphnode=trimnode_all(s,g,refstart,junction[njs]->start,graphnode,source,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);// do something to find intermediate nodes; alternatively, I could only do this for end nodes


	    		dropcov=true;
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
	    					//GStr je((int)junction[njs]->end);
	    					//GVec<int> *e=ends[je.chars()];
	    					GVec<int> *e=ends[junction[njs]->end];
	    					if(!e) {
	    						e = new GVec<int>();
	    						//ends.Add(je.chars(),e);
	    						ends.Add(junction[njs]->end, e);
	    					}
	    					e->Add(graphnode->nodeid);
	    				}
	    			}
	    			njs++;
	    		}

	    		if(pos<endbundle) { // there is still place for another node in this bundle (I might put a limit of length here for the graphnode -> because otherwise one can assume this is just a pre-mRNA fragment)

	    			// see if I should just skip node
	    			if(endbundle-pos<junctionsupport) {
	    				while(njs<njunctions && junction[njs]->strand+1 != 2*s) njs++;
	    				if(!mergeMode && (njs>=njunctions || junction[njs]->start > endbundle) && (nje>=njunctions || ejunction[nje]->end > endbundle)) { // there are no more junctions starting within this bundle
	    					float covleft=get_cov(1,2*pos-endbundle+1-refstart,pos-refstart,bpcov);
	    					float covright=get_cov(1,pos+1-refstart,endbundle-refstart,bpcov);
	    					if(covright<covleft*(1-ERROR_PERC)) { // adjust start here if needed
	    						completed=true;
	    					}
	    				}
	    			}

	    			if(!completed) {
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

	    						uint start=guideedge[nge].val;
	    						uint end=pos-1;
	    						bool sourceguide=false;
	    						if(guideedge[nge].val<guideedge[nge].endval) sourceguide=true;
	    						nge++;
	    						if(sourceguide) { if(guideedge[nge-1].endval>endbundle) continue;}
	    						else if(guideedge[nge-1].endval<currentstart) continue;

	    						while(nge<guideedge.Count() && guideedge[nge].strand!=s) nge++;
	    						if(nge<guideedge.Count() && guideedge[nge].val<pos-1) end=guideedge[nge].val;

	    						// I need to check there is no other trimming needed due to drops from longreads
	    						if(longreads && (lstart.Count() || lend.Count())) graphnode=longtrim(s,g,refstart,start,nls,nle,dropcov,!sourceguide,lstart,lend,
	    								graphnode,source,sink,futuretr,graphno,bpcov,bundlenode,bundle2graph,no2gnode,edgeno);

	    						if(sourceguide)	graphnode=source2guide(s,g,refstart,start,end,graphnode,source,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);
	    						else graphnode=guide2sink(s,g,refstart,start,end,graphnode,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);

	    					}
	    				}
	    			}
	    			//if(trim && !processguide && !mergeMode) graphnode=trimnode(s,g,refstart,pos-1,graphnode,source,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);// do something to find intermediate nodes; alternatively, I could only do this for end nodes
	    			else if(longreads && (lstart.Count() || lend.Count())) graphnode=longtrim(s,g,refstart,pos-1,nls,nle,dropcov,false,lstart,lend,
	    						graphnode,source,sink,futuretr,graphno,bpcov,bundlenode,bundle2graph,no2gnode,edgeno);
	    			if(trim && !longreads && !mergeMode) graphnode=trimnode_all(s,g,refstart,pos-1,graphnode,source,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);// do something to find intermediate nodes; alternatively, I could only do this for end nodes


	    			graphnode->end=pos-1; // set end of current graphnode here
	    			dropcov=false;
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

	    		//GStr spos((int)pos);
	    		//GVec<int> *e=ends[spos.chars()]; // WHY DOESN'T THIS REPEAT THE SAME THING IN CASE THE START HASN'T BEEN ADJUSTED? because nje is bigger now than the ones that end at the currentstart
	    		GVec<int> *e=ends[pos];
	    		if(e) for(int i=0;i<e->Count();i++) {
	    			CGraphnode *node=no2gnode[s][g][e->Get(i)];
	    			node->child.Add(graphnode->nodeid);  // this node is the child of previous node
	    			graphnode->parent.Add(node->nodeid); // this node has as parent the previous node
	    			// COUNT EDGE HERE
	    			edgeno++;
	    			//fprintf(stderr,"6 Edge %d-%d, edgeno=%d\n",node->nodeid,graphnode->nodeid,edgeno);
	    		}
	    	}

	    } while((nje<njunctions && (ejunction[nje]->end<=endbundle)) || (njs<njunctions && (junction[njs]->start<=endbundle)));


	    if(!completed) { // I did not finish node --> this will be an ending node

	    	// add guide starts/ends first
	    	if(processguide) {
	    		while(nge<guideedge.Count() && guideedge[nge].val<=graphnode->start) nge++;
	    		if(nge<guideedge.Count()) {

	    			while(true) {

	    				while(nge<guideedge.Count() && guideedge[nge].strand!=s) nge++;

	    				if(nge>=guideedge.Count() || guideedge[nge].val>=endbundle) break;

	    				uint start=guideedge[nge].val;
	    				uint end=endbundle;
	    				bool sourceguide=false;
	    				if(guideedge[nge].val<guideedge[nge].endval) sourceguide=true;
	    				nge++;
	    				if(sourceguide) { if(guideedge[nge-1].endval>endbundle) continue;}
	    				else if(guideedge[nge-1].endval<currentstart) continue;

	    				while(nge<guideedge.Count() && guideedge[nge].strand!=s) nge++;
	    				if(nge<guideedge.Count() && guideedge[nge].val<endbundle) end=guideedge[nge].val;

	    				// I need to check there is no other trimming needed due to drops from longreads
	    				if(longreads && (lstart.Count() || lend.Count())) graphnode=longtrim(s,g,refstart,start,nls,nle,dropcov,!sourceguide,lstart,lend,
	    						graphnode,source,sink,futuretr,graphno,bpcov,bundlenode,bundle2graph,no2gnode,edgeno);

	    				if(sourceguide)	graphnode=source2guide(s,g,refstart,start,end,graphnode,source,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);
	    				else graphnode=guide2sink(s,g,refstart,start,end,graphnode,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);

	    			}
	    		}
	    	}
	    	// if(trim && !processguide && !mergeMode) graphnode=trimnode(s,g,refstart,endbundle,graphnode,source,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno); // do something to find intermediate nodes; alternatively, I could only do this for end nodes
	    	else if(longreads && (lstart.Count() || lend.Count())) graphnode=longtrim(s,g,refstart,endbundle,nls,nle,dropcov,true,lstart,lend,
	    				graphnode,source,sink,futuretr,graphno,bpcov,bundlenode,bundle2graph,no2gnode,edgeno);
	    	if(trim && !longreads && !mergeMode) graphnode=trimnode_all(s,g,refstart,endbundle,graphnode,source,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno); // do something to find intermediate nodes; alternatively, I could only do this for end nodes

	    	graphnode->end=endbundle;
	    	// COUNT EDGE HERE (this is an edge to sink)
	    	edgeno++;
	    	//fprintf(stderr,"7 Edge to sink from %d, edgeno=%d\n",graphnode->nodeid,edgeno);
	    }

	    bundlenode=bundlenode->nextnode; // advance to next bundle
	} // end while(bundlenode!=NULL)

	//fprintf(stderr,"graphno=%d\n",graphno);

	// add source/sink links for very high coverage drops
	if(!mergeMode) for(int i=1;i<graphno;i++) {
		float icov=0;
		if(i>1 && no2gnode[s][g][i]->parent[0]) { // node i has parents, and does not come from source => might need to be linked to source
			/*icov=(get_cov(1,no2gnode[s][g][i]->start-refstart,no2gnode[s][g][i]->end-refstart,bpcov)-
					get_cov(2-2*s,no2gnode[s][g][i]->start-refstart,no2gnode[s][g][i]->end-refstart,bpcov))/no2gnode[s][g][i]->len();*/
			icov=get_cov_sign(2*s,no2gnode[s][g][i]->start-refstart,no2gnode[s][g][i]->end-refstart,bpcov)/no2gnode[s][g][i]->len();
			float parcov=0; // all parents' coverage
			bool addsource=true;
			for(int j=0;j<no2gnode[s][g][i]->parent.Count();j++) {
				if(!no2gnode[s][g][i]->parent[j]) { // parent is source
					addsource=false;break;
				}
				int p=no2gnode[s][g][i]->parent[j];
				/*parcov+=(get_cov(1,no2gnode[s][g][p]->start-refstart,no2gnode[s][g][p]->end-refstart,bpcov)-
						get_cov(2-2*s,no2gnode[s][g][p]->start-refstart,no2gnode[s][g][p]->end-refstart,bpcov))/no2gnode[s][g][p]->len();*/
				parcov+=get_cov_sign(2*s,no2gnode[s][g][p]->start-refstart,no2gnode[s][g][p]->end-refstart,bpcov)/no2gnode[s][g][p]->len();
			}
			if(addsource && parcov<icov*ERROR_PERC*DROP) {
				no2gnode[s][g][i]->parent.Insert(0,source->nodeid);
				source->child.Add(i);
				futuretr.cAdd(0.0);
				float tmp=i;
				futuretr.Add(tmp);
				tmp=(icov-parcov)/DROP;
				futuretr.Add(tmp);
				//fprintf(stderr,"Add link to source from node %d with abundance %f\n",i,tmp);
				edgeno++;
			}
		}

		if(i<graphno-1) {
			if(no2gnode[s][g][i]->child.Count()) { // might need to be linked to sink
				bool addsink=true;
				for(int j=0;j<sink->parent.Count();j++) {
					if(sink->parent[j]==i) {
						addsink=false;
						break;
					}
				}
				if(addsink) {
					/*if(!icov) icov=(get_cov(1,no2gnode[s][g][i]->start-refstart,no2gnode[s][g][i]->end-refstart,bpcov)-
							get_cov(2-2*s,no2gnode[s][g][i]->start-refstart,no2gnode[s][g][i]->end-refstart,bpcov))/no2gnode[s][g][i]->len();*/
					if(!icov) icov=get_cov_sign(2*s,no2gnode[s][g][i]->start-refstart,no2gnode[s][g][i]->end-refstart,bpcov)/no2gnode[s][g][i]->len();
					float chcov=0; // all children coverages
					for(int j=0;j<no2gnode[s][g][i]->child.Count();j++) {
						int c=no2gnode[s][g][i]->child[j];
						/*chcov+=(get_cov(1,no2gnode[s][g][c]->start-refstart,no2gnode[s][g][c]->end-refstart,bpcov)-
								get_cov(2-2*s,no2gnode[s][g][c]->start-refstart,no2gnode[s][g][c]->end-refstart,bpcov))/no2gnode[s][g][c]->len();*/
						chcov+=get_cov_sign(2*s,no2gnode[s][g][c]->start-refstart,no2gnode[s][g][c]->end-refstart,bpcov)/no2gnode[s][g][c]->len();
					}
					if(chcov<icov*ERROR_PERC*DROP) {
						sink->parent.Add(i); // prevnode is the parent of sink
						float tmp=i;
						futuretr.Add(tmp);
						futuretr.cAdd(-1.0);
						tmp=(icov-chcov)/DROP;
						futuretr.Add(tmp);
						edgeno++;
					}
				}
			}
			else edgeno++; // node has no children so it might get linked to sink in traverse_dfs function
		}
	}

	// here I know graphno => I can see if it's too big
	if(!mergeMode && graphno>allowed_nodes) { // TODO: define allowed_nodes as a default in stringtie.cpp that varies with the memory
		graphno=prune_graph_nodes(graphno,s,g,bundle2graph,bnode.Count(),no2gnode,junction,edgeno,futuretr,sink);
	}

	sink->nodeid=graphno;
	no2gnode[s][g].Add(sink);
	graphno++;

	if(mergeMode) { // I might have a bunch of sink's parents that are not linked to sink
		for(int i=0;i<sink->parent.Count();i++) {
			CGraphnode *node=no2gnode[s][g][sink->parent[i]];
			node->child.Add(sink->nodeid);
		}
	}

	//fprintf(stderr,"This graph has %d nodes and %d edges and starts at lastpos=%d\n",graphno,edgeno,graphno);
	lastgpos=graphno; // nodes are from 0 to graphno-1, so the first "available" position in GBitVec is graphno

	// now I can create the future transfrags because I know graphno
	for(int i=0;i<futuretr.Count();i+=3) {
		// add links between node and sink
		int n1=int(futuretr[i]);
		if(n1 >=0 ) {
			int n2=int(futuretr[i+1]);
			GBitVec trpat(graphno+edgeno);
			trpat[n1]=1;
			GVec<int> nodes;
			if(n2<0) {
				CGraphnode *node=no2gnode[s][g][n1];

				int n3=n1+1;
				if(n3<sink->nodeid) {
					CGraphnode *prevnode=node;
					CGraphnode *nextnode=no2gnode[s][g][n3]; // this should be okay
					uint dist=nextnode->len()-1;
					bool keeptr=true;
					while(nextnode->start==prevnode->end+1 && dist<longintronanchor) {
						for(int c=0;c<nextnode->child.Count();c++)
							if(nextnode->child[c]==sink->nodeid) {
								// do not continue with the futuretr
								keeptr=false;
								break;
							}
						n3++;
						if(n3==sink->nodeid) {
							// no link to sink found -> can continue with futuretr
							break;
						}
						prevnode=nextnode;
						nextnode=no2gnode[s][g][n3];
						dist+=nextnode->len();
					}
					if(!keeptr) continue; // found link to sink; no need to keep futuretr
				}

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
				fprintf(stderr,"Add future transfrag[%d][%d]= %d with %d nodes n1=%d n2=%d graphno=%d, abundance=%f and pattern",s,g,transfrag[s][g].Count(),tr->nodes.Count(),n1,n2,graphno,futuretr[i+2]);
				//printBitVec(trpat);
				fprintf(stderr,"\n");
			}
			*/

			/*if(mixedMode) {
				tr->abundance*=2;
			}*/

			transfrag[s][g].Add(tr);
			if(mixedMode) {
				CTransfrag *longtr=new CTransfrag(nodes,trpat,futuretr[i+2]);
				longtr->longread=true;
				transfrag[s][g].Add(longtr);
			}
			else
			if(longreads)// || mixedMode)
				transfrag[s][g].Last()->longread=true;
		}
	}

	// finished reading bundle -> now create the parents' and children's patterns
	GVec<bool> visit;
	visit.Resize(graphno);
	GBitVec parents(graphno+edgeno);

	//fprintf(stderr,"traverse graph[%d][%d] now with %d nodes, %d edges and lastgpos=%d....\n",s,g,graphno,edgeno,lastgpos);//edgeno=0;
	traverse_dfs(s,g,source,sink,parents,graphno,visit,no2gnode,transfrag,edgeno,gpos,lastgpos);
	//fprintf(stderr,"done traversing with edgeno=%d lastgpos=%d\n",edgeno,lastgpos);

	/*
	{ //DEBUG ONLY
		fprintf(stderr,"after traverse:\n");
		for(int i=0;i<graphno;i++) {
			fprintf(stderr,"Node %d with parents:",i);
			for(int p=0;p<no2gnode[s][g][i]->parent.Count();p++) fprintf(stderr," %d",no2gnode[s][g][i]->parent[p]);
			fprintf(stderr," and children:");
			for(int c=0;c<no2gnode[s][g][i]->child.Count();c++) fprintf(stderr," %d",no2gnode[s][g][i]->child[c]);
			fprintf(stderr,"\n");
		}
	}
	*/

	// delete variables created here, like e.g. ends; do I need to delete the GVec<int> elements created too?
	ends.Clear();

	return(graphno);

}


void get_read_pattern(int s, float readcov,GVec<int> &rgno, float rprop,GVec<int> *rnode,
		GList<CReadAln>& readlist,int n,GVec<int> *readgroup,GVec<int>& merge,GVec<int> *group2bundle,
		GVec<CGraphinfo> **bundle2graph,GPVec<CGraphnode> **no2gnode) {

	int lastgnode=-1;
	int lastngraph=-1;
	int ncoord=readlist[n]->segs.Count();
	int k=0; // need to keep track of coordinates already added to coverages of graphnodes
	int kmer=KMER-1; //f1

	GIntHash<bool> hashnode;

	for(int i=0;i<readgroup[n].Count();i++) { // how can a read be associated to multiple groups? ---> If it is spliced
		int gr=readgroup[n][i];
		while(merge[gr]!=gr) gr=merge[gr];
		readgroup[n][i]=gr;
		int bnode=group2bundle[2*s][gr]; // group was associated to bundle
		if(bnode>-1 && bundle2graph[s][bnode].Count()) { // group has a bundle node associated with it and bundle was processed
			if(!hashnode[bnode]) {
				hashnode.Add(bnode,true);
				int nbnode=bundle2graph[s][bnode].Count(); // number of nodes in bundle
				int j=0;

				//fprintf(stderr,"bundle2graph[%d] for read %d:",bnode,n);

				while(j<nbnode && k<ncoord) {
					int ngraph=bundle2graph[s][bnode][j].ngraph;
					int gnode=bundle2graph[s][bnode][j].nodeno;
					CGraphnode *node=no2gnode[s][ngraph].Get(gnode);
					// compute graphnode coverage by read here (assuming they come in order on the genomic line)
					bool intersect=false;
					while(k<ncoord) { // found how much of readlist[n] intersects node represented by bnode; here I assume that read is always included in bundle2graph

						if(readlist[n]->segs[k].end<node->start) k++;
						else {

							int bp = readlist[n]->segs[k].overlapLen(node);
							if(bp) {
								intersect=true;

								if(!readlist[n]->unitig)
										node->cov+=rprop*bp*readcov;
								if(readlist[n]->segs[k].end<=node->end) k++;
								else break;
							}
							else break;
						}
					}

					if(intersect) { // read intersects gnode
						int g=-1; // position of graph in list
						if(rgno.Count()) {
							if(ngraph!=rgno.Last()) {
								for(int r=0;r<rgno.Count()-1;r++) { // count should be small
									if(ngraph==rgno[r]) {
										g=r;
										break;
									}
								}
							}
							else g=rgno.Count()-1;
						}
						if(g<0) { // add ngraph to list of graphs if not seen before
							rgno.Add(ngraph);
							g=rgno.Count()-1;
						}

						if(lastngraph!=ngraph || lastgnode!=gnode) { // this is a new graph, or a new gnode

							bool trimu=false;
							if(readlist[n]->unitig) { //f1
								if(rnode[g].Count()==1) { // deal with first node trimming here
									if(no2gnode[s][ngraph][gnode]->parent.Count()>1) { // only if this node has more than one parent I might need to trim it
										if((int)(no2gnode[s][ngraph][rnode[g][0]]->end - readlist[n]->start) < kmer) { // I have to trim this part of the super-read because it might be wrong
											rnode[g].Pop();
										}
									}
								}
								if(rnode[g].Count() && no2gnode[s][ngraph][rnode[g].Last()]->child.Count()>1) {
									if((int)(readlist[n]->end-no2gnode[s][ngraph][gnode]->start) < kmer) { // super-read extends only a tiny bit into the node -> might be wrong
										trimu=true;
									}
								}
							}


							if(!trimu) rnode[g].Add(gnode);
						}

						lastgnode=gnode;
						lastngraph=ngraph;
					} // end if(intersect)

					j++;
				} // end while(j<nbnode && k<ncoord)
			} // end hashnode
		} // end if(bnode>-1 && bundle2graph[s][bnode].Count())
	} // end for
}


void get_read_pattern(float readcov,GBitVec& pattern0,GBitVec& pattern1,int *rgno, float *rprop,GVec<int> *rnode,
		GList<CReadAln>& readlist,int n,GVec<int> *readgroup,GVec<int>& merge,GVec<int> *group2bundle,
		GVec<CGraphinfo> **bundle2graph,GVec<int> *graphno,GVec<int> *edgeno,GIntHash<int> **gpos,
		GPVec<CGraphnode> **no2gnode) {


	int lastgnode[2]={-1,-1}; // lastgnode[0] is for - strand; [1] is for + strand -> I need these in order to add the edges to the read pattern; check this: if it's not correct than storage was wrong!
	int ncoord=readlist[n]->segs.Count();

	int k[2]={0,0}; // need to keep track of coordinates already added to coverages of graphnodes
    bool valid[2]={true,true};
    for(int s=0;s<2;s++) if(!rprop[s]) valid[s]=false;

    for(int i=0;i<readgroup[n].Count();i++) // how can a read be associated to multiple groups? ---> If it is spliced
    	if(valid[0] || valid[1]) { // there are still stranded bundles associated with the read
    		int gr=readgroup[n][i];
    		while(merge[gr]!=gr) gr=merge[gr];
    		for(int s=0;s<2;s++)
    			if(valid[s]) {
    				int bnode=group2bundle[2*s][gr]; // group was associated to bundle
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
    								if(!readlist[n]->unitig)
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
	//fprintf(stderr,"There are %d transfrags\n",transfrag.Count());
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
			if(pos && pattern[*pos]) { // there is an edge between node[n-1] and node[n]
				tree=tree->nextpat[gno-1-node[n-1]+node[n]-node[n-1]-1];
			}
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


CTransfrag *update_abundance(int s,int g,int gno,GIntHash<int>&gpos,GBitVec& pattern,float abundance,GVec<int>& node,
		GPVec<CTransfrag> **transfrag,CTreePat ***tr2no, GPVec<CGraphnode>& no2gnode,uint rstart,uint rend,bool is_sr=false,bool is_lr=false){

	/*
	{ // DEBUG ONLY
		fprintf(stderr,"Update transfrag[%d][%d] longread=%d:",s,g,is_lr);
		for(int i=0;i<node.Count();i++) fprintf(stderr," %d",node[i]);
		fprintf(stderr," with abundance=%f in graph[%d][%d]\n",abundance,s,g);
	}
	*/

	if(is_lr) { // see if I need to adjust nodes in pattern
		if(node.Count()>1) { // correct longreads that extend inside introns or past nodes that link back to source or sink

			bool is_source=(rstart>=no2gnode[node[0]]->start); // read start could be inside of the node (node could be extended due to other reads)
			bool is_sink=(rend<=no2gnode[node.Last()]->end); // read end could be inside of the node (node could be extended due to other reads
			if(no2gnode[node[0]]->hardstart) is_source=false; // do not do any trimming if this is an ending node
			if(no2gnode[node.Last()]->hardend) is_sink=false; // do not do any trimming if this is an ending node

			if(is_sink) { // check last node
				int i=node.Count()-1;
				while(i>0 && no2gnode[node[i-1]]->end+1==no2gnode[node[i]]->start) {
					uint dist=rend-no2gnode[node[i-1]]->end;
					if(dist<longintronanchor) { // I can potentially eliminate node[i] from transfrag due to proximity
						if(no2gnode[node[i-1]]->hardend || no2gnode[node[i-1]]->cov*DROP*no2gnode[node[i]]->len()>no2gnode[node[i]]->cov*no2gnode[node[i-1]]->len()) { // there is a significant drop
							int j=no2gnode[node[i-1]]->child.Count()-1;
							bool trim=false;
							while(j>=0) {
								if(no2gnode[node[i-1]]->child[j]!=node[i]) { // there is a different child and a significant drop within short distance
									// adjust pattern
									int k=i;
									while(k<node.Count()) {
										pattern[node[k]]=0;
										int *pos=gpos[edge(node[k-1],node[k],gno)];
										if(pos) pattern[*pos]=0;
										k++;
									}
									node.setCount(i); // adust number of nodes
									trim=true;
									break;
								}
								j--;
							}
							if(trim) break;
						}
					}
					else if(dist<CHI_WIN) { // check for sink connection
						if(no2gnode[node[i-1]]->hardend || no2gnode[node[i-1]]->cov*DROP*no2gnode[node[i]]->len()>no2gnode[node[i]]->cov*no2gnode[node[i-1]]->len()) { // there is a significant drop
							int j=no2gnode[node[i-1]]->child.Count()-1;
							bool trim=false;
							while(j>=0) {
								if(no2gnode[node[i-1]]->child[j]==gno-1) { // node[i-1] can connect to source
									// adjust pattern
									int k=i;
									while(k<node.Count()) {
										pattern[node[k]]=0;
										int *pos=gpos[edge(node[k-1],node[k],gno)];
										if(pos) pattern[*pos]=0;
										k++;
									}
									node.setCount(i); // adust number of nodes
									trim=true;
									break;
								}
								j--;
							}
							if(trim) break;
						}
					}
					else break;
					i--;
				}
			}

			if(is_source) { // check first node
				int i=0;
				while(i+1<node.Count() && no2gnode[node[i]]->end+1==no2gnode[node[i+1]]->start) {
					uint dist=no2gnode[node[i+1]]->start-rstart;
					if(dist<longintronanchor) { // I can potentially eliminate node[i] from transfrag due to proximity
						if(no2gnode[node[i+1]]->hardstart || no2gnode[node[i+1]]->cov*DROP*no2gnode[node[i]]->len()>no2gnode[node[i]]->cov*no2gnode[node[i+1]]->len()) { // there is a significant drop
							int j=0;
							bool trim=false;
							while(j<no2gnode[node[i+1]]->parent.Count()) {
								if(no2gnode[node[i+1]]->parent[j]!=node[i]) { // there is a different child and a significant drop within short distance
									// adjust pattern
									int k=i;
									while(k>=0) {
										pattern[node[k]]=0;
										int *pos=gpos[edge(node[k],node[k+1],gno)];
										if(pos) pattern[*pos]=0;
										k--;
									}
									node.Shift(i+1); // delete first i+1 nodes in transfrag
									trim=true;
									break;
								}
								j++;
							}
							if(trim) break;
						}
					}
					else if(dist<CHI_WIN) { // check for sink connection
						if(no2gnode[node[i+1]]->hardstart || no2gnode[node[i+1]]->cov*DROP*no2gnode[node[i]]->len()>no2gnode[node[i]]->cov*no2gnode[node[i+1]]->len()) { // there is a significant drop
							int j=0;
							bool trim=false;
							while(j<no2gnode[node[i+1]]->parent.Count()) {
								if(!no2gnode[node[i+1]]->parent[j]) { // there is a link to source from node[i+1]
									// adjust pattern
									int k=i;
									while(k>=0) {
										pattern[node[k]]=0;
										int *pos=gpos[edge(node[k],node[k+1],gno)];
										if(pos) pattern[*pos]=0;
										k--;
									}
									node.Shift(i+1); // delete first i+1 nodes in transfrag
									trim=true;
									break;
								}
								j++;
							}
							if(trim) break;
						}
					}
					else break;
					i++;
				}
			}
		}
	}
	else if(!mergeMode && node.Count()==1) return(NULL); // do the one node transfrag make any difference? CHECK IF YOU NEED TO KEEP THIS ONE

	CTransfrag *t=findtrf_in_treepat(gno,gpos,node,pattern,tr2no[s][g]);
	if(!t) { // t is NULL
	  t=new CTransfrag(node,pattern,0);

	  /*
		{ // DEBUG ONLY
			fprintf(stderr,"Add update transfrag[%d][%d]=%d and pattern",s,g,transfrag[s][g].Count());
			fprintf(stderr," (check nodes:");
			for(int i=0;i<node.Count();i++) fprintf(stderr," %d",node[i]);
			fprintf(stderr,")");
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
	//else if(!is_lr && t->longread) t->longread=false;
	if(is_lr) t->longread=true;
	else {
		//fprintf(stderr,"Set short read\n");
		t->shortread=true;
	}

	if(is_sr) t->srabund+=abundance;
	else t->abundance+=abundance;

	if(no2gnode[node[0]]->start<=rstart) if(!t->longstart || rstart<t->longstart) t->longstart=rstart; // value of longstart inside node; the other extensions could come from spliced nodes
	if(no2gnode[node.Last()]->end>=rend) if(!t->longend|| rend>t->longend) t->longend=rend; // value of longend inside node; the other extensions could come from spliced nodes

	return(t);
}

void get_fragment_pattern(GList<CReadAln>& readlist,int n, int np,float readcov,GVec<int> *readgroup,GVec<int>& merge, GVec<int> *group2bundle,
		GVec<CGraphinfo> **bundle2graph,GVec<int> *graphno,GVec<int> *edgeno, GIntHash<int> **gpos,GPVec<CGraphnode> **no2gnode,
		GPVec<CTransfrag> **transfrag,CTreePat ***tr2no,GPVec<CGroup> &group) {

	/*fprintf(stderr,"get fragment for read[%d]:%d-%d-%d-%d-%f with pair[%d] and longread=%d and exons: ",n,readlist[n]->start,readlist[n]->end,int(readlist[n]->strand),readlist[n]->nh,readlist[n]->read_count,np,readlist[n]->longread);
	for(int i=0;i<readlist[n]->segs.Count();i++) fprintf(stderr," %d-%d",readlist[n]->segs[i].start,readlist[n]->segs[i].end);
	if(np>-1) for(int i=0;i<readlist[np]->segs.Count();i++) fprintf(stderr," %d-%d",readlist[np]->segs[i].start,readlist[np]->segs[i].end);
	fprintf(stderr,"\n");*/
	uint rstart=readlist[n]->start; // this only works for unpaired long reads -> I will have to take into account the pair if I want to do it for all reads
	uint rend=readlist[n]->end;
	if(np>-1 && readlist[np]->end>rend) rend=readlist[np]->end;

	float rprop[2]={0,0}; // by default read does not belong to any strand
	// compute proportions of unstranded read associated to strands
	if(readlist[n]->nh && !readlist[n]->strand && np>-1 && readlist[np]->nh && !readlist[np]->strand) { // both reads are unstranded
		int ngroup=0; // due to errors read might not belong to any groups

		for(int i=0;i<readgroup[n].Count();i++) {
			ngroup++;
			int gr=readgroup[n][i]; // read is unstranded => it should belong to one group only -> not necessarily see above
			while(merge[gr]!=gr) gr=merge[gr];
			rprop[0]+=group[gr]->neg_prop;
		}
		for(int i=0;i<readgroup[np].Count();i++) {
			ngroup++;
			int gr=readgroup[np][i]; // read is unstranded => it should belong to one group only
			while(merge[gr]!=gr) gr=merge[gr];
			rprop[0]+=group[gr]->neg_prop;
		}

		if(ngroup) {
			rprop[0]/=ngroup;
			rprop[1]=1-rprop[0];
		}
	}
	else {
		if(readlist[n]->nh) {
			if(!readlist[n]->strand) { // the paired read is not present otherwise it would have the same strand from add_read_to_group

				int ngroup=0; // due to errors read might not belong to any groups
				for(int i=0;i<readgroup[n].Count();i++) {
					ngroup++;
					int gr=readgroup[n][i]; // read is unstranded => it should belong to one group only -> not necessarily see above
					while(merge[gr]!=gr) gr=merge[gr];
					rprop[0]+=group[gr]->neg_prop;
				}
				if(ngroup) {
					rprop[0]/=ngroup;
					rprop[1]=1-rprop[0];
				}

			}
			else {
				if(readlist[n]->strand==-1) rprop[0]=1;
				else rprop[1]=1;
			}
		}
		else if(np>-1 && readlist[np]->nh) { // readlist[n] is deleted
			if(!readlist[np]->strand) { // the paired read is not present otherwise it would have the same strand from add_read_to_group

				int ngroup=0; // due to errors read might not belong to any groups
				for(int i=0;i<readgroup[np].Count();i++) {
					ngroup++;
					int gr=readgroup[np][i]; // read is unstranded => it should belong to one group only -> not necessarily see above
					while(merge[gr]!=gr) gr=merge[gr];
					rprop[0]+=group[gr]->neg_prop;
				}
				if(ngroup) {
					rprop[0]/=ngroup;
					rprop[1]=1-rprop[0];
				}
			}
			else {
				if(readlist[np]->strand==-1) rprop[0]=1;
				else rprop[1]=1;
			}
		}
	}

	for(int s=0;s<2;s++) if(rprop[s]){

		GVec<int> rgno;
		GVec<int> *rnode=new GVec<int>[readgroup[n].Count()];
		if(readlist[n]->nh)
			get_read_pattern(s, readcov,rgno,rprop[s],rnode,readlist,n,readgroup,merge,group2bundle,bundle2graph,no2gnode);

		GVec<int> pgno;
		GVec<int> *pnode=NULL;
		// get pair pattern if pair exists and it hasn't been deleted
		if(np>-1 && readlist[np]->nh) {
			pnode=new GVec<int>[readgroup[n].Count()];
			get_read_pattern(s,readcov,pgno,rprop[s],pnode,readlist,np,readgroup,merge,group2bundle,bundle2graph,no2gnode);
		}


		GBitVec rpat;
		GBitVec ppat;
		int usedp=0;

		for(int r=0;r<rgno.Count();r++) {

			// set read pattern
			rpat.clear();
			rpat.resize(graphno[s][rgno[r]]+edgeno[s][rgno[r]]);
			int i=0;
			int clear_rnode=0;
			for(int j=0;j<rnode[r].Count();j++) {

				rpat[rnode[r][j]]=1;
				if(j) { // add edge pattern
					if(no2gnode[s][rgno[r]][rnode[r][j]]->start-1==no2gnode[s][rgno[r]][rnode[r][j-1]]->end) { // continuous node
						bool cont=true;
						if(readlist[n]->unitig) { // see if node j-1 has sink as child, or node j has source as parent -> this means trimming was involved
							CGraphnode *jnode=no2gnode[s][rgno[r]][rnode[r][j]];
							for(int p=0;p<jnode->parent.Count();p++) {
								if(!jnode->parent[p]) { cont=false; break;} // found source among parents
							}
							if(cont) {
								jnode=no2gnode[s][rgno[r]][rnode[r][j-1]];
								for(int c=0;c<jnode->child.Count();c++) {
									if(jnode->child[c]==graphno[s][rgno[r]]-1) { cont=false; break;} // found sink among children
								}
							}
							if(!cont) { // trimming was involved
								rpat[rnode[r][j]]=0;
								GVec<int> node;
								for(int o=clear_rnode;o<j;o++) {
									node.Add(rnode[r][o]);
								}
								update_abundance(s,rgno[r],graphno[s][rgno[r]],gpos[s][rgno[r]],rpat,rprop[s]*readcov,node,transfrag,tr2no,
										no2gnode[s][rgno[r]],rstart, rend,readlist[n]->unitig,readlist[n]->longread);
								rpat.reset();
								rpat[rnode[r][j]]=1; // restart the pattern
								clear_rnode=j;
							}
						}
						if(cont) {
							int *pos=gpos[s][rgno[r]][edge(rnode[r][j-1],rnode[r][j],graphno[s][rgno[r]])];
							if(pos!=NULL) rpat[*pos]=1;
						}
					}
					else { // non-continuous node
						while(i<readlist[n]->juncs.Count() && readlist[n]->juncs[i]->end<no2gnode[s][rgno[r]][rnode[r][j]]->start) i++; // end of junction should be after start of node
						if(i<readlist[n]->juncs.Count() && readlist[n]->juncs[i]->start<=no2gnode[s][rgno[r]][rnode[r][j-1]]->end) {
							if(readlist[n]->juncs[i]->strand && readlist[n]->juncs[i]->start==no2gnode[s][rgno[r]][rnode[r][j-1]]->end &&
									readlist[n]->juncs[i]->end==no2gnode[s][rgno[r]][rnode[r][j]]->start) { // valid junction
								int *pos=gpos[s][rgno[r]][edge(rnode[r][j-1],rnode[r][j],graphno[s][rgno[r]])];
								if(pos!=NULL) rpat[*pos]=1;
							}
							else { // update abundance for what I have so far
								rpat[rnode[r][j]]=0;
								GVec<int> node;
								for(int o=clear_rnode;o<j;o++) {
									node.Add(rnode[r][o]);
								}
								update_abundance(s,rgno[r],graphno[s][rgno[r]],gpos[s][rgno[r]],rpat,rprop[s]*readcov,node,transfrag,tr2no,
										no2gnode[s][rgno[r]],rstart, rend,readlist[n]->unitig,readlist[n]->longread);
								rpat.reset();
								rpat[rnode[r][j]]=1; // restart the pattern
								clear_rnode=j;
							}
						}
					}
				}
			}
			for(int j=0;j<clear_rnode;j++) rnode[r].Shift();

			int p=0;
			while(p<pgno.Count() && rgno[r]!=pgno[p]) p++;
			if(p<pgno.Count()) { // I found a node that is within the same graph -> can only happen if pair read is present;

				usedp++;

				// set pair pattern
				ppat.clear();
				ppat.resize(graphno[s][rgno[r]]+edgeno[s][rgno[r]]);
				i=0;
				for(int j=0;j<pnode[p].Count();j++) {
					ppat[pnode[p][j]]=1;
					if(j) { // add edge pattern
						if(no2gnode[s][rgno[r]][pnode[p][j]]->start-1==no2gnode[s][rgno[r]][pnode[p][j-1]]->end) { // continuous node
							int *pos=gpos[s][rgno[r]][edge(pnode[p][j-1],pnode[p][j],graphno[s][rgno[r]])];
							if(pos!=NULL) ppat[*pos]=1;
						}
						else { // non-continuous node
							while(i<readlist[np]->juncs.Count() && readlist[np]->juncs[i]->start<no2gnode[s][rgno[r]][pnode[p][j-1]]->end) i++;
							if(i<readlist[np]->juncs.Count() && readlist[np]->juncs[i]->strand && readlist[np]->juncs[i]->start==no2gnode[s][rgno[r]][pnode[p][j-1]]->end &&
									readlist[np]->juncs[i]->end==no2gnode[s][rgno[r]][pnode[p][j]]->start) {
								int *pos=gpos[s][rgno[r]][edge(pnode[p][j-1],pnode[p][j],graphno[s][rgno[r]])];
								if(pos!=NULL) ppat[*pos]=1;
							}
						}
					}
				}

				// check if there is a conflict of patterns
				CGraphnode *gnode=no2gnode[s][rgno[r]][pnode[p][0]];
				GBitVec conflictpattn=gnode->parentpat;
				conflictpattn[pnode[p][0]]=1;

				if((conflictpattn & rpat)==rpat) { // there isn't a conflict -> pair parents should contain read pattern
					i=0;
					if(pnode[p][0]==rnode[r].Last()) // read and pair share a node
						i++;
					while(i<pnode[p].Count()) { rnode[r].Add(pnode[p][i]);i++;}
					rpat=rpat|ppat;
					update_abundance(s,rgno[r],graphno[s][rgno[r]],gpos[s][rgno[r]],rpat,rprop[s]*readcov,rnode[r],transfrag,tr2no,
							no2gnode[s][rgno[r]],rstart, rend,readlist[n]->unitig,readlist[n]->longread);
				}
				else { // update both patterns separately
					update_abundance(s,rgno[r],graphno[s][rgno[r]],gpos[s][rgno[r]],rpat,rprop[s]*readcov,rnode[r],transfrag,tr2no,
							no2gnode[s][rgno[r]],rstart, rend,readlist[n]->unitig,readlist[n]->longread);
					update_abundance(s,pgno[p],graphno[s][pgno[p]],gpos[s][pgno[p]],ppat,rprop[s]*readcov,pnode[p],transfrag,tr2no,
							no2gnode[s][pgno[p]],rstart, rend,readlist[n]->unitig,readlist[n]->longread);
				}
				pgno[p]=-1;
			}
			else { // pair has no valid pattern in this graph
				update_abundance(s,rgno[r],graphno[s][rgno[r]],gpos[s][rgno[r]],rpat,rprop[s]*readcov,rnode[r],transfrag,tr2no,
						no2gnode[s][rgno[r]],rstart, rend,readlist[n]->unitig,readlist[n]->longread);
			}
		}

		if(usedp<pgno.Count()) for(int p=0;p<pgno.Count();p++)
			if(pgno[p]>-1) {

				// set pair pattern
				ppat.clear();
				ppat.resize(graphno[s][pgno[p]]+edgeno[s][pgno[p]]);
				int i=0;
				for(int j=0;j<pnode[p].Count();j++) {
					ppat[pnode[p][j]]=1;
					if(j) { // add edge pattern
						if(no2gnode[s][pgno[p]][pnode[p][j]]->start-1==no2gnode[s][pgno[p]][pnode[p][j-1]]->end) { // continuous node
							int *pos=gpos[s][pgno[p]][edge(pnode[p][j-1],pnode[p][j],graphno[s][pgno[p]])];
							if(pos!=NULL) ppat[*pos]=1;
						}
						else { // non-continuous node
							while(i<readlist[np]->juncs.Count() && readlist[np]->juncs[i]->start<no2gnode[s][pgno[p]][pnode[p][j-1]]->end) i++;
							if(i<readlist[np]->juncs.Count() && readlist[np]->juncs[i]->strand && readlist[np]->juncs[i]->start==no2gnode[s][pgno[p]][pnode[p][j-1]]->end &&
									readlist[np]->juncs[i]->end==no2gnode[s][pgno[p]][pnode[p][j]]->start) {
								int *pos=gpos[s][pgno[p]][edge(pnode[p][j-1],pnode[p][j],graphno[s][pgno[p]])];
								if(pos!=NULL) ppat[*pos]=1;
							}
						}
					}
				}
				update_abundance(s,pgno[p],graphno[s][pgno[p]],gpos[s][pgno[p]],ppat,rprop[s]*readcov,pnode[p],transfrag,tr2no,
						no2gnode[s][pgno[p]],rstart, rend,readlist[n]->unitig,readlist[n]->longread);
			}
		delete [] rnode;
		if(pnode) delete [] pnode;
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
			CTransfrag *t=update_abundance(s,rgno[s],graphno[s][rgno[s]],gpos[s][rgno[s]],rpat[s],rprop[s]*readlist[n]->read_count,rnode[s],transfrag,tr2no,no2gnode[s][rgno[s]],readlist[n]->start,readlist[n]->end);

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

bool eliminate_transfrags_under_thr(int gno,GIntHash<int>& gpos,GPVec<CTransfrag>& transfrag, CTreePat *tr2no,float threshold,GPVec<CTransfrag>& srfrag) {

	bool longret=false;

	//fprintf(stderr,"eliminate from %d transfrags\n",transfrag.Count());

	for(int t=transfrag.Count()-1;t>=0;t--)
		//if(transfrag[t]->srabund || (mixedMode && (transfrag[t]->longread || !transfrag[t]->nodes[0] || transfrag[t]->nodes.Last()==gno-1))) { // this is a super-read
			if(transfrag[t]->srabund || (mixedMode && transfrag[t]->longread)) { // this is a super-read
			srfrag.Add(transfrag[t]);
			if(mixedMode) {
				if(transfrag[t]->nodes[0] && transfrag[t]->nodes.Last()!=gno-1) // proper longread transfrag
					longret=true;
			}
		}
		else if(!transfrag[t]->guide && transfrag[t]->abundance<threshold &&
			  transfrag[t]->nodes[0] && transfrag[t]->nodes.Last()<gno-1) { // need to delete transfrag that doesn't come from source or ends at sink
			settrf_in_treepat(NULL,gno,gpos,transfrag[t]->nodes,transfrag[t]->pattern,tr2no); // this should be eliminated if I want to store transcripts from 0 node
			transfrag.Exchange(t,transfrag.Count()-1);
			transfrag.Delete(transfrag.Count()-1);
		}

	while(transfrag.Count()>max_trf_number) {
		threshold++;
		for(int t=transfrag.Count()-1;t>=0;t--)
			if(!transfrag[t]->guide && transfrag[t]->abundance<threshold &&
			    transfrag[t]->nodes[0] && transfrag[t]->nodes.Last()<gno-1) { // need to delete transfrag that doesn't come from source
				settrf_in_treepat(NULL,gno,gpos,transfrag[t]->nodes,transfrag[t]->pattern,tr2no); // this should be eliminated if I want to store transcripts from 0 node
				transfrag.Exchange(t,transfrag.Count()-1);
				transfrag.Delete(transfrag.Count()-1);
			}
	}
	return(longret);
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
		else { // found a good position to insert
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

				// I check if there are different paths to reach nextnode from node
				for(int j=0;j<node->child.Count();j++) {
					if(node->child[j]>nextnode) break; // if child of node is after nextnode -> I don't care
					if(node->child[j]==nextnode || no2gnode[node->child[j]]->childpat[nextnode]) { // I found a node on the path to nextnode
						CPath p(n,node->child[j]);
						int *pos=gpos[edge(n,node->child[j],gno)];
						if(pos) for(int f=0;f<ntrf;f++) {
							if(transfrag[node->trf[f]]->pattern[*pos]) { // this is a transfrag that goes through node and child; ATTENTION: it might not reach nextnode!! or be compatible with nextnode
								p.abundance+=transfrag[node->trf[f]]->abundance;
								totalcov+=transfrag[node->trf[f]]->abundance;
							}
						}
						if(p.abundance) {
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

float stitch_trf(CTransfrag *u,float available,int i_start, int i,GPVec<CTransfrag> &seltrfrag,GVec<int> *start,GVec<int> next) {
	if(seltrfrag[i]->abundance<available) available=seltrfrag[i]->abundance;
	int nl=seltrfrag[i]->nodes.Last();
	if(nl==u->nodes.Last()) nl=u->nodes[0];
	if(nl==seltrfrag[i_start]->nodes[0]) return(available);
	int nf=seltrfrag[i]->nodes[0];
	if(nf<seltrfrag[i_start]->nodes[0] && nl>seltrfrag[i_start]->nodes[0]) return(0);
	float maxabund=0;
	for(int j=0;j<start[nl].Count();j++) {
		if(start[nl][j]>i_start && seltrfrag[start[nl][j]]->abundance) maxabund=stitch_trf(u,available,i_start,start[nl][j],seltrfrag,start,next);
		if(maxabund) {
			next.Add(start[nl][j]);
			return(maxabund);
		}
	}
	return(0);
}


void process_srfrag(CTransfrag *u,GPVec<CTransfrag>& transfrag,GPVec<CGraphnode>& no2gnode,int gno, GIntHash<int> &gpos) { // v7

	//fprintf(stderr,"Process super-read with abundance=%f and srabundance=%f\n",u->abundance,u->srabund);

        float available=u->srabund*DROP*ERROR_PERC-u->abundance; // abundance of sr is higher than spliced read alignment
	if(available>epsilon) { // I have some abundance to explain
		GPVec<CTransfrag> seltrfrag(false);
		GVec<int> equalinc;
		GVec<int> equalext;
		int n=u->nodes[0];
		for(int j=0;j<no2gnode[n]->trf.Count();j++) if(transfrag[no2gnode[n]->trf[j]]!=u) {
			int t=no2gnode[n]->trf[j];
			if(transfrag[t]->nodes[0]<=n && transfrag[t]->nodes.Last()>=u->nodes.Last()) { // transcript includes u
				if((transfrag[t]->pattern & u->pattern) == u->pattern) { // t contains super-read
					available-=transfrag[t]->abundance;
					if(available<epsilon) {
						available=0;
						break;
					}
				}
				else { // this could still be compatible with u
					int mini=-1; // min node in t that appears in u
					int i=0;
					while(i<u->nodes.Count()) {
						if(transfrag[t]->pattern[u->nodes[i]]) { // first node in u that is in t
							mini=0;
							while(transfrag[t]->nodes[mini]!=u->nodes[i]) mini++;
							break;
						}
						i++;
					}
					if(mini>=0) { // found node in t -> check if it doesn't conflict with u
						bool conflict=false;
						if(i) { // only in this case I might have conflict
							if(no2gnode[u->nodes[i]]->parentpat[transfrag[t]->nodes[mini-1]]) {
								int *pos=gpos[edge(transfrag[t]->nodes[mini-1],transfrag[t]->nodes[mini],gno)];
								if(pos && transfrag[t]->pattern[*pos]) conflict=true;
							}
							else conflict=true;
						}
						if(!conflict) {
							int maxi=-1; // max node in t that appears in u
							i=u->nodes.Count()-1;
							while(i>=0) {
								if(transfrag[t]->pattern[u->nodes[i]]) { // last node in u that is in t
									maxi=transfrag[t]->nodes.Count()-1;
									while(transfrag[t]->nodes[maxi]!=u->nodes[i]) maxi--;
									break;
								}
								i--;
							}
							if(i<u->nodes.Count()-1) { // only in this case I might have conflict
								if(no2gnode[u->nodes[i]]->childpat[transfrag[t]->nodes[maxi+1]]) {
									int *pos=gpos[edge(transfrag[t]->nodes[maxi],transfrag[t]->nodes[maxi+1],gno)];
									if(pos && transfrag[t]->pattern[*pos]) conflict=true;
								}
								else conflict=true;
							}
							if(!conflict) { // check that the transcript is really not conflicting to u
								for(i=mini+1;i<=maxi;i++) {
									int *pos=gpos[edge(transfrag[t]->nodes[i-1],transfrag[t]->nodes[i],gno)];
									if(pos && transfrag[t]->pattern[*pos] && !u->pattern[*pos]) {
										conflict=true;
										break;
									}
								}
								if(!conflict) {
									if(transfrag[t]->nodes[0]==u->nodes[0] && transfrag[t]->nodes.Last()==u->nodes.Last())
										equalinc.Add(t);
									else equalext.Add(t);
								}
							}
						}
					}
				}
			}
			else if((transfrag[t]->pattern & u->pattern) == transfrag[t]->pattern) { // compatible transcripts; transcript is included in u
				seltrfrag.Add(transfrag[t]);
			}

		}
		if(available) {
			for(int i=0;i<equalinc.Count();i++) {
				int t=equalinc[i];
				float solveabund=transfrag[t]->abundance;
				if(available<solveabund) {
					solveabund=available;
				}
				transfrag[t]->abundance-=solveabund;
				transfrag[t]->srabund-=solveabund;
				if(transfrag[t]->abundance<epsilon) transfrag[t]->abundance=0;
				u->abundance+=solveabund;
				available-=solveabund;
				if(available<epsilon) {
					available=0;
					break;
				}
			}
			if(available) {
				for(int e=0;e<equalext.Count();e++) {
					int t=equalext[e];
					int j=0;
					for(int i=0;i<u->nodes.Count();i++) {
						if(!transfrag[t]->pattern[u->nodes[i]]) {
							while(transfrag[t]->nodes[j]<u->nodes[i]) j++;
							transfrag[t]->nodes.Insert(j,u->nodes[i]);
						}
					}
					transfrag[t]->pattern = transfrag[t]->pattern | u->pattern;
					available -= transfrag[t]->abundance;
					if(available<epsilon) {
						available=0;
						break;
					}
				}
				if(available) { // if I still have abundance left (very likely) -> stitch together transfrags in selfrags

					for(int i=1;i<u->nodes.Count();i++) { // compute abundance that needs to be explained at each node
						int n=u->nodes[i];
						for(int j=0;j<no2gnode[n]->trf.Count();j++) if(transfrag[no2gnode[n]->trf[j]]->nodes[0]==n) { // transfrag starts at node
							int t=no2gnode[n]->trf[j];
							if((transfrag[t]->pattern & u->pattern) == transfrag[t]->pattern) {
								seltrfrag.Add(transfrag[t]);
							}
						}
					}

					if(seltrfrag.Count()) {
						seltrfrag.Sort(trCmp);
						GVec<int> *start=new GVec<int>[gno]; // for each node remembers it's neighbours
						for(int i=0;i<seltrfrag.Count();i++) {
							start[seltrfrag[i]->nodes[0]].Add(i);
						}
						int i=0;
						while(available && i<seltrfrag.Count()) {
							if(seltrfrag[i]->abundance) {
								GVec<int> next; // this stores the augmenting path
								float maxabund=stitch_trf(u,available,i,i,seltrfrag,start,next);
								if(maxabund) {
									next.Add(i);
									for(int j=0;j<next.Count();j++) {
										int t=next[j];
										seltrfrag[t]->abundance-=maxabund;
										seltrfrag[t]->srabund-=maxabund;
										if(seltrfrag[t]->abundance<epsilon) seltrfrag[t]->abundance=0;
									}
									if(!seltrfrag[i]->abundance) i++; // I solved this transfrag
									available-=maxabund;
									if(available<epsilon) {
										available=0;
									}
									u->abundance+=maxabund;
								}
								else i++;
							}
							else i++;
						}
						delete [] start;
					}
				}
			}
		}
	}
}

int compatible_long(int* t,int *len,GPVec<CTransfrag>& transfrag,GPVec<CGraphnode>& no2gnode,int gno,GIntHash<int> &gpos) {
// return 0 if no compatibility, then 4*(rets)+rete where rets/rete = 1 if t[0] has extra intron, 2 if t[1] has extra intron, 3 if compatible starts/ends
// len[0] is how much t[0] extends past t[1] if positive, otherwise how much t[1] extends past t[0] at the start of transcripts
// len[1] positive: if t[0] has extra intron at start how much is t[1] inside of the intron -> default is 0; negative: same for t[1]
// len[2] is how much t[0] extends past t[1] if positive, otherwise how much t[1] extends past t[0]	at the end of transcripts
// len[3] positive: if t[0] has extra intron at end how much is t[1] inside of the intron -> default is 0; negative: same for t[1]


	if(transfrag[t[0]]->nodes[0]!=transfrag[t[1]]->nodes[0] &&
			no2gnode[transfrag[t[0]]->nodes[0]]->hardstart && no2gnode[transfrag[t[1]]->nodes[0]]->hardstart)
		return 0; // both starts can't be hard

	if(transfrag[t[0]]->nodes.Last()!=transfrag[t[1]]->nodes.Last() &&
			no2gnode[transfrag[t[0]]->nodes.Last()]->hardend && no2gnode[transfrag[t[1]]->nodes.Last()]->hardend)
		return 0; // both ends can't be hard

	uint tstart[2]={no2gnode[transfrag[t[0]]->nodes[0]]->start,no2gnode[transfrag[t[1]]->nodes[0]]->start};
	if(transfrag[t[0]]->longstart) tstart[0]=transfrag[t[0]]->longstart;
	if(transfrag[t[1]]->longstart) tstart[1]=transfrag[t[1]]->longstart;
	uint tend[2]={no2gnode[transfrag[t[0]]->nodes.Last()]->end,no2gnode[transfrag[t[1]]->nodes.Last()]->end};
	if(transfrag[t[0]]->longend) tend[0]=transfrag[t[0]]->longend;
	if(transfrag[t[1]]->longend) tend[1]=transfrag[t[1]]->longend;

	int f=0; // first starting transcript
	int s=1; // second starting transcript
	if(tstart[1]<tstart[0]) { // transfrag t1 starts before trf t0
		f=1;
		s=0;
	}

	if(transfrag[t[f]]->nodes.Last()<transfrag[t[s]]->nodes[0]) // no overlapping transcripts
		return 0;

	// check START
	int rets=3; // compatible starts by default
	int i[2]; // starting index in transfrags where they share a node
	i[0]=0;i[1]=0;
	len[0]=0;
	len[1]=0;
	while(transfrag[t[f]]->nodes[i[f]]<transfrag[t[s]]->nodes[0]) { // i[f]<transfrag[t[f]]->nodes.Count() always true due to condition of overlap before
		if(no2gnode[transfrag[t[f]]->nodes[i[f]]]->end+1<no2gnode[transfrag[t[f]]->nodes[i[f]+1]]->start) { // found intron in f
			rets=1+f;
			if(no2gnode[transfrag[t[s]]->nodes[0]]->hardstart) return 0; // there is no compatibility here
		}

		// len[0]+=no2gnode[transfrag[t[f]]->nodes[i[f]]]->len(); // this only approximates the length
		//  adj if I want real distances:
		if(i[f]) len[0]+=no2gnode[transfrag[t[f]]->nodes[i[f]]]->len();
		else len[0]+=no2gnode[transfrag[t[f]]->nodes[i[f]]]->end-tstart[f]+1;
		i[f]++;
	}

	if(transfrag[t[s]]->nodes.Last()<transfrag[t[f]]->nodes[i[f]]) return 0; // no overlap: t[s] is contained in t[f]

	// adj for real distance
	if(!i[f]) { // transfrag[t[f]]->nodes[0]==transfrag[t[s]]->nodes[0]
		len[0]=tstart[s]-tstart[f];
	}

	if(transfrag[t[s]]->nodes[0]<transfrag[t[f]]->nodes[i[f]]) { // i[f]>0 because f has to start first
		if(!no2gnode[transfrag[t[s]]->nodes[0]]->childpat[transfrag[t[f]]->nodes[i[f]]]) return 0; // I can not reach node i[f] from t[s][0]
		// to approximate:
		// len[1]=no2gnode[transfrag[t[f]]->nodes[i[f]]]->start-no2gnode[transfrag[t[s]]->nodes[0]]->start;
		// adj if I want real distances:
		len[1]=no2gnode[transfrag[t[f]]->nodes[i[f]]]->start-tstart[s];
		while(transfrag[t[s]]->nodes[i[s]]<transfrag[t[f]]->nodes[i[f]]) {
			i[s]++;
			if(no2gnode[transfrag[t[s]]->nodes[i[s]-1]]->end+1<no2gnode[transfrag[t[s]]->nodes[i[s]]]->start) { // not continuous
				return 0;
			}
		}
	}
	else if(i[f]){ // no intron -> adj lentgh for real distance
		len[0]+=tstart[s]-no2gnode[transfrag[t[s]]->nodes[0]]->start;
	}

	if(f) len[0]=-len[0]; // a negative value signals that t[1] extends past t[0]

	// check END
	int rete=3;
	f=0;
	s=1;
	if(tend[1]>tend[0]) { // f should be the longer transcript at the end
		f=1;
		s=0;
	}

	int j[2];
	j[0]=transfrag[t[0]]->nodes.Count()-1;
	j[1]=transfrag[t[1]]->nodes.Count()-1;

	len[2]=0;
	len[3]=0;
	while(transfrag[t[f]]->nodes[j[f]]>transfrag[t[s]]->nodes.Last()) { // s and f should have a node in common by now so fj doesn't have to go all the way back to 0
		if(no2gnode[transfrag[t[f]]->nodes[j[f]-1]]->end+1<no2gnode[transfrag[t[f]]->nodes[j[f]]]->start) { // found intron
			rete=1+f;
			if(rets<3 && rete!=rets) return 0;
			if(no2gnode[transfrag[t[s]]->nodes.Last()]->hardend) return 0; // there is no compatibility here
		}
		// to approximate:
		//len[2]+=no2gnode[transfrag[t[f]]->nodes[j[f]]]->len();
		// adj:
		if(j[f]<transfrag[t[f]]->nodes.Count()-1) len[2]+=no2gnode[transfrag[t[f]]->nodes[j[f]]]->len();
		else len[2]+=tend[f]-no2gnode[transfrag[t[f]]->nodes[j[f]]]->start+1;

		j[f]--;
	}
	if(transfrag[t[s]]->nodes[0]>transfrag[t[f]]->nodes[j[f]]) return 0; // no overlap (this shouldn't be the case since there is one node in common)

	if(j[f]==transfrag[t[f]]->nodes.Count()-1) len[2]+=tend[f]-tend[s];

	if(transfrag[t[s]]->nodes[j[s]]>transfrag[t[f]]->nodes[j[f]]) {
		if(!no2gnode[transfrag[t[f]]->nodes[j[f]]]->childpat[transfrag[t[s]]->nodes[j[s]]]) return 0; // I can not reach node j[s] from j[f]
		// len[3]=no2gnode[transfrag[t[s]]->nodes.Last()]->end-no2gnode[transfrag[t[f]]->nodes[j[f]]]->end; // to approximate dist
		//  adj if I want real distances:
		len[3]=tend[s]-no2gnode[transfrag[t[f]]->nodes[j[f]]]->end;
		while(transfrag[t[s]]->nodes[j[s]]>transfrag[t[f]]->nodes[j[f]]) {
			j[s]--;
			if(no2gnode[transfrag[t[s]]->nodes[j[s]]]->end+1<no2gnode[transfrag[t[s]]->nodes[j[s]+1]]->start) { // not continuous
				return 0;
			}
		}
	}
	else if(j[f]<transfrag[t[f]]->nodes.Count()-1) {
		len[2]+=no2gnode[transfrag[t[s]]->nodes.Last()]->end-tend[s];
	}

	if(i[f]>j[f] || i[s]>j[s]) return 0; // there should be i<=j

	if(f) len[2]=-len[2]; // a negative value signals that t[1] extends past t[0]

	// now if and is point at the same node in transcripts, and so do jf and js
	while(i[f]<=j[f] && i[s]<=j[s]) {
		if(transfrag[t[f]]->nodes[i[f]]==transfrag[t[s]]->nodes[i[s]]) { // skip equal nodes
			i[f]++;
			i[s]++;
		}
		else { // unequal
			if(transfrag[t[f]]->nodes[i[f]]>transfrag[t[s]]->nodes[i[s]]) { // make sure node f is always smallest
				int k=f;
				f=s;
				s=k;
			}
			while(transfrag[t[f]]->nodes[i[f]]<transfrag[t[s]]->nodes[i[s]] &&
					no2gnode[transfrag[t[f]]->nodes[i[f]-1]]->end+1==no2gnode[transfrag[t[f]]->nodes[i[f]]]->start) { // advance smallest node until I get to a gap
				i[f]++;
			}
			if(no2gnode[transfrag[t[f]]->nodes[i[f]-1]]->end+1==no2gnode[transfrag[t[f]]->nodes[i[f]]]->start) return 0; // gap in s is filled by transfrag f
			// see if gaps are both hard
			int *pos=gpos[edge(transfrag[t[f]]->nodes[i[f]-1],transfrag[t[f]]->nodes[i[f]],gno)];
			if(pos && transfrag[t[f]]->pattern[*pos]) { // this is a hard edge in t[f]
				pos=gpos[edge(transfrag[t[s]]->nodes[i[s]-1],transfrag[t[s]]->nodes[i[s]],gno)];
				if(pos && transfrag[t[s]]->pattern[*pos]) return 0; // edge is hard in both f and s
			}
			if(transfrag[t[f]]->nodes[i[f]]>transfrag[t[s]]->nodes[i[s]]) { // make sure node f is always smallest
				int k=f;
				f=s;
				s=k;
			}
			while(transfrag[t[f]]->nodes[i[f]]<transfrag[t[s]]->nodes[i[s]]) {
				if(no2gnode[transfrag[t[f]]->nodes[i[f]]]->end+1<no2gnode[transfrag[t[f]]->nodes[i[f]+1]]->start) return 0;
				i[f]++;
			}
		}
	}

	if(rets==1 || rete==1) return(1); // t[0] extends with introns past t[1] at at least one side
	if(rets==2 || rete==2) return(2); // t[1] extends with introns past t[0] at at least one side
	return(3); // t[0] and t[1] have compatible starts/ends
}

void process_transfrags(int s, int gno,int edgeno,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,CTreePat *tr2no,
		GIntHash<int> &gpos,GVec<CGuide>& guidetrf,GList<CPrediction>& pred,GVec<int>& trflong) {

	/*
	{ // DEBUG ONLY
		printTime(stderr);
		fprintf(stderr,"There are %d transfrags before clean up:\n",transfrag.Count());
		for(int i=0;i<transfrag.Count();i++) {
			fprintf(stderr,"transfrag[%d](%f,%f) long=%d short=%d usepath=%d:",i,transfrag[i]->abundance,transfrag[i]->srabund,transfrag[i]->longread,transfrag[i]->shortread,(int)transfrag[i]->usepath);
			for(int j=0;j<transfrag[i]->nodes.Count();j++) fprintf(stderr," %d",transfrag[i]->nodes[j]);
			fprintf(stderr,"\n");
		}
	}
	*/

	// add all guide patterns to the set of transfrags so that I can have a "backbone" for each guide
	// I need this because there might be an incompatible transfrag connecting the nodes in the guide
	//fprintf(stderr,"There are %d guides\n",guidetrf.Count());
	for(int i=0;i<guidetrf.Count();i++) if(guidetrf[i].trf->guide){

		CTransfrag *t=NULL;
		bool add=true;
		if(longreads || mixedMode) {
			/*guidetrf[i].trf->pattern[0]=0;
			guidetrf[i].trf->pattern[gno-1]=0;
			int *pos=gpos[edge(0,guidetrf[i].trf->nodes[1],gno)];
			if(pos) guidetrf[i].trf->pattern[*pos]=0;
			pos=gpos[edge(guidetrf[i].trf->nodes[guidetrf[i].trf->nodes.Count()-2],guidetrf[i].trf->nodes.Last(),gno)];
			if(pos) guidetrf[i].trf->pattern[*pos]=0;
			guidetrf[i].trf->nodes.Pop();
			guidetrf[i].trf->nodes.Shift();*/

			t=findtrf_in_treepat(gno,gpos,guidetrf[i].trf->nodes,guidetrf[i].trf->pattern,tr2no); // I need to adjust first/last node
			if(!t) { // t is NULL
				float abund=0;
				if(mixedMode) abund=trthr*ERROR_PERC;
				t=new CTransfrag(guidetrf[i].trf->nodes,guidetrf[i].trf->pattern,abund);
				t->longread=true;
			}
			else add=false;
		}
		else {
			t=new CTransfrag(guidetrf[i].trf->nodes,guidetrf[i].trf->pattern,trthr*ERROR_PERC);
		}

		if(!longreads) {
			if(includesource) {
				guidetrf[i].trf->nodes.Insert(0,0); // I need to comment this if I need path not to include the source
				guidetrf[i].trf->pattern[0]=1;
				int *pos=gpos[edge(0,guidetrf[i].trf->nodes[1],gno)];
				if(pos) guidetrf[i].trf->pattern[*pos]=1;
			}
			int sink=gno-1;
			guidetrf[i].trf->nodes.Add(sink);
			guidetrf[i].trf->pattern[sink]=1;
			int *pos=gpos[edge(guidetrf[i].trf->nodes[guidetrf[i].trf->nodes.Count()-2],guidetrf[i].trf->nodes.Last(),gno)];
			if(pos) guidetrf[i].trf->pattern[*pos]=1;
		}



		/*
		float abund=0;
		if(!longreads) {
			abund=trthr*ERROR_PERC;
		}
		CTransfrag *t=new CTransfrag(guidetrf[i].trf->nodes,guidetrf[i].trf->pattern,abund);
		*/
		/*
		{ // DEBUG ONLY
			fprintf(stderr,"Add guidetrf with nodes:");
			for(int j=0;j<guidetrf[i].trf->nodes.Count();j++) fprintf(stderr," %d",guidetrf[i].trf->nodes[j]);
			//fprintf(stderr," and pattern: ");
			//printBitVec(guidetrf[i].trf->pattern);
			fprintf(stderr,"\n");
			fprintf(stderr,"Found transfrag %p with nodes:",t);
			for(int j=0;j<t->nodes.Count();j++) fprintf(stderr," %d",t->nodes[j]);
			//fprintf(stderr," and pattern: ");
			//printBitVec(guidetrf[i].trf->pattern);
			fprintf(stderr,"\n");

		}
		*/

		/*if(!longreads) {
			t->pattern[0]=0;
			t->pattern[gno-1]=0;
			int *pos=gpos[edge(0,t->nodes[1],gno)];
			if(pos) t->pattern[*pos]=0;
			pos=gpos[edge(t->nodes[t->nodes.Count()-2],t->nodes.Last(),gno)];
			if(pos) t->pattern[*pos]=0;
			if(add) {
				t->nodes.Pop();
				t->nodes.Shift();
			}
		}*/
		t->guide=1+guidetrf[i].g;
		//fprintf(stderr,"t->guide set to=%d\n",1+guidetrf[i].g);
		t->longstart=no2gnode[t->nodes[0]]->start;
		t->longend=no2gnode[t->nodes.Last()]->end;
		//if(longreads) t->usepath=guidetrf[i].g; // guide index
		no2gnode[t->nodes[0]]->hardstart=1;  // I can always trust a guide's start
		no2gnode[t->nodes.Last()]->hardend=1; // I can always trust a guide's end
		if(add) transfrag.Add(t);
	}
	else if(!longreads) {
		if(includesource) {
			guidetrf[i].trf->nodes.Insert(0,0); // I need to comment this if I need path not to include the source
			guidetrf[i].trf->pattern[0]=1;
			int *pos=gpos[edge(0,guidetrf[i].trf->nodes[1],gno)];
			if(pos) guidetrf[i].trf->pattern[*pos]=1;
		}
		int sink=gno-1;
		guidetrf[i].trf->nodes.Add(sink);
		guidetrf[i].trf->pattern[sink]=1;
		int *pos=gpos[edge(guidetrf[i].trf->nodes[guidetrf[i].trf->nodes.Count()-2],guidetrf[i].trf->nodes.Last(),gno)];
		if(pos) guidetrf[i].trf->pattern[*pos]=1;
	}

	GPVec<CTransfrag> srfrag(false);
	// eliminate transfrags below threshold (they represent noise) if they don't come from source
	if(!eliminate_transfrags_under_thr(gno,gpos,transfrag,tr2no,trthr,srfrag) && srfrag.Count()) { // long transfrags but only to source/sink
		srfrag.Clear();
	}

	/*
	{ // DEBUG ONLY
		//printTime(stderr);
		fprintf(stderr,"\nThere are %d transfrags after clean up:\n",transfrag.Count());
		for(int i=0;i<transfrag.Count();i++) {
			fprintf(stderr,"transfrag[%d] abund=%f long:%d: short:%d usepath:%d",i,transfrag[i]->abundance,transfrag[i]->longread,transfrag[i]->shortread,(int)transfrag[i]->usepath);
			for(int j=0;j<transfrag[i]->nodes.Count();j++) fprintf(stderr," %d",transfrag[i]->nodes[j]);
			fprintf(stderr,"\n");
		}
	}
	*/

	GBitVec allpat(gno+edgeno);

	bool trsort=true;

	if(longreads) { // add source/sink links but only if they need to be added to explain the traversals in the graph
		transfrag.Sort(longtrCmp); // most abundant transfrag in the graph come first, then the ones with most nodes, then the ones more complete
		trsort=false;
		int source=0;
		int sink=gno-1;
		GVec<int> hassource(gno,-1); // remembers transcript number that links given node to source
		GVec<int> hassink(gno,-1); // remembers transcript number that links given node to sink
		GBitVec keepsource(gno); // if not set then I can remove link from node to source; keeps source link if it exists otherwise otherwise
		GBitVec keepsink(gno); // if not set then I can remove link from node to sink; keeps sink link if it exists otherwise otherwise
		GVec<CLongTrf> keeptrf; // keeps all potential transfrags that will be kept from most abundant to least, unassembled
		float zero=0;
		GVec<float> addsource(gno,zero);
		GVec<float> addsink(gno,zero);
		int edgedist=CHI_WIN; // I need to be consistent (if I change here then I need to change in update_abundance too)
		int ssdist=longintronanchor;

		/*
		{ // DEBUG ONLY
			//printTime(stderr);
			fprintf(stderr,"\nThere are %d transfrags after clean up:\n",transfrag.Count());
			for(int i=0;i<transfrag.Count();i++) {
				fprintf(stderr,"transfrag[%d] abund=%f:",i,transfrag[i]->abundance);
				for(int j=0;j<transfrag[i]->nodes.Count();j++) fprintf(stderr," %d",transfrag[i]->nodes[j]);
				if(transfrag[i]->guide) fprintf(stderr," guide");
				fprintf(stderr,"\n");
			}
		}
		*/

		for(int t1=0;t1<transfrag.Count();t1++) {
			/*fprintf(stderr,"Consider t=%d with abund=%f guide=%d and nodes:",t1,transfrag[t1]->abundance,transfrag[t1]->guide);
			for(int j=0;j<transfrag[t1]->nodes.Count();j++) {
				if(j) {
					int *pos=gpos[edge(transfrag[t1]->nodes[j-1],transfrag[t1]->nodes[j],gno)];
					if(pos && transfrag[t1]->pattern[*pos]) {
						fprintf(stderr,"-");
					}
					else fprintf(stderr," ");
				}
				fprintf(stderr,"%d",transfrag[t1]->nodes[j]);
			} fprintf(stderr,"\n");*/
			if(!transfrag[t1]->nodes[0]) {
				hassource[transfrag[t1]->nodes[1]]=t1;
				//fprintf(stderr,"Node %d in t=%d with cov=%f has source\n",transfrag[t1]->nodes[1],t1,transfrag[t1]->abundance);
			}
			else if(transfrag[t1]->nodes.Last()==gno-1) {
				hassink[transfrag[t1]->nodes[0]]=t1;
				//fprintf(stderr,"Node %d in t=%d with cov=%f has sink\n",transfrag[t1]->nodes[0],t1,transfrag[t1]->abundance);
			}
			else {
				if(eonly && !transfrag[t1]->guide) continue; // do not remember transfrags that are not guides
				if(!keepsource[transfrag[t1]->nodes[0]]) {
					if(transfrag[t1]->longstart) keepsource[transfrag[t1]->nodes[0]]=1; //fprintf(stderr,"keep source %d\n",transfrag[t1]->nodes[0]);}
					else if(no2gnode[transfrag[t1]->nodes[0]]->hardstart) keepsource[transfrag[t1]->nodes[0]]=1;
				}
				if(!keepsink[transfrag[t1]->nodes.Last()]) {
					if(transfrag[t1]->longend) keepsink[transfrag[t1]->nodes.Last()]=1;//fprintf(stderr,"keep sink %d\n",transfrag[t1]->nodes.Last());}
					else if(no2gnode[transfrag[t1]->nodes.Last()]) keepsink[transfrag[t1]->nodes.Last()]=1;
				}
				bool included=false;
				// a transfrag that starts at source and ends at sink can never be included in a kept transfrag, so I am safe to do next
				for(int t2=0; t2<keeptrf.Count();t2++) {
					int t[2]={t1,keeptrf[t2].t}; // t1 current, t2 the one I kept
					int len[4]={MAX_NODE,MAX_NODE,MAX_NODE,MAX_NODE};
					int ret=compatible_long(t,len,transfrag,no2gnode,gno,gpos);
					//fprintf(stderr,"  ret=%d t[0]=%d t[1]=%d len[0]=%d len[1]=%d len[2]=%d len[3]=%d\n",ret,t[0],t[1],len[0],len[1],len[2],len[3]);
					if(ret){
						switch(ret) {
						case 1: // t[0] includes t[1]: it extends with introns on on one or both sides -> keep unless I can eliminate a previous one
							if(!transfrag[t[1]]->guide && transfrag[t1]->longstart && transfrag[t1]->longend && // t[1] might be included in t[0] so I might eliminate if it doesn't pass threshold
									(!no2gnode[transfrag[t[1]]->nodes[0]]->hardstart || transfrag[t[0]]->nodes[0] == transfrag[t[1]]->nodes[0] ) &&
									(!no2gnode[transfrag[t[1]]->nodes.Last()]->hardend || transfrag[t[0]]->nodes.Last() == transfrag[t[1]]->nodes.Last())) {
								//if(transfrag[t[0]]->abundance>(1-ERROR_PERC/DROP)*transfrag[t[1]]->abundance) { // t[0] is within limits of t[1]
								if(transfrag[t[0]]->abundance>DROP*transfrag[t[1]]->abundance) { // t[0] is within limits of t[1]
									if(len[1]<ssdist && len[3]<ssdist) { // prefer t[0] instead of t[1]
										keeptrf[t2].t=t1;
										keeptrf[t2].cov+=transfrag[t1]->abundance;
										keeptrf[t2].group.Add(t1);
										included=true; // I do not want to store transcript
										//fprintf(stderr,"trf %d includes %d\n",t[0],t[1]);
									}
								}
							}
							break;
						case 2: // t[1] includes t[0]: extends with introns past ends of t[0] (t[1] possibly includes t[0]); t[1] is more abundant than t0
							//if(transfrag[t[1]]->guide || transfrag[t[1]]->abundance>(1-ERROR_PERC/DROP)*transfrag[t[0]]->abundance) {
							if(!transfrag[t[0]]->guide &&
									(!no2gnode[transfrag[t[0]]->nodes[0]]->hardstart || transfrag[t[0]]->nodes[0] == transfrag[t[1]]->nodes[0]) &&
									(!no2gnode[transfrag[t[0]]->nodes.Last()]->hardend || transfrag[t[0]]->nodes.Last() == transfrag[t[1]]->nodes.Last())) {
								if(len[1]<ssdist && len[3]<ssdist) {
									keeptrf[t2].cov+=transfrag[t1]->abundance;
									keeptrf[t2].group.Add(t1);
									included=true;
									//fprintf(stderr,"trf %d is intronic including %d\n",t[1],t[0]);
								}
							}
							//}
							break;
						case 3: // t1 and t0 are compatible --> just look for the edges; t1 goes further apart
							if(transfrag[t[0]]->nodes[0]!=transfrag[t[1]]->nodes[0] && transfrag[t[0]]->nodes.Last()!=transfrag[t[1]]->nodes.Last() &&
									((no2gnode[transfrag[t[0]]->nodes[0]]->hardstart && !no2gnode[transfrag[t[1]]->nodes[0]]->hardstart &&
											!no2gnode[transfrag[t[0]]->nodes.Last()]->hardend && no2gnode[transfrag[t[1]]->nodes.Last()]->hardend) ||
											(!no2gnode[transfrag[t[0]]->nodes[0]]->hardstart && no2gnode[transfrag[t[1]]->nodes[0]]->hardstart &&
													no2gnode[transfrag[t[0]]->nodes.Last()]->hardend && !no2gnode[transfrag[t[1]]->nodes.Last()]->hardend))) {
								// these two transcripts both have one good start and one different --> keep them both (different option would be to add them to another transfrag that is compatible and has both hardends but then it's more complicated
								break;
							}

							// I keep both if both are guides
							//if((!transfrag[t[0]]->guide || !transfrag[t[1]]->guide) && abs(len[0])<edgedist && abs(len[2])<edgedist) { // close by
							if((!transfrag[t[0]]->guide || !transfrag[t[1]]->guide) && len[0]<edgedist && len[2]<edgedist) { // close by
								if(transfrag[t[0]]->guide || (!transfrag[t[1]]->guide && no2gnode[transfrag[t[0]]->nodes[0]]->hardstart && no2gnode[transfrag[t[0]]->nodes.Last()]->hardend))
									keeptrf[t2].t=t1; // t[0] to replace t[1]
								keeptrf[t2].cov+=transfrag[t1]->abundance;
								keeptrf[t2].group.Add(t1);
								included=true;
								//fprintf(stderr,"trf %d %d equivalent start/ends\n",t[1],t[0]);
							}
							break;
						}
						if(included) break; // break from for loop
					}
				}

				if(!included){
					if(transfrag[t1]->guide || ((transfrag[t1]->longstart || no2gnode[transfrag[t1]->nodes[0]]->hardstart)&&
							(transfrag[t1]->longend || no2gnode[transfrag[t1]->nodes.Last()]->hardend))) { // if this is not included and has correct start/end
						CLongTrf kt(t1,transfrag[t1]->abundance);
						keeptrf.Add(kt);
						keeptrf.Last().group.Add(t1);
						//fprintf(stderr,"keep transfrag %d\n",t1);
					}
					else { // incomplete transcript, possibly wrong
						transfrag[t1]->weak=1;
						//fprintf(stderr,"Incomplete transcript %d\n",t1);
					}
				}
			}
		}


		char sign='-';
		if(s) { sign='+';}
		//GBitVec guidesource(gno);
		//GBitVec guidesink(gno);
		//for(int i=0;i<keeptrf.Count();i++) {
		for(int i=keeptrf.Count()-1;i>=0;i--) { // I add the kept transcripts to trflong from least significant to most in order to make deletion easier

			//fprintf(stderr,"Build source/sink for transfrag %d\n",keeptrf[i].t);
			int n1=transfrag[keeptrf[i].t]->nodes[0];
			int n2=transfrag[keeptrf[i].t]->nodes.Last();
			//if(hassource[n1]<0 || hassink[n2]<0) fprintf(stderr,"Build source/sink for transfrag %d\n",keeptrf[i].t);

			if(!rawreads && !transfrag[keeptrf[i].t]->guide && ((!no2gnode[n1]->hardstart && (hassource[n1]<0 || !keepsource[n1])) ||
					(!no2gnode[n2]->hardend && (hassink[n2]<0 || !keepsink[n2])))) {
			//if(!rawreads && (no2gnode[n1]->hardstart || (hassource[n1]>=0 && keepsource[n1])) && (no2gnode[n2]->hardend ||(hassink[n2]>=0 && keepsink[n2]))) {
				trflong.Add(keeptrf[i].t);
			}

			/******* previous implementation here

			addsource[n1]=keeptrf[i].cov;
			addsink[n2]=keeptrf[i].cov;
			 *******/

			/*
			if(!addsource[n1] && hassource[n1]<0) {
				int startpos=no2gnode[n1]->start-refstart;
				if(startpos-CHI_THR<0 || startpos+CHI_THR>bpcov->Count()) addsource[n1]=1;
				else {
					addsource[n1]=(get_cov(1,startpos,startpos+CHI_THR-1,bpcov)-get_cov(2-2*s,startpos,startpos+CHI_THR-1,bpcov)-
							get_cov(1,startpos-CHI_THR,startpos-1,bpcov)+get_cov(2-2*s,startpos-CHI_THR,startpos-1,bpcov))/(DROP*CHI_THR);
				}
			}

			if(!addsink[n2]  && hassink[n1]<0) {
				int endpos=no2gnode[n2]->end-refstart;
				if(endpos-CHI_THR<0 || endpos+CHI_THR>bpcov->Count()-1) addsink[n2]=1;
				else {
					addsink[n2]=(get_cov(1,endpos-CHI_THR+1,endpos,bpcov)-get_cov(2-2*s,endpos-CHI_THR+1,endpos,bpcov)-
							get_cov(1,endpos+1,endpos+CHI_THR,bpcov)+get_cov(2-2*s,endpos+1,endpos+CHI_THR,bpcov))/(DROP*CHI_THR);
				}
			}*/

			// all
			for(int j=0;j<keeptrf[i].group.Count();j++) {
				if(n1==transfrag[keeptrf[i].group[j]]->nodes[0]) {
					addsource[n1]+=transfrag[keeptrf[i].group[j]]->abundance;
					//fprintf(stderr,"Add source t[%d]->cov=%f to node %d = %f\n",keeptrf[i].group[j],transfrag[keeptrf[i].group[j]]->abundance,n1,addsource[n1]);
				}
				if(n2==transfrag[keeptrf[i].group[j]]->nodes.Last()) {
					addsink[n2]+=transfrag[keeptrf[i].group[j]]->abundance;
					//fprintf(stderr,"Add sink t[%d]->cov=%f to node %d = %f\n",keeptrf[i].group[j],transfrag[keeptrf[i].group[j]]->abundance,n2,addsink[n2]);
				}
			}


			if(rawreads) {
				GVec<GSeg> exons;
				int j=0;
				int len=0;
				int t=keeptrf[i].t;
				while(j<transfrag[t]->nodes.Count()) {
					int nodestart=no2gnode[transfrag[t]->nodes[j]]->start;
					int nodeend=no2gnode[transfrag[t]->nodes[j]]->end;
					len+=nodeend-nodestart+1;
					while(j+1<transfrag[t]->nodes.Count() && no2gnode[transfrag[t]->nodes[j]]->end+1==no2gnode[transfrag[t]->nodes[j+1]]->start) {
						j++;
						len+=no2gnode[transfrag[t]->nodes[j]]->len();
						nodeend=no2gnode[transfrag[t]->nodes[j]]->end;
					}
					GSeg exon(nodestart,nodeend);
					exons.Add(exon);
					j++;
				}
				uint tstart=exons[0].start;
				uint tend=exons.Last().end;
				if(transfrag[t]->longstart>tstart) {
					len-=transfrag[t]->longstart-tstart;
					tstart=transfrag[t]->longstart;
				}
				if(transfrag[t]->longend && transfrag[t]->longend<tend) {
					len-=tend-transfrag[t]->longend;
					tend=transfrag[t]->longend;
				}

				CPrediction *p=new CPrediction(s, NULL, tstart, tend, keeptrf[i].cov, sign, len);
				exons[0].start=tstart;
				exons.Last().end=tend;
				p->exons=exons;
				pred.Add(p);
			}
		}
		for(int i=keeptrf.Count()-1;i>=0;i--) {
			int n1=transfrag[keeptrf[i].t]->nodes[0];
			int n2=transfrag[keeptrf[i].t]->nodes.Last();

			if(transfrag[keeptrf[i].t]->guide || ((no2gnode[n1]->hardstart || (hassource[n1]>=0 && keepsource[n1])) && (no2gnode[n2]->hardend ||(hassink[n2]>=0 && keepsink[n2]))))
				trflong.Add(keeptrf[i].t);

		}

		/*
		{ // DEBUG ONLY
			fprintf(stderr,"%d keeptrf:\n",keeptrf.Count());
			for(int i=0;i<keeptrf.Count();i++) {
				fprintf(stderr,"(%d) %d abund=%f keepcov=%f",i,keeptrf[i].t,transfrag[keeptrf[i].t]->abundance,keeptrf[i].cov);
				for(int j=0;j<transfrag[keeptrf[i].t]->nodes.Count();j++) {
					fprintf(stderr," %d",transfrag[keeptrf[i].t]->nodes[j]);
				}
				fprintf(stderr,"\n");

			}
			fprintf(stderr,"%d trflong:",trflong.Count());
			for(int i=0;i<trflong.Count();i++)
				fprintf(stderr," %d",trflong[i]);
			fprintf(stderr,"\n");
		}
		*/


		// add source/sink connections
		for(int i=1;i<gno-1;i++) {
			//fprintf(stderr,"i=%d hassource=%d addsource=%f hassink=%d addsink=%f\n",i,hassource[i],addsource[i],hassink[i],addsink[i]);
			if(hassource[i]<0) {
				if(addsource[i] && no2gnode[i]->hardstart) {
					/*float abund=trthr;
					if(no2gnode[i]->hardstart) { // node i doesn't have source as parent but it should
						abund=addsource[i];
					}*/
					// else I am not that confident this is a start
					no2gnode[i]->parent.Insert(0,zero);
					no2gnode[0]->child.Add(i);
					GVec<int> nodes;
					nodes.Add(source);
					nodes.Add(i);
					//CTransfrag *t=new CTransfrag(nodes,allpat,abund);
					CTransfrag *t=new CTransfrag(nodes,allpat,addsource[i]);
					t->pattern[source]=1;
					t->pattern[i]=1;
					t->longread=true;
					transfrag.Add(t);
					//fprintf(stderr,"Add link to source: 0-%d with abundance=%f\n",i, addsource[i]);
				}
			}
			else {
				if(!keepsource[i]) { // it was a mistake to link this to source --> remove
					/*fprintf(stderr,"delete abundance of trf:");
					for(int j=0;j<transfrag[hassource[i]]->nodes.Count();j++) fprintf(stderr," %d",transfrag[hassource[i]]->nodes[j]);fprintf(stderr,"\n");*/
					transfrag[hassource[i]]->abundance=0;
				}
				else { // this one has source and it should be kept => update trf abundance to a more realistic value
					//fprintf(stderr,"source %d:%d abund=%f addsource=%f\n",i,no2gnode[i]->start,transfrag[hassource[i]]->abundance,addsource[i]);
					//if(guidesource[i] && transfrag[hassource[i]]->abundance<addsource[i])
					transfrag[hassource[i]]->abundance=addsource[i];
				}
			}

			if(hassink[i]<0) {
				if(addsink[i] && no2gnode[i]->hardend) {
					/*float abund=trthr;
					if(no2gnode[i]->hardend) {
						abund=addsink[i];
					}*/
					// else  not very confident about this end

					no2gnode[i]->child.Add(sink);
					no2gnode[sink]->parent.Add(i);
					GVec<int> nodes;
					nodes.Add(i);
					nodes.Add(sink);
					CTransfrag *t=new CTransfrag(nodes,allpat,addsink[i]);
					t->pattern[sink]=1;
					t->pattern[i]=1;
					t->longread=true;
					transfrag.Add(t);
					//fprintf(stderr,"Add link to sink: %d-%d with abundance=%f\n",i,sink,addsink[i]);
				}
			}
			else {
				if(!keepsink[i]) { // it was a mistake to link this to sink --> remove
					/*fprintf(stderr,"delete abundance of trf:");
					for(int j=0;j<transfrag[hassink[i]]->nodes.Count();j++) fprintf(stderr," %d",transfrag[hassink[i]]->nodes[j]);fprintf(stderr,"\n");*/
					transfrag[hassink[i]]->abundance=0;
				}
				else {
					//fprintf(stderr,"sink %d:%d abund=%f addsink=%f\n",i,no2gnode[i]->end,transfrag[hassink[i]]->abundance,addsink[i]);
					//if(guidesink[i] && transfrag[hassink[i]]->abundance<addsink[i])
					transfrag[hassink[i]]->abundance=addsink[i];
				}
			}
		}
	}
	else if(srfrag.Count() && mixedMode) { // add source/sink links but only if they need to be added to explain the traversals in the graph
			srfrag.Sort(longtrCmp); // most abundant transfrags in the graph come first, then the ones with most nodes, then the ones more complete
			int source=0;
			int sink=gno-1;
			GVec<int> hassource(gno,-1); // remembers transcript number that links given node to source
			GVec<int> hassink(gno,-1); // remembers transcript number that links given node to sink
			GBitVec keepsource(gno); // if not set then I can remove link from node to source; keeps source link if it exists otherwise otherwise
			GBitVec keepsink(gno); // if not set then I can remove link from node to sink; keeps sink link if it exists otherwise otherwise
			GVec<CLongTrf> keeptrf; // keeps all potential transfrags that will be kept from most abundant to least, unassembled
			float zero=0;
			GVec<float> addsource(gno,zero);
			GVec<float> addsink(gno,zero);
			int edgedist=CHI_WIN; // I need to be consistent (if I change here then I need to change in update_abundance too)
			int ssdist=longintronanchor;
			int ntrf=0; // number of trflong transcripts

			/*
			{ // DEBUG ONLY
				//printTime(stderr);
				fprintf(stderr,"\nThere are %d longsrtransfrags after clean up:\n",srfrag.Count());
				for(int i=0;i<srfrag.Count();i++) {
					fprintf(stderr,"srfrag[%d] longstart=%d longend=%d abund=%f t=%d:",i,srfrag[i]->longstart,srfrag[i]->longend,srfrag[i]->abundance,int(srfrag[i]->usepath));
					for(int j=0;j<srfrag[i]->nodes.Count();j++) fprintf(stderr," %d",srfrag[i]->nodes[j]);
					fprintf(stderr,"\n");
				}
			}
			*/

			for(int t1=0;t1<srfrag.Count();t1++) {
				/*fprintf(stderr,"Consider t=%d with abund=%f and nodes:",t1,srfrag[t1]->abundance);
				for(int j=0;j<srfrag[t1]->nodes.Count();j++) {
					if(j) {
						int *pos=gpos[edge(srfrag[t1]->nodes[j-1],srfrag[t1]->nodes[j],gno)];
						if(pos && srfrag[t1]->pattern[*pos]) {
							fprintf(stderr,"-");
						}
						else fprintf(stderr," ");
					}
					fprintf(stderr,"%d",srfrag[t1]->nodes[j]);
				} fprintf(stderr,"\n");*/
				if(!srfrag[t1]->nodes[0]) {
					hassource[srfrag[t1]->nodes[1]]=t1;
					//fprintf(stderr,"Node %d in t=%d with cov=%f has source\n",srfrag[t1]->nodes[1],t1,srfrag[t1]->abundance);
				}
				else if(srfrag[t1]->nodes.Last()==gno-1) {
					hassink[srfrag[t1]->nodes[0]]=t1;
					//fprintf(stderr,"Node %d in t=%d with cov=%f has sink\n",srfrag[t1]->nodes[0],t1,srfrag[t1]->abundance);
				}
				else {
					if(eonly && !srfrag[t1]->guide) continue; // do not remember transfrags that are not guides
					if(!keepsource[srfrag[t1]->nodes[0]]) {
						if(srfrag[t1]->longstart) keepsource[srfrag[t1]->nodes[0]]=1;
						else if(no2gnode[srfrag[t1]->nodes[0]]->hardstart) keepsource[srfrag[t1]->nodes[0]]=1;
					}
					if(!keepsink[srfrag[t1]->nodes.Last()]) {
						if(srfrag[t1]->longend) keepsink[srfrag[t1]->nodes.Last()]=1;
						else if(no2gnode[srfrag[t1]->nodes.Last()]) keepsink[srfrag[t1]->nodes.Last()]=1;
					}
					bool included=false;
					// a transfrag that starts at source and ends at sink can never be included in a kept transfrag, so I am safe to do next
					for(int t2=0; t2<keeptrf.Count();t2++) {
						int t[2]={t1,keeptrf[t2].t}; // t1 current, t2 the one I kept
						int len[4]={MAX_NODE,MAX_NODE,MAX_NODE,MAX_NODE};
						int ret=compatible_long(t,len,srfrag,no2gnode,gno,gpos);
						//fprintf(stderr,"  ret=%d t[0]=%d t[1]=%d len[0]=%d len[1]=%d len[2]=%d len[3]=%d\n",ret,t[0],t[1],len[0],len[1],len[2],len[3]);
						if(ret){
							switch(ret) {
							case 1: // t[0] includes t[1]: it extends with introns on on one or both sides -> keep unless I can eliminate a previous one
								if(!srfrag[t[1]]->guide && srfrag[t1]->longstart && srfrag[t1]->longend && // t[1] might be included in t[0] so I might eliminate if it doesn't pass threshold
										(!no2gnode[srfrag[t[1]]->nodes[0]]->hardstart || srfrag[t[0]]->nodes[0] == srfrag[t[1]]->nodes[0] ) &&
										(!no2gnode[srfrag[t[1]]->nodes.Last()]->hardend || srfrag[t[0]]->nodes.Last() == srfrag[t[1]]->nodes.Last())) {
									//if(srfrag[t[0]]->abundance>(1-ERROR_PERC/DROP)*srfrag[t[1]]->abundance) { // t[0] is within limits of t[1]
									if(srfrag[t[0]]->abundance>DROP*srfrag[t[1]]->abundance) { // t[0] is within limits of t[1]
										if(len[1]<ssdist && len[3]<ssdist) { // prefer t[0] instead of t[1]
											keeptrf[t2].t=t1;
											keeptrf[t2].cov+=srfrag[t1]->abundance;
											keeptrf[t2].group.Add(t1);
											included=true; // I do not want to store transcript
											//fprintf(stderr,"trf %d includes %d\n",t[0],t[1]);
										}
									}
								}
								break;
							case 2: // t[1] includes t[0]: extends with introns past ends of t[0] (t[1] possibly includes t[0]); t[1] is more abundant than t0
								//if(srfrag[t[1]]->guide || srfrag[t[1]]->abundance>(1-ERROR_PERC/DROP)*srfrag[t[0]]->abundance) {
								if(!srfrag[t[0]]->guide &&
										(!no2gnode[srfrag[t[0]]->nodes[0]]->hardstart || srfrag[t[0]]->nodes[0] == srfrag[t[1]]->nodes[0]) &&
										(!no2gnode[srfrag[t[0]]->nodes.Last()]->hardend || srfrag[t[0]]->nodes.Last() == srfrag[t[1]]->nodes.Last())) {
									if(len[1]<ssdist && len[3]<ssdist) {
										keeptrf[t2].cov+=srfrag[t1]->abundance;
										keeptrf[t2].group.Add(t1);
										included=true;
										//fprintf(stderr,"trf %d is intronic including %d\n",t[1],t[0]);
									}
								}
								//}
								break;
							case 3: // t1 and t0 are compatible --> just look for the edges; t1 goes further apart
								if(srfrag[t[0]]->nodes[0]!=srfrag[t[1]]->nodes[0] && srfrag[t[0]]->nodes.Last()!=srfrag[t[1]]->nodes.Last() &&
										((no2gnode[srfrag[t[0]]->nodes[0]]->hardstart && !no2gnode[srfrag[t[1]]->nodes[0]]->hardstart &&
												!no2gnode[srfrag[t[0]]->nodes.Last()]->hardend && no2gnode[srfrag[t[1]]->nodes.Last()]->hardend) ||
												(!no2gnode[srfrag[t[0]]->nodes[0]]->hardstart && no2gnode[srfrag[t[1]]->nodes[0]]->hardstart &&
														no2gnode[srfrag[t[0]]->nodes.Last()]->hardend && !no2gnode[srfrag[t[1]]->nodes.Last()]->hardend))) {
									// these two transcripts both have one good start and one different --> keep them both (different option would be to add them to another srfrag that is compatible and has both hardends but then it's more complicated
									break;
								}

								// I keep both if both are guides
								//if((!srfrag[t[0]]->guide || !srfrag[t[1]]->guide) && abs(len[0])<edgedist && abs(len[2])<edgedist) { // close by
								if((!srfrag[t[0]]->guide || !srfrag[t[1]]->guide) && len[0]<edgedist && len[2]<edgedist) { // close by
									if(srfrag[t[0]]->guide || (!srfrag[t[1]]->guide && no2gnode[srfrag[t[0]]->nodes[0]]->hardstart && no2gnode[srfrag[t[0]]->nodes.Last()]->hardend))
										keeptrf[t2].t=t1; // t[0] to replace t[1]
									keeptrf[t2].cov+=srfrag[t1]->abundance;
									keeptrf[t2].group.Add(t1);
									included=true;
									//fprintf(stderr,"trf %d %d equivalent start/ends\n",t[1],t[0]);
								}
								break;
							}
							if(included) break; // break from for loop
						}
					}

					if(!included){
						if(srfrag[t1]->guide || ((srfrag[t1]->longstart || no2gnode[srfrag[t1]->nodes[0]]->hardstart)&&
								(srfrag[t1]->longend || no2gnode[srfrag[t1]->nodes.Last()]->hardend))) { // if this is not included and has correct start/end
							CLongTrf kt(t1,srfrag[t1]->abundance);
							keeptrf.Add(kt);
							keeptrf.Last().group.Add(t1);
							//fprintf(stderr,"keep srfrag %d\n",t1);
						}
						else { // incomplete transcript, possibly wrong
							srfrag[t1]->weak=1;
							//fprintf(stderr,"Incomplete transcript %d\n",t1);
						}
					}
				}
			}


			//GBitVec guidesource(gno);
			//GBitVec guidesink(gno);
			//for(int i=0;i<keeptrf.Count();i++) {
			for(int i=keeptrf.Count()-1;i>=0;i--) { // I add the kept transcripts to trflong from least significant to most in order to make deletion easier

				//fprintf(stderr,"Build source/sink for srfrag %d\n",keeptrf[i].t);
				int n1=srfrag[keeptrf[i].t]->nodes[0];
				int n2=srfrag[keeptrf[i].t]->nodes.Last();
				//if(hassource[n1]<0 || hassink[n2]<0) fprintf(stderr,"Build source/sink for srfrag %d\n",keeptrf[i].t);

				if(!srfrag[keeptrf[i].t]->guide && ((!no2gnode[n1]->hardstart && (hassource[n1]<0 || !keepsource[n1])) ||
						(!no2gnode[n2]->hardend && (hassink[n2]<0 || !keepsink[n2])))) {
				//if((no2gnode[n1]->hardstart || (hassource[n1]>=0 && keepsource[n1])) && (no2gnode[n2]->hardend ||(hassink[n2]>=0 && keepsink[n2]))) {
					srfrag[keeptrf[i].t]->usepath=-2-ntrf; // order of transfrag
					//fprintf(stderr,"keeptrf[%d].t=%d srfrag.usepath=%d\n",i,keeptrf[i].t,t);
					ntrf++;
					//trflong.Add(keeptrf[i].t);
				}

				/******* previous implementation here

				addsource[n1]=keeptrf[i].cov;
				addsink[n2]=keeptrf[i].cov;
				 *******/

				/*
				if(!addsource[n1] && hassource[n1]<0) {
					int startpos=no2gnode[n1]->start-refstart;
					if(startpos-CHI_THR<0 || startpos+CHI_THR>bpcov->Count()) addsource[n1]=1;
					else {
						addsource[n1]=(get_cov(1,startpos,startpos+CHI_THR-1,bpcov)-get_cov(2-2*s,startpos,startpos+CHI_THR-1,bpcov)-
								get_cov(1,startpos-CHI_THR,startpos-1,bpcov)+get_cov(2-2*s,startpos-CHI_THR,startpos-1,bpcov))/(DROP*CHI_THR);
					}
				}

				if(!addsink[n2]  && hassink[n1]<0) {
					int endpos=no2gnode[n2]->end-refstart;
					if(endpos-CHI_THR<0 || endpos+CHI_THR>bpcov->Count()-1) addsink[n2]=1;
					else {
						addsink[n2]=(get_cov(1,endpos-CHI_THR+1,endpos,bpcov)-get_cov(2-2*s,endpos-CHI_THR+1,endpos,bpcov)-
								get_cov(1,endpos+1,endpos+CHI_THR,bpcov)+get_cov(2-2*s,endpos+1,endpos+CHI_THR,bpcov))/(DROP*CHI_THR);
					}
				}*/

				// all
				for(int j=0;j<keeptrf[i].group.Count();j++) {
					if(n1==srfrag[keeptrf[i].group[j]]->nodes[0]) {
						addsource[n1]+=srfrag[keeptrf[i].group[j]]->abundance;
						//fprintf(stderr,"Add source t[%d]->cov=%f to node %d = %f\n",keeptrf[i].group[j],srfrag[keeptrf[i].group[j]]->abundance,n1,addsource[n1]);
					}
					if(n2==srfrag[keeptrf[i].group[j]]->nodes.Last()) {
						addsink[n2]+=srfrag[keeptrf[i].group[j]]->abundance;
						//fprintf(stderr,"Add sink t[%d]->cov=%f to node %d = %f\n",keeptrf[i].group[j],srfrag[keeptrf[i].group[j]]->abundance,n2,addsink[n2]);
						}
					}

			}
			for(int i=keeptrf.Count()-1;i>=0;i--) {
				int n1=srfrag[keeptrf[i].t]->nodes[0];
				int n2=srfrag[keeptrf[i].t]->nodes.Last();

				if(srfrag[keeptrf[i].t]->guide || ((no2gnode[n1]->hardstart || (hassource[n1]>=0 && keepsource[n1])) &&
						(no2gnode[n2]->hardend ||(hassink[n2]>=0 && keepsink[n2])))) {
					//trflong.Add(keeptrf[i].t);
					srfrag[keeptrf[i].t]->usepath=-2-ntrf; // mark transfrag that it needs to be part of trflong
					ntrf++;
				}
			}

			/*
			{ // DEBUG ONLY
				fprintf(stderr,"%d keeptrf:\n",keeptrf.Count());
				for(int i=0;i<keeptrf.Count();i++) {
					fprintf(stderr,"(%d) %d abund=%f keepcov=%f",i,keeptrf[i].t,srfrag[keeptrf[i].t]->abundance,keeptrf[i].cov);
					for(int j=0;j<srfrag[keeptrf[i].t]->nodes.Count();j++) {
						fprintf(stderr," %d",srfrag[keeptrf[i].t]->nodes[j]);
					}
					fprintf(stderr,"\n");

				}
				fprintf(stderr,"%d trflong:",trflong.Count());
				for(int i=0;i<trflong.Count();i++)
					fprintf(stderr," %d",trflong[i]);
				fprintf(stderr,"\n");
			}
			*/


			// add source/sink connections
			for(int i=1;i<gno-1;i++) {
				//fprintf(stderr,"i=%d hassource=%d addsource=%f hassink=%d addsink=%f hardstart=%d hardend=%d\n",i,hassource[i],addsource[i],hassink[i],addsink[i],no2gnode[i]->hardstart,no2gnode[i]->hardend);
				if(hassource[i]<0) {
					if(addsource[i] && no2gnode[i]->hardstart) {
						/*float abund=trthr;
						if(no2gnode[i]->hardstart) { // node i doesn't have source as parent but it should
							abund=addsource[i];
						}*/
						// else I am not that confident this is a start
						no2gnode[i]->parent.Insert(0,zero);
						no2gnode[0]->child.Add(i);
						GVec<int> nodes;
						nodes.Add(source);
						nodes.Add(i);
						//CTransfrag *t=new CTransfrag(nodes,allpat,abund);
						CTransfrag *t=new CTransfrag(nodes,allpat,addsource[i]);
						t->pattern[source]=1;
						t->pattern[i]=1;
						t->longread=true; /// this is only true for longreads;
						transfrag.Add(t);
						//fprintf(stderr,"Add link to source: 0-%d with abundance=%f\n",i, addsource[i]);
					}
				}
				else {
					/*if(keepsource[i]) {
						srfrag[hassource[i]]->abundance+=addsource[i];
					}*/
					if(!keepsource[i]) { // it was a mistake to link this to source --> remove
						/*fprintf(stderr,"delete abundance of trf:");
						for(int j=0;j<srfrag[hassource[i]]->nodes.Count();j++) fprintf(stderr," %d",srfrag[hassource[i]]->nodes[j]);fprintf(stderr,"\n");*/
						srfrag[hassource[i]]->abundance=0;
					}
					else { // this one has source and it should be kept => update trf abundance to a more realistic value
						//fprintf(stderr,"source %d:%d abund=%f addsource=%f\n",i,no2gnode[i]->start,srfrag[hassource[i]]->abundance,addsource[i]);
						//if(guidesource[i] && srfrag[hassource[i]]->abundance<addsource[i])
						srfrag[hassource[i]]->abundance=addsource[i];
					}
				}

				if(hassink[i]<0) {
					if(addsink[i] && no2gnode[i]->hardend) {
						/*float abund=trthr;
						if(no2gnode[i]->hardend) {
							abund=addsink[i];
						}*/
						// else  not very confident about this end

						no2gnode[i]->child.Add(sink);
						no2gnode[sink]->parent.Add(i);
						GVec<int> nodes;
						nodes.Add(i);
						nodes.Add(sink);
						CTransfrag *t=new CTransfrag(nodes,allpat,addsink[i]);
						t->pattern[sink]=1;
						t->pattern[i]=1;
						t->longread=true; //// this is only true for longreads
						transfrag.Add(t);
						//fprintf(stderr,"Add link to sink: %d-%d with abundance=%f\n",i,sink,addsink[i]);
					}
				}
				else {
					if(!keepsink[i]) { // it was a mistake to link this to sink --> remove
						/*fprintf(stderr,"delete abundance of trf:");
						for(int j=0;j<srfrag[hassink[i]]->nodes.Count();j++) fprintf(stderr," %d",srfrag[hassink[i]]->nodes[j]);fprintf(stderr,"\n");*/
						srfrag[hassink[i]]->abundance=0;
					}
					else {
						//fprintf(stderr,"sink %d:%d abund=%f addsink=%f\n",i,no2gnode[i]->end,srfrag[hassink[i]]->abundance,addsink[i]);
						//if(guidesink[i] && srfrag[hassink[i]]->abundance<addsink[i])
						srfrag[hassink[i]]->abundance=addsink[i];
					}
					/*if(keepsink[i]) {
						//fprintf(stderr,"Add %f to longtransfrag=%d\n",addsink[i],hassink[i]);
						srfrag[hassink[i]]->abundance+=addsink[i];
					}*/
				}
			}
			//fprintf(stderr,"There should be %d trflong transcripts\n",ntrf);
			trflong.Resize(ntrf,-1);
	}


	// add edges between disconnected parent-child nodes
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
				if(longreads) t->longread=true;
				transfrag.Add(t);
			}
		}
	}

	// sort transfrag with smallest being the one that has the most nodes, and ties are decided by the abundance (largest abundance first); last transfrags all have 1 node
	if(trsort)
		transfrag.Sort(trCmp);

	/*
	{ // DEBUG ONLY
		printTime(stderr);
		fprintf(stderr,"There are %d transfrags that remained\n",transfrag.Count());
	}
	*/

	/*
	{ // DEBUG ONLY
		//printTime(stderr);
		fprintf(stderr,"\nThere are %d transfrags after sorting:\n",transfrag.Count());
		for(int i=0;i<transfrag.Count();i++) {
			fprintf(stderr,"transfrag[%d]:",i);
			for(int j=0;j<transfrag[i]->nodes.Count();j++) fprintf(stderr," %d",transfrag[i]->nodes[j]);
			fprintf(stderr," abundance=%f usepath=%.1f\n",transfrag[i]->abundance,transfrag[i]->usepath);
		}
	}
	*/


	GVec<int> incompletetrf; //remembers incomplete transfrags (the ones that don't have edges between two consecutive nodes


	// create compatibilities
	for(int t1=0;t1<transfrag.Count();t1++) { // transfrags are processed in increasing order -> important for the later considerations

		// update nodes
		int n1=transfrag[t1]->nodes.Count();

		if(mixedMode) {
			if(transfrag[t1]->usepath<-1) { // this is a mixedMode long transfrag (TODO: consider using this for longreads mode also)

				//fprintf(stderr,"Add transcript %d to trflong with value %d\n",t1,(int)transfrag[t1]->usepath);
				trflong[abs(2+(int)transfrag[t1]->usepath)]=t1;
			}
			//transfrag[t1]->usepath=-1; // restore order
			transfrag[t1]->usepath=transfrag[t1]->abundance; // store abundance for later
		}

		if(n1>1) { // add transfrag to nodes' in and out; if a transfrag only has one node then it is not added to a node; I might want to change this for the computation of fpkm's
			bool incomplete = false;
			bool nosplice=true; // try to give less priority to unspliced reads vs spliced reads
			for(int n=0;n<n1;n++) { // for all nodes in transfrag

				if(nosplice && n) { // reduce abundance of continuous transfrags
					if(transfrag[t1]->nodes[n]!=1+transfrag[t1]->nodes[n-1] || no2gnode[transfrag[t1]->nodes[n]]->start-1!=no2gnode[transfrag[t1]->nodes[n-1]]->end) {
						nosplice=false;
					}
				}

				if(n && n<transfrag[t1]->nodes.Count()-1) {// not first or last node
					// add t1 to in and out of node
					no2gnode[transfrag[t1]->nodes[n]]->trf.Add(t1);

					if(transfrag[t1]->nodes[n-1] && transfrag[t1]->nodes[n]<gno-1) {
						// check if transfrag t1 is incomplete between node[n-1] and node [n]
						int *pos=gpos[edge(transfrag[t1]->nodes[n-1],transfrag[t1]->nodes[n],gno)];
						if(!pos || !transfrag[t1]->pattern[*pos]) // there is no edge between node[n-1] and node[n]
							incomplete = assign_incomplete_trf_to_nodes(t1,transfrag[t1]->nodes[n-1],transfrag[t1]->nodes[n],no2gnode) or incomplete; 	// this ensures that I still have compatibilities when going through a certain path: !!! THIS IS NEW COMPARED TO PERL
					}
				}
				else if(n) { // last but not first node
					// add t1 to in of node
					no2gnode[transfrag[t1]->nodes[n]]->trf.Add(t1);

					if(transfrag[t1]->nodes[n-1] && transfrag[t1]->nodes[n]<gno-1) {
						// check if transfrag t1 is incomplete between node[n-1] and node [n]
						int *pos=gpos[edge(transfrag[t1]->nodes[n-1],transfrag[t1]->nodes[n],gno)];
						if(!pos || !transfrag[t1]->pattern[*pos]) // there is no edge between node[n-1] and node[n]
							incomplete = assign_incomplete_trf_to_nodes(t1,transfrag[t1]->nodes[n-1],transfrag[t1]->nodes[n],no2gnode) or incomplete; 	// this ensures that I still have compatibilities when going through a certain path: !!! THIS IS NEW COMPARED TO PERL
					}
				}
				else { // first node -> only add transfrag to out of node
					no2gnode[transfrag[t1]->nodes[n]]->trf.Add(t1);
				}
			}

			//if(nosplice) transfrag[t1]->abundance*=(1-isofrac);
			//transfrag[t1]->abundance*=0.5;


			if(incomplete) incompletetrf.Add(t1);
			else transfrag[t1]->real=true;

		}
		//else if(longreads) no2gnode[n1]->trf.Add(t1);
		/*
		else { // this transcript is included completely in node
			no2gnode[transfrag[t1]->nodes[0]]->frag+=transfrag[t1]->abundance;
		}
		*/
	} // end for(int t1=0;t1<transfrag.Count();t1++)


	if(srfrag.Count() && !mixedMode) {
		srfrag.Sort(trCmp); // always start with largest super-read to solve
		for(int u=0;u<srfrag.Count();u++) //process_srfrag(srfrag[u],transfrag,no2gnode,gno,gpos);
		  if(!srfrag[u]->abundance) srfrag[u]->abundance=srfrag[u]->srabund*ERROR_PERC;
	}


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
			if(t0>-1 && transfrag[t0]->abundance) { // found transfrag from source to node and the transfrag wasn't deleted
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
	GStr name("",128);

	// first compute coverage
	float cov=0;
	for(int i=0;i<alltr.Count();i++) {
		int a=alltr[i];
		cov+=mgt[alltr[i]]->transfrag->abundance;
		mgt[alltr[i]]->transfrag->abundance=0;
		if(enableNames || (alltr.Count()==1 && mgt[alltr[0]]->transfrag->real)) for(int j=0;j<mgt[a]->read.Count();j++) {
			int r=mgt[a]->read[j];
			if(readlist[r]->tinfo->g == -1) {
				int fidx=1+readlist[r]->tinfo->fileidx;
				//GStr fid(fidx);
				if(!name.is_empty()) {
					name+=", ";
				}
				name+=fidx;
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
	GStr name("", 128);

	if(enableNames) for(int i=0;i<alltr.Count();i++) {
		int a=alltr[i];
		for(int j=0;j<mgt[a]->read.Count();j++) {
			int r=mgt[a]->read[j];
			if(!name.is_empty()) {
				name+=", ";
			}
			if(readlist[r]->tinfo->g == -1) {
				int fidx=1+readlist[r]->tinfo->fileidx;
				//GStr fid(fidx);
				name+=fidx;
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


bool bfs(int n,GVec<float> *capacity,GVec<float> *flow,GVec<int> *link,GVec<int>& pred) {
	GVec<int> color;
	color.Resize(n+2);
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

		// add transfrag to node
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

				if(remove) {
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
			delete localpred[i].p;
			localpred[i].p=NULL;
		}
	}

	return(geneno);
}



bool is_compatible(int t1,int t2, int n,GBitVec& compatible) {
	if(t1<t2) return(compatible[comptbl_pos(t1,t2,n)]);
	else return(compatible[comptbl_pos(t2,t1,n)]);
}

GVec<int> *max_compon_size_with_penalty(int trnumber,float &maxsize,GVec<CTrInfo>& set,GBitVec& compatible,
		GBitVec& mark,GBitVec& removable, GHash<CComponent*>& computed) {

	// this max_compon presumes the set is always sorted according to the set.trno

	GVec<int> *result=NULL;
	float penalty=0;

	for(int i=0;i<set.Count();i++) {
		float size=set[i].abundance-penalty;
		float maxagreesize=0;
		GVec<CTrInfo> agreeset;
		GStr s("",64);
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

bool onpath_long(GBitVec& trpattern,GVec<int>& trnode,GBitVec& pathpattern,int minp,int maxp,GPVec<CGraphnode>& no2gnode,int gno,
		GIntHash<int>& gpos) {

	int j=0;
	int *edgep=NULL;
	int tn=trnode.Count();
	edgep=gpos[edge(0,minp,gno)];
	if(edgep && pathpattern[*edgep]) { // minp links to source
		if(!trnode[0] && trnode[1]!=minp) return false;
		else { // transfrag does not link to source
			if(trnode[0]<minp) return false;
		}
		edgep=NULL;
	}
	edgep=gpos[edge(maxp,gno-1,gno)];
	if(edgep && pathpattern[*edgep]) { // maxp links to sink
		if(trnode.Last()==gno-1 && trnode[tn-2]!=maxp) return false;
		else { // transfrag does not link to sink
			if(trnode.Last()>maxp) return false;
		}
		edgep=NULL;
	}
	int prevp=-1;
	int p=minp;
	while(1) {
		while(j<tn && trnode[j]<p) {
			if(!prevp || (edgep && pathpattern[*edgep]) || p==gno-1) return false; // there is an edge between prevp<trnode[j] and p>trnode[j]; when prevp==0/gno-1 I might not have an edge there
			if(!no2gnode[trnode[j]]->childpat[p]) {
			        //fprintf(stderr,"Node p=%d cannot be reached from node=%d\n",p,trnode[j]);
				return false; // I can not reach p from trnode[j]
			}
			j++;
		}
		if(j==tn) return true; // went through transcript
		// now trnode[j]>=p and trnode[j-1]<p
		if(trnode[j]>p) {
			if(!no2gnode[p]->childpat[trnode[j]]) return false; // I can not reach trnode[j] from p
			int *edget=NULL;
			if(j) edget=gpos[edge(trnode[j-1],trnode[j],gno)];
			if(edget && trpattern[*edget]) return false; // there is an edge between trnode[j-1]<p and trnode[j]>p
		}
		else j++; // now trnode[j]>p and trnode[j-1]==prevp

		if(p==maxp) return true; // went through full path
		prevp=p;
		p=pathpattern.find_next(prevp); // p is not -1 since prevp<maxp
		edgep=gpos[edge(prevp,p,gno)];
	}

	return true;
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

	if(pattern[path.Last()] && pathincov.Last()-abundance<epsilon) {
		penalty=abundance;
		return(false); // transcript entering the last added node to path is the only one
	}

	if(path.Last()!=firstnode) { // transcript doesn't start at last node added to path -> this shouldn't happen since the transcript has already been added
		for(int i=path.Count()-2;i>=0;i--) {
			if(pattern[path[i]] && path[i]==firstnode) { // transcript starts at this node
				penalty=abundance;
				if(pathoutcov[i]-abundance<epsilon) {
					return(false);
				}
				break;
			}
			else // middle of the path
				if(pattern[path[i]]) {
					penalty=abundance;
					if((pathincov[i]-abundance<epsilon) || (pathoutcov[i]-abundance<epsilon)) {
						return(false);
					}
				}
		}
	}

  return(true);
}

void compute_weak(GPVec<CTransfrag>& transfrag,int t,GVec<float>& nodecov,GPVec<CGraphnode>& no2gnode) {
	transfrag[t]->weak=0;
	int n=1;
	while(n<transfrag[t]->nodes.Count()) {
		//if(inode->child[c]==i+1 && i<gno-2 && inode->end+1==cnode->start &&
		//					nodecov[i+1]/cnode->len() <1000 && nodecov[i]*(DROP+ERROR_PERC)>nodecov[i+1])  { // adjacent to child
		if(transfrag[t]->nodes[n]==transfrag[t]->nodes[n-1]+1 &&
				no2gnode[transfrag[t]->nodes[n-1]]->end+1==no2gnode[transfrag[t]->nodes[n]]->start &&
				(nodecov[transfrag[t]->nodes[n-1]]*(DROP+ERROR_PERC)>nodecov[transfrag[t]->nodes[n]] ||
						nodecov[transfrag[t]->nodes[n]]*(DROP+ERROR_PERC)>nodecov[transfrag[t]->nodes[n-1]])) {
			transfrag[t]->weak=1;
			break;
		}
		n++;
	}
}

void replace_transfrag(int t,int& tmax,GPVec<CTransfrag>& transfrag,GVec<float>& nodecov,GPVec<CGraphnode>& no2gnode) {

	if(tmax==-1) {
		tmax=t;
		if(transfrag[tmax]->weak<0) {
			compute_weak(transfrag,tmax,nodecov,no2gnode);
		}
	}
	else {
		if(transfrag[t]->weak<0) {
			compute_weak(transfrag,t,nodecov,no2gnode);
		}

		if(!transfrag[tmax]->weak) {
			if(!transfrag[t]->weak && transfrag[t]->abundance>transfrag[tmax]->abundance) {
				tmax=t;
			}
		}
		else { // tmax is weak
			if(!transfrag[t]->weak || transfrag[t]->abundance>transfrag[tmax]->abundance) tmax=t;
		}

	}

}

bool fwd_to_sink_fast_long(int i,GVec<int>& path,int& minpath,int& maxpath,GBitVec& pathpat,GPVec<CTransfrag>& transfrag,GPVec<CGraphnode>& no2gnode,
		GVec<float>& nodecov,int gno,GIntHash<int>& gpos){
	// find all parents -> if parent is source then go back
	CGraphnode *inode=no2gnode[i];

	if(i<maxpath && !inode->childpat[maxpath]) return false; // I can not reach maxpath from node

	int nchildren=inode->child.Count(); // number of children

	/*
	fprintf(stderr,"Children of node %d with maxpath=%d are:",i,maxpath);
	for(int c=0;c<nchildren;c++) fprintf(stderr," %d",inode->child[c]);
	fprintf(stderr," pathpat=");
	printBitVec(pathpat);
	fprintf(stderr,"\n");
	*/

	if(!nchildren) return true; // node is sink
	int maxc=-1;

	/* I can not do this because then I don't check for continuity
	if(nchildren==1) {
		maxc=inode->child[0]; // only one child to consider
	}
	else {
	*/
	float maxcov=0;
	int tmax=-1;
	bool exclude=false;

	int nextnode=0; // nextnode on path
	bool reach=false;
	if(maxpath<=i) reach=true;

	//int maxchild=inode->child[0];
	//float maxchildcov=-1;
	if(pathpat[i+1]) {
		maxc=i+1;
		tmax=-1;
		reach=true;
	}
	else {
		for(int c=0;c<nchildren;c++) {

			bool childonpath=false;
			if(pathpat[inode->child[c]]) childonpath=true;
			// check if child is already on the path
			int *pos=gpos[edge(i,inode->child[c],gno)];
			if(pos && pathpat[*pos]) {
				maxc=inode->child[c];
				tmax=-1;
				reach=true;
				break;
			}

			// check if I can reach nextnode from c
			if(maxpath>i) { // actually maxpath is >= i+2 otherwise I wouldn't get here
				if(!nextnode) {
					int j=i+2;
					while(j<=maxpath) {
						if(pathpat[j]) {
							nextnode=j;
							break;
						}
						j++;
					}
				}
				if(inode->child[c]!=nextnode && !no2gnode[inode->child[c]]->childpat[nextnode]) {
					continue;
				}
				reach=true;
			}

			/*if(no2gnode[maxpath]->hardend && pathpat[gno-1]) { // this is a trflong transcript
				if(!pathpat[c]) continue;
			}*/

			float childcov=0;
			int tchild=-1;
			int endpath=maxpath;
			if(inode->child[c]>endpath) endpath=inode->child[c];
			CGraphnode *cnode=no2gnode[inode->child[c]];
			/*
				if(nodecov[inode->child[c]]>maxchildcov) {
					maxchildcov=nodecov[inode->child[c]];
					maxchild=inode->child[c];
				}
			 */

			if(inode->child[c]==i+1 && i<gno-2 && inode->end+1==cnode->start &&
					nodecov[i+1]/cnode->len() <1000 && nodecov[i]*(DROP+ERROR_PERC)>nodecov[i+1])  { // adjacent to child
				// cnode->len() this is redundant because i<gno-2
				//((nodecov[i+1]/cnode->len() <1000 && nodecov[i]*DROP>nodecov[i+1]) || (cnode->len()<longintronanchor && cnode->child.Count()==1 && cnode->child[0]==gno-1 ) ))  { // adjacent to child
				exclude=true;
			}
			else {
				pathpat[inode->child[c]]=1;
				if(pos) pathpat[*pos]=1;
				//else GError("Found parent-child: %d-%d not linked by edge!\n",i,inode->child[c]);
				for(int j=0;j<cnode->trf.Count();j++) { // for all transfrags going through child
					int t=cnode->trf[j];
					if(transfrag[t]->abundance<epsilon) { // this transfrag was used before -> needs to be deleted
						if(!mixedMode) { cnode->trf.Delete(j); j--;}
						else transfrag[t]->abundance=0;
					}
					else if(transfrag[t]->longread && transfrag[t]->nodes[0]) {
						if(inode->child[c]==gno-1) { // child is sink
							if(transfrag[t]->nodes[0]==i && maxpath<=i) { // need to check it doesn't violate path
								childcov+=transfrag[t]->abundance;
								if(tchild==-1 || transfrag[t]->abundance>transfrag[tchild]->abundance) tchild=t;
							}
						}
						else if(transfrag[t]->nodes[0]<=i && transfrag[t]->nodes.Last()>=inode->child[c] &&    // transfrag goes from i to c
								//(transfrag[t]->pattern[inode->child[c]] || transfrag[t]->pattern[i]) &&    // transfrag is not incomplete through these nodes
								onpath_long(transfrag[t]->pattern,transfrag[t]->nodes,pathpat,minpath,endpath,no2gnode,gno,gpos)) { // transfrag is compatible with path
							//fprintf(stderr,"add transfr[%d]->abund=%f\n",t,transfrag[t]->abundance);
							childcov+=transfrag[t]->abundance;
							replace_transfrag(t,tchild,transfrag,nodecov,no2gnode);
							//if(tchild==-1 || transfrag[t]->abundance>transfrag[tchild]->abundance) tchild=t;
						}
					}
				}

				if(childcov>maxcov) {
					if(tmax==-1 || transfrag[tmax]->weak>0 || !transfrag[tchild]->weak) {
						maxcov=childcov;
						maxc=inode->child[c];
						tmax=tchild;
					}
				}
				else if(maxc!=-1 && childcov>maxcov-epsilon) {
					if(transfrag[tmax]->weak>0 || !transfrag[tchild]->weak) {
						if(nodecov[maxc]<nodecov[inode->child[c]]) {
							maxc=inode->child[c];
							tmax=tchild;
						}
					}
				}

				if(pos) pathpat[*pos]=0;
				if(childonpath) break;
				pathpat[inode->child[c]]=0;


			}
		}
	}

	if(!reach) return false;

	if(maxc==-1) {
		if(exclude && nodecov[i+1]) {
			CGraphnode *cnode=no2gnode[i+1];
			tmax=-1;
			float childcov=0;
			pathpat[i+1]=1;
			int *pos=gpos[edge(i,i+1,gno)];
			if(pos) pathpat[*pos]=1;
			int endpath=maxpath;
			if(i+1>endpath) endpath=i+1;
			//else GError("Found parent-child: %d-%d not linked by edge\n",i,i+1);
			for(int j=0;j<cnode->trf.Count();j++) { // for all transfrags going through child
				int t=cnode->trf[j];
				if(transfrag[t]->abundance<epsilon) { // this transfrag was used before -> needs to be deleted
					if(!mixedMode) { cnode->trf.Delete(j); j--;}
					else transfrag[t]->abundance=0;
				}
				else if(transfrag[t]->longread && transfrag[t]->nodes[0]<=i && transfrag[t]->nodes.Last()>=i+1 &&   // transfrag goes from i to c
						//(transfrag[t]->pattern[i+1] || transfrag[t]->pattern[i]) &&    // transfrag is not incomplete through these nodes
						onpath_long(transfrag[t]->pattern,transfrag[t]->nodes,pathpat,minpath,endpath,no2gnode,gno,gpos)) { // transfrag is compatible with path
					childcov+=transfrag[t]->abundance;
					//if(tmax==-1 || transfrag[t]->abundance>transfrag[tmax]->abundance) tmax=t;
					replace_transfrag(t,tmax,transfrag,nodecov,no2gnode);
				}
			}
			pathpat[i+1]=0; // these were 1 in the previous code, but it doesn't make sense -> I hope it was wrong
			if(pos) pathpat[*pos]=0;

			if(childcov) {
				maxc=i+1;
			}
			else {
				//goback=true;
				return false;
			}
		}
		else {
			//goback=true;
			return false; //maxc=maxchild;
		}

	}

	/*
	{ // DEBUG ONLY
		//printTime(stderr);
		fprintf(stderr,"maxc=%d maxcov=%f tmax=%d\n",maxc,maxcov,tmax);
	}
	*/

	// add maxp to path
	path.Add(maxc);
	pathpat[maxc]=1;

	int *pos=gpos[edge(i,maxc,gno)];
	if(pos) pathpat[*pos]=1;
	//else GError("Found parent-child %d-%d not linked by edge\n",i,maxc);
	if(tmax>-1) {
		pathpat=pathpat | transfrag[tmax]->pattern;
		if(transfrag[tmax]->nodes[0]<minpath) minpath=transfrag[tmax]->nodes[0];
		if(transfrag[tmax]->nodes.Last()>maxpath) maxpath=transfrag[tmax]->nodes.Last();
	}

	return fwd_to_sink_fast_long(maxc,path,minpath,maxpath,pathpat,transfrag,no2gnode,nodecov,gno,gpos);
}

bool back_to_source_fast_long(int i,GVec<int>& path,int& minpath,int& maxpath,GBitVec& pathpat,GPVec<CTransfrag>& transfrag,GPVec<CGraphnode>& no2gnode,
		GVec<float>& nodecov,int gno,GIntHash<int>& gpos){

	// find all parents -> if parent is source then go back
	CGraphnode *inode=no2gnode[i];

	if(minpath<i && !inode->parentpat[minpath]) return false; // I can not reach maxpath from node

	int nparents=inode->parent.Count(); // number of parents

	/*
	fprintf(stderr,"Parents of node %d are:",i);
	for(int p=0;p<nparents;p++) fprintf(stderr," %d",inode->parent[p]);
	fprintf(stderr," pathpat=");
	printBitVec(pathpat);
	fprintf(stderr,"\n");
	*/

	if(!nparents) return true; // node is source
	int maxp=-1;


	/* I can not do this because then I don't check for continuity
	if(nparents==1) {
		maxp=inode->parent[0];
	}
	else {
	*/
	float maxcov=0;
	int tmax=-1;
	bool exclude=false;

	int nextnode=gno; // nextnode on path
	bool reach=false;
	if(minpath>=i) reach=true;

	// check if adjacent parent is on the path
	if(pathpat[i-1]) {
		maxp=i-1;
		tmax=-1;
		reach=true;
	}
	else {
		//int maxparent=inode->parent[0];
		//float maxparentcov=-1;
		for(int p=0;p<nparents;p++) {

			bool parentonpath=false;
			if(pathpat[inode->parent[p]]) parentonpath=true;

			// check if parent is already on the path
			int *pos=gpos[edge(inode->parent[p],i,gno)];
			if(pos && pathpat[*pos]) { // this is next child on the path
				maxp=inode->parent[p];
				tmax=-1;
				reach=true;
				break;
			}

			// check if I can reach nextnode from p
			if(minpath<i) { // actually minpath is <= i-2 otherwise I wouldn't get here
				if(nextnode==gno) {
					int j=i-2;
					while(j>=minpath) {
						if(pathpat[j]) {
							nextnode=j;
							break;
						}
						j++;
					}
				}
				if(inode->parent[p]!=nextnode && !no2gnode[inode->parent[p]]->parentpat[nextnode]) {
					continue;
				}
				reach=true;
			}

			/*if(no2gnode[minpath]->hardstart && pathpat[gno-1]) { // this is a trflong transcript
				if(!pathpat[p]) continue;
			}*/

			float parentcov=0;
			int tpar=-1;
			int startpath=minpath;
			if(inode->parent[p]<startpath) startpath=inode->parent[p];
			CGraphnode *pnode=no2gnode[inode->parent[p]];

			/*
				if(nodecov[inode->parent[p]]>maxparentcov) {
					maxparentcov=nodecov[inode->parent[p]];
					maxparent=inode->parent[p];
				}
			 */

			if(inode->parent[p]==i-1 && i>1 && inode->start==pnode->end+1 &&
					nodecov[i-1]/pnode->len() <1000 && nodecov[i]*(DROP+ERROR_PERC)>nodecov[i-1]) { // adjacent to parent
				exclude=true;
			}
			else {
				pathpat[inode->parent[p]]=1;
				if(pos) pathpat[*pos]=1;
				//else GError("Found parent-child %d-%d not linked by edge\n",inode->parent[p],i);

				for(int j=0;j<pnode->trf.Count();j++) { // for all transfrags going through parent
					int t=pnode->trf[j];
					if(transfrag[t]->abundance<epsilon) { // this transfrag was used before -> needs to be deleted
						if(!mixedMode) { pnode->trf.Delete(j);j--;}
						else transfrag[t]->abundance=0;

					}
					else if(transfrag[t]->longread && transfrag[t]->nodes.Last()!=gno-1) {
						if(!inode->parent[p]) { // parent is source
							if(transfrag[t]->nodes.Last()==i && minpath>=i) {
								parentcov+=transfrag[t]->abundance;
								if(tpar==-1 || transfrag[t]->abundance>transfrag[tpar]->abundance) tpar=t;
							}
						}
						else if(transfrag[t]->nodes[0]<=inode->parent[p] && transfrag[t]->nodes.Last()>=i &&   // transfrag goes from p to i
								//(transfrag[t]->pattern[inode->parent[p]]||transfrag[t]->pattern[i]) &&    // transfrag is not incomplete through these nodes
								onpath_long(transfrag[t]->pattern,transfrag[t]->nodes,pathpat,startpath,maxpath,no2gnode,gno,gpos)) { // transfrag is compatible with path
							// comment: I need to check for the parent to make sure I am not going through another node!!
							//fprintf(stderr,"parent=%d add transfr[%d]->abund=%f\n",inode->parent[p],t,transfrag[t]->abundance);
							parentcov+=transfrag[t]->abundance;
							//if(tpar==-1 || transfrag[t]->abundance>transfrag[tpar]->abundance) tpar=t;
							replace_transfrag(t,tpar,transfrag,nodecov,no2gnode);
						}
					}
				}

				if(parentcov>maxcov) {
					if(tmax==-1 || transfrag[tmax]->weak>0 || !transfrag[tpar]->weak) {
						maxcov=parentcov;
						maxp=inode->parent[p];
						tmax=tpar;
					}
				}
				else if(maxp!=-1 && parentcov>maxcov-epsilon) {
					if(transfrag[tmax]->weak>0 || !transfrag[tpar]->weak) {
						if(nodecov[maxp]<nodecov[inode->parent[p]]) {
							maxp=inode->parent[p];
							tmax=tpar;
						}
					}
				}

				if(pos) pathpat[*pos]=0;
				if(parentonpath) break;
				pathpat[inode->parent[p]]=0;

			}
		}
	}

	if(!reach) return false;

	if(maxp==-1) {
		if(exclude && nodecov[i-1]) {
			CGraphnode *pnode=no2gnode[i-1];
			tmax=-1;
			float parentcov=0;
			pathpat[i-1]=1;
			int *pos=gpos[edge(i-1,i,gno)];
			if(pos) pathpat[*pos]=1;
			//else GError("Found parent-child %d-%d not linked by edge\n",i-1,i);
			int startpath=minpath;
			if(i-1<startpath) startpath=i-1;
			for(int j=0;j<pnode->trf.Count();j++) { // for all transfrags going through parent
				int t=pnode->trf[j];
				if(transfrag[t]->abundance<epsilon) { // this transfrag was used before -> needs to be deleted
					if(!mixedMode) { pnode->trf.Delete(j); j--;}
					else transfrag[t]->abundance=0;
				}
				else if(transfrag[t]->longread && transfrag[t]->nodes[0]<=i-1 && transfrag[t]->nodes.Last()>=i &&   // transfrag goes from p to i
						//(transfrag[t]->pattern[i-1]||transfrag[t]->pattern[i]) &&    // transfrag is not incomplete through these nodes
						onpath_long(transfrag[t]->pattern,transfrag[t]->nodes,pathpat,startpath,maxpath,no2gnode,gno,gpos)) { // transfrag is compatible with path
					// comment: I need to check for the parent to make sure I am not going through another node!!
					//fprintf(stderr,"parent=%d add transfr[%d]->abund=%f\n",i-1,t,transfrag[t]->abundance);
					parentcov+=transfrag[t]->abundance;
					//if(tmax==-1 || transfrag[t]->abundance>transfrag[tmax]->abundance) tmax=t;
					replace_transfrag(t,tmax,transfrag,nodecov,no2gnode);
				}
			}
			pathpat[i-1]=0;
			if(pos) pathpat[*pos]=0;

			if(parentcov) {
				maxp=i-1;
			}
			else {
				//goback=true;
				return false;
			}
		}
		else {
			//goback=true;
			return false; //maxp=maxparent;
		}
	}

	/*
	{ // DEBUG ONLY
		//printTime(stderr);
		fprintf(stderr,"maxp=%d maxcov=%f tmax=%d\n",maxp,maxcov,tmax);
	}
	*/

	if(maxp) { // add maxp to path only if not source
		path.Add(maxp);
	}

	pathpat[maxp]=1;                 // if maxp is source I added it in the pathpat
	int *pos=gpos[edge(maxp,i,gno)];
	if(pos) pathpat[*pos]=1;
	//else GError("Found parent-child %d-%d not linked by edge\n",maxp,i);
	if(tmax>-1) {
		pathpat=pathpat | transfrag[tmax]->pattern;
		if(transfrag[tmax]->nodes[0]<minpath) minpath=transfrag[tmax]->nodes[0];
		if(transfrag[tmax]->nodes.Last()>maxpath) maxpath=transfrag[tmax]->nodes.Last();
	}

	return back_to_source_fast_long(maxp,path,minpath,maxpath,pathpat,transfrag,no2gnode,nodecov,gno,gpos);
}

bool fwd_to_sink_fast(int i,GVec<int>& path,GBitVec& pathpat,GPVec<CTransfrag>& transfrag,GPVec<CGraphnode>& no2gnode,
		GVec<float>& nodecov,int gno,GIntHash<int>& gpos){
	// find all parents -> if parent is source then go back
	CGraphnode *inode=no2gnode[i];

	int nchildren=inode->child.Count(); // number of children

	/*
	fprintf(stderr,"Children of node %d are:",i);
	for(int c=0;c<nchildren;c++) fprintf(stderr," %d",inode->child[c]);
	fprintf(stderr,"\n");
	*/

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

		if(inode->child[c]==i+1 && i<gno-2 && inode->end+1==cnode->start &&
				nodecov[i+1]/cnode->len() <1000 && nodecov[i]*(DROP+ERROR_PERC)>nodecov[i+1])  { // adjacent to child
				// cnode->len() this is redundant because i<gno-2
				//((nodecov[i+1]/cnode->len() <1000 && nodecov[i]*DROP>nodecov[i+1]) || (cnode->len()<longintronanchor && cnode->child.Count()==1 && cnode->child[0]==gno-1 ) ))  { // adjacent to child
			exclude=true;
		}
		else {
			pathpat[inode->child[c]]=1;
	    	int *pos=gpos[edge(i,inode->child[c],gno)];
	    	if(pos) pathpat[*pos]=1;
	    	else GError("1 Found parent-child: %d-%d not linked by edge!\n",i,inode->child[c]);
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
			else if(childcov==maxcov && maxc!=-1) {
				if(nodecov[maxc]<nodecov[inode->child[c]]) {
					maxc=inode->child[c];
				}
			}

			pathpat[inode->child[c]]=0;
			if(pos) pathpat[*pos]=0;
		}
	}
	if(maxc==-1) {
		//bool goback=false;
		if(exclude && nodecov[i+1]) {
			CGraphnode *cnode=no2gnode[i+1];
			float childcov=0;
			pathpat[i+1]=1;
	    	int *pos=gpos[edge(i,i+1,gno)];
	    	if(pos) pathpat[*pos]=1;
	    	else GError("2 Found parent-child: %d-%d not linked by edge\n",i,i+1);
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
			else {
				//goback=true;
				//fprintf(stderr,"fwd: return false\n");
				return false;
			}
		}
		else {
			//goback=true;
			//fprintf(stderr,"fwd: return false maxc=maxchild=%d\n",maxc);
			return false; //maxc=maxchild;
		}

		/*if(goback) { // no continuation to source was found
			// check if I have a better path before returning false; might want to restrict to long reads?
			while(path.Count()>1) {
				int c=path.Pop();
				pathpat[c]=0;
				int p=path.Last();
				int *pos=gpos[edge(p,c,gno)];
				if(pos) pathpat[*pos]=0;
				if(no2gnode[p]->child.Last()==gno-1) { // has sink child --> check if there is any abundance left
					for(int j=0;j<no2gnode[p]->trf.Count();j++) { // for all transfrags going through parent
						int t=no2gnode[p]->trf[j];
						if(transfrag[t]->abundance<epsilon) { // this transfrag was used before -> needs to be deleted
							no2gnode[p]->trf.Delete(j);
							j--;
						}
						else if(transfrag[t]->nodes.Last()==gno-1) { // transfrag goes through sink
							pathpat[gno-1]=1;
							int *pos=gpos[edge(p,gno-1,gno)];
							if(pos) pathpat[*pos]=1;
							c=gno-1;
							path.Add(c);
							return true;
						}
					}
				}
			}
			return false;
		}*/

	}

	/*
	{ // DEBUG ONLY
		//printTime(stderr);
		fprintf(stderr,"maxc=%d maxcov=%f\n",maxc,maxcov);
	}
	*/

	// add maxp to path
	path.Add(maxc);
	pathpat[maxc]=1;

	int *pos=gpos[edge(i,maxc,gno)];
	if(pos) pathpat[*pos]=1;
	else GError("3 Found parent-child %d-%d not linked by edge\n",i,maxc);

	return fwd_to_sink_fast(maxc,path,pathpat,transfrag,no2gnode,nodecov,gno,gpos);
}

bool back_to_source_fast(int i,GVec<int>& path,GBitVec& pathpat,GPVec<CTransfrag>& transfrag,GPVec<CGraphnode>& no2gnode,
		GVec<float>& nodecov,int gno,GIntHash<int>& gpos){

	// find all parents -> if parent is source then go back
	CGraphnode *inode=no2gnode[i];

	int nparents=inode->parent.Count(); // number of parents

	/*
	fprintf(stderr,"Parents of node %d are:",i);
	for(int p=0;p<nparents;p++) fprintf(stderr," %d",inode->parent[p]);
	fprintf(stderr,"\n");
	*/

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

		if(inode->parent[p]==i-1 && i>1 && inode->start==pnode->end+1 &&
				nodecov[i-1]/pnode->len() <1000 && nodecov[i]*(DROP+ERROR_PERC)>nodecov[i-1]) { // adjacent to parent
				// pnode->len() && this is redundant because source should always be 0
				//((nodecov[i-1]/pnode->len() <1000 && nodecov[i]*DROP>nodecov[i-1]) || (pnode->len()<longintronanchor && pnode->parent.Count()==1 && !pnode->parent[0])))  { // adjacent to parent
			exclude=true;
		}
		else {
			pathpat[inode->parent[p]]=1;
			int *pos=gpos[edge(inode->parent[p],i,gno)];
			if(pos) pathpat[*pos]=1;
			else GError("4 Found parent-child %d-%d not linked by edge\n",inode->parent[p],i);

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
			else if(parentcov==maxcov && maxp!=-1) {
				if(nodecov[maxp]<nodecov[inode->parent[p]]) {
					maxp=inode->parent[p];
				}
			}

			pathpat[inode->parent[p]]=0;
			if(pos) pathpat[*pos]=0;
		}
	}
	if(maxp==-1) {
		//bool goback=false;
		if(exclude && nodecov[i-1]) {
			CGraphnode *pnode=no2gnode[i-1];
			float parentcov=0;
			pathpat[i-1]=1;
			int *pos=gpos[edge(i-1,i,gno)];
			if(pos) pathpat[*pos]=1;
			else GError("5 Found parent-child %d-%d not linked by edge\n",i-1,i);

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
			else {
				//goback=true;
				//fprintf(stderr,"return false\n");
				return false;
			}
		}
		else {
			//goback=true;
			//fprintf(stderr,"return false maxp=maxparent=%d\n",maxp);
			return false; //maxp=maxparent;
		}

	}

	/*
	{ // DEBUG ONLY
		//printTime(stderr);
		fprintf(stderr,"maxp=%d maxcov=%f\n",maxp,maxcov);
	}
	*/

	if(maxp) { // add maxp to path only if not source
		path.Add(maxp);
	}

	pathpat[maxp]=1;                 // if maxp is source I added it in the pathpat
	int *pos=gpos[edge(maxp,i,gno)];
	if(pos) pathpat[*pos]=1;
	else GError("6 Found parent-child %d-%d not linked by edge\n",maxp,i);

	return back_to_source_fast(maxp,path,pathpat,transfrag,no2gnode,nodecov,gno,gpos);
}


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
	color.Resize(n);
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
		capacity[i].Resize(n);
		flow[i].Resize(n);
	}

	// establish capacities in the network
	for(int i=0;i<n;i++) {
		int nt=no2gnode[path[i]]->trf.Count();
		for(int j=0;j<nt;j++) {
			int t=no2gnode[path[i]]->trf[j];
			if(transfrag[t]->abundance && (istranscript[t] || ((pathpat & transfrag[t]->pattern)==transfrag[t]->pattern))) { // this needs to be adjusted
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

float long_max_flow(int gno,GVec<int>& path,GBitVec& istranscript,GPVec<CTransfrag>& transfrag,GPVec<CGraphnode>& no2gnode,
		GVec<float>& nodecapacity,GBitVec& pathpat) {

	float flux=0;
	int n=path.Count();
	GVec<float> *capacity=new GVec<float>[n]; // capacity of edges in network
	GVec<float> *flow=new GVec<float>[n]; // flow in network
	GVec<int> *link=new GVec<int>[n]; // for each node remembers it's neighbours
	GVec<int> pred; // this stores the augmenting path
	pred.Resize(n,-1);
	GVec<int> node2path;
	GVec<float> noderate;
	node2path.Resize(gno,-1);

	/*
	{ // DEBUG ONLY
		printTime(stderr);
		fprintf(stderr,"Start max flow algorithm for gno=%d and path ",gno);
		printBitVec(pathpat);
		fprintf(stderr," :");
		for(int i=0;i<n;i++) fprintf(stderr," %d:%d",i,path[i]);
		fprintf(stderr,"\n");
		//fprintf(stderr,"Used transcripts:");
		//for(int i=0;i<transfrag.Count();i++) if(istranscript[i]) fprintf(stderr," %d",i);
		//fprintf(stderr,"\n");
	}
	*/

	for(int i=0;i<n;i++) {
		node2path[path[i]]=i;
		nodecapacity.cAdd(0.0);
		noderate.cAdd(1.0); // I set up all rates to be 1 for now
		capacity[i].Resize(n);
		flow[i].Resize(n);
	}

	float max_fl=0;

	// establish capacities in the network
	for(int i=0;i<n;i++) {
		int nt=no2gnode[path[i]]->trf.Count();
		float sumleft=0;
		float sumright=0;
		for(int j=0;j<nt;j++) {
			int t=no2gnode[path[i]]->trf[j];
			if(transfrag[t]->longread && transfrag[t]->nodes[0]==path[i] && transfrag[t]->abundance && (istranscript[t] || ((pathpat & transfrag[t]->pattern)==transfrag[t]->pattern))) {
			//if(transfrag[t]->nodes[0]==path[i] && transfrag[t]->abundance && (istranscript[t] || ((pathpat & transfrag[t]->pattern)==transfrag[t]->pattern))) {
				bool keeptr=true;
				if(i==0) max_fl=transfrag[t]->abundance; // this is the flow from source on this path
				else {
					int ti=1;
					int pi=i+1;
					int lenp=0;
					while(ti<transfrag[t]->nodes.Count()) {
						if(path[pi]!=transfrag[t]->nodes[ti]) { // I found a gap in transfrag => I need to check if it's not too big and can not support path
							// node on path is coming before transfrag; otherwise I wouldn't have the match above
							lenp+=no2gnode[path[pi]]->len();
							if(lenp>CHI_WIN) {
								keeptr=false;
								break;
							}
						}
						else ti++; // I can only advance path as I reach it
						pi++;
					}
				}
				if(keeptr) {
					istranscript[t]=1;
					int n1=i;
					int n2=node2path[transfrag[t]->nodes.Last()];
					//fprintf(stderr,"t=%d n1=%d n2=%d(%d) ",t,n1,n2,transfrag[t]->nodes.Last());
					if(!capacity[n1][n2]) { // haven't seen this link before
						link[n1].Add(n2);
						link[n2].Add(n1);
					}
					capacity[n1][n2]+=transfrag[t]->abundance;
					//fprintf(stderr,"capacity[%d][%d]=%f after adding transfrag[%d]->abundance=%f\n",n1,n2,capacity[n1][n2],t,transfrag[t]->abundance);
				}
			}
			if(i>1 && i<n-2) if(transfrag[t]->longread){ // rate needs to stay 1 for nodes next to source and sink
				if(path[i]==transfrag[t]->nodes[0]) sumright+=transfrag[t]->abundance; // how many transfrags exit node
				if(path[i]==transfrag[t]->nodes.Last()) sumleft+=transfrag[t]->abundance; // how many transfrags enter node
			}
		}
		if(sumright && sumleft) {
			noderate[i]=sumright/sumleft;
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
		float increment=max_fl;
		rate[r++]=1;
		for(int u=n-1;pred[u]>=0;u=pred[u]) {
			float adjflux=(capacity[pred[u]][u]-flow[pred[u]][u])*rate[r-1];
			increment = increment < adjflux ? increment : adjflux; // minimum flux increment on the path
			if(pred[pred[u]]>=0) {
				if(pred[u]<u) {
					if(pred[pred[u]]<pred[u]) rate[r]=rate[r-1]*noderate[pred[u]]; //no2gnode[path[pred[u]]]->rate;
					else rate[r]=rate[r-1];
				}
				else {
					if(pred[pred[u]]<pred[u]) rate[r]=rate[r-1];
					else rate[r]=rate[r-1]/noderate[pred[u]];  //no2gnode[path[pred[u]]]->rate;
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
					if(flow[n1][n2]>0) {
						if(flow[n1][n2]<transfrag[t]->abundance) {
							if(!i) sumout+=flow[n1][n2];
							//fprintf(stderr,"Update capacity of transfrag[%d] with value %f to %f\n",t,transfrag[t]->abundance, transfrag[t]->abundance-flow[n1][n2]);
							update_capacity(0,transfrag[t],flow[n1][n2],nodecapacity,node2path);
							//if(path[i] && transfrag[t]->nodes.Last()!=gno-1) fragno+=flow[n1][n2];
							if(transfrag[t]->abundance<DBL_ERROR) transfrag[t]->abundance=0; // stricter threshold for long transfrags
							flow[n1][n2]=0;
						}
						else {
							if(!i) sumout+=transfrag[t]->abundance;
							flow[n1][n2]-=transfrag[t]->abundance;
							//if(path[i] && transfrag[t]->nodes.Last()!=gno-1) fragno+=transfrag[t]->abundance;
							//fprintf(stderr,"Update capacity of transfrag[%d] with value=%f to 0\n",t,transfrag[t]->abundance);
							update_capacity(0,transfrag[t],transfrag[t]->abundance,nodecapacity,node2path);
							if(transfrag[t]->abundance<DBL_ERROR) transfrag[t]->abundance=0; // stricter threshold for long transfrags
						}
					}
				}
				else if(!i && transfrag[t]->nodes.Last()==path[i]) pos=j;
			}
		}
		if(!i && pos>-1) { // this is first node -> adjust entering transfrag
			int t=no2gnode[path[i]]->trf[pos];
			float val=sumout/noderate[i]; //no2gnode[path[i]]->rate;
			transfrag[t]->abundance-=val;
			if(transfrag[t]->abundance<DBL_ERROR) transfrag[t]->abundance=0;
		}
	}

	// clean up
	delete [] capacity;
	delete [] flow;
	delete [] link;

	return(flux);
}


// version of push_max_flow where I weight the incomplete transfrags
float push_max_flow(int gno,GVec<int>& path,GBitVec& istranscript,GPVec<CTransfrag>& transfrag,GPVec<CGraphnode>& no2gnode,
		GVec<float>& nodeflux,GBitVec& pathpat, GIntHash<int> &gpos, bool &full) {

	int n=path.Count();
	GVec<int> node2path;
	node2path.Resize(gno,-1);
	for(int i=0;i<n;i++) {
		node2path[path[i]]=i;
		nodeflux.cAdd(0.0);
	}
	GVec<float> capacityleft;   // how many transcripts compatible to path enter node
	GVec<float> capacityright;  // how many transcripts compatible to path exit node
	capacityleft.Resize(n);
	capacityright.Resize(n);
	GVec<float> sumleft;        // how many transcripts enter node
	GVec<float> sumright;       // how many transcripts exit node
	sumleft.Resize(n);
	sumright.Resize(n);

	//bool full=true;
	//if(longreads && path.Count()>3) full=false;

	/*
	{ // DEBUG ONLY
		//printTime(stderr);
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
				bool keeptr=false;
				if(istranscript[t]) keeptr=true;
				else if(!transfrag[t]->nodes[0]) {
					if(transfrag[t]->nodes.Last()==path[1]) keeptr=true;
				}
				else if(transfrag[t]->nodes.Last()==path.Last()) {
					if(transfrag[t]->nodes[0]==path[n-2]) keeptr=true;
				}
				else if(transfrag[t]->nodes[0]==path[i] && ((pathpat & transfrag[t]->pattern)==transfrag[t]->pattern)) { // only need to check transfrag the first time I encounter it
					keeptr=true;

					if(longreads) { // an extremely gapped transcript should not be considered to support path (I am doing this for longreads but it might work for paired reads too
						int ti=1;
						int pi=i+1;
						int lenp=0;
						while(ti<transfrag[t]->nodes.Count()) {
							if(path[pi]!=transfrag[t]->nodes[ti]) { // I found a gap in transfrag => I need to check if it's not too big and can not support path
								// node on path is coming before transfrag node; otherwise I wouldn't have the match above
								lenp+=no2gnode[path[pi]]->len();
								if(lenp>CHI_WIN) {
									keeptr=false;
									break;
								}
							}
							else ti++;
							pi++;
						}
					}

					if(keeptr && !full) { // check if transcript fully supports path (full is false means I have not found any transcript to fully support path)
						full=true;
						int p=1;
						if(!transfrag[t]->longstart || !transfrag[t]->longend) full=false;
						if(full) while(path[p]<transfrag[t]->nodes[0]) {
							if(no2gnode[path[p]]->end+1<no2gnode[path[p+1]]->start) {
								full=false;
								break;
							}
							p++;
						}
						if(full) {
							p=path.Count()-2;
							while(path[p]>transfrag[t]->nodes.Last()) {
								if(no2gnode[path[p-1]]->end+1<no2gnode[path[p]]->start) {
									full=false;
									break;
								}
								p--;
							}
						}
						if(full) {
							for(p=2;p<path.Count()-2;p++) {
								if(!transfrag[t]->pattern[path[p]] && no2gnode[path[p]]->len()>longintronanchor) {
									full=false;
									break;
								}
							}
						}
					}
				}

				if(keeptr) { // transcript on path
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
						else if(transfrag[t]->usepath>-1 && int(transfrag[t]->usepath)<transfrag[t]->path.Count()) { //TODO: this crashes with intv.gtf guides -> to fix
							capacityleft[i]+=transfrag[t]->abundance*transfrag[t]->path[int(transfrag[t]->usepath)].abundance;
						}

						//fprintf(stderr,"add transfrag t=%d i=%d sumleft=%f capacityleft=%f\n",t,i,sumleft[i],capacityleft[i]);
					}
					if(transfrag[t]->nodes.Last()>path[i]) { // transfrag ends after this node
						sumright[i]+=transfrag[t]->abundance;
						if(transfrag[t]->real) capacityright[i]+=transfrag[t]->abundance;
						else if(transfrag[t]->usepath>-1 && int(transfrag[t]->usepath)<transfrag[t]->path.Count())
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

	//if(!full) return(0);

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
					if(transfrag[t]->usepath>-1 && int(transfrag[t]->usepath)<transfrag[t]->path.Count())
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
	}

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
	capacityleft.Resize(n);
	capacityright.Resize(n);
	GVec<float> sumleft;        // how many transcripts enter node
	GVec<float> sumright;       // how many transcripts exit node
	sumleft.Resize(n);
	sumright.Resize(n);

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

	bool marginal=false;
	if(longreads) marginal=true;

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

					if(marginal && transfrag[t]->nodes[0] && transfrag[t]->nodes.Last()!=gno-1) marginal=false;

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

	if(marginal) return(0);

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

	for(int i=0;i<n;i++) {
		node2path[path[i]]=i;
		capacity[i].Resize(n);
		flow[i].Resize(n);
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

	//fprintf(stderr,"Process guide %d with pattern: ",g);
	//printBitVec(guidetrf[g].trf->pattern);

	for(int i=0;i<n;i++) {
		node2path[guidetrf[g].trf->nodes[i]]=i;
		nodeflux.cAdd(0.0);
	}

	GVec<float> capacityleft;	// how many transcripts compatible to path enter node
	GVec<float> capacityright;  // how many transcripts compatible to path exit node
	capacityleft.Resize(n);
	capacityright.Resize(n);
	GVec<float> sumleft;        // how many transcripts enter node
	GVec<float> sumright;       // how many transcripts exit node
	sumleft.Resize(n);
	sumright.Resize(n);

	// compute capacities and sums for all nodes
	for(int i=1;i<n-1;i++) {
		int pathi=guidetrf[g].trf->nodes[i];

		int nt=no2gnode[pathi]->trf.Count();
		for(int j=0;j<nt;j++) {
			int t=no2gnode[pathi]->trf[j];
			if(transfrag[t]->abundance) {
				if(istranscript[t] || ((guidetrf[g].trf->pattern & transfrag[t]->pattern)==transfrag[t]->pattern)) { // transcript on path
					istranscript[t]=1;
					// check if there are other guides sharing this transcript so that I can allocate proportionally to guide abundances
					float totalcov=guidetrf[g].trf->abundance;
					for(int r=g-1;r>=0;r--) if((guidetrf[r].trf->pattern & transfrag[t]->pattern)==transfrag[t]->pattern) {
						totalcov+=guidetrf[r].trf->abundance;
					}
					float prop=1;
					if(totalcov>guidetrf[g].trf->abundance) prop=guidetrf[g].trf->abundance/totalcov;

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
	through.Resize(n);

	for(int i=0;i<m;i++) {
		if(i<n) node2path[path[i]]=i;
		if(i<n) nodecapacity.cAdd(0.0);
		capacity[i].Resize(m);
		flow[i].Resize(m);
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
	GIntHash<float> tabund;

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
								//GStr tid(t);
								tabund.Add(t,flow[n1][n2]);
								flow[n1][n2]=0;
							}
							else {
								flow[n1][n2]-=transfrag[t]->abundance;
								//GStr tid(t);
								tabund.Add(t,transfrag[t]->abundance);
							}
						}
					}
					else if(transfrag[t]->nodes[0]<path[i] && transfrag[t]->nodes.Last()>path[i] && transfrag[t]->pattern[path[i]]) { // through transfrag
						//GStr tid(t);
						const float *abund=tabund[t];
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
			for(int i=0;i<m;i++) {
				flow[i].Clear();
				flow[i].Resize(m);
			}

		iterations++;

	}

	// adjust transfrag abundances
	for(int i=0;i<n;i++) {
		int nt=no2gnode[path[i]]->trf.Count();
		for(int j=0;j<nt;j++) {
			int t=no2gnode[path[i]]->trf[j];
			if(istranscript[t] && transfrag[t]->abundance && transfrag[t]->nodes[0]==path[i]) {
				//GStr tid(t);
				const float *abund=tabund[t];
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
		capacity[i].Resize(n);
		flow[i].Resize(n);
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


float store_transcript(GList<CPrediction>& pred,GVec<int>& path,GVec<float>& nodeflux,GVec<float>& nodecov,
		GPVec<CGraphnode>& no2gnode,int& geneno,bool& first,int strand,int gno,GIntHash<int>& gpos, bool& included,
		GBitVec& prevpath, bool full=false,BundleData *bdata=NULL, //float fragno, char* id=NULL) {
		   GffObj* t=NULL) {


	float cov=0;
	int len=0;
	CGraphnode *prevnode=NULL;
	GVec<GSeg> exons;
	GVec<float> exoncov;
	float excov=0;

	uint refstart=0;
	if(bdata) refstart=(uint)bdata->start;

	/*
	fprintf(stderr,"store transcript path[0]=%d",path[0]);
	if(t) fprintf(stderr," with id=%s",t->getID());
	fprintf(stderr,"\n");
	*/

	int s=0;
	if(!path[0]) s=1;

	bool firstex=true;
	bool lastex=false;

	for(int i=s;i<path.Count()-1;i++) {
		if(!prevpath[path[i]]) { // if I can find a node that was not included previously in any path then this is a new path
			included=false;
			prevpath[path[i]]=1;
		}
		int *pos=gpos[edge(path[i-1],path[i],gno)]; // if I can find an edge that was not included in any previous path then this is a new path
		if(i && pos && !prevpath[*pos]) {
			included=false;
			if(pos) prevpath[*pos]=1;
		}

		CGraphnode *node=no2gnode[path[i]];

	    // moved this one before the usedcov computation since that one wasn't used
		if(t && (node->end<t->start || lastex)) { // I am skipping the nodes that do not overlap the guide so that I don't add them up to coverage
			prevnode=node; continue;
		}

		// push
		float usedcov=nodecov[path[i]]*nodeflux[i]*(node->end-node->start+1);
		//fprintf(stderr,"usedcov=%f for nodecov[path[%d]]=%f nodeflux[%d]=%f node->end=%d node->start=%d\n",usedcov,i,nodecov[path[i]],i,nodeflux[i],node->end,node->start);

		uint nodestart=node->start;
		uint nodeend=node->end;

		if(t) { // I am adjusting the start/end of the exon but shouldn't I also adjust the coverage? -> I added two ifs below
			float firstprop=0;
			float lastprop=0;
			if(firstex) {
				if(nodestart<t->start && bdata) {

					float rightcov=0;

					// cummulative bpcov
					float leftcov=get_cov(2*strand,node->start-refstart,t->start-1-refstart,bdata->bpcov);
					if(node->end>t->start) rightcov=get_cov(2*strand,t->start-refstart,node->end-refstart,bdata->bpcov);

					if(leftcov) firstprop=leftcov/(leftcov+rightcov);
				}
				nodestart=t->start;
			}
			if(t->end<=nodeend || i==path.Count()-2) {
				lastex=true;
				if(t->end<node->end && bdata) {
					//usedcov*=(t->end-nodestart+1)/(nodeend-nodestart+1); // this way I am keeping coverage proportions right
					float leftcov=0;

					// cummulative bpcov
					if(node->start<t->end) {
						leftcov=get_cov(2*strand,node->start-refstart,t->end-refstart,bdata->bpcov);
					}
					float rightcov=get_cov(2*strand,t->end+1-refstart,node->end-refstart,bdata->bpcov);

					if(rightcov) lastprop=rightcov/(leftcov+rightcov);
				}
				nodeend=t->end;
			}
			if(firstprop || lastprop) {
				usedcov-=(firstprop +lastprop)*usedcov;
				nodeflux[i]*=(1-firstprop-lastprop);
			}
		}

		nodecov[path[i]]*=(1-nodeflux[i]); // don't allow this to be less than 0

		if(!prevnode || firstex || node->start>prevnode->end+1) { // this is a new exon
			if(prevnode && !firstex) { // compute exon coverage
				excov/=exons.Last().end-exons.Last().start+1;
				exoncov.Add(excov);
				excov=0;
			}
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
		fprintf(stderr,"Predicted transcript cov=%f usedcov=%f len=%d path.count=%d ",cov/len, cov,len,path.Count());
		fprintf(stderr,"and exons cov:");
		for(int e=0;e<exons.Count();e++) fprintf(stderr," %d-%d",exons[e].start,exons[e].end);
		fprintf(stderr,"\n");
		if(t) fprintf(stderr,"Ref_id=%s\n",t->getID());
	}
	*/

	// Add last exon coverage
	if(prevnode) { // compute exon coverage
		excov/=exons.Last().end-exons.Last().start+1;
		exoncov.Add(excov);
	}
	if(len) cov/=len;

	//if(t || (cov>=readthr && len>=mintranscriptlen)) { // store transcript here; also accept some coverage fuzziness that would get eliminated later
	//if(t || (cov>=1 && len>=mintranscriptlen)) { // store transcript here; also accept some coverage fuzziness that would get eliminated later
	// sensitive mode:
	if(t || (cov && len>=mintranscriptlen)) { // store transcript here;
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
			fprintf(stderr,"Store transcript in prediction %d with coverage %f \n",pred.Count(),gcov);
			fprintf(stderr,"And exons cov:");
			for(int e=0;e<exons.Count();e++) fprintf(stderr," %g",exoncov[e]);
			fprintf(stderr,"\n");
		}
		*/

		CPrediction *p=new CPrediction(geneno-1, t, exons[0].start, exons.Last().end, gcov, sign, len);
		p->exons=exons;
		if(t && t->exons.Count()==1) exoncov[0]=gcov;
		p->exoncov=exoncov;
		if(full) p->mergename+='.';
		if(longreads) p->tlen=-p->tlen;
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
				pred[np]->cov+=excov/abs(pred[np]->tlen);
				e++;
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
	if(longreads) p->tlen=-p->tlen;
	pred.Add(p);

	if (t && t->uptr) {
		RC_TData &td = *(RC_TData*) (t->uptr);
		td.in_bundle=3;
		//fprintf(stderr,"sg guide %s is stored\n",t->getID());
	}

	update_guide_pred(pred,np,path,nodeflux,nodecov,no2gnode,gno,update);

	return(np);
}

/*
int best_trf_match(CTransfrag *t,GVec<CTransfrag>& keeptrf,GPVec<CGraphnode>& no2gnode,int gno) {
	int mineditdist=no2gnode[gno-2]->end-no2gnode[1]->start+1; // bp different between the 2
	int mininternaldist=mineditdist; // internal bp that are different -> preferably equal to 0
	int maxintersect=0; // nodes in common
	int mink=-1;
	for(int k=0;k<keeptrf.Count();k++) if(t->nodes[0]<=keeptrf[k].nodes[keeptrf[k].nodes.Count()-2] && keeptrf[k].nodes[1]<=t->nodes.Last()){
		int editdist=0;
		int internaldist=0;
		int intersect=0;
		int i=0; // t index
		int j=1; // k index

		while(j<keeptrf[k].nodes.Count()-1 && t->nodes[i]>keeptrf[k].nodes[j]) {
			editdist+=no2gnode[keeptrf[k].nodes[j]]->len();
			j++;
		}

		while(i<t->nodes.Count() && t->nodes[i]<keeptrf[k].nodes.Last() && j<keeptrf[k].nodes.Count()-1) {
			if(t->nodes[i]<keeptrf[k].nodes[j]) {
				editdist+=no2gnode[t->nodes[i]]->len();
				if(j>1) internaldist+=no2gnode[t->nodes[i]]->len();
				i++;
			}
			else if(t->nodes[i]>keeptrf[k].nodes[j]) {
				editdist+=no2gnode[keeptrf[k].nodes[j]]->len();
				internaldist+=no2gnode[keeptrf[k].nodes[j]]->len();
				j++;
			}
			else { // nodes are equal
				i++;
				j++;
				intersect++;
			}
		}

		if(!intersect || internaldist>CHI_THR) continue;

		if(intersect>maxintersect) { // more nodes in common
			mink=k;
			while(i<t->nodes.Count()) { // i>=1 because t intersects k
				editdist+=no2gnode[t->nodes[i]]->len();
				if(no2gnode[t->nodes[i-1]]->end+1<no2gnode[t->nodes[i]]->start) internaldist+=no2gnode[t->nodes[i]]->len();
				i++;
			}
			while(j<keeptrf[k].nodes.Count()-1) {
				editdist+=no2gnode[keeptrf[k].nodes[j]]->len();
				j++;
			}
			mineditdist=editdist;
			mininternaldist=internaldist;
			maxintersect=intersect;
		}
		else if(intersect==maxintersect) {
			if(internaldist<mininternaldist) {
				while(i<t->nodes.Count()) { // i>=1 because t intersects k
					editdist+=no2gnode[t->nodes[i]]->len();
					if(no2gnode[t->nodes[i-1]]->end+1<no2gnode[t->nodes[i]]->start) internaldist+=no2gnode[t->nodes[i]]->len();
					i++;
				}
				while(j<keeptrf[k].nodes.Count()-1) {
					editdist+=no2gnode[keeptrf[k].nodes[j]]->len();
					j++;
				}
				if(internaldist<mininternaldist) {
					mink=k;
					mineditdist=editdist;
					mininternaldist=internaldist;
				}
			}
			else if(internaldist==mininternaldist) {
				if(mink>=0) {
					if(editdist>=mineditdist) { // the only time I can replace mink here is if a(mink)==0 && a(k)>0
						if(keeptrf[mink].abundance || !keeptrf[k].abundance) continue;
					}
					else if(keeptrf[mink].abundance && !keeptrf[k].abundance) continue; // there is no way I can replace mink in this case
				}

				// finish computing editdist
				while(i<t->nodes.Count()) {
					editdist+=no2gnode[t->nodes[i]]->len();
					i++;
				}
				while(j<keeptrf[k].nodes.Count()-1) {
					editdist+=no2gnode[keeptrf[k].nodes[j]]->len();
					j++;
				}

				if(mink<0 || (!keeptrf[mink].abundance && keeptrf[k].abundance)) { // I always replace in this case
					mink=k;
					mineditdist=editdist;
				}
				else if(editdist<mineditdist && (!keeptrf[mink].abundance || keeptrf[k].abundance)) {
						mink=k;
						mineditdist=editdist;
				}

			}
		}
	}
	return(mink);
}
*/

float best_trf_match(CTransfrag *t,GVec<CTransfrag>& keeptrf,GPVec<CGraphnode>& no2gnode,int gno,GVec<int>& tmatch) {
	int mininternaldist=no2gnode[gno-2]->end-no2gnode[1]->start+1; // bp different between the 2 internally (in between nodes) -> preferably equal to 0
	int maxintersect=0; // nodes in common
	float abundancesum=0;
	for(int k=0;k<keeptrf.Count();k++) if(t->nodes[0]<=keeptrf[k].nodes[keeptrf[k].nodes.Count()-2] && keeptrf[k].nodes[1]<=t->nodes.Last()){
		int internaldist=0;
		int intersect=0;
		int i=0; // t index
		int j=1; // k index

		while(j<keeptrf[k].nodes.Count()-1 && t->nodes[i]>keeptrf[k].nodes[j]) {
			j++;
		}

		while(i<t->nodes.Count() && t->nodes[i]<keeptrf[k].nodes.Last() && j<keeptrf[k].nodes.Count()-1) {
			if(t->nodes[i]<keeptrf[k].nodes[j]) {
				if(j>1) internaldist+=no2gnode[t->nodes[i]]->len();
				i++;
			}
			else if(t->nodes[i]>keeptrf[k].nodes[j]) {
				internaldist+=no2gnode[keeptrf[k].nodes[j]]->len();
				j++;
			}
			else { // nodes are equal
				i++;
				j++;
				intersect++;
			}
		}

		if(!intersect || internaldist>CHI_THR) continue;

		if(intersect>maxintersect) { // more nodes in common
			tmatch.Clear();
			tmatch.Add(k);
			abundancesum=keeptrf[k].abundance;
			while(i<t->nodes.Count()) { // i>=1 because t intersects k
				if(no2gnode[t->nodes[i-1]]->end+1<no2gnode[t->nodes[i]]->start) internaldist+=no2gnode[t->nodes[i]]->len();
				i++;
			}
			while(j<keeptrf[k].nodes.Count()-1) {
				j++;
			}
			mininternaldist=internaldist;
			maxintersect=intersect;
		}
		else if(intersect==maxintersect) { // same number of nodes that intersect
			if(internaldist<mininternaldist) {
				while(i<t->nodes.Count()) { // i>=1 because t intersects k
					if(no2gnode[t->nodes[i-1]]->end+1<no2gnode[t->nodes[i]]->start) internaldist+=no2gnode[t->nodes[i]]->len();
					i++;
				}
				while(j<keeptrf[k].nodes.Count()-1) {
					j++;
				}
				if(internaldist<mininternaldist) { // better internal distance
					tmatch.Clear();
					tmatch.Add(k);
					mininternaldist=internaldist;
					abundancesum=keeptrf[k].abundance;
				}
			}
			else if(internaldist==mininternaldist) { // same internal distance --> keep them both
				tmatch.Add(k);
				abundancesum+=keeptrf[k].abundance;
			}
		}
	}
	return(abundancesum);
}

/*
void get_trf_long(int gno,int edgeno, GIntHash<int> &gpos,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,
		int& geneno,int strand,GList<CPrediction>& pred,GVec<int>& trflong,BundleData *bdata) {

	GPVec<GffObj>& guides = bdata->keepguides;

	GVec<float> nodecov; // the coverage of all transfrags entering a node
	GVec<float> noderate;
	for(int i=0;i<gno;i++) {
		CGraphnode *inode=no2gnode[i]; // this is here only because of the DEBUG option below
		nodecov.cAdd(0.0); // this are all transfrags that link nodes together
		float rate=1;
		if(i && i<gno-1) {
			for(int j=0;j<inode->trf.Count();j++) { // for all transfrags going through node
				int t=inode->trf[j];
				if(transfrag[t]->nodes[0]<i) { // entering transfrags:
					nodecov[i]+=transfrag[t]->abundance;
				}
			}

			if(nodecov[i]) rate=nodecov[i];
			if(rate<=0) rate=1; // this shouldn't happen
			rate=inode->cov/rate;
		}
		noderate.Add(rate);
	}

	GBitVec istranscript(transfrag.Count());

	 GVec<int> path;
	 GBitVec pathpat(gno+edgeno);
	 int minp;
	 int maxp;
	 int maxi;

	 char sign='-';
	 if(strand) { sign='+';}
	 int npred=pred.Count();

	 GVec<CTransfrag> keeptrf;
	 GVec<int> checktrf;
	 for(int f=trflong.Count()-1;f>=0;f--) { // if this is a guide it should be reflected in the prediction downstream
		 path.Clear();
		 int t=trflong[f];
		 if(t<0) GError("Stored long transcript is negative!\n");
		 pathpat=transfrag[t]->pattern;
		 minp=transfrag[t]->nodes[0];
		 maxp=transfrag[t]->nodes.Last();

		 //if(no2gnode[transfrag[t]->nodes[0]]->hardstart) {
			 int *pos=gpos[edge(0,minp,gno)];
			 if(pos) pathpat[*pos]=1;
			 //minp=0;
		 //}
		 //if(no2gnode[transfrag[t]->nodes.Last()]->hardend) {
			 pos=gpos[edge(maxp,gno-1,gno)];
			 if(pos) pathpat[*pos]=1;
			 //maxp=gno-1;
		 //}

		 maxi=minp;
		 path.Add(maxi);
		 pathpat[maxi]=1;

		 istranscript.reset();

		 float flux=0;
		 //float fragno=0;
		 GVec<float> nodeflux;

		 bool tocheck=true;
		 if(back_to_source_fast_long(maxi,path,minp,maxp,pathpat,transfrag,no2gnode,nodecov,gno,gpos)) {
			 path.cAdd(0);
			 path.Reverse(); // back to source adds the nodes at the end to avoid pushing the list all the time

			 if(fwd_to_sink_fast_long(maxi,path,minp,maxp,pathpat,transfrag,no2gnode,nodecov,gno,gpos)) {

				 flux=long_max_flow(gno,path,istranscript,transfrag,no2gnode,nodeflux,pathpat);


				 { // DEBUG ONLY
					 //printTime(stderr);
					 fprintf(stderr,"flux=%g Path:",flux);
					 for(int i=0;i<path.Count();i++) fprintf(stderr," %d",path[i]);
					 fprintf(stderr,"\n");
					 fprintf(stderr,"Nodecapacities:");
					 for(int i=0;i<path.Count();i++) fprintf(stderr," %f",nodeflux[i]);
					 fprintf(stderr,"***\n");
				 }


				 if(flux) { // these are not valid paths in the graph
					 tocheck=false;

					 GVec<GSeg> exons;
					 GVec<float> exoncov;
					 int j=1;
					 int len=0;
					 float cov=0;
					 int startnode=j;
					 int lastnode=path.Count()-2;
					 uint startpoint=no2gnode[path[lastnode]]->end;
					 uint endpoint=no2gnode[path[1]]->start;

					 if(mixedMode) { // establish start/end point of path
						 while(j<path.Count()-1) {
							 CGraphnode *jnode=no2gnode[path[j]];
							 for(int i=0;i<jnode->trf.Count();i++) { // for all transfrags going through node
								 int t=jnode->trf[i];
								 if(istranscript[t] && transfrag[t]->longread && transfrag[t]->nodes[0] && transfrag[t]->nodes.Last()!=gno-1) {
									 if(transfrag[t]->nodes[0]==path[j]) {
										 if(transfrag[t]->longstart) {
											 if(transfrag[t]->longstart<startpoint) { startpoint=transfrag[t]->longstart; startnode=j;}
										 }
										 else if(jnode->start<startpoint) { startpoint=jnode->start; startnode=j;}
									 }
									 if(transfrag[t]->nodes.Last()==path[j]) {
										 if(transfrag[t]->longend) {
											 if(transfrag[t]->longend>endpoint) { endpoint=transfrag[t]->longend; lastnode=j;}
										 }
										 else if(jnode->end>endpoint) {endpoint=jnode->end; lastnode=j;}
									 }
								 }
							 }
							 j++;
						 }
						 if(startpoint>endpoint) {
							 j=1;
							 lastnode=path.Count()-2;
							 startpoint=no2gnode[path[1]]->start;
							 endpoint=no2gnode[path[lastnode]]->end;
						 }
						 else j=startnode;
					 }

					 while(j<=lastnode) {
					 //while(j<path.Count()-1) {
						 int nodestart=no2gnode[path[j]]->start;
						 if(mixedMode && j==startnode) nodestart=startpoint;
						 int nodeend=no2gnode[path[j]]->end;
						 if(mixedMode && j==lastnode) nodeend=endpoint;
						 nodecov[path[j]]-=nodeflux[j];
						 len+=nodeend-nodestart+1;
						 float ecov=nodeflux[j]*noderate[path[j]];
						 float excov=ecov;
						 //if(mixedMode) {
							// no2gnode[path[j]]->cov-=ecov;
							 //if(no2gnode[path[j]]->cov<epsilon) no2gnode[path[j]]->cov=0;
						 //}
						 //float excov=nodeflux[j]*noderate[path[j]];
						 while(j+1<=lastnode && no2gnode[path[j]]->end+1==no2gnode[path[j+1]]->start) {
						 //while(j+1<path.Count()-1 && no2gnode[path[j]]->end+1==no2gnode[path[j+1]]->start) {
							 j++;
							 nodeend=no2gnode[path[j]]->end;
							 ecov=nodeflux[j]*noderate[path[j]];
							 excov+=ecov;
							 if(mixedMode) {
								 //no2gnode[path[j]]->cov-=ecov;
								 //if(no2gnode[path[j]]->cov<epsilon) no2gnode[path[j]]->cov=0;
								 if(j==lastnode) nodeend=endpoint;
							 }
							 len+=nodeend-no2gnode[path[j]]->start+1;

						 }
						 GSeg exon(nodestart,nodeend);
						 exons.Add(exon);
						 cov+=excov;
						 exoncov.Add(excov);
						 j++;
					 }
					 if(transfrag[t]->nodes.Count()==1) transfrag[t]->abundance=0;
					 GffObj *g=NULL;
					 if(transfrag[t]->guide && longreads) {
						 g=guides[int(transfrag[t]->usepath)];
						 if (g && g->uptr) {
							 RC_TData &td = *(RC_TData*) (g->uptr);
							 td.in_bundle=3;
						 }
					 }
					 if(!eonly || g) {
						 CPrediction *p=new CPrediction(geneno, g,exons[0].start , exons.Last().end, cov, sign, len);
						 p->exons=exons;
						 p->exoncov=exoncov;
						 p->mergename='.'; // I should not delete this prediction
						 p->tlen=-p->tlen; // negative transcript length signifies assembly is from a long read
						 pred.Add(p);

						 CTransfrag u(path,pathpat,cov/len);
						 keeptrf.Add(u);
					 }
				 }
				 else if(transfrag[t]->guide) {
					 checktrf.Add(t);
				 }
			 }
		 }

		 if(tocheck)  { // try to see if you can rescue transfrag -> they are stored from more abundant to least -> not if using mixedMode
			 checktrf.Add(t);
		 }
	 }

	 //keeptrf.Sort(longtrCmp); // most abundant transfrag in the graph come first, then the ones with most nodes, then the ones more complete

	 for(int c=0;c<checktrf.Count();c++) if(transfrag[checktrf[c]]->guide || transfrag[checktrf[c]]->abundance>=readthr) { // only in this case it is worth considering it as a potential prediction
		 int t=checktrf[c];

		 int mink=best_trf_match(transfrag[t],keeptrf,no2gnode,gno);
		 if(mink>=0) { // found good match
			 int p=0;
			 int i=0;
			 int np=npred+mink; // how do I know that keeptrf lead to a prediction -> because keeptrf represent the predictions that were added
			 while(i<transfrag[t]->nodes.Count() && p<pred[np]->exons.Count()) {
				 if(no2gnode[transfrag[t]->nodes[i]]->end<pred[np]->exons[p].start) i++;
				 else if(pred[np]->exons[p].end<no2gnode[transfrag[t]->nodes[i]]->start) p++;
				 else { // the two intersect (I can only have the full node included in exon)
					 if(nodecov[transfrag[t]->nodes[i]]>epsilon){
						 //float addcov=transfrag[t]->abundance*noderate[transfrag[t]->nodes[i]];
						 float addcov=transfrag[t]->abundance*no2gnode[transfrag[t]->nodes[i]]->len();
						 float newnodecov=nodecov[transfrag[t]->nodes[i]]-addcov/noderate[transfrag[t]->nodes[i]];
						 if(newnodecov<0) {
							 addcov=nodecov[transfrag[t]->nodes[i]]*noderate[transfrag[t]->nodes[i]];
							 newnodecov=0;
						 }
						 nodecov[transfrag[t]->nodes[i]]=newnodecov;

						 pred[np]->exoncov[p]+=addcov;
						 pred[np]->cov+=addcov;
					 }
					 i++;
				 }
			 }
		 }
		 else if(!eonly || transfrag[t]->guide) { // store it as an independent prediction
			 pathpat=transfrag[t]->pattern; // not used right now but maybe in the future?
			 path.Clear();
			 path.cAdd(0);
			 path.Add(transfrag[t]->nodes[0]);
			 for(int j=1;j<transfrag[t]->nodes.Count();j++) {
				 if(transfrag[t]->nodes[j]!=1+transfrag[t]->nodes[j-1] ||
						 no2gnode[transfrag[t]->nodes[j]]->start-1!=no2gnode[transfrag[t]->nodes[j-1]]->end) {
					 // check if transfrag t1 is incomplete between node[n-1] and node [n]
					 int *pos=gpos[edge(transfrag[t]->nodes[j-1],transfrag[t]->nodes[j],gno)];
					 if(!pos || !transfrag[t]->pattern[*pos]) { // incomplete transfrag
						 break;
					 }
					 if(pos) pathpat[*pos]=1;
					 path.Add(transfrag[t]->nodes[j]);
				 }
				 else path.Add(transfrag[t]->nodes[j]);
			 }
			 if(path.Last()==transfrag[t]->nodes.Last()) { // this transfrag is complete, might be worth rescuing
				 int sink=gno-1;
				 path.Add(sink);

				 GVec<GSeg> exons;
				 GVec<float> exoncov;
				 int j=1;
				 int len=0;
				 float cov=0;
				 while(j<path.Count()-1) {
					 int nodestart=no2gnode[path[j]]->start;
					 int nodeend=no2gnode[path[j]]->end;
					 // nodecov[path[j]]-=transfrag[t]->abundance; // do not need this here anymore
					 len+=nodeend-nodestart+1;

					 float excov=0;
					 if(nodecov[path[j]]>epsilon) {
						 excov=transfrag[t]->abundance*no2gnode[path[j]]->len();
						 float newnodecov=nodecov[path[j]]-excov/noderate[path[j]];
						 if(newnodecov<0) {
							 excov=nodecov[path[j]]*noderate[path[j]];
							 newnodecov=0;
						 }
						 nodecov[path[j]]=newnodecov;
					 }

					 while(j+1<path.Count()-1 && no2gnode[path[j]]->end+1==no2gnode[path[j+1]]->start) {
						 j++;
						 len+=no2gnode[path[j]]->len();
						 nodeend=no2gnode[path[j]]->end;

						 float addcov=0;
						 if(nodecov[path[j]]>epsilon) {
							 addcov=transfrag[t]->abundance*no2gnode[path[j]]->len();
							 float newnodecov=nodecov[path[j]]-addcov/noderate[path[j]];
							 if(newnodecov<0) {
								 addcov=nodecov[path[j]]*noderate[path[j]];
								 newnodecov=0;
							 }
							 nodecov[path[j]]=newnodecov;
						 }
						 excov+=addcov;
						 //excov+=transfrag[t]->abundance*noderate[path[j]];
					 }
					 GSeg exon(nodestart,nodeend);
					 exons.Add(exon);
					 cov+=excov;
					 exoncov.Add(excov);
					 j++;
				 }
				 GffObj *g=NULL;
				 if(transfrag[t]->guide && longreads) {
					 g=guides[int(transfrag[t]->usepath)];
					 if (g && g->uptr) {
						 RC_TData &td = *(RC_TData*) (g->uptr);
						 td.in_bundle=3;
					 }
				 }
				 CPrediction *p=new CPrediction(geneno, g,exons[0].start , exons.Last().end, cov, sign, len);
				 p->exons=exons;
				 p->exoncov=exoncov;
				 p->tlen=-p->tlen; // negative transcript length signifies assembly is from a long read
				 pred.Add(p);

				 CTransfrag u(path,pathpat,cov/len);
				 keeptrf.Add(u);
			 }
		 }
		 transfrag[t]->abundance=0;
	 }

	 if(pred.Count()>npred) {
		 if(mixedMode) {
			 for(int t=0;t<transfrag.Count();t++) {
			   if(transfrag[t]->longread && (!transfrag[t]->nodes[0] || transfrag[t]->nodes.Last()==gno-1) )
			     transfrag[t]->abundance=0;
			   else
			     transfrag[t]->abundance=transfrag[t]->usepath;
			   transfrag[t]->usepath=-1;
			 }
		 }
		 else for(int t=0;t<transfrag.Count();t++) { // longreads mode --> tries to add all transfrags to predictions
			 if(transfrag[t]->longread && transfrag[t]->abundance>epsilon && transfrag[t]->nodes[0] && transfrag[t]->nodes.Last()!=gno-1 ) {
				 int mink=best_trf_match(transfrag[t],keeptrf,no2gnode,gno);

				 if(mink>=0) { // mink gives the prediction's position too
					 int p=0;
					 int i=0;
					 int np=npred+mink;
					 while(i<transfrag[t]->nodes.Count() && p<pred[np]->exons.Count()) {
						 if(no2gnode[transfrag[t]->nodes[i]]->end<pred[np]->exons[p].start) i++;
						 else if(pred[np]->exons[p].end<no2gnode[transfrag[t]->nodes[i]]->start) p++;
						 else { // the two intersect (I can only have the full node included in exon)
							 if(nodecov[transfrag[t]->nodes[i]]>epsilon){
								 //float addcov=transfrag[t]->abundance*noderate[transfrag[t]->nodes[i]];
								 float addcov=transfrag[t]->abundance*no2gnode[transfrag[t]->nodes[i]]->len();
								 float newnodecov=nodecov[transfrag[t]->nodes[i]]-addcov/noderate[transfrag[t]->nodes[i]];
								 if(newnodecov<0) {
									 addcov=nodecov[transfrag[t]->nodes[i]]*noderate[transfrag[t]->nodes[i]];
									 newnodecov=0;
								 }
								 nodecov[transfrag[t]->nodes[i]]=newnodecov;

								 pred[np]->exoncov[p]+=addcov;
								 pred[np]->cov+=addcov;
							 }
							 i++;
						 }
					 }
					 transfrag[t]->abundance=0; // delete abundance in order not to use it in short reads
				 }
			 }
		 }

		 int p=npred;
		 while(p<pred.Count()) {
			 if(pred[p]->cov) {
				 pred[p]->cov/=abs(pred[p]->tlen);
				 for(int i=0;i<pred[p]->exons.Count();i++)
					 pred[p]->exoncov[i]/=pred[p]->exons[i].len();
				 p++;
			 }
			 else if(!eonly) { // || !pred[p]->t_eq) {
				 pred.Delete(p); // I delete all predictions that have 0 coverage unless it's eonly mode
			 }
			 else p++;
		 }
	 }
}
*/

void get_trf_long_mix(int gno,int edgeno, GIntHash<int> &gpos,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,int& geneno,int strand,
		GList<CPrediction>& pred,GVec<int>& trflong,GVec<float>& nodecovall,GBitVec& istranscript,GBitVec& prevpath,BundleData *bdata,bool &first) {

	GPVec<GffObj>& guides = bdata->keepguides;
	GVec<float> nodecov; // the coverage of all transfrags entering a node
	GVec<float> noderate;
	for(int i=0;i<gno;i++) {
		CGraphnode *inode=no2gnode[i]; // this is here only because of the DEBUG option below
		nodecov.cAdd(0.0); // this are all transfrags that link nodes together
		float rate=1;
		if(i && i<gno-1) {
			for(int j=0;j<inode->trf.Count();j++) { // for all transfrags going through node
				int t=inode->trf[j];
				//if(transfrag[t]->longread && transfrag[t]->nodes[0]<i) { // entering transfrags:
				if(transfrag[t]->longread && transfrag[t]->nodes.Last()>i) { // exiting transfrags: this is more consistent with the nodeflux computation
					nodecov[i]+=transfrag[t]->abundance;
				}
			}

			if(nodecov[i]) rate=nodecov[i];
			if(rate<=0) rate=1; // this shouldn't happen
			//fprintf(stderr,"rate=%f\n",rate);
			rate=inode->cov/rate;
		}
		noderate.Add(rate);
		//fprintf(stderr,"Node[%d]:%d-%d no2gnode->cov=%f nodecov=%f noderate=%f\n",i,inode->start,inode->end,inode->cov,nodecov[i],noderate[i]);
	}

	GVec<int> path;
	GBitVec pathpat(gno+edgeno);
	int minp;
	int maxp;
	int maxi;


	char sign='-';
	if(strand) { sign='+';}
	int npred=pred.Count();

	GVec<CTransfrag> keeptrf;
	GVec<int> checktrf;
	for(int f=trflong.Count()-1;f>=0;f--) { // if this is a guide it should be reflected in the prediction downstream
		path.Clear();
		 int t=trflong[f];
		 if(t<0) GError("Stored long transcript is negative!\n");
		 pathpat=transfrag[t]->pattern;
		 minp=transfrag[t]->nodes[0];
		 maxp=transfrag[t]->nodes.Last();

		 int *pos=gpos[edge(0,minp,gno)];
		 if(pos) pathpat[*pos]=1;

		 pos=gpos[edge(maxp,gno-1,gno)];
		 if(pos) pathpat[*pos]=1;

		 maxi=minp;
		 path.Add(maxi);
		 pathpat[maxi]=1;

		 istranscript.reset();

		 float flux=0;
		 GVec<float> nodeflux;

		 /*
	 	 { // DEBUG ONLY
	 	 fprintf(stderr,"\n\n***Start get_trf_long_mix with maxi=%d minp=%d maxp=%d guide=%d and transcript:",maxi,minp,maxp,transfrag[t]->guide);
	 	 for(int i=0;i<transfrag[t]->nodes.Count();i++) fprintf(stderr," %d",transfrag[t]->nodes[i]);
	 	 fprintf(stderr," pathpat=");
	 	 //printBitVec(pathpat);
		 fprintf(stderr,"\n");

#ifdef GMEMTRACE
	 	 double vm,rsm;
	 	 get_mem_usage(vm, rsm);
	 	 GMessage("\t\tM(s):parse_trf memory usage: rsm=%6.1fMB vm=%6.1fMB\n",rsm/1024,vm/1024);
#endif
	 	 }
	 	 */


		 bool tocheck=true;
		 if(back_to_source_fast_long(maxi,path,minp,maxp,pathpat,transfrag,no2gnode,nodecov,gno,gpos)) {
			 path.cAdd(0);
			 path.Reverse(); // back to source adds the nodes at the end to avoid pushing the list all the time

			 if(fwd_to_sink_fast_long(maxi,path,minp,maxp,pathpat,transfrag,no2gnode,nodecov,gno,gpos)) {

				 flux=long_max_flow(gno,path,istranscript,transfrag,no2gnode,nodeflux,pathpat);

				 /*
				 { // DEBUG ONLY
					 //printTime(stderr);
					 fprintf(stderr,"flux=%g Path:",flux);
					 for(int i=0;i<path.Count();i++) fprintf(stderr," %d",path[i]);
					 fprintf(stderr,"\n");
					 fprintf(stderr,"Nodecapacities:");
					 for(int i=0;i<path.Count();i++) fprintf(stderr," %f",nodeflux[i]);
					 fprintf(stderr,"***\n");
				 }
				 */

				 if(flux) { // these are not valid paths in the graph


					 tocheck=false;

					 GVec<GSeg> exons;
					 GVec<float> exoncov;
					 int j=1;
					 int len=0;
					 float cov=0;
					 int startnode=1;
					 int lastnode=path.Count()-2;

					 uint startpoint=no2gnode[path[startnode]]->end;
					 uint endpoint=no2gnode[path[lastnode]]->start;
					 CGraphnode *jnode=no2gnode[path[startnode]];
					 for(int i=0;i<jnode->trf.Count();i++) { // for all transfrags going through startnode
						 int t=jnode->trf[i];
						 if(istranscript[t] && transfrag[t]->longread && transfrag[t]->nodes[0] && transfrag[t]->nodes.Last()!=gno-1) {
							 if(transfrag[t]->nodes[0]==path[startnode]) {
								 if(transfrag[t]->longstart) {
									 if(transfrag[t]->longstart<startpoint) { startpoint=transfrag[t]->longstart;}
								 }
							 }
						 }
					 }
					 jnode=no2gnode[path[lastnode]];
					 for(int i=0;i<jnode->trf.Count();i++) { // for all transfrags going through lastnode
						 int t=jnode->trf[i];
						 if(istranscript[t] && transfrag[t]->longread && transfrag[t]->nodes[0] && transfrag[t]->nodes.Last()!=gno-1) {
							 if(transfrag[t]->nodes.Last()==path[lastnode]) {
								 if(transfrag[t]->longend) {
									 if(transfrag[t]->longend>endpoint) { endpoint=transfrag[t]->longend;}
								 }
							 }
						 }
					 }
					 if(startpoint==no2gnode[path[startnode]]->end) startpoint=no2gnode[path[1]]->start;
					 if(endpoint==no2gnode[path[lastnode]]->start) endpoint=no2gnode[path[lastnode]]->end;


					 while(j<=lastnode) {
						 int nodestart=no2gnode[path[j]]->start;
						 int nodeend=no2gnode[path[j]]->end;
						 len+=nodeend-nodestart+1;
						 nodecov[path[j]]-=nodeflux[j];
						 float ecov=nodeflux[j]*noderate[path[j]];
						 float excov=ecov;
						 /*if(mixedMode) {
							 no2gnode[path[j]]->cov-=ecov;
							 if(no2gnode[path[j]]->cov<epsilon) no2gnode[path[j]]->cov=0;
						 }*/
						 //float excov=nodeflux[j]*noderate[path[j]];
						 while(j+1<=lastnode && no2gnode[path[j]]->end+1==no2gnode[path[j+1]]->start) {
						 //while(j+1<path.Count()-1 && no2gnode[path[j]]->end+1==no2gnode[path[j+1]]->start) {
							 j++;
							 nodeend=no2gnode[path[j]]->end;
							 ecov=nodeflux[j]*noderate[path[j]];
							 nodecov[path[j]]-=nodeflux[j];
							 excov+=ecov;
							 len+=nodeend-no2gnode[path[j]]->start+1;
						 }
						 GSeg exon(nodestart,nodeend);
						 exons.Add(exon);
						 cov+=excov;
						 exoncov.Add(excov);
						 j++;
					 }
					 if(transfrag[t]->nodes.Count()==1) transfrag[t]->abundance=0;
					 //fprintf(stderr,"Store prediction %d  with abundance=%f len=%d\n",pred.Count(),cov/len,len);
					 //GffObj *g=NULL;

					 if(len>=mintranscriptlen) {
						 if(first) { geneno++; first=false;}
						 /*fprintf(stderr,"1 Store prediction %d  with abundance=%f totalabundance=%f len=%d startpoint=%d endpoint=%d and exons:",pred.Count(),cov/len,cov,len,startpoint,endpoint);
						 for(int i=0;i<exons.Count();i++) fprintf(stderr," %d-%d",exons[i].start,exons[i].end);
						 fprintf(stderr,"\n");*/
						 GffObj *g=NULL;
						 if(transfrag[t]->guide) {
							 g=guides[int(transfrag[t]->guide-1)];
							 if (g && g->uptr) {
								 RC_TData &td = *(RC_TData*) (g->uptr);
								 td.in_bundle=3;
								 //fprintf(stderr,"sg guide %s is stored\n",g->getID());
							 }
						 }

						 CPrediction *p=new CPrediction(geneno, g,startpoint , endpoint, cov, sign, len);
						 p->exons=exons;
						 p->exoncov=exoncov;
						 p->mergename='.'; // I should not delete this prediction
						 p->tlen=-p->tlen; // negative transcript length signifies assembly is from a long read
						 pred.Add(p);

						 //fprintf(stderr,"Added prediction=%d with totalcov=%.1f\n",pred.Count()-1,pred.Last()->cov);

						 CTransfrag u(path,pathpat,cov/len);
						 keeptrf.Add(u);
					 }
				 }
			 }
		 }

		 if(tocheck)  { // try to see if you can rescue transfrag
			if(!guided || transfrag[t]->guide || (no2gnode[transfrag[t]->nodes[0]]->parent[0]==0 &&
				   no2gnode[transfrag[t]->nodes.Last()]->child.Last()==gno-1) )
				// only accept long transfrags that are linked to source and sink
			 checktrf.Add(t);
		 }
	 }

	//keeptrf.Sort(longtrCmp); // most abundant transfrag in the graph come first, then the ones with most nodes, then the ones more complete

	 for(int c=0;c<checktrf.Count();c++) if(transfrag[checktrf[c]]->guide || transfrag[checktrf[c]]->abundance>=readthr) { // only in this case it is worth considering it as a potential prediction
		 int t=checktrf[c];

		 /*fprintf(stderr,"checktrf[%d]=%d with abundance=%f with start=%d end=%d nodes:",c,t,transfrag[t]->abundance,transfrag[t]->longstart,transfrag[t]->longend);
		 for(int i=0;i<transfrag[t]->nodes.Count();i++) fprintf(stderr," %d",transfrag[t]->nodes[i]);
		 fprintf(stderr,"\n");*/

		 GVec<int> tmatch;
		 float abundancesum=best_trf_match(transfrag[t],keeptrf,no2gnode,gno,tmatch);
		 if(abundancesum>0) {
			 if(!transfrag[t]->shortread && transfrag[t]->nodes.Count()>1) {
				 for(int j=0;j<tmatch.Count();j++) { // found good match(es) but transcript has more than one node for this to make sense
					 float abundprop=transfrag[t]->abundance*keeptrf[tmatch[j]].abundance/abundancesum;
					 int p=0;
					 int i=0;
					 int np=npred+tmatch[j]; // how do I know that keeptrf lead to a prediction -> because keeptrf represent the predictions that were added
					 //fprintf(stderr,"Add %.1f to prediction %d with cov=%.1f\n",abundprop,np,pred[np]->cov);
					 while(i<transfrag[t]->nodes.Count() && p<pred[np]->exons.Count()) {
						 if(no2gnode[transfrag[t]->nodes[i]]->end<pred[np]->exons[p].start) i++;
						 else if(pred[np]->exons[p].end<no2gnode[transfrag[t]->nodes[i]]->start) p++;
						 else { // the two intersect (I can only have the full node included in exon)
							 if(nodecov[transfrag[t]->nodes[i]]>epsilon){
								 //float addcov=transfrag[t]->abundance*noderate[transfrag[t]->nodes[i]];
								 float addcov=abundprop*no2gnode[transfrag[t]->nodes[i]]->len();
								 float newnodecov=nodecov[transfrag[t]->nodes[i]]-addcov/noderate[transfrag[t]->nodes[i]];
								 //fprintf(stderr,"newnodecov=%.1f\n",newnodecov);
								 if(newnodecov<0) {
									 newnodecov=nodecov[transfrag[t]->nodes[i]]*keeptrf[tmatch[j]].abundance/abundancesum;
									 addcov=newnodecov*noderate[transfrag[t]->nodes[i]];
									 newnodecov=nodecov[transfrag[t]->nodes[i]]-newnodecov;
								 }
								 nodecov[transfrag[t]->nodes[i]]=newnodecov;

								 pred[np]->exoncov[p]+=addcov;
								 pred[np]->cov+=addcov;
								 //fprintf(stderr,"...add cov=%1.f (totalcov=%.1f) to exon[%d] and pred[%d]->cov=%f\n",addcov/no2gnode[transfrag[t]->nodes[i]]->len(),addcov,p,np,pred[np]->cov);
							 }
							 i++;
						 }
					 }
				 }
				 transfrag[t]->abundance=0;
			 }
		 }
		 else {
			 if(!eonly || transfrag[t]->guide) { // store it as an independent prediction
				 pathpat=transfrag[t]->pattern; // not used right now but maybe in the future?
				 path.Clear();
				 path.cAdd(0);
				 path.Add(transfrag[t]->nodes[0]);
				 for(int j=1;j<transfrag[t]->nodes.Count();j++) {
					 if(transfrag[t]->nodes[j]!=1+transfrag[t]->nodes[j-1] ||
							 no2gnode[transfrag[t]->nodes[j]]->start-1!=no2gnode[transfrag[t]->nodes[j-1]]->end) {
						 // check if transfrag t1 is incomplete between node[n-1] and node [n]
						 int *pos=gpos[edge(transfrag[t]->nodes[j-1],transfrag[t]->nodes[j],gno)];
						 if(!pos || !transfrag[t]->pattern[*pos]) { // incomplete transfrag
							 break;
						 }
						 if(pos) pathpat[*pos]=1;
						 path.Add(transfrag[t]->nodes[j]);
					 }
					 else path.Add(transfrag[t]->nodes[j]);
				 }
				 if(path.Last()==transfrag[t]->nodes.Last()) { // this transfrag is complete, might be worth rescuing

					 uint startpoint=no2gnode[path[1]]->end;
					 uint endpoint=no2gnode[path.Last()]->start;
					 if(transfrag[t]->longstart) startpoint=transfrag[t]->longstart;
					 if(transfrag[t]->longend) endpoint=transfrag[t]->longend;

					 int sink=gno-1;
					 path.Add(sink);

					 GVec<GSeg> exons;
					 GVec<float> exoncov;
					 int j=1;
					 int len=0;
					 float cov=0;
					 while(j<path.Count()-1) {
						 int nodestart=no2gnode[path[j]]->start;
						 int nodeend=no2gnode[path[j]]->end;
						 // nodecov[path[j]]-=transfrag[t]->abundance; // do not need this here anymore
						 len+=nodeend-nodestart+1;

						 float excov=0;
						 if(nodecov[path[j]]>epsilon) {
							 excov=transfrag[t]->abundance*no2gnode[path[j]]->len();
							 float newnodecov=nodecov[path[j]]-excov/noderate[path[j]];
							 if(newnodecov<0) {
								 excov=nodecov[path[j]]*noderate[path[j]];
								 newnodecov=0;
							 }
							 nodecov[path[j]]=newnodecov;
						 }

						 while(j+1<path.Count()-1 && no2gnode[path[j]]->end+1==no2gnode[path[j+1]]->start) {
							 j++;
							 len+=no2gnode[path[j]]->len();
							 nodeend=no2gnode[path[j]]->end;

							 float addcov=0;
							 if(nodecov[path[j]]>epsilon) {
								 addcov=transfrag[t]->abundance*no2gnode[path[j]]->len();
								 float newnodecov=nodecov[path[j]]-addcov/noderate[path[j]];
								 if(newnodecov<0) {
									 addcov=nodecov[path[j]]*noderate[path[j]];
									 newnodecov=0;
								 }
								 nodecov[path[j]]=newnodecov;
							 }
							 excov+=addcov;
							 //excov+=transfrag[t]->abundance*noderate[path[j]];
						 }
						 GSeg exon(nodestart,nodeend);
						 exons.Add(exon);
						 cov+=excov;
						 exoncov.Add(excov);
						 j++;
					 }
					 //GffObj *g=NULL;
					 if(len>=mintranscriptlen) {
						 if(first) { geneno++; first=false;}
						 GffObj *g=NULL;
						 if(transfrag[t]->guide) {
							 g=guides[int(transfrag[t]->guide-1)];
							 if (g && g->uptr) {
								 RC_TData &td = *(RC_TData*) (g->uptr);
								 td.in_bundle=3;
								 //fprintf(stderr,"sg guide %s is stored\n",g->getID());
							 }
						 }
						 //fprintf(stderr,"2 Store prediction %d:%d-%d  with len=%d and abundance=%f startpoint=%d endpoint=%d\n",pred.Count(),exons[0].start ,exons.Last().end,len,cov/len,startpoint,endpoint);
						 CPrediction *p=new CPrediction(geneno, g,startpoint , endpoint, cov, sign, len);
						 p->exons=exons;
						 p->exoncov=exoncov;
						 p->tlen=-p->tlen; // negative transcript length signifies assembly is from a long read
						 pred.Add(p);

						 CTransfrag u(path,pathpat,cov/len);
						 keeptrf.Add(u);
					 }
				 }
			 }
			 transfrag[t]->abundance=0;
		 }
	 }


	 if(pred.Count()>npred) {
		 for(int t=0;t<transfrag.Count();t++) if(transfrag[t]->longread) {
			 if(!transfrag[t]->shortread && transfrag[t]->nodes.Count()>1 && transfrag[t]->abundance>epsilon && transfrag[t]->nodes[0] && transfrag[t]->nodes.Last()!=gno-1 ) {
				 /*fprintf(stderr,"Consider transfrag[%d]->abundance=%f with start=%d end=%d nodes:",t,transfrag[t]->abundance,transfrag[t]->longstart,transfrag[t]->longend);
				 for(int i=0;i<transfrag[t]->nodes.Count();i++) fprintf(stderr," %d",transfrag[t]->nodes[i]);
				 fprintf(stderr,"\n");*/

				 GVec<int> tmatch;
				 float abundancesum=best_trf_match(transfrag[t],keeptrf,no2gnode,gno,tmatch); // abundancesum is the sum of all matching transcripts
				 if(abundancesum>0) for(int j=0;j<tmatch.Count();j++){ // found good match(es)
					 float abundprop=transfrag[t]->abundance*keeptrf[tmatch[j]].abundance/abundancesum; // proportion of transcript that will be allocated to this matching one
					 int p=0;
					 int i=0;
					 int np=npred+tmatch[j];
					 //fprintf(stderr,"Add %.1f to prediction %d with cov=%.1f",abundprop,np,pred[np]->cov);
					 //for(int i=0;i<keeptrf[tmatch[j]].nodes.Count();i++) fprintf(stderr," %d",keeptrf[tmatch[j]].nodes[i]);
					 //fprintf(stderr,"\n");
					 while(i<transfrag[t]->nodes.Count() && p<pred[np]->exons.Count()) {
						 if(no2gnode[transfrag[t]->nodes[i]]->end<pred[np]->exons[p].start) i++;
						 else if(pred[np]->exons[p].end<no2gnode[transfrag[t]->nodes[i]]->start) p++;
						 else { // the two intersect (I can only have the full node included in exon)
							 if(nodecov[transfrag[t]->nodes[i]]>epsilon){
								 //float addcov=transfrag[t]->abundance*noderate[transfrag[t]->nodes[i]];
								 float addcov=abundprop*no2gnode[transfrag[t]->nodes[i]]->len();
								 float newnodecov=nodecov[transfrag[t]->nodes[i]]-addcov/noderate[transfrag[t]->nodes[i]];
								 //fprintf(stderr,"newnodecov=%.1f\n",newnodecov);
								 if(newnodecov<0) {
									 newnodecov=nodecov[transfrag[t]->nodes[i]]*keeptrf[tmatch[j]].abundance/abundancesum;
									 addcov=newnodecov*noderate[transfrag[t]->nodes[i]];
									 newnodecov=nodecov[transfrag[t]->nodes[i]]-newnodecov;
								 }
								 nodecov[transfrag[t]->nodes[i]]=newnodecov;
								 pred[np]->exoncov[p]+=addcov;
								 pred[np]->cov+=addcov;
								 //fprintf(stderr,"...add cov=%1.f (totalcov=%.1f) to exon[%d] and pred[%d]->cov=%f\n",addcov/no2gnode[transfrag[t]->nodes[i]]->len(),addcov,p,np,pred[np]->cov);
							 }
							 i++;
						 }
					 }
				 }
			 }
			 if(!transfrag[t]->nodes[0] || transfrag[t]->nodes.Last()==gno-1) // deletes starting transfrangs
				 transfrag[t]->abundance=0;
			 else transfrag[t]->abundance=transfrag[t]->usepath; // this restores tranfrag[t] abundance to what it was before so when I do the short read flux I can utilizes fully

			 transfrag[t]->usepath=-1;
		 }

		 int p=npred;
		 while(p<pred.Count()) {
			 if(pred[p]->cov) {
				 //fprintf(stderr,"Pred[%d] has coverage=%.1f (totalcov=%.1f)\n",p,pred[p]->cov,pred[p]->cov/abs(pred[p]->tlen));
				 pred[p]->cov/=abs(pred[p]->tlen);
				 for(int i=0;i<pred[p]->exons.Count();i++)
					 pred[p]->exoncov[i]/=pred[p]->exons[i].len();
				 // adjust start/endpoints
				 if(pred[p]->start!=pred[p]->exons[0].start) {
					 pred[p]->tlen+=pred[p]->start-pred[p]->exons[0].start;
					 pred[p]->exons[0].start=pred[p]->start;
				 }
				 if(pred[p]->end!=pred[p]->exons.Last().end) {
					 pred[p]->tlen+=pred[p]->exons.Last().end-pred[p]->end;
					 pred[p]->exons.Last().end=pred[p]->end;
				 }
				 p++;
			 }
			 else if(!eonly) { // || !pred[p]->t_eq) {
				 //fprintf(stderr,"delete prediction %d\n",p);
				 pred.Delete(p); // I delete all predictions that have 0 coverage unless it's eonly mode
			 }
			 else p++;
		 }
	 }

	 int nkeep=keeptrf.Count();
	 /*for(int i=nkeep-1;i>nkept-1;i--) { // delete all kept transfrags we are not confident in
			 keeptrf.Delete(i);
	 }*/
	 for(int i=0;i<nkeep;i++) keeptrf[i].weak=i+npred; // I need to remember the prediction it represents
	 keeptrf.Sort(longtrCmp); // most abundant transfrag in the graph comes first, then the one with most nodes, then the one more complete

	 float flux=0;
	 GVec<float> nodeflux;
	 path.Clear();


	 for(int i=0;i<nkeep;i++) { // compute flux from short read data here
		 istranscript.reset();
		 nodeflux.Clear();
		 bool full=true;
		 // make sure source and sink are present and they are connected in the pattern of the prediction
		 if(keeptrf[i].nodes[0]) {
			 keeptrf[i].nodes.Insert(0,0);
			 keeptrf[i].pattern[0]=1;
		 }
		 int sink=gno-1;
		 if(keeptrf[i].nodes.Last()!=gno-1) {
			 keeptrf[i].nodes.Add(sink);
			 keeptrf[i].pattern[sink]=1;
		 }
		 int key=edge(0,keeptrf[i].nodes[1],gno);
		 int *pos=NULL;
		 if(key<(int)gpos.Count()) pos=gpos[key];
		 if(pos!=NULL) keeptrf[i].pattern[*pos]=1;
		 pos=NULL;
		 int n=keeptrf[i].nodes.Count();
		 key=edge(keeptrf[i].nodes[n-2],sink,gno);
		 if(pos!=NULL) keeptrf[i].pattern[*pos]=1;

		 flux=push_max_flow(gno,keeptrf[i].nodes,istranscript,transfrag,no2gnode,nodeflux,keeptrf[i].pattern,gpos,full);

		 /*
		 { // DEBUG ONLY
			 //printTime(stderr);
			 fprintf(stderr,"flux=%g Path:",flux);
			 for(int j=0;j<keeptrf[i].nodes.Count();j++) fprintf(stderr," %d",keeptrf[i].nodes[j]);
			 fprintf(stderr,"***\n");
		 }
		 */

		 if(flux>epsilon) {

			 bool included=true;
			 float cov=store_transcript(pred,keeptrf[i].nodes,nodeflux,nodecovall,no2gnode,geneno,first,strand,gno,gpos,included,prevpath);
			 if(cov) {

				 //fprintf(stderr,"Coverage of long pred[%d]=%f vs. short pred=%f\n",keeptrf[i].weak,pred[keeptrf[i].weak]->cov,pred.Last()->cov);

				 if(pred.Last()->cov>pred[keeptrf[i].weak]->cov) { // new prediction is better than the previous one -- shouldn't I add this to previous prediction?
					 int p=keeptrf[i].weak;
					 pred[p]->cov=pred.Last()->cov;
					 for(int k=0;k<pred[p]->exons.Count();k++) {
						 pred[p]->exoncov[k]=pred.Last()->exoncov[k];
					 }
					 /*pred[p]->start=pred.Last()->start;
					 pred[p]->end=pred.Last()->end;
					 pred[p]->exons[0].start=pred.Last()->exons[0].start;
					 pred[p]->exons.Last().end=pred.Last()->exons.Last().end;
					 pred[p]->tlen=pred.Last()->tlen;*/
					 //pred[p]->tlen=abs(pred[p]->tlen);
				 }
				 delete pred.Pop(); //prevent memory leak, popped element is otherwise "forgotten"
			 }
		 }

	 }

}

void get_trf_long(int gno,int edgeno, GIntHash<int> &gpos,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,
		int& geneno,int strand,GList<CPrediction>& pred,GVec<int>& trflong,BundleData *bdata) {

	GPVec<GffObj>& guides = bdata->keepguides;

	GVec<float> nodecov; // the coverage of all transfrags entering a node
	GVec<float> noderate;
	for(int i=0;i<gno;i++) {
		CGraphnode *inode=no2gnode[i]; // this is here only because of the DEBUG option below
		nodecov.cAdd(0.0); // this are all transfrags that link nodes together
		float rate=1;
		if(i && i<gno-1) {
			for(int j=0;j<inode->trf.Count();j++) { // for all transfrags going through node
				int t=inode->trf[j];
				//if(transfrag[t]->nodes[0]<i) { // entering transfrags:
				if(transfrag[t]->nodes.Last()>i) { // exiting transfrags: this is more consistent with the nodeflux computation
					nodecov[i]+=transfrag[t]->abundance;
				}
			}

			if(nodecov[i]) rate=nodecov[i];
			if(rate<=0) rate=1; // this shouldn't happen
			//fprintf(stderr,"rate=%f\n",rate);
			rate=inode->cov/rate;
		}
		noderate.Add(rate);
		//fprintf(stderr,"Node[%d]:%d-%d no2gnode->cov=%f nodecov=%f noderate=%f\n",i,inode->start,inode->end,inode->cov,nodecov[i],noderate[i]);
	}

	GBitVec istranscript(transfrag.Count());

	 GVec<int> path;
	 GBitVec pathpat(gno+edgeno);
	 int minp;
	 int maxp;
	 int maxi;

	 char sign='-';
	 if(strand) { sign='+';}
	 int npred=pred.Count();

	 GVec<CTransfrag> keeptrf;
	 GVec<int> checktrf;
	 for(int f=trflong.Count()-1;f>=0;f--) { // if this is a guide it should be reflected in the prediction downstream
		 path.Clear();
		 int t=trflong[f];
		 if(t<0) GError("Stored long transcript is negative!\n");
		 pathpat=transfrag[t]->pattern;
		 minp=transfrag[t]->nodes[0];
		 maxp=transfrag[t]->nodes.Last();

		 //if(no2gnode[transfrag[t]->nodes[0]]->hardstart) {
			 int *pos=gpos[edge(0,minp,gno)];
			 if(pos) pathpat[*pos]=1;
			 //minp=0;
		 //}
		 //if(no2gnode[transfrag[t]->nodes.Last()]->hardend) {
			 pos=gpos[edge(maxp,gno-1,gno)];
			 if(pos) pathpat[*pos]=1;
			 //maxp=gno-1;
		 //}

		 maxi=minp;
		 path.Add(maxi);
		 pathpat[maxi]=1;

		 istranscript.reset();

		 float flux=0;
		 //float fragno=0;
		 GVec<float> nodeflux;

		 /*
	 	 { // DEBUG ONLY
	 	 fprintf(stderr,"\n\n***Start get_trf_long with maxi=%d minp=%d maxp=%d guide=%d and transcript:",maxi,minp,maxp,transfrag[t]->guide);
	 	 for(int i=0;i<transfrag[t]->nodes.Count();i++) {
	 		 if(i) {
	 			 pos=gpos[edge(transfrag[t]->nodes[i-1],transfrag[t]->nodes[i],gno)];
	 			 if(pos && pathpat[*pos])
	 				 fprintf(stderr,"-");
	 		 }
	 		 fprintf(stderr," %d",transfrag[t]->nodes[i]);

	 	 }
	 	 fprintf(stderr," pathpat=");
	 	 //printBitVec(pathpat);
		 fprintf(stderr,"\n");

#ifdef GMEMTRACE
	 	 double vm,rsm;
	 	 get_mem_usage(vm, rsm);
	 	 GMessage("\t\tM(s):parse_trf memory usage: rsm=%6.1fMB vm=%6.1fMB\n",rsm/1024,vm/1024);
#endif
	 	 }
	 	 */

		 bool tocheck=true;
		 if(back_to_source_fast_long(maxi,path,minp,maxp,pathpat,transfrag,no2gnode,nodecov,gno,gpos)) {
			 path.cAdd(0);
			 path.Reverse(); // back to source adds the nodes at the end to avoid pushing the list all the time

			 if(fwd_to_sink_fast_long(maxi,path,minp,maxp,pathpat,transfrag,no2gnode,nodecov,gno,gpos)) {

				 flux=long_max_flow(gno,path,istranscript,transfrag,no2gnode,nodeflux,pathpat);

				 /*
				 { // DEBUG ONLY
					 //printTime(stderr);
					 fprintf(stderr,"flux=%g Path:",flux);
					 for(int i=0;i<path.Count();i++) fprintf(stderr," %d",path[i]);
					 fprintf(stderr,"\n");
					 fprintf(stderr,"Nodecapacities:");
					 for(int i=0;i<path.Count();i++) fprintf(stderr," %f",nodeflux[i]);
					 fprintf(stderr,"***\n");
				 }
				 */

				 if(flux) { // these are not valid paths in the graph
					 tocheck=false;

					 GVec<GSeg> exons;
					 GVec<float> exoncov;
					 int j=1;
					 int len=0;
					 float cov=0;
					 //int startnode=j;
					 int lastnode=path.Count()-2;
					 //uint startpoint=no2gnode[path[lastnode]]->end;
					 //uint endpoint=no2gnode[path[1]]->start;

					 /* if(mixedMode) { // establish start/end point of path
						 while(j<path.Count()-1) {
							 CGraphnode *jnode=no2gnode[path[j]];
							 for(int i=0;i<jnode->trf.Count();i++) { // for all transfrags going through node
								 int t=jnode->trf[i];
								 if(istranscript[t] && transfrag[t]->longread && transfrag[t]->nodes[0] && transfrag[t]->nodes.Last()!=gno-1) {
									 if(transfrag[t]->nodes[0]==path[j]) {
										 if(transfrag[t]->longstart) {
											 if(transfrag[t]->longstart<startpoint) { startpoint=transfrag[t]->longstart; startnode=j;}
										 }
										 else if(jnode->start<startpoint) { startpoint=jnode->start; startnode=j;}
									 }
									 if(transfrag[t]->nodes.Last()==path[j]) {
										 if(transfrag[t]->longend) {
											 if(transfrag[t]->longend>endpoint) { endpoint=transfrag[t]->longend; lastnode=j;}
										 }
										 else if(jnode->end>endpoint) {endpoint=jnode->end; lastnode=j;}
									 }
								 }
							 }
							 j++;
						 }
						 if(startpoint>endpoint) {
							 j=1;
							 lastnode=path.Count()-2;
							 startpoint=no2gnode[path[1]]->start;
							 endpoint=no2gnode[path[lastnode]]->end;
						 }
						 else j=startnode;
						 //fprintf(stderr,"startnode=%d lastnode=%d startpoint=%d endpoint=%d\n",path[startnode],path[lastnode],startpoint,endpoint);
					 }*/

					 while(j<=lastnode) {
					 //while(j<path.Count()-1) {
						 int nodestart=no2gnode[path[j]]->start;
						 int nodeend=no2gnode[path[j]]->end;
						 nodecov[path[j]]-=nodeflux[j];
						 len+=nodeend-nodestart+1;
						 float ecov=nodeflux[j]*noderate[path[j]];
						 float excov=ecov;
						 //fprintf(stderr,"excov+=%f * %f = %f\n",nodeflux[j],noderate[path[j]],excov);
						 /*if(mixedMode) {
							 no2gnode[path[j]]->cov-=ecov;
							 if(no2gnode[path[j]]->cov<epsilon) no2gnode[path[j]]->cov=0;
						 }*/
						 //float excov=nodeflux[j]*noderate[path[j]];
						 while(j+1<=lastnode && no2gnode[path[j]]->end+1==no2gnode[path[j+1]]->start) {
						 //while(j+1<path.Count()-1 && no2gnode[path[j]]->end+1==no2gnode[path[j+1]]->start) {
							 j++;
							 nodecov[path[j]]-=nodeflux[j];
							 nodeend=no2gnode[path[j]]->end;
							 ecov=nodeflux[j]*noderate[path[j]];
							 excov+=ecov;
							 //fprintf(stderr,"excov+=%f * %f = %f\n",nodeflux[j],noderate[path[j]],excov);
							 len+=nodeend-no2gnode[path[j]]->start+1;

						 }
						 GSeg exon(nodestart,nodeend);
						 exons.Add(exon);
						 //fprintf(stderr,"excov=%f\n",excov/(exon.end-exon.start+1));
						 cov+=excov;
						 //fprintf(stderr,"cov+=%f=%f\n",excov,cov);
						 exoncov.Add(excov);
						 j++;
					 }
					 if(transfrag[t]->nodes.Count()==1) transfrag[t]->abundance=0;
					 //fprintf(stderr,"Store prediction %d  with abundance=%f len=%d\n",pred.Count(),cov/len,len);
					 GffObj *g=NULL;
					 if(transfrag[t]->guide) {
						 g=guides[int(transfrag[t]->guide-1)];
						 if (g && g->uptr) {
							 RC_TData &td = *(RC_TData*) (g->uptr);
							 td.in_bundle=3;
							 //fprintf(stderr,"sg guide %s is stored\n",g->getID());
						 }
					 }
					 if(!eonly || g) {
						 /*fprintf(stderr,"1 Store prediction %d  with abundance=%f len=%d and exons:",pred.Count(),cov/len,len);
						 for(int i=0;i<exons.Count();i++) fprintf(stderr," %d-%d",exons[i].start,exons[i].end);
						 fprintf(stderr,"\n");*/
						 CPrediction *p=new CPrediction(geneno, g,exons[0].start , exons.Last().end, cov, sign, len);
						 p->exons=exons;
						 p->exoncov=exoncov;
						 p->mergename='.'; // I should not delete this prediction
						 p->tlen=-p->tlen; // negative transcript length signifies assembly is from a long read
						 pred.Add(p);

						 CTransfrag u(path,pathpat,cov/len);
						 keeptrf.Add(u);
					 }
				 }
				 else if(transfrag[t]->guide) {
					 checktrf.Add(t);
				 }
			 }
		 }

		 if(tocheck)  { // try to see if you can rescue transfrag -> they are stored from more abundant to least -> not if using mixedMode
			 checktrf.Add(t);
		 }
	 }

	 //keeptrf.Sort(longtrCmp); // most abundant transfrag in the graph come first, then the ones with most nodes, then the ones more complete

	 for(int c=0;c<checktrf.Count();c++) if(transfrag[checktrf[c]]->guide || transfrag[checktrf[c]]->abundance>=readthr) { // only in this case it is worth considering it as a potential prediction
		 int t=checktrf[c];

		 //fprintf(stderr,"checktrf[%d]=%d with abundance=%f\n",c,t,transfrag[t]->abundance);

		 GVec<int> tmatch;
		 float abundancesum=best_trf_match(transfrag[t],keeptrf,no2gnode,gno,tmatch);
		 if(abundancesum>0) for(int j=0;j<tmatch.Count();j++){ // found good match(es)
			 float abundprop=transfrag[t]->abundance*keeptrf[tmatch[j]].abundance/abundancesum;
			 int p=0;
			 int i=0;
			 int np=npred+tmatch[j]; // how do I know that keeptrf lead to a prediction -> because keeptrf represent the predictions that were added
			 while(i<transfrag[t]->nodes.Count() && p<pred[np]->exons.Count()) {
				 if(no2gnode[transfrag[t]->nodes[i]]->end<pred[np]->exons[p].start) i++;
				 else if(pred[np]->exons[p].end<no2gnode[transfrag[t]->nodes[i]]->start) p++;
				 else { // the two intersect (I can only have the full node included in exon)
					 if(nodecov[transfrag[t]->nodes[i]]>epsilon){
						 //float addcov=transfrag[t]->abundance*noderate[transfrag[t]->nodes[i]];
						 float addcov=abundprop*no2gnode[transfrag[t]->nodes[i]]->len();
						 float newnodecov=nodecov[transfrag[t]->nodes[i]]-addcov/noderate[transfrag[t]->nodes[i]];
						 if(newnodecov<0) {
							 newnodecov=nodecov[transfrag[t]->nodes[i]]*keeptrf[tmatch[j]].abundance/abundancesum;
							 addcov=newnodecov*noderate[transfrag[t]->nodes[i]];
							 newnodecov=nodecov[transfrag[t]->nodes[i]]-newnodecov;
						 }
						 nodecov[transfrag[t]->nodes[i]]=newnodecov;

						 pred[np]->exoncov[p]+=addcov;
						 pred[np]->cov+=addcov;
					 }
					 i++;
				 }
			 }
		 }
		 else if(!eonly || transfrag[t]->guide) { // store it as an independent prediction
			 pathpat=transfrag[t]->pattern; // not used right now but maybe in the future?
			 path.Clear();
			 path.cAdd(0);
			 path.Add(transfrag[t]->nodes[0]);
			 for(int j=1;j<transfrag[t]->nodes.Count();j++) {
				 if(transfrag[t]->nodes[j]!=1+transfrag[t]->nodes[j-1] ||
						 no2gnode[transfrag[t]->nodes[j]]->start-1!=no2gnode[transfrag[t]->nodes[j-1]]->end) {
					 // check if transfrag t1 is incomplete between node[n-1] and node [n]
					 int *pos=gpos[edge(transfrag[t]->nodes[j-1],transfrag[t]->nodes[j],gno)];
					 if(!pos || !transfrag[t]->pattern[*pos]) { // incomplete transfrag
						 break;
					 }
					 if(pos) pathpat[*pos]=1;
					 path.Add(transfrag[t]->nodes[j]);
				 }
				 else path.Add(transfrag[t]->nodes[j]);
			 }
			 if(path.Last()==transfrag[t]->nodes.Last()) { // this transfrag is complete, might be worth rescuing
				 int sink=gno-1;
				 path.Add(sink);

				 GVec<GSeg> exons;
				 GVec<float> exoncov;
				 int j=1;
				 int len=0;
				 float cov=0;
				 while(j<path.Count()-1) {
					 int nodestart=no2gnode[path[j]]->start;
					 int nodeend=no2gnode[path[j]]->end;
					 // nodecov[path[j]]-=transfrag[t]->abundance; // do not need this here anymore
					 len+=nodeend-nodestart+1;

					 float excov=0;
					 if(nodecov[path[j]]>epsilon) {
						 excov=transfrag[t]->abundance*no2gnode[path[j]]->len();
						 float newnodecov=nodecov[path[j]]-excov/noderate[path[j]];
						 if(newnodecov<0) {
							 excov=nodecov[path[j]]*noderate[path[j]];
							 newnodecov=0;
						 }
						 nodecov[path[j]]=newnodecov;
					 }

					 while(j+1<path.Count()-1 && no2gnode[path[j]]->end+1==no2gnode[path[j+1]]->start) {
						 j++;
						 len+=no2gnode[path[j]]->len();
						 nodeend=no2gnode[path[j]]->end;

						 float addcov=0;
						 if(nodecov[path[j]]>epsilon) {
							 addcov=transfrag[t]->abundance*no2gnode[path[j]]->len();
							 float newnodecov=nodecov[path[j]]-addcov/noderate[path[j]];
							 if(newnodecov<0) {
								 addcov=nodecov[path[j]]*noderate[path[j]];
								 newnodecov=0;
							 }
							 nodecov[path[j]]=newnodecov;
						 }
						 excov+=addcov;
						 //excov+=transfrag[t]->abundance*noderate[path[j]];
					 }
					 GSeg exon(nodestart,nodeend);
					 exons.Add(exon);
					 cov+=excov;
					 exoncov.Add(excov);
					 j++;
				 }
				 GffObj *g=NULL;
				 if(transfrag[t]->guide) {
					 g=guides[int(transfrag[t]->guide-1)];
					 if (g && g->uptr) {
						 RC_TData &td = *(RC_TData*) (g->uptr);
						 td.in_bundle=3;
						 //fprintf(stderr,"sg guide %s is stored\n",g->getID());
					 }
				 }
				 //fprintf(stderr,"2 Store prediction %d:%d-%d  with len=%d and abundance=%f\n",pred.Count(),exons[0].start ,exons.Last().end,len,cov/len);
				 CPrediction *p=new CPrediction(geneno, g,exons[0].start , exons.Last().end, cov, sign, len);
				 p->exons=exons;
				 p->exoncov=exoncov;
				 p->tlen=-p->tlen; // negative transcript length signifies assembly is from a long read
				 pred.Add(p);

				 CTransfrag u(path,pathpat,cov/len);
				 keeptrf.Add(u);
			 }
		 }
		 transfrag[t]->abundance=0;
	 }

	 if(pred.Count()>npred) {
		 /*if(mixedMode) {
			 for(int t=0;t<transfrag.Count();t++) {
			   if(transfrag[t]->longread && (!transfrag[t]->nodes[0] || transfrag[t]->nodes.Last()==gno-1) )
			     transfrag[t]->abundance=0;
			   else
			     transfrag[t]->abundance=transfrag[t]->usepath;
			   transfrag[t]->usepath=-1;
			 }
		 }
		 else */for(int t=0;t<transfrag.Count();t++) if(transfrag[t]->longread) { // longreads mode --> tries to add all transfrags to predictions
			 if(transfrag[t]->abundance>epsilon && transfrag[t]->nodes[0] && transfrag[t]->nodes.Last()!=gno-1 ) {
			   /*fprintf(stderr,"Consider transfrag[%d]->abundance=%f with start=%d end=%d nodes:",t,transfrag[t]->abundance,transfrag[t]->longstart,transfrag[t]->longend);
				 for(int i=0;i<transfrag[t]->nodes.Count();i++) fprintf(stderr," %d",transfrag[t]->nodes[i]);
				 fprintf(stderr,"\n");*/

				 GVec<int> tmatch;
				 float abundancesum=best_trf_match(transfrag[t],keeptrf,no2gnode,gno,tmatch);
				 if(abundancesum>0) for(int j=0;j<tmatch.Count();j++){ // found good match(es)
					 float abundprop=transfrag[t]->abundance*keeptrf[tmatch[j]].abundance/abundancesum;
					 //fprintf(stderr,"Added %.1f to prediction:",abundprop);
					 //for(int i=0;i<keeptrf[tmatch[j]].nodes.Count();i++) fprintf(stderr," %d",keeptrf[tmatch[j]].nodes[i]);
					 int p=0;
					 int i=0;
					 int np=npred+tmatch[j];
					 while(i<transfrag[t]->nodes.Count() && p<pred[np]->exons.Count()) {
						 if(no2gnode[transfrag[t]->nodes[i]]->end<pred[np]->exons[p].start) i++;
						 else if(pred[np]->exons[p].end<no2gnode[transfrag[t]->nodes[i]]->start) p++;
						 else { // the two intersect (I can only have the full node included in exon)
							 if(nodecov[transfrag[t]->nodes[i]]>epsilon){
								 //float addcov=transfrag[t]->abundance*noderate[transfrag[t]->nodes[i]];
								 float addcov=abundprop*no2gnode[transfrag[t]->nodes[i]]->len();
								 float newnodecov=nodecov[transfrag[t]->nodes[i]]-addcov/noderate[transfrag[t]->nodes[i]];
								 if(newnodecov<0) {
									 newnodecov=nodecov[transfrag[t]->nodes[i]]*keeptrf[tmatch[j]].abundance/abundancesum;
									 addcov=newnodecov*noderate[transfrag[t]->nodes[i]];
									 newnodecov=nodecov[transfrag[t]->nodes[i]]-newnodecov;
								 }
								 nodecov[transfrag[t]->nodes[i]]=newnodecov;
								 //fprintf(stderr," newnodecov=%f\n",newnodecov);
								 pred[np]->exoncov[p]+=addcov;
								 pred[np]->cov+=addcov;
							 }
							 i++;
						 }
					 }
					 //fprintf(stderr," new abundance=%f\n",pred[np]->cov/pred[np]->tlen);
				 }
			 }
			 transfrag[t]->abundance=0; // delete abundance in order not to use it in short reads
		 }

		 int p=npred;
		 while(p<pred.Count()) {
			 if(pred[p]->cov) {
				 pred[p]->cov/=abs(pred[p]->tlen);
				 for(int i=0;i<pred[p]->exons.Count();i++)
					 pred[p]->exoncov[i]/=pred[p]->exons[i].len();
				 p++;
			 }
			 else if(!eonly) { // || !pred[p]->t_eq) {
				 //fprintf(stderr,"delete prediction %d\n",p);
				 pred.Delete(p); // I delete all predictions that have 0 coverage unless it's eonly mode
			 }
			 else p++;
		 }
	 }
}

/*
void parse_trf_long(int maxi,int gno,int edgeno, GIntHash<int> &gpos,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,
		int& geneno,bool first,int strand,GList<CPrediction>& pred,GVec<float>& nodecov,
		GBitVec& istranscript,GBitVec& removable,GBitVec& usednode,float maxcov,GBitVec& prevpath,GVec<int>& trflong) {

	 GVec<int> path;
	 GBitVec pathpat(gno+edgeno);
	 int minp=maxi;
	 int maxp=maxi;

	 bool usetrflong=false;

	 if(trflong.Count()) {
		 int t=trflong.Pop();
		 usetrflong=true;
		 pathpat=pathpat | transfrag[t]->pattern;
		 minp=transfrag[t]->nodes[0];
		 maxp=transfrag[t]->nodes.Last();
		 int *pos=gpos[edge(0,minp,gno)];
		 if(pos) pathpat[*pos]=1;
		 pos=gpos[edge(maxp,gno-1,gno)];
		 if(pos) pathpat[*pos]=1;
		 maxi=minp;
	 }

	 path.Add(maxi);
	 pathpat[maxi]=1;


	 bool full=false;
	 istranscript.reset();

	 float flux=0;
	 //float fragno=0;
	 GVec<float> nodeflux;


	 if(back_to_source_fast_long(maxi,path,minp,maxp,pathpat,transfrag,no2gnode,nodecov,gno,gpos)) {
		 if(includesource) path.cAdd(0);
		 path.Reverse(); // back to source adds the nodes at the end to avoid pushing the list all the time

		 if(fwd_to_sink_fast_long(maxi,path,minp,maxp,pathpat,transfrag,no2gnode,nodecov,gno,gpos)) {

			 flux=push_max_flow(gno,path,istranscript,transfrag,no2gnode,nodeflux,pathpat,gpos,full);

	 		}

	 }


	 if(flux>epsilon) {
		 bool included=true;
		 float cov=store_transcript(pred,path,nodeflux,nodecov,no2gnode,geneno,first,strand,gno,gpos,included,prevpath,full);


		 if(!usetrflong && (included || cov<isofrac*maxcov)) { // maybe not use it for long reads if
			 usednode[maxi]=1;
			 maxi=0;
			 maxcov=0;
			 //cont=false;

		 }
		 else if(cov>maxcov) maxcov=cov;
	 }
	 else {

		 usednode[maxi]=1;

		 maxi=0;
		 maxcov=0;
		 //cont=false;
	 }

	 // Node coverages:
	 for(int i=1;i<gno;i++) {
		 if(!usednode[i] && nodecov[i]>nodecov[maxi]) maxi=i;
		 if(i<gno-1 && no2gnode[i]->end==no2gnode[i+1]->start) {
			 for(int t=0;t<no2gnode[i]->trf.Count();t++) if(transfrag[t]->weak>0 && transfrag[t]->pattern[i+1]) transfrag[t]->weak=-1;
		 }
	 }


	 //if(nodecov[maxi]>=readthr && (!specific || cont)) { // if I still have nodes that are above coverage threshold
	 if(nodecov[maxi]>=1) { // if I still have nodes to use

		 path.Clear();
		 nodeflux.Clear();
		 parse_trf_long(maxi,gno,edgeno,gpos,no2gnode,transfrag,geneno,first,strand,pred,nodecov,istranscript,removable,usednode,maxcov,prevpath,trflong);
	 }

}
*/

void parse_trf(int maxi,int gno,int edgeno, GIntHash<int> &gpos,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,
		int& geneno,bool first,int strand,GList<CPrediction>& pred,GVec<float>& nodecov,
		GBitVec& istranscript,GBitVec& usednode,float maxcov,GBitVec& prevpath) {

	 GVec<int> path;
	 path.Add(maxi);
	 GBitVec pathpat(gno+edgeno);
	 pathpat[maxi]=1;
	 istranscript.reset();

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


	if(back_to_source_fast(maxi,path,pathpat,transfrag,no2gnode,nodecov,gno,gpos)) {
		 	 if(includesource) path.cAdd(0);
	 		 path.Reverse(); // back to source adds the nodes at the end to avoid pushing the list all the time

			if(fwd_to_sink_fast(maxi,path,pathpat,transfrag,no2gnode,nodecov,gno,gpos)) {
	 			 bool full=true;

				 flux=push_max_flow(gno,path,istranscript,transfrag,no2gnode,nodeflux,pathpat,gpos,full);

				 /*
	 			 { // DEBUG ONLY
	 				 //printTime(stderr);
	 				 fprintf(stderr,"flux=%g Path:",flux);
	 				 for(int i=0;i<path.Count();i++) fprintf(stderr," %d",path[i]);
	 				 fprintf(stderr,"***\n");
	 			 }
	 			 */
	 		}
			/*else {
	 			//pathpat.reset();
	 			//pathpat[maxi]=1;
	 		}*/
	 }
	 /*else {
		 //pathpat.reset();
		 //pathpat[maxi]=1;
	 }*/

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

		 float frac=isofrac;
		 if(mixedMode && frac<ERROR_PERC) frac=ERROR_PERC;

		 if(included || cov<frac*maxcov) {
			 /*
			 if(sensitivitylevel) usednode[maxi]=1;
			 else usednode = usednode | prevpath;
			 */
			 //fprintf(stderr,"included\n");
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

	 //fprintf(stderr,"maxcov=%g gno=%d\n",maxcov,gno);

	 // Node coverages:
	 for(int i=1;i<gno;i++)
		 if(!usednode[i] && nodecov[i]>nodecov[maxi]) maxi=i;

	 //fprintf(stderr," maxi=%d nodecov=%f\n",maxi,nodecov[maxi]);

	 //if(nodecov[maxi]>=readthr && (!specific || cont)) { // if I still have nodes that are above coverage threshold
	 if(nodecov[maxi]>=1) { // if I still have nodes to use

		 /*
		 { // DEBUG ONLY
			 //printTime(stderr);
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
		 parse_trf(maxi,gno,edgeno,gpos,no2gnode,transfrag,geneno,first,strand,pred,nodecov,istranscript,usednode,maxcov,prevpath);
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
				//else GError("Found parent-child %d-%d not linked by edge\n",i,inode->child[k]);
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
						//else GError("Found parent-child %d-%d not linked by edge\n",i,i+1);
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
		GPVec<CTransfrag>& transfrag,int s,GVec<CGuide>& guidetrf,BundleData *bdata) {

	GPVec<GffObj>& guides = bdata->keepguides;

	//fprintf(stderr,"In refguides with %d guides\n",guides.Count());

	char strand='-';
	if(s) strand='+';

	// find guides' patterns
	for(int g=0;g<guides.Count();g++) {
		//fprintf(stderr,"Consider guide[%d out of %d] %s in_bundle=%d\n",g,guides.Count(),guides[g]->getID(),((RC_TData*)(guides[g]->uptr))->in_bundle);
		if((guides[g]->strand==strand || guides[g]->strand=='.') && ((RC_TData*)(guides[g]->uptr))->in_bundle>=2 && (guides[g]->overlap(no2gnode[1]->start,no2gnode[gno-2]->end))) {
			CTransfrag *trguide=find_guide_pat(guides[g],no2gnode,gno,edgeno,gpos);
			if(trguide) { // the guide can be found among the graph nodes
				/*if(longreads) { // do not allow guide to be too far away from start/end of the transfrag
					if(abs((int)(no2gnode[trguide->nodes[0]]->start-guides[g]->start))<CHI_WIN+CHI_THR &&  abs((int)(guides[g]->end-no2gnode[trguide->nodes.Last()]->end))<CHI_WIN+CHI_THR) {
						CGuide newguide(trguide,g);
						guidetrf.Add(newguide);
						trguide->guide=true;
					}
					else delete trguide;
				}
				else {*/
					//CGuide newguide(trguide,guides[g]);
					CGuide newguide(trguide,g);
					//if(!longreads)
					trguide->real=true;
					trguide->guide=1;
					guidetrf.Add(newguide);
				//}

				/*
				{ // DEBUG ONLY
					if(trguide) {
						fprintf(stderr,"Added guidetrf %d (%d) with ID=%s overlapping transcript interval %d - %d with nodes:",g,guidetrf.Count()-1,guides[g]->getID(),no2gnode[1]->start,no2gnode[gno-2]->end);
						for(int i=0;i<trguide->nodes.Count();i++) fprintf(stderr," %d",trguide->nodes[i]);
						fprintf(stderr,"\n");
						//fprintf(stderr,"s=%d strand = %c trguide[%d]=",s,strand,g);
						//printBitVec(trguide->pattern);
						//fprintf(stderr,"\n");
					}
				}
				*/
			}
		}
	} //for g

	/*if(longreads && guidetrf.Count()) { // if a guide is not matching any transcript the best do not include TODO: rewrite this!!
		int ng=guidetrf.Count();
		GVec<bool> covered(ng,false);
		GVec<int> gcount; // count number of bits set in each guide's pattern
		for(int g=0;g<guidetrf.Count();g++) {
			int gc=guidetrf[g].trf->pattern.count();
			gcount.Add(gc);
		}
		for(int t=0;t<transfrag.Count();t++) if(transfrag[t]->nodes[0] && transfrag[t]->nodes.Last()!=gno-1){
			int maxg=-1;
			int maxcovscore=0;
			//int mingscore=0;
			for(int g=0;g<guidetrf.Count();g++) {
				GBitVec intersect=transfrag[t]->pattern & guidetrf[g].trf->pattern;
				int covscore=intersect.count();
				if(covscore>maxcovscore) { // only allows guides that have something in common
					maxg=g;
					maxcovscore=covscore;
					//mingscore=gcount[g]-covscore;
				}
				else if(covscore==maxcovscore){
					int gcov=gcount[g]-covscore;
					if(gcov<mingscore) {
						mingscore=gcov;
						maxg=g;
					}
					else if(gcov==mingscore && !covered[maxg] && covered[g]) {
						maxg=g;
					}
				}
			}
			if(maxg>=0 && !covered[maxg]) {
				covered[maxg]=true;
				ng--;
				if(!ng) break;
			}
		}
		if(ng) {
			for(int g=0;g<guidetrf.Count();g++) if(!covered[g]) {
				delete guidetrf[g].trf;
				guidetrf.Delete(g);
			}
		}
	}*/

	// print guide coverages
	if(c_out) {
		for(int g=0;g<guidetrf.Count();g++) {
			GBitVec guidepattern;
			float guideabundance=0;
			for(int t=0;t<transfrag.Count();t++) { // check compatibility with guides
				int tstart=0;
				GBitVec trpattern(transfrag[t]->pattern);
				bool edgecompatible=true;
				if(guidetrf[g].trf->nodes[0]>transfrag[t]->nodes[0]) { // transfrag starts before guide -> see if I can adjust trpattern to the guide limit
					tstart=1;
					while(tstart<transfrag[t]->nodes.Count() && transfrag[t]->nodes[tstart-1]<guidetrf[g].trf->nodes[0]) {
						if(transfrag[t]->nodes[tstart]==transfrag[t]->nodes[tstart-1]+1 && 1+no2gnode[transfrag[t]->nodes[tstart-1]]->end==no2gnode[transfrag[t]->nodes[tstart]]->start) {
							trpattern[transfrag[t]->nodes[tstart-1]]=0;
							int *pos=gpos[edge(transfrag[t]->nodes[tstart-1],transfrag[t]->nodes[tstart],gno)];
							if(pos) trpattern[*pos]=0;
						}
						else {
							edgecompatible=false;
							break;
						}
						tstart++;
					}
				}
				int tend=transfrag[t]->nodes.Count()-1;
				if(edgecompatible && guidetrf[g].trf->nodes.Last()<transfrag[t]->nodes.Last()) { // transfrag starts before guide -> see if I can adjust trpattern to the guide limit
					tend--;
					while(tend>=0 && transfrag[t]->nodes[tend+1]>guidetrf[g].trf->nodes.Last()) {
						if(transfrag[t]->nodes[tend]==transfrag[t]->nodes[tend+1]-1 && 1+no2gnode[transfrag[t]->nodes[tend]]->end==no2gnode[transfrag[t]->nodes[tend+1]]->start) {
							trpattern[transfrag[t]->nodes[tend+1]]=0;
							int *pos=gpos[edge(transfrag[t]->nodes[tend],transfrag[t]->nodes[tend+1],gno)];
							if(pos) trpattern[*pos]=0;
						}
						else {
							edgecompatible=false;
							break;
						}
						tend--;
					}
				}

				if(edgecompatible && ((trpattern & guidetrf[g].trf->pattern) == trpattern)) { // transfrag is completely included in guide
					guidepattern = guidepattern | trpattern;
					// compute guideabundance here
					for(int i=tstart; i<=tend;i++) guideabundance+=transfrag[t]->abundance*no2gnode[transfrag[t]->nodes[i]]->len();
				}
			}
			if(guidepattern==guidetrf[g].trf->pattern) {
				int len=0;
				for(int i=0;i<guidetrf[g].trf->nodes.Count();i++) len+=no2gnode[guidetrf[g].trf->nodes[i]]->len();
				guideabundance/=len;
				GStr guidecov;
				guidecov.appendfmt("%.2f",guideabundance);
				guides[guidetrf[g].g]->addAttr("coverage",guidecov.chars());
				printCovMutex.lock();
				guides[guidetrf[g].g]->printTranscriptGff(c_out);
				printCovMutex.unlock();
			}
		}
	}


	// **** NEXT PORTION NEEDS TO BE CHECKED!!!

	guidetrf.Sort(guideCmp); // guide with most bits set in pattern comes first
	uint refstart=(uint)bdata->start;
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
				guidetrf[g].trf->guide=0;  // this marks a guide that is included in another one -> I might want to swap them if they have the same abundance
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

			/*if(includesource && !longreads) {
				guidetrf[g].trf->nodes.Insert(0,0); // I need to comment this if I need path not to include the source
				guidetrf[g].trf->pattern[0]=1;
				guidetrf[g].trf->pattern[*pos]=1;
			}*/

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

					// cummulative bpcov
					leftcov+=get_cov(1,inode->start-refstart,guides[guidetrf[g].g]->start-1-refstart,bdata->bpcov);
					leftcov/=(guides[guidetrf[g].g]->start-inode->start);
					uint thisend=inode->end+1;
					if(guides[guidetrf[g].g]->end<inode->end) thisend=guides[guidetrf[g].g]->end+1;

					// cummulative bpcov
					if(thisend>guides[guidetrf[g].g]->start) {
						rightcov+=get_cov(1,guides[guidetrf[g].g]->start-refstart,thisend-1-refstart,bdata->bpcov);
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
				//fprintf(stderr,"introduce node from source to %d wih abundance=%g where pos=%d\n",nodei,maxabund,*pos);
				if(longreads) tr->longread=true;
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
			int *pos=gpos[edge(nodei,sink,gno)];
			if(!pos) {
				//fprintf(stderr,"Add sink link from position %d for guide=%d\n",no2gnode[nodei]->end,guidetrf[g].g);
				int key=edge(nodei,sink,gno);
				gpos.Add(key,lastgpos);
				lastgpos++;
				pos=gpos[key];
			}
			/*if(!longreads) {
				guidetrf[g].trf->nodes.Add(sink);
				guidetrf[g].trf->pattern[sink]=1;
				guidetrf[g].trf->pattern[*pos]=1;
			}*/

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

					// cummulative bpcov
					leftcov+=get_cov(1,thisstart-refstart,guides[guidetrf[g].g]->end-refstart,bdata->bpcov);
					leftcov/=(guides[guidetrf[g].g]->end-thisstart+1);

					// cummulative bpcov
					rightcov+=get_cov(1,guides[guidetrf[g].g]->end-refstart,inode->end-refstart,bdata->bpcov);
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
				if(longreads) tr->longread=true;
				transfrag.Add(tr);

				// add sink among maxnode children
				inode->child.Add(sink);
				inode->childpat[*pos]=1; // source should be already among parents of maxnode but not the edge to source
			}

		}

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

int guides_pushmaxflow(int gno,int edgeno,GIntHash<int>& gpos,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,GVec<CGuide>& guidetrf,int& geneno,
		int s,GList<CPrediction>& pred,GVec<float>& nodecov,GBitVec& istranscript,GBitVec& pathpat,bool &first,GPVec<GffObj>& guides,GVec<int> &guidepred, BundleData *bdata) {

	int maxi=1;
	int ng=guidetrf.Count();

	if(ng==1) { // if only one guide I do not need to do the 2 pass
		GVec<float> nodeflux;
		//float fragno=0;
		bool full=true;
		float flux= push_max_flow(gno,guidetrf[0].trf->nodes,istranscript,transfrag,no2gnode,nodeflux,guidetrf[0].trf->pattern,gpos,full);
		istranscript.reset();

		/*
		{ // DEBUG ONLY
			fprintf(stderr,"guide=%s flux[0]=%g\n",guides[guidetrf[0].g]->getID(),flux);
		}
		*/

		if(flux>epsilon) {
			bool include=true;
			if(guidepred[guidetrf[0].g]==-1) {

				store_transcript(pred,guidetrf[0].trf->nodes,nodeflux,nodecov,no2gnode,geneno,first,s,gno,gpos,include,pathpat,false,bdata,guides[guidetrf[0].g]);
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
					store_transcript(pred,guidetrf[g].trf->nodes,nodeflux,nodecov,no2gnode,geneno,first,s,gno,gpos,include,pathpat,false,bdata,guides[guidetrf[g].g]);
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
							store_transcript(pred,guidetrf[g].trf->nodes,nodeflux,nodecov,no2gnode,geneno,first,s,gno,gpos,include,pathpat,false,bdata,guides[guidetrf[g].g]);
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
					initgcov[i]=pred[np]->cov*abs(pred[np]->tlen);
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


int find_transcripts(int gno,int edgeno, GIntHash<int> &gpos,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,int geneno,int strand,
		GVec<CGuide>& guidetrf,GPVec<GffObj>& guides,GVec<int>& guidepred,BundleData* bdata,GVec<int>& trflong) {

	GList<CPrediction>& pred = bdata->pred;

	/*
	if(trflong.Count()) get_trf_long(gno,edgeno, gpos,no2gnode,transfrag,geneno,strand,pred,trflong,bdata);
	if(longreads) return(geneno);*/

	if(longreads) {
		if(trflong.Count()) get_trf_long(gno,edgeno, gpos,no2gnode,transfrag,geneno,strand,pred,trflong,bdata);
		return(geneno);
	}


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

	if (eonly)
		guides_pushmaxflow(gno,edgeno,gpos,no2gnode,transfrag,guidetrf,geneno,strand,pred,nodecov,istranscript,pathpat,first,guides,guidepred,bdata);
	else { //if(!eonly) {
		if(mixedMode && trflong.Count()) get_trf_long_mix(gno,edgeno, gpos,no2gnode,transfrag,geneno,strand,pred,trflong,nodecov,istranscript,pathpat,bdata,first);

		if(!mixedMode && guidetrf.Count()) maxi=guides_pushmaxflow(gno,edgeno,gpos,no2gnode,transfrag,guidetrf,geneno,strand,pred,nodecov,istranscript,pathpat,first,guides,guidepred,bdata);
		/*
		{ // DEBUG ONLY
			if(mixedMode) {
				fprintf(stderr,"After get_trf_long:\n");
				for(int i=0;i<gno;i++) {
					CGraphnode *inode=no2gnode[i];
					printTime(stderr);
					fprintf(stderr,"Node %d: cov=%f capacity=%f rate=%f ",i,inode->cov/(inode->end-inode->start+1),inode->capacity,inode->rate);
					fprintf(stderr,"trf=");
					for(int t=0;t<inode->trf.Count();t++) fprintf(stderr," %d(%f)",inode->trf[t],transfrag[inode->trf[t]]->abundance);
					fprintf(stderr," maxi=%d maxcov=%f\n",maxi,nodecov[maxi]);
				}
				fprintf(stderr,"There are %d transfrags:\n",transfrag.Count());
				for(int t=0;t<transfrag.Count();t++) {
					fprintf(stderr,"%d: ",t);
					//printBitVec(transfrag[s][b][t]->pattern);
					fprintf(stderr," %f(%f,%d,%d) long=%d nodes=%d",transfrag[t]->abundance,transfrag[t]->srabund, transfrag[t]->longstart,transfrag[t]->longend,transfrag[t]->longread,transfrag[t]->nodes.Count());
					for(int i=0;i<transfrag[t]->nodes.Count();i++) fprintf(stderr," %d",transfrag[t]->nodes[i]);
					if(!transfrag[t]->abundance) fprintf(stderr," *");
					fprintf(stderr,"\n");
				}
			}
		}
		*/


		if(nodecov[maxi]>=1) { // sensitive mode only; otherwise >=readthr

			// 1:
			// parse_trf_weight_max_flow(gno,no2gnode,transfrag,geneno,strand,pred,nodecov,pathpat);
			// 2:
			GBitVec usednode(gno+edgeno);
			parse_trf(maxi,gno,edgeno,gpos,no2gnode,transfrag,geneno,first,strand,pred,nodecov,istranscript,usednode,0,pathpat);

		}

	}
	return(geneno);
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

		//if(covered && maxguidelen>=(uint)mintranscriptlen) { if(c_out) guide->printTranscriptGff(c_out);}
		//else covered=false;
		if(maxguidelen<(uint)mintranscriptlen) covered=false;

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

	//fprintf(stderr,"consider junction:%d-%d:%d support=%f,%f nreads=%f\n",jd.start,jd.end,jd.strand,jd.leftsupport,jd.rightsupport,jd.nreads);

	/*** 1. if only we do estimation and it does not match guide do not consider ***/
	if(eonly && !jd.guide_match) { // this way I am using only reads that are compatible to annotated transcripts
		jd.strand=0;
		return false;
	}

	/*** 2. keep if it's in the annotation ***/
	if(jd.guide_match) return true; // this junction is covered by at least one read: the one that calls good_junc
													  // ^ KEEP THIS IN MIND IF WE CHANGE HOW WE USE GOOD_JUNC

	/*** 3. don't keep if it's below threshold ***/
	if (jd.nreads_good<junctionthr) { // junctionthr is always positive
	//if (jd.leftsupport<junctionthr || jd.rightsupport<junctionthr || (jd.end-jd.start>longintron && jd.nreads_good<junctionthr)) {
		jd.strand=0;
		//fprintf(stderr,"nosupport: left=%f right=%f good=%f\n",jd.leftsupport,jd.rightsupport,jd.nreads_good);
		return false;
	}

	// if(!jd.strand) return false; the strand has to be non-zero when we call good_junc -> KEEP THIS IN MIND IF WE CHANGE HOW WE USE GOOD_JUNC

	// don't trust spliced reads that have a very low coverage:
	//int sno=(int)jd.strand+1;
	bool mismatch=false;

	/*** 5. don't keep if it's a long intron and all junctions are from bad reads ***/
	if (jd.nm && round(jd.nm)==round(jd.nreads)) {
		mismatch=true;
		if(jd.end-jd.start>longintron && jd.nreads<CHI_WIN*ERROR_PERC) { // don't believe long intron if all junctions are from bad reads unless highly covered
			jd.strand=0;
			//jd.nm=0; // not sure why this was set here, as it doesn't seem to be used again
			//fprintf(stderr,"bad longintron\n");
			return false;
		}
	}


	/*** 4. don't keep if coverage remains relatively constant and the junction is not a high fraction of coverage ***/
	// cummulative bpcov on bw=3 (I should be looking in a 'codon' window
	float mult=1/ERROR_PERC; // test if this is needed for longreads
	//if(!longreads) {
	if(longreads) mult/=ERROR_PERC; // higher tolerance for long reads
	int bw=5;
	int j=jd.start-refstart;

	float lleftcov=0;
	float lrightcov=0;
	int sno=(int)jd.strand+1;
	if(j+bw+1<bpcov[1].Count() && j>=bw) {
		lleftcov=get_cov(1,j-bw+1,j,bpcov);
		lrightcov=get_cov(1,j+1,j+bw,bpcov);
		if(longreads) {
			lleftcov-=get_cov(2-sno,j-bw+1,j,bpcov);
			lrightcov-=get_cov(2-sno,j+1,j+bw,bpcov);
		}
	}

	// shouldn't leftcov be normalized by bw below when comparing to leftsupport?
	if(lleftcov>1/ERROR_PERC && jd.leftsupport*mult<ERROR_PERC*lleftcov && (mismatch || lrightcov>lleftcov*(1-ERROR_PERC))) { // gave some boost (*10) to junctions here assumming not all of them are captured
		jd.strand=0;
		//fprintf(stderr," left j=%d leftcov=%f leftsupprt=%f rightcov=%f\n",j+refstart,lleftcov,jd.leftsupport,lrightcov);
		return false;
	}

	j=jd.end-refstart-1;

	/* strand based */
	float rleftcov=0;
	float rrightcov=0;
	if(j+bw+1<bpcov[1].Count() && j>=bw) {
		rleftcov=get_cov(1,j-bw+1,j,bpcov);
		rrightcov=get_cov(1,j+1,j+bw,bpcov);
		if(longreads) {
			rleftcov-=get_cov(2-sno,j-bw+1,j,bpcov);
			rrightcov-=get_cov(2-sno,j+1,j+bw,bpcov);
		}
	}

	if(rrightcov>1/ERROR_PERC && jd.rightsupport*mult<rrightcov*ERROR_PERC && (mismatch || rleftcov>rrightcov*(1-ERROR_PERC))) { // gave some boost (*10) to junctions here assumming not all of them are captured
		jd.strand=0;
		//fprintf(stderr," right j=%d leftcov=%f rightsupprt=%f rightcov=%f\n",j+refstart,rleftcov,jd.rightsupport,rrightcov);
		return false;
	}

	/*if(lleftcov<rrightcov*ERROR_PERC*DROP || lleftcov*ERROR_PERC*DROP>rrightcov) { // significant drop between exons
			jd.strand=0;
			return false;
	}*/
	//}

	/*** 6. don't keep it if it's not above 1% of all junctions -> might be spliceosome error ***/
	//if(jd.nreads<1/ERROR_PERC && (jd.nreads*10<ERROR_PERC*jd.leftsupport || jd.nreads*10<ERROR_PERC*jd.rightsupport)) {
	if(jd.nreads*10<ERROR_PERC*jd.leftsupport || jd.nreads*10<ERROR_PERC*jd.rightsupport) {
		//fprintf(stderr,"nreads=%f leftsuport=%f rightsupport=%f\n",jd.nreads,jd.leftsupport,jd.rightsupport);
		jd.strand=0;
		return false;
	}


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
	}
	j=nj+1;
	while(j<junction.Count() && junction[j]->start==jd.start) {
		sumjunc+=junction[j]->nreads_good; j++;
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
	/*GPVec<GPtFeature>& feature = bdata->ptfs; // these are point features (confirmed starts/stops)

	for(int i=0;i<feature.Count();i++) {
		if(feature[i]->ftype==GPFT_TSS)
			fprintf(stderr,"TSS at position %d on strand %d\n",feature[i]->coord,feature[i]->strand);
		if(feature[i]->ftype==GPFT_CPAS)
			fprintf(stderr,"CPAS at position %d on strand %d\n",feature[i]->coord,feature[i]->strand);
	}*/

	//fprintf(stderr,"build_graphs with %d guides\n",guides.Count());


	if(guides.Count()) {

		guideedge.setSorted(true);
		guideedge.setUnique(true);

		//if(eonly)
		for(int g=0;g<guides.Count();g++) {
			guidepred.cAdd(-1);
			bool covered=true;
			RC_TData* tdata=(RC_TData*)(guides[g]->uptr);
			if(longreads || mixedMode) {
				for(int i=1;i<guides[g]->exons.Count();i++) {
					char s=0; // unknown strand
					if(guides[g]->strand=='+') s=1; // guide on positive strand
					else if(guides[g]->strand=='-') s=-1; // guide on negative strand
					CJunction jn(guides[g]->exons[i-1]->end,guides[g]->exons[i]->start,s);
					int oidx=-1;
					if (!junction.Found(&jn, oidx)) {
						covered=false;
						break;
					}
				}
			}
			else {
				for(int i=0;i<tdata->t_introns.Count();i++) {
					if(!tdata->t_introns[i]->rcount) {
						covered=false;
						break;
					}
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
	double leftsupport[2]={0,0};    // should be strand based
	double rightsupport[2]={0,0};
	bool higherr=false;
	char leftcons=-1;
	char rightcons=-1;
	for(int i=0;i<junction.Count();i++) {

		//fprintf(stderr,"check junction:%d-%d:%d leftsupport=%f rightsupport=%f nm=%f nreads=%f\n",junction[i]->start,junction[i]->end,junction[i]->strand,junction[i]->leftsupport,junction[i]->rightsupport,junction[i]->nm,junction[i]->nreads);

		if((!higherr || mixedMode) && junction[i]->strand && junction[i]->nm==junction[i]->nreads && !junction[i]->guide_match) {
			higherr=true;
			if(mixedMode) {
				int j=i-1;
				while(j>=0 && junction[j]->start+sserror>junction[i]->start) {
					if(junction[j]->strand && junction[i]->strand!=junction[j]->strand && abs((int)(junction[j]->end-junction[i]->end))<(int)sserror && junction[j]->nm<junction[j]->nreads) {
						junction[i]->strand = 0;
						break;
					}
					j--;
				}
				if(junction[i]->strand) {
					j=i+1;
					while(j<junction.Count() && junction[j]->start-sserror<junction[i]->start) {
						if(junction[j]->strand && junction[i]->strand!=junction[j]->strand && abs((int)(junction[j]->end-junction[i]->end))<(int)sserror && junction[j]->nm<junction[j]->nreads) {
							junction[i]->strand = 0;
							break;
						}
						j++;
					}
				}
			}
		}

		if(junction[i]->start!=start) {  // new junction starting at i

			int j=i-1;
			while(j>=0 && junction[j]->start==start) {
				junction[j]->consleft=leftcons;
				if(!junction[j]->guide_match && !leftcons && junction[j]->nreads_good<DROP/ERROR_PERC) junction[j]->strand=0;
				if(junction[j]->strand) junction[j]->leftsupport=leftsupport[(1+junction[j]->strand)/2];
				j--;
			}


			j=i+1; // check if there is the same junction further ahead first
			if(junction[i]->strand) while(j<junction.Count() && junction[j]->start==junction[i]->start && junction[j]->end==junction[i]->end) {
				if(junction[j]->strand && junction[i]->strand!=junction[j]->strand) {
					//possible missaligned junction --> delete junction strand
					if(junction[i]->nreads>junction[j]->nreads && junction[i]->nreads_good>junction[j]->nreads_good) junction[j]->strand=0;
					if(junction[i]->nreads<junction[j]->nreads && junction[i]->nreads_good<junction[j]->nreads_good) junction[i]->strand=0;
				}
				j++;
			}

			leftsupport[0]=0;
			leftsupport[1]=0;
			if(junction[i]->strand) leftsupport[(1+junction[i]->strand)/2]=junction[i]->leftsupport; // I might have deleted the i junction above
			start=junction[i]->start;

			if(bdata->gseq) {
				if(junction[i]->strand>0) {
					if((bdata->gseq[junction[i]->start+1-refstart]=='g'||bdata->gseq[junction[i]->start+1-refstart]=='G') &&
							(bdata->gseq[junction[i]->start+2-refstart]=='T'||bdata->gseq[junction[i]->start+2-refstart]=='t'||
									bdata->gseq[junction[i]->start+2-refstart]=='C'||bdata->gseq[junction[i]->start+2-refstart]=='c')) {
						junction[i]->consleft=1;
						leftcons=1;
					}
					else {
						junction[i]->consleft=0;
						leftcons=0;
					}
				}
				else if(junction[i]->strand<0) {
					if((bdata->gseq[junction[i]->start+1-refstart]=='c'||bdata->gseq[junction[i]->start+1-refstart]=='C') &&
							(bdata->gseq[junction[i]->start+2-refstart]=='T'||bdata->gseq[junction[i]->start+2-refstart]=='t')) {
						junction[i]->consleft=1;
						leftcons=1;
					}
					else {
						junction[i]->consleft=0;
						leftcons=0;
					}
				}
				//fprintf(stderr,"junction:%d-%d:%d is %c%c-%c%c\n",junction[i]->start,junction[i]->end,junction[i]->strand,bdata->gseq[junction[i]->start+1-refstart],bdata->gseq[junction[i]->start+2-refstart],bdata->gseq[junction[i]->end-2-refstart],bdata->gseq[junction[i]->end-1-refstart]);
			}
		}
		else if(junction[i]->strand) leftsupport[(1+junction[i]->strand)/2]+=junction[i]->leftsupport;
		//fprintf(stderr,"leftsupport[%d]=%f\n",(1+junction[i]->strand)/2,leftsupport[(1+junction[i]->strand)/2]);


		//fprintf(stderr,"check ejunction:%d-%d:%d leftsupport=%f rightsupport=%f nm=%f nreads=%f\n",ejunction[i]->start,ejunction[i]->end,ejunction[i]->strand,ejunction[i]->leftsupport,ejunction[i]->rightsupport,ejunction[i]->nm,ejunction[i]->nreads);
		if(ejunction[i]->end!=end) { // I do not check if I deleted the junction here for support

			int j=i-1;
			while(j>=0 && ejunction[j]->end==end) {
				ejunction[j]->consright=rightcons;
				if(!ejunction[j]->guide_match && !rightcons && ejunction[j]->nreads_good<DROP/ERROR_PERC) ejunction[j]->strand=0;
				if(ejunction[j]->strand) ejunction[j]->rightsupport=rightsupport[(1+ejunction[j]->strand)/2];
				j--;
			}
			// I do not check here for possible missalignments
			rightsupport[0]=0;
			rightsupport[1]=0;
			if(ejunction[i]->strand) rightsupport[(1+ejunction[i]->strand)/2]=ejunction[i]->rightsupport;
			end=ejunction[i]->end;

			if(bdata->gseq) {
				if(ejunction[i]->strand>0) {
					if((bdata->gseq[ejunction[i]->end-2-refstart]=='a'||bdata->gseq[ejunction[i]->end-2-refstart]=='A') &&
							(bdata->gseq[ejunction[i]->end-1-refstart]=='G'||bdata->gseq[ejunction[i]->end-1-refstart]=='g')) {
						ejunction[i]->consright=1;
						rightcons=1;
					}
					else {
						ejunction[i]->consright=0;
						rightcons=0;
					}
				}
				else if(ejunction[i]->strand<0) {
					if((bdata->gseq[ejunction[i]->end-1-refstart]=='C'||bdata->gseq[ejunction[i]->end-1-refstart]=='c') &&
							(bdata->gseq[ejunction[i]->end-2-refstart]=='A'||bdata->gseq[ejunction[i]->end-2-refstart]=='a'||
									bdata->gseq[ejunction[i]->end-2-refstart]=='G'||bdata->gseq[ejunction[i]->end-2-refstart]=='g')) {
						ejunction[i]->consright=1;
						rightcons=1;
					}
					else {
						ejunction[i]->consright=0;
						rightcons=0;
					}
				}
				//fprintf(stderr,"ejunction:%d-%d:%d is %c%c-%c%c\n",ejunction[i]->start,ejunction[i]->end,ejunction[i]->strand,bdata->gseq[ejunction[i]->start+1-refstart],bdata->gseq[ejunction[i]->start+2-refstart],bdata->gseq[ejunction[i]->end-2-refstart],bdata->gseq[ejunction[i]->end-1-refstart]);
			}

		}
		else if(ejunction[i]->strand) rightsupport[(1+ejunction[i]->strand)/2]+=ejunction[i]->rightsupport;
		//fprintf(stderr,"rightsupport[%d]=%f\n",(1+ejunction[i]->strand)/2,rightsupport[(1+ejunction[i]->strand)/2]);
	}
	// end adjusting leftsupport and rightsupport

	//fprintf(stderr,"junction support computed\n");

	if(higherr) { // there are some reads that contain very bad junctions -> need to find better closest junctions

		uint juncsupport=junctionsupport;
		if(longreads)
			juncsupport=sserror;
		else if(mixedMode) juncsupport=sserror/DROP;
		//fprintf(stderr,"In higherr!\n");
		GVec<int> jstarts; // good starting junctions
		GVec<int> jends; // good ending junctions
		if(viral) {
			for(int i=1;i<junction.Count();i++) { // junction is sorted based on start

				if(junction[i]->strand && junction[i]->nm && !junction[i]->guide_match && junction[i]->nm>=junction[i]->nreads) { // this is a bad junction -> check if it's maximal;
					if(junction[i]->nreads_good>=0 && (junction[i]->nreads_good<1.25*junctionthr || !good_junc(*junction[i],refstart,bpcov))) { // threshold for bad junctions is higher; (should I also add that too short junctions not to be accepted?)
						//junction[i]->strand=0; // just delete junction if it's low count
						junction[i]->mm=-1;
						//fprintf(stderr,"...delete due to being under threshold\n");
					}

					int j=i-1;
					while(j>0 && junction[i]->start-junction[j]->start<juncsupport) {
						if(junction[j]->strand==junction[i]->strand && junction[j]->nreads_good>=0 && junction[i]->nreads<junction[j]->nreads &&
								abs((int)junction[i]->end-(int)junction[j]->end)<(int)juncsupport) { // j was not elminated and is better
							if(junction[i]->nreads_good<0) { // i was eliminated before
								int k=-junction[i]->nreads_good;
								if(junction[k]->nreads<junction[j]->nreads) junction[i]->nreads_good=-j;
							}
							else {
								junction[i]->nreads_good=-j;
							}
						}
						j--;
					}
				}
			}
			for(int i=junction.Count()-1;i>0;i--) {
				if(junction[i]->strand && junction[i]->nm && !junction[i]->guide_match && junction[i]->nm>=junction[i]->nreads) { // this is a bad junction -> check if it's maximal;
					int j=i+1;
					while(j<junction.Count() && junction[j]->start-junction[i]->start<juncsupport) {
						if(junction[j]->strand==junction[i]->strand && junction[j]->nreads_good>=0 && junction[i]->nreads<junction[j]->nreads &&
								abs((int)junction[i]->end-(int)junction[j]->end)<(int)juncsupport) { // j was not elminated and is better
							if(junction[i]->nreads_good<0) { // i was eliminated before
								int k=-junction[i]->nreads_good;
								if(junction[k]->nreads<junction[j]->nreads) junction[i]->nreads_good=-j;
							}
							else {
								junction[i]->nreads_good=-j;
							}
						}
						j++;
					}
				}
			}
		}
		else {

			float tolerance=1-ERROR_PERC;

			// strand based version
			for(int i=1;i<junction.Count();i++) { // junction is sorted based on start

				//fprintf(stderr,"junct[%d]:%d-%d:%d lefttsupport=%f nm=%f mm=%f nreads=%f nreads_good=%f\n",i,junction[i]->start,junction[i]->end,junction[i]->strand,junction[i]->leftsupport,junction[i]->nm,junction[i]->mm,junction[i]->nreads,junction[i]->nreads_good);
				if(junction[i]->strand) {
					if(junction[i]->nm && !junction[i]->guide_match && junction[i]->nm>=junction[i]->nreads) { // this is a bad junction -> check if it's maximal;
					if(junction[i]->nreads_good>=0 && junction[i]->nreads_good<1.25*junctionthr) { // threshold for bad junctions is higher; (should I also add that too short junctions not to be accepted?)
						//junction[i]->strand=0; // just delete junction if it's low count
						junction[i]->mm=-1;
						//fprintf(stderr,"...delete due to being under threshold\n");
					}

					int j=i-1;
					float support=0;
					bool searchjunc=true;
					bool reliable=false;
					//if(j>=0) fprintf(stderr,"...start at junct:%d-%d:%d leftsupport=%f dist=%d\n",junction[j]->start,junction[j]->end,junction[j]->strand,junction[j]->leftsupport,junction[i]->start-junction[j]->start);
					while(j>0 && junction[i]->start-junction[j]->start<juncsupport) {
						if(junction[j]->strand==junction[i]->strand) {
							if(junction[j]->start==junction[i]->start) { // found a junction with the same start -> I have already searched it if it's bad
								if(junction[j]->nreads<0) {
									junction[i]->nreads=junction[j]->nreads;
									searchjunc=false;
								}
								break;
							}
							else if(junction[j]->guide_match || junction[j]->nm<junction[j]->nreads) { // nearby junction is much more reliable
								if(junction[j]->leftsupport>junction[i]->leftsupport*tolerance) { // the good junction is close enough
									reliable=true;
									junction[i]->nreads=-j;
									support=junction[j]->leftsupport;
									break;
								}
							}
							else if(junction[j]->leftsupport>support && junction[i]->start-junction[j]->start<sserror && junction[j]->leftsupport*tolerance>junction[i]->leftsupport) {
								//fprintf(stderr,"...1 compare to [%d]:%d-%d:%d leftsupport=%f\n",j,junction[j]->start,junction[j]->end,junction[j]->strand,junction[j]->leftsupport);
								junction[i]->nreads=-j;
								support=junction[j]->leftsupport;
							}
						}
						j--;
					}
					if(searchjunc) {
						j=i+1;
						int dist=juncsupport;
						if(junction[i]->nreads<0) {
							dist=junction[i]->start-junction[abs((int)junction[i]->nreads)]->start;
						}
						//if(j<junction.Count()) fprintf(stderr,"...start at junct:%d-%d:%d leftsupport=%f dist=%d\n",junction[j]->start,junction[j]->end,junction[j]->strand,junction[j]->leftsupport,junction[j]->start-junction[i]->start);
						while(j<junction.Count() && junction[j]->start-junction[i]->start<juncsupport) {
							if(junction[j]->strand==junction[i]->strand && junction[j]->start!=junction[i]->start) {
								int d=(int)(junction[j]->start-junction[i]->start);
								if(junction[j]->guide_match || junction[j]->nm<junction[j]->nreads) {
									if((d<dist || (d==dist && junction[j]->leftsupport>support)) && junction[j]->leftsupport>junction[i]->leftsupport*tolerance) {
										junction[i]->nreads=-j;
										support=junction[j]->leftsupport;
										break;
									}
								}
								else if(!reliable && junction[j]->leftsupport>support && (uint)d<sserror && junction[j]->leftsupport*tolerance>junction[i]->leftsupport) { // junction is not best within window
									//fprintf(stderr,"...2 compare to [%d]:%d-%d:%d leftsupport=%f\n",j,junction[j]->start,junction[j]->end,junction[j]->strand,junction[j]->leftsupport);
									junction[i]->nreads=-j;
									support=junction[j]->leftsupport;
								}
							}
							j++;
						}
					}
				}
					else if(mixedMode && junction[i]->nm<junction[i]->nreads) {
						jstarts.Add(i);
					}
				}
				//fprintf(stderr,"ejunct[%d]:%d-%d:%d rightsupport=%f nm=%f nreads=%f\n",i,ejunction[i]->start,ejunction[i]->end,ejunction[i]->strand,ejunction[i]->rightsupport,ejunction[i]->nm,ejunction[i]->nreads);

				if(ejunction[i]->strand) {
					if(ejunction[i]->nm && !ejunction[i]->guide_match && ejunction[i]->nm>=ejunction[i]->nreads) { // this is a bad junction -> check if it's maximal
					if(ejunction[i]->nreads_good>=0 && ejunction[i]->nreads_good<1.25*junctionthr) { // threshold for bad junctions is higher
						//ejunction[i]->strand=0;
						ejunction[i]->mm=-1;
						//fprintf(stderr,"...delete due to being under threshold\n");
					}
					int j=i-1;
					float support=0;
					bool searchjunc=true;
					bool reliable=false;
					//if(j>=0) fprintf(stderr,"...start at junct:%d-%d:%d rightsupport=%f dist=%d\n",ejunction[j]->start,ejunction[j]->end,ejunction[j]->strand,ejunction[j]->rightsupport,ejunction[i]->end-ejunction[j]->end);
					while(j>0 && ejunction[i]->end-ejunction[j]->end<juncsupport) {
						if(ejunction[j]->strand==ejunction[i]->strand) {
							if(ejunction[j]->end==ejunction[i]->end) {
								if(ejunction[j]->nreads_good<0) {
									ejunction[i]->nreads_good=ejunction[j]->nreads_good;
									searchjunc=false;
								}
								break;
							}
							else if(ejunction[j]->guide_match || ejunction[j]->nm<ejunction[j]->nreads) { // nearby junction is much more reliable
								if(ejunction[j]->rightsupport>ejunction[i]->rightsupport*tolerance) { // the good junction is close enough
									reliable=true;
									ejunction[i]->nreads_good=-j;
									support=ejunction[j]->rightsupport;
									break;
								}
							}
							else if(ejunction[j]->rightsupport>support && ejunction[i]->end-ejunction[j]->end < sserror && ejunction[j]->rightsupport*tolerance>ejunction[i]->rightsupport) {
								//fprintf(stderr,"...1 compare to [%d]:%d-%d:%d rightsupport=%f\n",j,ejunction[j]->start,ejunction[j]->end,ejunction[j]->strand,ejunction[j]->rightsupport);
								ejunction[i]->nreads_good=-j;
								support=ejunction[j]->rightsupport;
							}
						}
						j--;
					}
					if(searchjunc) {
						j=i+1;
						int dist=juncsupport;
						if(ejunction[i]->nreads_good<0) {
							dist=ejunction[i]->end-ejunction[abs((int)ejunction[i]->nreads_good)]->end;
						}
						//if(j<junction.Count()) fprintf(stderr,"...start at junct:%d-%d:%d rightsupport=%f dist=%d\n",ejunction[j]->start,ejunction[j]->end,ejunction[j]->strand,ejunction[j]->rightsupport,ejunction[j]->end-ejunction[i]->end);
						while(j<junction.Count() && ejunction[j]->end-ejunction[i]->end<juncsupport) {
							if(ejunction[j]->strand==ejunction[i]->strand && ejunction[j]->end!=ejunction[i]->end) {
								int d=ejunction[j]->end-ejunction[i]->end;
								if(ejunction[j]->guide_match || ejunction[j]->nm<ejunction[j]->nreads) {
									if((d<dist || (d==dist && ejunction[j]->rightsupport>support)) && ejunction[j]->rightsupport>ejunction[i]->rightsupport*tolerance) {
										ejunction[i]->nreads_good=-j;
										support=ejunction[j]->rightsupport;
										break;
									}
								}
								else if((!reliable && ejunction[j]->rightsupport>support && ejunction[j]->end-ejunction[i]->end < sserror && ejunction[j]->rightsupport*tolerance>ejunction[i]->rightsupport) || ((int)(ejunction[j]->end-ejunction[i]->end)<dist && (ejunction[j]->guide_match || ejunction[j]->nm<ejunction[j]->nreads))) {
									//fprintf(stderr,"...2 compare to [%d]:%d-%d:%d rightsupport=%f\n",j,ejunction[j]->start,ejunction[j]->end,ejunction[j]->strand,ejunction[j]->rightsupport);
									ejunction[i]->nreads_good=-j;
									support=ejunction[j]->rightsupport;
								}
							}
							j++;
						}
					}
				}
					else if(mixedMode && ejunction[i]->nm<ejunction[i]->nreads) {
						jends.Add(i);
					}
				}
			}
		}
		if(mixedMode) { // check if there are junctions inside bigger junctions that can form small exons
			int si=0;
			for(int ei=0;ei<jends.Count();ei++) {
				while(si<jstarts.Count() && junction[jstarts[si]]->start<=ejunction[jends[ei]]->end) si++;
				int k=si;
				while(k<jstarts.Count() && junction[jstarts[k]]->start-ejunction[jends[ei]]->end<SMALL_EXON) {
					if(junction[jstarts[k]]->strand == ejunction[jends[ei]]->strand) {
						CJunction jn(ejunction[jends[ei]]->start,junction[jstarts[k]]->end,junction[jstarts[k]]->strand);
						int oidx=-1;
						if (junction.Found(&jn, oidx) && junction[oidx]->nm>=junction[oidx]->nreads) { // candidate junction for deletion
							junction[oidx]->strand=0;
							break;
						}
					}
					k++;
				}
			}
		}
	} //if higherr

	/*
	{ // DEBUG ONLY
		for(int i=0;i<junction.Count();i++) {
			if(junction[i]->guide_match) fprintf(stderr,"G");
			if(junction[i]->strand && junction[i]->nreads>0 && junction[i]->nreads_good>0 && junction[i]->mm>=0) fprintf(stderr,"***");
			//if(junction[i]->strand && junction[i]->mm>=0)
				fprintf(stderr,"Junction[%d] %d-%d:%d has nm=%f mm=%f nreads=%f nreads_good=%f leftsupport=%f and rightsupport=%f\n",i,junction[i]->start,junction[i]->end,junction[i]->strand,
					junction[i]->nm,junction[i]->mm,junction[i]->nreads,
					junction[i]->nreads_good,junction[i]->leftsupport,junction[i]->rightsupport);
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

	//GHash<bool> boundaryleft;
	//GHash<bool> boundaryright;
	GIntHash<bool> boundaryleft;
	GIntHash<bool> boundaryright;

	bool resort=false;
	int njunc=junction.Count();
	for (int n=0;n<readlist.Count();n++) {
		CReadAln & rd=*(readlist[n]);

		// version that does not adjust read but discards it completely
		bool keep=true;
		int i=0;

		/*fprintf(stderr,"check read[%d]:%d-%d:%d refstart=%d w/exons:",n,readlist[n]->start,readlist[n]->end,readlist[n]->strand,refstart);
		for(i=0;i<rd.juncs.Count();i++) { fprintf(stderr," %d-%d:%d",rd.segs[i].start,rd.segs[i].end,rd.juncs[i]->strand);}
		fprintf(stderr," %d-%d\n",rd.segs[i].start,rd.segs[i].end);
		i=0;*/


		while(i<rd.juncs.Count()) {
			CJunction& jd=*(rd.juncs[i]);
			//fprintf(stderr, " read junc %d-%d:%d nreads=%f nreads_good=%f nm=%f support=%f,%f between exons:%d-%d and %d-%d\n", jd.start, jd.end, jd.strand,jd.nreads,jd.nreads_good,jd.nm,jd.leftsupport,jd.rightsupport,rd.segs[i].start,rd.segs[i].end,rd.segs[i+1].start,rd.segs[i+1].end);

			bool changeright=jd.nreads_good<0;
			bool changeleft=jd.nreads<0;
			if(viral && changeright) {
				changeleft=true;
				jd.nreads=jd.nreads_good;
			}

			if(jd.strand) {
				if(!good_junc(jd,refstart,bpcov)) jd.mm=-1; // bad junction
				if(changeleft || changeright) jd.strand=0;
				if(jd.mm<0) jd.strand=0; // I can't believe the junction because it's too low
			}

			if(!jd.strand) { // found a bad junction
				//instead of deleting read -> just delete junction and make sure that add_read_to_group adds read to multiple groups
				//jd.nreads-=rd.read_count;

				//if(jd.mm>=0 && (changeleft || changeright)) { // not sure why jd.mm should be positive here since I am searching the junction against the others anyway
				if(changeleft || changeright) {
					uint newstart=rd.segs[i].end;
					uint newend=rd.segs[i+1].start;
					int jk=-1;
					int ek=-1;
					if(changeleft) {
						jk=abs(int(jd.nreads));
						if(junction[jk]->nreads<0) { // not a valid left support -> delete junction
							newstart=rd.segs[i].start-1;
						}
						else {
							newstart=junction[jk]->start; // junction jk is good
							//fprintf(stderr,"junction has newstart=%d from junction[jk=%d]\n",newstart,jk);
							if(junction[jk]->strand) jk=-1;
						}
					}
					if(changeright) {
						ek=abs(int(jd.nreads_good));
						if(viral) {
							if(junction[ek]->nreads_good<0) {
								newend=rd.segs[i+1].end+1;
							}
							else {
								newend=junction[ek]->end;
								//fprintf(stderr,"junction has newend=%d from junction[ek=%d]\n",newend,ek);
							}
						}
						else {
							if(ejunction[ek]->nreads_good<0) {
								newend=rd.segs[i+1].end+1;
							}
							else {
								newend=ejunction[ek]->end;
								//fprintf(stderr,"junction has newend=%d from junction[ek=%d]\n",newend,ek);
								if(ejunction[ek]->strand) ek=-1;
							}
						}
					}
					//if(jd.mm>=0 && newstart>=rd.segs[i].start && newend<=rd.segs[i+1].end) { // junction inside read boundaries
					if(newstart>=rd.segs[i].start && newend<=rd.segs[i+1].end && newstart<=newend) { // junction inside read boundaries
						bool searchjunc=true;
						bool addjunction=true;
						rd.segs[i].end=newstart;   // adjust start
						rd.segs[i+1].start=newend; // adjust end
						if(jd.mm>=0) {
							if(jk>0) { // search junctions
								int k=jk;
								while(k>0 && junction[k]->start==newstart) {
								   if(mixedMode && junction[k]->nm<junction[k]->nreads) addjunction=false;
								   if(rd.strand==junction[k]->strand && junction[k]->end==newend) {
									rd.juncs.Put(i,junction[k]);
										searchjunc=false;
										break;
									}
									k--;
								}
								if(searchjunc) {
									k=jk+1;
									while(k<junction.Count() && junction[k]->start==newstart) {
										if(mixedMode && junction[k]->nm<junction[k]->nreads) addjunction=false;
										if(rd.strand==junction[k]->strand && junction[k]->end==newend) {
											rd.juncs.Put(i,junction[k]);
											searchjunc=false;
											break;
										}
										k++;
									}
								}
								if(searchjunc) { // I did not find junction -> I need to create a new one
									// first check if I already added such a junction
									for(k=njunc;k<junction.Count();k++) {
										if(rd.strand==junction[k]->strand && junction[k]->start==newstart && junction[k]->end==newend) {
											rd.juncs.Put(i,junction[k]);
											searchjunc=false;
											break;
										}
									}
									if(searchjunc && addjunction) {
										if(!resort) {
											junction.setSorted(false);
											ejunction.setSorted(false);
											resort=true;
										}
										//fprintf(stderr,"Add new junction:%d-%d:%d at position %d njunc=%d\n",newstart,newend,rd.strand,junction.Count(),njunc);
										CJunction *junc=new CJunction(newstart,newend,rd.strand);
										junction.Add(junc);
										ejunction.Add(junc);
										rd.juncs.Put(i,junc);
									}
								}
							}
							else if(ek>0) { // ek>0 => search ejunctions; ek>0 because I have either changeleft or changeright
								int k=ek;
								while(k>0 && ejunction[k]->end==newend) {
									if(mixedMode && junction[k]->nm<junction[k]->nreads) addjunction=false;
									if(rd.strand==ejunction[k]->strand && ejunction[k]->start==newstart) {
										rd.juncs.Put(i,ejunction[k]);
										searchjunc=false;
										break;
									}
									k--;
								}
								if(searchjunc) {
									k=ek+1;
									while(k<ejunction.Count() && ejunction[k]->end==newend) {
										if(mixedMode && junction[k]->nm<junction[k]->nreads) addjunction=false;
										if(rd.strand==ejunction[k]->strand && ejunction[k]->start==newstart) {
											rd.juncs.Put(i,ejunction[k]);
											searchjunc=false;
											break;
										}
										k++;
									}
								}
								if(searchjunc && addjunction) { // I did not find junction -> I need to create a new one
									// first check if I already added such a junction
									for(k=njunc;k<ejunction.Count();k++) {
										if(rd.strand==ejunction[k]->strand && ejunction[k]->start==newstart && ejunction[k]->end==newend) {
											rd.juncs.Put(i,ejunction[k]);
											searchjunc=false;
											break;
										}
									}
									if(searchjunc) {
										if(!resort) {
											junction.setSorted(false);
											ejunction.setSorted(false);
											resort=true;
										}
										//fprintf(stderr,"Add new junction:%d-%d:%d at position %d njunc=%d\n",newstart,newend,rd.strand,junction.Count(),njunc);
										CJunction *junc=new CJunction(newstart,newend,rd.strand);
										junction.Add(junc);
										ejunction.Add(junc);
										rd.juncs.Put(i,junc);
									}
								}
							}
						}
					}
					//fprintf(stderr, "read[%d] adjusted to junction:%d-%d\n",n,rd.segs[i].end,rd.segs[i+1].start);
				}

				// because read might be poorly mapped I also have to unpair it
				if(keep) { // this is the first time I unpair read -> remove strand if pair has single exon?
					for(int p=0;p<rd.pair_idx.Count();p++) {
						int np=rd.pair_idx[p];
						if(np>-1) {
							rd.pair_idx[p]=-1;
							for(int j=0;j<readlist[np]->pair_idx.Count();j++) // also unpair read np to n
								if(readlist[np]->pair_idx[j]==n) {
									readlist[np]->pair_idx[j]=-1;
									break;
								}
							// remove strand for pair read if single exon ? is strand assigned before here ?
							if(!readlist[np]->juncs.Count()) readlist[np]->strand=0;
						}
					}
					keep=false;
				}
			}
			else if(guides.Count()){ // need to remember boundary
				bool exist=true;
		    	//GStr bs((int)jd.start);
		    	if(!boundaryleft[jd.start]) boundaryleft.Add(jd.start,exist);
		    	//GStr be((int)jd.end);
		    	if(!boundaryright[jd.end]) boundaryright.Add(jd.end,exist);
			}
			i++;

		}

		if(!keep) { // read has bad junctions -> check if I need to unstrand it
			if(rd.strand) {
				bool keepstrand=false;
				for(int j=0;j<rd.juncs.Count();j++) if(rd.juncs[j]->strand) { keepstrand=true; break; }
				if(!keepstrand) rd.strand=0;
			}
		}
		else if(!rd.juncs.Count() && rd.strand) { // check if I need to unstrand current read due to future mapping
			if(rd.pair_idx.Count()) {
				bool keepstrand=false;
				for(int p=0;p<rd.pair_idx.Count();p++) {
					int np=rd.pair_idx[p];
					if(np>-1 && n<np) {
						if(!readlist[np]->juncs.Count()) {
							if(readlist[np]->strand==rd.strand) {
								keepstrand=true;
								break;
							}
						}
						else {
							for(int j=0;j<readlist[np]->juncs.Count();j++) {
								CJunction& jd=*(readlist[np]->juncs[j]);
								if(jd.strand && good_junc(jd,refstart,bpcov)) {
									keepstrand=true;
									break;
								}
							}
							if(!keepstrand) readlist[np]->strand=0;
						}
					}
					if(keepstrand) break;
				}
				if(!keepstrand) rd.strand=0;
			}
		}


		//if(rd.juncs.Count()) fprintf(stderr,"read[%d] keep=%d\n",n,keep);
		//if(rd.strand) fprintf(stderr,"read[%d] has strand %d\n",n,rd.strand);


		//if(keep) { // if it's a good read that needs to be kept


			/*fprintf(stderr,"add read %d:%d-%d w/count=%g for color=%d with npairs=%d\n",n,readlist[n]->start,readlist[n]->end,readlist[n]->read_count,color,readlist[n]->pair_idx.Count());
			fprintf(stderr,"add read[%d]:%d-%d:%d w/count=%g w/exons:",n,readlist[n]->start,readlist[n]->end,readlist[n]->strand,readlist[n]->read_count);
			for(i=0;i<rd.juncs.Count();i++) { fprintf(stderr," %d-%d:%d",rd.segs[i].start,rd.segs[i].end,rd.juncs[i]->strand);}
			fprintf(stderr," %d-%d\n",rd.segs[i].start,rd.segs[i].end);*/

			color=add_read_to_group(n,readlist,color,group,currgroup,startgroup,readgroup,equalcolor,merge);

			// count fragments
			if(!rd.unitig)
				bdata->frag_len+=rd.len*rd.read_count; // TODO: adjust this to work with FPKM for super-reads and Pacbio
			double single_count=rd.read_count;
			if(keep) for(int i=0;i<rd.pair_idx.Count();i++) {
				// I am not counting the fragment if I saw the pair before and it wasn't deleted
				if(rd.pair_idx[i]!=-1 && n>rd.pair_idx[i] && readlist[rd.pair_idx[i]]->nh) {// only if read is paired and it comes first in the pair I count the fragments
					single_count-=rd.pair_count[i];
				}
			}
			if(!rd.unitig && single_count>epsilon) {
				bdata->num_fragments+=single_count; // TODO: FPKM will not work for super-reads here because I have multiple fragments in
												    // a super-read -> I might want to re-estimate this from coverage and have some input for read length; or I might only use TPM
			}

			//fprintf(stderr,"now color=%d\n",color);
		//}
		//else { fprintf(stderr,"read[%d] is not kept\n",n);}
		//else clean_read_junctions(readlist[n]);
	}

	if(resort) {
		junction.setSorted(true);
		ejunction.setSorted(juncCmpEnd);
	}


	//fprintf(stderr,"fragno=%d fraglen=%g\n",fragno,fraglen);

	//if(fragno) fraglen/=fragno;

	// merge groups that are close together or __groups that are within the same exon of a reference gene__
	if(bundledist || (guides.Count() && !longreads)) {
		for(int sno=0;sno<3;sno++) {
			CGroup *lastgroup=NULL;
			CGroup *procgroup=startgroup[sno];
			while(procgroup!=NULL) {

				if(lastgroup) {

					//fprintf(stderr,"sno=%d lastgroup->end=%d procgroup->start=%d procgroup->end=%d\n",sno,lastgroup->end,procgroup->start,procgroup->end);

					//GStr bstart((int)lastgroup->end);
					//GStr bend((int)procgroup->start);
					if(!boundaryleft[lastgroup->end] && !boundaryright[procgroup->start] && (procgroup->start-lastgroup->end<=bundledist ||
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


	// each unstranded group needs to remember what proportion of stranded group it overlaps so that it can distribute reads later on -> maybe I can do this in the following while?

	while(currgroup[0]!=NULL || currgroup[1]!=NULL || currgroup[2]!=NULL) { // there are still groups to process

		int nextgr=get_min_start(currgroup); // gets the index of currgroup with the left most begining

		int grcol = currgroup[nextgr]->color;    // set smallest color for currgroup[$nextgr]

		while(equalcolor[grcol]!=grcol) {
			grcol=equalcolor[grcol];
		}
		currgroup[nextgr]->color=grcol;

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
    	//exit(0);
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
				//fprintf(stderr,"Add group=%d to bundle[%d][%d]\n",currgroup[nextgr]->grid,nextgr,bno);
				add_group_to_bundle(currgroup[nextgr],bundle[nextgr][bno],bnode[nextgr],bundledist);
			}
			else { // create new bundle
				bno=create_bundle(bundle[nextgr],currgroup[nextgr],bnode[nextgr]);
				//fprintf(stderr,"Add group=%d to new bundle[%d][%d]\n",currgroup[nextgr]->grid,nextgr,bno);
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
					//fprintf(stderr,"Add group=%d to bundle[%d:0][%d]\n",currgroup[nextgr]->grid,nextgr,bno);
					add_group_to_bundle(currgroup[nextgr],bundle[0][bno],bnode[0],bundledist); // this needs bundledist
				}
				else { // create new bundle
					bno=create_bundle(bundle[0],currgroup[nextgr],bnode[0]);
					//fprintf(stderr,"Add group=%d to new bundle[%d:0][%d]\n",currgroup[nextgr]->grid,nextgr,bno);
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
					//fprintf(stderr,"Add group=%d to new bundle[%d][%d:2]\n",currgroup[nextgr]->grid,nextgr,bno);
					add_group_to_bundle(currgroup[nextgr],bundle[2][bno],bnode[2],bundledist);
				}
				else { // create new bundle
					bno=create_bundle(bundle[2],currgroup[nextgr],bnode[2]);
					//fprintf(stderr,"Add group=%d to new bundle[%d:2][%d]\n",currgroup[nextgr]->grid,nextgr,bno);
					bundlecol[poscol]=bno;
				}
				group2bundle[2][currgroup[nextgr]->grid]=bundle[2][bno]->lastnodeid;
			}
		}

		currgroup[nextgr]=currgroup[nextgr]->next_gr;

	} // while(currgroup[0]!=NULL || currgroup[1]!=NULL || currgroup[2]!=NULL)


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

	int g=0;

	for(int b=0;b<bundle[1].Count();b++) { // these are neutral bundles that do not overlap any signed reads

		// I need to address features here too -> TODO

		bool guide_ovlp=false;
		while(g<guides.Count() && (guides[g]->exons.Count()>1 || guides[g]->end<bnode[1][bundle[1][b]->startnode]->start)) {
			g++;
		}
		// now guides[g]->end>=bnode[1][bundle[1][b]->startnode]->start
		if(g<guides.Count() && guides[g]->start<=bnode[1][bundle[1][b]->startnode]->end) guide_ovlp=true;

		if(bundle[1][b]->cov && ((bundle[1][b]->multi/bundle[1][b]->cov)<=mcov*(1-ERROR_PERC) || guide_ovlp || rawreads)) { // && (guides.Count() || adaptive || bundle[1][b]->len >= mintranscriptlen)) { // there might be small transfrags that are worth showing, but here I am ignoring them
    		// bundle might contain multiple fragments of a transcript but since we don't know the complete structure -> print only the pieces that are well represented
    		CBundlenode *currbnode=bnode[1][bundle[1][b]->startnode];
    		int t=1;
    		while(currbnode!=NULL) {
    			//int len=currbnode->end-currbnode->start+1;
    			//float cov=currbnode->cov/(currbnode->end-currbnode->start+1);

    			bool printguides=false;

    			if(!rawreads) for(int i=0;i<bnodeguides[currbnode->bid].Count();i++) {
    				int g=bnodeguides[currbnode->bid][i];
    				geneno++;
    				int glen=guides[g]->end-guides[g]->start+1;
    				if(glen && guides[g]->exons.Count()==1) {
    					RC_TData* tdata=(RC_TData*)(guides[g]->uptr);
    					tdata->in_bundle=3;
    					float gcov=(tdata->t_exons[0])->movlcount/glen;
    					// if(cov<gcov) gcov=cov; WHY DO I DO THIS?? CHECK!!!
    					CPrediction *p=new CPrediction(geneno-1, guides[g], guides[g]->start, guides[g]->end, gcov, guides[g]->strand, glen);
    					if(c_out) {
    						GStr guidecov;
    						guidecov.appendfmt("%.2f",gcov);
    						guides[g]->addAttr("coverage",guidecov.chars());
    						printCovMutex.lock();
    						guides[g]->printTranscriptGff(c_out);
    						printCovMutex.unlock();
    					}
    					GSeg exon(guides[g]->start, guides[g]->end);
    					p->exons.Add(exon);
    					p->exoncov.Add(gcov);
    					if(longreads) p->tlen=-p->tlen;
    					pred.Add(p);
    					printguides=true;
    					guidepred[g]=pred.Count()-1;
    				}
    			}

    			if(!printguides) { // && (adaptive || (cov>=readthr && len>=mintranscriptlen))) {
    				if(t==1) { geneno++;}
    				char sign='.';

    				GVec<CTrimPoint> trimpoint;
    				find_all_trims(refstart,0,currbnode->start, currbnode->end,bpcov,trimpoint); // sign should not matter as I am in a totally neutral zone

    				uint predstart=currbnode->start;
    				uint predend=currbnode->end;

    				for(int i=0;i<trimpoint.Count();i++) {
    					if(trimpoint[i].pos) {
    						if(trimpoint[i].start) {
    							int len=trimpoint[i].pos-CHI_WIN-predstart+1;
    							if(len>mintranscriptlen) {
    								float cov=get_cov(1,predstart-refstart,trimpoint[i].pos-CHI_WIN-refstart,bpcov)/len;

    								//fprintf(stderr,"Store single prediction:%d - %d with cov=%f\n",predstart, trimpoint[i].pos-CHI_WIN, cov);

    								CPrediction *p=new CPrediction(geneno-1, NULL, predstart, trimpoint[i].pos-CHI_WIN, cov, sign, len);
    								GSeg exon(predstart, trimpoint[i].pos-CHI_WIN);

    								p->exons.Add(exon);
    								if(!rawreads) p->exoncov.Add(cov);
    								if(longreads) p->tlen=-p->tlen;
    								pred.Add(p);
    								t++;
    							}
    							predstart=trimpoint[i].pos;
    						}
    						else {
    							int len=trimpoint[i].pos-predstart+1;
    							if(len>mintranscriptlen) {
    								float cov=get_cov(1,predstart-refstart,trimpoint[i].pos-refstart,bpcov)/len;

    			    				//fprintf(stderr,"Store single prediction:%d - %d with cov=%f\n",predstart, trimpoint[i].pos, cov);

    								CPrediction *p=new CPrediction(geneno-1, NULL, predstart, trimpoint[i].pos, cov, sign, len);
    								GSeg exon(predstart, trimpoint[i].pos);

    								p->exons.Add(exon);
    								if(!rawreads) p->exoncov.Add(cov);
    								if(longreads) p->tlen=-p->tlen;
    								pred.Add(p);
    								t++;
    							}
    							predstart=trimpoint[i].pos+CHI_WIN;

    						}
    					}
    				}

    				int len=predend-predstart+1;
    				if(len>mintranscriptlen) {
    					float cov=get_cov(1,predstart-refstart,predend-refstart,bpcov)/len;

    					//fprintf(stderr,"Store single prediction:%d - %d with cov=%f\n",predstart, predend, cov);

    					CPrediction *p=new CPrediction(geneno-1, NULL, predstart, predend, cov, sign, len);
    					GSeg exon(predstart, predend);

    					p->exons.Add(exon);
    					if(!rawreads) p->exoncov.Add(cov);
    					if(longreads) p->tlen=-p->tlen;
    					pred.Add(p);
    					t++;
    				}
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
    					if(bundle[sno][b]->cov) {
    						fprintf(stderr,"proc bundle[%d][%d] %f/%f is %f len=%d and %d guides\n",sno,b,
    								bundle[sno][b]->multi,bundle[sno][b]->cov,(float)bundle[sno][b]->multi/bundle[sno][b]->cov,
    								bundle[sno][b]->len,nolap);
    					}
    				}
    				*/

    				// here I can add something in stringtie to lower the mintranscript len if there are guides?

    				if(bundle[sno][b]->cov &&
    						(((bundle[sno][b]->multi/bundle[sno][b]->cov)<=mcov && bundle[sno][b]->len >= mintranscriptlen)
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
    							bundle2graph,no2gnode,transfrag,gpos,bdata,edgeno[s][b],lastgpos[s][b],guideedge); // also I need to remember graph coverages somewhere -> probably in the create_graph procedure

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


    	// because of this going throu
    	// compute probabilities for stranded bundles

    	for (int n=0;n<readlist.Count();n++) {

	  /*if(readlist[n]->unitig) { // super-reads are unpaired
    			float srcov=0;
    			for(int i=0;i<readlist[n]->segs.Count();i++)
    				srcov+=get_cov(1,readlist[n]->segs[i].start-refstart,readlist[n]->segs[i].end-refstart,bpcov)/readlist[n]->segs[i].len();
    			if(srcov) get_fragment_pattern(readlist,n,-1,readlist[n]->read_count/srcov,readgroup,merge,group2bundle,bundle2graph,graphno,edgeno,gpos,no2gnode,transfrag,tr2no,group);

    		}
    		else {*/
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
			//}
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
    			//fprintf(stderr,"Process graph[%d][%d] with %d nodes\n",s,b,graphno[s][b]);
    			if(graphno[s][b]) {

    				// include source to guide starts links
    				GVec<CGuide> guidetrf;

    				/*
    				{ // DEBUG ONLY
    					fprintf(stderr,"process refguides for s=%d b=%d edgeno=%d gno=%d lastgpos=%d guidescount=%d\n",s,b,edgeno[s][b],graphno[s][b],lastgpos[s][b],guides.Count());
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

    				if(guides.Count()) process_refguides(graphno[s][b],edgeno[s][b],gpos[s][b],lastgpos[s][b],no2gnode[s][b],transfrag[s][b],s,guidetrf,bdata);

    				GVec<int> trflong; // non-redundant long transfrags that I can use to guide the long read assemblies
    				//process transfrags to eliminate noise, and set compatibilities, and node memberships
    				process_transfrags(s,graphno[s][b],edgeno[s][b],no2gnode[s][b],transfrag[s][b],tr2no[s][b],gpos[s][b],guidetrf,pred,trflong);
    				//get_trf_long(graphno[s][b],edgeno[s][b], gpos[s][b],no2gnode[s][b],transfrag[s][b],geneno,s,pred,trflong);

    				/*
    				{ //DEBUG ONLY
    					//printTime(stderr);
    					fprintf(stderr,"There are %d nodes for graph[%d][%d]:\n",graphno[s][b],s,b);
    					for(int i=0;i<graphno[s][b];i++) {
    						fprintf(stderr,"%d (%d-%d): %f len=%d cov=%f",i,no2gnode[s][b][i]->start,no2gnode[s][b][i]->end,no2gnode[s][b][i]->cov,no2gnode[s][b][i]->len(),no2gnode[s][b][i]->cov/no2gnode[s][b][i]->len());
    						fprintf(stderr," parents:");
    						for(int j=0;j<no2gnode[s][b][i]->parent.Count();j++) fprintf(stderr," %d",no2gnode[s][b][i]->parent[j]);
    						fprintf(stderr," children:");
    						for(int j=0;j<no2gnode[s][b][i]->child.Count();j++) fprintf(stderr," %d",no2gnode[s][b][i]->child[j]);
    						fprintf(stderr," trf=");
    						for(int j=0;j<no2gnode[s][b][i]->trf.Count();j++) fprintf(stderr," %d",no2gnode[s][b][i]->trf[j]);
    						fprintf(stderr,"\n");
    					}
    					fprintf(stderr,"There are %d transfrags[%d][%d]:\n",transfrag[s][b].Count(),s,b);
    					for(int t=0;t<transfrag[s][b].Count();t++) {
    						fprintf(stderr,"%d: ",t);
    						//printBitVec(transfrag[s][b][t]->pattern);
    						fprintf(stderr," %f(%f) long=%d short=%d nodes=%d",transfrag[s][b][t]->abundance,transfrag[s][b][t]->srabund, transfrag[s][b][t]->longread,transfrag[s][b][t]->shortread,transfrag[s][b][t]->nodes.Count());
    						for(int i=0;i<transfrag[s][b][t]->nodes.Count();i++) fprintf(stderr," %d",transfrag[s][b][t]->nodes[i]);
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

    				//if(!longreads) {
    				// find transcripts now
    				if(!rawreads) geneno=find_transcripts(graphno[s][b],edgeno[s][b],gpos[s][b],no2gnode[s][b],transfrag[s][b],
    						geneno,s,guidetrf,guides,guidepred,bdata,trflong);
    				//}
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
    		fprintf(stderr,"pred[%d] (cov=%f,strand=%c):",i,pred[i]->cov,pred[i]->strand);
    		for(int j=0;j<pred[i]->exons.Count();j++) fprintf(stderr," %d-%d",pred[i]->exons[j].start,pred[i]->exons[j].end);
    		fprintf(stderr,"\n");
    	}
    }
    */

    // don't forget to clean up the allocated data here
    return(geneno);
}

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

			//if(rd.juncs.Count()) { // why is this commented?
			int m=n+1;
			while(m<nreads && readlist[m]->start<=rd.segs[0].end) {
				if(readlist[m]->nh && readlist[m]->tinfo->g==-1 && readlist[m]->juncs.Count()==rd.juncs.Count()) {
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
					add_group_to_bundle(currgroup[nextgr],bundle[0][bno],bnode[0],0); // last parameter was bundledist before
				}
				else { // create new bundle
					bno=create_bundle(bundle[0],currgroup[nextgr],bnode[0]);
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
					add_group_to_bundle(currgroup[nextgr],bundle[2][bno],bnode[2],0); // last parameter was bundledist before
				}
				else { // create new bundle
					bno=create_bundle(bundle[2],currgroup[nextgr],bnode[2]);
					bundlecol[poscol]=bno;
				}
				group2bundle[2][currgroup[nextgr]->grid]=bundle[2][bno]->lastnodeid;
			}
		}

		currgroup[nextgr]=currgroup[nextgr]->next_gr;

	} // while(currgroup[0]!=NULL || currgroup[1]!=NULL || currgroup[2]!=NULL)


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

    	if(bundle[1][b]->cov) {
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

// compute coverages and junction support here
void count_good_junctions(BundleData* bdata) {

	GList<CReadAln>& readlist = bdata->readlist;
	GList<CJunction>& junction = bdata->junction;
	GVec<float>* bpcov = bdata->bpcov;
	int refstart=bdata->start;
	int refend=bdata->end;
	bool modified=false;
	//GHash<CJunction*> jhash(false);
	GHashMap<CJunction*, CJunction*> jhash(false); //hash of pointers
	//char sbuf[20];

	if((longreads||mixedMode) && bdata->keepguides.Count()) { // there are guides to consider with longreads -> I might need to adjust the splice sites
		GPVec<GffObj>& guides = bdata->keepguides;
		GList<CJunction> gjunc(true, true, true);
		for(int g=0;g<guides.Count();g++) {
			char s=0; // unknown strand
			if(guides[g]->strand=='+') s=1; // guide on positive strand
			else if(guides[g]->strand=='-') s=-1; // guide on negative strand
			for(int i=1;i<guides[g]->exons.Count();i++) { //add_junction((int)guides[g]->exons[i-1]->end,(int)guides[g]->exons[i]->start,gjunc,s);

				int gidx=gjunc.AddedIfNew(new CJunction(guides[g]->exons[i-1]->end,guides[g]->exons[i]->start,s));
				if(gidx>-1) { // item is new
					gjunc[gidx]->guide_match=1;
					if(mixedMode) { // this should favor junctions already covered by short reads
						CJunction jn(guides[g]->exons[i-1]->end,guides[g]->exons[i]->start,s);
						int oidx=-1;
						if (!junction.Found(&jn, oidx) || junction[oidx]->nm>=junction[oidx]->nreads) {
							gjunc[gidx]->guide_match=0;
							break;
						}
					}
				}
				//fprintf(stderr,"Add guide junction=%d-%d:%d from reference %s\n",guides[g]->exons[i-1]->end,guides[g]->exons[i]->start,s,guides[g]->getID());
			}
		}

		if(gjunc.Count()) {

			//for(int i=0;i<gjunc.Count();i++) fprintf(stderr,"Guide junction=%d-%d:%d\n",gjunc[i]->start,gjunc[i]->end,gjunc[i]->strand);

			//int ssdist=1/ERROR_PERC;
			//int ssdist=longintronanchor;

			GVec<int> smodjunc; // keeps a list of all modified start junctions' indices (might keep the same junction twice)
			GVec<int> emodjunc; // keeps a list of all modified end junctions' indices (might keep the same junction twice)
			GList<CJunction> ejunction(junction);
			ejunction.setFreeItem(false);
			if(ejunction.Count()) ejunction.setSorted(juncCmpEnd);
			GList<CJunction> egjunc(gjunc);
			egjunc.setFreeItem(false);
			egjunc.setSorted(juncCmpEnd);
			int s=0;
			int e=0;
			for(int i=1;i<junction.Count();i++) {

				if(junction[i]->nm>=junction[i]->nreads){ // for all junctions -> try to see if I can correct them

					//fprintf(stderr,"check junction[%d]:%d-%d:%d rightsupport=%f nm=%f nreads=%f\n",i,junction[i]->start,junction[i]->end,junction[i]->strand,junction[i]->rightsupport,junction[i]->nm,junction[i]->nreads);

					// check start junction
					while(s<gjunc.Count() && gjunc[s]->start+sserror<junction[i]->start) s++;
					int k=s;
					int c=-1;
					int dist=1+sserror;
					while(k<gjunc.Count() && gjunc[k]->start<=junction[i]->start+sserror) {
						if(!junction[i]->strand || gjunc[k]->strand==junction[i]->strand) {
							if(gjunc[k]->start==junction[i]->start && gjunc[k]->guide_match) { // perfect match --> no need to change anything
								c=-1;
								break;
							}
							int d=dist;
							if(c<0 || gjunc[c]->guide_match==gjunc[k]->guide_match) d=abs((int)gjunc[k]->start-(int)junction[i]->start);
							if(d<dist || (!gjunc[c]->guide_match && gjunc[k]->guide_match)) {
								dist=d;
								c=k;
								smodjunc.Add(i);
							}
						}
						k++;
					}
					if(c>=0) {
						//fprintf(stderr,"...correct start of junction[%d] to %d\n",i,gjunc[c]->start);
						junction[i]->start=gjunc[c]->start;
						junction[i]->strand=gjunc[c]->strand;
						while(c<gjunc.Count() && gjunc[c]->start==junction[i]->start) {
							if(junction[i]->end==gjunc[c]->end && junction[i]->strand==gjunc[c]->strand) {
								junction[i]->guide_match=true;
								break;
							}
							c++;
						}
					}
				}
				if(ejunction[i]->nm>=ejunction[i]->nreads){ // for all junctions -> try to see if I can correct them

					//fprintf(stderr,"check ejunction[%d]:%d-%d:%d rightsupport=%f nm=%f nreads=%f\n",i,ejunction[i]->start,ejunction[i]->end,ejunction[i]->strand,ejunction[i]->rightsupport,ejunction[i]->nm,ejunction[i]->nreads);

					// check end junction
					while(e<egjunc.Count() && egjunc[e]->end+sserror<ejunction[i]->end) e++;
					int k=e;
					int c=-1;
					int dist=1+sserror;
					while(k<egjunc.Count() && egjunc[k]->end<=ejunction[i]->end+sserror) {
						if(!ejunction[i]->strand || egjunc[k]->strand==ejunction[i]->strand) {
							if(egjunc[k]->end==ejunction[i]->end  && egjunc[k]->guide_match) { // perfect match --> no need to change anything
								c=-1;
								break;
							}

							int d=dist;
							if(c<0 || egjunc[c]->guide_match==egjunc[k]->guide_match) d=abs((int)egjunc[k]->end-(int)ejunction[i]->end);
							if(d<dist || (!egjunc[c]->guide_match && egjunc[k]->guide_match)) {
								dist=d;
								c=k;
								emodjunc.Add(i);
							}
						}
						k++;
					}
					if(c>=0) {
						//fprintf(stderr,"...correct end of ejunction[%d] to %d\n",i,egjunc[c]->end);
						ejunction[i]->end=egjunc[c]->end;
						ejunction[i]->strand=egjunc[c]->strand;
						while(c<egjunc.Count() && egjunc[c]->end==ejunction[i]->end) {
							if(ejunction[i]->start==egjunc[c]->start && ejunction[i]->strand==egjunc[c]->strand) {
								ejunction[i]->guide_match=true; break;
							}
							c++;
						}
					}
					if(s==gjunc.Count() && e==egjunc.Count()) break;
				}
			}

			if(smodjunc.Count()||emodjunc.Count()) {
				modified=true;
				for(int i=0;i<smodjunc.Count();i++) {
					int j=smodjunc[i]; // junction that I modified --> try to see if there are other equal junctions
					//sprintf(sbuf, "%p", junction[j]);
					const CJunction* jp=jhash[junction[j]];
					if(!jp) { // did not process junction before
						//fprintf(stderr,"smodified junction[%d]:%d-%d:%d %p nreads=%f\n",j,junction[j]->start,junction[j]->end,junction[j]->strand,junction[j],junction[j]->nreads);
						s=j-1;
						GVec<int> equal;
						while(s>=0 && junction[s]->start==junction[j]->start) {
							if(junction[s]->end==junction[j]->end && junction[s]->strand==junction[j]->strand) equal.Add(s);
							s--;
						}
						s=j+1;
						while(s<junction.Count() && junction[s]->start==junction[j]->start) {
							if(junction[s]->end==junction[j]->end && junction[s]->strand==junction[j]->strand) equal.Add(s);
							s++;
						}
						if(equal.Count()) { // junction j is equal to other junctions
							for(s=0;s<equal.Count();s++) {
								//sprintf(sbuf, "%p", junction[equal[s]]);
								//jhash.Add(sbuf,junction[j]);
								CJunction* jct=junction[equal[s]];
								//fprintf(stderr,"...equal to junction[%d]:%d-%d:%d %p nreads=%f\n",equal[s],jct->start,jct->end,jct->strand,jct,jct->nreads);
								jhash.Add(jct, junction[j]);
								if(mixedMode) { // I can trust the reads coming from mixed data but not so much from long reads
									junction[j]->nreads+=jct->nreads;
									junction[j]->nm+=jct->nm;
									jct->nreads=0;
								}
								jct->strand=0;
								jct->guide_match=false;
							}
						}
					}
				}
				for(int i=0;i<emodjunc.Count();i++) {
					int j=emodjunc[i];
					//sprintf(sbuf, "%p", ejunction[j]);
					CJunction* jp=jhash[ejunction[j]];
					if(!jp) { // did not process junction before
						//fprintf(stderr,"emodified ejunction[%d]:%d-%d:%d %p nreads=%f\n",j,ejunction[j]->start,ejunction[j]->end,ejunction[j]->strand,ejunction[j],ejunction[j]->nreads);
						s=j-1;
						GVec<int> equal;
						while(s>=0 && ejunction[s]->end==ejunction[j]->end) {
							if(ejunction[s]->start==ejunction[j]->start && ejunction[s]->strand==ejunction[j]->strand) equal.Add(s);
							s--;
						}
						s=j+1;
						while(s<ejunction.Count() && ejunction[s]->end==ejunction[j]->end) {
							if(ejunction[s]->start==ejunction[j]->start && ejunction[s]->strand==ejunction[j]->strand) equal.Add(s);
							s++;
						}
						if(equal.Count()) { // junction j is equal to other junctions
							for(s=0;s<equal.Count();s++) {
								//sprintf(sbuf, "%p", ejunction[equal[s]]);
								//jhash.Add(sbuf,ejunction[j]);
								//ejunction[equal[s]]->strand=0;
								//ejunction[equal[s]]->guide_match=false;
								CJunction* ej=ejunction[equal[s]];
								jhash.Add(ej, ejunction[j]);
								//fprintf(stderr,"...equal to ejunction[%d]:%d-%d:%d %p nreads=%f\n",equal[s],ej->start,ej->end,ej->strand,ej,ej->nreads);
								if(mixedMode) { // I can trust the reads coming from mixed data but not so much from long reads
									ejunction[j]->nreads+=ej->nreads;
									ejunction[j]->nm+=ej->nm;
									ej->nreads=0;
								}
								ej->strand=0;
								ej->guide_match=false;
							}
						}
					}
				}
			}
			junction.Sort();
		}
		gjunc.Clear();
	}

	for(int s=0;s<3;s++) bpcov[s].Resize(refend-refstart+3);

	GVec<int> unstranded; // remembers unstranded reads

	for(int n=0;n<readlist.Count();n++) {
		CReadAln & rd=*(readlist[n]);
		float rdcount=rd.read_count;

		int nex=rd.segs.Count();

		if(!rd.unitig) add_read_to_cov(readlist,n,bpcov,refstart);
		else if(rdcount>1) rdcount=1;

		GVec<uint> leftsup;
		GVec<uint> rightsup;
		uint maxleftsupport=0;
		uint maxrightsupport=0;

		//int sno=(int)rd.strand+1; // 0(-),1(.),2(+)
		//if(nex>1) fprintf(stderr,"Process spliced read[%d] with cov=%f and sno=%d: ",n,rdcount,sno);
		for(int i=0;i<nex;i++) {
			//if(nex>1) fprintf(stderr," %d-%d",rd.segs[i].start,rd.segs[i].end);
			if(i) {
				//fprintf(stderr,":%d",rd.juncs[i-1]->strand);

				if(modified) { // see if read uses modified junction -> correct it
					//sprintf(sbuf, "%p", rd.juncs[i-1]);
					CJunction* jp=jhash[rd.juncs[i-1]];
					if(jp) {
						if(rd.segs[i-1].start>jp->start || rd.segs[i].end<jp->end) {

							if(rd.segs[i-1].start<=jp->start) rd.segs[i-1].end=jp->start;
							if(rd.segs[i].end>=jp->end) rd.segs[i].start=jp->end;

							rd.juncs[i-1]->nreads-=rd.read_count;
							//if(rd.juncs[i-1]->nreads<ERROR_PERC) rd.juncs[i-1]->strand=0; // this approach removes a perfectly valid (guide) junction which might not be the desired effect
							rd.juncs[i-1]=junction[0];

						}
						else {
							rd.juncs[i-1]=jp;
							if(!rd.strand) rd.strand=jp->strand;
							//fprintf(stderr," [correct rd from %d-%d to %d-%d]",rd.segs[i-1].end,rd.segs[i].start,jp->start,jp->end);
							if(rd.segs[i-1].start<=jp->start) rd.segs[i-1].end=jp->start;
							if(rd.segs[i].end>=jp->end) rd.segs[i].start=jp->end;
						}

					}
					else {
						if(rd.segs[i-1].start>rd.juncs[i-1]->start || rd.segs[i].end<rd.juncs[i-1]->end) {
							if(rd.segs[i-1].end!=rd.juncs[i-1]->start && rd.segs[i-1].start<=rd.juncs[i-1]->start) rd.segs[i-1].end=rd.juncs[i-1]->start;
							if(rd.segs[i].start!=rd.juncs[i-1]->end && rd.segs[i].end>=rd.juncs[i-1]->end) rd.segs[i].start=rd.juncs[i-1]->end;
							rd.juncs[i-1]->nreads-=rd.read_count;
							//if(rd.juncs[i-1]->nreads<ERROR_PERC) rd.juncs[i-1]->strand=0; // this approach removes a perfectly valid (guide) junction which might not be the desired effect
							rd.juncs[i-1]=junction[0];
						}
						else {
							//if(rd.segs[i-1].end!=rd.juncs[i-1]->start || rd.segs[i].start!=rd.juncs[i-1]->end) fprintf(stderr," [chg rd from %d-%d to %d-%d]",rd.segs[i-1].end,rd.segs[i].start,rd.juncs[i-1]->start,rd.juncs[i-1]->end);
							if(rd.segs[i-1].end!=rd.juncs[i-1]->start && rd.segs[i-1].start<=rd.juncs[i-1]->start) rd.segs[i-1].end=rd.juncs[i-1]->start;
							if(rd.segs[i].start!=rd.juncs[i-1]->end && rd.segs[i].end>=rd.juncs[i-1]->end) rd.segs[i].start=rd.juncs[i-1]->end;
						}
					}
				}
				/* DEL AWARE*/
				else if(!rd.juncs[i-1]->strand && (longreads || mixedMode) && (rd.segs[i-1].end!=rd.juncs[i-1]->start || rd.segs[i].start!=rd.juncs[i-1]->end)){ // see if I need to adjust read start/ends due to junction having deletions around it

					CJunction *nj=NULL;
					CJunction jn(rd.juncs[i-1]->start, rd.juncs[i-1]->end, rd.strand); // if the deletion makes sense then keep it
					int oidx=-1;
					if (junction.Found(&jn, oidx)) {
						nj=junction.Get(oidx);
						nj->nreads+=rd.juncs[i-1]->nreads;
						nj->mm+=rd.juncs[i-1]->mm;
						nj->nm+=rd.juncs[i-1]->nm;
						rd.juncs[i-1]=nj;
						// adjust start/end of read
						rd.segs[i-1].end=rd.juncs[i-1]->start;
						rd.segs[i].start=rd.juncs[i-1]->end;
					}
					else { // restore the correct junction
						CJunction jn(rd.segs[i-1].end, rd.segs[i].start, rd.strand);
						if (junction.Found(&jn, oidx)) {
							nj=junction.Get(oidx);
							nj->nreads+=rd.juncs[i-1]->nreads;
							nj->mm+=rd.juncs[i-1]->mm;
							nj->nm+=rd.juncs[i-1]->nm;
							rd.juncs[i-1]=nj;
						}
						else { // new junction
							rd.juncs[i-1]->strand=rd.strand;
							rd.juncs[i-1]->start=rd.segs[i-1].end;
							rd.juncs[i-1]->end=rd.segs[i].start;
						}
					}
				}


				if(rd.segs[i-1].len()>maxleftsupport) maxleftsupport=rd.segs[i-1].len();
				if(rd.segs[nex-i].len()>maxrightsupport) maxrightsupport=rd.segs[nex-i].len();
				leftsup.Add(maxleftsupport);
				rightsup.Add(maxrightsupport);
				//rd.juncs[i-1]->nreads+=rd.read_count;
				//if(rd.unitig) rd.juncs[i-1]->guide_match=true; // v7 this might be a little too much!
			}
			//if(!rd.unitig) cov_edge_add(bpcov,sno,rd.segs[i].start-refstart,rd.segs[i].end+1-refstart,rd.read_count);

		}
		//if(nex>1) fprintf(stderr," With anchors: ");
		for(int i=1;i<nex;i++) {
			//if(!rd.juncs[i-1]->strand) { rd.juncs[i-1]->strand = rd.strand; }
			uint anchor=junctionsupport;
			//if(rd.juncs[i-1]->len()>longintron || rd.juncs[i-1]->nreads-rd.juncs[i-1]->mm<junctionthr || rd.juncs[i-1]->nreads-rd.juncs[i-1]->nm<junctionthr) {
			if(rd.juncs[i-1]->len()>longintron || (!longreads && rd.juncs[i-1]->nreads-rd.juncs[i-1]->nm<junctionthr && !rd.juncs[i-1]->mm)) {
				//if(rd.juncs[i-1]->len()>verylongintron) { if(anchor<verylongintronanchor) anchor=verylongintronanchor; } // I want to use a longer anchor for long introns to believe them
				//else
					if(anchor<longintronanchor) anchor=longintronanchor;

			}

			//if(rd.unitig) anchor=1; // v7 also **** if not trimming involved in super-read creation then comment this; unitigs should be trimmed so I will accept them as anchors

			//if(leftsup[i-1]>=anchor && rightsup[nex-i-1]>=anchor) rd.juncs[i-1]->nreads_good+=rd.read_count;
			if(leftsup[i-1]>=anchor) { // support only comes from spliced reads that are bigger than the anchor
				rd.juncs[i-1]->leftsupport+=rdcount;
				if(rightsup[nex-i-1]>=anchor) {
					rd.juncs[i-1]->rightsupport+=rdcount;
					rd.juncs[i-1]->nreads_good+=rdcount;
					//rd.juncs[i-1]->nreads_good+=rd.read_count*rd.nh; // try to see if this helps
					//if(leftsup[i-1]>=anchor+1/ERROR_PERC && rightsup[nex-i-1]>=anchor+1/ERROR_PERC) rd.juncs[i-1]->strong=true;
				}
			}
			else if(rightsup[nex-i-1]>=anchor) {
				rd.juncs[i-1]->rightsupport+=rdcount;
			}
			//fprintf(stderr," %d(%d-%d)[%f-%f][%d][%f][%f]",anchor,leftsup[i-1],rightsup[nex-i-1],rd.juncs[i-1]->leftsupport,rd.juncs[i-1]->rightsupport,rd.juncs[i-1]->strand,rd.juncs[i-1]->nm,rd.juncs[i-1]->nreads);
		}

		if((mixedMode || longreads) && !rd.strand && nex>1) unstranded.Add(n);

		//if(nex>1) fprintf(stderr,"\n");
	}


	//fprintf(stderr,"Compute coverages:\n");

	// this code enssures that I can find interval coverages very fast
	int m=int((bpcov[1].Count()-1)/BSIZE);
	//fprintf(stderr,"m=%d refstart=%d refend=%d bpcount=%d\n",m,refstart,refend,bpcov[1].Count());
	int k=0;
	float prev_val[3]={0,0,0};
	float prev_sum[3]={0,0,0};
	while(k<m) {
		int end=(k+1)*BSIZE;


		//fprintf(stderr,"end=%d\n",end);
		for(int i=k*BSIZE;i<end;i++) {

			for(int s=0;s<3;s++) {
				//fprintf(stderr,"(1)bpcov[%d][%d]=%f ",s,i,bpcov[s][i]);
				bpcov[s][i]+=prev_val[s];
				if(bpcov[s][i]<0) bpcov[s][i]=0;
				prev_val[s]=bpcov[s][i];
				//fprintf(stderr,"(2)bpcov[%d][%d]=%f ",s,i,bpcov[s][i]);
				bpcov[s][i]+=prev_sum[s];
				prev_sum[s]=bpcov[s][i];
				//fprintf(stderr,"(2)bpcov[%d][%d]=%f ",s,i,bpcov[s][i]);
				//fprintf(stderr,"bpPos[%d][%d]: prev_val=%f prev_sum=%f\n",s,i,prev_val[s],prev_sum[s]);
			}
		}
		for(int s=0;s<3;s++) prev_sum[s]=0;
		k++;
	}

	for(int i=k*BSIZE;i<bpcov[1].Count();i++)
		for(int s=0;s<3;s++) {
			bpcov[s][i]+=prev_val[s];
			if(bpcov[s][i]<0) bpcov[s][i]=0;
			prev_val[s]=bpcov[s][i];
			bpcov[s][i]+=prev_sum[s];
			prev_sum[s]=bpcov[s][i];
			//fprintf(stderr,"bpPos[%d][%d]: cov=%f sumbpcov=%f\n",s,i,prev_val[s],prev_sum[s]);
		}

	if(longreads || mixedMode) {
		for(int n=0;n<unstranded.Count();n++) {
			CReadAln & rd=*(readlist[unstranded[n]]);
			if(rd.start>=(uint)refstart && (int)(rd.end-(uint)refstart)<bpcov[1].Count()) {
				int nex=rd.segs.Count();
				float pos=0;
				float neg=0;
				if(longreads) {
					pos=get_cov(2,rd.segs[0].start-refstart,rd.segs[0].end-refstart,bpcov);
					neg=get_cov(0,rd.segs[0].start-refstart,rd.segs[0].end-refstart,bpcov);
				}
				bool posjunc=false;
				bool negjunc=false;
				for(int i=1;i<nex;i++) { // some junctions don't get an assigned strand by buggy aligners
					if(longreads) {
						pos+=get_cov(2,rd.segs[i].start-refstart,rd.segs[i].end-refstart,bpcov);
						neg+=get_cov(0,rd.segs[i].start-refstart,rd.segs[i].end-refstart,bpcov);
					}
					if(rd.juncs[i-1]->strand>0) posjunc=true;
					else if(rd.juncs[i-1]->strand<0) negjunc=true;
					else {
						int oidx=-1;
						CJunction jn(rd.juncs[i-1]->start, rd.juncs[i-1]->end, 1);
						if (junction.Found(&jn, oidx)) {
							posjunc=true;
						}
						jn.strand=-1;
						if (junction.Found(&jn, oidx)) {
							negjunc=true;
						}
					}
				}
				if(posjunc && !negjunc) {
					rd.strand=1;
				}
				else if(negjunc && !posjunc) {
					rd.strand=-1;
				}
				else {
					if(neg<1) neg=0;
					if(pos<1) pos=0;
					if(neg>pos) {rd.strand=-1;}
					else if(pos>neg) { rd.strand=1;}
				}
				//fprintf(stderr,"read strand is:%d for pos=%f neg=%f\n",rd.strand,pos,neg);
				if(rd.strand) {
					for(int i=1;i<nex;i++)
						if(!rd.juncs[i-1]->strand) {
							// first check to see if the junction already exists
							int oidx=-1;
							CJunction *nj=NULL;
							CJunction jn(rd.juncs[i-1]->start, rd.juncs[i-1]->end, rd.strand);
							if (junction.Found(&jn, oidx)) {
								nj=junction.Get(oidx);
								//fprintf(stderr,"found strand junction at %p",nj);
								nj->nreads+=rd.juncs[i-1]->nreads;
								nj->nreads_good+=rd.juncs[i-1]->nreads_good;
								nj->nm+=rd.juncs[i-1]->nm;
								nj->mm+=rd.juncs[i-1]->mm;
								nj->leftsupport+=rd.juncs[i-1]->leftsupport;
								nj->rightsupport+=rd.juncs[i-1]->rightsupport;
								rd.juncs[i-1]=nj;
							}
							else rd.juncs[i-1]->strand = rd.strand;
					}
				}
			}
		}
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

		count_good_junctions(bundle);

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
	if(abs(a->tlen) < abs(b->tlen)) return -1;
	if(abs(a->tlen) > abs(b->tlen)) return 1;
	return 0;
}

int predcovCmp(const pointer p1, const pointer p2) { // sort from highest to lowest coverage
	CPrediction *a=(CPrediction*)p1;
	CPrediction *b=(CPrediction*)p2;
	if(a->cov < b->cov) return 1;
	if(a->cov > b->cov) return -1;
	return 0;
}


int predordCmp(const pointer p1, const pointer p2) { // sort from highest to lowest coverage
	CPred *a=(CPred*)p1;
	CPred *b=(CPred*)p2;
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

//bool is_exon_above_frac(CInterval *maxcov,uint start,uint end,float cov,int exonno) {
bool is_exon_above_frac(CInterval *maxcov,uint start,uint end,float cov) {

	CInterval *lastinterv=NULL;
	while(maxcov && start>=maxcov->pos) {
		lastinterv=maxcov;
		maxcov=maxcov->next;
	}
	//if(lastinterv && ((exonno==1 && cov<lastinterv->val) || cov<isofrac*lastinterv->val)) return(false); // I need to deal with single exons too here
	if(lastinterv &&  cov<isofrac*lastinterv->val) return(false); // I need to deal with single exons too here
	while(maxcov && end>=maxcov->pos) {
		//if((exonno==1 && cov<maxcov->val) || cov<isofrac*maxcov->val) return(false);
		if(cov<isofrac*maxcov->val) return(false);
		maxcov=maxcov->next;
	}
	return(true);
}

void delete_Cinterval(CInterval *interv){
	if(interv) {
		if(interv->next) delete_Cinterval(interv->next);
		delete interv;
	}
}


uint min(uint n1,uint n2) {
	if(n1<n2) return(n1);
	return(n2);
}

uint max(uint n1,uint n2) {
	if(n1<n2) return(n2);
	return n1;
}

/*
bool retainedintron(GList<CPrediction>& pred,int n1,int n2,GVec<GBitVec>& lowintron) {
	int j=1;
	for(int i=1;i<pred[n1]->exons.Count();i++) {
		if(j>pred[n2]->exons.Count()-2) return false;
		if(lowintron[n1][i-1]) {
			while(j<pred[n2]->exons.Count()-1 && pred[n2]->exons[j].end<pred[n1]->exons[i].start) j++; // now pred[n2]->exons[j].end>=pred[n1]->exons[i].start
			if(j<pred[n2]->exons.Count()-1 && pred[n2]->exons[j].start<=pred[n1]->exons[i-1].end) return true;
		}
	}
	return false;
}
*/

int retainedintron(GList<CPrediction>& pred,int n1,int n2,GVec<GBitVec>& lowintron) {

  float frac=ERROR_PERC;
  if(mixedMode && isofrac<frac && pred[n2]->tlen<0 && pred[n2]->cov>DROP/ERROR_PERC) frac=isofrac;

	int j=0;
	for(int i=1;i<pred[n1]->exons.Count();i++) {
		if(j>pred[n2]->exons.Count()-1) return(0);
		if(lowintron[n1][i-1]) { // pred[n1] has lowintron i-1
			if(j==pred[n2]->exons.Count()-1 && pred[n2]->cov<frac*pred[n1]->cov && pred[n2]->exons[j].start<=pred[n1]->exons[i-1].end)
				return(1);
			while(j<pred[n2]->exons.Count() && pred[n2]->exons[j].end<pred[n1]->exons[i].start) j++; // now pred[n2]->exons[j].end>=pred[n1]->exons[i].start
			if(!j && pred[n2]->cov<frac*pred[n1]->cov) return(1);
			if(j<pred[n2]->exons.Count() && pred[n2]->exons[j].start<=pred[n1]->exons[i-1].end) {
				if(j && j<pred[n2]->exons.Count()-1) {
					/*fprintf(stderr,"middle exon i=%d j=%d pred[n2=%d]=%f < pred[n1=%d]=%f for intron %d-%d in n1 and exon %d-%d in n2\n",i,j,
							n2,pred[n2]->cov,n1,pred[n1]->cov,pred[n1]->exons[i-1].end,pred[n1]->exons[i].start,pred[n2]->exons[j].start,pred[n2]->exons[j].end);*/
					return(2); // middle exon
				}
				else if(pred[n2]->cov<frac*pred[n1]->cov) {
					/*fprintf(stderr,"i=%d j=%d pred[n2=%d]=%f < pred[n1=%d]=%f for intron %d-%d in n1 and exon %d-%d in n2\n",i,j,
							n2,pred[n2]->cov,n1,pred[n1]->cov,pred[n1]->exons[i-1].end,pred[n1]->exons[i].start,pred[n2]->exons[j].start,pred[n2]->exons[j].end);*/
					return(1);
				}
			}
		}
	}
	return(0);
}

bool included_pred(GList<CPrediction>& pred,int n1,int n2,uint refstart=0,GVec<float>* bpcov = NULL,bool reference=true) { // check if the small prediction is included in the larger prediction

	int intronextension=2*longintronanchor;
	if(longreads) intronextension=CHI_WIN; // more aggressive inclusion rules for long reads

	if(pred[n1]->start > pred[n2]->end || pred[n2]->start>pred[n1]->end) return false;

	int big=n1;
	int small=n2;

	if(pred[n1]->exons.Count()<pred[n2]->exons.Count()) {
		big=n2;
		small=n1;
	}

	if(reference && pred[small]->t_eq) { // this is a reference guide included into another transcript -> only keep one of them
		if(pred[big]->t_eq  && pred[small]->t_eq!=pred[big]->t_eq) return false;
		if(pred[small]->exons.Count()!=pred[big]->exons.Count()) return false;
	}


	int bex=0;
	while(bex<pred[big]->exons.Count()) {
		if(pred[small]->exons[0].start>pred[big]->exons[bex].end) bex++;
		else { // now pred[small]->exons[0].start <= pred[big]->exons[bex].end
			if(pred[small]->exons[0].end<pred[big]->exons[bex].start) return false; // no overlap
			// now: pred[small]->exons[0].end>=pred[big]->exons[bex].start
			int sex=0;

			//if(pred[small]->exons.Count()>1) while(bex<pred[big]->exons.Count() && pred[big]->exons[bex].end<pred[small]->exons[0].end) bex++;
			//if(bex==pred[big]->exons.Count()) return false;

			while(sex<pred[small]->exons.Count() && bex<pred[big]->exons.Count()) {
				if(sex==pred[small]->exons.Count()-1) { // I am at end of small pred

					if(bex==pred[big]->exons.Count()-1) return true; // same intron structure and there is overlap; if small pred is single exon then I didn't check if it extends past the end of bex which is problematic since above I do

					if(sex && pred[small]->exons[sex].end>pred[big]->exons[bex+1].start) return false; // might be able to eliminate some transcripts here if the drop is big -> I should probably just comment this because below I take care of it; or make it as an else to the if below

					//if(sex && pred[small]->exons[sex].end>pred[big]->exons[bex].end+longintronanchor) { // small prediction extends well into intron of big prediction
					if(sex && pred[small]->exons[sex].end>pred[big]->exons[bex].end+intronextension) { // small prediction extends well into intron of big prediction
						uint endin=pred[small]->exons[sex].end > pred[big]->exons[bex+1].start ? pred[big]->exons[bex+1].start : pred[small]->exons[sex].end;
						float covintron=get_cov(1,pred[big]->exons[bex].end+1-refstart,endin-refstart,bpcov)/(endin-pred[big]->exons[bex].end);
						if(covintron>singlethr) {
							uint startex=pred[small]->exons[sex].start > pred[big]->exons[bex].start ? pred[small]->exons[sex].start : pred[big]->exons[bex].start;
							float covexon=get_cov(1,startex-refstart,pred[big]->exons[bex].end-refstart,bpcov)/(pred[big]->exons[bex].end-startex+1);
							if(covintron > ERROR_PERC*covexon && pred[big]->exons[bex].end-longintronanchor>refstart) { // return false
								float exondrop=get_cov(1,pred[big]->exons[bex].end-longintronanchor+1-refstart,pred[big]->exons[bex].end-refstart,bpcov);
								float introndrop=get_cov(1,pred[big]->exons[bex].end+1-refstart,pred[big]->exons[bex].end+longintronanchor-refstart,bpcov);
								if(introndrop>ERROR_PERC*exondrop)
									return false;
							}
						}
					}

					return true;
				}

				// sex is not last exon in small prediction but overlaps bex
				if(bex==pred[big]->exons.Count()-1) return false; // small pred extends past big pred
				if(pred[small]->exons[sex].end != pred[big]->exons[bex].end) return false;

				//if(!sex && bex && pred[small]->exons[sex].start<pred[big]->exons[bex].start-longintronanchor) {  // first exon of small prediction extends within intron of big prediction
				if(!sex && bex && pred[small]->exons[sex].start<pred[big]->exons[bex].start-intronextension) {  // first exon of small prediction extends within intron of big prediction
					uint startin=pred[small]->exons[sex].start > pred[big]->exons[bex-1].end ? pred[small]->exons[sex].start : pred[big]->exons[bex-1].end;
					float covintron=get_cov(1,startin-refstart,pred[big]->exons[bex].start-1-refstart,bpcov)/(pred[big]->exons[bex].start-startin);
					if(covintron>singlethr) {
						uint endex=pred[small]->exons[sex].end > pred[big]->exons[bex].end ? pred[big]->exons[bex].end : pred[small]->exons[sex].end;
						float covexon=get_cov(1,pred[big]->exons[bex].start-refstart,endex-refstart,bpcov)/(endex-pred[big]->exons[bex].start+1);
						if(covintron > ERROR_PERC*covexon && pred[big]->exons[bex].start-longintronanchor>=refstart) { // return false
							float introndrop=get_cov(1,pred[big]->exons[bex].start-longintronanchor-refstart,pred[big]->exons[bex].start-1-refstart,bpcov);
							float exondrop=get_cov(1,pred[big]->exons[bex].start-refstart,pred[big]->exons[bex].start+longintronanchor-1-refstart,bpcov);
							if(introndrop>ERROR_PERC*exondrop)
								return false;
						}
					}
				}


				bex++;
				sex++;
				if(pred[small]->exons[sex].start != pred[big]->exons[bex].start) return false;
			}
			return false;
		}
	}

	return false;
}

void update_cov(GList<CPrediction>& pred,int big,int small,float frac=1) { // small gets included into big

	if(pred[big]->strand=='.') pred[big]->strand=pred[small]->strand;

	if(pred[small]->t_eq && !pred[big]->t_eq) { // this can only happen if big pred has the same number of exons
		pred[big]->tlen=pred[small]->tlen;
		pred[big]->exons[0].start=pred[small]->exons[0].start;
		pred[big]->start=pred[small]->start;
		pred[big]->exons.Last().end=pred[small]->exons.Last().end;
		pred[big]->end=pred[small]->end;
		pred[big]->t_eq=pred[small]->t_eq;
	}

	int bex=0; // first exon in big pred that overlaps small pred
	while(pred[small]->exons[0].start>pred[big]->exons[bex].end) bex++;

	if(pred[big]->cov<pred[small]->cov && !pred[big]->t_eq && !bex && pred[small]->exons.Count()>1 && pred[big]->exons[0].start<pred[small]->exons[0].start) { // adjust start to account for trimming -> this messes up coverage of big pred though
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


void update_overlap(GList<CPrediction>& pred,int p,int e,GVec<CExon>& node,GVec<bool>& overlap) {

	int n=pred.Count();
	for(int i=0;i<node.Count();i++) if(!overlap[n*node[i].predno+p]) { // there is no overlap detected yet

		// exon e in prediction p overlaps exon ei in prediction pi
		int pi=node[i].predno;
		int ei=node[i].exonno;

		overlap[n*pi+p]=true;

		// overlap is not significant in the following cases
		if(!e && pred[p]->exons[0].end>=pred[pi]->end) { // I am at the beginning of prediction p
			int len=pred[pi]->end-pred[p]->start+1;
			if(len<ERROR_PERC*abs(pred[p]->tlen)) overlap[n*pi+p]=false; // prediction p doesn't overlap prediction pi by a considerable amount
		}

		if(overlap[n*pi+p] && !ei && pred[pi]->exons[0].end>=pred[p]->end) { // I am at the beginning of prediction pi
			int len=pred[p]->end-pred[pi]->start+1;
			if(len<ERROR_PERC*abs(pred[pi]->tlen)) overlap[n*pi+p]=false; // prediction pi doesn't overlap prediction p by a considerable amount
		}

		if(overlap[n*pi+p] && e==pred[p]->exons.Count()-1 && pred[p]->exons[e].start<=pred[pi]->start) { // I am at the end of prediction p
			int len=pred[p]->end-pred[pi]->start+1;
			if(len<ERROR_PERC*abs(pred[p]->tlen)) overlap[n*pi+p]=false; // prediction p doesn't overlap prediction pi by a considerable amount
		}

		if(overlap[n*pi+p] && ei==pred[pi]->exons.Count()-1 && pred[pi]->exons[ei].start<=pred[p]->start) { // I am at the end of prediction pi
			int len=pred[pi]->end-pred[p]->start+1;
			if(len<ERROR_PERC*abs(pred[pi]->tlen)) overlap[n*pi+p]=false; // prediction pi doesn't overlap prediction p by a considerable amount
		}

	}
}


CMaxIntv *add_exon_to_maxint(CMaxIntv *maxint,uint start,uint end,int p,int e,float c,GList<CPrediction>& pred,GVec<bool>& overlap) {

	CMaxIntv *prevmaxint=NULL;
	while(maxint && start>maxint->end) {
		prevmaxint=maxint;
		maxint=maxint->next; // skip maxint until I get to the one that overlaps my exon or after it
	}

	CExon ex=CExon(p,e,c);

	CMaxIntv *nextmaxint;
	if(maxint) { // start <= maxint->end, start>prexmaxint->end ; there might be overlap with maxint

		if(end<maxint->start) { // --> no overlap -> exon is in between nodes of maxint -> I need to keep previous maxint for this
			if(!prevmaxint) { // this shouldn't happen due to maxint construction
				GError("No overlap with exon in printing results!\n");
			}
			nextmaxint=new CMaxIntv(start,end);
			nextmaxint->node.Add(ex);
			nextmaxint->next=maxint;
			prevmaxint->next=nextmaxint;
			return nextmaxint;
		}
		else { // end>=maxint->start --> there is overlap

			if(maxint->start<start) { // some region of maxint comes before exon
				nextmaxint=new CMaxIntv(maxint->node,start,maxint->end,maxint->next);
				maxint->next=nextmaxint;
				maxint->end=start-1;
			}
			else if(start<maxint->start) { // this does nothing to maxint
				nextmaxint=new CMaxIntv(start,maxint->start-1);
				nextmaxint->next=maxint;
				if(prevmaxint) prevmaxint->next=nextmaxint;
			}
			else { // start==maxint->start and this does nothing to maxint
				nextmaxint=maxint;
			}

			update_overlap(pred,p,e,nextmaxint->node,overlap);
			nextmaxint->node.Add(ex);

			if(end<nextmaxint->end) { // nextmaxint now ends where previous maxint ended
				maxint=new CMaxIntv(nextmaxint->node,end+1,nextmaxint->end,nextmaxint->next);
				nextmaxint->end=end;
				maxint->node.Pop();
				nextmaxint->next=maxint;
			}
			else if(end>nextmaxint->end) {
				maxint=nextmaxint;
				while(maxint->next && end>=maxint->next->start) {
					if(maxint->end+1<maxint->next->start) {
						CMaxIntv *newintv=new CMaxIntv(maxint->end+1,maxint->next->start-1);
						newintv->node.Add(ex);
						newintv->next=maxint->next;
						maxint->next=newintv;
						maxint=newintv;
					}
					update_overlap(pred,p,e,maxint->next->node,overlap);
					maxint->next->node.Add(ex);
					maxint=maxint->next;
				}

				// end<maxint->next->start or maxint->next is NULL

				if(end>maxint->end) {
					CMaxIntv *newintv=new CMaxIntv(maxint->end+1,end);
					newintv->node.Add(ex);
					newintv->next=maxint->next;
					maxint->next=newintv;
				}
				else if(end<maxint->end) {
					CMaxIntv *newintv=new CMaxIntv(maxint->node,end+1,maxint->end,maxint->next);
					newintv->node.Pop();
					maxint->next=newintv;
					maxint->end=end;
				}
			}
		}
	}
	else { // start > prevmaxint->end and maxint is NULL --> no overlap here
		nextmaxint=new CMaxIntv(start,end);
		nextmaxint->node.Add(ex);
		prevmaxint->next=nextmaxint;
	}

	return nextmaxint;
}

bool intronic(GList<CPrediction>& pred,int m,int M) {
	int i=0;
	while(i<pred[M]->exons.Count() && pred[m]->start>=pred[M]->exons[i].start) i++;
	if(i==pred[M]->exons.Count() || !i) return false;

	// now pred[m]->start<pred[M]->exons[i].start
	if(pred[m]->exons[0].end<pred[M]->exons[i-1].end || pred[m]->exons.Last().start>pred[M]->exons[i].start) return false; // full first|last exon inside exon of big prediction

	if(pred[m]->start<=pred[M]->exons[i-1].end) {
		int len=pred[M]->exons[i-1].end-pred[m]->start+1;
		if(len>ERROR_PERC*abs(pred[m]->tlen)) return false; // overlap too big for prediction to be considered intronic
	}
	if(pred[m]->end>=pred[M]->exons[i].start) {
		int len=pred[m]->end-pred[M]->exons[i].start+1;
		if(len>ERROR_PERC*abs(pred[m]->tlen)) return false;
	}

	return true;
}

bool transcript_overlap(GList<CPrediction>& pred,int n1,int n2) {
	if(pred[n1]->end>=pred[n2]->start && pred[n2]->end>=pred[n1]->start) { // there is overlap --> check if it's substantial by looking at first/last exons as in update_overlap
		if(pred[n1]->exons.Last().start<=pred[n2]->exons[0].end || pred[n2]->exons.Last().start<=pred[n1]->exons[0].end) return false;
		return true;
	}
	return false;
}

int print_predcluster(GList<CPrediction>& pred,int geneno,GStr& refname,
		GVec<CGene>& refgene, GHash<int>& hashgene, GVec<CGene>& predgene, BundleData* bundleData,bool checkincomplete) {

	uint runoffdist=200;
	if(longreads) runoffdist=0;
	if(bundledist>runoffdist) runoffdist=bundledist;

	int npred=pred.Count();
	GVec<bool> overlap;
	overlap.Resize(npred*npred-npred);

	GVec<CPred> predord;
	//CPred p(0,pred[0]->cov);
	CPred p(0,abs(pred[0]->tlen)*pred[0]->cov); // priority based on number of bases covered -> this should a give a boost to predictions that are longer

	predord.Add(p);
	GVec<int> color;
	color.cAdd(0);

	GVec<int> incomplete;
	GHash<bool> guideintron;

	GVec<GBitVec> lowintron;
	GBitVec intron;

	if(pred[0]->exons.Count()>1) {
		intron.resize(pred[0]->exons.Count()-1,false);
		if(checkincomplete && !pred[0]->t_eq && pred[0]->mergename.is_empty()) { // incomplete transcript
			incomplete.cAdd(0);
		}
	}

	GVec<float>* bpcov = bundleData->bpcov;

	float intronfrac=DROP-ERROR_PERC;
	if(longreads && isofrac<intronfrac) intronfrac=isofrac;

	// build maxIntv
	CMaxIntv *maxint=new CMaxIntv(pred[0]->exons[0].start,pred[0]->exons[0].end);

	if(!pred[0]->t_eq) {
		if(eonly) pred[0]->flag=false;
		else if(pred[0]->exons.Count()==1 && pred[0]->strand!='.') { // neutral strand if most reads are neutral
			int s=0;
			if(pred[0]->strand=='+') s=2;
			float totalcov=get_cov(1,pred[0]->start-bundleData->start,pred[0]->end-bundleData->start,bpcov);
			float strandcov=get_cov(s,pred[0]->start-bundleData->start,pred[0]->end-bundleData->start,bpcov);
			if(strandcov<ERROR_PERC*totalcov) pred[0]->strand='.';
		}
	}

	float excov=pred[0]->cov;
	if(!longreads) excov*=abs(pred[0]->tlen); // it was len() before but that doesn't make any sense
	CExon ex(0,0,excov); // this keeps the exon flow based on per bp coverage
	//CExon ex(0,0,pred[0]->exoncov[0]*pred[0]->exons[0].len()); // this keeps the exon flow based on read coverage (elen factor)
	//pred[0]->exoncov[0]=0;
	maxint->node.Add(ex);
	CMaxIntv *nextmaxint=maxint;
	bool exist=true;
	GStr id("", 32);
	for(int j=1;j<pred[0]->exons.Count();j++) {
		if(checkincomplete && pred[0]->t_eq) {
			id.assign((int)pred[0]->exons[j-1].end);
			id+=pred[0]->strand;
			id+=(int)pred[0]->exons[j].start;
			bool *gi=guideintron[id.chars()];
			if(!gi) guideintron.Add(id.chars(), exist);
		}
		float introncov=get_cov(1,pred[0]->exons[j-1].end+1-bundleData->start,pred[0]->exons[j].start-1-bundleData->start,bpcov)/(pred[0]->exons[j].start-pred[0]->exons[j-1].end-1);
		if(introncov) {
		  if(introncov<1 || (!longreads && introncov<singlethr)) intron[j-1]=1;
		  else {
				  float exoncov=(get_cov(1,pred[0]->exons[j-1].start-bundleData->start,pred[0]->exons[j-1].end-bundleData->start,bpcov)+
						  get_cov(1,pred[0]->exons[j].start-bundleData->start,pred[0]->exons[j].end-bundleData->start,bpcov))/(pred[0]->exons[j-1].len()+pred[0]->exons[j].len());
		    if(introncov<exoncov*intronfrac) intron[j-1]=1;
				  else {
						  // left drop
						  int start=pred[0]->exons[j-1].end-longintronanchor+1-bundleData->start;
						  if(start<0) start=0;
						  int end=pred[0]->exons[j-1].end+longintronanchor-bundleData->start;
						  if(end>bpcov->Count()-2) end=bpcov->Count()-2;
						  exoncov=get_cov(1,start,pred[0]->exons[j-1].end-bundleData->start,bpcov);
						  introncov=get_cov(1,pred[0]->exons[j-1].end+1-bundleData->start,end,bpcov);
							  if(introncov<exoncov*intronfrac) intron[j-1]=1;
		    	else { // right drop
								  int start=pred[0]->exons[j].start-longintronanchor-bundleData->start;
								  if(start<0) start=0;
								  int end=pred[0]->exons[j].start+longintronanchor-1-bundleData->start;
								  if(end>bpcov->Count()-2) end=bpcov->Count()-2;
								  introncov=get_cov(1,start,pred[0]->exons[j].start-1-bundleData->start,bpcov);
								  exoncov=get_cov(1,pred[0]->exons[j].start-bundleData->start,end,bpcov);
		    		if(introncov<exoncov*intronfrac) intron[j-1]=1;
							  }
						  }
			//fprintf(stderr,"pred[%d] intron[%d]:%d-%d introncov=%f exoncov=%f\n",0,j-1,pred[0]->exons[j-1].end,pred[0]->exons[j].start,introncov,exoncov);
				  }
			  }
		nextmaxint=add_exon_to_maxint(nextmaxint,pred[0]->exons[j].start,pred[0]->exons[j].end,0,j,excov,pred,overlap); // per bp coverage
		//nextmaxint=add_exon_to_maxint(nextmaxint,pred[0]->exons[j].start,pred[0]->exons[j].end,0,j,pred[0]->exoncov[j]*pred[0]->exons[j].len(),pred,overlap); // per read coverage
		//pred[0]->exoncov[j]=0;
	}

	lowintron.Add(intron);

	nextmaxint=maxint;
	for(int n=1;n<npred;n++) {

		if(!pred[n]->t_eq) {
			if(eonly) pred[n]->flag=false;
			else if(pred[n]->exons.Count()==1 && pred[n]->strand!='.') { // neutral strand if most reads are neutral
				int s=0;
				if(pred[n]->strand=='+') s=2;
				float totalcov=get_cov(1,pred[n]->start-bundleData->start,pred[n]->end-bundleData->start,bpcov);
				float strandcov=get_cov(s,pred[n]->start-bundleData->start,pred[n]->end-bundleData->start,bpcov);
				if(strandcov<ERROR_PERC*totalcov) pred[n]->strand='.';
			}
		}

		color.Add(n);

		excov=pred[n]->cov;
		if(!longreads) excov*=abs(pred[n]->tlen);
		//CPred p(n,pred[n]->cov); // priority based on cov/bp
		CPred p(n,abs(pred[n]->tlen)*pred[n]->cov); // priority based on number of bases covered
		predord.Add(p);
		nextmaxint=add_exon_to_maxint(nextmaxint,pred[n]->exons[0].start,pred[n]->exons[0].end,n,0,excov,pred,overlap); // per bp coverage
		//nextmaxint=add_exon_to_maxint(nextmaxint,pred[n]->exons[0].start,pred[n]->exons[0].end,n,0,pred[n]->exoncov[0]*pred[n]->exons[0].len(),pred,overlap); // read coverage
		//pred[n]->exoncov[0]=0;
		intron.clear();
		intron.resize(pred[n]->exons.Count()-1,false);
		if(checkincomplete && !pred[n]->t_eq && pred[n]->exons.Count()>1 && pred[n]->mergename.is_empty()) { // incomplete transcript
			incomplete.Add(n);
		}

		CMaxIntv *nextintv=nextmaxint;
		for(int j=1;j<pred[n]->exons.Count();j++) {
			if(checkincomplete && pred[n]->t_eq) {
				//GStr id((int)pred[n]->exons[j-1].end);
				id.assign((int)pred[n]->exons[j-1].end);
				id+=pred[n]->strand;
				id+=(int)pred[n]->exons[j].start;
				bool *gi=guideintron[id.chars()];
				if(!gi) guideintron.Add(id.chars(),exist);
			}
			float introncov=get_cov(1,pred[n]->exons[j-1].end+1-bundleData->start,pred[n]->exons[j].start-1-bundleData->start,bpcov)/(pred[n]->exons[j].start-pred[n]->exons[j-1].end-1);
			if(introncov) {
			  if(introncov<1 || (!longreads && introncov<singlethr)) intron[j-1]=1;
			  else {
					  float exoncov=(get_cov(1,pred[n]->exons[j-1].start-bundleData->start,pred[n]->exons[j-1].end-bundleData->start,bpcov)+
							  get_cov(1,pred[n]->exons[j].start-bundleData->start,pred[n]->exons[j].end-bundleData->start,bpcov))/(pred[n]->exons[j-1].len()+pred[n]->exons[j].len());

					if(introncov<exoncov*intronfrac) intron[j-1]=1;
					  else {
							  // left drop
							  int start=pred[n]->exons[j-1].end-longintronanchor+1-bundleData->start;
							  if(start<0) start=0;
							  int end=pred[n]->exons[j-1].end+longintronanchor-bundleData->start;
							  if(end>bpcov->Count()-2) end=bpcov->Count()-2;
							  exoncov=get_cov(1,start,pred[n]->exons[j-1].end-bundleData->start,bpcov);
							  introncov=get_cov(1,pred[n]->exons[j-1].end+1-bundleData->start,end,bpcov);
								  if(introncov<exoncov*intronfrac) intron[j-1]=1;
						else { // right drop
									  int start=pred[n]->exons[j].start-longintronanchor-bundleData->start;
									  if(start<0) start=0;
									  int end=pred[n]->exons[j].start+longintronanchor-1-bundleData->start;
									  if(end>bpcov->Count()-2) end=bpcov->Count()-2;
									  introncov=get_cov(1,start,pred[n]->exons[j].start-1-bundleData->start,bpcov);
									  exoncov=get_cov(1,pred[n]->exons[j].start-bundleData->start,end,bpcov);
							if(introncov<exoncov*intronfrac) intron[j-1]=1;
								  }
							  }
						  }
					  }
			nextintv=add_exon_to_maxint(nextintv,pred[n]->exons[j].start,pred[n]->exons[j].end,n,j,excov,pred,overlap); // per bp coverage
			//nextintv=add_exon_to_maxint(nextintv,pred[n]->exons[j].start,pred[n]->exons[j].end,n,j,pred[n]->exoncov[j]*pred[n]->exons[j].len(),pred,overlap); // read coverage
			//pred[n]->exoncov[j]=0;
		}

		lowintron.Add(intron);

		/*
		{ // DEBUG ONLY
			fprintf(stderr,"After adding pred %d\n",n);
			CMaxIntv *intv=maxint;
			while(intv) {
				fprintf(stderr,"...average coverage for interval:%d-%d=%f (%f perbp)\n",intv->start,intv->end,bpcov[1][intv->end-bundleData->start]-bpcov[1][intv->start-1-bundleData->start],(bpcov[1][intv->end-bundleData->start]-bpcov[1][intv->start-1-bundleData->start])/(intv->end-intv->start+1));
				intv=intv->next;
			}
		}
		*/

	}

	if(checkincomplete) {
		for(int i=0;i<incomplete.Count();i++) {
			int n=incomplete[i];
			bool eliminate=true;
			for(int j=1;j<pred[n]->exons.Count();j++) {
				//GStr id((int)pred[n]->exons[j-1].end);
				id.assign((int)pred[n]->exons[j-1].end);
				id+=pred[n]->strand;
				id+=pred[n]->exons[j].start;
				bool *gi=guideintron[id.chars()];
				if(!gi) {
					eliminate=false; break;
				}
			}
			if(eliminate) {
				pred[n]->flag=false;
				//fprintf(stderr,"falseflag: incomplete pred[%d]\n",n);
			}
		}
	}


	//GVec<int> bettercov[npred];
	GVec<int> ivec;
	GVec< GVec<int> > bettercov(npred, ivec);

	predord.Sort(predordCmp); // sort predictions from the most highest coverage to lowest
	for(int i=0;i<predord.Count()-1;i++) {

		if(pred[predord[i].predno]->flag) {
			int n1=predord[i].predno;
			//fprintf(stderr,"first consideration of pred[%d] with cov=%f and strand=%c\n",n1,pred[n1]->cov,pred[n1]->strand);
			for(int j=i+1;j<predord.Count();j++) if(pred[predord[j].predno]->flag) {
				int n2=predord[j].predno; // always predord.cov(n1)>=predord.cov(n2)
				int m = n1;
				int M = n2;
				//fprintf(stderr,"...vs prediction n2=%d with cov=%f and strand=%c\n",n2,pred[n2]->cov,pred[n2]->strand);
				if(m>M) { m=n2; M=n1;}
				if(overlap[npred*m+M]) { // n1 could get eliminated by some other genes for instance if it has fewer exons and it's included in a bigger one
				  //fprintf(stderr,"overlap\n");
					if(!pred[n2]->t_eq) { // this is not a known gene -> only then I can eliminate it
						if(pred[n1]->t_eq && pred[n2]->cov<ERROR_PERC*pred[n1]->cov) { // more strict about novel predictions if annotation is available
							pred[n2]->flag=false;
							continue;
						}
						uint anchor=longintronanchor;
						if(pred[n2]->tlen<0 && pred[n2]->cov>readthr) anchor=junctionsupport;
						if(pred[n2]->exons[0].len()<anchor || pred[n2]->exons.Last().len()<anchor) { // VAR8: leave this one out
							//fprintf(stderr,"falseflag: pred[%d] n2=%d has low first/last exon\n",n2,n2);
							pred[n2]->flag=false;
						}
						else if(retainedintron(pred,n1,n2,lowintron)) {
							//if(ret>1 || pred[n2]->cov<ERROR_PERC*pred[n1]->cov) {
							//fprintf(stderr,"falseflag: pred[%d] n2=%d has low intron coverage\n",n2,n2);
							pred[n2]->flag=false;
							//}
						}
						else if(pred[n1]->strand != pred[n2]->strand) {

							//fprintf(stderr,"diff strand\n");
							if(pred[n2]->exons.Count()==1 || (pred[n1]->tlen<0 && pred[n2]->tlen>0 && pred[n2]->cov<1/ERROR_PERC)) { // why was I restricting it to single exons?
								//fprintf(stderr,"falseflag: ...strand elmination of pred[%d] n2=%d by n1=%d\n",n2,n2,n1);
								pred[n2]->flag=false;
							}
							else if(!pred[n1]->t_eq && pred[n1]->exons.Count()==1) {
								if(pred[n1]->cov < pred[n2]->cov+singlethr) {
									//fprintf(stderr,"flaseflag: pred[%d] n1=%d with 1 exon eliminated by n2=%d\n",n1,n1,n2);
									pred[n1]->flag=false;
									break;
								}
								/*else {
								  float totalcov=get_cov(1,pred[n1]->start-bundleData->start,pred[n1]->end-bundleData->start,bundleData->bpcov);
								  if(totalcov*DROP>pred[n1]->cov*pred[n1]->len())
								    pred[n1]->flag=false; break;
								    }*/
							}
							else if(pred[n2]->exons.Count()<=pred[n1]->exons.Count() && included_pred(pred,n1,n2,(uint)bundleData->start,bpcov)) {
								//fprintf(stderr,"falseflag: pred[%d] n2=%d is included in n1=%d\n",n2,n2,n1);
								pred[n2]->flag=false;
							}
							else if(pred[n2]->tlen>0 && !mixedMode){
	    						bool lowcov=true;
	    						for(int e=0;e<pred[n2]->exons.Count();e++) {
	    							float totalcov=get_cov(1,pred[n2]->exons[e].start-bundleData->start,pred[n2]->exons[e].end-bundleData->start,bundleData->bpcov);
	    							if(pred[n2]->exoncov[e]*pred[n2]->exons[e].len()>totalcov/2) {
	    								lowcov=false;
	    								break;
	    							}
	    						}
	    						if(lowcov) {
	    							pred[n2]->flag=false;
	    							//fprintf(stderr,"falseflag: pred[%d] n2=%d has lowcov\n",n2,n2);
	    						}
	    					}

						} // end pred[n1]->strand != pred[n2]->strand
	    				//else if(pred[n2]->cov<isofrac*pred[n1]->cov) pred[n2]->flag=false;  // I am only considering inclusions here since this might change allocations
	    				else if(pred[n2]->exons.Count()<=pred[n1]->exons.Count() && pred[n1]->cov>pred[n2]->cov*DROP && included_pred(pred,n1,n2,(uint)bundleData->start,bpcov)) {
	    					pred[n2]->flag=false;
	    					//fprintf(stderr,"falseflag: ...included elmination of pred[%d] n2=%d by n1=%d\n",n2,n2,n1);
	    				}
	    				else if(!pred[n1]->t_eq && n1 && ((pred[n1]->tlen>0 && pred[n1]->exons.Count()<=2 && pred[n2]->cov>lowisofrac*pred[n1]->cov) || pred[n1]->cov < pred[n2]->cov+singlethr) && pred[n1]->exons.Count()<pred[n2]->exons.Count() && included_pred(pred,n2,n1,(uint)bundleData->start,bpcov)) {
	    					pred[n1]->flag=false;
	    					//fprintf(stderr,"falseflag: ...included elmination of pred[%d] n1=%d by n2=%d\n",n1,n1,n2);
	    					break;
	    				}
	    			}
	    			else if(!pred[n1]->t_eq && n1 && ((pred[n1]->tlen>0 && pred[n1]->exons.Count()<=2 && pred[n2]->cov>lowisofrac*pred[n1]->cov) || pred[n1]->cov < pred[n2]->cov+singlethr) && pred[n1]->exons.Count()<=pred[n2]->exons.Count() && included_pred(pred,n2,n1,(uint)bundleData->start,bpcov)) {
	    				//fprintf(stderr,"falseflag: ...included elmination of pred[%d] n1=%d(%f) by n2=%d(%f)\n",n1,n1,pred[n1]->cov,n2,pred[n2]->cov);
	    				pred[n1]->flag=false;
	    				break;
	    			}
	    			else if(pred[n2]->tlen<0 && !pred[n2]->t_eq && pred[n2]->cov < DROP && pred[n2]->exons.Count()<=pred[n1]->exons.Count() && included_pred(pred,n1,n2,(uint)bundleData->start,bpcov,false)) {
	    				pred[n2]->flag=false;
	    				//fprintf(stderr,"falseflag: pred[%d] n2=%d is included in n1=%d\n",n2,n2,n1);
	    			}
	    		}
	    		//else if(!pred[n2]->t_eq && pred[n2]->exons.Count()==1 && pred[n2]->cov<DROP*pred[n1]->cov && pred[n1]->start<pred[n2]->start && pred[n2]->end<pred[n1]->end) {
	    		else if(!pred[n2]->t_eq && pred[n2]->exons.Count()==1 && ((pred[n2]->tlen>0 && pred[n2]->cov<pred[n1]->cov) || pred[n2]->cov<isofrac*pred[n1]->cov) && pred[n1]->start<=pred[n2]->start && pred[n2]->end<=pred[n1]->end) {
	    			//fprintf(stderr,"falseflag: ...single exon elmination of pred[%d] n2=%d by n1=%d\n",n2,n2,n1);
	    			pred[n2]->flag=false; // delete intronic prediction if single exon and at realtively low coverage
	    		}
	    		else if(pred[n1]->tlen>0) {
	    			if(!pred[n1]->t_eq && pred[n1]->exons.Count()==1 && pred[n1]->cov<pred[n2]->cov && pred[n2]->start<=pred[n1]->start && pred[n1]->end<=pred[n2]->end) {
	    				//fprintf(stderr,"falseflag: ...single exon elmination of pred[%d] n1=%d by n2=%d\n",n1, n1,n2);
	    				pred[n1]->flag=false;
	    				break;
	    			}
	    			else if(pred[n2]->tlen>0 && !pred[n2]->t_eq) {
	    				if(transcript_overlap(pred,n1,n2)) {
	    				for(int p=0;p<bettercov[n2].Count();p++) {
	    					int m = n1;
	    					int M = bettercov[n2][p];
	    					if(m>M) { m=bettercov[n2][p]; M=n1;}
	    					if(!overlap[npred*m+M]) {
	    						pred[n2]->flag=false;
	    						//fprintf(stderr,"falseflag: pred[%d] eliminated because of overlap to %d and %d\n",n2,m,M);
	    						break;
	    					}
	    				}
	    				if(pred[n2]->flag) {
	    						if((mixedMode || pred[n1]->strand==pred[n2]->strand) && ((mixedMode && pred[n2]->cov<singlethr) || pred[n2]->cov<pred[n1]->cov*ERROR_PERC) && pred[n1]->start<=pred[n2]->start && pred[n2]->end<=pred[n1]->end && intronic(pred,n2,n1)) { // n2 is an intronic prediction to n1
							  //fprintf(stderr,"falseflag: eliminate pred[%d] is intronic into pred[%d]\n",n2,n1);
	    						pred[n2]->flag=false;
	    					}
	    					else bettercov[n2].Add(n1);
	    				}
	    			}
	    				else if(mixedMode && pred[n2]->strand!=pred[n1]->strand && pred[n2]->cov<singlethr) {
	    					if(pred[n1]->start<=pred[n2]->end && pred[n2]->start<=pred[n1]->end) pred[n2]->flag=false; // antisense short read overlap
	    				}
	    		}
	    		}
	    		else if(pred[n2]->tlen>0  && !pred[n2]->t_eq ) { // && pred[n2]->strand!=pred[n1]->strand){ // now pred[n1]->tlen<0; n1 & n2 are antisense
	    			if(pred[n1]->start<=pred[n2]->end && pred[n2]->start<=pred[n1]->end) {// n1 & n2 overlap but not within the exons
	    				if(pred[n2]->cov<1/ERROR_PERC || (pred[n2]->strand==pred[n1]->strand && pred[n2]->cov<ERROR_PERC*pred[n1]->cov)) pred[n2]->flag=false; // mixed mode doesn't accept low covered anti-senses
	    			}
	    		}
	    		//if(!pred[n1]->flag) fprintf(stderr,"pred[%d] eliminated\n",n1);
	    		//if(!pred[n2]->flag) fprintf(stderr,"pred[%d] eliminated\n",n2);
	    	}
	    }
	}

	/*
	// eliminate included predictions
	GVec<bool> keep;
	keep.Resize(npred,false);
	nextmaxint=maxint;
	while(nextmaxint) {
		GVec<CPred> exord; // only remember those predictions that have an intron
		for(int i=0;i<nextmaxint->node.Count();i++) if(pred[nextmaxint->node[i].predno]->flag) {
			int n=nextmaxint->node[i].predno;
			int e=nextmaxint->node[i].exonno;
			if(e<pred[n]->exons.Count()-1 && pred[n]->exons[e].end==nextmaxint->end) { // intron coming here
				CPred ex(i,pred[nextmaxint->node[i].predno]->tlen*pred[nextmaxint->node[i].predno]->cov);
				exord.Add(ex); //exord only remembers the non-false predictions that have introns starting here
			}
		}
		if(exord.Count()) { // there are introns starting here
			GVec<int> introns; // remembers end of introns starting here
			exord.Sort(predordCmp);
			for(int j=0;j<exord.Count();j++) {
				int i=exord[j].predno;
				int n=nextmaxint->node[i].predno;
				int e=nextmaxint->node[i].exonno;
				if(e<pred[n]->exons.Count()-1 && pred[n]->exons[e].end==nextmaxint->end) { // intron coming here
					// check if new intron
					bool newintron=true;
					int end=pred[n]->exons[e+1].start;
					for(int k=0;k<introns.Count();k++) {
						if(introns[k]==end) {
							newintron=false;
							break;
						}
					}
					if(newintron) {
						introns.Add(end); // add new intron end here
						keep[n]=true;
					}
				}
			}
		}
		nextmaxint=nextmaxint->next;
	}
	*/

	if(longreads) { // eliminate predictions under threshold without doing the coverage adjustment
		nextmaxint=maxint;
		while(nextmaxint) {
			if((int)nextmaxint->start-bundleData->start<bpcov[1].Count()) {
				GVec<CPred> exord;
				for(int i=0;i<nextmaxint->node.Count();i++) if(pred[nextmaxint->node[i].predno]->flag) {
					CPred ex(i,abs(pred[nextmaxint->node[i].predno]->tlen)*pred[nextmaxint->node[i].predno]->cov);
					exord.Add(ex); //exord only remembers the non-false predictions
				}
				if(exord.Count()) {
					exord.Sort(predordCmp); // sort exons from the most highest coverage to lowest
					int i=exord[0].predno;
					int n=nextmaxint->node[i].predno;
					float usedcov[2]={0,0};
					float multicov[2]={0,0};
					int s=0;
					if(pred[n]->strand=='+') s=1;
					usedcov[s]=pred[n]->cov;
					if(pred[n]->exons.Count()>1) multicov[s]=usedcov[s];
					//fprintf(stderr,"maxint=%d-%d has first pred[%d] with cov=%f usedcov=%f\n",nextmaxint->start,nextmaxint->end,n,pred[n]->cov,usedcov[s]);

					for(int j=1;j<exord.Count();j++) {
						i=exord[j].predno;
						n=nextmaxint->node[i].predno;
						bool longunder=false;
						if(pred[n]->strand=='+') s=1;
						else s=0;
						if(!pred[n]->t_eq) {
							float isofraclong=isofrac;
							if(pred[n]->cov>CHI_WIN) isofraclong=isofrac*ERROR_PERC*DROP; // more abundant transcript have higher error tollerance
							else if(pred[n]->cov>CHI_THR) isofraclong=isofrac*DROP;
							//else if(pred[n]->cov>1/ERROR_PERC) isofraclong=isofrac*DROP;
							if(pred[n]->exons.Count()>1) {
								if((!multicov[s] && pred[n]->cov<isofraclong*usedcov[s] && pred[n]->cov<DROP/ERROR_PERC) ||
										pred[n]->cov<isofraclong*multicov[s]) longunder=true;
							}
							else if(pred[n]->cov<isofraclong*usedcov[s]) longunder=true;
						}
						if(longunder) {
							//fprintf(stderr,"flaseflag: pred[%d] has low cov=%f percentage from usedcoverage=%f\n",n,pred[n]->cov,usedcov[s]);
							pred[n]->flag=false;
						}
						else {
							usedcov[s]+=pred[n]->cov;
							if(pred[n]->exons.Count()>1) multicov[s]+=pred[n]->cov;
							//fprintf(stderr,"maxint=%d-%d add pred[%d]->cov=%f to usedcov=%f\n",nextmaxint->start,nextmaxint->end,n,pred[n]->cov,usedcov[s]);
						}
					}
				}
			}
			nextmaxint=nextmaxint->next;
		}
	}
	else { // recompute coverages -> mostly for short reads; long reads should be more reliable in the way coverages are estimated
		bool adjust=true;
		while(adjust) { // by here I eliminated inclusions and intronic single exons
			adjust=false;

			for(int n=0;n<npred;n++) if(pred[n]->flag) {
				/*if(!pred[n]->t_eq && !keep[n] && pred[n]->cov<lowcov) pred[n]->flag=false; // this was commented and it was almost there Pr 29.8 for 534291
		      	  else*/
				for(int i=0;i<pred[n]->exons.Count();i++) if(pred[n]->tlen>0)
					pred[n]->exoncov[i]=0;
			}

			nextmaxint=maxint;
			while(nextmaxint) {
				if((int)nextmaxint->start-bundleData->start<bpcov[1].Count()) {
					float covsum=0;
					GVec<CPred> exord;
					for(int i=0;i<nextmaxint->node.Count();i++) if(pred[nextmaxint->node[i].predno]->flag) {
				        //CPred ex(i,pred[nextmaxint->node[i].predno]->cov);
						CPred ex(i,abs(pred[nextmaxint->node[i].predno]->tlen)*pred[nextmaxint->node[i].predno]->cov);
						exord.Add(ex); //exord only remembers the non-false predictions

						//covsum+=pred[nextmaxint->node[i].predno]->cov; // priority based on cov/bp
						//covsum+=pred[nextmaxint->node[i].predno]->tlen*pred[nextmaxint->node[i].predno]->cov; // priority based on number of bases covered
						covsum+=nextmaxint->node[i].exoncov; // priority based on exoncov
					}

					if(covsum) {
						float totalcov=get_cov(1,nextmaxint->start-bundleData->start,nextmaxint->end-bundleData->start,bpcov);
						exord.Sort(predordCmp); // sort exons from the most highest coverage to lowest
						int i=exord[0].predno;
						int n=nextmaxint->node[i].predno;
						int e=nextmaxint->node[i].exonno;
						//float exoncov=totalcov*pred[n]->tlen*pred[n]->cov/covsum; // priority based on number of bases covered
						//float exoncov=totalcov*pred[n]->cov/covsum; // priority based on per bp cov
						float exoncov=totalcov*nextmaxint->node[i].exoncov/covsum; // priority based on per bp cov
						bool shortmax=false;
						if(pred[n]->tlen>0) {
						  pred[n]->exoncov[e]+=exoncov;
						  if(mixedMode && exord.Count()==1) {
						    if(pred[n]->exons.Count()==1 && pred[n]->cov<singlethr/DROP) pred[n]->flag=false; // stricter criteria for random single exon short-read genes
						    else if(pred[n]->exons.Count()==2 && pred[n]->cov<singlethr*DROP) pred[n]->flag=false; // stricter criteria for random two-exon short-read genes
						  }
						  shortmax=true;
						}
						//fprintf(stderr,"pred[%d]->exoncov[%d]=%f exoncov=%f totalcov=%f covsum=%f\n",n,e,pred[n]->exoncov[e],exoncov,totalcov,covsum);
						float usedcov=pred[n]->cov;

						//fprintf(stderr,"maxint=%d-%d has first pred[%d] with cov=%f usedcov=%f\n",nextmaxint->start,nextmaxint->end,n,pred[n]->cov,usedcov);

						for(int j=1;j<exord.Count();j++) {
							i=exord[j].predno;
							n=nextmaxint->node[i].predno;
							float allowedfrac=isofrac;
							if(mixedMode) {
								if(pred[n]->cov>1/ERROR_PERC) allowedfrac*=DROP;
								else if(pred[n]->tlen>0) {
								  allowedfrac/=DROP; // stricter criteria for short read data
								  if(shortmax) allowedfrac=DROP; // very strict criteria if short reads are maximal
								}
							}
							if(!pred[n]->t_eq && (pred[n]->cov<allowedfrac*usedcov || (pred[n]->cov<lowcov && pred[n]->cov<lowisofrac*usedcov) ||
									(pred[n]->exons.Count()==1 && pred[n]->cov<usedcov) || pred[n]->cov<1 ||
									(pred[n]->exons.Count()==2 && (pred[n]->cov<lowcov || (pred[n]->geneno>=0 && pred[n]->cov<ERROR_PERC*usedcov))))) {
								adjust=true;
								//fprintf(stderr,"flaseflag: pred[%d] has low cov=%f percentage from usedcoverage=%f < %f\n",n,pred[n]->cov,usedcov,allowedfrac*usedcov);
								pred[n]->flag=false;
							}
							else {
								e=nextmaxint->node[i].exonno;
								//exoncov=totalcov*pred[n]->cov/covsum; // priority based on per bp cov
								//float exoncov=totalcov*pred[n]->tlen*pred[n]->cov/covsum; // priority based on number of bases covered
								exoncov=totalcov*nextmaxint->node[i].exoncov/covsum; // priority based on per bp cov
								if(pred[n]->tlen>0) pred[n]->exoncov[e]+=exoncov;
								usedcov+=pred[n]->cov;
								//fprintf(stderr,"maxint=%d-%d add pred[%d]->cov=%f to usedcov=%f\n",nextmaxint->start,nextmaxint->end,n,pred[n]->cov,usedcov);
							}
						}
					}
				}
				nextmaxint=nextmaxint->next;
			}

			//for(int p=0;p<npred;p++) {
			//	int n=predord[p].predno;
			for(int n=0;n<npred;n++) {
				if(pred[n]->flag && pred[n]->tlen>0) {

					// compute prediction coverage first
					pred[n]->cov=0;
					for(int i=0;i<pred[n]->exons.Count();i++) {
						pred[n]->cov+=pred[n]->exoncov[i];
						pred[n]->exoncov[i]/=pred[n]->exons[i].len();
					}
					pred[n]->cov/=pred[n]->tlen;

					if(!pred[n]->t_eq) {
						if(pred[n]->cov<1 || (pred[n]->exons.Count()==1 && pred[n]->cov<singlethr)) {
							//if(pred[n]->cov<readthr || (pred[n]->exons.Count()==1 && pred[n]->cov<singlethr) ) {
							pred[n]->flag=false;
							//fprintf(stderr,"falseflag: pred[%d]->cov=%f eliminated due to low coverage\n",n,pred[n]->cov);
							adjust=true;
						}
					}

					//predord[p].cov=pred[n]->tlen*pred[n]->cov; ** use this
					//predord[p].cov=pred[n]->cov;
					//fprintf(stderr,"coverage of pred[%d]=%f\n",n,pred[n]->cov);
				}
			}
		}
	}

	// eliminate possible polymerase run-offs
	for(int i=0;i<npred;i++) if(pred[i]->flag && pred[i]->exons.Count()==1 && !pred[i]->t_eq) { // only check prediction that can be deleted
		int p=i-1;
		while(p>-1 && !pred[p]->flag) p--;
		int c=i+1;
		while(c<npred && !pred[c]->flag) c++;

		uint runoff=runoffdist;
		//if(mixedMode && pred[i]->tlen<0) runoff=0;

		if((p<0 || pred[p]->end<pred[i]->start) && (c==npred || pred[c]->start>pred[i]->end)) { // pred[i] alone between two predictions
			//if(pred[i]->exons.Count()==1) { // single exon prediction -> check if I can delete it
				if(p>-1 && pred[p]->exons.Count()>1 && pred[p]->end+runoff>=pred[i]->start) { // previous prediction is nearby
					float exoncov=get_cov(1,pred[p]->exons.Last().start-bundleData->start,pred[p]->exons.Last().end-bundleData->start,bpcov);
					if(pred[i]->cov<exoncov+singlethr) {
						pred[i]->flag=false;
						//fprintf(stderr,"falseflag: elim pred[%d] as possible polymerase runoff\n",i);
					}
				}
				if(pred[i]->flag && c<npred && pred[c]->exons.Count()>1 && pred[i]->end+runoff>=pred[c]->start) { // next prediction is nearby
					float exoncov=get_cov(1,pred[c]->exons[0].start-bundleData->start,pred[c]->exons[0].end-bundleData->start,bpcov);
					if(pred[i]->cov<exoncov+singlethr) {
						pred[i]->flag=false;
						//fprintf(stderr,"falseflag: elim pred[%d] as possible polymerase runoff\n",i);
					}
				}
			//}
			// else if pred[i]->cov<singlethr -> here I would do the stiching of partial genes
		}
	}

	// assign gene numbers
	GVec<int> genes; // for each prediction remembers its geneno
	genes.Resize(npred,-1);
	GVec<int> transcripts; // for each gene remembers how many transcripts were printed

	//pred.Sort();
	for(int i=0;i<npred;i++)
	  if(pred[i]->flag) {
		 //TODO: this eliminates e.g. 0.98 cov prediction based on long read that otherwise covers all the junctions!
		 //  =>  implement a jcov metric (junction coverage) which should supersede base coverage ?
		if ( pred[i]->cov<1 ||
				(!pred[i]->t_eq && (pred[i]->cov<readthr || (mixedMode && guided && pred[i]->cov<singlethr)))) {
		//if ( !pred[i]->t_eq && (pred[i]->cov<readthr || (mixedMode && guided && pred[i]->cov<singlethr)) ) {
			pred[i]->flag=false;
			//fprintf(stderr,"falseflag: elim pred[%d] due to low cov=%f\n",i,pred[i]->cov);
			continue;
		}
		//fprintf(stderr,"check pred i=%d with end=%d and next start=%d\n",i,pred[i]->end,pred[i+1]->start);
		int ci=color[i];
		while(ci!=color[ci]) { ci=color[ci];color[i]=ci;}
		int j=i+1;
		while(j<npred && pred[i]->end>=pred[j]->start) {
			//fprintf(stderr,"... check pred j=%d\n",j);
			if(pred[j]->flag && pred[i]->strand==pred[j]->strand && overlap[npred*i+j]) {
				if(!pred[j]->t_eq && pred[j]->cov<readthr) {
					pred[j]->flag=false;
					//fprintf(stderr,"falseflag: elim pred[%d] due to low cov=%f\n",j,pred[j]->cov);
					j++;
					continue;
				}
				//fprintf(stderr,"pred %d overlaps pred %d\n",j,i);
				int cj=color[j];
				while(cj!=color[cj]) { cj=color[cj];color[j]=cj;}
				if(cj<ci) {
					color[ci]=cj;
					color[i]=cj;
				}
				else if(ci<cj){
					color[cj]=ci;
					color[j]=ci;
				}
			}
			j++;
		}
	}

	for(int n=0;n<npred;n++) if(pred[n]->flag && genes[n]==-1) { // no gene was assigned to this prediction
		int cn=color[n];
		genes[n]=transcripts.Count();
		transcripts.cAdd(0);
		while(cn!=color[cn]) { cn=color[cn];color[n]=cn;}
		uint currentend=pred[n]->end;
		int m=n+1;
		while(m<npred && pred[m]->start<=currentend) {
			if(pred[m]->flag) {
				int cm=color[m];
				while(cm!=color[cm]) { cm=color[cm];color[m]=cm;}
				if(cn==cm) {
					if(genes[m]>-1) GError("Overlapping predictions with different gene number!\n");
					genes[m]=genes[n];
				}
				if(pred[m]->end>currentend) currentend=pred[m]->end;
			}
			m++;
		}
	}

	for(int n=0;n<npred;n++) if(pred[n]->flag){

		/*
		{ // DEBUG ONLY
			  fprintf(stderr,"\nprint prediction %d with cov=%f len=%d",n,pred[n]->cov,pred[n]->tlen);
			  if(pred[n]->flag) fprintf(stderr," with true flag");
			  fprintf(stderr," with geneno=%d and exons:",pred[n]->geneno);
			  for(int i=0;i<pred[n]->exons.Count();i++) fprintf(stderr," cov=%f len=%d",pred[n]->exoncov[i],pred[n]->exons[i].len());
			  fprintf(stderr,"\n");
		}
		*/


		transcripts[genes[n]]++;
		pred[n]->geneno=genes[n]+geneno+1;
		uint t_id=0;
		if (pred[n]->t_eq && pred[n]->t_eq->uptr) {
			t_id = ((RC_TData*)pred[n]->t_eq->uptr)->t_id;
		}

		pred[n]->tlen=abs(pred[n]->tlen);

		//fprintf(f_out,"%d %d %d %.6f %.6f\n",pred[n]->exons.Count()+1,pred[n]->tlen, t_id, pred[n]->frag,pred[n]->cov);
		fprintf(f_out,"1 %d %d %d %.6f\n",pred[n]->exons.Count()+1,pred[n]->tlen, t_id,pred[n]->cov);

		GStr geneid(label, 10);geneid+='.';geneid+=pred[n]->geneno;
		GStr trid(geneid.chars(), 10); trid+='.';trid+=transcripts[genes[n]];
		if(eonly && pred[n]->t_eq && pred[n]->t_eq->getGeneID()) geneid=pred[n]->t_eq->getGeneID();
		if(eonly && pred[n]->t_eq) trid=pred[n]->t_eq->getID();

		fprintf(f_out,"%s\tStringTie\ttranscript\t%d\t%d\t1000\t%c\t.\tgene_id \"%s\"; transcript_id \"%s\"; ",
				refname.chars(),pred[n]->start,pred[n]->end,pred[n]->strand,geneid.chars(),trid.chars());

		//fprintf(stderr,"print pred[%d] gene_id %s transcript_id %s\n",n,geneid.chars(),trid.chars());

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

		// now deal with the genes
		// predicted:
		int gno=genes[n];
		//fprintf(stderr,"gno=%d\n",gno);
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
					bundleData->sum_cov+=pred[n]->cov; // I am adding to cov_sum the reference transcript's coverage --> isn't this double counting?
				}
			}
			else {
				//fprintf(stderr,"2:add pred[%d]'s coverage=%g\n",n,pred[n]->cov);
				bundleData->sum_cov+=pred[n]->cov;
			}
		}
	}

	// clean-up maxint
	while(maxint) {
		nextmaxint=maxint->next;
		delete maxint;
		maxint=nextmaxint;
	}


	geneno+=transcripts.Count();

	return(geneno);
}


void add_pred(GList<CPrediction>& pred,int x,int y, float cov) { // add single exon prediction y to prediction x

	if(pred[y]->start<pred[x]->exons[0].end) { // add y to first exon in x
		int addlen=0;
		if(pred[y]->end<pred[x]->start)  // predictions do not overlap
			addlen=pred[x]->start-pred[y]->end-1;
		pred[x]->cov=(pred[x]->cov*abs(pred[x]->tlen)+cov*abs(pred[y]->tlen))/(abs(pred[x]->tlen)+addlen+abs(pred[y]->tlen));
		pred[x]->exoncov[0]= (pred[x]->exoncov[0]*(pred[x]->exons[0].end-pred[x]->exons[0].start+1)+
						cov*abs(pred[y]->tlen))/(pred[x]->exons[0].end-pred[x]->exons[0].start+1+addlen+abs(pred[y]->tlen));
		if(pred[y]->start<pred[x]->start) {
			pred[x]->start=pred[y]->start;
			pred[x]->exons[0].start=pred[y]->start;
		}
		if(pred[x]->tlen<0) pred[x]->tlen-=addlen;
		else pred[x]->tlen+=addlen;
	}
	else { // add y to last exon in x
		int addlen=0;
		if(pred[x]->end<pred[y]->start) // predictions do not overlap
			addlen=pred[y]->start-pred[x]->end-1;
		pred[x]->cov=(pred[x]->cov*abs(pred[x]->tlen)+cov*abs(pred[y]->tlen))/(abs(pred[x]->tlen)+addlen+abs(pred[y]->tlen));
		pred[x]->exoncov.Last()= (pred[x]->exoncov.Last()*(pred[x]->exons.Last().end-pred[x]->exons.Last().start+1)+
						cov*abs(pred[y]->tlen))/(pred[x]->exons.Last().end-pred[x]->exons.Last().start+1+addlen+abs(pred[y]->tlen));
		if(pred[x]->end<pred[y]->end) {
			pred[x]->end=pred[y]->end;
			pred[x]->exons.Last().end=pred[y]->end;
		}
		if(pred[x]->tlen<0) pred[x]->tlen-=addlen;
		else pred[x]->tlen+=addlen;
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
					if((pred[n]->exons.Count()==1 || pred[m]->exons.Count()==1)) {

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
					GStr trid("", 64);
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


void find_endpoints(int refstart,uint start,uint end,GVec<float>* bpcov, uint &newstart, uint &newend) {


	newstart=start+longintronanchor-1-(uint)refstart;
	newend=end-longintronanchor+1-(uint)refstart;
	start-=(uint)refstart;
	end-=(uint)refstart;

	float maxcov=get_cov(1,start,newstart,bpcov)/(newstart-start+1)-get_cov(1,newstart+1,newstart+longintronanchor,bpcov)/longintronanchor+
			get_cov(1,newend,end,bpcov)/(end-newend+1)-get_cov(1,newend-longintronanchor,newend-1,bpcov)/longintronanchor;

	uint istart=newstart+1;
	uint iend=newend-1;
	while(istart<iend) {
		float icov=get_cov(1,start,istart,bpcov)/(istart-start+1)-get_cov(1,istart+1,istart+longintronanchor,bpcov)/longintronanchor+
				get_cov(1,iend,end,bpcov)/(end-iend+1)-get_cov(1,iend-longintronanchor,iend-1,bpcov)/longintronanchor;

		if(icov>maxcov) {
			newstart=istart;
			newend=iend;
			maxcov=icov;
		}
		istart++;
		iend--;
	}
	newstart+=(uint)refstart;
	newend+=(uint)refstart;

}

/*
uint find_midhash(int refstart,int start,int end,GVec<float>* bpcov) {

	//if(end-refstart>bpcov[1].Count()) return(0);
	uint midpoint=start+longintronanchor-refstart;
	float mincov=get_cov(1,midpoint,midpoint,bpcov);
	start-=refstart;
	end-=refstart;
	float maxcov=fabs(get_cov(1,start,midpoint,bpcov)/(midpoint-start+1)-get_cov(1,midpoint+1,end,bpcov)/(end-midpoint));

	int i=midpoint+1;
	int endcheck=end-longintronanchor;
	while(i<endcheck) {
		float newcov=get_cov(1,i,i,bpcov);
		if(newcov<mincov) {
			midpoint=(uint)i;
			mincov=newcov;
			maxcov=fabs(get_cov(1,start,midpoint,bpcov)/(midpoint-start+1)-get_cov(1,midpoint+1,end,bpcov)/(end-midpoint));
			if(!newcov) break;
		}
		else if(newcov==mincov) { // check if there is a better place to split; longreads modification
			float icov=fabs(get_cov(1,start,i,bpcov)/(i-start+1)-get_cov(1,i+1,end,bpcov)/(end-i));
			if(icov>maxcov) {
				midpoint=(uint)i;
				maxcov=icov;
			}
		}
		i++;
	}
	midpoint+=refstart;

	return(midpoint);
}

uint find_midhash(int refstart,int start,int end,GVec<float>* bpcov) {

	if(end-refstart>bpcov[1].Count()) return(0);
	uint midpoint=start+longintronanchor;
	float mincov=get_cov(1,midpoint-refstart,midpoint-refstart,bpcov);
	end-=longintronanchor;
	int i=midpoint+1;
	while(i<end) {
		float newcov=get_cov(1,i-refstart,i-refstart,bpcov);
		if(newcov<mincov) {
			midpoint=(uint)i;
			mincov=newcov;
			if(!newcov) break;
		}
		i++;
	}

	return(midpoint);
}
*/


int printResults(BundleData* bundleData, int geneno, GStr& refname) {

	uint runoffdist=200;
	if(longreads) runoffdist=0;
	if(bundledist>runoffdist) runoffdist=bundledist;

	// print transcripts including the necessary isoform fraction cleanings
	GList<CPrediction>& pred = bundleData->pred;

	/*
	{ // DEBUG ONLY
		fprintf(stderr,"Pred set before sorting:\n");
		for(int i=0;i<pred.Count();i++) {
			if(pred[i]->t_eq) fprintf(stderr,"%s ",pred[i]->t_eq->getID());
			fprintf(stderr,"pred[%d]:%d-%d (cov=%f, readcov=%f, strand=%c):",i,pred[i]->start,pred[i]->end,pred[i]->cov,pred[i]->tlen*pred[i]->cov,pred[i]->strand);
			for(int j=0;j<pred[i]->exons.Count();j++) fprintf(stderr," %d-%d",pred[i]->exons[j].start,pred[i]->exons[j].end);
			fprintf(stderr,"\n");
		}
		fprintf(stderr,"\n");
	}
	*/

	int npred=pred.Count();

	pred.setSorted(predCmp);

	GVec<CGene> predgene;
	int refstart=bundleData->start;
	GVec<float>* bpcov = bundleData->bpcov;
	int startgno=geneno+1;

	if(rawreads) { //-R mode
		// assign gene numbers
		GVec<int> genes; // for each prediction remembers it's geneno
		genes.Resize(npred,-1);
		GVec<int> transcripts; // for each gene remembers how many transcripts were printed
		GVec<int> color;
		for(int i=0;i<npred;i++) color.Add(i);

		for(int i=0;i<npred-1;i++) if(pred[i]->flag) {
			if(pred[i]->cov<readthr) {
				pred[i]->flag=false;
				//fprintf(stderr,"falseflag: elim pred[%d] due to low cov=%f\n",i,pred[i]->cov);
				continue;
			}
			//fprintf(stderr,"check pred i=%d with end=%d and next start=%d\n",i,pred[i]->end,pred[i+1]->start);
			int ci=color[i];
			while(ci!=color[ci]) { ci=color[ci];color[i]=ci;}
			int j=i+1;
			while(j<npred && pred[i]->end>=pred[j]->start) {
				//fprintf(stderr,"... check pred j=%d\n",j);
				if(pred[j]->cov<readthr) {
					pred[j]->flag=false;
					//fprintf(stderr,"falseflag: elim pred[%d] due to low cov=%f\n",j,pred[j]->cov);
					j++;
					continue;
				}
				if(pred[j]->flag && pred[i]->strand==pred[j]->strand) {
					// check overlap betweeg genes
					bool overlap=false;
					int m=0;
					for(int k=0; k<pred[i]->exons.Count();k++) {
						while(m<pred[j]->exons.Count() && pred[j]->exons[m].end<pred[i]->exons[k].start) m++;
						if(m==pred[j]->exons.Count()) break;
						if(pred[j]->exons[m].overlap(pred[i]->exons[k])) {
							overlap=true;
							break;
						}
						m++;
					}
					if(overlap) {
						//fprintf(stderr,"pred %d overlaps pred %d\n",j,i);
						int cj=color[j];
						while(cj!=color[cj]) { cj=color[cj];color[j]=cj;}
						if(cj<ci) {
							color[ci]=cj;
							color[i]=cj;
						}
						else if(ci<cj){
							color[cj]=ci;
							color[j]=ci;
						}
					}
				}
				j++;
			}
		}

		for(int n=0;n<npred;n++) if(pred[n]->flag && genes[n]==-1) { // no gene was assigned to this prediction
			int cn=color[n];
			genes[n]=transcripts.Count();
			transcripts.cAdd(0);
			while(cn!=color[cn]) { cn=color[cn];color[n]=cn;}
			uint currentend=pred[n]->end;
			int m=n+1;
			while(m<npred && pred[m]->start<=currentend) {
				if(pred[m]->flag) {
					int cm=color[m];
					while(cm!=color[cm]) { cm=color[cm];color[m]=cm;}
					if(cn==cm) {
						if(genes[m]>-1) GError("Overlapping predictions with different gene number!\n");
						genes[m]=genes[n];
					}
					if(pred[m]->end>currentend) currentend=pred[m]->end;
				}
				m++;
			}
		}

		for(int n=0;n<npred;n++) if(pred[n]->flag){

			/*
			{ // DEBUG ONLY
				fprintf(stderr,"print prediction %d with cov=%f len=%d",n,pred[n]->cov,pred[n]->tlen);
				if(pred[n]->flag) fprintf(stderr," with true flag");
				fprintf(stderr," with geneno=%d and exons:",pred[n]->geneno);
				for(int i=0;i<pred[n]->exons.Count();i++) fprintf(stderr," len=%d",pred[n]->exons[i].len());
				fprintf(stderr,"\n");
			}
			*/

			transcripts[genes[n]]++;
			pred[n]->geneno=genes[n]+geneno+1;

			//fprintf(f_out,"%d %d %d %.6f %.6f\n",pred[n]->exons.Count()+1,pred[n]->tlen, t_id, pred[n]->frag,pred[n]->cov);
			fprintf(f_out,"1 %d %d 0 %.6f\n",pred[n]->exons.Count()+1,pred[n]->tlen,pred[n]->cov);

			GStr geneid(label);geneid+='.';geneid+=pred[n]->geneno;
			GStr trid(geneid.chars()); trid+='.';trid+=transcripts[genes[n]];

			fprintf(f_out,"%s\tStringTie\ttranscript\t%d\t%d\t1000\t%c\t.\tgene_id \"%s\"; transcript_id \"%s\"; ",
					refname.chars(),pred[n]->start,pred[n]->end,pred[n]->strand,geneid.chars(),trid.chars());

			//fprintf(stderr,"print pred[%d] gene_id %s transcript_id %s\n",n,geneid.chars(),trid.chars());

			fprintf(f_out,"cov \"%.6f\";\n",pred[n]->cov);
			for(int j=0;j<pred[n]->exons.Count();j++) {
				fprintf(f_out,"%s\tStringTie\texon\t%d\t%d\t1000\t%c\t.\tgene_id \"%s\"; transcript_id \"%s\"; exon_number \"%d\"; ",
						refname.chars(),pred[n]->exons[j].start,pred[n]->exons[j].end,pred[n]->strand,geneid.chars(),
						trid.chars(),j+1); // maybe add exon coverage here
				fprintf(f_out,"\n");
			}

			// now deal with the genes
			// predicted:
			int gno=genes[n];
			//fprintf(stderr,"gno=%d\n",gno);
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

		geneno+=transcripts.Count();

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

				// cummulative bpcov
				float localcov[3]={0,0,0};
				for(int i=0;i<3;i++) {
					if(end) localcov[i]=get_cov(i,start,end-1,bpcov);
				}

				switch(s) {
				case 0:
					if(localcov[2]) cov+=localcov[0]+(localcov[1]-localcov[0]-localcov[2])*localcov[0]/(localcov[0]+localcov[2]);
					else cov+=localcov[1];
					break;
				case 1: cov+=localcov[1]-localcov[2]-localcov[0];break;
				case 2:
					if(localcov[0]) cov+=localcov[2]+(localcov[1]-localcov[0]-localcov[2])*localcov[2]/(localcov[0]+localcov[2]);
					else cov+=localcov[1];
					break;
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

		return(geneno);
	} // ^^^^ rawreads processing mode (-R)

	/*
	{ // DEBUG ONLY
		fprintf(stderr,"Initial set:\n");
		for(int i=0;i<pred.Count();i++) {
			if(pred[i]->t_eq) fprintf(stderr,"%s ",pred[i]->t_eq->getID());
			fprintf(stderr,"pred[%d]:%d-%d (cov=%f, readcov=%f, strand=%c):",i,pred[i]->start,pred[i]->end,pred[i]->cov,pred[i]->tlen*pred[i]->cov,pred[i]->strand);
			for(int j=0;j<pred[i]->exons.Count();j++) fprintf(stderr," %d-%d",pred[i]->exons[j].start,pred[i]->exons[j].end);
			fprintf(stderr," (");
			for(int j=0;j<pred[i]->exons.Count();j++) fprintf(stderr," %f",pred[i]->exoncov[j]);
			fprintf(stderr,")\n");
		}
		fprintf(stderr,"\n");
	}
	*/

	bool preddel=false;

	// this are needed for gene abundance estimations
	GVec<CGene> refgene;
	GHash<int> hashgene;

	bool isintronguide=false;

	// process predictions that equal the same single exon guide and stich them together
	GPVec<GffObj>& guides = bundleData->keepguides;
	if(guides.Count()) {

		//fprintf(stderr,"guides count=%d\n",guides.Count());

		// first create reference genes
		int gno=0;
		for(int i=0;i<guides.Count();i++) {

			if(guides[i]->exons.Count()>1 && longreads) isintronguide=true;

			if (eonly) { // if eonly I need to print all guides that were not printed yet
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
					hashgene.Add(gid.chars(),gno);
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


		// this version is more inclusive by stiching together single exons to reference guides that overlap them, but doesn't print them -> this is still done later on
		GVec< GVec<int> > reflink(npred); reflink.Resize(npred);
		for(int n=0;n<npred;n++) {

			if(pred[n] && pred[n]->t_eq) {

				if(mixedMode && pred[n]->cov<DROP) { // need to be more strict with mixed data since we introduced the guides by default
					pred[n]->flag=false;
					continue;
				}

				//fprintf(stderr,"pred[%d]: start=%d end=%d ID=%s strand=%c refstrand=%c refstart=%d\n",n,pred[n]->start,pred[n]->end,pred[n]->t_eq->getID(),pred[n]->strand,pred[n]->t_eq->strand,pred[n]->t_eq->start);

				// check if there are single exon on the left of the predicted transcript
				int i=n-1;
				while(i>=0 && pred[i]->start>=pred[n]->t_eq->start) {
				  //fprintf(stderr,"%d: excnot=%d end=%d start=%d str=%c\n",i,pred[i]->exons.Count(), pred[i]->end,pred[n]->start, pred[i]->strand);
					if(pred[i]->exons.Count()==1 && pred[i]->end<pred[n]->start && (pred[i]->strand==pred[n]->strand || pred[i]->strand=='.')) {
						if(!pred[i]->t_eq || !strcmp(pred[i]->t_eq->getID(),pred[n]->t_eq->getID())) reflink[i].Add(n);
					}
					if(!pred[i]->t_eq && pred[i]->end>pred[n]->start && pred[i]->cov<pred[n]->cov*ERROR_PERC) {
						//fprintf(stderr,"1 set pred[%d] to false\n",i);
						pred[i]->flag=false;
					}
					i--;
				}
				while(i>=0 && pred[i]->end>pred[n]->start){
					if(!pred[i]->t_eq && pred[i]->cov<pred[n]->cov*ERROR_PERC)
						pred[i]->flag=false;
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
					if(!pred[i]->t_eq && pred[i]->start<pred[n]->end && pred[i]->cov<pred[n]->cov*ERROR_PERC) {
						//fprintf(stderr,"2 set pred[%d] to false\n",i);
						pred[i]->flag=false;
					}
					i++;
				}
				while(i<npred && pred[i]->start<pred[n]->end) {
					if(!pred[i]->t_eq && pred[i]->cov<pred[n]->cov*ERROR_PERC)
						pred[i]->flag=false;
					i++;
			}
		}
		}

		for(int n=0;n<npred;n++) {

		  //fprintf(stderr,"pred[%d] falseflag=%d reflink_cnt=%d\n",n,pred[n]->flag,reflink[n].Count());
			if(pred[n]) {
				if(reflink[n].Count()) {
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
						//fprintf(stderr,"Replace pred[%d] with null\n",n);
					CPrediction *p=pred[n];
					pred.Forget(n);
					delete p;
				}
			}
				/*if(!del && !pred[n]->flag) {
					CPrediction *p=pred[n];
					pred.Forget(n);
					delete p;
				}*/
			}
		}
		pred.Pack();
		pred.Sort();
		npred=pred.Count();

	}

	//preddel=false;
	/* adaptive mode: stitch together nearby single predictions if not within trimming parameters */
	//if(adaptive) { // only in adaptive mode I am storing all single transcripts
	GHash<uint> starthash;
	GHash<uint> endhash;

	bool incomplete=false;
	GStr id("", 32);
	for(int n=0;n<npred-1;n++) {
		//fprintf(stderr,"check pred[%d]:%d-%d:%c with noexon=%d and cov=%f\n",n,pred[n]->start,pred[n]->end,pred[n]->strand,pred[n]->exons.Count(),pred[n]->cov);
		bool ndel=false;
		if(!pred[n]->t_eq && (abs(pred[n]->tlen)<mintranscriptlen)) { ndel=true;}
		else {
			int m=n+1;
			while(!ndel && m<npred && pred[m]->start<=pred[n]->end) {
				if(pred[n]->strand==pred[m]->strand) {
					if(equal_pred(pred,n,m)) {

						if(pred[n]->t_eq && pred[m]->t_eq && pred[n]->t_eq!=pred[m]->t_eq) { m++; continue;} // both are equal but represent different transcripts

						if(mixedMode && pred[n]->tlen*pred[m]->tlen<0) { // choose the larger one
							if(pred[n]->t_eq) {
								pred[m]->t_eq=pred[n]->t_eq;
								pred[m]->flag=true; // make sure I do not delete the prediction
								pred[n]->cov+=pred[m]->cov; // add coverage to annotated gene in this case
							}
							if(pred[n]->cov>pred[m]->cov) { // replace with higher coverage
								pred[m]->start=pred[n]->start;
								pred[m]->end=pred[n]->end;
								pred[m]->cov=pred[n]->cov;
								pred[m]->tlen=-abs(pred[n]->tlen); // here is to say structure is supported by long read --> might decrease precision
								pred[m]->exons[0].start=pred[n]->exons[0].start;
								pred[m]->exons.Last().end=pred[n]->exons.Last().end;
								for(int k=0;k<pred[m]->exons.Count();k++) {
									pred[m]->exoncov[k]=pred[n]->exoncov[k];
								}
								pred[m]->geneno=-abs(pred[m]->geneno);
							}
							//pred[m]->tlen=abs(pred[m]->tlen); // make sure the prediction is conserved as it is
							ndel=true;
							break;
						}

						if(pred[n]->exons.Count()==1) {
							if(pred[m]->cov>ERROR_PERC*pred[n]->cov && pred[n]->cov>ERROR_PERC*pred[m]->cov) { // predictions close to each other in abundance -> otherwise I just delete the less abundant overlapping single exon gene

								if(!pred[m]->t_eq && pred[n]->end>pred[m]->end) pred[m]->end=pred[n]->end;


								//fprintf(stderr,"pred[%d] start=%d end=%d\n",m,pred[m]->start,pred[m]->end);

								if(!pred[m]->t_eq) pred[m]->start=pred[n]->start; // pred[n]->start comes before pred[m]->start -> adjust if pred[m] is not reference
								if(pred[n]->t_eq) {
									pred[m]->t_eq=pred[n]->t_eq;
									pred[m]->end=pred[n]->end;
									pred[m]->flag=true;
								}

								//fprintf(stderr,"pred[%d] start=%d end=%d\n",m,pred[m]->start,pred[m]->end);

								pred[m]->cov=(pred[n]->cov*abs(pred[n]->tlen)+pred[m]->cov*abs(pred[m]->tlen))/(pred[m]->end-pred[m]->start+1);
								pred[m]->tlen=pred[m]->end-pred[m]->start+1;
								if(longreads) pred[m]->tlen=-pred[m]->tlen;
								else if(mixedMode){
									if(pred[m]->cov>pred[n]->cov) {
										if(pred[m]->tlen<0) pred[m]->tlen=-pred[m]->tlen;
									}
									else if(pred[n]->tlen<0) pred[m]->tlen=-pred[m]->tlen;
								}
								pred[m]->exoncov[0]=pred[m]->cov;
								pred[m]->exons[0].start=pred[m]->start;
								pred[m]->exons[0].end = pred[m]->end;

								if(pred[m]->strand=='.') pred[m]->strand=pred[n]->strand;
							}
							else if(pred[n]->t_eq) { m++; continue;}
							else { // prefer prediction with higher coverage
								if(!pred[m]->t_eq && (pred[n]->t_eq || pred[n]->cov>pred[m]->cov)) { // prefer prediction with higher coverage here
									pred[m]->start=pred[n]->start;
									pred[m]->exons[0].start=pred[n]->exons[0].start;
									pred[m]->end=pred[n]->end;
									pred[m]->exons.Last().end=pred[n]->exons.Last().end;
									pred[m]->tlen=pred[n]->tlen;
									pred[m]->t_eq=pred[n]->t_eq;
									if(pred[n]->t_eq) pred[m]->flag=true;
								}
								pred[m]->cov+=pred[n]->cov;
								//fprintf(stderr,"--ndel Prediction m=%d (new cov=%f) is equal\n",m,pred[m]->cov);
								for(int k=0;k<pred[m]->exons.Count();k++) {
									pred[m]->exoncov[k]+=pred[n]->exoncov[k];
								}
							}
							//fprintf(stderr,"--ndel pred[%d] start=%d end=%d\n",m,pred[m]->start,pred[m]->end);
						}
						else {
							if(!pred[m]->t_eq && (pred[n]->t_eq || pred[n]->cov>pred[m]->cov)) { // prefer prediction with higher coverage here
								pred[m]->start=pred[n]->start;
								pred[m]->exons[0].start=pred[n]->exons[0].start;
								pred[m]->end=pred[n]->end;
								pred[m]->exons.Last().end=pred[n]->exons.Last().end;
								pred[m]->tlen=pred[n]->tlen;
								pred[m]->t_eq=pred[n]->t_eq;
								if(pred[n]->t_eq) pred[m]->flag=true;
							}
							pred[m]->cov+=pred[n]->cov;
							//fprintf(stderr,"--ndel Prediction m=%d (new cov=%f) is equal\n",m,pred[m]->cov);
							for(int k=0;k<pred[m]->exons.Count();k++) {
								pred[m]->exoncov[k]+=pred[n]->exoncov[k];
							}
						}
						ndel=true;
						break;
					}
				}
				else if(pred[n]->exons.Count()>1) { // check if overlap on different strands
					// first check for split read alignment errors
					if(m==n+1 && !pred[n]->t_eq && pred[n]->exons.Count()==2) {
						int p=n-1;
						bool err=true;
						while(p>=0) {
							if(pred[p]!=NULL && pred[p]->end>pred[n]->start) {
								err=false;
								break;
							}
							p--;
						}
						if(err) {
							/* version before longreads
							float firstexoncov=get_cov(1,pred[n]->start-bundleData->start,pred[n]->exons[0].end-bundleData->start,bundleData->bpcov)/pred[n]->exons[0].len();
							float lastexoncov=get_cov(1,pred[n]->exons.Last().start-bundleData->start,pred[n]->end-bundleData->start,bundleData->bpcov)/pred[n]->exons.Last().len();*/
							int len1=pred[n]->exons[0].len();
							int len2=pred[n]->exons.Last().len();
							if(len1>CHI_WIN) len1=CHI_WIN;
							if(len2>CHI_WIN) len2=CHI_WIN;
							float firstexoncov=get_cov(1,pred[n]->exons[0].end-len1+1-bundleData->start,pred[n]->exons[0].end-bundleData->start,bundleData->bpcov)/len1;
							float lastexoncov=get_cov(1,pred[n]->exons.Last().start-bundleData->start,pred[n]->exons.Last().start+len2-1-bundleData->start,bundleData->bpcov)/len2;
							if(firstexoncov<lastexoncov*ERROR_PERC) {
								//fprintf(stderr,"--ndel pred[%d] has firstexoncov:%d-%d=%f and lastexoncov:%d-%d=%f\n",n,pred[n]->start,pred[n]->exons[0].end,firstexoncov,pred[n]->exons.Last().start,pred[n]->end,lastexoncov);
								ndel=true;
							}
						}
					}

					if(!ndel && !longreads && pred[m]->exons.Count()>1 && (!pred[n]->t_eq || !pred[m]->t_eq) ) { // at least one prediction can be modified
						int l=0;
						int f=0;
						if(pred[n]->exons.Last().start<pred[m]->exons[0].end) {
							l=n;f=m;
						}
						else if(pred[m]->exons.Last().start<pred[n]->exons[0].end) {
							l=m;f=n;
						}
						if(f || l) {

							int lastval=(int)pred[l]->exons.Last().start;
							if(pred[l]->t_eq) lastval=pred[l]->end;
							int firstval=(int)pred[f]->exons[0].end;
							if(pred[f]->t_eq) firstval=pred[f]->start;

							if(firstval-lastval>(int)(2*longintronanchor)) { // far apart so I can store hash
								//GStr id(lastval);
								id.assign(lastval);
								id+=':';
								id+=firstval;
								uint startval=0;
								uint endval=0;
								uint *val=starthash[id];
								if(!val) {
									find_endpoints(bundleData->start,(uint)lastval,(uint)firstval,bundleData->bpcov,startval,endval);
									starthash.Add(id.chars(),startval);
									endhash.Add(id.chars(),endval);
								}
								else {
									startval=*val;
									val=endhash[id];
									endval=*val;
								}
								// adjust pred[l] and pred[f] here
								if(startval) {
									int prednlen=abs(pred[l]->tlen)-pred[l]->end+startval;
									int predmlen=abs(pred[f]->tlen)-endval+pred[f]->start;

									//fprintf(stderr,"l=%d f=%d prednlen=%d predmlen=%d startval=%d endval=%d predlend=%d predstart=%d\n",l,f,prednlen,predmlen,startval,endval,pred[l]->end,pred[f]->start);

									if(prednlen>mintranscriptlen && predmlen>mintranscriptlen) {

										if(!pred[l]->t_eq && pred[l]->end>startval) {
											float totalcov=get_cov(1,pred[l]->exons.Last().start-bundleData->start,pred[l]->end-bundleData->start,bundleData->bpcov);
											float ratio=0;
											float exoncov=pred[l]->exons.Last().len()*pred[l]->exoncov.Last();
											if(totalcov) ratio=exoncov/totalcov;
											pred[l]->cov=pred[l]->cov*abs(pred[l]->tlen)-exoncov;
											//pred[l]->tlen-=pred[l]->end-midpoint;
											if(pred[l]->tlen<0) pred[l]->tlen=-prednlen;
											else pred[l]->tlen=prednlen;
											pred[l]->exoncov.Last()=ratio*get_cov(1,pred[l]->exons.Last().start-bundleData->start,startval-bundleData->start,bundleData->bpcov);
											pred[l]->cov=(pred[l]->cov+pred[l]->exoncov.Last())/abs(pred[l]->tlen);
											pred[l]->exons.Last().end=startval;
											pred[l]->end=startval;
											pred[l]->exoncov.Last()/=pred[l]->exons.Last().len();
											preddel=true;
											//if(pred[l]->tlen<mintranscriptlen) ndel=true;
										}
										if(!pred[f]->t_eq && pred[f]->start<endval) { // adjust pred[f] coverage
											float totalcov=get_cov(1,pred[f]->start-bundleData->start,pred[f]->exons[0].end-bundleData->start,bundleData->bpcov);
											float ratio=0;
											float exoncov=pred[f]->exons[0].len()*pred[f]->exoncov[0];
											if(totalcov) ratio=exoncov/totalcov;
											pred[f]->cov=pred[f]->cov*abs(pred[f]->tlen)-exoncov;
											//pred[f]->tlen-=midpoint-pred[f]->start;
											if(pred[f]->tlen<0) pred[f]->tlen=-predmlen;
											else pred[f]->tlen=predmlen;
											pred[f]->exoncov[0]=ratio*get_cov(1,endval-bundleData->start,pred[f]->exons[0].end-bundleData->start,bundleData->bpcov);
											pred[f]->cov=(pred[f]->cov+pred[f]->exoncov[0])/abs(pred[f]->tlen);
											pred[f]->exons[0].start=endval;
											pred[f]->start=endval;
											pred[f]->exoncov[0]/=pred[f]->exons[0].len();
											preddel=true;
										}
										//fprintf(stderr,"midpoint=%d pred[n=%d]->cov=%f pred[m=%d]->cov=%f\n",midpoint,n,pred[l]->cov,m,pred[f]->cov);
									}
								}
							}
						}
					}
				}
				m++;
			}
			// now pred[m]->start > pred[n]->end
			uint runoff=runoffdist;
			if(mixedMode && pred[n]->tlen<0) runoff=0;
			if(m<npred && !ndel && !pred[n]->t_eq) { // pred[n] can be deleted
				if (pred[n]->exons.Count()==1) {
					if(pred[m]->exons.Count()==1 && pred[m]->start<pred[n]->end+runoff &&
							(pred[n]->strand == pred[m]->strand || pred[n]->strand == '.' || pred[m]->strand == '.') &&
							pred[m]->cov>ERROR_PERC*pred[n]->cov && pred[n]->cov>ERROR_PERC*pred[m]->cov && !pred[m]->t_eq) {
						// if both are single predictions and are within bundledist on the same strand and not in reference and are not within error_perc from each other
						//fprintf(stderr,"Stich predictions %d and %d with cov=%f and %f\n",n,m,pred[n]->cov,pred[m]->cov);
						pred[m]->start=pred[n]->start;
						pred[m]->cov=(pred[n]->cov*abs(pred[n]->tlen)+pred[m]->cov*abs(pred[m]->tlen))/(pred[m]->end-pred[m]->start+1);
						int mult=1;
						if(longreads) mult=-1;
						else if(mixedMode){
							if(pred[m]->cov>pred[n]->cov) {
								if(pred[m]->tlen<0) mult=-1;
							}
							else if(pred[n]->tlen<0) mult=-1;
						}
						pred[m]->tlen=mult*(pred[m]->end-pred[m]->start+1);
						pred[m]->exoncov[0]=pred[m]->cov;
						pred[m]->exons[0].start=pred[m]->start;
						if(pred[m]->strand=='.') pred[m]->strand=pred[n]->strand;
						ndel=true;
					}
					else if(pred[n]->cov<singlethr) { /****** single exon different threshold ******/
						ndel=true;
					}
					else if(m==n+1){
						int m=n-1;
						while(m>=0 && pred[m]==NULL) m--;
						//if(m>=0 && pred[m]->end>=pred[n]->start) fprintf(stderr,"pred[%d] overlaps pred[%d]\n",n,m);
						if(m<0 || pred[m]->end<pred[n]->start) { // pred[n] does not overlap previous prediction --> check if I can delete it
							if(pred[n+1]->start-pred[n]->end<runoff) {
								float exoncov=get_cov(1,pred[n+1]->start-bundleData->start,pred[n+1]->exons[0].end-bundleData->start,bundleData->bpcov)/pred[n+1]->exons[0].len();
								//fprintf(stderr,"Pred %d to be deleted because under first exoncov of pred %d\n",n,n+1);
								if(pred[n]->cov<exoncov+singlethr) ndel=true;
							}
							if(!ndel && m>=0 && pred[n]->start-pred[m]->end<runoff) {
								float exoncov=get_cov(1,pred[m]->exons.Last().start-bundleData->start,pred[m]->end-bundleData->start,bundleData->bpcov)/pred[m]->exons.Last().len();
								//fprintf(stderr,"Pred %d to be deleted because under last exoncov of pred %d\n",n,m);
								if(pred[n]->cov<exoncov+singlethr) ndel=true;
							}
						}
					}
				}
				else if(m==n+1 && pred[n]->exons.Count()==2) { // check for splice read alignment error
					int p=n-1;
					bool err=true;
					while(p>=0) {
						if(pred[p]!=NULL && (pred[p]->end>=pred[n]->exons.Last().start || (pred[p]->strand==pred[n]->strand && pred[p]->end>=pred[n]->start))) {
							err=false;
							break;
						}
						p--;
					}
					if(err) {
						/* version before longreads
						float firstexoncov=get_cov(1,pred[n]->start-bundleData->start,pred[n]->exons[0].end-bundleData->start,bundleData->bpcov)/pred[n]->exons[0].len();
						float lastexoncov=get_cov(1,pred[n]->exons.Last().start-bundleData->start,pred[n]->end-bundleData->start,bundleData->bpcov)/pred[n]->exons.Last().len();*/
						int len1=pred[n]->exons[0].len();
						int len2=pred[n]->exons.Last().len();
						if(len1>CHI_WIN) len1=CHI_WIN;
						if(len2>CHI_WIN) len2=CHI_WIN;
						float firstexoncov=get_cov(1,pred[n]->exons[0].end-len1+1-bundleData->start,pred[n]->exons[0].end-bundleData->start,bundleData->bpcov)/len1;
						float lastexoncov=get_cov(1,pred[n]->exons.Last().start-bundleData->start,pred[n]->exons.Last().start+len2-1-bundleData->start,bundleData->bpcov)/len2;
						if(lastexoncov<firstexoncov*ERROR_PERC) {
							ndel=true;
							//fprintf(stderr,"pred[%d] has firstexoncov:%d-%d=%f and lastexoncov:%d-%d=%f\n",n,pred[n]->start,pred[n]->exons[0].end,firstexoncov,pred[n]->exons.Last().start,pred[n]->end,lastexoncov);
						}
					}
				}
			}
		}
		if(ndel) {
			// now replace pred[n] with null
			CPrediction *p=pred[n];
			pred.Forget(n);
			delete p;
			preddel=true;
		}
		else if(isintronguide && !incomplete && !pred[n]->t_eq && pred[n]->exons.Count()>1 && pred[n]->mergename.is_empty()){
			incomplete=true;
		}
	}

	// stich last prediction
	if(npred && !pred[npred-1]->t_eq) {
		bool check=true;
		if(abs(pred[npred-1]->tlen)<mintranscriptlen || (pred[npred-1]->exons.Count()==1 && pred[npred-1]->cov<singlethr)) { 	/****** single exon different threshold ******/
			// now replace pred[n] with null
			CPrediction *p=pred[npred-1];
			pred.Forget(npred-1);
			delete p;
			preddel=true;
			check=false;
		}
		else if(pred[npred-1]->exons.Count()==1) {
			if(npred>1) {
				int m=npred-2;
				while(m>=0 && pred[m]==NULL) m--;
				uint runoff=runoffdist;
				if(mixedMode && pred[npred-1]->tlen<0) runoff=0;
				if(m>=0 && pred[m]->end<pred[npred-1]->start && pred[npred-1]->start-pred[m]->end<runoff) {
					float exoncov=get_cov(1,pred[m]->exons.Last().start-bundleData->start,pred[m]->end-bundleData->start,bundleData->bpcov)/pred[m]->exons.Last().len();
					if(pred[npred-1]->cov<exoncov) {
						CPrediction *p=pred[npred-1];
						pred.Forget(npred-1);
						delete p;
						preddel=true;
						check=false;
					}
				}
			}
		}
		else if(pred[npred-1]->exons.Count()==2) { // check for spliced read alignment error
			int k=npred-2;
			bool err=true;
			while(k>=0) {
				if(pred[k]!=NULL && (pred[k]->end>pred[npred-1]->exons.Last().start || (pred[k]->strand==pred[npred-1]->strand && pred[k]->end>=pred[npred-1]->start))) {
					err=false;
					break;
				}
				k--;
			}
			if(err) {
				/* version before longreads
				float firstexoncov=get_cov(1,pred[npred-1]->start-bundleData->start,pred[npred-1]->exons[0].end-bundleData->start,bundleData->bpcov)/pred[npred-1]->exons[0].len();
				float lastexoncov=get_cov(1,pred[npred-1]->exons.Last().start-bundleData->start,pred[npred-1]->end-bundleData->start,bundleData->bpcov)/pred[npred-1]->exons.Last().len();*/
				int len1=pred[npred-1]->exons[0].len();
				int len2=pred[npred-1]->exons.Last().len();
				if(len1>CHI_WIN) len1=CHI_WIN;
				if(len2>CHI_WIN) len2=CHI_WIN;
				float firstexoncov=get_cov(1,pred[npred-1]->exons[0].end-len1+1-bundleData->start,pred[npred-1]->exons[0].end-bundleData->start,bundleData->bpcov)/len1;
				float lastexoncov=get_cov(1,pred[npred-1]->exons.Last().start-bundleData->start,pred[npred-1]->exons.Last().start+len2-1-bundleData->start,bundleData->bpcov)/len2;
				if(lastexoncov<firstexoncov*ERROR_PERC) {
					CPrediction *p=pred[npred-1];
					pred.Forget(npred-1);
					delete p;
					check=false;
					preddel=true;
				}
			}
		}
		if(check && isintronguide && !incomplete && !pred[npred-1]->t_eq && pred[npred-1]->exons.Count()>1 && pred[npred-1]->mergename.is_empty()){
			incomplete=true;
		}
	}
	if(preddel) {
		pred.Pack();
		pred.Sort();
		npred=pred.Count();
	}
	//}

	/*
    { // DEBUG ONLY
    	fprintf(stderr,"Before predcluster:\n");
		for(int i=0;i<pred.Count();i++) {
    		if(pred[i]->t_eq) fprintf(stderr,"%s ",pred[i]->t_eq->getID());
    		fprintf(stderr,"pred[%d]:%d-%d (cov=%f, readcov=%f, strand=%c falseflag=%d):",i,pred[i]->start,pred[i]->end,pred[i]->cov,pred[i]->tlen*pred[i]->cov,pred[i]->strand,pred[i]->flag);
    		//for(int j=0;j<pred[i]->exons.Count();j++) fprintf(stderr," %d-%d(%f)",pred[i]->exons[j].start,pred[i]->exons[j].end,pred[i]->exoncov[j]);
			for(int j=0;j<pred[i]->exons.Count();j++) fprintf(stderr," %d-%d",pred[i]->exons[j].start,pred[i]->exons[j].end);
			fprintf(stderr," (");
			for(int j=0;j<pred[i]->exons.Count();j++) fprintf(stderr," %f",pred[i]->exoncov[j]);
			fprintf(stderr,")\n");
    	}
		fprintf(stderr,"\n");
    }
    */

	if(npred) geneno=print_predcluster(pred,geneno,refname,refgene,hashgene,predgene,bundleData,incomplete);

	hashgene.Clear();
	// I am done printing all transcripts, now evaluate/print the gene abundances
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

					if(start<bpcov[1].Count()) {
						if(start<0) start=0;
						if(end>=bpcov[1].Count()) end=bpcov[1].Count()-1;

						// cummulative bpcov
						float localcov[3]={0,0,0};
						for(int i=0;i<3;i++) {
							if(end) localcov[i]=get_cov(i,start,end-1,bpcov);
						}

						switch(s) {
						case 0:
							if(localcov[2]) cov+=localcov[0]+(localcov[1]-localcov[0]-localcov[2])*localcov[0]/(localcov[0]+localcov[2]);
							else cov+=localcov[1];
							break;
						case 1: cov+=localcov[1]-localcov[2]-localcov[0];break;
						case 2:
							if(localcov[0]) cov+=localcov[2]+(localcov[1]-localcov[0]-localcov[2])*localcov[2]/(localcov[0]+localcov[2]);
							else cov+=localcov[1];
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

				// cummulative bpcov
				float localcov[3]={0,0,0};
				for(int i=0;i<3;i++) {
					if(end) localcov[i]=get_cov(i,start,end-1,bpcov);
				}

				switch(s) {
				case 0:
					if(localcov[2]) cov+=localcov[0]+(localcov[1]-localcov[0]-localcov[2])*localcov[0]/(localcov[0]+localcov[2]);
					else cov+=localcov[1];
					break;
				case 1: cov+=localcov[1]-localcov[2]-localcov[0];break;
				case 2:
					if(localcov[0]) cov+=localcov[2]+(localcov[1]-localcov[0]-localcov[2])*localcov[2]/(localcov[0]+localcov[2]);
					else cov+=localcov[1];
					break;
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


	/*** I might resque this if I want to write a non-conflicting (among threads) version of the covered references
	if (c_out) {
		for (int i=0;i<bundleData->covguides.Count();i++)
			bundleData->covguides[i]->print(c_out);
	}
	***/


	//rc_write_counts(refname.chars(), *bundleData);
	return(geneno);
}

