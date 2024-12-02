/*
 * print_results.cpp
 *
 *  Created on: Oct 2, 2023
 *      Author: mpertea
 */


#include "rlink.h"


extern bool longreads;
extern uint bundledist;  // reads at what distance should be considered part of separate bundles
                        // <- this is not addressed everywhere, e.g. in infer_transcripts -> look into this
extern bool geneabundance; // need to compute the gene abundance
extern bool rawreads;
extern float readthr;     // read coverage per bundle bp to accept it; otherwise considered noise
extern FILE* f_out;
extern GStr label;
extern bool eonly;
extern bool mixedMode;
extern bool guided;
extern int mintranscriptlen; // minimum number for a transcript to be printed
extern float singlethr;     // read coverage per bundle bp to accept it; otherwise considered noise
extern float isofrac;
extern uint junctionsupport; // anchor length for junction to be considered well supported <- consider shorter??
extern bool includecov;
extern bool enableNames;

int predordCmp(const pointer p1, const pointer p2) { // sort from highest to lowest coverage
	CPred *a=(CPred*)p1;
	CPred *b=(CPred*)p2;
	if(a->cov < b->cov) return 1;
	if(a->cov > b->cov) return -1;
	return 0;
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

/*
 *

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

int predcovCmp(const pointer p1, const pointer p2) { // sort from highest to lowest coverage
	CPrediction *a=(CPrediction*)p1;
	CPrediction *b=(CPrediction*)p2;
	if(a->cov < b->cov) return 1;
	if(a->cov > b->cov) return -1;
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

 *
 */

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
*/

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

/*void update_cov(GList<CPrediction>& pred,int big,int small,float frac=1) { // small gets included into big

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

}*/


/*
uint min(uint n1,uint n2) {
	if(n1<n2) return(n1);
	return(n2);
}

uint max(uint n1,uint n2) {
	if(n1<n2) return(n2);
	return n1;
}
*/



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
				nextmaxint=new CMaxIntv(maxint->node,start,maxint->end, maxint->cov, maxint->next);
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
				maxint=new CMaxIntv(nextmaxint->node,end+1,nextmaxint->end, nextmaxint->cov, nextmaxint->next);
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
					CMaxIntv *newintv=new CMaxIntv(maxint->node,end+1,maxint->end, maxint->cov, maxint->next);
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
				if (guides[i]->uptr && getGuideStatus(guides[i])<GBST_STORED) {
				                       //((RC_TData*)guides[i]->uptr)->in_bundle<3) {
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

/*
void printGff3Header(FILE* f, GArgs& args) {
  fprintf(f, "# ");
  args.printCmdLine(f);
  fprintf(f, "##gff-version 3\n");
}
*/


