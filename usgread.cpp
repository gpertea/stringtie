//#define GMEMTRACE 1  //debugging memory allocation
#ifdef GMEMTRACE
#include "proc_mem.h"
#endif

#include "usgread.h"

extern uint sserror;
extern int mintranscriptlen; // minimum number for a transcript to be printed
extern bool longreads;
extern uint bundledist;  // reads at what distance should be considered part of separate bundles

void mark_usgstart(GList<CJunction>& junction,int &js,int &jend,uint start,uint end,int val) {
	int je=jend;
	while(js<=je) {
		int mid= js +(je-js)/2;
		uint midval=junction[mid]->start;
		if(midval==start) { // found junction in array
			int j=mid;
			while(j>=0 && junction[j]->start==start) {
				junction[j]->usg_start=val;
				j--;
			}
			j=mid+1;
			while(j<=jend && junction[j]->start==start) {
				junction[j]->usg_start=val;
				j++;
			}
			js=j;
			return;
		}
		else if(midval<start) { //Narrow search to the right half
			js=mid+1;
		}
		else { // Narrow search to the left half
			je=mid-1;
			if(junction[mid]->start>end) jend=mid; // junction at mid starts after the read end
		}
	}
}

// fund junction start returns the node that comes right at or right before the end of the junction; if later then the return will be -2 - u_position
int find_junction_start_in_USG(GList<CJunction>& junction,int &js,int &je,CJunction *rdjunc,uint readend,SGBundle* usgbundle,int u,int nu) {
// u is starting at position >= than read; junction makes sure that I do not search jstart more than once -> once I find it once I can update all junctions nearby

	if(rdjunc->usg_start==-2) return(-2); // this was searched before


	if(rdjunc->start<usgbundle->nodes[u]->position || rdjunc->start>usgbundle->nodes.Last()->position) {
		// mark all junctions starting at this position with -2 so I won't search them again
		mark_usgstart(junction,js,je,rdjunc->start,readend,-2);
		return(-2);
	}

	nu--;

	while(u<=nu) {
		int mid = u + (nu - u) / 2; // Find the middle index
		uint midvalue=usgbundle->nodes[mid]->position; // Value at the middle index
		if(midvalue==rdjunc->start) { // Target found at index mid
			// mark all junctions starting at this position with mid so I won't search them again
			mark_usgstart(junction,js,je,rdjunc->start,readend,mid);

			return(mid);
		}
		else if(midvalue < rdjunc->start)
			u=mid +1; //Narrow search to the right half
		else nu=mid-1; // Narrow search to the left half
	}


	// mark all junctions starting at this position with -2 so I won't search them again
	mark_usgstart(junction,js,je,rdjunc->start,readend,-2);
	return(-2);

}

int find_USG_predecessor(SGBundle* usgbundle,int low,int high,uint val) {
	int predecessor = low;  // Default value if no smaller element is found

	while(low<=high) {
		int mid = low + (high - low) / 2; // Find the middle index
		uint midvalue=usgbundle->nodes[mid]->position; // Value at the middle index
		if(midvalue==val) { // Target found at index mid -> only can happen if junction start is not in USG
			return(mid);
		}
		else if(midvalue < val) {
			predecessor=mid;  // Update predecessor since arr[mid] is smaller than v
			low=mid +1; // Narrow search to the right half
		}
		else high=mid-1; // Narrow search to the left half
	}

	return(-3-predecessor);
}

int find_node_in_array(int n,GArray<int> &result,int low, int high) { // binary search of n node position in result

	while(low<=high) {
		int mid = low + (high - low) / 2; // Find the middle index
		if(result[mid]==n) {
			return(mid);
		}
		else if(result[mid]<n) {
			low=mid+1;
		}
		else high=mid-1;
	}

	return(-1); // did not find the node
}

// finds a USG node of a certain type that is within the required distance; assumes node u->position is <= than val; but u+1 is > val
int find_nearby_type(uint val,SGNodeType type,SGBundle* usgbundle,int u, int nu,uint dist) {

	int near=u;

	int cu=u;
	bool found=false;
	while(cu>=0 && val-usgbundle->nodes[cu]->position<dist) {
		if(usgbundle->nodes[cu]->type==type) {
			near=cu;
			found=true;
			break;
		}
		cu--;
	}

	cu=u+1;
	while(cu<nu && usgbundle->nodes[cu]->position-val<dist) {
		if(usgbundle->nodes[cu]->type==type) {
			if(!found || usgbundle->nodes[cu]->position-val<val-usgbundle->nodes[near]->position) near=cu; // closer than the previous junction
			break;
		}
		cu++;
	}

	return(near);
}


// find junction start returns the node that comes right at or right before the end of the junction;
// assumptions: if junction start is within close proximity (sserror) of a JSTART it will adjust the start at that JSTART
// if junction end is within close proximity of a JEND it will adjust the start at that JEND
// if USG node is not found then the return will be -3 - u_position, where u is node right before
// if start is before u[0]->position then returns -2
// if start is after u.Last->position then returns -3-nu+1
// it also assumes no two type of nodes in the USG boundle can start at the same position
void find_junction_in_USG(GList<CJunction>& junction,int &js,int &je,CReadAln & rd,int j,SGBundle* usgbundle,int u,int nu) {
// u is starting at position >= than read; junction makes sure that I do not search jstart more than once -> once I find it once I can update all junctions nearby

	if(rd.juncs[j]->usg_start==-1) { // only if this start wasn't searched before

		// first search start of junction
		if(rd.juncs[j]->start<usgbundle->nodes[0]->position) { // junction starts before the USG bundle

			int cu=0;
			if(rd.juncs[j]->len()>sserror) // only if junction end is far enough
				while(cu<nu && usgbundle->nodes[cu]->position-rd.juncs[j]->start<sserror) {
					if(usgbundle->nodes[cu]->type==JSTART) {
						rd.juncs[j]->usg_start=-3-cu; // mark predicted start
						break;
					}
					cu++;
				}

			if(rd.juncs[j]->usg_start==-1) rd.juncs[j]->usg_start=-2;

			mark_usgstart(junction,js,je,rd.juncs[j]->start,rd.end,rd.juncs[j]->usg_start);
		}
		else { // junction start is within the bundle or after it
			if(rd.juncs[j]->start>usgbundle->nodes.Last()->position) { // junction after the bundle -> just mark it
				rd.juncs[j]->usg_start=-3-nu+1;
				mark_usgstart(junction,js,je,rd.juncs[j]->start,rd.end,rd.juncs[j]->usg_start);
				rd.juncs[j]->usg_end=-3-nu+1; // the junction end also has to be after the end of bundle
				return; // I do not need to search anymore
			}

			// search for junction start in usg; because junction start is within close proximity of the read start or the previous junction end I do not need to do a binary search
			if(j) { // second juction in read -> adjust u
				if(rd.juncs[j-1]->usg_end>0) u=rd.juncs[j-1]->usg_end;
				else {
					u=-3-rd.juncs[j-1]->usg_end;
				}
			}

			if(u<0) u++; // junction has to start after the USG bundle start according to the first if above

			uint juncstart=rd.juncs[j]->start+1; // TODO: check that USG node always includes the first base in the node

			while(u<nu && usgbundle->nodes[u]->position<=juncstart) u++; // while nodes in sgbundle start before read advance u
			u--;

			if(usgbundle->nodes[u]->position==juncstart) { // make sure that I do not allow USG to have the same position twice; otherwise I need to make sure this USG node is a jstart
				if(usgbundle->nodes[u]->type==JSTART) rd.juncs[j]->usg_start=u;
				else  { // try to find nearby start
					if(rd.juncs[j]->len()>sserror) u=find_nearby_type(juncstart,JSTART,usgbundle,u,nu,sserror);
					rd.juncs[j]->usg_start=-3-u;
				}
			}
			else {
				rd.juncs[j]->usg_start=-3-u;
			}
			mark_usgstart(junction,js,je,rd.juncs[j]->start,rd.end,rd.juncs[j]->usg_start);
		}
	}

	// search for junction end now
	if(rd.juncs[j]->usg_end==-1) { // TODO: check that end doesn't need to be adjusted because it starts at the junction end
		if(rd.juncs[j]->usg_start>0) u=rd.juncs[j]->usg_start;
		else {
			u=-3-rd.juncs[j]->usg_start;
		}
		if(u<0) { // junction start comes before USG bundle start
			if(rd.juncs[j]->end<usgbundle->nodes[0]->position) { // junction starts before the USG bundle
				rd.juncs[j]->usg_end=-2;
				return;
			}
			u++; // make sure I have always a positive u
		}

		if(rd.juncs[j]->usg_start>0) { // the junction might be in the USG already -> no need to search for it
			int k=0;
			while(k<usgbundle->nodes[u]->jxLinks.Count() && usgbundle->nodes[u]->jxLinks[k].pos < rd.juncs[j]->end) k++;
			if(usgbundle->nodes[u]->jxLinks[k].pos == rd.juncs[j]->end) { // matching end -> I assume there is only such end that can match -> no two junctions with different strands ending here

				char s=0; // unknown strand
				if(usgbundle->nodes[u]->jxLinks[j].strand=='+') s=1; // junction on positive strand
				else if(usgbundle->nodes[u]->jxLinks[j].strand=='-') s=-1; // guide on negative strand

				/* do this if I want to mark junction strand not in graph
				if(s==rd.juncs[j]->strand || !rd.juncs[j]->strand) { // matching strand
					rd.juncs[j]->usg_end=usgbundle->nodes[u]->jxLinks[j].jx->bidx;
				}
				else {
					rd.juncs[j]->usg_end=-3-usgbundle->nodes[u]->jxLinks[j].jx->bidx; // junction start might match but junction end might differ if the strand doesn't match
					rd.juncs[j]->usg_start=-3-rd.juncs[j]->usg_start; // mark junction start to signal that I have not found a start
				}*/

				if(s!=rd.juncs[j]->strand && rd.juncs[j]->strand)  // non-matching strand -> unstrand it
					rd.juncs[j]=0;

				rd.juncs[j]->usg_end=usgbundle->nodes[u]->jxLinks[j].jx->bidx;

				return;
			}
		}
		// junction start doesn't match or the end is not in USG -> I need to do a binary search for the junction end

		/*
		int pos=find_USG_predecessor(usgbundle,u,nu-1,rd.juncs[j]->end);
		if(pos>0) pos=-3-pos; // do this if I want to mark junction not in usg
		*/
		u=find_USG_predecessor(usgbundle,u,nu-1,rd.juncs[j]->end);

		if(rd.juncs[j]->end==usgbundle->nodes[u]->position && usgbundle->nodes[u]->type==JEND) {
			rd.juncs[j]->usg_end=u;
		}
		else {
			if(rd.juncs[j]->len()>sserror) u=find_nearby_type(rd.juncs[j]->end,JEND,usgbundle,u,nu,sserror);
			rd.juncs[j]->usg_end=-3-u;
		}
	}
}

int add2child(UCNode *node,int s,UCNode *child) {

	int low=0;
	int high=node->childlink[s].Count()-1;

	while(low<=high) {
		int mid = low + (high - low) / 2; // Find the middle index
		int midvalue=node->childlink[s][mid]->prevSGnode; // Value at the middle index
		if(midvalue==child->prevSGnode) { // Target found at index mid -> only can happen if junction start is not in USG
			return(mid);
		}
		else if(midvalue < child->prevSGnode) {
			low=mid +1; // Narrow search to the right half
		}
		else high=mid-1; // Narrow search to the left half
	}

	node->childlink[s].Insert(low,child);

	return(low);

}

int add2parent(UCNode *node,int s,UCNode *parent) {

	int low=0;
	int high=node->parentlink[s].Count()-1;

	while(low<=high) {
		int mid = low + (high - low) / 2; // Find the middle index
		int midvalue=node->parentlink[s][mid]->prevSGnode; // Value at the middle index
		if(midvalue==parent->prevSGnode) { // Target found at index mid -> only can happen if junction start is not in USG
			return(mid);
		}
		else if(midvalue < parent->prevSGnode) {
			low=mid +1; // Narrow search to the right half
		}
		else high=mid-1; // Narrow search to the left half
	}

	node->parentlink[s].Insert(low,parent);

	return(low);

}


void insert_covlink(GPVec<UCNode> &SG2UCnode,int n,int nx,int s,float count) {
	if(!SG2UCnode[n]) {
		SG2UCnode[n]= new UCNode();
		SG2UCnode[n]->prevSGnode=n;
	}
	if(!SG2UCnode[nx]) {
		SG2UCnode[nx]= new UCNode();
		SG2UCnode[n]->prevSGnode=nx;
	}
	// search if SG2UCnode[nx] was inserted before
	int p=add2child(SG2UCnode[n],s,SG2UCnode[nx]); // and child nx to n
	SG2UCnode[n]->childlink[s][p]=SG2UCnode[nx];
	SG2UCnode[n]->childcov[s][p]+=count;
	add2parent(SG2UCnode[nx],s,SG2UCnode[n]); // add parent n to nx
}

void merge_UCNodecov(GPVec<UCNode> &SG2UCnode,int n,uint end,uint start) {

	// TODO: write this --> merge all coverages in SG2UCnode that end at >=end and start >=start
	// first: extend the coverage that ends at  >=end to have an end at start if not already goes beyond

}

// do intervals that come right after one another need to be merged? -> no if I employ bundledist -> TODO: maybe use a smaller bundledist for USG? and merge stranded ones no matter what?
void add_cov2UCNode(GPVec<UCNode> &SG2UCnode,int n,uint start,uint end,int s,float count) {

	if(!SG2UCnode[n]) { // node is NULL
		SG2UCnode[n]= new UCNode();
		SG2UCnode[n]->prevSGnode=n;
	}

	if(!SG2UCnode[n]->covintv[s]) { // node has no coverage yet
		UCov *cov=new UCov(start,end,count*(end-start+1));
		SG2UCnode[n]->covintv[s]=cov;
		return;
	}

	UCov *cov=SG2UCnode[n]->covintv[s];
	if(end<cov->start) { // insert new coverage before
		UCov *prevcov=new UCov(start,end,count*(end-start+1));
		prevcov->next=cov;
		SG2UCnode[n]->covintv[s]=prevcov;
		return;
	}

	UCov *prevcov=NULL;
	while(cov && start>cov->end) { // skip intervals that new coverage does not overlap
		prevcov=cov;
		cov=cov->next; // skip cov until I get to the one that overlaps my coverage interval or after it
	}

	if(!cov || end<cov->start) { // no overlap: start > prevcov->end end<cov->start or cov is null

		UCov *newcov=new UCov(start,end,count*(end-start+1));
		newcov->next=cov;
		prevcov->next=newcov;
	}
	else { // end>=cov->start --> there is overlap

		if(start<cov->start) {
			cov->start=start;
		}
		cov->cov+=(end-start+1)*count;
		UCov *nextcov=cov->next;

		while(nextcov && end>=nextcov->start) { // merge while the new interval intersect the existing coverages
			cov->end=nextcov->end;
			cov->cov+=nextcov->cov;
			cov->next=nextcov->next;
			delete nextcov;
			nextcov=cov->next;
		}

		if(end>cov->end) cov->end=end;

	}

}


void add_to_UCNode(int s,GList<CReadAln>& readlist,int n,int np,float count,GVec<int> &read2unode,SGBundle* usgbundle,int nu,GPVec<UCNode> &SG2UCnode) {

	CReadAln & rd=*(readlist[n]);
	int u=read2unode[n]; // previous node in USG bundle from read
	if(u<0) u=0; // start unode

	int nex=rd.segs.Count(); // number of exons in read

	int i=0;
	int lu=u; // last USG node spanned by read

	while(i<nex) { // process exon i

		if(rd.segs[i].end<usgbundle->nodes[0]->position) continue; // skip exons before the start of bundle
		if(u==nu-1) return; // I do not need to create links after this node which is the last in the USG (terminal node)

		uint start=rd.segs[i].start; // start of exon
		uint end=rd.segs[i].end; // end of exon

		int su=u; // remember u

		// adjust start to be after the start USG node -> start might be before in case of a junction nearby
		if(start<usgbundle->nodes[su]->position) start=usgbundle->nodes[su]->position;

		if(i+1<nex) { // there is another exon after this one

			// how far would this read exon go; lu=rd.juncs[i]->usg_start; if rd.juncs[i]->usg_start>-1
			if(rd.juncs[i]->usg_start<0) { // correct end of exon is it passes junction start by a bit
				lu=-3-rd.juncs[i]->usg_start;
				if(end-usgbundle->nodes[lu]->position<sserror && usgbundle->nodes[lu]->type==JSTART)
					end=usgbundle->nodes[lu]->position-1;
			}

			// calculate next exon's start USG node
			if(rd.juncs[i]->usg_end>-1) u=rd.juncs[i]->usg_end;
			else u=-3-rd.juncs[i]->usg_end;

		}

		if(start>end) continue; // this is an empty exon -> would this ever happen?

		while(start<=end && su<nu) { // start>=usgbundle->nodes[su]->position
			uint minend = end;
			if(su+1<nu && end>=usgbundle->nodes[su+1]->position) {
				minend=usgbundle->nodes[su+1]->position-1;
				if(!SG2UCnode[su+1]) {
					SG2UCnode[su+1]= new UCNode();
				}
				insert_covlink(SG2UCnode,su,su+1,s,count);
			}
			add_cov2UCNode(SG2UCnode,su,start,end,s,count);
			start=minend+1;
			su++;
		}
		su--;

		lu=su; // last node spanned by exon

		// link su to u
		if(su<u) {
			if(!SG2UCnode[u]) {
				SG2UCnode[u]= new UCNode();
			}
			insert_covlink(SG2UCnode,su,u,s,count);
		}

		i++;
	}

	// add link to pair if present
	if(np>-1) {
		if(lu<read2unode[np]) insert_covlink(SG2UCnode,lu,read2unode[np],s,count);
		else if(rd.end<readlist[np]->start) { // if paired read is within the same USG node -> I need to merge coverages in node u up to this point
			merge_UCNodecov(SG2UCnode,lu,rd.end,readlist[np]->start);
		}
	}

}


void populate_UCNodes(BundleData* bdata, GPVec<UCNode> &SG2UCnode,GVec<int> &read2unode) {
	GList<CReadAln>& readlist = bdata->readlist;
	GList<CJunction>& junction = bdata->junction;
	GVec<float>* bpcov = bdata->bpcov;
	int refstart=bdata->start;
	int refend=bdata->end;
	SGBundle* usgbundle = bdata->usgbundle;

	for(int s=0;s<3;s++) bpcov[s].Resize(refend-refstart+3);

	int nu=usgbundle->nodes.Count(); // number of nodes in usg bundle
	for(int i=0;i<nu;i++) SG2UCnode[i]=NULL; // do I need this?


	// if useUSG for a bundle between position1 to position 2 in USG, processRead only loads reads from the BAM file that overlap the region
	// from position1 to position2. Reads outside this range in the sample are ignored -> means that bdata->usgbundle should always be present if useUSG

	int u=0; // next usg node to process
	int jst=0; // for USG only keeps track of junctions that were seen already

	for(int n=0;n<readlist.Count();n++) {

		while(u<nu && usgbundle->nodes[u]->position<=readlist[n]->start) u++; // while nodes in sgbundle start before read advance u
		u--;
		read2unode[n]=u;

		while(jst<junction.Count() && junction[jst]->start<readlist[n]->start) jst++; // advance jst to next junction that makes sense to search

		int js=jst;                // defines junction search space
		int je=junction.Count()-1; // defines junction search space

		CReadAln & rd=*(readlist[n]);
		float rdcount=rd.read_count;
		int nex=rd.segs.Count();

		if(nex>1) {
			char rdstrand=0;
			bool juncstrand=false; // did not find a junction strand

			for(int i=1;i<nex;i++) {
				// search junction
				find_junction_in_USG(junction,js,je,rd,i-1,usgbundle,u,nu); // u is starting at position >= than read; junction makes sure that I do not search jstart more than once -> once I find it once I can update all junctions nearby

				if(rd.juncs[i-1]->usg_start>-1 && rd.juncs[i-1]->usg_end>-1 && rd.juncs[i-1]->strand) { // found a stranded junction
					if(juncstrand && rdstrand!=rd.juncs[i-1]->strand) { // there is conflict in assigning the read strand
						rd.strand=0;
					}
					else {
						rdstrand=rd.juncs[i-1]->strand;
						rd.strand=rdstrand;
						juncstrand=true;
					}
				}

			}
			if(!juncstrand) rd.strand=0; // delete strand of read if no junction found
		}

		if(!rd.unitig) add_read_to_cov(readlist,n,bpcov,refstart); // unitig special case; do I need this here?
		else if(rdcount>1) rdcount=1;
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


	// strandness of reads was established -> second pass to assign reads to UCnodes
	for(int n=0;n<readlist.Count();n++) {

		float single_count=readlist[n]->read_count; // need to compute read single count

		for(int p=0;p<readlist[n]->pair_idx.Count();p++) { // for all pairs of read

			// for each pair I have to pretend strand is independent of the previous one since I grouped together several reads which might be independent
			int sno=(int)readlist[n]->strand+1; // 0: negative strand; 1: zero strand; 2: positive strand (starting from -1,0,1) // I have to reset sno every time because it might change due to snop

			// check if I've seen read's pair and if yes get its readcol; at the least get read's pair strand if available
			int np=readlist[n]->pair_idx[p]; // pair read number

			if(np>-1) { // read pair still exists  -> need to update counts even if pair came before

				// see if I have the correct read strand
				int snop=(int)readlist[np]->strand+1;  // snop is the strand of pair read
				if(sno!=snop) { // different strands for read and pair

					if(sno==1) { // read n is on zero (neutral) strand, but pair has strand
						sno=snop;  // assign strand of pair to read
					}
					else if(snop!=1 && n<np) { // conflicting strands -> un-pair reads in the hope that one is right -> if it was already processed I do not need to unpair
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

					if(n<np) { // process the two reads together
						add_to_UCNode(sno,readlist,n,np,readlist[n]->pair_count[p],read2unode,usgbundle,nu,SG2UCnode); // populate UCNodes
						add_to_UCNode(sno,readlist,np,-1,readlist[n]->pair_count[p],read2unode,usgbundle,nu,SG2UCnode); // populate UCNodes
					}

				} // this ends if(np>-1) : reads didn't get split
			} // if(np>-1) : read pair exists and it wasn't deleted
		} // for(int p=0,...)

		// now I need to deal with the single count
		if(single_count>epsilon) { // my way of controlling for rounding errors

			add_to_UCNode(readlist[n]->strand+1,readlist,n,-1,single_count,read2unode,usgbundle,nu,SG2UCnode); // populate UCNodes

		}
	}

}

void add_neutral_pred(int start,int end,float pcov,int n,BundleData* bdata,int &geneno) {

	int len=end-start+1;

	if(len>mintranscriptlen) {
		//fprintf(stderr,"1 Store single prediction:%d - %d with cov=%f\n",cov->start, cov->end, cov->cov/cov->len());
		geneno++;
		char sign='.';
		CPrediction *p=new CPrediction(geneno-1, NULL, start, end, pcov, sign, len,n);

		GSeg exon(start, end);
		p->exons.Add(exon);
		p->exoncov.Add(pcov);
		if(longreads) {
			p->tlen=-p->tlen;
			p->longcov=pcov*len;
		}

		GList<CPrediction>& pred = bdata->pred;
		pred.Add(p);
	}
}

// s is strand; if neutral is 1;
void find_connected_cov(UCov *cov,GPVec<UCNode> &SG2UCnode,int i,int s,GArray<int> &result,GBitVec &visited) {

	if(cov->next) {
		int v=2*i;
		result.Add(v); // first in graph
		return;
	}

	// visited.reset(); // No need to do this if I take care of reset in a smarter way
	// visited keeps two elements for each node i: 2i and 2i+1 for first in cov and last in cov
	GVec <int> queue;

	int p=0;
	int v=2*i+1; // last in graph
	visited[v]=1;
	queue.Add(v);

	while(p<queue.Count()) {
		int q=queue[p];
		p++;
		result.Add(q);

		int r=q % 2; // if r then coverage is last; otherwise it is first
		q=int(q/2);

		if(r || !SG2UCnode[q]->covintv[s]->next) // if last in coverage or coverage is unique -> look at the children
			for(int j=0;j<SG2UCnode[q]->childlink[s];j++) {
				int u=SG2UCnode[q]->childlink[s][j]->prevSGnode;
				v=2*u; // first in coverage
				if(!visited[v] && SG2UCnode[u]->covintv[s]) {
					queue.Add(v);
					visited[v]=1;
				}
			}

		if(!r || !SG2UCnode[q]->covintv[s]->next)  // if first in coverage or the coverage node is unique -> also look at the parents
			for(int j=0;j<SG2UCnode[q]->parentlink[s];j++) {
				int u=SG2UCnode[q]->parentlink[s][j]->prevSGnode;
				v=2*u+1; // last in coverage
				if(!visited[v] && SG2UCnode[u]->covintv[s]) {
					queue.Add(v);
					visited[v]=1;
				}
			}
	}

	// now result contains all UCnodes that are connected to last interval cov in SG2UCnode[i] -> make sure it is sorted
	// result.Sort(); // I used a sorted GArray -> everytime I add it should be in order
}


// distributes connected components and coverages and return next node to procecss from SG2UCnode
int get_neutral_proportions(UCov *cov,BundleData* bdata,GPVec<UCNode> &SG2UCnode,int i,GBitVec &visited,int &geneno) {

	GArray<int> result(true,true);	// contains all the nodes that are connected to cov making sure they are sorted
	find_connected_cov(cov,SG2UCnode,i,1,result,visited);

	int n=result.Count(); // n should at least have one componnent
	GVec <float> perbp; // perbp coverages of the result intervals
 	GVec<float> pos; // positive proportions of the result intervals
	GVec<float> neg; // negative proportions of the result intervals
	GVec<int> covend;
 	float zero=0;
	GVec<float> poscovxlink(n,zero); // sum of all positive coverage_proportions x links entering the nodes in result
	GVec<float> negcovxlink(n,zero); // sum of all positive coverage_proportions x links entering the nodes in result
	GVec<float> sumlink(n,zero); // sum of all links entering the nodes in result

	// compute the c's
	int j=0;
	uint start=cov->start;
	uint end=cov->end;
	float c=cov->cov;
	while(j<n) {

		visited[j]=0;
		int k=j+1;
		while(k<n && result[k]==result[k-1]+1) { // consecutive node

			int q=result[k]/2;
			UCov *qcov=SG2UCnode[q]->covintv[1];
			if(result[k]%2) { // k is the last node in cov interval
				while(qcov->next) { // only add coverages that I did not see before
					qcov=qcov->next; // get to last coverage interval
					c+=qcov->cov;
				}
			}
			else { // first in interval
				c+=qcov->cov;
			}
			end=qcov->end;

			visited[k]=0;
			k++;
		}

		covend.Add(k);

		c/=(end-start+1); // perbp coverage of interval
		float posprop=get_cov(2, start-bdata->start, end-bdata->start, bdata->bpcov);
		float negprop=get_cov(0, start-bdata->start, end-bdata->start, bdata->bpcov);

		if(posprop) {
			posprop=posprop/(posprop+negprop);
			if(posprop<epsilon) { posprop=0; negprop=1;}
			else {
				negprop=1-posprop;
				if(negprop<epsilon) negprop=0;
			}
		}
		else if(negprop) {
			negprop=negprop/(posprop+negprop);
			if(negprop<epsilon) { posprop=1; negprop=0;}
			else {
				posprop=1-negprop;
				if(posprop<epsilon) posprop=0;
			}
		}


		// I know links coming here from before

		perbp.Add(c);

		for(int u=j;u<k;u++) { // find all connections that leave nodes in interval

			// compute covxlink and sumlink
			if(result[u]%2>1 && (posprop || negprop)) { // if at last node
				int x=result[u]/2; // I am at node x
				for(int l=0;l<SG2UCnode[x]->childlink[1].Count();l++) {
					int lx=2*SG2UCnode[x]->childlink[1][l]->prevSGnode;
					if(visited[lx]) { // if lx is outside of my j to k interval
						int p=find_node_in_array(lx,result,k,n-1); // binary search of lx node position in result
						if(p<0) GError("Node %d in bundle starting at %d not among connected components\n",lx/2,bdata->start); // this should never happen
						if(posprop || negprop) {
							poscovxlink[p]+=posprop*SG2UCnode[x]->childcov[1][l];
							negcovxlink[p]+=negprop*SG2UCnode[x]->childcov[1][l];
							sumlink[p]+=SG2UCnode[x]->childcov[1][l];
						}
					}
				}
			}
			pos.Add(posprop);
			neg.Add(negprop);
		}

		if(k<n) {
			int q=result[k]/2;
			UCov *qcov=SG2UCnode[q]->covintv[1];
			if(result[k]%2) {
				while(qcov->next) qcov=qcov->next;
			}
			start=qcov->start;
			end=qcov->end;
			c=qcov->cov;
		}

		j=k;

	} // visited should be reset now

	// now I have to compute for each node the new coverages and covlinks
	j=0;
	int jj=0;
	start=cov->start;
	end=cov->end;
	while(j<n) {

		float jcov=perbp[jj];
		int jk=covend[jj];
		jj++;

		float posprop=0;
		float negprop=0;
		float sumbp=0;
		if(pos[j] || neg[j]) {
			posprop+=pos[j]*jcov;
			negprop+=neg[j]*jcov;
			sumbp+=jcov;
		}

		for(int u=j;u<jk;u++) { // consecutive coverages

			if(poscovxlink[u] || negcovxlink[u]) {
				posprop+=poscovxlink[u];
				negprop+=negcovxlink[u];
				sumbp+=sumlink[u];
			}

			if(result[u]%2>1) { // if at last node -> see if it connects to stranded interval
				int x=result[u]/2; // I am at node x
				for(int l=0;l<SG2UCnode[x]->childlink[1].Count();l++) {
					int lx=2*SG2UCnode[x]->childlink[1][l]->prevSGnode;
					if(lx>result[jk-1]) { // if lx is outside of my j to k interval
						int p=find_node_in_array(lx,result,jk,n-1); // binary search of lx node position in result
						if(p<0) GError("Node %d in bundle starting at %d not among connected components\n",lx/2,bdata->start); // this should never happen
						if(pos[p] || neg[p]) {
							posprop+=pos[p]*SG2UCnode[x]->childcov[1][l];
							negprop+=neg[p]*SG2UCnode[x]->childcov[1][l];
							sumbp+=sumlink[p];
						}
					}
				}
			}
		}

		// distribute coverages for nodes from j to k
		if(posprop || negprop) { // distribute coverages among the strands



			for(int k=j;k<jk;k++) {
				int q=result[k]/2;
				UCov *qcov=SG2UCnode[q]->covintv[1];
				if(qcov) {
					if(k>j && result[k]%2) { // k is the last node in cov interval
						while(qcov->next) { // add all coverages up to this point

							if(posprop) add_cov2UCNode(SG2UCnode,q,qcov->start,qcov->end,2,posprop*qcov->cov);
							if(negprop) add_cov2UCNode(SG2UCnode,q,qcov->start,qcov->end,0,negprop*qcov->cov);
							qcov=qcov->next; // get to last coverage interval
						}
					}
					else { // first in interval


						if(posprop) add_cov2UCNode(SG2UCnode,q,qcov->start,qcov->end,2,posprop*qcov->cov);
						if(negprop) add_cov2UCNode(SG2UCnode,q,qcov->start,qcov->end,0,negprop*qcov->cov);
					}

				// also distribute the coveragelinks
				for(int l=0;l<SG2UCnode[x]->childlink[1].Count();l++) {
					int lx=SG2UCnode[x]->childlink[1][l]->prevSGnode;

					if(posprop) insert_covlink(SG2UCnode,x,lx,2,posprop*SG2UCnode[x]->childcov[1]);
					if(negprop) insert_covlink(SG2UCnode,x,lx,0,negprop*SG2UCnode[x]->childcov[1]);

				}

				// delete coverages from node?




		}
		else { // neutral prediction
			if(jk-2>j && result[jk-1]%2) {
				int q=result[jk-1]/2;
				UCov *qcov=SG2UCnode[q]->covintv[1];
				while(qcov->next) qcov=qcov->next; // get to last coverage interval
				end=qcov->end;
			}

			int x=result[j]/2; // I am at node x
			add_neutral_pred(start,end,jcov,x,bdata,geneno);
			}

		if(jk<n) {

			int q=result[jk]/2;
			UCov *qcov=SG2UCnode[q]->covintv[1];
			if(result[jk]%2) {
				while(qcov->next) qcov=qcov->next;
			}
			start=qcov->start;
			end=qcov->end;
			}

		j=jk;

	}
}


void predict_neutral(BundleData* bdata,GPVec<UCNode> &SG2UCnode, int &geneno) {

	//int refstart = bdata->start;  // reference start

	int n=SG2UCnode.Count();
	float zero=0;
	GVec<float> pospropprior(n,zero); // remember positive prior proportions from previous nodes if priorcov are also present
	GVec<float> perbpprior(n,zero); // if this is non-negative then p-=1-p+

	GBitVec visited(2*n); // set it here so I don't have to recreate and reset it all the time

	for(int i=0;i<SG2UCnode.Count();i++) {
		if(SG2UCnode[i] && SG2UCnode[i]->covintv[1]) { // there are nodes covered in interval

			UCov *cov=SG2UCnode[i]->covintv[1];
			UCov *prevcov=NULL;

			/**** to do with stranded coverages, but should I do this with neutral ones too? */
			UCov *nextcov=NULL;

			// merge nearby coverages together
			while(nextcov) {
				if(nextcov->start-cov->end<bundledist) {
					float covleft=cov->cov/cov->len();
					float covright=nextcov->cov/nextcov->len();
					if(covleft>DROP*covright || covright>DROP*covleft) { // comparable coverages -> I might comment this out
						cov->cov+=nextcov->cov;
						cov->end=nextcov->end;
						cov->next=nextcov->next;
						delete nextcov;
						continue;
					}
				}
				cov=nextcov;
				nextcov=nextcov->next;
			}
			cov=SG2UCnode[i][1];
			/**** end merging */

			while(cov) { // this is not the last coverage for node
				bool toprint=true;

				float posprop=0;
				float negprop=0;
				get_strand_proportions(cov,bdata,SG2UCnode,i,posprop,negprop);

				float posprop=get_cov(2, cov->start-refstart, cov->end-refstart, bdata->bpcov);
				float negprop=get_cov(0, cov->start-refstart, cov->end-refstart, bdata->bpcov);

				float perbp=0;
				if(posprop||negprop) {
					posprop=posprop/(posprop+negprop);
					if(posprop<epsilon) { posprop=0; negprop=1;}
					else {
						negprop=1-posprop;
						if(negprop<epsilon) negprop=0;
					}
					perbp=cov->cov/cov->len;
				}

				if(perbpprior[i]) { // there is coverage coming in
					posprop=(posprop*perbp+pospropprior[i])/(perbp+perbpprior[i]);
					if(posprop<epsilon) { posprop=0; negprop=1;}
					else {
						negprop=1-posprop;
						if(negprop<epsilon) negprop=0;
					}
					perbpprior[i]=0; // used the priors
				}


				if(!cov->next) { // this is the last coverage in interval -> I need to compute future priors

					// now cov is the last covered interval in UCnode

					for(int i=0;i<SG2UCnode[i]->childlink[1].Count();i++) {

					}

				}

				if(posprop) {
					add_cov2UCNode(SG2UCnode,i,cov->start,cov->end,2,posprop*cov->cov);
					toprint=false;
				}
				if(negprop) {
					add_cov2UCNode(SG2UCnode,i,cov->start,cov->end,0,negprop*cov->cov);
					toprint=false;
				}

				if(toprint) { // add neutral bundle to predictions if needed
					add_neutral_pred(cov,i,bdata,geneno);
				}
				else { // cov was distributed among the +/-
					if(prevcov) prevcov->next=cov->next;
					else SG2UCnode[i]->covintv[1]=cov->next;
				}

				// I delete processed coverages
				UCNode* tmpcov=cov;
				cov=cov->next;
				delete tmpcov;

			}







		}
	}

}

// WIP FIXME
int build_usg(BundleData* bdata) {

	int refstart = bdata->start;  // reference start
	GList<CReadAln>& readlist = bdata->readlist; // all reads in bundle
	GList<CJunction>& junction = bdata->junction; // all junction in bundle
	GPVec<GffObj>& guides = bdata->keepguides; // all guides in bundle
	GVec<float>* bpcov = bdata->bpcov; // I might want to use a different type of data for bpcov to save memory in the case of very long bundles
	GList<CPrediction>& pred = bdata->pred;

	SGBundle* usgbundle = bdata->usgbundle;
	if(!usgbundle) return(0);  // no USG bundle to process

	GPVec<UCNode> SG2UCnode; // links SGnodes to a UCnode;
	GVec<int> read2unode; // remembers position of read in USG
	populate_UCNodes(bdata,SG2UCnode,read2unode);

	int geneno=0;

	// merge and predict neutral groups
	predict_neutral(bdata,SG2UCnode,geneno);

	// build stranded groups


	for(int n=0;n<readlist.Count();n++) { // for all reads in bundle


	}


    //TODO: don't forget to clean up the allocated data here
}
