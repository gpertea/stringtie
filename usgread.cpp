//#define GMEMTRACE 1  //debugging memory allocation
#ifdef GMEMTRACE
#include "proc_mem.h"
#endif

#include "usgread.h"

extern uint sserror;

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

int find_USG_predecessor(SGBundle* usgbundle,int low,int high,int val) {
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

int add2child(UCNode *node,char s,UCNode *child) {

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


void insert_covlink(GPVec<UCNode> &SG2UCnode,int n,int nx,char s,float count) {
	if(!SG2UCnode[n]) {
		SG2UCnode[n]= new UCNode();
		SG2UCnode[n]->prevSGnode=n;
	}
	if(!SG2UCnode[nx]) {
		SG2UCnode[nx]= new UCNode();
		SG2UCnode[n]->prevSGnode=nx;
	}
	// search if SG2UCnode[nx] was inserted before
	int p=add2child(SG2UCnode[n],s,SG2UCnode[nx]);
	SG2UCnode[n]->childlink[s]=SG2UCnode[nx];
	SG2UCnode[n]->childcov[s][p]+=count;
}

void add_cov2UCNode(GPVec<UCNode> &SG2UCnode,int n,uint start,uint end,char s,float count) {
	if(!SG2UCnode[n]) {
		SG2UCnode[n]= new UCNode();
		SG2UCnode[n]->prevSGnode=n;
		//SG2UCnode[n]->
	}

	/*
	UCov *prevcov=NULL;
	UCov *cov=SG2UCnode[n]->covintv[s];
	while(cov && start>cov->end) {
		prevcov=cov;
		cov=cov->next; // skip cov until I get to the one that overlaps my coverage interval or after it
	}

	UCov *nextcov;
	if(cov) { // start <= cov->end, start > prevcov->end ; there might be overlap with maxint

		if(end<cov->start) { // --> no overlap -> covintv is in between covg of node -> I need to keep previous cov for this
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
					nextmaxint=new CMaxIntv(maxint->node,start,maxint->end,0,maxint->next);
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

				update_overlap(pred,p,e,nextmaxint->node, overlaps);
				nextmaxint->node.Add(ex);

				if(end<nextmaxint->end) { // nextmaxint now ends where previous maxint ended
					maxint=new CMaxIntv(nextmaxint->node,end+1,nextmaxint->end,0,nextmaxint->next);
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
						update_overlap(pred,p,e,maxint->next->node, overlaps);
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
						CMaxIntv *newintv=new CMaxIntv(maxint->node,end+1,maxint->end,0,maxint->next);
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
	}*/

}

void add_to_UCNode(char s,GList<CReadAln>& readlist,int n,int np,float count,GVec<int> &read2unode,SGBundle* usgbundle,int nu,GPVec<UCNode> &SG2UCnode) {

	CReadAln & rd=*(readlist[n]);
	int u=read2unode[n]; // previous node in USG bundle from read
	if(u<0) u=0; // start unode

	int nex=rd.segs.Count();

	int i=0;
	while(i<nex) {

		if(rd.segs[i].end<usgbundle->nodes[0]->position) continue; // skip exons before the start of bundle
		if(u==nu-1) return; // I do not need to create links after this node

		uint start=rd.segs[i].start;
		uint end=rd.segs[i].end;

		int su=u;
		uint pos=usgbundle->nodes[0]->position;

		//if(usgbundle->nodes[0]->type==JSTART || usgbundle->nodes[0]->type==TEND) pos++; // TODO: check -> assuming that I need my nodes to start at jstart

		if(start<usgbundle->nodes[su]->position) start=usgbundle->nodes[su]->position;
		int lu=nu;
		if(i+1<nex) { // there is another exon after this one
			// how far would this read exon go
			if(rd.juncs[i]->usg_start>-1) lu=rd.juncs[i]->usg_start;
			else {
				lu=-3-rd.juncs[i]->usg_start;
				if(end-usgbundle->nodes[lu]->position<sserror && usgbundle->nodes[lu]->type==JSTART) end=usgbundle->nodes[lu]->position-1;
			}

			// calculate next exon in read
			if(rd.juncs[i]->usg_end>-1) u=rd.juncs[i]->usg_end;
			else u=-3-rd.juncs[i]->usg_end;

		}

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
	if(np>-1 && u<read2unode[np]) insert_covlink(SG2UCnode,u,read2unode[np],s,count);

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

		char rdstrand=rd.strand; // -1,0,1: I need to populate nodes but I can't do this until I know that the strand of the read is confirmed -> need to do 2 passes when nex=1

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
			char sno=readlist[n]->strand+1; // 0: negative strand; 1: zero strand; 2: positive strand (starting from -1,0,1) // I have to reset sno every time because it might change due to snop

			// check if I've seen read's pair and if yes get its readcol; at the least get read's pair strand if available
			int np=readlist[n]->pair_idx[p]; // pair read number

			if(np>-1) { // read pair still exists  -> need to update counts even if pair came before

				// see if I have the correct read strand
				char snop=readlist[np]->strand+1;  // snop is the strand of pair read
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


	for(int i=0;i<usgbundle->nodes.Count();i++) { // mark junctions that are in the USG
		if(usgbundle->nodes[i]->type==JSTART) { // junction
			for(int j=0;j<usgbundle->nodes[i]->jxLinks.Count();j++) {
				char s=0; // unknown strand
				if(usgbundle->nodes[i]->jxLinks[j].strand=='+') s=1; // junction on positive strand
				else if(usgbundle->nodes[i]->jxLinks[j].strand=='-') s=-1; // guide on negative strand
				CJunction jn(usgbundle->nodes[i]->position,usgbundle->nodes[i]->jxLinks[j].pos,s);
				int oidx=-1;
				if (junction.Found(&jn, oidx)) { // junction present in data
					junction[oidx]->usg_start=i;
					junction[oidx]->usg_end=usgbundle->nodes[i]->jxLinks[j].jx->bidx;
				}
			}

		}
	}


	for(int n=0;n<readlist.Count();n++) { // for all reads in bundle


	}


    //TODO: don't forget to clean up the allocated data here
}
