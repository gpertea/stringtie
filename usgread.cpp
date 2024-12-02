//#define GMEMTRACE 1  //debugging memory allocation
#ifdef GMEMTRACE
#include "proc_mem.h"
#endif

#include "rlink.h"

// this functions merges a read to ugroup and returns the ugroup position in vector
//WIP FIXME
void merge_read_to_ugroup(int n, float readcov, int sno, GList<CReadAln>& readlist, UGroup **ugroup, UGroup **currugroup,
		                  SGBundle* ubundle) {

	if(currugroup[sno]) { // group exists

	}
	else { // this is the first time I see a read from this group

	}
}

void merge_ugroups(UGroup *grp1, UGroup *grp2) {
	// grp2 is added to grp1 and then it is deleted from the links
    //WIP FIXME
	//int i=
	//if(grp1->unode.Last()->end<grp2->unode)

}


void merge_pair_to_ugroup(int n, int np, float readcov, int sno,GList<CReadAln>& readlist, UGroup **ugroup, UGroup **currugroup,
		                  SGBundle* ubundle, UGroup **pairstart) {


	merge_read_to_ugroup(n,readcov,sno,readlist,ugroup,currugroup,ubundle);
	merge_read_to_ugroup(np,readcov,sno,readlist,ugroup,pairstart,ubundle);

	if(currugroup[sno]!=pairstart[sno]) merge_ugroups(currugroup[sno],pairstart[sno]);

}

void add_read_to_ugroup(int n, GList<CReadAln>& readlist, SGBundle* ubundle, GVec<int> &read2unode, GPVec<UGroup> *unode2ugroup) {


	float single_count=readlist[n]->read_count; // need to compute read single count
	int sno=readlist[n]->strand+1; // 0: negative strand; 1: zero strand; 2: positive strand (starting from -1,0,1) // I have to reset sno every time because it might change due to snop

	for(int p=0;p<readlist[n]->pair_idx.Count();p++) { // for all pairs of read


		// check if I've seen read's pair and if yes get its readcol; at the least get read's pair strand if available
		int np=readlist[n]->pair_idx[p]; // pair read number

		if(np>-1 && readlist[np]->nh) { // read pair still exists and it wasn't deleted

			// see if I have the correct read strand
			char snop=readlist[np]->strand+1;  // snop is the strand of pair read

			if(sno!=snop) { // different strands for read and pair

				if(snop==1) snop=sno; // pair gets strand from read
				else if(sno!=1) { // conflicting strands -> un-pair reads in the hope that one is right
					readlist[n]->pair_idx[p]=-1;
					for(int j=0;j<readlist[np]->pair_idx.Count();j++) // also unpair read np to n
						if(readlist[np]->pair_idx[j]==n) {
							readlist[np]->pair_idx[j]=-1;
							break;
						}
					np=-1;
				}
			}

			if(np>-1) { // reads didn't get split -> treat them together here

				single_count-=readlist[n]->pair_count[p]; // update single count of read here

				if(n<np) { // read pair comes after read n
					float readcov=readlist[n]->pair_count[p];
					//WIP FIXME - ugroup, currugroup, pair_start undefined
					//merge_pair_to_ugroup(n,np,readcov,snop,readlist,ugroup,currugroup,ubundle,pair_start);
				}
			} // this ends if(np>-1)
		} // if(np>-1 && readlist[np]->nh) : read pair exists and it wasn't deleted
	} // for(int p=0,...)




	// now I need to deal with the single count
	if(single_count>epsilon) { // my way of controlling for rounding errors
	    //WIP FIXME - ugroup, currugroup, undefined
		//merge_read_to_ugroup(n,single_count,sno,readlist,ugroup,currugroup,ubundle);
	}


}

// WIP FIXME
void build_usg(BundleData* bdata, GVec<int> &read2unode) {

	int refstart = bdata->start;  // reference start
	GList<CReadAln>& readlist = bdata->readlist; // all reads in bundle
	GList<CJunction>& junction = bdata->junction; // all junction in bundle
	GPVec<GffObj>& guides = bdata->keepguides; // all guides in bundle
	GVec<float>* bpcov = bdata->bpcov; // I might want to use a different type of data for bpcov to save memory in the case of very long bundles
	GList<CPrediction>& pred = bdata->pred;

	SGBundle* usgbundle = bdata->usgbundle;
	GPVec<UGroup> unode2ugroup[usgbundle->nodes.Count()];


	for(int n=0;n<readlist.Count();n++) { // for all reads in bundle
		if(readlist[n]->nh) add_read_to_ugroup(n,readlist,usgbundle,read2unode,unode2ugroup);
	}

	for(int i=0;i<usgbundle->nodes.Count();i++) {
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



    //TODO: don't forget to clean up the allocated data here
}
