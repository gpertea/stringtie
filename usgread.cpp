#include "rlink.h"


//#define GMEMTRACE 1  //debugging memory allocation
#ifdef GMEMTRACE
#include "proc_mem.h"
#endif


// this functions merges a read to ugroup and returns the ugroup position in vector
void merge_read_to_ugroup(int n,float readcov, int sno,GList<CReadAln>& readlist,UGroup **ugroup,UGroup **currugroup,
		SGBundle* ubundle) {

	if(currugroup[sno]) { // group exists

	}
	else { // this is the first time I see a read from this group

	}
}

void merge_ugroups(UGroup *grp1, UGroup *grp2) {

	// grp2 is added to grp1 and then it is deleted from the links

	int i=
	if(grp1->unode.Last()->end<grp2->unode)

}


void merge_pair_to_ugroup(int n,int np, float readcov, int sno,GList<CReadAln>& readlist,UGroup **ugroup,UGroup **currugroup,
		SGBundle* ubundle,UGroup **pairstart) {


	merge_read_to_ugroup(n,readcov,sno,readlist,ugroup,currugroup,ubundle);
	merge_read_to_ugroup(np,readcov,sno,readlist,ugroup,pairstart,ubundle);

	if(currugroup[sno]!=pairstart[sno]) merge_ugroups(currugroup[sno],pairstart[sno]);

}

void add_read_to_ugroup(int n,GList<CReadAln>& readlist,SGBundle* ubundle,GVec<int> &read2unode,GPVec<UGroup> *unode2ugroup) {


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
					merge_pair_to_ugroup(n,np,readcov,snop,readlist,ugroup,currugroup,ubundle,pair_start);
				}
			} // this ends if(np>-1)
		} // if(np>-1 && readlist[np]->nh) : read pair exists and it wasn't deleted
	} // for(int p=0,...)




	// now I need to deal with the single count
	if(single_count>epsilon) { // my way of controlling for rounding errors
		merge_read_to_ugroup(n,single_count,sno,readlist,ugroup,currugroup,ubundle);
	}


}

