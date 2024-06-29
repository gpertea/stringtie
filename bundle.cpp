#include "bundle.h"

const byte txSTATUS_MASK=0x03; // 2 bits mask for bundle status, former RC_TData::in_bundle
           // 1 if used in a bundle (guide added to keepguides, default value)
	         // 2 if all introns are covered by at least one read 
					 // 3 if it is stored to be printed
GuideBundleStatus getGuideStatus(GffObj* guide) {
  return (GuideBundleStatus)guide->getUserFlags(txSTATUS_MASK);
}
void setGuideStatus(GffObj* guide, GuideBundleStatus status) {
  if (status>3) {
    GError("Error: invalid status value %d set for transcript %s\n", status, guide->getID());
    exit(1);
  }
  guide->setUserFlags(txSTATUS_MASK & status); // can be 1, 2 or 0 (clear)
}

void BundleData::keepGuide(GffObj* guide, Ref_RC_Data& ref_rc) {
	if (rc_data==NULL) {
	  rc_init(guide, ref_rc.rc_tdata, ref_rc.rc_edata, 
               ref_rc.rc_idata);
	}
	keepguides.Add(guide);
	guide->udata=(int)rc_data->addTranscript(*guide); //this also adds exon/intron info
}

bool BundleData::evalReadAln(GReadAlnData& alndata, char& xstrand) {
	            //GSamRecord& brec, char& strand, int nh) {
 if (rc_data==NULL) {
	  return false; //no ref transcripts available for this reads' region
 }
 GSamRecord& brec=*(alndata.brec);
 int mate_pos=brec.mate_start();
 int nh=alndata.nh;
 if ((int)brec.end<rc_data->lmin || (int)brec.start>rc_data->rmax) {
	 return false; //hit outside coverage area
 }
 if (rc_data->g_exons.Count()==0 || rc_data->g_tdata.Count()==0)
	 return false; //nothing to do without transcripts
 //check this read alignment against ref exons and introns
 char strandbits=0;
 bool result=false;
 bool is_in_guide=true; // exons and junctions are in reference transcripts but they might be in different guides
 for (int i=0;i<brec.exons.Count();i++) {
	 if (ballgown)
		 rc_data->updateCov(xstrand, nh, brec.exons[i].start, brec.exons[i].len());
	 GArray<RC_ExonOvl> exonOverlaps(true, true); //overlaps sorted by decreasing length
	 if (rc_data->findOvlExons(exonOverlaps, brec.exons[i].start,
			 brec.exons[i].end, xstrand, mate_pos)) {
		 result=true;
		 int max_ovl=exonOverlaps[0].ovlen;
		 if(is_in_guide && (uint)max_ovl<brec.exons[i].len()) is_in_guide=false;
		 //alndata.g_exonovls.Add(new GVec<RC_ExonOvl>(exonOverlaps));
			 for (int k=0;k<exonOverlaps.Count();++k) {
				 //if (exonOverlaps[k].ovlen < 5) break; //ignore very short overlaps
				 if (k && (exonOverlaps[k].mate_ovl < exonOverlaps[0].mate_ovl
						 || exonOverlaps[k].ovlen+5<max_ovl) )
					 break; //ignore further overlaps after a mate matched or if they are shorter than max_overlap-5
				 if (exonOverlaps[k].feature->strand=='+') strandbits |= 0x01;
				 else if (exonOverlaps[k].feature->strand=='-') strandbits |= 0x02;
				 //TODO: perhaps we could use a better approach for non-overlapping ref exons
				 //      spanned by this same read alignment
				 //counting this overlap for multiple exons if it's similarly large
				 //(in the shared region of overlapping exons)
				 rc_updateExonCounts(exonOverlaps[k], nh);
			 }
	 } //ref exon overlaps
	 if (i>0) { //intron processing
		 int j_l=brec.exons[i-1].end+1;
		 int j_r=brec.exons[i].start-1;
		 RC_Feature* ri=rc_data->findIntron(j_l, j_r, xstrand);
		 alndata.juncs.Add(new CJunction(j_l, j_r)); //don'guide set strand, etc. for now
		 if (ri) { //update guide intron counts
			 ri->rcount++;
			 ri->mrcount += (nh > 1) ? (1.0/nh) : 1;
			 if (nh==1)  ri->ucount++;
			 alndata.juncs.Last()->guide_match=1;
		 }
		 else is_in_guide=false;

	 } //intron processing
 }
 if (xstrand=='.' && strandbits && strandbits<3) {
	xstrand = (strandbits==1) ? '+' : '-';
 }

 if(!mergeMode && is_in_guide) alndata.in_guide=true;

 return result;
}
