#include "bundle.h"

const byte txNASCENT_MASK=0x0C; // 2 bit flag for synthetic nascent transcripts
//0 = not nascent, 1 = nascent, 2 = nascent replacing a guide

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

// set/get synthetic nascent flag for a transcript
void setNascent(GffObj* guide, byte nasc_st) { // 0, 1 or 2
      if (nasc_st>2) {
        GError("Error: invalid nascent status value %d set for transcript %s\n", nasc_st, guide->getID());
        exit(1);
      }
      if (nasc_st) guide->setUserFlags((nasc_st<<2) & txNASCENT_MASK);
      else guide->clearUserFlags(txNASCENT_MASK);
}

byte isNascent(GffObj* guide) {
      byte v=guide->getUserFlags(txNASCENT_MASK);
      return (v>>2);
}



// genNascent is included from rlink.h
void BundleData::keepGuide(GffObj* guide, Ref_RC_Data& ref_rc) {
	if (rc_data==NULL) {
	  rc_init(guide, ref_rc.rc_tdata, ref_rc.rc_edata, 
               ref_rc.rc_idata);
	}
	keepguides.Add(guide);
  if (isNascent(guide)) {
    // add it to corresponding refdata.synrnas
    RC_TData* tdata=new RC_TData(*guide, ref_rc.rc_tdata->Count());
	  guide->uptr=tdata;
	  ref_rc.rc_tdata->Add(tdata);
    ref_rc.refdata->synrnas.Add(guide);
  }
	guide->udata=(int)rc_data->addTranscript(*guide); //this also adds exon/intron info
}


//generate all nascent transcripts for a guide
//keep nascents in 2 lists: tnlist for multi-exon, setnlist for single-exon
void genTxNascents(GffObj &guide, GList<GffObj>& tnlist, GPVec<GffObj>& guides) {
    // DEBUG only: keeping track of the number of nascent transcripts generated vs redundant
    //GMessage("Generating %d nascents for %s(%c)\n", guide.exons.Count()-1, guide.getID(), guide.strand);
    int nxr=0; //number of redundant nascent transcripts discarded
    int n = guide.exons.Count() - 1; // number of introns
    for (int i = 1; i <= n; ++i) {
        GffObj* nt = new GffObj(true, nullptr, guide.gseq_id, guide.strand); // Create a new transcript
        if (guide.strand == '-') { // on the reverse strand, we need to add the exons in reverse order
          for (int j=n; j>n-i;--j) {
              nt->addExon(guide.exons[j]->start, guide.exons[j]->end);
          }
          // Extend the first exon of nt through the previous intron up to the start of the previous exon in guide
          nt->exons[0]->start = guide.exons[n-i]->end + 1; 
          nt->start = nt->exons[0]->start;
        } else { //forward strand, add exons in the same order as in guide
          for (int j = 0; j < i; ++j) {
              // add i exons from guide to nt
              nt->addExon(guide.exons[j]->start, guide.exons[j]->end);
          }
          // Extend the last exon of nt through the next intron up to the start of the next exon in guide
          nt->exons.Last()->end = guide.exons[i]->start-1;
          nt->end = nt->exons.Last()->end;
        }
          
        bool keepNascent = true;
         //   set keepNascent to false if a matching previous nascent is found
        int ridx=-1;
        if (tnlist.Found(nt, ridx)) { //exact match found with a previous nascent
          keepNascent=false;
          nxr++;
        } else { // no exact match in tnlist
          // might still be a structural match with a previous nascent
          // backward search --------------
          // set a limit to the backward search, to avoid checking too distant nascents
          int64_t nt_bck_limit=((int64_t)nt->start - (int64_t)(nt->end - nt->start)/2);
          
          for (int i = ridx - 1; i >= 0; --i) {
              // break condition for the backward search is not as straightforward as the forward search
              // because the end coordinates do not influence the sort order.
              if ((int64_t)tnlist[i]->start < nt_bck_limit) {
                  // if the current transcript's start is before the middle of nt, no previous transcripts
                  // can structurally match nt, because they will start even earlier.
                  break;
              }
              if (txStructureMatch(*nt, *tnlist[i])) {
                  keepNascent = false;
                  nxr++;
                  break;
              }
              // Note: An optimization here might consider the spatial distance between transcripts, but 
              // without a clear rule on how much overlap is needed for a match, it's safer not to break early.
          } // backward search
          
          if (keepNascent)
            //forward search
            for (int i = ridx; i < tnlist.Count(); ++i) {
                 if (tnlist[i]->start > nt->end) {
                        // If the current transcript's start is beyond the end of guide, no subsequent transcripts
                        // can structurally match guide, because they will start even later.
                        break;
                 }
                 if (txStructureMatch(*nt, *tnlist[i])) {
                      keepNascent = false;
                      nxr++;
                      break;
                 }
            }
        } 
        
        if (keepNascent) {
          // insert it in the right location in the list
          tnlist.sortInsert(ridx, nt);
          GStr nid(guide.getID());
          nid.appendfmt(".nasc%03d", i);
          nt->setGeneID(guide.getGeneID());
          nt->setGeneName(guide.getGeneName());
          nt->setID(nid.chars());
          setNascent(nt, 1); // Set the nascent flag
        } else 
           delete nt; //discard this nascent transcript (matching guide/other nascent found)
    } //for i 
    //DEBUG only: keeping track of the number of nascent transcripts generated vs redundant
    //if (nxr>0)  GMessage("   %d redundant nascents discarded\n", nxr);
}

void BundleData::generateAllNascents(int from_guide_idx, Ref_RC_Data& ref_rc) {
  //GPVec<GffObj> oguides(keepguides); // copy list of pointers to guides
  // tnlist is a sorted list of nascent transcripts to be added to keepguides
  GList<GffObj> tnlist(true, true); //list to collect all nascents in this bundle
  tnlist.setSorted(txCmpByExons); //sort and match by exon coordinates
  for (int i=from_guide_idx;i<keepguides.Count();++i) {
    GffObj* guide=keepguides[i];
    if (guide->exons.Count()<2) continue;
    genTxNascents(*guide, tnlist, keepguides);
  }
  //add all nascents from tnlist to keepguides
  numNascents+=tnlist.Count();
  for (int i=0;i<tnlist.Count();++i) {
    keepGuide(tnlist[i], ref_rc);
    tnlist.Put(i, NULL); //don't delete nascent, it's now owned by keepguides
  } 
  keepguides.Sort();
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
