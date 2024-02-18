#include "bundle.h"

const byte txNASCENT_MASK=0x0C; // 2 bit flag for synthetic nascent transcripts
//0 = not nascent, 1 = nascent, 2 = nascent replacing a guide

const byte txSTATUS_MASK=0x03; // 2 bits mask for bundle status, former RC_TData::in_bundle
           // 1 if used in a bundle (guide added to keepguides, default value)
	         // 2 if all introns are covered by at least one read 
					 // 3 if it is stored to be printed
GuideBundleStatus getGuideStatus(GffObj* t) {
  return (GuideBundleStatus)t->getUserFlags(txSTATUS_MASK);
}
void setGuideStatus(GffObj* t, GuideBundleStatus status) {
  if (status>3) {
    GError("Error: invalid status value %d set for transcript %s\n", status, t->getID());
    exit(1);
  }
  t->setUserFlags(txSTATUS_MASK & status); // can be 1, 2 or 0 (clear)
}

// set/get synthetic nascent flag for a transcript
void setNascent(GffObj* t, byte nasc_st) { // 0, 1 or 2
      if (nasc_st>2) {
        GError("Error: invalid nascent status value %d set for transcript %s\n", nasc_st, t->getID());
        exit(1);
      }
      if (nasc_st) t->setUserFlags((nasc_st<<2) & txNASCENT_MASK);
      else t->clearUserFlags(txNASCENT_MASK);
}

byte isNascent(GffObj* t) {
      byte v=t->getUserFlags(txNASCENT_MASK);
      return (v>>2);
}



// genNascent is included from rlink.h
void BundleData::keepGuide(GffObj* t, Ref_RC_Data& ref_rc) {
	if (rc_data==NULL) {
	  rc_init(t, ref_rc.rc_tdata, ref_rc.rc_edata, 
               ref_rc.rc_idata);
	}
	keepguides.Add(t);
  if (isNascent(t)) {
    // add it to corresponding refdata.synrnas
    RC_TData* tdata=new RC_TData(*t, ref_rc.rc_tdata->Count());
	  t->uptr=tdata;
	  ref_rc.rc_tdata->Add(tdata);
    ref_rc.refdata->synrnas.Add(t);
  }
	t->udata=(int)rc_data->addTranscript(*t);
}

//multi-exon transcripts basic comparison/ordering function
// based on intron chain comparison
int txCompareProc(const pointer p1, const pointer p2) {
  //const GffObj* a = static_cast<const GffObj*>(p1);
  //const GffObj* b = static_cast<const GffObj*>(p2); 
  GffObj* a = (GffObj*)p1;
  GffObj* b = (GffObj*)p2;
  // check chromosome and strand
  if (a->gseq_id != b->gseq_id) return a->gseq_id < b->gseq_id ? -1 : 1;
  if (a->strand != b->strand) return a->strand < b->strand ? -1 : 1;
  // compare intron chains by their start and end, regardless of exon count
  int min_exon_count = GMIN(a->exons.Count(), b->exons.Count());
  for (int i = 0; i < min_exon_count-1; ++i) {
        int a_istart = a->exons[i]->end+1;
        int a_iend = a->exons[i+1]->start-1;
        int b_istart = b->exons[i]->end+1;
        int b_iend = b->exons[i+1]->start-1;

        if (a_istart != b_istart) return a_istart < b_istart ? -1 : 1;
        if (a_iend != b_iend) return a_iend < b_iend ? -1 : 1;
  }
  // after comparing shared introns, tx with more exons is "greater"
  if (a->exons.Count() != b->exons.Count()) {
      return a->exons.Count() < b->exons.Count() ? -1 : 1;
  }
  // intron chains are identical, they are considered "equal"
  return 0;
}

//single-exon transcripts basic comparison/ordering function
int seTxCompareProc(pointer* p1, pointer* p2) {
  // check chromosome and strand
  GffObj* a = (GffObj*)p1;
  GffObj* b = (GffObj*)p2;
  if (a->gseq_id != b->gseq_id) return a->gseq_id < b->gseq_id ? -1 : 1;
  if (a->strand != b->strand) return a->strand < b->strand ? -1 : 1;
  //first by start positions
  if (a->start != b->start) return a->start < b->start ? -1 : 1;
  // then by end positions
  if (a->end != b->end) return a->end < b->end ? -1 : 1;
  return 0; // they are equal (both start and end match)
}

//generate all nascent transcripts for a guide
//keep nascents in 2 lists: tnlist for multi-exon, setnlist for single-exon
void genTxNascents(GffObj &t, GList<GffObj>& tnlist, GList<GffObj>& setnlist, GPVec<GffObj>& guides) {
    int n = t.exons.Count() - 1; // number of introns
    for (int i = 1; i <= n; ++i) {
        GffObj* nt = new GffObj(); // Create a new transcript
        nt->isTranscript(true);
        nt->ftype_id=gff_fid_transcript;
        nt->subftype_id=gff_fid_exon;
        nt->gseq_id = t.gseq_id; // copying reference ID, strand, gene_id, gene_name, etc.
        nt->strand = t.strand;
        nt->isFinalized(true);
        if (t.strand == '-') { // on the reverse strand, we need to add the exons in reverse order
          for (int j=n; j>n-i;--j) {
              nt->addExon(t.exons[j]->start, t.exons[j]->end);
          }
          // Extend the first exon of nt through the previous intron up to the start of the previous exon in t
          nt->exons[0]->start = t.exons[n-i]->end + 1; 
        } else { //forward strand, add exons in the same order as in t
          for (int j = 0; j < i; ++j) {
              // add i exons from t to nt
              nt->addExon(t.exons[j]->start, t.exons[j]->end);
          }
          // Extend the last exon of nt through the next intron up to the start of the next exon in t
          nt->exons.Last()->end = t.exons[i]->start-1;
        }
          
        bool keepNascent = true;
        //TODO: check here for txStructureMatch() with any tx in guides 
        //   set keepNascent to false if a matching guide is found
        //   (or even replace the guide if we decide to keep the nascent instead)
        
        if (nt->exons.Count() > 1) {
            /* 
            int fidx=-1;
            if (tnlist.Found(nt, fidx)) { //uses special intron-chain comparison
              if (tnlist[fidx]->covlen < nt->covlen) {
                //TODO check this
                delete tnlist[fidx]; //discard the old one
                tnlist.Put(fidx, nt); //replace with the new one
              } else { //TODO: check this
                 keepNascent = false; // we already have a nascent with the same structure (intron chain)
              }
            }
            else */ 
               tnlist.Add(nt); // Add the new transcript to tnlist
          } else { // single-exon nascent
            /*
            int fidx=-1;
            if (setnlist.Found(nt, fidx)) {
              if (setnlist[fidx]->covlen < nt->covlen) {
                //TODO check this
                delete setnlist[fidx]; //discard the old one
                setnlist.Put(fidx, nt); //replace with the new one
              } else {
                keepNascent = false; // we already have a nascent with the same structure (intron chain)
              } 
            }
            else {*/

              //TODO: check for matching SE nascent in setnlist from fidx
              setnlist.Add(nt); // Add the new transcript to setnlist
              /*
              for (int p = fidx; p < setnlist.Count(); ++p) {
                   GffObj* ctx = setnlist[p];
                   if (ctx->start > nt->end) 
                        break; // further transcripts cannot overlap, so break the loop
                   if (has80PercentOverlap(nt, ctx)) {
                    // Found a transcript that overlaps with at least 80%
                    // TODO Process this match as needed
                   }
              } //for p
            }*/
          } //single-exon nascent
        if (keepNascent) {
          GStr nid(t.getID());
          nid.appendfmt(".nasc%03d", i);
          nt->setGeneID(t.getGeneID());
          nt->setGeneName(t.getGeneName());
          nt->setID(nid.chars());
          setNascent(nt, 1); // Set the nascent flag
        } else 
           delete nt; //discard this nascent transcript (matching guide/other nascent found)
    } //for i 
}

void BundleData::generateAllNascents(int from_guide_idx, Ref_RC_Data& ref_rc) {
  //GPVec<GffObj> oguides(keepguides); // copy list of pointers to guides
  // tnlist is a sorted list of nascent transcripts to be added to keepguides
  GList<GffObj> tnlist(true, true); //list to collect all nascents in this bundle
  tnlist.setSorted(txCompareProc);
  GList<GffObj> setnlist(true, true); //list to collect all nascents in this bundle
  setnlist.setSorted((GCompareProc*)seTxCompareProc);
  for (int i=from_guide_idx;i<keepguides.Count();++i) {
    GffObj* t=keepguides[i];
    if (t->exons.Count()<2) continue;
    GMessage("Generating %d nascents for %s(%c)\n", t->exons.Count()-1, t->getID(), t->strand);
    genTxNascents(*t, tnlist, setnlist, keepguides);
  }
  //TODO add all nascents from tnlist and setnlist to keepguides
  for (int i=0;i<tnlist.Count();++i) {
    //keepguides.Add(tnlist[i]);
    keepGuide(tnlist[i], ref_rc);
    tnlist.Put(i, NULL); //don't delete the transcript, it's now owned by keepguides
  } //TODO: CHECK THIS
  for (int i=0;i<setnlist.Count();++i) {
    //keepguides.Add(setnlist[i]);
    keepGuide(setnlist[i], ref_rc);
    setnlist.Put(i, NULL); //don't delete the transcript, it's now owned by keepguides
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
		 alndata.juncs.Add(new CJunction(j_l, j_r)); //don't set strand, etc. for now
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
