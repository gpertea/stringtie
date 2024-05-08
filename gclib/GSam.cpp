#include "GSam.h"
#include <ctype.h>

#define _cigOp(c) ((c)&BAM_CIGAR_MASK)
#define _cigLen(c) ((c)>>BAM_CIGAR_SHIFT)

GSamRecord::GSamRecord(const char* qname, int32_t gseq_tid,
                 int pos, bool reverse, GDynArray<uint32_t>& cigar,
                 const char* qseq, const char* quals):iflags(0), exons(1),
				 juncsdel(1),clipL(0), clipR(0), mapped_len(0) {
   uint16_t flag=0;
   novel=true;
   b_hdr=NULL;
   b=bam_init1();
   size_t l_qname=strlen(qname);
   size_t l_seq=0;
   if (qseq!=NULL) {
	   l_seq=strlen(qseq);
   }
   //-- passed pos is 1-based
   if (pos<=0 || gseq_tid<0) {
	   flag |= BAM_FUNMAP;
	   pos=-1;
   }
   pos--;
   if (reverse) flag|=BAM_FREVERSE;

   bam_set1(b, l_qname, qname, flag,
		   gseq_tid, pos, 60, cigar.Count(), cigar(),
		   -1, 0, 0, l_seq, qseq, quals, 0);
}

void GSamRecord::add_aux(const char* str) {
     //requires: being called AFTER add_quals() for built-from-scratch records
     //--check the "// aux" section in sam_parse1() htslib/sam.c
     char tag[2];
     uint8_t abuf[512];

     int strl=strlen(str);
     //int doff = b->core.l_qname + b->core.n_cigar*4 + (b->core.l_qseq+1)/2 + b->core.l_qseq + b->l_aux;
     //int doff0=doff;
     if (strl < 6 || str[2] != ':' || str[4] != ':')
         parse_error("missing colon in auxiliary data");
     tag[0] = str[0]; tag[1] = str[1];
     uint8_t atype = str[3];
     uint8_t* adata=abuf;
     int alen=0;
     if (atype == 'A' || atype == 'a' || atype == 'c' || atype == 'C') { // c and C for backward compatibility
         atype='A';
         alen=1;
         adata=(uint8_t*)&str[5];
         }
      else if (atype == 'I' || atype == 'i') {

		 long long x=strtoll(str+5, NULL, 10); //(long long)atoll(str + 5);
		 //long x=(long)atol(str + 5);
		 if (x < 0) {
             if (x >= -127) {
                 atype='c';
                 abuf[0] =  (int8_t)x;
                 alen=1;
                 }
             else if (x >= -32767) {
                 atype = 's';
                 *(int16_t*)abuf = (int16_t)x;
                 alen=2;
                 }
             else {
                 atype='i';
                 *(int32_t*)abuf = (int32_t)x;
                 alen=4;
                 if (x < -2147483648ll)
                     GMessage("Parse warning: integer %lld is out of range.", x);
                 }
             } else { //x >=0
             if (x <= 255) {
                 atype = 'C';
                 abuf[0] = (uint8_t)x;
                 alen=1;
                 }
             else if (x <= 65535) {
                 atype='S';
                 *(uint16_t*)abuf = (uint16_t)x;
                 alen=2;
                 }
             else {
                 atype='I';
                 *(uint32_t*)abuf = (uint32_t)x;
                 alen=4;
                 if (x > 4294967295ll)
                     GMessage("Parse warning: integer %lld is out of range.", x);
                 }
             }
         } //integer type
         else if (atype == 'f') {
             *(float*)abuf = (float)atof(str + 5);
             alen = sizeof(float);
             }
         else if (atype == 'd') { //?
             *(float*)abuf = (float)atof(str + 9);
             alen=8;
             }
         else if (atype == 'Z' || atype == 'H') {
             if (atype == 'H') { // check whether the hex string is valid
                 if ((strl - 5) % 2 == 1) parse_error("length of the hex string not even");
                 for (int i = 0; i < strl - 5; ++i) {
                     int c = toupper(str[5 + i]);
                     if (!((c >= '0' && c <= '9') || (c >= 'A' && c <= 'F')))
                         parse_error("invalid hex character");
                     }
                 }
             memcpy(abuf, str + 5, strl - 5);
             abuf[strl-5] = 0;
             alen=strl-4;
         }
         else if (atype == 'B') { //FIXME
            // -- see sam_parse_B_vals() function in htslib/sam.c
            //Integer or numeric array -- too messy to parse again here
            GMessage("Warning: sorry, B tags not supported yet.\n");
         }
         else parse_error("unrecognized aux type");
  //this->add_aux(tag, atype, alen, adata);
  bam_aux_append(b, tag, atype, alen, adata);
 }//add_aux()

int interpret_CIGAR(char cop, int cl, int aln_start) {
// returns the number of bases "aligned" (matches or mismatches) from the read
// gpos = current genomic position (will end up as right coordinate on the genome)
// rpos = read position (will end up as the length of the read)
// cop = CIGAR operation, cl = operation length
int mbases = 0; //count "aligned" bases (includes mismatches)
int rpos = 0;
int gpos = aln_start;
int num_mismatches=0; //NM tag value = edit distance
switch (cop) {
 case BAM_CDIFF:  // X
      num_mismatches+=cl;
 	 //fall-through
 case BAM_CMATCH: // M
     //have to actually check for mismatches: num_mismatches+=count_mismatches;
	 //fall-through
 case BAM_CEQUAL: // =
      //printf("[%d-%d]", gpos, gpos + cl - 1);
      gpos+=cl;
      rpos+=cl;
      ++mbases;
      break;

 case BAM_CPAD:
      // printf("[%d-%d]", pos, pos + cl - 1); // Spans positions, No Coverage
      gpos+=cl;
      break;

 case BAM_CHARD_CLIP:
      // printf("[%d]", pos);  // No coverage
      // gpos is not advanced by this operation
      break;

 case BAM_CSOFT_CLIP: // S
      //soft clipped bases, present in SEQ
      rpos+=cl;
      break;

 case BAM_CINS: // I
      // No Coverage
      // adds cl bases "throughput" but not genomic position "coverage" (gap in genomic seq)
      // should also add cl to the number of "mismatches" (unaligned read bases)
      num_mismatches+=cl;
      // How you handle this is application dependent
      // gpos is not advanced by this operation
      rpos+=cl;
      break;

 case BAM_CDEL: // D
      //deletion in reference sequence relative to the read (gap in read sequence)
      // printf("[%d-%d]", pos, pos + cl - 1);
      // Spans positions
      num_mismatches+=cl;
      gpos += cl;
      break;

 case BAM_CREF_SKIP: // N
      // intron
      //special skip operation, not contributing to "edit distance",
      // printf("[%d-%d]", pos, pos + cl - 1); // Spans positions, No Coverage
      //   so num_mismatches is not updated
      gpos+=cl;
      break;

 default:
      fprintf(stderr, "Unhandled cigar_op %d:%d\n", cop, cl);
      //printf("?");
  }
 return mbases;
} // interpret_CIGAR(), just a reference of CIGAR operations interpretation

 void GSamRecord::setupCoordinates() {
	if (!b) return;
	const bam1_core_t *c = &b->core;
	if (c->flag & BAM_FUNMAP) return; /* skip unmapped reads */
	uint32_t *cigar = bam_get_cigar(b);
	int l=0;
	mapped_len=0;
	clipL=0;
	clipR=0;
	start=c->pos+1; //genomic start coordinate, 1-based (BAM core.pos is 0-based)
	int exstart=c->pos;
	GSeg exon;
	bool exonStarted=false; //keep track of exon start (M/X/= detected in this alignment)
	bool intron=false;
	bool ins=false;
	uint del=0;
	uint prevdel=0;

	for (uint i = 0; i < c->n_cigar; ++i) {
		unsigned char op = _cigOp(cigar[i]);
		switch(op) {
		  case BAM_CEQUAL:    // =
		  case BAM_CDIFF:     // X
		  case BAM_CMATCH:    // M
			exonStarted=true;
		    l+=_cigLen(cigar[i]);
			if(intron) { // op comes after intron --> update juncsdel
				GSeg deljunc(prevdel,0);
				juncsdel.Add(deljunc);
			}
		    intron=false;
		    ins=false;
		    del=0;
		    break;
		  case BAM_CDEL:
			del=_cigLen(cigar[i]);
			l+=del;
			if (intron) { // deletion after intron --> update juncsdel
				GSeg deljunc(prevdel,del);
				juncsdel.Add(deljunc);
			}
			ins=false;
			break;
		  case BAM_CINS:      // I
		    //rpos+=cl; //gpos not advanced
		    //take care of cases where there is an ins within an intron
		    ins=true;
		    break;
		  case BAM_CREF_SKIP: // N
			if (!exonStarted) break; // anomalous alignment starting with an intron before exon (uLTRA)
		    if(!ins || !intron) { // not preceded by insertion in the middle of an intron
		      //so we can close the current (preceding) exon
		      exon.end=c->pos+l;
		      exon.start=exstart+1;
		      exons.Add( exon ); //keep the preceding exon
		      mapped_len+=exon.len();
		    }
		    has_Introns=true;
		    l += _cigLen(cigar[i]);
		    exstart=c->pos+l;
			prevdel=del;
		    intron=true;
		    del=0;
		    break;
		  case BAM_CSOFT_CLIP: // S
		    soft_Clipped=true;
		    if (l) clipR=_cigLen(cigar[i]);
		      else clipL=_cigLen(cigar[i]);
		    //intron=false;
		    ins=false;
		    del=0;
		    break;
		  case BAM_CHARD_CLIP:
		    hard_Clipped=true;
		    //intron=false;
		    ins=false;
		    del=0;
		    break;
		  case BAM_CPAD: //?
		    //gpos+=cl;
		    break;
		  default:
		    int cl=_cigLen(cigar[i]);
		    fprintf(stderr, "Unhandled CIGAR operation %d:%d\n", op, cl);
		}
	}
	if (!intron) {
	  exon.start=exstart+1;
	  exon.end=c->pos+l;
	  exons.Add(exon);
	  mapped_len+=exon.len();
	}
	if (exon.end) end=exon.end; //genomic end coordinate
	if (end==0) GError("Error: invalid CIGAR record for %s !\n", bam_get_qname(b));
 }

 uint8_t* GSamRecord::find_tag(const char tag[2]) {
   return bam_aux_get(this->b, tag);
 }

 int GSamRecord::remove_tag(const char tag[2]) {
   uint8_t* p=bam_aux_get(this->b, tag);
   if (p!=NULL) return bam_aux_del(this->b, p);
   return 0;
 }


 char GSamRecord::tag_char(const char tag[2]) { //retrieve tag data as single char
   uint8_t* s=find_tag(tag);
   if (s) return ( bam_aux2A(s) );
   return 0;
 }

 char GSamRecord::tag_char1(const char tag[2]) { //just the first char from Z type tags
	uint8_t* s=bam_aux_get(this->b, tag);
	if (s==NULL) return 0;
 	int type;
 	type = *s++;
 	if (s == 0) return 0;
 	if (type == 'A' || type == 'Z') return *(char*)s;
 	else return 0;
 }

 int64_t GSamRecord::tag_int(const char tag[2], int nfval) { //get the numeric value of tag
   uint8_t *s=find_tag(tag);
   if (s) return ( bam_aux2i(s) );
   return nfval;
 }

 double GSamRecord::tag_float(const char tag[2]) { //get the float value of tag
    uint8_t *s=bam_aux_get(this->b, tag);;
    if (s) return ( bam_aux2f(s) );
    return 0;
 }

 char* GSamRecord::tag_str(const char tag[2]) { //return string value for a tag
   uint8_t *s=find_tag(tag);
   if (s) return ( bam_aux2Z(s) );
   return NULL;
 }

 char GSamRecord::spliceStrand() { // '+', '-' from the XS tag, or 0 if no XS tag
   char c=tag_char1("XS");
   if (c==0) {
    //try minimap2's "ts" tag
    char m=tag_char1("ts");
    if (m=='+' || m=='-') {
       if ((this->b->core.flag & BAM_FREVERSE) != 0) c=((m=='+') ? '-' : '+');
         else c=m;
    }
   }
   return ((c=='+' || c=='-') ? c : '.');
 }

 char* GSamRecord::sequence() { //user must free this after use
   char *s = (char*)bam_get_seq(b);
   char* qseq=NULL;
   GMALLOC(qseq, (b->core.l_qseq+1));
   int i;
   for (i=0;i<(b->core.l_qseq);i++) {
     int8_t v = bam_seqi(s,i);
     qseq[i] = seq_nt16_str[v];
     }
   qseq[i] = 0;
   return qseq;
 }

 char* GSamRecord::qualities() {//user must free this after use
   char *qual  = (char*)bam_get_qual(b);
   char* qv=NULL;
   GMALLOC(qv, (b->core.l_qseq+1) );
   int i;
   for(i=0;i<(b->core.l_qseq);i++) {
     qv[i]=qual[i]+33;
     }
   qv[i]=0;
   return qv;
 }

 char* GSamRecord::cigar() { //returns text version of the CIGAR string; must be freed by user
   kstring_t str = KS_INITIALIZE;
   if (b->core.n_cigar == 0) kputc('*', &str);
    else {
      for (uint i = 0; i < b->core.n_cigar; ++i) {
         kputw(bam_get_cigar(b)[i]>>BAM_CIGAR_SHIFT, &str);
         kputc(BAM_CIGAR_STR[bam_get_cigar(b)[i]&BAM_CIGAR_MASK], &str);
         }
      }
   return str.s;
 }
