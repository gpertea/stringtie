#include "GSam.h"
#include <ctype.h>
/*
//for bam1_t (re)allocation functions:
// sam_realloc_bam_data(), realloc_bam_data(), possibly_expand_bam_data()
#include "sam_internal.h"

//for parsing functions hts_str2uint() etc.:
#include "textutils_internal.h"

#define  _get_bmem(type_t, x, b, l) if (possibly_expand_bam_data((b), (l)) < 0) \
 GError("Error: cannot allocate SAM data\n"); \
 *(x) = (type_t*)((b)->data + (b)->l_data); (b)->l_data += (l)
*/
#define _parse_err(cond, msg) if (cond) GError("Error [SAM parsing]: %s\n",msg);
#define _parse_warn(cond, msg) if (cond) GMessage("Warning [SAM parsing]: %s\n",msg);
#define _parse_mem_err() GError("Error [SAM parsing]: memory allocation problem!\n");

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
/*
GSamRecord::GSamRecord(const char* qname, int32_t samflags, int32_t g_tid,
             int pos, int map_qual, const char* cigar, int32_t mg_tid, int mate_pos,
             int insert_size, const char* qseq, const char* quals,
             GVec<char*>* aux_strings):iflags(0), exons(1)  {
  novel=true;
  b_hdr=NULL;
  b=bam_init1();
  b->core.tid=g_tid;
  b->core.pos = (pos<=0) ? -1 : pos-1; //BAM is 0-based
  b->core.qual=map_qual;

  int l_qname=strlen(qname);
  //from sam_parse1() in sam.c:
  _parse_warn(l_qname <= 1, "empty query name");
  _parse_err(l_qname > 255, "query name too long");
  // resize large enough for name + extranul
  if (possibly_expand_bam_data(b, l_qname + 4) < 0)
     _parse_mem_err();
  memcpy(b->data + b->l_data, qname, l_qname);
  b->l_data += l_qname;

  b->core.l_extranul = (4 - (b->l_data & 3)) & 3;
  memcpy(b->data + b->l_data, "\0\0\0\0", b->core.l_extranul);
  b->l_data += b->core.l_extranul;

  b->core.l_qname = l_qname + b->core.l_extranul;

  set_cigar(cigar); //this will also set core.bin
  if (qseq!=NULL) add_sequence(qseq, strlen(qseq));
  if (quals!=NULL) add_quals(quals); //quals must be given as Phred33
  set_flags(samflags);
  set_mdata(mg_tid, (int32_t)(mate_pos-1), (int32_t)insert_size);
  if (aux_strings!=NULL) {
    for (int i=0;i<aux_strings->Count();i++) {
       add_aux(aux_strings->Get(i));
       }
    }
}
*/

/*
void GSamRecord::set_cigar(const char* str) {
   //requires b->core.pos and b->core.flag to have been set properly PRIOR to this call
  // also expects the b record memory to not be allocated already (fresh record creation)
  if (b==NULL) GError("Error: invalid call to ::set_cigar() (b is NULL)\n");
  //SAM header ptr is in b_hdr
  char *p = const_cast<char*>(str);
  bam1_core_t *c = &b->core;
  int overflow = 0;
  hts_pos_t cigreflen;
  if (*p == '*') {
      _parse_warn(!(c->flag&BAM_FUNMAP),
    		  "mapped query must have a CIGAR; treated as unmapped");
      c->flag |= BAM_FUNMAP;
      cigreflen = 1;
  }
  else {
	  char* q;
	  uint i=0;
      uint32_t *cigar;
      size_t n_cigar = 0;
      for (q = p; *p!=0; ++p)
          if (!isdigit(*p)) ++n_cigar;
      _parse_err( n_cigar == 0, "no CIGAR operations");
      _parse_err(n_cigar >= 2147483647, "too many CIGAR operations");
      c->n_cigar = n_cigar;
      _get_bmem(uint32_t, &cigar, b, c->n_cigar * sizeof(uint32_t));
      for (i = 0; i < c->n_cigar; ++i) {
          int op;
          cigar[i] = hts_str2uint(q, &q, 28, &overflow)<<BAM_CIGAR_SHIFT;
          op = bam_cigar_table[(unsigned char)*q++];
          if (op<0) parse_error("unrecognized CIGAR operator");
          cigar[i] |= op;
      }
      // can't use bam_endpos() directly as some fields not yet set up
      cigreflen = (!(c->flag&BAM_FUNMAP))? bam_cigar2rlen(c->n_cigar, cigar) : 1;
  }
  _parse_err(HTS_POS_MAX - cigreflen <= c->pos,
             "read ends beyond highest supported position");
  c->bin = hts_reg2bin(c->pos, c->pos + cigreflen, 14, 5);

}


 void GSamRecord::add_sequence(const char* qseq, int slen) {
    // ---- see sam_parse1() in htslib/sam.c for details
    //must be called AFTER set_cigar (cannot replace existing sequence for now)
   if (qseq==NULL) return; //should we ever care about this?
   if (slen<0) slen=strlen(qseq);
   if (strcmp(qseq, "*")!=0) {
       b->core.l_qseq=slen;
       hts_pos_t ql = bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b));
      _parse_err(b->core.n_cigar && ql != b->core.l_qseq,
    		  "CIGAR and query sequence are of different length");

       if (b->core.n_cigar && b->core.l_qseq != ql)
           GError("Error: CIGAR and sequence length are inconsistent!(%s:%s)\n",
        		   bam_get_qname(b), qseq);
       int v = (b->core.l_qseq + 1) >> 1;
       uint8_t *t;
       _get_bmem(uint8_t, &t, b, v);

       unsigned int lqs2 = b->core.l_qseq&~1, i;
       for (i = 0; i < lqs2; i+=2)
           t[i>>1] = (seq_nt16_table[(unsigned char)qseq[i]] << 4) | seq_nt16_table[(unsigned char)qseq[i+1]];
       for (; i < (unsigned int)b->core.l_qseq; ++i)
           t[i>>1] = seq_nt16_table[(unsigned char)qseq[i]] << ((~i&1)<<2);
   }
     else b->core.l_qseq = 0;
 }

 void GSamRecord::add_quals(const char* quals) {
   //must be called immediately AFTER add_sequence()
   uint8_t *t;
   //this will just append newly allocated mem to b->data:
   _get_bmem(uint8_t, &t, b, b->core.l_qseq);
   if (quals==NULL || strcmp(quals, "*") == 0) {
      memset(t, 0xff, b->core.l_qseq);
      return;
   }
   int qlen=strlen(quals);
   _parse_err(qlen!=b->core.l_qseq, "SEQ and QUAL are of different length!" );
   //uint8_t* p=bam_get_qual(b);
   for (int i=0;i < b->core.l_qseq; i++) t[i] = quals[i]-33;
 }
*/
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
		    l+=_cigLen(cigar[i]);
			if(intron) { // op comes after intron --> update juncdel
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
		    // exon ends here
		    if(!ins || !intron) { // insertion in the middle of an intron --> adjust last exon
		      exon.end=c->pos+l;
		      exon.start=exstart+1;
		      exons.Add( exon );
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
		    intron=false; ins=false;
		    del=0;
		    break;
		  case BAM_CHARD_CLIP:
		    hard_Clipped=true;
		    intron=false; ins=false;
		    del=0;
		    break;
		  case BAM_CPAD:
		    //gpos+=cl;
		    break;
		  default:
		    int cl=_cigLen(cigar[i]);
		    fprintf(stderr, "Unhandled CIGAR operation %d:%d\n", op, cl);
		}
	}
	exon.start=exstart+1;
	exon.end=c->pos+l;
	exons.Add(exon);
	mapped_len+=exon.len();
	end=exon.end; //genomic end coordinate
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
