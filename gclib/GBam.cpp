#include "GBam.h"
#include <ctype.h>
#include "kstring.h"

//auxiliary functions for low level BAM record creation
uint8_t* realloc_bdata(bam1_t *b, int size) {
  if (b->m_data < size) {
        b->m_data = size;
        kroundup32(b->m_data);
        b->data = (uint8_t*)realloc(b->data, b->m_data);
        }
  if (b->data_len<size) b->data_len=size;
  return b->data;
}

uint8_t* dupalloc_bdata(bam1_t *b, int size) {
  //same as realloc_bdata, but does not free previous data
  //but returns it instead
  //it ALWAYS duplicates data
  b->m_data = size;
  kroundup32(b->m_data);
  uint8_t* odata=b->data;
  b->data = (uint8_t*)malloc(b->m_data);
  memcpy((void*)b->data, (void*)odata, b->data_len);
  b->data_len=size;
  return odata; //user must FREE this after
}

GBamRecord::GBamRecord(const char* qname, int32_t gseq_tid,
                 int pos, bool reverse, const char* qseq,
                 const char* cigar, const char* quals):iflags(0), exons(1),
                		 clipL(0), clipR(0), mapped_len(0) {
   novel=true;
   bam_header=NULL;
   b=bam_init1();
   if (pos<=0 || gseq_tid<0) {
               b->core.pos=-1; //unmapped
               b->core.flag |= BAM_FUNMAP;
               gseq_tid=-1;
               }
          else b->core.pos=pos-1; //BAM is 0-based
   b->core.tid=gseq_tid;
   b->core.qual=255;
   b->core.mtid=-1;
   b->core.mpos=-1;
   int l_qseq=strlen(qseq);
   //this may not be accurate, setting CIGAR is the correct way
   //b->core.bin = bam_reg2bin(b->core.pos, b->core.pos+l_qseq-1);
   b->core.l_qname=strlen(qname)+1; //includes the \0 at the end
   memcpy(realloc_bdata(b, b->core.l_qname), qname, b->core.l_qname);
   set_cigar(cigar); //this will also set core.bin
   add_sequence(qseq, l_qseq);
   add_quals(quals); //quals must be given as Phred33
   if (reverse) { b->core.flag |= BAM_FREVERSE ; }
   }

GBamRecord::GBamRecord(const char* qname, int32_t samflags, int32_t g_tid,
             int pos, int map_qual, const char* cigar, int32_t mg_tid, int mate_pos,
             int insert_size, const char* qseq, const char* quals,
             GVec<char*>* aux_strings):iflags(0), exons(1)  {
  novel=true;
  bam_header=NULL;
  b=bam_init1();
  b->core.tid=g_tid;
  b->core.pos = (pos<=0) ? -1 : pos-1; //BAM is 0-based
  b->core.qual=map_qual;
  int l_qseq=strlen(qseq);
  b->core.l_qname=strlen(qname)+1; //includes the \0 at the end
  memcpy(realloc_bdata(b, b->core.l_qname), qname, b->core.l_qname);
  set_cigar(cigar); //this will also set core.bin
  add_sequence(qseq, l_qseq);
  add_quals(quals); //quals must be given as Phred33
  set_flags(samflags);
  set_mdata(mg_tid, (int32_t)(mate_pos-1), (int32_t)insert_size);
  if (aux_strings!=NULL) {
    for (int i=0;i<aux_strings->Count();i++) {
       add_aux(aux_strings->Get(i));
       }
    }
}

 void GBamRecord::set_cigar(const char* cigar) {
   //requires b->core.pos and b->core.flag to have been set properly PRIOR to this call
   int doff=b->core.l_qname;
   uint8_t* after_cigar=NULL;
   int after_cigar_len=0;
   uint8_t* prev_bdata=NULL;
   if (b->data_len>doff) {
      //cigar string already allocated, replace it
      int d=b->core.l_qname + b->core.n_cigar * 4;//offset of after-cigar data
      after_cigar=b->data+d;
      after_cigar_len=b->data_len-d;
      }
   const char *s;
   char *t;
   int i, op;
   long x;
   b->core.n_cigar = 0;
   if (cigar != NULL && strcmp(cigar, "*") != 0) {
        for (s = cigar; *s; ++s) {
            if (isalpha(*s)) b->core.n_cigar++;
            else if (!isdigit(*s)) {
                 GError("Error: invalid CIGAR character (%s)\n",cigar);
                 }
            }
        if (after_cigar_len>0) { //replace/insert into existing full data
             prev_bdata=dupalloc_bdata(b, doff + b->core.n_cigar * 4 + after_cigar_len);
             memcpy((void*)(b->data+doff+b->core.n_cigar*4),(void*)after_cigar, after_cigar_len);
             free(prev_bdata);
             }
           else {
             realloc_bdata(b, doff + b->core.n_cigar * 4);
             }
        for (i = 0, s = cigar; i != b->core.n_cigar; ++i) {
            x = strtol(s, &t, 10);
            op = toupper(*t);
            if (op == 'M' || op == '=' || op == 'X') op = BAM_CMATCH;
            else if (op == 'I') op = BAM_CINS;
            else if (op == 'D') op = BAM_CDEL;
            else if (op == 'N') op = BAM_CREF_SKIP; //has_Introns=true;
            else if (op == 'S') op = BAM_CSOFT_CLIP; //soft_Clipped=true;
            else if (op == 'H') op = BAM_CHARD_CLIP; //hard_Clipped=true;
            else if (op == 'P') op = BAM_CPAD;
            else GError("Error: invalid CIGAR operation (%s)\n",cigar);
            s = t + 1;
            bam1_cigar(b)[i] = x << BAM_CIGAR_SHIFT | op;
        }
        if (*s) GError("Error: unmatched CIGAR operation (%s)\n",cigar);
        b->core.bin = bam_reg2bin(b->core.pos, bam_calend(&b->core, bam1_cigar(b)));
    } else {//no CIGAR string given
        if (!(b->core.flag&BAM_FUNMAP)) {
            GMessage("Warning: mapped sequence without CIGAR (%s)\n", (char*)b->data);
            b->core.flag |= BAM_FUNMAP;
        }
        b->core.bin = bam_reg2bin(b->core.pos, b->core.pos + 1);
    }
   setupCoordinates();
   } //set_cigar()

 void GBamRecord::add_sequence(const char* qseq, int slen) {
   //must be called AFTER set_cigar (cannot replace existing sequence for now)
   if (qseq==NULL) return; //should we ever care about this?
   if (slen<0) slen=strlen(qseq);
   int doff = b->core.l_qname + b->core.n_cigar * 4;
   if (strcmp(qseq, "*")!=0) {
       b->core.l_qseq=slen;
       if (b->core.n_cigar && b->core.l_qseq != (int32_t)bam_cigar2qlen(&b->core, bam1_cigar(b)))
           GError("Error: CIGAR and sequence length are inconsistent!(%s)\n",
                  qseq);
       uint8_t* p = (uint8_t*)realloc_bdata(b, doff + (b->core.l_qseq+1)/2 + b->core.l_qseq) + doff;
       //also allocated quals memory
       memset(p, 0, (b->core.l_qseq+1)/2);
       for (int i = 0; i < b->core.l_qseq; ++i)
           p[i/2] |= bam_nt16_table[(int)qseq[i]] << 4*(1-i%2);
       } else b->core.l_qseq = 0;
   }

 void GBamRecord::add_quals(const char* quals) {
   //requires core.l_qseq already set
   //and must be called AFTER add_sequence(), which also allocates the memory for quals
   uint8_t* p = b->data+(b->core.l_qname + b->core.n_cigar * 4 + (b->core.l_qseq+1)/2);
   if (quals==NULL || strcmp(quals, "*") == 0) {
      for (int i=0;i < b->core.l_qseq; i++) p[i] = 0xff;
      return;
      }
   for (int i=0;i < b->core.l_qseq; i++) p[i] = quals[i]-33;
   }

 void GBamRecord::add_aux(const char* str) {
     //requires: being called AFTER add_quals()
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
             } else parse_error("unrecognized aux type");
  this->add_aux(tag, atype, alen, adata);
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
 case BAM_CMATCH: // M
      //have to actually check for mismatches: num_mismatches+=count_mismatches;
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

void GBamRecord::setupCoordinates() {
	const bam1_core_t *c = &b->core;
	if (c->flag & BAM_FUNMAP) return; /* skip unmapped reads */
	uint32_t *p = bam1_cigar(b);
	//--- prevent alignment error here (reported by UB-sanitazer):
	uint32_t *cigar= new uint32_t[c->n_cigar];
	memcpy(cigar, p, c->n_cigar * sizeof(uint32_t));
	//--- UBsan protection end
	int l=0;
	mapped_len=0;
	clipL=0;
	clipR=0;
	start=c->pos+1; //genomic start coordinate, 1-based (BAM core.pos is 0-based)
	int exstart=c->pos;
	for (int i = 0; i < c->n_cigar; ++i) {
		int op = cigar[i]&0xf;
		if (op == BAM_CMATCH || op==BAM_CEQUAL ||
				op == BAM_CDIFF || op == BAM_CDEL) {
			l += cigar[i]>>4;
		}
		else if (op == BAM_CREF_SKIP) { //N
			//intron starts
			//exon ends here
			has_Introns=true;
			GSeg exon(exstart+1,c->pos+l);
			exons.Add(exon);
			mapped_len+=exon.len();
			l += cigar[i]>>4;
			exstart=c->pos+l;
		}
		else if (op == BAM_CSOFT_CLIP) {
			soft_Clipped=true;
			if (l) clipR=(cigar[i]>>4);
			else clipL=(cigar[i]>>4);
		}
		else if (op == BAM_CHARD_CLIP) {
			hard_Clipped=true;
		}
	}
	GSeg exon(exstart+1,c->pos+l);
	exons.Add(exon);
	mapped_len+=exon.len();
	end=c->pos+l; //genomic end coordinate
	delete[] cigar; //UBsan protection
}


 uint8_t* GBamRecord::find_tag(const char tag[2]) {
   return bam_aux_get(this->b, tag);
 }

 char GBamRecord::tag_char(const char tag[2]) { //retrieve tag data as single char
   uint8_t* s=find_tag(tag);
   if (s) return ( bam_aux2A(s) );
   return 0;
  }

 int GBamRecord::tag_int(const char tag[2]) { //get the numeric value of tag
   uint8_t *s=find_tag(tag);
   if (s) return ( bam_aux2i(s) );
   return 0;
   }

 float GBamRecord::tag_float(const char tag[2]) { //get the float value of tag
    uint8_t *s=bam_aux_get(this->b, tag);;
    if (s) return ( bam_aux2f(s) );
    return 0;
    }

 char* GBamRecord::tag_str(const char tag[2]) { //return string value for a tag
   uint8_t *s=find_tag(tag);
   if (s) return ( bam_aux2Z(s) );
   return NULL;
   }

 char GBamRecord::spliceStrand() { // '+', '-' from the XS tag, or 0 if no XS tag
   char c=tag_char("XS");
   if (c==0) {
    //try minimap2's "ts" tag
    char m=tag_char("ts");
    if (m=='+' || m=='-') {
       if ((this->b->core.flag & BAM_FREVERSE) != 0) c=((m=='+') ? '-' : '+');
         else c=m;
    }
   }
   return ((c=='+' || c=='-') ? c : '.');
 }

 char* GBamRecord::sequence() { //user must free this after use
   char *s = (char*)bam1_seq(b);
   char* qseq=NULL;
   GMALLOC(qseq, (b->core.l_qseq+1));
   int i;
   for (i=0;i<(b->core.l_qseq);i++) {
     int8_t v = bam1_seqi(s,i);
     qseq[i] = bam_nt16_rev_table[v];
     }
   qseq[i] = 0;
   return qseq;
   }

 char* GBamRecord::qualities() {//user must free this after use
   char *qual  = (char*)bam1_qual(b);
   char* qv=NULL;
   GMALLOC(qv, (b->core.l_qseq+1) );
   int i;
   for(i=0;i<(b->core.l_qseq);i++) {
     qv[i]=qual[i]+33;
     }
   qv[i]=0;
   return qv;
   }

 char* GBamRecord::cigar() { //returns text version of the CIGAR string; must be freed by user
   kstring_t str;
   str.l = str.m = 0; str.s = 0;
   if (b->core.n_cigar == 0) kputc('*', &str);
    else {
      for (int i = 0; i < b->core.n_cigar; ++i) {
         kputw(bam1_cigar(b)[i]>>BAM_CIGAR_SHIFT, &str);
         kputc("MIDNSHP=X"[bam1_cigar(b)[i]&BAM_CIGAR_MASK], &str);
         }
      }
   return str.s;
   }
