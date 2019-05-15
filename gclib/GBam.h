#ifndef _G_BAM_H
#define _G_BAM_H
#include "GBase.h"
#include "GList.hh"
#include "bam.h"
#include "sam.h"

class GBamReader;
class GBamWriter;

class GBamRecord: public GSeg {
   friend class GBamReader;
   friend class GBamWriter;
   bam1_t* b;
   // b->data has the following strings concatenated:
   //  qname (including the terminal \0)
   //  +cigar (each event encoded on 32 bits)
   //   +seq  (4bit-encoded)
   //    +qual
   //     +aux

   union {
	  uint16_t iflags;
      struct {
    	  bool novel         :1; //if set, the destructor must free b
    	  bool hard_Clipped  :1;
    	  bool soft_Clipped  :1;
    	  bool has_Introns   :1;
      };
   };
   bam_header_t* bam_header;
   char tag[2];
   uint8_t abuf[512];
 public:
   GVec<GSeg> exons; //coordinates will be 1-based
   int clipL; //soft clipping data, as seen in the CIGAR string
   int clipR;
   int mapped_len; //sum of exon lengths
   bool isHardClipped() { return hard_Clipped; }
   bool isSoftClipped() { return soft_Clipped; }
   bool hasIntrons() { return has_Introns; }
   //created from a reader:
   void bfree_on_delete(bool b_free=true) { novel=b_free; }
   GBamRecord(bam1_t* from_b=NULL, bam_header_t* b_header=NULL, bool b_free=true):iflags(0), exons(1),
		   clipL(0), clipR(0), mapped_len(0) {
      bam_header=NULL;
      if (from_b==NULL) {
           b=bam_init1();
           novel=true;
      }
      else {
           b=from_b; //it'll take over from_b
           novel=b_free;
      }

      bam_header=b_header;
      setupCoordinates();//set 1-based coordinates (start, end and exons array)
   }

   GBamRecord(GBamRecord& r):GSeg(r.start, r.end), iflags(0), exons(r.exons),
		   clipL(r.clipL), clipR(r.clipR), mapped_len(r.mapped_len) { //copy constructor
	      //makes a new copy of the bam1_t record etc.
	      clear();
	      b=bam_dup1(r.b);
	      novel=true; //will also free b when destroyed
   }

   const GBamRecord& operator=(GBamRecord& r) {
	  //copy operator
      //makes a new copy of the bam1_t record etc.
      clear();
      b=bam_dup1(r.b);
      iflags=r.iflags;
      novel=true; //will also free b when destroyed
      start=r.start;
      end=r.end;
      exons = r.exons;
      clipL = r.clipL;
      clipR = r.clipR;
      mapped_len=r.mapped_len;
      return *this;
      }

     void setupCoordinates();
     void clear() {
        if (novel) {
           bam_destroy1(b);
           //novel=false;
        }
        b=NULL;
        exons.Clear();
        mapped_len=0;
        bam_header=NULL;
        iflags=0;
    }

    ~GBamRecord() {
       clear();
    }

    void parse_error(const char* s) {
      GError("BAM parsing error: %s\n", s);
      }

    bam1_t* get_b() { return b; }

    void set_mdata(int32_t mtid, int32_t m0pos, //0-based coordinate, -1 if not available
                     int32_t isize=0) { //mate info for current record
      b->core.mtid=mtid;
      b->core.mpos=m0pos; // should be -1 if '*'
      b->core.isize=isize; //should be 0 if not available
    }

    void set_flags(uint16_t samflags) {
      b->core.flag=samflags;
    }

    //creates a new record from 1-based alignment coordinate
    //quals should be given as Phred33
    //Warning: pos and mate_pos must be given 1-based!
    GBamRecord(const char* qname, int32_t gseq_tid,
                    int pos, bool reverse, const char* qseq, const char* cigar=NULL, const char* quals=NULL);
    GBamRecord(const char* qname, int32_t samflags, int32_t g_tid,
             int pos, int map_qual, const char* cigar, int32_t mg_tid, int mate_pos,
             int insert_size, const char* qseq, const char* quals=NULL,
             GVec<char*>* aux_strings=NULL);
             //const std::vector<std::string>* aux_strings=NULL);
    void set_cigar(const char* cigar); //converts and adds CIGAR string given in plain SAM text format
    void add_sequence(const char* qseq, int slen=-1); //adds the DNA sequence given in plain text format
    void add_quals(const char* quals); //quality values string in Phred33 format
    void add_aux(const char* str); //adds one aux field in plain SAM text format (e.g. "NM:i:1")
    void add_aux(const char tag[2], char atype, int len, uint8_t *data) {
      //IMPORTANT:  strings (Z,H) should include the terminal \0
      int addz=0;
      if ((atype=='Z' || atype=='H') && data[len-1]!=0) {
        addz=1;
        }
      int ori_len = b->data_len;
      b->data_len += 3 + len + addz;
      b->l_aux += 3 + len + addz;
      if (b->m_data < b->data_len) {
        b->m_data = b->data_len;
        kroundup32(b->m_data);
        b->data = (uint8_t*)realloc(b->data, b->m_data);
      }
      b->data[ori_len] = tag[0]; b->data[ori_len + 1] = tag[1];
      b->data[ori_len + 2] = atype;
      if (addz) {
        b->data[ori_len+len+4]=0;
        }
      memcpy(b->data + ori_len + 3, data, len);
      }

    void add_tag(const char tag[2], char atype, int len, uint8_t *data) {
      //same with add_aux()
      add_aux(tag,atype,len,data);
      }
 //--query methods:
 uint32_t flags() { return b->core.flag; } //return SAM flags
 bool isUnmapped() { return ((b->core.flag & BAM_FUNMAP) != 0); }
 bool isMapped() { return ((b->core.flag & BAM_FUNMAP) == 0); }
 bool isPaired() { return ((b->core.flag & BAM_FPAIRED) != 0); }
 const char* name() { return bam1_qname(b); }
 int pairOrder() {
    //which read in the pair: 0 = unpaired, 1=first read, 2=second read
    int r=0;
    if ((b->core.flag & BAM_FREAD1) != 0) r=1;
    else if ((b->core.flag & BAM_FREAD2) != 0) r=2;
    return r;
    }
 bool revStrand() {
   //this is the raw alignment strand, NOT the transcription/splice strand
   return ((b->core.flag & BAM_FREVERSE) != 0);
 }
 const char* refName() {
   return (bam_header!=NULL) ?
         ((b->core.tid<0) ? "*" : bam_header->target_name[b->core.tid]) : NULL;
   }
 int32_t refId() { return b->core.tid; }
 int32_t mate_refId() { return b->core.mtid; }
 const char* mate_refName() {
    return (bam_header!=NULL) ?
       ((b->core.mtid<0) ? "*" : bam_header->target_name[b->core.mtid]) : NULL;
    }
 int32_t insertSize() { return b->core.isize; }
 int32_t mate_start() { return b->core.mpos<0? 0 : b->core.mpos+1; }

 //int find_tag(const char tag[2], uint8_t* & s, char& tag_type);
 uint8_t* find_tag(const char tag[2]);
 //position s at the beginning of tag data, tag_type is set to the found tag type
 //returns length of tag data, or 0 if tag not found

 char* tag_str(const char tag[2]); //return tag value for tag type 'Z'
 int tag_int(const char tag[2]); //return numeric value of tag (for numeric types)
 float tag_float(const char tag[2]); //return float value of tag (for float types)
 char tag_char(const char tag[2]); //return char value of tag (for type 'A')
 char spliceStrand(); // '+', '-' from the XS tag, or '.' if no XS tag

 char* sequence(); //user should free after use
 char* qualities();//user should free after use
 char* cigar(); //returns text version of the CIGAR string; user must free
};

// from sam.c:
#define FTYPE_BAM  1
#define FTYPE_READ 2

class GBamReader {
   samfile_t* bam_file;
   char* fname;
   // from bam_import.c:
   struct samtools_tamFile_t {
   	gzFile fp;
   	void *ks;
   	void *str;
   	uint64_t n_lines;
   	int is_first;
   };

 public:
   void bopen(const char* filename, bool forceBAM=false) {
      if (strcmp(filename, "-")==0) {
        //if stdin was given, we assume it's text SAM, unless forceBAM was given
        if (forceBAM) bam_file=samopen(filename, "rb", 0);
        else bam_file=samopen(filename, "r", 0);
        }
      else {
        FILE* f=Gfopen(filename);
        if (f==NULL) {
           GError("Error opening SAM/BAM file %s!\n", filename);
        }
        if (forceBAM) {
           //directed to open this as a BAM file
            if (forceBAM) bam_file=samopen(filename, "rb", 0);
        }
        else {
          //try to guess if it's BAM or SAM
          //BAM files have the zlib signature bytes at the beginning: 1F 8B 08
          //if that's not present then we assume text SAM
          byte fsig[3];
          size_t rd=fread(fsig, 1, 3, f);
          fclose(f);
          if (rd<3) GError("Error reading from file %s!\n",filename);
          if ((fsig[0]==0x1F && fsig[1]==0x8B && fsig[2]==0x08) ||
            (fsig[0]=='B' && fsig[1]=='A' && fsig[2]=='M')) {
            bam_file=samopen(filename, "rb", 0); //BAM or uncompressed BAM
          }
          else { //assume text SAM file
            bam_file=samopen(filename, "r", 0);
          }
        }
      }
      if (bam_file==NULL)
         GError("Error: could not open SAM file %s!\n",filename);
      fname=Gstrdup(filename);
   }

   GBamReader(const char* fn, bool forceBAM=false) {
      bam_file=NULL;
      fname=NULL;
      bopen(fn, forceBAM);
      }

   bam_header_t* header() {
      return bam_file? bam_file->header : NULL;
      }
   void bclose() {
      if (bam_file) {
        samclose(bam_file);
        bam_file=NULL;
        }
      }
   ~GBamReader() {
      bclose();
      GFREE(fname);
      }
   int64_t fpos() {
  	 if ( bam_file->type & FTYPE_BAM )
  	   return bgzf_tell(bam_file->x.bam);
  	 else
  		 return (int64_t)gztell(((samtools_tamFile_t*)(bam_file->x.tamr))->fp);
   }
   int64_t fseek(int64_t offs) {
  	 if ( bam_file->type & FTYPE_BAM )
  		 return bgzf_seek(bam_file->x.bam, offs, SEEK_SET);
  	 else
  		 return (int64_t)gzseek(((samtools_tamFile_t*)(bam_file->x.tamr))->fp, offs, SEEK_SET);
   }
   void rewind() {
     if (fname==NULL) {
       GMessage("Warning: GBamReader::rewind() called without a file name.\n");
       return;
       }
     bclose();
     char* ifname=fname;
     bopen(ifname);
     GFREE(ifname);
     }

   GBamRecord* next() {
      if (bam_file==NULL)
        GError("Warning: GBamReader::next() called with no open file.\n");
      bam1_t* b = bam_init1();
      if (samread(bam_file, b) >= 0) {
        GBamRecord* bamrec=new GBamRecord(b, bam_file->header, true);
        return bamrec;
        }
      else {
        bam_destroy1(b);
        return NULL;
        }
      }
};


class GBamWriter {
   samfile_t* bam_file;
   bam_header_t* bam_header;
 public:
   void create(const char* fname, bool uncompressed=false) {
      if (bam_header==NULL)
         GError("Error: no bam_header for GBamWriter::create()!\n");
      if (uncompressed) {
         bam_file=samopen(fname, "wbu", bam_header);
         }
        else {
         bam_file=samopen(fname, "wb", bam_header);
         }
      if (bam_file==NULL)
         GError("Error: could not create BAM file %s!\n",fname);
      }
   void create(const char* fname, bam_header_t* bh, bool uncompressed=false) {
     bam_header=bh;
     create(fname,uncompressed);
     }

   GBamWriter(const char* fname, bam_header_t* bh, bool uncompressed=false) {
      create(fname, bh, uncompressed);
      }

   GBamWriter(const char* fname, const char* samfname, bool uncompressed=false) {
      tamFile samf_in=sam_open(samfname);
      if (samf_in==NULL)
         GError("Error: could not open SAM file %s\n", samfname);
      bam_header=sam_header_read(samf_in);
      if (bam_header==NULL)
         GError("Error: could not read SAM header from %s!\n",samfname);
      sam_close(samf_in);
      create(fname, uncompressed);
      }

    ~GBamWriter() {
      samclose(bam_file);
      bam_header_destroy(bam_header);
      }
   bam_header_t* get_header() { return bam_header; }
   int32_t get_tid(const char *seq_name) {
      if (bam_header==NULL)
         GError("Error: missing SAM header (get_tid())\n");
      return bam_get_tid(bam_header, seq_name);
      }

   //just a convenience function for creating a new record, but it's NOT written
   //given pos must be 1-based (so it'll be stored as pos-1 because BAM is 0-based)
   GBamRecord* new_record(const char* qname, const char* gseqname,
            int pos, bool reverse, const char* qseq, const char* cigar=NULL, const char* qual=NULL) {
      int32_t gseq_tid=get_tid(gseqname);
      if (gseq_tid < 0 && strcmp(gseqname, "*")) {
            if (bam_header->n_targets == 0) {
               GError("Error: missing/invalid SAM header\n");
               } else
                   GMessage("Warning: reference '%s' not found in header, will consider it '*'.\n",
                                   gseqname);
            }

      return (new GBamRecord(qname, gseq_tid, pos, reverse, qseq, cigar, qual));
      }

   GBamRecord* new_record(const char* qname, int32_t samflags, const char* gseqname,
         int pos, int map_qual, const char* cigar, const char* mgseqname, int mate_pos,
         int insert_size, const char* qseq, const char* quals=NULL,
                          GVec<char*>* aux_strings=NULL) {
      int32_t gseq_tid=get_tid(gseqname);
      if (gseq_tid < 0 && strcmp(gseqname, "*")) {
            if (bam_header->n_targets == 0) {
               GError("Error: missing/invalid SAM header\n");
               } else
                   GMessage("Warning: reference '%s' not found in header, will consider it '*'.\n",
                                   gseqname);
            }
      int32_t mgseq_tid=-1;
      if (mgseqname!=NULL) {
         if (strcmp(mgseqname, "=")==0) {
            mgseq_tid=gseq_tid;
            }
          else {
            mgseq_tid=get_tid(mgseqname);
            if (mgseq_tid < 0 && strcmp(mgseqname, "*")) {
                GMessage("Warning: reference '%s' not found in header, will consider it '*'.\n",
                                   mgseqname);
                }
            }
          }
      return (new GBamRecord(qname, samflags, gseq_tid, pos, map_qual, cigar,
              mgseq_tid, mate_pos, insert_size, qseq, quals, aux_strings));
      }

   void write(GBamRecord* brec) {
      if (brec!=NULL)
          samwrite(this->bam_file,brec->get_b());
      }
   void write(bam1_t* b) {
      samwrite(this->bam_file, b);
      }
};

#endif
