//#define GFF_DEBUG 1
#include "rlink.h"
#ifndef NOTHREADS
#include "GThreads.h"
#endif

//#define GMEMTRACE 1  //debugging mem allocation

#ifdef GMEMTRACE
#include "proc_mem.h"
#endif

#define VERSION "1.1.0"

//uncomment this to show DBGPRINT messages (for threads)
//#define DEBUGPRINT 1

#ifdef DEBUGPRINT
#define DBGPRINT(x) GMessage(x)
#define DBGPRINT2(a,b) GMessage(a,b)
#define DBGPRINT3(a,b,c) GMessage(a,b,c)
#define DBGPRINT4(a,b,c,d) GMessage(a,b,c,d)
#define DBGPRINT5(a,b,c,d,e) GMessage(a,b,c,d,e)
#else
#define DBGPRINT(x) 
#define DBGPRINT2(a,b) 
#define DBGPRINT3(a,b,c) 
#define DBGPRINT4(a,b,c,d) 
#define DBGPRINT5(a,b,c,d,e) 
#endif

#define USAGE "StringTie v"VERSION" usage:\n\
 stringtie <input.bam> [-G <guide_gff>] [-l <label>] [-o <out_gtf>] [-p <cpus>]\n\
  [-v] [-a <min_anchor_len>] [-m <min_tlen>] [-j <min_anchor_cov>] \n\
  [-C <coverage_file_name>] [-c <min_bundle_cov>] [-g <bdist>]\n\
  [-e] [-x <seqid,..>] [-A <gene_file_name>] {-B | -b <dir_path>} \n\
\nAssemble RNA-Seq alignments into potential transcripts.\n\
 \n\
 Options:\n\
 -G reference annotation to use for guiding the assembly process (GTF/GFF3)\n\
 -l name prefix for output transcripts (default: STRG)\n\
 -f minimum isoform fraction (default: 0.1)\n\
 -m minimum assembled transcript length (default: 200)\n\
 -o output path/file name for the assembled transcripts GTF (default: stdout)\n\
 -a minimum anchor length for junctions (default: 10)\n\
 -j minimum junction coverage (default: 1)\n\
 -t disable trimming of predicted transcripts based on coverage\n\
    (default: coverage trimming is enabled)\n\
 -c minimum reads per bp coverage to consider for transcript assembly (default: 2.5)\n\
 -v verbose (log bundle processing details)\n\
 -g gap between read mappings triggering a new bundle (default: 50)\n\
 -C output file with reference transcripts that are covered by reads\n\
 -M fraction of bundle allowed to be covered by multi-hit reads (default:0.95)\n\
 -p number of threads (CPUs) to use (default: 1)\n\
 -B enable output of Ballgown table files which will be created in the\n\
    same directory as the output GTF (requires -G, -o recommended)\n\
 -b enable output of Ballgown table files but these files will be \n\
    created under the directory path given as <dir_path>\n\
 -e only estimates the abundance of given reference transcripts (requires -G)\n\
 -A gene abundance estimation\n\
 -x do not assemble any transcripts on these reference sequence(s)\n\
 "
/* 
 -n sensitivity level: 0,1, or 2, 3, with 3 the most sensitive level (default 0)\n\ \\ deprecated for now
 -O disable the coverage saturation limit and use a slower two-pass approach\n\
    to process the input alignments, collapsing redundant reads\n\
 -F disable fast computing for transcript path; default: yes\n\
 -i the reference annotation contains partial transcripts\n\
 -w weight the maximum flow algorithm towards the transcript with higher rate (abundance); default: no\n\
 -y include EM algorithm in max flow estimation; default: no\n\
 -z don't include source in the max flow algorithm\n\
 -P output file with all transcripts in reference that are partially covered by reads
 -M fraction of bundle allowed to be covered by multi-hit reads (paper uses default: 1)\n\
 -c minimum bundle reads per bp coverage to consider for assembly (paper uses default: 3)\n\
 -S more sensitive run (default: no) disabled for now \n\
 -s coverage saturation threshold; further read alignments will be\n\
    ignored in a region where a local coverage depth of <maxcov> \n\
    is reached (default: 1,000,000);\n\ \\ deprecated
*/
//---- globals

FILE* f_out=NULL;
FILE* c_out=NULL;
//#define B_DEBUG 1
#ifdef B_DEBUG
 FILE* dbg_out=NULL;
#endif

GStr outfname;
GStr out_dir;
GStr tmpfname;
GStr genefname;
bool guided=false;
bool trim=true;
bool fast=true;
bool eonly=false; // parameter -e
bool specific=false;
bool complete=true;
bool geneabundance=false;
//bool partialcov=false;
int num_cpus=1;
int mintranscriptlen=200; // minimum length for a transcript to be printed
int sensitivitylevel=1;
uint junctionsupport=10; // anchor length for junction to be considered well supported <- consider shorter??
int junctionthr=1; // number of reads needed to support a particular junction
float readthr=2.5;     // read coverage per bundle bp to accept it; otherwise considered noise; paper uses 3
uint bundledist=50;  // reads at what distance should be considered part of separate bundles
float mcov=0.95; // fraction of bundle allowed to be covered by multi-hit reads paper uses 1

// different options of implementation reflected with the next three options
bool includesource=true;
bool EM=false;
bool weight=false;

float isofrac=0.1;
GStr label("STRG");
GStr ballgown_dir;

GStr guidegff;

bool debugMode=false;
bool verbose=false;
bool ballgown=false;

//int maxReadCov=1000000; //max local read coverage (changed with -s option)
//no more reads will be considered for a bundle if the local coverage exceeds this value
//(each exon is checked for this)

bool singlePass=true; //-O will set this to False

int GeneNo=0; //-- global "gene" counter
double Num_Fragments=0; //global fragment counter (aligned pairs)
double Frag_Len=0;
double Cov_Sum=0;
//bool firstPrint=true; //just for writing the GFF header before the first transcript is printed

/*
double SumReads=0;
double SumFrag=0;
double NumCov=0;
double NumFrag=0;
double NumReads=0;
double NumFrag3=0;
double SumFrag3=0;
*/

GffNames* gseqNames=NULL; //used as a dictionary for genomic sequence names and ids

int refseqCount=0;


#ifdef GMEMTRACE
 double maxMemRS=0;
 double maxMemVM=0;
 GStr maxMemBundle;
#endif


#ifndef NOTHREADS

GMutex dataMutex; //manage availability of data records ready to be loaded by main thread
GVec<int> dataClear; //indexes of data bundles cleared for loading by main thread (clear data pool)
GConditionVar haveBundles; //will notify a thread that a bundle was loaded in the ready queue
                           //(or that no more bundles are coming)
int bundleWork=1; // bit 0 set if bundles are still being prepared (BAM file not exhausted yet)
                  // bit 1 set if there are Bundles ready in the queue


//GFastMutex waitMutex;
GMutex waitMutex; // controls threadsWaiting (idle threads counter)

int threadsWaiting; // idle worker threads
GConditionVar haveThreads; //will notify the bundle loader when a thread
                          //is available to process the currently loaded bundle

GConditionVar haveClear; //will notify when bundle buf space available

GMutex queueMutex; //controls bundleQueue and bundles access

GFastMutex printMutex; //for writing the output to file

GFastMutex logMutex; //only when verbose - to avoid mangling the log output

GFastMutex bamReadingMutex;

GFastMutex countMutex;

#endif

GHash<int> excludeGseqs; //hash of chromosomes/contigs to exclude (e.g. chrM)

bool NoMoreBundles=false;
bool moreBundles(); //thread-safe retrieves NoMoreBundles
void noMoreBundles(); //sets NoMoreBundles to true
//--
GStr Process_Options(GArgs* args);
char* sprintTime();

void processBundle(BundleData* bundle);
//void processBundle1stPass(BundleData* bundle); //two-pass testing

#ifndef NOTHREADS

bool waitForThreads(); //wait for at least 1 worker thread to enter "ready" state

bool noThreadsWaiting();

void workerThread(GThreadData& td); // Thread function

//check if a worker thread popped the data queue:
bool queuePopped(GPVec<BundleData>& bundleQueue, int prevCount); 

//prepare the next free bundle for loading
int waitForData(BundleData* bundles);
#endif

int main(int argc, char * const argv[]) {

 // == Process arguments.
 GArgs args(argc, argv, 
   //"debug;help;fast;xhvntj:D:G:C:l:m:o:a:j:c:f:p:g:");
   "debug;help;exclude=yzwFShvtiex:n:j:s:D:G:C:l:m:o:a:j:c:f:p:g:P:M:Bb:A:");
 args.printError(USAGE, true);

 GStr bamfname=Process_Options(&args);
 // == Done argument processing.

 GVec<GRefData> refguides; // plain vector with transcripts for each chromosome

 //table indexes for Ballgown Raw Counts data (-B/-b option)
 GPVec<RC_TData> guides_RC_tdata(true); //raw count data for all guide transcripts
 GPVec<RC_Feature> guides_RC_exons(true); //raw count data for all guide exons
 GPVec<RC_Feature> guides_RC_introns(true);//raw count data for all guide introns

 GVec<int> alncounts(30,0); //keep track of the number of read alignments per chromosome [gseq_id]

#ifdef DEBUGPRINT
  verbose=true;
#endif

const char* ERR_BAM_SORT="\nError: the input alignment file is not sorted!\n";

 if(guided) { // read guiding transcripts from input gff file
	 if (verbose) {
		 printTime(stderr);
		 GMessage(" Loading reference annotation (guides)..\n");
	 }
   FILE* f=fopen(guidegff.chars(),"r");
   if (f==NULL) GError("Error: could not open reference annotation file (%s)!\n",
       guidegff.chars());
   //                transcripts_only    sort gffr->gfflst by loc?
   GffReader gffr(f,       true,                   true); //loading only recognizable transcript features
   gffr.showWarnings(verbose);
   //        keepAttrs    mergeCloseExons   noExonAttrs
   gffr.readAll(false,          true,        true);
   //the list of GffObj is in gffr.gflst, sorted by chromosome and start-end coordinates
   //collect them in other data structures, if it's kept for later call gffobj->isUsed(true)
   // (otherwise it'll be deallocated when gffr is destroyed due to going out of scope)
   refseqCount=gffr.gseqStats.Count();
   if (refseqCount==0 || gffr.gflst.Count()==0) {
	   GError("Error: could not read reference annotation transcripts from GTF/GFF %s - invalid file?\n",
			   guidegff.chars());
   }
   refguides.setCount(refseqCount); //maximum gseqid
   uint c_tid=0;
   uint c_exon_id=0;
   uint c_intron_id=0;
   GList<RC_Feature> uexons(true, false, true); //sorted, free items, unique
   GList<RC_Feature> uintrons(true, false, true);
   //assign unique transcript IDs based on the sorted order
   int last_refid=-1;
   bool skipGseq=false;
   for (int i=0;i<gffr.gflst.Count();i++) {
	   GffObj* m=gffr.gflst[i];
	   if (last_refid!=m->gseq_id) {
		   //chromosome switch
		   if (ballgown) { //prepare memory storage/tables for all guides on this chromosome
			   uexons.Clear();
			   uintrons.Clear();
		   }
		   last_refid=m->gseq_id;
		   skipGseq=excludeGseqs.hasKey(m->getGSeqName());
	   }
	   //sanity check: make sure there are no exonless "genes" or other
	   if (skipGseq) continue;
	   if (m->exons.Count()==0) {
		    if (verbose)
		    	GMessage("Warning: exonless GFF object %s found, added implicit exon %d-%d.\n",
		    			m->getID(),m->start, m->end);
		    m->exons.Add(new GffExon(m->start, m->end));
	   }
	   if (ballgown) {
		   RC_TData* tdata=new RC_TData(*m, ++c_tid);
		   m->uptr=tdata;
		   guides_RC_tdata.Add(tdata);
		   tdata->rc_addFeatures(c_exon_id, uexons, guides_RC_exons,
		          c_intron_id, uintrons, guides_RC_introns);
	   }
	   GRefData& grefdata = refguides[m->gseq_id];
	   grefdata.add(&gffr, m); //transcripts already sorted by location
   }
	 if (verbose) {
		 printTime(stderr);
		 GMessage(" %d reference transcripts loaded.\n", gffr.gflst.Count());
	 }
	fclose(f);
 }

#ifdef GFF_DEBUG
  for (int r=0;r<refguides.Count();++r) {
	  GRefData& grefdata = refguides[r];
      for (int k=0;k<grefdata.rnas.Count();++k) {
    	  GMessage("#transcript #%d : %s (%d exons)\n", k, grefdata.rnas[k]->getID(), grefdata.rnas[k]->exons.Count());
    	  grefdata.rnas[k]->printGff(stderr);
      }
  }
  GMessage("GFF Debug mode, exiting...\n");
  exit(0);
#endif

 // --- here we do the input processing
 gseqNames=GffObj::names; //might have been populated already by gff data
 gffnames_ref(gseqNames);  //initialize the names collection if not guided

 GBamReader bamreader(bamfname.chars());

 GHash<int> hashread;      //read_name:pos:hit_index => readlist index

 //set of annotation transcript for the current locus
 GList<GffObj>* guides=NULL; //list of transcripts on a specific chromosome

 int currentstart=0, currentend=0;
 int ng_start=0;
 int ng_end=-1;
 int ng=0;
 GStr lastref;
 bool no_ref_used=true;
 int lastref_id=-1; //last seen gseq_id
 // int ncluster=0; used it for debug purposes only

 //Ballgown files
 FILE* f_tdata=NULL;
 FILE* f_edata=NULL;
 FILE* f_idata=NULL;
 FILE* f_e2t=NULL;
 FILE* f_i2t=NULL;
if (ballgown)
 Ballgown_setupFiles(f_tdata, f_edata, f_idata, f_e2t, f_i2t);
#ifndef NOTHREADS
 GThread* threads=new GThread[num_cpus]; //bundle processing threads

 GPVec<BundleData> bundleQueue(false); //queue of loaded bundles

 BundleData* bundles=new BundleData[num_cpus+1];
 //bundles[0..num_cpus-1] are processed by threads, loading bundles[num_cpus] first

 dataClear.setCapacity(num_cpus+1);
 for (int b=0;b<num_cpus;b++) {
	 threads[b].kickStart(workerThread, (void*) &bundleQueue);
	 bundles[b+1].idx=b+1;
	 dataClear.Push(b);
   }
 BundleData* bundle = &(bundles[num_cpus]);
#else
 BundleData bundles[1];
 BundleData* bundle = &(bundles[0]);
#endif
 GBamRecord* brec=NULL;
 bool more_alns=true;
 int prev_pos=0;
 bool skipGseq=false;
 while (more_alns) {
	 bool chr_changed=false;
	 int pos=0;
	 const char* refseqName=NULL;
	 //char strand=0;
	 char xstrand=0;
	 int nh=1;
	 int hi=0;
	 int gseq_id=lastref_id;  //current chr id
	 bool new_bundle=false;
	 delete brec;
	 if ((brec=bamreader.next())!=NULL) {
		 if (brec->isUnmapped()) continue;
		 refseqName=brec->refName();
		 xstrand=brec->spliceStrand();
		 if (xstrand=='.' && brec->exons.Count()>1)
			 continue; //skip spliced alignments lacking XS tag (e.g. HISAT alignments)
		 if (refseqName==NULL) GError("Error: cannot retrieve target seq name from BAM record!\n");
		 pos=brec->start; //BAM is 0 based, but GBamRecord makes it 1-based
		 chr_changed=(lastref.is_empty() || lastref!=refseqName);
		 if (chr_changed) {
			 skipGseq=excludeGseqs.hasKey(refseqName);
			 gseq_id=gseqNames->gseqs.addName(refseqName);
			 if (guided) {
				 if (gseq_id>=refseqCount) {
					 if (verbose)
						 GMessage("WARNING: no reference transcripts found for genomic sequence \"%s\"! (different naming convention?)\n",
								refseqName);
				 }
				 else no_ref_used=false;
			 }

			 if (alncounts.Count()<=gseq_id) {
				 alncounts.Resize(gseq_id+1, 0);
			 }
			 else if (alncounts[gseq_id]>0) GError(ERR_BAM_SORT);
			 prev_pos=0;
		 }
		 if (pos<prev_pos) GError(ERR_BAM_SORT);
		 prev_pos=pos;
		 if (skipGseq) continue;
		 alncounts[gseq_id]++;
		 nh=brec->tag_int("NH");
		 if (nh==0) nh=1;
		 hi=brec->tag_int("HI");
		 if (!chr_changed && currentend>0 && pos>currentend+(int)bundledist)
			   new_bundle=true;
	 }
	 else { //no more alignments
		 more_alns=false;
		 new_bundle=true; //fake a new start (end of last bundle)
	 }

	 if (new_bundle || chr_changed) {
		 hashread.Clear();
		 if (bundle->readlist.Count()>0) { // process reads in previous bundle
			// (readthr, junctionthr, mintranscriptlen are globals)
			bundle->getReady(currentstart, currentend);
#ifndef NOTHREADS
			//push this in the bundle queue where it'll be picked up by the threads
			DBGPRINT2("##> Locking queueMutex to push loaded bundle into the queue (bundle.start=%d)\n", bundle->start);
			int qCount=0;
			queueMutex.lock();
			bundleQueue.Push(bundle);
			bundleWork |= 0x02; //set bit 1
			qCount=bundleQueue.Count();
			queueMutex.unlock();
			DBGPRINT2("##> bundleQueue.Count()=%d)\n", qCount);
			//wait for a thread to pop this bundle from the queue
			waitMutex.lock();
			DBGPRINT("##> waiting for a thread to become available..\n");
			while (threadsWaiting==0) {
				haveThreads.wait(waitMutex);
			}
			waitMutex.unlock();
			haveBundles.notify_one();
			this_thread::yield();
			queueMutex.lock();
			while (bundleQueue.Count()==qCount) {
				queueMutex.unlock();
				haveBundles.notify_one();
				this_thread::yield();
				queueMutex.lock();
			}
			queueMutex.unlock();

#else //no threads
			//Num_Fragments+=bundle->num_fragments;
			//Frag_Len+=bundle->frag_len;
			processBundle(bundle);
#endif
			// ncluster++; used it for debug purposes only
		 } //have alignments to process
		 else { //no read alignments in this bundle?
			bundle->Clear();
#ifndef NOTHREADS
			dataMutex.lock();
			dataClear.Push(bundle->idx);
			dataMutex.unlock();
#endif
		 } //nothing to do with this bundle

		 if (chr_changed) {
			 if (guided) {
				 ng=0;
				 guides=NULL;
				 ng_start=0;
				 ng_end=-1;
				 if (refguides.Count()>gseq_id && refguides[gseq_id].rnas.Count()>0) {
					 guides=&(refguides[gseq_id].rnas);
					 ng=guides->Count();
				 }
			 }
			 lastref=refseqName;
			 lastref_id=gseq_id;
			 currentend=0;
		 }

		 if (!more_alns) {
				if (verbose) {
#ifndef NOTHREADS
					GLockGuard<GFastMutex> lock(logMutex);
#endif
					printTime(stderr);
					if(Num_Fragments) GMessage(" %g aligned fragments found.\n", Num_Fragments);
					//GMessage(" Done reading alignments.\n");
				}
			 noMoreBundles();
			 break;
		 }
#ifndef NOTHREADS

		 int new_bidx=waitForData(bundles);
		 if (new_bidx<0) {
			 //should never happen!
			 GError("Error: waitForData() returned invalid bundle index(%d)!\n",new_bidx);
			 break;
		 }
		 bundle=&(bundles[new_bidx]);
#endif
		 currentstart=pos;
		 currentend=brec->end;
		 if (guides) { //guided and guides!=NULL
			 ng_start=ng_end+1;
			 while (ng_start<ng && (int)(*guides)[ng_start]->end < pos) {
				 // skip guides that have no read coverage
				 ng_start++;
			 }
			 int ng_ovlstart=ng_start;
			 //add all guides overlapping the current read and other guides that overlap them
			 while (ng_ovlstart<ng && (int)(*guides)[ng_ovlstart]->start<=currentend) {
				 if (currentstart>(int)(*guides)[ng_ovlstart]->start)
					 currentstart=(*guides)[ng_ovlstart]->start;
				 if (currentend<(int)(*guides)[ng_ovlstart]->end)
					 currentend=(*guides)[ng_ovlstart]->end;
				  bundle->keepGuide((*guides)[ng_ovlstart],
						   &guides_RC_tdata, &guides_RC_exons, &guides_RC_introns);
				 ng_ovlstart++;
			 }
			 ng_end=ng_ovlstart-1; //MUST update ng_end here, even if no overlaps were found
		 } //guides present on the current chromosome
		bundle->refseq=lastref;
		bundle->start=currentstart;
		bundle->end=currentend;
	 } //<---- new bundle just started

	 if (currentend<(int)brec->end) {
		 //current read extends the bundle
		 //this might not happen if a longer guide had already been added to the bundle
		 currentend=brec->end;
		 if (guides) { //add any newly overlapping guides to bundle
			 bool cend_changed;
			 do {
				 cend_changed=false;
				 while (ng_end+1<ng && (int)(*guides)[ng_end+1]->start<=currentend) {
					 ++ng_end;
					 //more transcripts overlapping this bundle?
					 if ((int)(*guides)[ng_end]->end>=currentstart) {
						 //it should really overlap the bundle
						 bundle->keepGuide((*guides)[ng_end],
								  &guides_RC_tdata, &guides_RC_exons, &guides_RC_introns);
						 if(currentend<(int)(*guides)[ng_end]->end) {
							 currentend=(*guides)[ng_end]->end;
							 cend_changed=true;
						 }
				 	 }
				 }
			 } while (cend_changed);
		 }
	 } //adjusted currentend and checked for overlapping reference transcripts
	 GReadAlnData alndata(brec, 0, nh, hi);
     bundle->evalReadAln(alndata, xstrand); //xstrand, nh);
     if (xstrand=='+') alndata.strand=1;
		else if (xstrand=='-') alndata.strand=-1;
     //GMessage("%s\t%c\t%d\thi=%d\n",brec->name(), xstrand, alndata.strand,hi);
	 //countFragment(*bundle, *brec, hi,nh); // we count this in build_graphs to only include mapped fragments that we consider correctly mapped
	 //fprintf(stderr,"fragno=%d fraglen=%lu\n",bundle->num_fragments,bundle->frag_len);if(bundle->num_fragments==100) exit(0);
	 //if (!ballgown || ref_overlap)
	 processRead(currentstart, currentend, *bundle, hashread, alndata);
			  // *brec, strand, nh, hi);
 } //for each read alignment
 if (guided && no_ref_used)
    GMessage("WARNING: no reference transcripts were found for the genomic sequences where reads were mapped!\n"
    		"Please make sure the -G annotation file uses the same naming convention for the genome sequences.\n");
 //cleaning up
 delete brec;
 bamreader.bclose();
#ifndef NOTHREADS
 for (int t=0;t<num_cpus;t++)
	 threads[t].join();
 if (verbose) {
   printTime(stderr);
   GMessage(" All threads finished.\n");
 }
 delete[] threads;
 delete[] bundles;
#else
 if (verbose) {
    printTime(stderr);
    GMessage(" Done.\n");
 }
#endif

#ifdef B_DEBUG
 fclose(dbg_out);
#endif
 //if (f_out && f_out!=stdout) fclose(f_out);
 fclose(f_out);
 if (c_out && c_out!=stdout) fclose(c_out);


 if(verbose) {
	 GMessage("Total count of aligned fragments: %g\n",Num_Fragments);
	 GMessage("Fragment coverage length:%g\n",Frag_Len/Num_Fragments);
 }

 f_out=stdout;
 if(outfname!="stdout") {
	 f_out=fopen(outfname.chars(), "w");
	 if (f_out==NULL) GError("Error creating output file %s\n", outfname.chars());
 }

 fprintf(f_out,"# ");
 args.printCmdLine(f_out);
 fprintf(f_out,"# StringTie version %s\n",VERSION);

 FILE *g_out;
 if(geneabundance) {
	 g_out=fopen(genefname.chars(),"w");
	 if (g_out==NULL) GError("Error creating output file %s\n", genefname.chars());
	 fprintf(g_out,"Gene ID\tGene Name\tStrand\tStart\tEnd\tLength\tCoverage\tFPKM\tTPM\n");
 }

 FILE* t_out=fopen(tmpfname.chars(),"rt");
 if (t_out!=NULL) {
	 char* linebuf=NULL;
	 int linebuflen=5000;
     GMALLOC(linebuf, linebuflen);
	 int nl;
	 int istr;
	 int tlen;
	 float tcov;
	 //float fpkm;
	 float calc_fpkm;
	 float calc_tpm;
	 int t_id;
	 while(fgetline(linebuf,linebuflen,t_out)) {
		 //sscanf(linebuf,"%d %d %d %g %g", &nl, &tlen, &t_id, &fpkm, &tcov);
		 sscanf(linebuf,"%d %d %d %d %g", &istr, &nl, &tlen, &t_id, &tcov);
		 calc_fpkm=tcov*1000000000/Frag_Len;
		 calc_tpm=tcov*1000000/Cov_Sum;
		 /*
		 double eff_len=1;
		 if((float)tlen>AvgFrag) eff_len=tlen-AvgFrag+1;
		 calc_fpkm2=calc_fpkm*tlen/eff_len;
		 */

		 //fprintf(stderr,"tid=%d tlen=%d fpkm=%g calc_fpkm=%g tcov=%g new_fpkm=%g newcalc_fpkm=%g calc_fpkm3=%g\n",t_id,tlen,fpkm,calc_fpkm,tcov,tcov*tlen*1000000000/(efflen*SumReads),calc_fpkm2,calc_fpkm3);
		 if(istr) { // this is a transcript
			 if (ballgown && t_id>0) {
				 guides_RC_tdata[t_id-1]->fpkm=calc_fpkm;
				 guides_RC_tdata[t_id-1]->cov=tcov;
			 }
			 for(int i=0;i<nl;i++) {
				 fgetline(linebuf,linebuflen,t_out);
				 if(!i) {
					 //linebuf[strlen(line)-1]='\0';
					 fprintf(f_out,"%s",linebuf);
					 fprintf(f_out," FPKM \"%.6f\";",calc_fpkm);
					 fprintf(f_out," TPM \"%.6f\";",calc_tpm);
					 //fprintf(f_out,"FPKM \"%.6f\"; calculated_FPKM \"%.6f\";",tcov*1000000000/Frag_Len,fpkm*1000000000/(Num_Fragments*tlen));
					 //fprintf(f_out,"flen \"%.6f\"; FPKM \"%.6f\";",fpkm,fpkm*1000000000/Num_Fragments);
					 fprintf(f_out,"\n");
				 }
				 else fprintf(f_out,"%s\n",linebuf);
			 }
		 }
		 else { // this is a gene -> different file pointer
			 fgetline(linebuf,linebuflen,t_out);
			 fprintf(g_out,"%s\t%.6f\t%.6f\n",linebuf,calc_fpkm,calc_tpm);
		 }
	 }
	 fclose(f_out);
	 fclose(t_out);
	 if(geneabundance) fclose(g_out);
	 GFREE(linebuf);
	 remove(tmpfname.chars());
 }
 else {
	 fclose(f_out);
	 GError("No temporary file %s present!\n",tmpfname.chars());
 }

 //lastly, for ballgown, rewrite the tdata file with updated cov and fpkm
 if (ballgown) {
	 rc_writeRC(guides_RC_tdata, guides_RC_exons, guides_RC_introns,
			 f_tdata, f_edata, f_idata, f_e2t, f_i2t);
 }

 gffnames_unref(gseqNames); //deallocate names collection


#ifdef GMEMTRACE
 if(verbose) GMessage(" Max bundle memory: %6.1fMB for bundle %s\n", maxMemRS/1024, maxMemBundle.chars());
#endif
} // -- END main

//----------------------------------------
char* sprintTime() {
	static char sbuf[32];
	time_t ltime; /* calendar time */
	ltime=time(NULL);
	struct tm *t=localtime(&ltime);
	sprintf(sbuf, "%02d_%02d_%02d:%02d:%02d",t->tm_mon+1, t->tm_mday,
			t->tm_hour, t->tm_min, t->tm_sec);
	return(sbuf);
}


GStr Process_Options(GArgs* args) {

	if (args->getOpt('h') || args->getOpt("help")) {
		GMessage("%s",USAGE);
	    exit(1);
	}

	 debugMode=(args->getOpt("debug")!=NULL || args->getOpt('D')!=NULL);
	 fast=!(args->getOpt('F')!=NULL);
	 verbose=(args->getOpt('v')!=NULL);
	 if (verbose) {
	     fprintf(stderr, "Command line was:\n");
	     args->printCmdLine(stderr);
	 }
	 complete=!(args->getOpt('i')!=NULL);
	 trim=!(args->getOpt('t')!=NULL);
	 includesource=!(args->getOpt('z')!=NULL);
	 EM=(args->getOpt('y')!=NULL);
	 weight=(args->getOpt('w')!=NULL);
	 GStr s=args->getOpt('m');
	 if (!s.is_empty()) {
	   mintranscriptlen=s.asInt();
	   if (mintranscriptlen<30) 
	     GError("Error: invalid -m value, must be >=30)\n");
	 }
     s=args->getOpt('x');
     if (!s.is_empty()) {
    	 //split by comma and populate excludeGSeqs
    	 s.startTokenize(" ,\t");
    	 GStr chrname;
    	 while (s.nextToken(chrname)) {
    		 excludeGseqs.Add(chrname.chars(),new int(0));
    	 }
     }
	 s=args->getOpt('n');
	 if (!s.is_empty()) {
		 sensitivitylevel=s.asInt();
		 if(sensitivitylevel<0) {
			 sensitivitylevel=0;
			 GMessage("sensitivity level out of range: setting sensitivity level at 0\n");
		 }
		 if(sensitivitylevel>3) {
			 sensitivitylevel=3;
			 GMessage("sensitivity level out of range: setting sensitivity level at 2\n");
		 }
	 }

	 s=args->getOpt('g');
	 if (!s.is_empty()) bundledist=s.asInt();
	 s=args->getOpt('p');
	 if (!s.is_empty()) {
		   num_cpus=s.asInt();
		   if (num_cpus<=0) num_cpus=1;
	 }

	 s=args->getOpt('a');
	 if (!s.is_empty()) junctionsupport=(uint)s.asInt();
	 s=args->getOpt('j');
	 if (!s.is_empty()) junctionthr=s.asInt();
	 s=args->getOpt('c');
	 if (!s.is_empty()) {
		 readthr=(float)s.asDouble();
		 if (readthr<0.001) {
			 GError("Error: invalid -c value, must be >=0.001)\n");
		 }
	 }
	 s=args->getOpt('l');
	 if (!s.is_empty()) label=s;
	 s=args->getOpt('f');
	 if (!s.is_empty()) {
		 isofrac=(float)s.asDouble();
		 if(isofrac>=1) GError("Miminum isoform fraction (-f coefficient: %f) needs to be less than 1\n",isofrac);
	 }
	 s=args->getOpt('M');
	 if (!s.is_empty()) {
		 mcov=(float)s.asDouble();
	 }

	 genefname=args->getOpt('A');
	 if(!genefname.is_empty()) {
		 geneabundance=true;
	 }

	 //f_out=stdout;
	 tmpfname=args->getOpt('o');
	 outfname="stdout";
	 out_dir="./";
	 if (!tmpfname.is_empty() && tmpfname!="-") {
		 outfname=tmpfname;
		 int pidx=outfname.rindex('/');
		 if (pidx>=0) //path given
			 out_dir=outfname.substr(0,pidx+1);
	 }
	 else { // stdout
		tmpfname=outfname;
		char *stime=sprintTime();
		tmpfname+='.';
		tmpfname+=stime;
	 }
	 if (out_dir!="./") {
		 if (fileExists(out_dir.chars())==0) {
			//directory does not exist, create it
			Gmkdir(out_dir.chars());
		 }
	 }
#ifdef B_DEBUG
	 GStr dbgfname(tmpfname);
	 dbgfname+=".dbg";
	 dbg_out=fopen(dbgfname.chars(), "w");
	 if (dbg_out==NULL) GError("Error creating debug output file %s\n", dbgfname.chars());
#endif

	 tmpfname+=".tmp";
	 f_out=fopen(tmpfname.chars(), "w");
	 if (f_out==NULL) GError("Error creating output file %s\n", tmpfname.chars());
     /*
	 if (args->getOpt('O')) {
		 singlePass=false;
		 maxReadCov=0;
	 }
     */

	 if (args->getOpt('S')) {
		 // sensitivitylevel=2; no longer supported from version 1.0.3
		 sensitivitylevel=1;
	 }

	 // coverage saturation no longer used after version 1.0.4; left here for compatibility with previous versions
	 s=args->getOpt('s');
	 if (!s.is_empty()) {
		 GMessage("Coverage saturation parameter is deprecated starting at version 1.0.5");
		 /*
		 int r=s.asInt();
		 if (r<2) {
			 GMessage("Warning: invalid -s value, setting coverage saturation threshold, using default (%d)\n", maxReadCov);
		 }
		 else maxReadCov=r;
		 */
	 }

	 if (args->getOpt('G')) {
	   guidegff=args->getOpt('G');
	   if (fileExists(guidegff.chars())>1)
	        guided=true;
	   else GError("Error: reference annotation file (%s) not found.\n",
	             guidegff.chars());
	 }

	 eonly=(args->getOpt('e')!=NULL);
	 if(eonly && !guided)
		 GError("Error: invalid -e usage, GFF reference not given (-G option required).\n");

	 ballgown_dir=args->getOpt('b');
	 ballgown=(args->getOpt('B')!=NULL);
	 if (ballgown && !ballgown_dir.is_empty()) {
		 GError("Error: please use either -B or -b <path> options, not both.");
	 }
	 if (ballgown) ballgown_dir=out_dir;
	 else if (!ballgown_dir.is_empty()) {
		    ballgown=true;
		    ballgown_dir.chomp('/');ballgown_dir+='/';
			if (fileExists(ballgown_dir.chars())==0) {
				//directory does not exist, create it
				Gmkdir(ballgown_dir.chars());
			}
	 	  }
	 if (ballgown && !guided)
		 GError("Error: invalid -B/-b usage, GFF reference not given (-G option required).\n");

	 /* s=args->getOpt('P');
	 if (!s.is_empty()) {
		 if(!guided) GError("Error: option -G with reference annotation file has to be specified.\n");
		 c_out=fopen(s.chars(), "w");
		 if (c_out==NULL) GError("Error creating output file %s\n", s.chars());
		 partialcov=true;
	 }
	 else { */
		 s=args->getOpt('C');
		 if (!s.is_empty()) {
			 if(!guided) GError("Error: invalid -C usage, GFF reference not given (-G option required).\n");
			 c_out=fopen(s.chars(), "w");
			 if (c_out==NULL) GError("Error creating output file %s\n", s.chars());
		 }
	 //}
	 int numbam=args->startNonOpt();
#ifndef GFF_DEBUG
	 if (numbam==0 || numbam>1) {
	 	 GMessage("%s\nError: no BAM input file provided!\n",USAGE);
	 	 exit(1);
	 }
#endif

	 s=args->nextNonOpt();
	 return(s);
}

//---------------
bool moreBundles() { //getter (interogation)
	bool v=true;
#ifndef NOTHREADS
  GLockGuard<GFastMutex> lock(bamReadingMutex);
#endif
  v = ! NoMoreBundles;
  return v;
}

void noMoreBundles() {
#ifndef NOTHREADS
		bamReadingMutex.lock();
		NoMoreBundles=true;
		bamReadingMutex.unlock();
		queueMutex.lock();
		bundleWork &= ~(int)0x01; //clear bit 0;
		queueMutex.unlock();
		bool areThreadsWaiting=true;
		do {
		  waitMutex.lock();
		   areThreadsWaiting=(threadsWaiting>0);
		  waitMutex.unlock();
		  if (areThreadsWaiting) {
		    DBGPRINT("##> NOTIFY ALL workers: no more data!\n");
		    haveBundles.notify_all();
		    this_thread::sleep_for(chrono::milliseconds(10));
		    waitMutex.lock();
		     areThreadsWaiting=(threadsWaiting>0);
		    waitMutex.unlock();
		    this_thread::sleep_for(chrono::milliseconds(10));
		  }
		} while (areThreadsWaiting); //paranoid check that all threads stopped waiting
#else
	  NoMoreBundles=true;
#endif
}

void processBundle(BundleData* bundle) {
	if (verbose) {
	#ifndef NOTHREADS
			GLockGuard<GFastMutex> lock(logMutex);
	#endif
		printTime(stderr);
		/*GMessage(">bundle %s:%d-%d(%d) (%djs, %d guides) loaded, begins processing...\n",
				bundle->refseq.chars(), bundle->start, bundle->end, bundle->numreads,
                bundle->junction.Count(), bundle->keepguides.Count());*/
		GMessage(">bundle %s:%d-%d(%lu) (%d guides) loaded, begins processing...\n",
				bundle->refseq.chars(), bundle->start, bundle->end, bundle->numreads,
                bundle->keepguides.Count());
#ifdef GMEMTRACE
		double vm,rsm;
		get_mem_usage(vm, rsm);
		GMessage("\t\tstart memory usage: %6.1fMB\n",rsm/1024);
		if (rsm>maxMemRS) {
			maxMemRS=rsm;
			maxMemVM=vm;
			maxMemBundle.format("%s:%d-%d(%d)", bundle->refseq.chars(), bundle->start, bundle->end, bundle->readlist.Count());
		}
#endif
	}
#ifdef B_DEBUG
	for (int i=0;i<bundle->keepguides.Count();++i) {
		GffObj& t=*(bundle->keepguides[i]);
		RC_TData* tdata=(RC_TData*)(t.uptr);
		fprintf(dbg_out, ">%s (t_id=%d) %s%c %d %d\n", t.getID(), tdata->t_id, t.getGSeqName(), t.strand, t.start, t.end );
		for (int fe=0;fe < tdata->t_exons.Count(); ++fe) {
			RC_Feature& exoninfo = *(tdata->t_exons[fe]);
			fprintf(dbg_out, "%d\texon\t%d\t%d\t%c\t%d\t%d\n", exoninfo.id, exoninfo.l, exoninfo.r,
					    exoninfo.strand, exoninfo.rcount, exoninfo.ucount);
			if (! (exoninfo==*(bundle->rc_data->guides_RC_exons->Get(exoninfo.id-1))))
				 GError("exoninfo with id (%d) not matching!\n", exoninfo.id);
		}
		for (int fi=0;fi < tdata->t_introns.Count(); ++fi) {
			RC_Feature& introninfo = *(tdata->t_introns[fi]);
			fprintf(dbg_out, "%d\tintron\t%d\t%d\t%c\t%d\t%d\n", introninfo.id, introninfo.l, introninfo.r,
					introninfo.strand, introninfo.rcount, introninfo.ucount);
			if (! (introninfo==*(bundle->rc_data->guides_RC_introns->Get(introninfo.id-1))))
				 GError("introninfo with id (%d) not matching!\n", introninfo.id);
		}
		//check that IDs are properly assigned
		if (tdata->t_id!=(uint)t.udata) GError("tdata->t_id(%d) not matching t.udata(%d)!\n",tdata->t_id, t.udata);
		if (tdata->t_id!=bundle->rc_data->guides_RC_tdata->Get(tdata->t_id-1)->t_id)
			 GError("tdata->t_id(%d) not matching rc_data[t_id-1]->t_id (%d)\n", tdata->t_id, bundle->rc_data->g_tdata[tdata->t_id-1]->t_id);

	}
#endif
	//int ngenes=infer_transcripts(bundle, fast | bundle->covSaturated);
	int ngenes=infer_transcripts(bundle, fast);
	if (ballgown && bundle->rc_data) {
		rc_update_exons(*(bundle->rc_data));
	}
	if (bundle->pred.Count()>0) {
#ifndef NOTHREADS
		GLockGuard<GFastMutex> lock(printMutex);
#endif
		GeneNo=printResults(bundle, ngenes, GeneNo, bundle->refseq);
	}

	if (bundle->num_fragments) {
		#ifndef NOTHREADS
				GLockGuard<GFastMutex> lock(countMutex);
		#endif
		Num_Fragments+=bundle->num_fragments;
		Frag_Len+=bundle->frag_len;
		Cov_Sum+=bundle->sum_cov;
	}

	if (verbose) {
		#ifndef NOTHREADS
				GLockGuard<GFastMutex> lock(logMutex);
		#endif
	  printTime(stderr);
	  /*
	  SumReads+=bundle->sumreads;
	  SumFrag+=bundle->sumfrag;
	  NumCov+=bundle->num_cov;
	  NumReads+=bundle->num_reads;
	  NumFrag+=bundle->num_frag;
	  NumFrag3+=bundle->num_fragments3;
	  SumFrag3+=bundle->sum_fragments3;
	  fprintf(stderr,"Number of fragments in bundle: %g with length %g\n",bundle->num_fragments,bundle->frag_len);
	  */
	  fprintf(stderr,"Number of fragments in bundle: %g with sum %g\n",bundle->num_fragments,bundle->frag_len);
	  GMessage("^bundle %s:%d-%d(%d) done (%d processed potential transcripts).\n",bundle->refseq.chars(),
	  		bundle->start, bundle->end, bundle->readlist.Count(), bundle->pred.Count());
	#ifdef GMEMTRACE
		    double vm,rsm;
		    get_mem_usage(vm, rsm);
		    GMessage("\t\tfinal memory usage: %6.1fMB\n",rsm/1024);
		    if (rsm>maxMemRS) {
			    maxMemRS=rsm;
			    maxMemVM=vm;
			    maxMemBundle.format("%s:%d-%d(%d)", bundle->refseq.chars(), bundle->start, bundle->end, bundle->readlist.Count());
		    }
	#endif
	    }
	bundle->Clear();
}

#ifndef NOTHREADS


bool noThreadsWaiting() {
	waitMutex.lock();
	int v=threadsWaiting;
	waitMutex.unlock();
	return (v<1);
}

bool waitForThreads() {
	bool noneWaiting=true;
	DBGPRINT("##> waiting for a thread to become available..\n");
	while (noneWaiting) {
	  waitMutex.lock();
	  noneWaiting=(threadsWaiting<1);
	  waitMutex.unlock();
	  if (noneWaiting) {
		DBGPRINT("##>all threads busy, sleep_for 2ms\n");
	    this_thread::sleep_for(chrono::milliseconds(2));
	  }
	}
 DBGPRINT("##> there are workers ready now.\n");
 return(!noneWaiting);
}

void workerThread(GThreadData& td) {
	GPVec<BundleData>* bundleQueue = (GPVec<BundleData>*)td.udata;
	//wait for a ready bundle in the queue, until there is no hope for incoming bundles
	DBGPRINT2("---->> Thread%d starting..\n",td.thread->get_id());
	DBGPRINT2("---->> Thread%d locking queueMutex..\n",td.thread->get_id());
	queueMutex.lock(); //enter wait-for-notification loop
	while (bundleWork) {
		DBGPRINT3("---->> Thread%d: waiting.. (queue len=%d)\n",td.thread->get_id(), bundleQueue->Count());
		waitMutex.lock();
		 threadsWaiting++;
		queueMutex.unlock();
		waitMutex.unlock();
		haveThreads.notify_one(); //in case main thread is waiting
		this_thread::yield();
		queueMutex.lock();
		while (bundleWork && bundleQueue->Count()==0) {
		    haveBundles.wait(queueMutex);//unlocks queueMutex and wait until notified
		               //when notified, locks queueMutex and resume
		}
		waitMutex.lock();
		if (threadsWaiting>0) threadsWaiting--;
		waitMutex.unlock();
		DBGPRINT3("---->> Thread%d: awakened! (queue len=%d)\n",td.thread->get_id(),bundleQueue->Count());
		BundleData* readyBundle=NULL;
		if ((bundleWork & 0x02)!=0 && (readyBundle=bundleQueue->Pop())!=NULL) { //is bit 1 set?
				if (bundleQueue->Count()==0)
					 bundleWork &= ~(int)0x02; //clear bit 1 (queue is empty)
				//Num_Fragments+=readyBundle->num_fragments;
				//Frag_Len+=readyBundle->frag_len;
				queueMutex.unlock();
				processBundle(readyBundle);
				DBGPRINT2("---->> Thread%d processed bundle, now locking back queueMutex\n", td.thread->get_id());
				dataMutex.lock();
				dataClear.Push(readyBundle->idx);
				dataMutex.unlock();
				haveClear.notify_one(); //inform main thread
				this_thread::yield();
				queueMutex.lock();
				DBGPRINT2("---->> Thread%d locked back queueMutex\n", td.thread->get_id());

		}
		//haveThreads.notify_one();
	} //while there is reason to live
	queueMutex.unlock();
	DBGPRINT2("---->> Thread%d DONE.\n", td.thread->get_id());
}

bool queuePopped(GPVec<BundleData>& bundleQueue, int prevCount) {
  int c;
  queueMutex.lock();
   c=bundleQueue.Count();
  queueMutex.unlock();
  DBGPRINT3("##> post-notification check: qlen is now %d (was %d)\n", c, prevCount);
  return (c==0 || c<prevCount);
}

//prepare the next available bundle slot for loading
int waitForData(BundleData* bundles) {
	int bidx=-1;
	/*
	while (bidx<0) {
	  dataMutex.lock();
	  if (dataClear.Count()>0) {
	    bidx=dataClear.Pop();
	    bundles[bidx].status=BUNDLE_STATUS_LOADING;
	    dataMutex.unlock();
	    return bidx;
	    }
	  dataMutex.unlock();
	  //this_thread::sleep_for(chrono::milliseconds(20));
	}
	return -1; // should NEVER happen
	*/
	dataMutex.lock();
	while (dataClear.Count()==0) {
		haveClear.wait(dataMutex);
	}
	bidx=dataClear.Pop();
	if (bidx>=0) {
	  bundles[bidx].status=BUNDLE_STATUS_LOADING;
	}
	dataMutex.unlock();
	return bidx;
}

#endif
