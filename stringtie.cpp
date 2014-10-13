#include "rlink.h"
#ifndef NOTHREADS
#include "GThreads.h"
#endif

//comment this out if you don't want memory tracing in log
#ifdef GDEBUG
#define GMEMTRACE 1
#else
#undef GMEMTRACE
#endif

//#undef GMEMTRACE //-- comment out to track memory use for GDEBUG

#ifdef GMEMTRACE
#include "proc_mem.h"
#endif

#define VERSION "0.97"
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
 stringtie <input.bam> [-G <guide_gff>] [-l <label>] [-o <out_gff>] [-p <cpus>]\n\
  [-v] [-a <min_anchor_len>] [-m <min_tlen>] [-j <min_anchor_cov>] [-n sens]\n\
  [-C <coverage_file_name>] [-s <maxcov>] [-c <min_bundle_cov>] [-g <bdist>]\n\
\nAssemble RNA-Seq alignments into potential transcripts.\n\
 \n\
 Options:\n\
 -G  reference annotation to use for guiding the assembly process (GTF/GFF3)\n\
 -l  name prefix for output transcripts (default: STRG)\n\
 -f  minimum isoform fraction (default: 0.15)\n\
 -m  minimum assembled transcript length to report (default 200bp)\n\
 -o  output file with the assembled transcripts (default: stdout)\n\
 -a  minimum anchor length for junctions (default: 10)\n\
 -j  minimum junction coverage (default: 1)\n\
 -t  disable trimming of predicted transcripts based on coverage (default: trimming enabled)\n\
 -c  minimum bundle reads per bp coverage to consider for assembly (default: 2.5)\n\
 -s  coverage saturation threshold; further read alignments will be\n\
     ignored in a region where a local coverage depth of <maxcov> \n\
     is reached (default: 1,000,000);\n\
 -v  verbose (log bundle processing details)\n\
 -e  abundance estimation only of input transcripts (for -G option)\n\
 -g  gap between read mappings triggering a new bundle (default: 50)\n\
 -S  more sensitive run (default: no)\n\
 -C  output file with all transcripts in reference that are fully\n\
     covered by reads\n\
 -M  fraction of bundle allowed to be covered by multi-hit reads (default:0.95)\n\
 -p  number of threads (CPUs) to use (default: 1)\n\
 "
/* 
 -n  sensitivity level: 0,1, or 2, 3, with 3 the most sensitive level (default 0)\n\
 -O  disable the coverage saturation limit and use a slower two-pass approach\n\
     to process the input alignments, collapsing redundant reads\n\
 -x  disable fast computing for transcript path; default: yes\n\
 -i  the reference annotation contains partial transcripts\n\
 -w  weight the maximum flow algorithm towards the transcript with higher rate (abundance); default: no\n\
 -y  include EM algorithm in max flow estimation; default: no\n\
 -z don't include source in the max flow algorithm\n\
 -P output file with all transcripts in reference that are partially covered by reads
 -M  fraction of bundle allowed to be covered by multi-hit reads (paper uses default: 1)\n\
 -c  minimum bundle reads per bp coverage to consider for assembly (paper uses default: 3)\n\
*/
//---- globals

FILE* f_out=NULL;
FILE* c_out=NULL;
GStr outfname;
GStr tmpfname;
bool guided=false;
bool trim=true;
bool fast=true;
bool eonly=false; // parameter -e
bool specific=false;
bool complete=true;
bool partialcov=false;
int num_cpus=1;
int mintranscriptlen=200; // minimum length for a transcript to be printed
int sensitivitylevel=1;
int junctionsupport=10; // anchor length for junction to be considered well supported <- consider shorter??
int junctionthr=1; // number of reads needed to support a particular junction
float readthr=2.5;     // read coverage per bundle bp to accept it; otherwise considered noise; paper uses 3
uint bundledist=50;  // reads at what distance should be considered part of separate bundles
float mcov=0.95; // fraction of bundle allowed to be covered by multi-hit readsl paper uses 1

// different options of implementation reflected with the next three options
bool includesource=true;
bool EM=false;
bool weight=false;

float isofrac=0.15;
GStr label("STRG");
GStr guidegff;

bool debugMode=false;
bool verbose=false;

int maxReadCov=1000000; //max local read coverage (changed with -s option)
//no more reads will be considered for a bundle if the local coverage exceeds this value
//(each exon is checked for this)

bool singlePass=true; //-O will set this to False


int GeneNo=0; //-- global "gene" counter
unsigned long long int Num_Fragments=0; //global fragment counter (aligned pairs)
unsigned long long int Frag_Len=0;
//bool firstPrint=true; //just for writing the GFF header before the first transcript is printed

GffNames* gseqNames=NULL; //used as a dictionary for genomic sequence names and ids

#ifdef GMEMTRACE
 double maxMemRS=0;
 double maxMemVM=0;
 GStr maxMemBundle;
#endif


#ifndef NOTHREADS

GFastMutex dataMutex; //to manage availability of data records ready to be loaded by main thread
GVec<int> dataClear; //indexes of data bundles cleared for loading by main thread (clear data pool)

GFastMutex waitMutex; // for main program to make sure there are threads ready/waiting
int threadsWaiting=0; // how many worker threads are waiting

GFastMutex printMutex; //for writing the output to file
GFastMutex logMutex; //only when verbose - to avoid mangling the log output
GMutex queueMutex; //whenever bundleQueue is updated
GFastMutex bamReadingMutex;
GConditionVar haveBundles; //will notify all threads when bundles are pushed in the ready queue
                           //or no more bundles are coming

int bundleWork=1; // bit 0 set if bundles are still being prepared (BAM file not exhausted yet)
                  // bit 1 set if there are Bundles ready in the queue
#endif

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
   "debug;help;xyzwShvtin:j:s:D:G:C:l:m:o:a:j:c:f:p:g:P:M:");
 args.printError(USAGE, true);

 GStr bamfname=Process_Options(&args);
 // == Done argument processing.

 GVec<GRefData> refguides; // plain vector with transcripts for each chromosome

#ifdef DEBUGPRINT
  verbose=true;
#endif

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
   refguides.setCount(gffr.gseqStats.Count()); //maximum gseqid
   for (int i=0;i<gffr.gflst.Count();i++) {
     GffObj* m=gffr.gflst[i];
     GRefData& grefdata = refguides[m->gseq_id];
     grefdata.add(&gffr, m); //transcripts already sorted by location
   }
	 if (verbose) {
		 printTime(stderr);
		 GMessage(" %d reference transcripts loaded.\n", gffr.gflst.Count());
	 }

 }


 // ---OK, here we do the input processing:

 gseqNames=GffObj::names; //might have been populated already by gff data
 gffnames_ref(gseqNames);  //initialize the names collection if not guided

 GBamReader bamreader(bamfname.chars());

 GHash<int> hashread;      //read_name:pos:hit_index => readlist index
 //my %hashjunction;  //junction coords and strand => junction index
 // we won't need this because we can quick-search in junction GList directly
 //my @guides=(); //set of annotation transcript for the current locus
 GList<GffObj>* guides=NULL; //list of transcripts on a specific chromosome

 int currentstart=0, currentend=0;
 int ng_start=0;
 int ng_end=-1;
 int ng=0;
 GStr lastref;
 int lastref_id=-1; //last seen gseq_id
 // int ncluster=0; used it for debug purposes only

#ifndef NOTHREADS
 GThread* threads=new GThread[num_cpus];
 GPVec<BundleData> bundleQueue(false);
 BundleData* bundles=new BundleData[num_cpus+1]; //extra one being prepared while all others are processed
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
 while (more_alns) {
	 bool chr_changed=false;
	 int pos=0;
	 const char* rname=NULL;
	 char strand=0;
	 char xstrand=0;
	 int nh=1;
	 int hi=0;
	 int gseq_id=lastref_id;  //current chr id
	 bool new_bundle=false;
	 delete brec;
	 if ((brec=bamreader.next())!=NULL) {
		 if (brec->isUnmapped()) continue;
		 rname=brec->refName();
		 if (rname==NULL) GError("Error: cannot retrieve target seq name from BAM record!\n");
		 pos=brec->start; //BAM is 0 based, but GBamRecord makes it 1-based
		 chr_changed=(lastref.is_empty() || lastref!=rname);
		 if (chr_changed) {
			 gseq_id=gseqNames->gseqs.addName(rname);
		 }
		 xstrand=brec->spliceStrand();
		 if (xstrand=='+') strand=1;
		 else if (xstrand=='-') strand=-1;
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
			 if (guides && ng_end>=ng_start) {
				 for (int gi=ng_start;gi<=ng_end;gi++)
					 bundle->keepguides.Add((*guides)[gi]);
			 }
			// geneno=infer_transcripts(geneno, lastref, $label,\@readlist,$readthr,\@junction,$junctionthr,$mintranscriptlen,\@keepguides);
			// (readthr, junctionthr, mintranscriptlen are globals)
			bundle->getReady(currentstart, currentend);
#ifndef NOTHREADS
			//push this in the bundle queue, where it'll be picked up by the threads
			DBGPRINT2("##> Locking queueMutex to push loaded bundle into the queue (bundle.start=%d)\n", bundle->start);
			queueMutex.lock();
			bundleQueue.Push(bundle);
			bundleWork |= 0x02; //set bit 1
			int qCount=bundleQueue.Count();
			queueMutex.unlock();
			do {
			     waitForThreads();
			     DBGPRINT("##> NOTIFY any thread...\n");
			     haveBundles.notify_one();
			     //this_thread::sleep_for(chrono::milliseconds(1));
			     sleep(0);
			} while (!queuePopped(bundleQueue, qCount));
#else //no threads
			Num_Fragments+=bundle->num_fragments;
			Frag_Len+=bundle->frag_len;
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
		 }

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
			 lastref=rname;
			 lastref_id=gseq_id;
			 currentend=0;
		 }

		 if (!more_alns) {
				if (verbose) {
#ifndef NOTHREADS
						GLockGuard<GFastMutex> lock(logMutex);
#endif
					printTime(stderr);
					GMessage(" %llu aligned fragments found.\n", Num_Fragments);
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
		 if (guides) { //guided and guides!=NULL
			 ng_start=ng_end+1;
			 while (ng_start<ng && (int)(*guides)[ng_start]->end < pos) { ng_start++; } // skip guides that have no read coverage
			 if(ng_start<ng && (int)(*guides)[ng_start]->start<pos) {
				 currentstart=(*guides)[ng_start]->start;
				 currentend=(*guides)[ng_start]->end;
				 ng_end=ng_start;
				 while(ng_end+1<ng && (int)(*guides)[ng_end+1]->start<=pos) {
					 ng_end++;
					 if(currentend<(int)(*guides)[ng_end]->end) {
						 currentend=(*guides)[ng_end]->end;
					 }
				 }
			 }
		 } //guides present on the current chromosome
		bundle->refseq=lastref;
		bundle->start=currentstart;
		bundle->end=currentend;
	 } //<---- new bundle
	 //currentend=process_read(currentstart, currentend, bundle->readlist, hashread,
		//	 bundle->junction, *brec, strand, nh, hi, bundle->bpcov);
     currentend=processRead(currentstart, currentend, *bundle, hashread, *brec, strand, nh, hi);
   //update current end to be at least as big as the start of the read pair in the fragment?? -> maybe not because then I could introduce some false positives with paired reads mapped badly

	 if(guides) { // I need to adjust end according to guides
		 while( ng_end+1 < ng && (int)(*guides)[ng_end+1]->start<=currentend) {
			 ng_end++;
			 if(currentend < (int)(*guides)[ng_end]->end) {
				 currentend=(*guides)[ng_end]->end;
			 }
		 }
	 }
 } //for each read alignment

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

 gffnames_unref(gseqNames); //deallocate names collection
 //if (f_out && f_out!=stdout) fclose(f_out);
 fclose(f_out);
 if (c_out) fclose(c_out);

 // write the FPKMs

 if(verbose) {
	 GMessage("Total count of aligned fragments: %llu\n",Num_Fragments);
	 //GMessage("Fragment length:%llu\n",Frag_Len);
	 GMessage("Average fragment length:%g\n",(float)Frag_Len/Num_Fragments);
 }

 f_out=stdout;
 if(outfname!="stdout") {
	 f_out=fopen(outfname.chars(), "w");
	 if (f_out==NULL) GError("Error creating output file %s\n", outfname.chars());
 }
 c_out=fopen(tmpfname.chars(),"rt");
 if (c_out!=NULL) {
	 char line[5001];
	 int nl;
	 int tlen;
	 float tcov;
	 float fpkm;
	 while(fgets(line,5000,c_out)) {
		 sscanf(line,"%d %d %g %g",&nl,&tlen,&fpkm,&tcov);
		 for(int i=0;i<nl;i++) {
			 fgets(line,5000,c_out);
			 if(!i) {
				 line[strlen(line)-1]='\0';
				 fprintf(f_out,"%s",line);
				 fprintf(f_out,"FPKM \"%.6f\";",tcov*1000000000/Frag_Len);
				 //fprintf(f_out,"FPKM \"%.6f\"; calculated_FPKM \"%.6f\";",tcov*1000000000/Frag_Len,fpkm*1000000000/(Num_Fragments*tlen));
				 //fprintf(f_out,"flen \"%.6f\"; FPKM \"%.6f\";",fpkm,fpkm*1000000000/Num_Fragments);
				 fprintf(f_out,"\n");
			 }
			 else fprintf(f_out,"%s",line);
		 }
	 }
	 fclose(f_out);
	 fclose(c_out);
	 remove(tmpfname.chars());
 }
 else {
	 fclose(f_out);
	 GError("No temporary file %s present!\n",tmpfname.chars());
 }

#ifdef GMEMTRACE
 if(verbose) GMessage(" Max bundle memory: %6.1fMB for bundle %s\n", maxMemRS/1024, maxMemBundle.chars());
#endif
} // -- END main

//----------------------------------------

int predCmp(const pointer p1, const pointer p2) {
	CPrediction *a=(CPrediction*)p1;
	CPrediction *b=(CPrediction*)p2;
	if(a->start < b->start) return 1;
	if(a->start > b->start) return -1;
	if(sensitivitylevel!=1 || sensitivitylevel!=2) {
		if(a->exons.Count() < b->exons.Count()) return 1;
		if(a->exons.Count() > b->exons.Count()) return -1;
		int i=0;
		uint a1=0;
		uint b1=0;
		while(i<a->exons.Count()) {
			if(a->exons[i].start<b->exons[i].start) {
				a1=a->exons[i].start;
				b1=b->exons[i].start;
				break;
			}
			if(a->exons[i].end<b->exons[i].end) {
				a1=a->exons[i].end;
				b1=b->exons[i].end;
				break;
			}
			i++;
		}
		if(a1) {
			if(a1<b1) return 1;
			if(a1>b1) return -1;
		}
	}
	return 0;
}

int predexCmp(const pointer p1, const pointer p2) {
	CPrediction *a=(CPrediction*)p1;
	CPrediction *b=(CPrediction*)p2;
	if(a->exons.Count() < b->exons.Count()) return -1;
	if(a->exons.Count() > b->exons.Count()) return 1;
	if(a->tlen < b->tlen) return -1;
	if(a->tlen > b->tlen) return 1;
	return 0;
}

int predcovCmp(const pointer p1, const pointer p2) {
	CPrediction *a=(CPrediction*)p1;
	CPrediction *b=(CPrediction*)p2;
	if(a->cov < b->cov) return 1;
	if(a->cov > b->cov) return -1;
	return 0;
}


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
	 fast=!(args->getOpt('x')!=NULL);
	 eonly=(args->getOpt('e')!=NULL);
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
	 if (!s.is_empty()) mintranscriptlen=s.asInt();

	 s=args->getOpt('n');
	 if (!s.is_empty()) {
		 sensitivitylevel=s.asInt();
		 if(sensitivitylevel<0) {
			 sensitivitylevel=0;
			 fprintf(stderr,"sensitivity level out of range: setting sensitivity level at 0\n");
		 }
		 if(sensitivitylevel>3) {
			 sensitivitylevel=3;
			 fprintf(stderr,"sensitivity level out of range: setting sensitivity level at 2\n");
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
	 if (!s.is_empty()) junctionsupport=s.asInt();
	 s=args->getOpt('j');
	 if (!s.is_empty()) junctionthr=s.asInt();
	 s=args->getOpt('c');
	 if (!s.is_empty()) readthr=(float)s.asDouble();
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

	 //f_out=stdout;
	 tmpfname=args->getOpt('o');
	 outfname="stdout";
	 if (!tmpfname.is_empty() && tmpfname!="-") {
		 outfname=tmpfname;
	 }
	 else { // s is stdout
		tmpfname=outfname;
		char *stime=sprintTime();
		tmpfname+='.';
		tmpfname+=stime;
	 }
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
		 sensitivitylevel=2;
	 }

	 s=args->getOpt('s');
	 if (!s.is_empty()) {
		 int r=s.asInt();
		 if (r<2) {
			 GMessage("Warning: invalid -s value, setting coverage saturation threshold, using default (%d)\n", maxReadCov);
		 }
		 else maxReadCov=r;
	 }
	 /*
	 {//DEBUG ONLY:
		 GStr fname(outfname);
		 fname+=".reads";
		 unlink(fname.chars());
	 }
	 */

	 if (args->getOpt('G')) {
	   guidegff=args->getOpt('G');
	   if (fileExists(guidegff.chars())>1)
	        guided=true;
	   else GError("Error: reference annotation file (%s) not found.\n",
	             guidegff.chars());
	 }


	 s=args->getOpt('P');
	 if (!s.is_empty()) {
		 if(!guided) GError("Error: option -G with reference annotation file has to be specified.\n");
		 c_out=fopen(s.chars(), "w");
		 if (c_out==NULL) GError("Error creating output file %s\n", s.chars());
		 partialcov=true;
	 }
	 else {
		 s=args->getOpt('C');
		 if (!s.is_empty()) {
			 c_out=fopen(s.chars(), "w");
			 if (c_out==NULL) GError("Error creating output file %s\n", s.chars());
		 }
	 }



	 int numbam=args->startNonOpt();
	 if (numbam==0 || numbam>1) {
	 	 GMessage("%s\nError: no BAM input file provided!\n",USAGE);
	 	 exit(1);
	 }

	 s=args->nextNonOpt();

	 return(s);
}

bool equal_pred(GList<CPrediction>& pred,int n1,int n2){

	if(pred[n1]->strand!=pred[n2]->strand) return false;

	//if(pred[n1]->start!=pred[n2]->start) return(false); // this allows genes with different start/ends to be merged together
	//if(pred[n1]->end!=pred[n2]->end) return(false);     // but I need to check if they overlap in case of single exons

	if((pred[n1]->end < pred[n2]->start) || (pred[n2]->end<pred[n1]->start)) return(false); // genes don't overlap

	int nex=pred[n1]->exons.Count();
	if(nex!=pred[n2]->exons.Count()) return(false);
	for(int i=0;i<nex;i++) {
		//if(pred[n1]->exons[i].start!=pred[n2]->exons[i].start) { fprintf(stderr,"ret false start[%d]: %d vs %d\n",i,pred[n1]->exons[i].start,pred[n2]->exons[i].start); return(false);}
		//if(pred[n1]->exons[i].end!=pred[n2]->exons[i].end) { fprintf(stderr,"ret false end[%d]: %d vs %d\n",i,pred[n1]->exons[i].end,pred[n2]->exons[i].end); return(false);}
		if(i>0 && (pred[n1]->exons[i].start!=pred[n2]->exons[i].start)) return(false);
		if(i<nex-1 && (pred[n1]->exons[i].end!=pred[n2]->exons[i].end)) return(false);
	}

	return(true);
}


CInterval *add_pred_to_cov(CInterval *maxcov, CPrediction* pred, bool *abundant=NULL) { // maybe I can eliminate some genes here


	//fprintf(stderr,"add pred: %d-%d with cov=%f to maxpos\n",pred->start,pred->end,pred->cov);

	if(maxcov==NULL || pred->end<maxcov->pos) { // prediction before current intervals or no current intervals
		uint linkstart=pred->end+1;
		CInterval *link=maxcov;
		if(maxcov && linkstart==maxcov->pos) { // next interval starts immediately after this one ends
			maxcov = new CInterval(pred->start,pred->cov,link);
		}
		else {
			CInterval *interval=new CInterval(linkstart,0,maxcov);
			maxcov = new CInterval(pred->start,pred->cov,interval);
		}
	}
	else { // I need to place current prediction
		CInterval *lastinterv=NULL;
		CInterval *interv=maxcov;
		while(interv && pred->start>=interv->pos) {
			lastinterv=interv;
			interv=interv->next;
		}
		if(interv){ // pred->start < interv->pos
			if(lastinterv) { // pred->start >= lastinterv->pos

				if(pred->end<interv->pos) {
					if(lastinterv->val<pred->cov) { // only in this case I am interested to do something, otherwise I am done
						if(pred->end+1<interv->pos) {
							lastinterv->next=new CInterval(pred->end+1,lastinterv->val,interv);
						}
						if(lastinterv->pos==pred->start) lastinterv->val=pred->cov;
						else {
							interv=lastinterv->next;
							lastinterv->next=new CInterval(pred->start,pred->cov,interv);
						}
					}
					else if(abundant && pred->cov<isofrac*lastinterv->val) *abundant=false;
				}
				else {
					float lastcov=lastinterv->val; // lastcov might be 0 here
					if(pred->cov>lastinterv->val) { // create new interval unless the position is in the same place as lastinterv
						if(pred->start==lastinterv->pos) lastinterv->val=pred->cov;
						else {
							lastinterv->next=new CInterval(pred->start,pred->cov,interv);
							lastinterv=lastinterv->next;
							lastcov=lastinterv->val;
						}
					}
					else if(abundant && pred->cov<isofrac*lastinterv->val) *abundant=false;
					while(interv && pred->end>=interv->pos) {
						if(interv->val<pred->cov) {
							if(lastinterv->val==pred->cov) { // I have to skip current interval
								lastinterv->next=interv->next;
								free(interv);
								interv=lastinterv;
							}
							else { lastcov=interv->val; interv->val=pred->cov;lastinterv=interv;}
						}
						else {
							if(abundant && pred->cov<isofrac*interv->val) *abundant=false;
							lastinterv=interv;
							lastcov=lastinterv->val;
						}
						interv=interv->next;
					}
					if(interv) { // pred->end<interv->pos: I might need to create new interval
						if(lastinterv->val<pred->cov) lastinterv->val=pred->cov;
						else if(abundant && pred->cov<isofrac*lastinterv->val) *abundant=false;
						uint newstart=pred->end+1;
						if(newstart<interv->pos) {
							lastinterv->next=new CInterval(newstart,lastcov,interv);
						}
					}
					else { // pred->end >= lastinterv->pos
						if(lastinterv->val<pred->cov) lastinterv->val=pred->cov;
						else if(abundant && pred->cov<isofrac*lastinterv->val) *abundant=false;
						lastinterv->next=new CInterval(pred->end+1);
					}
				}
			}
			else { // lastinterv == NULL
				maxcov=new CInterval(pred->start,pred->cov,interv);
				if(pred->end<interv->pos) {
					if(pred->end+1<interv->pos) maxcov->next=new CInterval(pred->end+1,0,interv);
				}
				else { // pred->end >= interv->pos
					lastinterv=maxcov;
					while(interv && pred->end>=interv->pos) {
						if(interv->val<pred->cov) {
							if(lastinterv->val==pred->cov) { // I have to skip current interval
								lastinterv->next=interv->next;
								free(interv);
							}
							else { interv->val=pred->cov;lastinterv=interv;}
						}
						else {
							if(abundant && pred->cov<isofrac*interv->val) *abundant=false;
							lastinterv=interv;
						}
						interv=interv->next;
					}
					if(interv) { // pred->end<interv->pos I might need to create new interval
						float lastcov=lastinterv->val;
						if(lastinterv->val<pred->cov) lastinterv->val=pred->cov;
						else if(abundant && pred->cov<isofrac*lastinterv->val) *abundant=false;
						uint newstart=pred->end+1;
						if(newstart<interv->pos) lastinterv->next=new CInterval(newstart,lastcov,interv);
					}
					else {
						if(lastinterv->val<pred->cov) lastinterv->val=pred->cov;
						else if(abundant && pred->cov<isofrac*lastinterv->val) *abundant=false;
						lastinterv->next=new CInterval(pred->end+1);
					}
				}
			}
		}
		else { // prediction is after the end of maxcov, or is equal to lastinterv->pos
			interv=new CInterval(pred->end+1);
			if(pred->start==lastinterv->pos) {
				if(pred->cov>=lastinterv->val) {
					lastinterv->val=pred->cov;
					lastinterv->next=interv;
				}
				else {
					if(abundant && pred->cov<isofrac*lastinterv->val) *abundant=false;
					lastinterv->next= new CInterval(pred->start+1,pred->cov,interv);
				}
			}
			else lastinterv->next = new CInterval(pred->start,pred->cov,interv);
		}
	}

	return(maxcov);
}

bool is_pred_above_frac(CInterval *maxcov,CPrediction* pred) {

	CInterval *lastinterv=NULL;
	while(maxcov && pred->start>=maxcov->pos) {
		lastinterv=maxcov;
		maxcov=maxcov->next;
	}
	if(lastinterv && ((pred->exons.Count()==1 && pred->cov<lastinterv->val) || pred->cov<isofrac*lastinterv->val)) return(false); // I need to deal with single exons too here
	while(maxcov && pred->end>=maxcov->pos) {
		if((pred->exons.Count()==1 && pred->cov<maxcov->val) || pred->cov<isofrac*maxcov->val) return(false);
		maxcov=maxcov->next;
	}
	return(true);
}

void delete_interval(CInterval *interv){
	if(interv) {
		if(interv->next) delete_interval(interv->next);
		delete interv;
	}
}


/*
// I don't use this
int print_cluster(GList<CPrediction>& pred,GVec<int>& genes,GVec<int>& transcripts,
		int nstart, int nend, int geneno,GStr& refname) {

  GVec<int> keep;

  CInterval *maxposcov=NULL; //remembers intervals of maximum positive coverage
  CInterval *maxnegcov=NULL; //remembers intervals of maximum positive coverage

  bool pos;
  bool neg;
  int lastadded=0;
  for(int n=nstart;n<=nend;n++) {
	  if(pred[n]->strand=='+') { pos=true;neg=false;}
	  else if(pred[n]->strand=='-') { pos=false;neg=true;}
	  else { pos=true;neg=true;}
	  if(n>nstart) {
		  if(equal_pred(pred,lastadded,n)) {
			  pred[lastadded]->cov+=pred[n]->cov;
			  if(pos) maxposcov=add_pred_to_cov(maxposcov,pred[lastadded]);
			  if(neg) maxnegcov=add_pred_to_cov(maxnegcov,pred[lastadded]);
			  continue;
		  }
	  }
	  if(pos) maxposcov=add_pred_to_cov(maxposcov,pred[n]);
	  if(neg) maxnegcov=add_pred_to_cov(maxnegcov,pred[n]);
	  keep.Add(n);
	  lastadded=n;
  }

  for(int i=0;i<keep.Count();i++) {
	  int n=keep[i];
	  if(pred[n]->strand=='+') { pos=true;neg=false;}
	  else if(pred[n]->strand=='-') { pos=false;neg=true;}
	  else { pos=false;neg=false;}
	  if((pos && is_pred_above_frac(maxposcov,pred[n])) || (neg && is_pred_above_frac(maxnegcov,pred[n])) ||
			  (!pos && !neg && is_pred_above_frac(maxnegcov,pred[n]))) {
		  if(genes[pred[n]->geneno]==-1) genes[pred[n]->geneno]=++geneno;
		  transcripts[pred[n]->geneno]++;
		  fprintf(f_out,"%s\tStringTie\ttranscript\t%d\t%d\t1000\t%c\t.\tgene_id \"%s.%d\"; transcript_id \"%s.%d.%d\"; cov \"%.6f\";\n",
				  refname.chars(),pred[n]->start,pred[n]->end,pred[n]->strand,label.chars(),genes[pred[n]->geneno],
				  label.chars(),genes[pred[n]->geneno],transcripts[pred[n]->geneno],pred[n]->cov);
		  for(int j=0;j<pred[n]->exons.Count();j++)
			  fprintf(f_out,"%s\tStringTie\texon\t%d\t%d\t1000\t%c\t.\tgene_id \"%s.%d\"; transcript_id \"%s.%d.%d\"; exon_number \"%d\"; cov \"%.6f\";\n",
			  		 refname.chars(),pred[n]->exons[j].start,pred[n]->exons[j].end,pred[n]->strand,label.chars(),genes[pred[n]->geneno],
			  		 label.chars(),genes[pred[n]->geneno],transcripts[pred[n]->geneno],j+1,pred[n]->cov); // maybe add exon coverage here
	  }
  }

  delete_interval(maxposcov);
  delete_interval(maxnegcov);

  return(geneno);
}
*/

/*
// I don't use this
int print_transcript_cluster(GList<CPrediction>& pred,GVec<int>& genes,GVec<int>& transcripts,
		int nstart, int nend, int geneno,GStr& refname) {

  GVec<int> keep;

  float maxcovpos=0;
  float maxcovneg=0;
  bool pos;
  bool neg;
  int lastadded=0;
  for(int n=nstart;n<=nend;n++) {
	  if(pred[n]->strand=='+') { pos=true;neg=false;}
	  else if(pred[n]->strand=='-') { pos=false;neg=true;}
	  else { pos=true;neg=true;}
	  if(n>nstart) {
		  if(equal_pred(pred,lastadded,n)) {
			  pred[lastadded]->cov+=pred[n]->cov;
			  if(pos && pred[lastadded]->cov > maxcovpos) {
				  maxcovpos=pred[lastadded]->cov;
			  }
			  if(neg && pred[lastadded]->cov > maxcovneg) {
				  maxcovneg=pred[lastadded]->cov;
			  }
			  continue;
		  }
	  }
	  if(pos && pred[n]->cov>maxcovpos) {
		  maxcovpos=pred[n]->cov;
	  }
	  if(neg && pred[n]->cov>maxcovneg) {
		  maxcovneg=pred[n]->cov;
	  }
	  keep.Add(n);
	  lastadded=n;
  }

  for(int i=0;i<keep.Count();i++) {
	  int n=keep[i];
	  if(pred[n]->strand=='+') { pos=true;neg=false;}
	  else if(pred[n]->strand=='-') { pos=false;neg=true;}
	  else { pos=false;neg=false;}
	  if((pos && ((pred[n]->exons.Count()==1 && pred[n]->cov>=maxcovpos) || (pred[n]->exons.Count()>1 && pred[n]->cov/maxcovpos>=isofrac))) ||
			  (neg && ((pred[n]->exons.Count()==1 && pred[n]->cov>=maxcovneg) || (pred[n]->exons.Count()>1 && pred[n]->cov/maxcovneg>=isofrac))) ||
			  (!pos && !neg && ((pred[n]->exons.Count()==1 && pred[n]->cov>=maxcovneg && pred[n]->cov>=maxcovpos) ||
					  (pred[n]->exons.Count()>1 && pred[n]->cov/maxcovpos>=isofrac && pred[n]->cov/maxcovneg>=isofrac)))) { // print this transcript
		  if(genes[pred[n]->geneno]==-1) genes[pred[n]->geneno]=++geneno;
		  transcripts[pred[n]->geneno]++;
		  fprintf(f_out,"%s\tStringTie\ttranscript\t%d\t%d\t1000\t%c\t.\tgene_id \"%s.%d\"; transcript_id \"%s.%d.%d\"; cov \"%.6f\";\n",
				  refname.chars(),pred[n]->start,pred[n]->end,pred[n]->strand,label.chars(),genes[pred[n]->geneno],
				  label.chars(),genes[pred[n]->geneno],transcripts[pred[n]->geneno],pred[n]->cov);
		  for(int j=0;j<pred[n]->exons.Count();j++)
			  fprintf(f_out,"%s\tStringTie\texon\t%d\t%d\t1000\t%c\t.\tgene_id \"%s.%d\"; transcript_id \"%s.%d.%d\"; exon_number \"%d\"; cov \"%.6f\";\n",
			  		 refname.chars(),pred[n]->exons[j].start,pred[n]->exons[j].end,pred[n]->strand,label.chars(),genes[pred[n]->geneno],
			  		 label.chars(),genes[pred[n]->geneno],transcripts[pred[n]->geneno],j+1,pred[n]->cov); // maybe add exon coverage here
	  }
  }

  return(geneno);
}
*/

int print_signcluster(char strand,GList<CPrediction>& pred,GVec<int>& genes,GVec<int>& transcripts,
		int nstart, int nend, int geneno,GStr& refname) {

  GVec<int> keep;

  CInterval *maxpos=NULL; //remembers intervals of maximum coverage

  int lastadded=0;
  for(int n=nstart;n<=nend;n++) if(pred[n]->strand==strand || pred[n]->strand=='.'){
	  if(n>nstart) {
		  if(equal_pred(pred,lastadded,n)) {
			  if(pred[n]->cov>pred[lastadded]->cov){
				  pred[lastadded]->flag=pred[n]->flag;
				  if(pred[lastadded]->exons[0].start != pred[n]->exons[0].start ||
						  pred[lastadded]->exons.Last().end!=pred[n]->exons.Last().end) { // new prediction has to replace old one
					  pred[lastadded]->tlen=pred[n]->tlen;
					  pred[lastadded]->exons[0].start=pred[n]->exons[0].start;
					  pred[lastadded]->exons.Last().end=pred[n]->exons.Last().end;
				  }
			  }
			  for(int j=0;j<pred[n]->exons.Count();j++) {
				  pred[lastadded]->exoncov[j]+=pred[n]->exoncov[j];
			  }
			  pred[lastadded]->cov+=pred[n]->cov;
			  pred[lastadded]->frag+=pred[n]->frag;

			  maxpos=add_pred_to_cov(maxpos,pred[lastadded]);
			  if(pred[n]->id && !pred[lastadded]->id) { pred[lastadded]->id=Gstrdup(pred[n]->id);}
			  continue;
		  }
	  }
	  bool abundant=true;
	  maxpos=add_pred_to_cov(maxpos,pred[n],&abundant);
	  //maxpos=add_pred_to_cov(maxpos,pred[n]);
	  if(pred[n]->id || abundant) {
		  keep.Add(n);
		  lastadded=n;
	  }
  }

  for(int i=0;i<keep.Count();i++) {
	  int n=keep[i];
	  if(is_pred_above_frac(maxpos,pred[n])) { // print this transcript
		  if(genes[pred[n]->geneno]==-1) genes[pred[n]->geneno]=++geneno;
		  transcripts[pred[n]->geneno]++;
		  if(pred[n]->flag) {
			  fprintf(f_out,"%d %d %.6f %.6f\n",pred[n]->exons.Count()+1,pred[n]->tlen,pred[n]->frag,pred[n]->cov);
			  fprintf(f_out,"%s\tStringTie\ttranscript\t%d\t%d\t1000\t%c\t.\tgene_id \"%s.%d\"; transcript_id \"%s.%d.%d\"; ",
					  refname.chars(),pred[n]->start,pred[n]->end,pred[n]->strand,label.chars(),genes[pred[n]->geneno],
					  label.chars(),genes[pred[n]->geneno],transcripts[pred[n]->geneno]);
			  if(pred[n]->id) fprintf(f_out,"reference_id \"%s\"; ",pred[n]->id);
			  fprintf(f_out,"cov \"%.6f\";\n",pred[n]->cov);
			  for(int j=0;j<pred[n]->exons.Count();j++) {
				  fprintf(f_out,"%s\tStringTie\texon\t%d\t%d\t1000\t%c\t.\tgene_id \"%s.%d\"; transcript_id \"%s.%d.%d\"; exon_number \"%d\"; ",
						  refname.chars(),pred[n]->exons[j].start,pred[n]->exons[j].end,pred[n]->strand,label.chars(),genes[pred[n]->geneno],
						  label.chars(),genes[pred[n]->geneno],transcripts[pred[n]->geneno],j+1); // maybe add exon coverage here
				  if(pred[n]->id) fprintf(f_out,"reference_id \"%s\"; ",pred[n]->id);
				  fprintf(f_out,"cov \"%.6f\";\n",pred[n]->exoncov[j]);
			  }
		  }
		  else pred[n]->flag=true;
	  }
	  else pred[n]->flag=false;
  }

  delete_interval(maxpos);

  return(geneno);
}

uint min(uint n1,uint n2) {
	if(n1<n2) return(n1);
	return(n2);
}

uint max(uint n1,uint n2) {
	if(n1<n2) return(n2);
	return n1;
}

bool included_pred(GPVec<CPrediction>& pred,int n1,int n2) { // check if the small prediction is included in the larger prediction

	if(pred[n1]->start > pred[n2]->end || pred[n2]->start>pred[n1]->end) return false;

	int big=n1;
	int small=n2;
	if(pred[n1]->exons.Count()<pred[n2]->exons.Count()) {
		big=n2;
		small=n1;
	}

	if(pred[small]->id) {
		if(pred[big]->id) if(strcmp(pred[small]->id,pred[big]->id)) return false;
		if(pred[small]->exons.Count()!=pred[big]->exons.Count()) return false;
	}

	int bex=0;
	while(bex<pred[big]->exons.Count()) {
		if(pred[small]->exons[0].start>pred[big]->exons[bex].end) bex++;
		else { // now pred[small]->exons[0].start <= pred[big]->exons[bex].end
			if(pred[small]->exons[0].end<pred[big]->exons[bex].start) return false; // no overlap
			int sex=0;
			while(sex<pred[small]->exons.Count() && bex<pred[big]->exons.Count()) {
				if(sex==pred[small]->exons.Count()-1) { // I am at end of small pred
					if(bex==pred[big]->exons.Count()-1) return true; // same intron structure and there is overlap
					// here not at last exon in big prediction
					if(pred[small]->exons[sex].end>=pred[big]->exons[bex+1].start) return false;
					return true;
				}
				// sex is not last exon in small prediction but overlaps bex
				if(bex==pred[big]->exons.Count()-1) return false; // small pred extends past big pred
				if(pred[small]->exons[sex].end != pred[big]->exons[bex].end) return false;
				bex++;
				sex++;
				if(pred[small]->exons[sex].start != pred[big]->exons[bex].start) return false;
			}
			return false;
		}
	}

	return false;
}

void update_cov(GPVec<CPrediction>& pred,int big,int small,float frac=1) {

	int bex=0;
	while(pred[small]->exons[0].start>pred[big]->exons[bex].end) bex++;

	if(!bex && pred[small]->exons.Count()>1 && pred[big]->exons[0].start<pred[small]->exons[0].start) { // adjust start to account for trimming
		pred[big]->tlen-=pred[small]->exons[0].start-pred[big]->exons[0].start;
		pred[big]->exons[0].start=pred[small]->exons[0].start;
	}

	int sex=0;
	int overlap=0;
	while(sex<pred[small]->exons.Count()) {
		int exovlp=(pred[small]->exons[sex].end<pred[big]->exons[bex].end ? pred[small]->exons[sex].end : pred[big]->exons[bex].end)-
				(pred[small]->exons[sex].start>pred[big]->exons[bex].start ? pred[small]->exons[sex].start : pred[big]->exons[bex].start)+1;

		if(bex==pred[big]->exons.Count()-1 && sex>=1 && pred[big]->exons[bex].end>pred[small]->exons[sex].end) { // adjust end
			pred[big]->tlen-=pred[big]->exons[bex].end-pred[small]->exons[sex].end;
			pred[big]->exons[bex].end=pred[small]->exons[sex].end;
		}

		pred[big]->exoncov[bex]=(pred[big]->exoncov[bex]*pred[big]->exons[bex].len()+frac*pred[small]->exoncov[sex]*exovlp)/pred[big]->exons[bex].len();
		overlap+=exovlp;
		sex++;bex++;
	}

	pred[big]->cov=(pred[big]->tlen*pred[big]->cov+overlap*frac*pred[small]->cov)/pred[big]->tlen;

}

int print_cluster(GPVec<CPrediction>& pred,GVec<int>& genes,GVec<int>& transcripts, int geneno,GStr& refname) {

	//fprintf(stderr,"start print cluster...\n");
	// sort predictions from the most abundant to the least:
	pred.Sort(predcovCmp);
	GVec<int> keep;

	CInterval *maxpos=NULL; //remembers intervals of maximum coverage

	for(int n=0;n<pred.Count();n++) { // don't need this anymore since I already took care of it before: if(pred[n]->strand==strand || pred[n]->strand=='.'){

		/*
		{ // DEBUG ONLY
			fprintf(stderr,"Consider prediction[%d] %c cov=%f:",n,pred[n]->strand,pred[n]->cov);
			for(int i=0;i<pred[n]->exons.Count();i++) fprintf(stderr," %d-%d",pred[n]->exons[i].start,pred[n]->exons[i].end);
			fprintf(stderr,"\n");
		}
		*/

		int k=0;
		bool included=false;
		while(!included && k<keep.Count() && keep[k]<n) {
			//if(pred[keep[k]]->exons.Count()>=pred[n]->exons.Count() && included_pred(pred,keep[k],n)) {

			if(included_pred(pred,keep[k],n)) {
				if(pred[keep[k]]->exons.Count()<pred[n]->exons.Count()) {
					//if(pred[keep[k]]->cov>pred[n]->cov) break; // this is new and improves specificity but I loose some things -> TO CHECK WHAT IT ACT: also this should always happen because of the sort procedure
					if(pred[keep[k]]->exons.Count()>2) break;
					update_cov(pred,n,keep[k]);
					pred[keep[k]]->cov=pred[n]->cov;
					pred[keep[k]]->exons.Clear();
					pred[keep[k]]->exons.Add(pred[n]->exons);
					pred[keep[k]]->exoncov.Clear();
					pred[keep[k]]->exoncov.Add(pred[n]->exoncov);
					pred[keep[k]]->flag=true;
					pred[keep[k]]->start=pred[n]->start;
					pred[keep[k]]->end=pred[n]->end;
					pred[keep[k]]->tlen=pred[n]->tlen;
				}
				else update_cov(pred,keep[k],n);

				if(pred[n]->id && !pred[keep[k]]->id) pred[keep[k]]->id=Gstrdup(pred[n]->id);
				pred[keep[k]]->frag+=pred[n]->frag;

				//fprintf(stderr,"...included in prediction[%d] with cov=%f\n",keep[k],pred[keep[k]]->cov);
				included=true;
				break;
			}
			k++;
		}
		if(included) {
			maxpos=add_pred_to_cov(maxpos,pred[keep[k]]);

			/*
			{ // DEBUG ONLY
				fprintf(stderr,"Maxpos is:");
				CInterval *interval=maxpos;
				while(interval!=NULL) {
					fprintf(stderr," pos=%d val=%f",interval->pos,interval->val);
					interval=interval->next;
				}
				fprintf(stderr,"\n");
			}
			*/

			continue;
		}
		bool abundant=true;
		maxpos=add_pred_to_cov(maxpos,pred[n],&abundant);
		//maxpos=add_pred_to_cov(maxpos,pred[n]);

		/*
		{ // DEBUG ONLY
			fprintf(stderr,"Maxpos is:");
			CInterval *interval=maxpos;
			while(interval!=NULL) {
				fprintf(stderr," pos=%d val=%f",interval->pos,interval->val);
				interval=interval->next;
			}
			fprintf(stderr,"\n");
		}
		*/

		if(pred[n]->id || abundant) {
			keep.Add(n);
			//fprintf(stderr,"...keep prediction\n");
		}
	}

  for(int i=0;i<keep.Count();i++) {
	  int n=keep[i];
	  if(is_pred_above_frac(maxpos,pred[n])) { // print this transcript

		  /*
		  { // DEBUG ONLY
			  fprintf(stderr,"print prediction %d",n);
			  if(pred[n]->flag) fprintf(stderr," with true flag");
			  fprintf(stderr,"\n");
		  }
		  */

		  if(pred[n]->flag) {
			  if(genes[pred[n]->geneno]==-1) genes[pred[n]->geneno]=++geneno;
			  transcripts[pred[n]->geneno]++;

			  fprintf(f_out,"%d %d %.6f %.6f\n",pred[n]->exons.Count()+1,pred[n]->tlen,pred[n]->frag,pred[n]->cov);
			  fprintf(f_out,"%s\tStringTie\ttranscript\t%d\t%d\t1000\t%c\t.\tgene_id \"%s.%d\"; transcript_id \"%s.%d.%d\"; ",
					  refname.chars(),pred[n]->start,pred[n]->end,pred[n]->strand,label.chars(),genes[pred[n]->geneno],
					  label.chars(),genes[pred[n]->geneno],transcripts[pred[n]->geneno]);
			  if(pred[n]->id) fprintf(f_out,"reference_id \"%s\"; ",pred[n]->id);
			  fprintf(f_out,"cov \"%.6f\";\n",pred[n]->cov);
			  for(int j=0;j<pred[n]->exons.Count();j++) {
				  fprintf(f_out,"%s\tStringTie\texon\t%d\t%d\t1000\t%c\t.\tgene_id \"%s.%d\"; transcript_id \"%s.%d.%d\"; exon_number \"%d\"; ",
						  refname.chars(),pred[n]->exons[j].start,pred[n]->exons[j].end,pred[n]->strand,label.chars(),genes[pred[n]->geneno],
						  label.chars(),genes[pred[n]->geneno],transcripts[pred[n]->geneno],j+1); // maybe add exon coverage here
				  if(pred[n]->id) fprintf(f_out,"reference_id \"%s\"; ",pred[n]->id);
				  fprintf(f_out,"cov \"%.6f\";\n",pred[n]->exoncov[j]);
			  }
		  }
		  else pred[n]->flag=true;
	  }
	  else pred[n]->flag=false;
  }

  delete_interval(maxpos);

  return(geneno);
}

int print_cluster_inclusion(GPVec<CPrediction>& pred,GVec<int>& genes,GVec<int>& transcripts, int geneno,GStr& refname, int limit=3) {

	//fprintf(stderr,"start print cluster...\n");
	// sort predictions from the one with the most exons to the one with the least:
	pred.Sort(predexCmp);

	GVec<int> included[pred.Count()];
	GVec<float> maxcov;
	GVec<float> totalcov;

	for(int n1=0;n1<pred.Count()-1;n1++) {
		float elem=0;
		maxcov.Add(elem);
		totalcov.Add(elem);
		bool equal=false;
		for(int n2=n1+1;n2<pred.Count();n2++)
			if(included_pred(pred,n1,n2)) {
				if(equal && pred[n1]->exons.Count()<pred[n2]->exons.Count()) break;
				if(pred[n1]->exons.Count()==pred[n2]->exons.Count()) {
					if(maxcov[n1]<pred[n1]->cov) maxcov[n1]=pred[n1]->cov;
					equal=true;
				}
				included[n1].Add(n2);
				if(pred[n2]->cov>maxcov[n1]) maxcov[n1]=pred[n2]->cov;
				totalcov[n1]+=pred[n2]->cov;
			}
	}

	GVec<int> keep;

	CInterval *maxpos=NULL; //remembers intervals of maximum coverage

	for(int n=0;n<pred.Count();n++) { // don't need this anymore since I already took care of it before: if(pred[n]->strand==strand || pred[n]->strand=='.'){

		/*
		{ // DEBUG ONLY
			fprintf(stderr,"Consider prediction[%d] %c cov=%f:",n,pred[n]->strand,pred[n]->cov);
			for(int i=0;i<pred[n]->exons.Count();i++) fprintf(stderr," %d-%d",pred[n]->exons[i].start,pred[n]->exons[i].end);
			fprintf(stderr," included in predictions:");
			for(int i=0;i<included[n].Count();i++) fprintf(stderr," %d",included[n][i]);
			if(n<pred.Count()-1) fprintf(stderr," maxcov=%f totalcov=%f",maxcov[n],totalcov[n]);
			fprintf(stderr,"\n");
		}
		*/

		if(included[n].Count() && (maxcov[n]>=pred[n]->cov || pred[n]->exons.Count()<limit)) { // this prediction is included in others, and is less abundant or has very few exons: make it <=2 if two exon genes also are to be ignored
			for(int k=0;k<included[n].Count();k++) {
				update_cov(pred,included[n][k],n,pred[included[n][k]]->cov/totalcov[n]);
			}
		}
		else {

			bool abundant=true;
			maxpos=add_pred_to_cov(maxpos,pred[n],&abundant);
			//maxpos=add_pred_to_cov(maxpos,pred[n]);

			/*
			{ // DEBUG ONLY
				fprintf(stderr,"Maxpos is:");
				CInterval *interval=maxpos;
				while(interval!=NULL) {
					fprintf(stderr," pos=%d val=%f",interval->pos,interval->val);
					interval=interval->next;
				}
				fprintf(stderr,"\n");
			}
			*/

			if(pred[n]->id || abundant) {
				keep.Add(n);
				//fprintf(stderr,"...keep prediction\n");
			}
		}
	}

	for(int i=0;i<keep.Count();i++) {
	  int n=keep[i];
	  if(is_pred_above_frac(maxpos,pred[n])) { // print this transcript

		  //fprintf(stderr,"print prediction %d",n);
		  //if(pred[n]->flag) fprintf(stderr," with true flag");
		  //fprintf(stderr,"\n");

		  if(pred[n]->flag) {
			  if(genes[pred[n]->geneno]==-1) genes[pred[n]->geneno]=++geneno;
			  transcripts[pred[n]->geneno]++;

			  fprintf(f_out,"%d %d %.6f %.6f\n",pred[n]->exons.Count()+1,pred[n]->tlen,pred[n]->frag,pred[n]->cov);
			  fprintf(f_out,"%s\tStringTie\ttranscript\t%d\t%d\t1000\t%c\t.\tgene_id \"%s.%d\"; transcript_id \"%s.%d.%d\"; ",
					  refname.chars(),pred[n]->start,pred[n]->end,pred[n]->strand,label.chars(),genes[pred[n]->geneno],
					  label.chars(),genes[pred[n]->geneno],transcripts[pred[n]->geneno]);
			  if(pred[n]->id) fprintf(f_out,"reference_id \"%s\"; ",pred[n]->id);
			  fprintf(f_out,"cov \"%.6f\";\n",pred[n]->cov);
			  for(int j=0;j<pred[n]->exons.Count();j++) {
				  fprintf(f_out,"%s\tStringTie\texon\t%d\t%d\t1000\t%c\t.\tgene_id \"%s.%d\"; transcript_id \"%s.%d.%d\"; exon_number \"%d\"; ",
						  refname.chars(),pred[n]->exons[j].start,pred[n]->exons[j].end,pred[n]->strand,label.chars(),genes[pred[n]->geneno],
						  label.chars(),genes[pred[n]->geneno],transcripts[pred[n]->geneno],j+1); // maybe add exon coverage here
				  if(pred[n]->id) fprintf(f_out,"reference_id \"%s\"; ",pred[n]->id);
				  fprintf(f_out,"cov \"%.6f\";\n",pred[n]->exoncov[j]);
			  }
		  }
		  else pred[n]->flag=true;
	  }
	  else pred[n]->flag=false;
  }

  delete_interval(maxpos);

  return(geneno);
}


int print_transcript_signcluster(char strand,GList<CPrediction>& pred,GVec<int>& genes,GVec<int>& transcripts,
		int nstart, int nend, int geneno,GStr& refname) {

  GVec<int> keep;

  float maxcov=0;
  int lastadded=0;
  for(int n=nstart;n<=nend;n++) if(pred[n]->strand==strand || pred[n]->strand=='.'){
	  if(n>nstart) {
		  if(equal_pred(pred,lastadded,n)) {
			  /*
			  fprintf(stderr,"pred %d is equal to pred %d btw %d-%d\n",lastadded,n,pred[lastadded]->start,pred[lastadded]->end);
			  for(int j=0;j<pred[n]->exons.Count();j++) {
				  fprintf(stderr,"%d-%d vs %d-%d\n",pred[lastadded]->exons[j].start,pred[lastadded]->exons[j].end,pred[n]->exons[j].start,pred[n]->exons[j].end);
			  }
			  */
		    if(pred[n]->cov>pred[lastadded]->cov){
		      pred[lastadded]->flag=pred[n]->flag;
		      if(pred[lastadded]->exons[0].start != pred[n]->exons[0].start ||
		    		  pred[lastadded]->exons.Last().end!=pred[n]->exons.Last().end) { // new prediction has to replace old one
		    	  pred[lastadded]->tlen=pred[n]->tlen;
		    	  pred[lastadded]->exons[0].start=pred[n]->exons[0].start;
		    	  pred[lastadded]->exons.Last().end=pred[n]->exons.Last().end;

		      }
		    }
		    for(int j=0;j<pred[n]->exons.Count();j++) {
		      pred[lastadded]->exoncov[j]+=pred[n]->exoncov[j];
		    }
		    pred[lastadded]->cov+=pred[n]->cov;
		    pred[lastadded]->frag+=pred[n]->frag;
		    if(pred[lastadded]->cov > maxcov) {
		      maxcov=pred[lastadded]->cov;
		    }
		    if(pred[n]->id && !pred[lastadded]->id) { pred[lastadded]->id=Gstrdup(pred[n]->id);}
		    continue;
		  }
	  }
	  if(pred[n]->cov>maxcov) {
		  maxcov=pred[n]->cov;
	  }
	  keep.Add(n);
	  lastadded=n;
  }

  //fprintf(stderr,"keepcount=%d\n",keep.Count());

  for(int i=0;i<keep.Count();i++) {
	  int n=keep[i];
	  if((pred[n]->exons.Count()==1 && pred[n]->cov>=maxcov) || (pred[n]->exons.Count()>1 && pred[n]->cov/maxcov>=isofrac))	{ // print this transcript
		  if(genes[pred[n]->geneno]==-1) genes[pred[n]->geneno]=++geneno;
		  transcripts[pred[n]->geneno]++;
		  if(pred[n]->flag) {
			  fprintf(f_out,"%d %d %.6f %.6f\n",pred[n]->exons.Count()+1,pred[n]->tlen,pred[n]->frag,pred[n]->cov);
			  fprintf(f_out,"%s\tStringTie\ttranscript\t%d\t%d\t1000\t%c\t.\tgene_id \"%s.%d\"; transcript_id \"%s.%d.%d\"; ",
				  refname.chars(),pred[n]->start,pred[n]->end,pred[n]->strand,label.chars(),genes[pred[n]->geneno],
				  label.chars(),genes[pred[n]->geneno],transcripts[pred[n]->geneno]);
			  if(pred[n]->id) fprintf(f_out,"reference_id \"%s\"; ",pred[n]->id);
			  fprintf(f_out,"cov \"%.6f\";\n",pred[n]->cov);
			  for(int j=0;j<pred[n]->exons.Count();j++) {
				  fprintf(f_out,"%s\tStringTie\texon\t%d\t%d\t1000\t%c\t.\tgene_id \"%s.%d\"; transcript_id \"%s.%d.%d\"; exon_number \"%d\"; ",
			  		 refname.chars(),pred[n]->exons[j].start,pred[n]->exons[j].end,pred[n]->strand,label.chars(),genes[pred[n]->geneno],
			  		 label.chars(),genes[pred[n]->geneno],transcripts[pred[n]->geneno],j+1); // maybe add exon coverage here
				  if(pred[n]->id) fprintf(f_out,"reference_id \"%s\"; ",pred[n]->id);
				  fprintf(f_out,"cov \"%.6f\";\n",pred[n]->exoncov[j]);
			  }
		  }
		  else pred[n]->flag=true;
	  }
	  else pred[n]->flag=false;
  }

  return(geneno);
}


int print_transcripts(GList<CPrediction>& pred,int ngenes, int geneno, GStr& refname) {

	// print transcripts including the necessary isoform fraction cleanings
	pred.setSorted(predCmp);
	//pred.setSorted(true);

	int npred=pred.Count();

	int currentstartpos=-1;
	uint currentendpos=0;
	int nstartpos=0;
	int nendpos=0;
	int currentstartneg=-1;
	uint currentendneg=0;
	int nstartneg=0;
	int nendneg=0;
	GVec<int> genes(true); // for each gene remembers it's geneno
	genes.Resize(ngenes,-1);
	GVec<int> transcripts(true); // for each gene remembers how many transcripts were printed
	transcripts.Resize(ngenes,0);

	GPVec<CPrediction> pospred(false);
	GPVec<CPrediction> negpred(false);

	for(int n=0;n<npred;n++) {

		if(pred[n]->strand=='.') pred[n]->flag=false; // only let it print if it passes threshold for printing

		if(pred[n]->strand=='+' || pred[n]->strand=='.') {
			if(pred[n]->start > currentendpos) { // begin new cluster
				// first print predictions I've seen so far
				if(currentstartpos>-1) { // I've seen a cluster before
					switch (sensitivitylevel) {
					case 0: geneno=print_transcript_signcluster('+',pred,genes,transcripts,nstartpos,nendpos,geneno,refname);break;
					case 1: geneno=print_cluster(pospred,genes,transcripts,geneno,refname);break;
					case 2: geneno=print_cluster_inclusion(pospred,genes,transcripts,geneno,refname);break;
					case 3: geneno=print_signcluster('+',pred,genes,transcripts,nstartpos,nendpos,geneno,refname);
					}
					pospred.Clear();
				}

				pospred.Add(pred[n]);

				currentstartpos=pred[n]->start;
				currentendpos=pred[n]->end;
				nstartpos=n;
				nendpos=n;
			}
			else {
				if(pred[n]->end > currentendpos) currentendpos=pred[n]->end;
				nendpos=n;
				pospred.Add(pred[n]);

			}
		}

		if(pred[n]->strand=='-' || pred[n]->strand=='.') {
			if(pred[n]->start > currentendneg) { // begin new cluster

				// first print predictions I've seen so far
				if(currentstartneg>-1) { // I've seen a cluster before
					switch (sensitivitylevel) {
					case 0: geneno=print_transcript_signcluster('-',pred,genes,transcripts,nstartneg,nendneg,geneno,refname);break;
					case 1: geneno=print_cluster(negpred,genes,transcripts,geneno,refname);break;
					case 2: geneno=print_cluster_inclusion(negpred,genes,transcripts,geneno,refname);break;
					case 3: geneno=print_signcluster('-',pred,genes,transcripts,nstartneg,nendneg,geneno,refname);break;
					}
					negpred.Clear();
				}

				negpred.Add(pred[n]);

				currentstartneg=pred[n]->start;
				currentendneg=pred[n]->end;
				nstartneg=n;
				nendneg=n;
			}
			else {
				negpred.Add(pred[n]);
				if(pred[n]->end > currentendneg) currentendneg=pred[n]->end;
				nendneg=n;
			}
		}
	}

	if(currentstartpos>-1) { // I've seen a cluster before
		switch (sensitivitylevel) {
		case 0: geneno=print_transcript_signcluster('+',pred,genes,transcripts,nstartpos,nendpos,geneno,refname);break;
		case 1: geneno=print_cluster(pospred,genes,transcripts,geneno,refname);break;
		case 2: geneno=print_cluster_inclusion(pospred,genes,transcripts,geneno,refname);break;
		case 3: geneno=print_signcluster('+',pred,genes,transcripts,nstartpos,nendpos,geneno,refname);
		}
		pospred.Clear();
	}

	if(currentstartneg>-1) { // I've seen a cluster before
		switch (sensitivitylevel) {
		case 0: geneno=print_transcript_signcluster('-',pred,genes,transcripts,nstartneg,nendneg,geneno,refname);break;
		case 1: geneno=print_cluster(negpred,genes,transcripts,geneno,refname);break;
		case 2: geneno=print_cluster_inclusion(negpred,genes,transcripts,geneno,refname);break;
		case 3: geneno=print_signcluster('-',pred,genes,transcripts,nstartneg,nendneg,geneno,refname);
		}
		negpred.Clear();
	}

	return(geneno);
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

void noMoreBundles() { //setter
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
		    this_thread::sleep_for(chrono::milliseconds(30));
		    waitMutex.lock();
		     areThreadsWaiting=(threadsWaiting>0);
		    waitMutex.unlock();
		    this_thread::sleep_for(chrono::milliseconds(30));
		  }
		} while (areThreadsWaiting); //paranoid check that all threads stopped waiting
#else
	  NoMoreBundles=true;
#endif
}

/*
void processBundle1stPass(BundleData* bundle) {
	// code executed on bundle data after 1st pass
	//bundle->readlist() is empty here also no pairs data are available
	//just splice sites and coverage info for the bundle
	//TODO: build the groups here (bundle->groups) so processRead() can
	//collapse reads efficiently in 2nd pass
	if (verbose) {
		printTime(stderr);
		GMessage(">bundle %s:%d-%d(%d) (%djs) begins 1st pass processing...\n",
				bundle->refseq.chars(), bundle->start, bundle->end,
				bundle->numreads, bundle->junction.Count());
	}
	// generate groups here, storing them in some bundle->groups data structure
	if (verbose) {
//	#ifndef NOTHREADS
//			GLockGuard<GFastMutex> lock(logMutex);
//	#endif
	  printTime(stderr);
	  GMessage("^bundle %s:%d-%d(%d) 1st pass done.\n",bundle->refseq.chars(),
	  		bundle->start, bundle->end, bundle->numreads);
	}
}

*/

/*
void processBundle(BundleData* bundle) {
	if (verbose) {
	#ifndef NOTHREADS
			GLockGuard<GFastMutex> lock(logMutex);
	#endif
		printTime(stderr);
		GMessage(">bundle %s:%d-%d(%d) (%djs) begins processing...\n",
				bundle->refseq.chars(), bundle->start, bundle->end, bundle->numreads, bundle->junction.Count());
#ifdef GMEMTRACE
		double vm,rsm;
		get_mem_usage(vm, rsm);
		GMessage(" memory usage: %6.1fMB\n",rsm/1024);
		if (rsm>maxMemRS) {
			maxMemRS=rsm;
			maxMemVM=vm;
			maxMemBundle.format("%s:%d-%d(%d)", bundle->refseq.chars(), bundle->start, bundle->end, bundle->numreads);
		}
#endif
	}
	int ngenes=infer_transcripts(bundle->start,bundle->readlist,
	        bundle->junction, bundle->keepguides, bundle->bpcov, bundle->pred);
	if (bundle->pred.Count()>0) {
#ifndef NOTHREADS
		GLockGuard<GFastMutex> lock(printMutex);
#endif
		GeneNo=print_transcripts(bundle->pred, ngenes, GeneNo, bundle->refseq);
	}
	if (verbose) {
	#ifndef NOTHREADS
			GLockGuard<GFastMutex> lock(logMutex);
	#endif
	  printTime(stderr);
	  GMessage("^bundle %s:%d-%d(%d) done (%d predicted transcripts).\n",bundle->refseq.chars(),
	  		bundle->start, bundle->end, bundle->numreads, bundle->pred.Count());
	}
	bundle->Clear(); //full clear (after the 2nd pass unless singlePass was requested)
#ifndef NOTHREADS
	dataMutex.lock();
	dataClear.Push(bundle->idx);
	dataMutex.unlock();
#endif
}

*/

void processBundle(BundleData* bundle) {
	if (verbose) {
	#ifndef NOTHREADS
			GLockGuard<GFastMutex> lock(logMutex);
	#endif
		printTime(stderr);
		GMessage(">bundle %s:%d-%d(%d) (%djs) loaded, begins processing...\n",
				bundle->refseq.chars(), bundle->start, bundle->end, bundle->numreads, bundle->junction.Count());
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
	int ngenes=infer_transcripts(bundle->start,bundle->readlist,
			bundle->junction, bundle->keepguides, bundle->bpcov, bundle->pred, fast | bundle->covSaturated);
	if (bundle->pred.Count()>0) {
#ifndef NOTHREADS
		GLockGuard<GFastMutex> lock(printMutex);
#endif
		GeneNo=print_transcripts(bundle->pred, ngenes, GeneNo, bundle->refseq);
	}
	if (verbose) {
	#ifndef NOTHREADS
			GLockGuard<GFastMutex> lock(logMutex);
	#endif
	  printTime(stderr);
	  GMessage("^bundle %s:%d-%d(%d) done (%d predicted transcripts).\n",bundle->refseq.chars(),
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
#ifndef NOTHREADS
	dataMutex.lock();
	dataClear.Push(bundle->idx);
	dataMutex.unlock();
#endif
}

#ifndef NOTHREADS

bool waitForThreads() {
	bool noneWaiting=true;
	DBGPRINT("##> waiting for workers to enter wait state..\n");
	while (noneWaiting) {
	  waitMutex.lock();
	  noneWaiting=(threadsWaiting<1);
	  waitMutex.unlock();
	  if (noneWaiting)
	    this_thread::sleep_for(chrono::milliseconds(30));
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
		waitMutex.unlock();
		haveBundles.wait(queueMutex);
		waitMutex.lock();
		 if (threadsWaiting>0) threadsWaiting--;
		waitMutex.unlock();
		DBGPRINT3("---->> Thread%d: awakened! (queue len=%d)\n",td.thread->get_id(),bundleQueue->Count());
		BundleData* readyBundle=NULL;
		if ((bundleWork & 0x02)!=0 && (readyBundle=bundleQueue->Pop())!=NULL) { //is bit 1 set?
			 //while ()!=NULL) {
				if (bundleQueue->Count()==0)
					 bundleWork &= ~(int)0x02; //clear bit 1 (queue is empty)
				Num_Fragments+=readyBundle->num_fragments;
				Frag_Len+=readyBundle->frag_len;
				queueMutex.unlock();
				processBundle(readyBundle);
				DBGPRINT2("---->> Thread%d processed bundle, now locking back queueMutex\n", td.thread->get_id());
				queueMutex.lock();
				DBGPRINT2("---->> Thread%d locked back queueMutex\n", td.thread->get_id());
			// }
		}
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

//prepare the next free bundle for loading
int waitForData(BundleData* bundles) {
	int bidx=-1;
	while (bidx<0) {
	  dataMutex.lock();
	  if (dataClear.Count()>0) {
	    bidx=dataClear.Pop();
	    bundles[bidx].status=BUNDLE_STATUS_LOADING;
	    dataMutex.unlock();
	    return bidx;
	    }
	  dataMutex.unlock();
	  this_thread::sleep_for(chrono::milliseconds(20));
	}
	return -1; // should NEVER happen
}

#endif
