#ifndef SNOWTOOLS_BAM_WALKER_H__
#define SNOWTOOLS_BAM_WALKER_H__

#include <cassert>

#include "MiniRules.h"
#include "GenomicRegion.h"

#include "reads.h"
#include "SnowUtils.h"

// Phred score transformations
inline int char2phred(char b) {
  uint8_t v = b;
  assert(v >= 33);
  return v - 33;
}

/////////////// 
// Hold read counts
//////////////
struct ReadCount {

  int keep = 0;
  int total = 0;
  
  int percent () const {
    int perc  = SnowUtils::percentCalc<int>(keep, total); 
    return perc;
  }

  string totalString() const {
    return SnowUtils::AddCommas<int>(total);
  }

  string keepString() const {
    return SnowUtils::AddCommas<int>(keep);
  }

};

/** Walk along a BAM or along BAM regions and stream in/out reads
 */
class BamWalker {

 public:

  /** Construct a new BamWalker for streaming data in and streaming
   * out to a new BAM file.
   *
   * This constructor will open the BAM file and read the header. 
   * It will also open the out BAM file.
   */
  BamWalker(const std::string& in, const std::string& out);

  /** Construct a new BamWalker for reading a BAM
   */
  BamWalker(const std::string& in);

  ~BamWalker() {
#ifdef HAVE_BAMTOOLS
    m_writer->Close();
    m_reader->Close();
    delete m_writer;
    delete m_reader;
#endif

#ifdef HAVE_HTSLIB
    bgzf_close(fp);
    bam_hdr_destroy(br);
    hts_itr_destroy(hts_itr);
    hts_idx_destroy(idx);
    sam_close(fop);
#endif
  }

  void saveAlignment(Read &r);

  struct timespec start;

  void printRuleCounts(unordered_map<string, size_t> &rm) const;
  
  /** Set a part of the BAM to walk.
   *
   * This will set the BAM pointer to the given region.
   * @param gp Location to point the BAM to
   */
  void setBamWalkerRegion(const GenomicRegion& gp);

  /** Set up multiple regions. Overwrites current regions. 
   * 
   * This will set the BAM pointer to the first element of the
   * input list.
   * @param grv Set of location to point BAM to
   */
  void setBamWalkerRegions(const GenomicRegionVector& grv);

  // create the index file for the output bam
  void MakeIndex();

  // print to stdout
  void printMessage(const ReadCount &rc_main, const Read &r) const;

  void rewind() { 
    #ifdef HAVE_BAMTOOLS
    m_reader->Rewind(); 
    #endif
  }
  
  /** Print out some basic info about this walker, 
   * including MiniRules
   */
  friend std::ostream& operator<<(std::ostream& out, const BamWalker& b);

  /** Open a BAM file for streaming in
   */
  bool OpenInBam(const std::string& bam);

  /** Open a BAM file for streaming out
   */
  bool OpenOutBam(const std::string& bam);

  /** Pass a set of MiniRules to the BAM so it
   * can decide which reads to keep when streaming in
   */
  void SetMiniRulesCollection(const std::string& rules);

  /** Retrieve the next read from the BAM.
   *
   * If a MiniRulesCollection is defined for this BAM
   * will grab the next valid read.
   * r Read to fill with data
   * ref String identifying which rule this read passed. Empty if not passed
   * @return true if the next read is available
   */
  bool GetNextRead(Read &r, std::string& ref);

  void WriteAlignment(Read &r);

  

 private:

  // point index to this region of bam
  bool __set_region(const GenomicRegion& gp);

  //open bam, true if success
  bool __open_BAM_for_reading();

  // open m_out, true if success
  bool __open_BAM_for_writing();
  
  std::string m_in;
  std::string m_out;

  MiniRulesCollection m_mr;

#ifdef HAVE_BAMTOOLS
  BamTools::BamReader * m_reader;
  BamTools::BamWriter * m_writer;
#endif

  GenomicRegionVector m_region;

  int m_verbose = 1;

  size_t m_reads_seen = 0;
  
  size_t m_reads_seen_valid = 0;

#ifdef HAVE_HTSLIB
  // hts
  BGZF * fp = 0;
  hts_idx_t * idx = 0; // hts_idx_load(bamfile.c_str(), HTS_FMT_BAI);
  hts_itr_t * hts_itr = 0; // sam_itr_queryi(idx, 3, 60000, 80000);
  bam_hdr_t * br = 0;

  samFile* fop = 0;
  //fp = sam_open(fn, mode);
#endif

};


#endif 


