#include "SnowTools/BamWalker.h"

// set the bam region
bool BamWalker::__set_region(const GenomicRegion& gp) {

  assert(m_in.length());

#ifdef HAVE_BAMTOOLS
  // set the region
  if (!m_reader->SetRegion(m_region.chr, m_region.pos1, m_region.chr, m_region.pos2)) {
    std::cerr << "Error: Failed to set region: " << gp << endl; 
    return false;
    //exit(EXIT_FAILURE);
  }
#endif 

#ifdef HAVE_HTSLIB
  //HTS set region
  if (!idx)
    idx = hts_idx_load(m_in.c_str(), HTS_FMT_BAI);
  hts_itr = sam_itr_queryi(idx, gp.chr, gp.pos1, gp.pos2);
  if (!hts_itr) {
    std::cerr << "Error: Failed to set region: " << gp << endl; 
    return false;
    //exit(EXIT_FAILURE);
  }
#endif

  return true;
}

void BamWalker::setBamWalkerRegion(const GenomicRegion& g)
{
  m_region.clear();
  m_region.add(g);
  __set_region(g);
}

void BamWalker::setBamWalkerRegions(const GenomicRegionVector& grv) 
{
  m_region = grv;
  m_region.rewind();
  __set_region(grv.at(0));
}

// closes the BamWriter and makes an index file
void BamWalker::MakeIndex() {

#ifdef HAVE_BAMTOOLS
  m_writer->Close();
  
  // open the file 
  BamReader reader;
  if (!reader.Open(m_out)) {
    cerr << "Error: Could not open the output BAM to create index " << m_out << endl;
    exit(EXIT_FAILURE);
  }

  // create the index
  if (!reader.CreateIndex()) {
    cerr << "Error: Could not create the output BAM index for " << m_out << endl;
    exit(EXIT_FAILURE);
  }

  reader.Close();
#endif
}

bool BamWalker::OpenInBam(const std::string& bam) 
{
  m_in = bam;

  return __open_BAM_for_reading();
}

bool BamWalker::OpenOutBam(const std::string& bam) 
{
  m_in = bam;

  return __open_BAM_for_writing();
}


BamWalker::BamWalker(const std::string& in) : m_in(in)
{
  // open for reading
  if (!__open_BAM_for_reading())
    throw 20;

}

bool BamWalker::__open_BAM_for_reading()
{

  assert(m_in.length());
  
#ifdef HAVE_BAMTOOLS
  // open the reader
  m_reader = new BamTools::BamReader();
  if (!m_reader->Open(in)) {
    cerr << "Error: Cannot open " << in << " for reading" << endl;
    exit(EXIT_FAILURE);
  }
#endif
  
#ifdef HAVE_HTSLIB
  // HTS open the reader
  const char rflag = 'r';
  fp = bgzf_open(m_in.c_str(), &rflag); 
  
  if (!fp) {
    cerr << "Error using HTS reader on opening " << m_in << endl;
    exit(EXIT_FAILURE);
  }
  br = bam_hdr_read(fp);
  
  if (!br) {
    cerr << "Error using HTS reader on opening " << m_in << endl;
    return false;
    //exit(EXIT_FAILURE);
  }
#endif
  
  return true;

}

bool BamWalker::__open_BAM_for_writing() 
{

  assert(m_out.length());

#ifdef HAVE_BAMTOOLS
  // open the writer
  m_writer = new BamTools::BamWriter();
  if (!m_writer->Open(m_out, m_reader->GetHeaderText(), m_reader->GetReferenceData())) {
    cerr << "Error: Cannot open BAM for writing " << m_out << endl;
    return false;
    //exit(EXIT_FAILURE);
  }

#endif

#ifdef HAVE_HTSLIB
  // hts open the writer
  fop = sam_open(m_out.c_str(), "wb");
  if (!fop) {
    cerr << "Error: Cannot open BAM for writing " << m_out << endl;
    return false;
    //exit(EXIT_FAILURE);
  }

  // hts write the header
  sam_hdr_write(fop, br);
#endif

  return true;

}

// this is also opens the files for reading/writing and checks
// that they are readable/writable
BamWalker::BamWalker(const std::string& in, const std::string& out) : m_in(in), m_out(out)
{

  // open for reading
  if (!__open_BAM_for_reading())
    throw 20;

  // open for writing
  if (!__open_BAM_for_writing())
    throw 20;

}

void BamWalker::SetMiniRulesCollection(const std::string& rules)
{
  // construct the minirules
  m_mr = MiniRulesCollection(rules);

  // check that it worked
  if (!m_mr.size()) {
    std::cerr << "No MiniRules were successfully parsed" << std::endl;
    throw 20;
  }
}

void BamWalker::printMessage(const ReadCount &rc_main, const Read &r) const {

  char buffer[100];
  string posstring = SnowUtils::AddCommas<int>(r_pos(r));
  sprintf (buffer, "Reading read %11s at position %2s:%-11s. Kept %11s (%2d%%) [running count across whole BAM]",  
	   rc_main.totalString().c_str(), GenomicRegion::chrToString(r_id(r)).c_str(), posstring.c_str(),  
	   rc_main.keepString().c_str(), rc_main.percent());
  
  printf ("%s | ",buffer);
  SnowUtils::displayRuntime(start);
  cout << endl;
  
}

void BamWalker::printRuleCounts(unordered_map<string, size_t> &rm) const {

  size_t total = 0;
  for (auto& i : rm)
    total += i.second;
  for (auto& i : rm) {
    cout << "  " << i.first << ":" << i.second << "(" << SnowUtils::percentCalc<size_t>(i.second, total) << "%)" << endl;
  }
  
}

bool BamWalker::GetNextRead(Read& r, std::string& rule)
{
  
#ifdef HAVE_HTSLIB

  void* dum = 0;
  bam1_t* b = bam_init1(); 
  if (hts_itr == 0) { 
    if (bam_read1(fp, b) < 0) { 
      bam_destroy1(b); 
      return false;
    } 
  } 
  
  bool valid = hts_itr_next(fp, hts_itr, b, dum) <= 0;
  if (!valid) { // read not found
    do {
      // try next region, return if no others to try
      GenomicRegion next_region; 
      if (!m_region.getNextGenomicRegion(next_region)) {
	bam_destroy1(b);
	return false;
      }
      // next region exists, try it
      __set_region(next_region);
      valid = hts_itr_next(fp, hts_itr, b, dum) <= 0;
    } while (!valid); // keep trying regions until works
  }
  
  r = std::shared_ptr<bam1_t> (b, free_delete());
#endif

#ifdef HAVE_BAMTOOLS
  BamTools::BamAlignment a; 
  if (!reader->GetNextAlignment(a)) 
    return false; 
  r = std::shared_ptr<BamTools::BamAlignment>(new BamTools::BamAlignment(a)) ;
#endif

  // check if it passed the rules
  rule = m_mr.isValid(r);

  return true;
}

void BamWalker::WriteAlignment(Read &r)
{
#ifdef HAVE_HTSLIB
  sam_write1(fop, br, r.get());
#endif
}

// save alignment
void BamWalker::saveAlignment(Read &r) {

#ifdef HAVE_HTSLIB
    sam_write1(fop, br, r.get()); 
#endif

#ifdef HAVE_BAMTOOLS
    m_writer->SaveAlignment(*(r.get()));
#endif

}

std::ostream& operator<<(std::ostream& out, const BamWalker& b)
{
  out << " -- In Bam:  " << b.m_in << std::endl;
  out << " -- Out Bam: " << b.m_out << std::endl;
  out << b.m_mr << std::endl;
  return out;
}
