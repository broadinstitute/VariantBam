#include "SnowTools/GenomicRegionVector.h"

#include <iostream>
#include "SnowTools/gzstream.h"
#include <sstream>
#include <cassert>

void GenomicRegionVector::readMuTect(const std::string &file, int pad) {
  
  assert(pad >= 0);

  std::cout << "Reading MuTect CallStats"  << std::endl;
  std::string curr_chr = "dum";
  
  igzstream iss(file.c_str());
  if (!iss || file.length() == 0) { 
    std::cerr << "MuTect call-stats file does not exist: " << file << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string line;
  while (std::getline(iss, line, '\n')) {
    size_t counter = 0;
      std::string chr, pos, judge;
      std::istringstream iss_line(line);
      std::string val;
      if (line.find("KEEP") != std::string::npos) {
	while(std::getline(iss_line, val, '\t')) {
	  switch (counter) { 
	  case 0 : chr = val; break; 
	  case 1 : pos = val; break;
	  }
	  if (counter >= 1)
	    break;
	  counter++;
	  
	  if (curr_chr != chr) {
	    std::cout << "...reading MuTect call-stats -- chr" << chr << std::endl;
	    curr_chr = chr;
	  }

	}
	if (GenomicRegion::chrToNumber(chr) >= 0) {
	  GenomicRegion gr(chr, pos, pos);
	  gr.pad(pad);
	  m_grv.push_back(gr);
	}
      } // end "keep" conditional
    } // end main while

  createTreeMap();
}

void GenomicRegionVector::readBEDfile(const std::string & file, int pad) {

  assert(pad >= 0);

  igzstream iss(file.c_str());
  if (!iss || file.length() == 0) { 
    std::cerr << "BED file does not exist: " << file << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string line;
  std::cout << "Reading normal BED" << std::endl;
  std::string curr_chr = "-1";
  while (std::getline(iss, line, '\n')) {

    size_t counter = 0;
    std::string chr, pos1, pos2;
    std::istringstream iss_line(line);
    std::string val;
    
    if (line.find("#") == std::string::npos) {
      while(std::getline(iss_line, val, '\t')) {
	switch (counter) { 
	case 0 : chr = SnowUtils::scrubString(val, "chr"); break; 
	case 1 : pos1 = val; break;
	case 2 : pos2 = val; break;
	}
	if (counter >= 2)
	  break;
	counter++;
	
	if (chr != curr_chr) {
	  //std::cout << "...reading from BED - chr" << chr << std::endl;
	  curr_chr = chr;
	}
	
      }
      if (GenomicRegion::chrToNumber(chr) >= 0) {
	GenomicRegion gr(chr, pos1, pos2);
	gr.pad(pad);
	m_grv.push_back(gr);
      }
    } // end "keep" conditional
  } // end main while

  createTreeMap();
  
}

void GenomicRegionVector::readVCFfile(const std::string & file, int pad) {

  assert(pad >= 0);

  igzstream iss(file.c_str());
  if (!iss || file.length() == 0) { 
    std::cerr << "VCF file does not exist: " << file << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string line;
  
  std::cout << "Parsing VCF file "  << std::endl;
  while (std::getline(iss, line, '\n')) {
    if (line.length() > 0) {
      if (line.at(0) != '#') { // its a valid line
	std::istringstream iss_this(line);
	int count = 0;
	std::string val, chr, pos;
	
	while (std::getline(iss_this, val, '\t')) {
	  switch (count) {
	  case 0 : chr = val;
	  case 1 : pos = val;
	  }
	  count++;
	}
	if (count < 3) {
	  std::cerr << "Didn't parse VCF line properly: " << line << std::endl;
	} else {
	  GenomicRegion gr(chr, pos, pos);
	  gr.pad(pad);
	  m_grv.push_back(gr);
	}
	
      }
    }
  }

  createTreeMap();
  
}

void GenomicRegionVector::regionFileToGRV(const std::string &file, int pad) {

  igzstream iss(file.c_str());
  if (!iss || file.length() == 0) { 
    std::cerr << "Region file does not exist: " << file << std::endl;
    exit(EXIT_FAILURE);
  }

  // get the header line to check format
  std::string header;
  if (!std::getline(iss, header, '\n'))
    std::cerr << "Region file is empty: " << file << std::endl;
  iss.close();

  GenomicRegionVector grv;

  // MUTECT CALL STATS
  if (header.find("MuTect") != std::string::npos)
    readMuTect(file, pad);
  // BED file
  else if (file.find(".bed") != std::string::npos)
    readBEDfile(file, pad);
  // VCF file
  else if (file.find(".vcf") != std::string::npos) 
    readVCFfile(file, pad);

}

// reduce a set of GenomicRegions into the minium overlapping set (same as GenomicRanges "reduce")
void GenomicRegionVector::mergeOverlappingIntervals() {

  // make the list
  std::list<GenomicRegion> intervals(m_grv.begin(), m_grv.end());

  intervals.sort();
  std::list<GenomicRegion>::iterator inext(intervals.begin());
  ++inext;
  for (std::list<GenomicRegion>::iterator i(intervals.begin()), iend(intervals.end()); inext != iend;) {
    if((i->pos2 > inext->pos1) && (i->chr == inext->chr))
      {
	if(i->pos2 >= inext->pos2) intervals.erase(inext++);
	else if(i->pos2 < inext->pos2)
	  { i->pos2 = inext->pos2; intervals.erase(inext++); }
      }
    else { ++i; ++inext; }
  }

  // move it over to a grv
  std::vector<GenomicRegion> v{ std::make_move_iterator(std::begin(intervals)), 
      std::make_move_iterator(std::end(intervals)) };
  m_grv = v;

  // recreate the tree map
  createTreeMap();

}

void GenomicRegionVector::createTreeMap() {

  // sort the 
  sort(m_grv.begin(), m_grv.end());

  GenomicIntervalMap map;
  for (auto it : m_grv) {
    map[it.chr].push_back(GenomicInterval(it.pos1, it.pos2, it));
  }

  for (auto it : map) 
    m_tree[it.first] = GenomicIntervalTree(it.second);

}

void GenomicRegionVector::sendToBED(const std::string file) {
  
  if (m_grv.size() ==  0) {
    std::cerr << "sendToBED: GenomicRegionVector is empty" << std::endl;
    return; 
  }

  std::ofstream ofile(file.c_str(), ios::out);
  for (auto it : m_grv)
    ofile << GenomicRegion::chrToString(it.chr) << "\t" << it.pos1 << "\t" << it.pos2 << std::endl;
  ofile.close();

}

// divide a region into pieces of width and overlaps
GenomicRegionVector::GenomicRegionVector(int width, int ovlp, const GenomicRegion &gr) {

  // undefined otherwise
  assert(width > ovlp);

  uint32_t start = gr.pos1;
  uint32_t end = gr.pos1 + width;

  // region is smaller than width
  if ( end >= gr.pos2 ) {
    std::cerr << "GenomicRegionVector constructor: GenomicRegion is smaller than bin width" << std::endl;
    return; 
  }

  // loop through the sizes until done
  while (end <= gr.pos2) {
    m_grv.push_back(GenomicRegion(gr.chr, start, end));
    end += width - ovlp; // make the new one
    start += width - ovlp;
  }
  assert(m_grv.size() > 0);
  
  // finish the last one
  start = m_grv.back().pos2 - width;
  end = gr.pos2;
  m_grv.push_back(GenomicRegion(gr.chr, start, end));

  createTreeMap();
}


size_t GenomicRegionVector::findOverlapping(const GenomicRegion &gr) {

  GenomicIntervalVector giv;
  GenomicIntervalTreeMap::iterator ff = m_tree.find(gr.chr);
  if (ff == m_tree.end())
    return 0;
  ff->second.findOverlapping(gr.pos1, gr.pos2, giv);
  return (giv.size());


}

string GenomicRegionVector::sendToBED() const {
  
  if (m_grv.size() ==  0)
    return ""; 

  stringstream ss;
  for (auto& i : m_grv)
    ss << i.chr << "\t" << i.pos1 << "\t" << i.pos2 << endl;

  return ss.str();

}

void GenomicRegionVector::concat(const GenomicRegionVector& g)
{
  m_grv.insert(m_grv.begin(), g.m_grv.begin(), g.m_grv.end());
  createTreeMap();
}

bool GenomicRegionVector::getNextGenomicRegion(GenomicRegion& gr)
{
  if (idx >= m_grv.size())
    return false;

  gr = m_grv[++idx];
  return true;
  
}

GenomicRegionVector::GenomicRegionVector(std::vector<GenomicRegion>& vec) : m_grv(vec)
{
  createTreeMap();
}

const GenomicRegion& GenomicRegionVector::at(size_t i) const
{ 
  if (i >= m_grv.size()) 
    throw 20;
  return m_grv[i]; 
}
