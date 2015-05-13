#ifndef SWAP_GENOMIC_REGION_VECTOR_H__
#define SWAP_GENOMIC_REGION_VECTOR_H__

#include <vector>
#include <string>
#include <cstdlib>
#include <list>
#include <unordered_map>

#include "IntervalTree.h"
#include "GenomicRegion.h"

typedef EInterval<GenomicRegion> GenomicInterval;
typedef unordered_map<int, vector<GenomicInterval> > GenomicIntervalMap;
typedef EIntervalTree<GenomicRegion> GenomicIntervalTree;
typedef unordered_map<int, GenomicIntervalTree> GenomicIntervalTreeMap;
typedef vector<GenomicInterval> GenomicIntervalVector;

/** Class to store vector of intervals on the genome
 */
class GenomicRegionVector {

 public:

  /** Construct an empty GenomicRegionVector 
   */
  GenomicRegionVector() {};

  /** Construct from a plain vector of GenomicRegion objects
   */
  GenomicRegionVector(std::vector<GenomicRegion>& vec);

  /** Construct a GenomicRegionVector with overlapping intervals
   * 
   * @param width Desired bin width
   * @param ovlp Amount that the bins should overlap
   * @param gr GenomicRegion to divide into smaller overlapping bins
   */
  GenomicRegionVector(int width, int ovlp, const GenomicRegion &gr);

  /** Read in a MuTect call-stats file and adds to GenomicRegionVector object.
   *
   * Reads a MuTect call-stats file and imports only
   * lines with KEEP marked. 
   * @param file Path to call-stats file
   * @param pad Amount to pad intervals by
   */
  void readMuTect(const std::string &file, int pad = 0);

  /** Read in a BED file and adds to GenomicRegionVector object
   * @param file Path to BED file
   * @param pad Amount to pad intervals by
   */
  void readBEDfile(const std::string &file, int pad = 0);

  /** Read in a VCF file and adds to GenomicRegionVector object
   * @param file Path to BED file
   * @param pad Amount to pad intervals by
   */
  void readVCFfile(const std::string &file, int pad = 0);

  /** Read in a text file (can be gzipped) and add to GenomicRegionVector
   *
   * This function will automatically detect which file type is being input:
   * -- ends in .vcf -> readVCFfile
   * -- ends in .bed -> readBEDfile
   * -- header contains "MuTect" -> readMuTect
   * The values are appended to existing vector of GenomicRegion objects
   * @param file Text file to read and store intervals
   * @param pad Amount to pad the intervals by (calls GenomicRegion::pad)
   */
  void regionFileToGRV(const std::string &file, int pad = 0);

  /** Fill in the GenomicIntervalTreeMap stored in this object. 
   *
   * A GenomicIntervalTreeMap is an unordered_map of GenomicIntervalTrees for 
   * each chromosome. A GenomicIntervalTree is an interval tree on the ranges
   * defined by the genomic interval, with cargo set at the same GenomicRegion object.
   */
  void createTreeMap();
  
  /** Send the GenomicRegionVector to a BED file
  * @param file Output file
  */
  void sendToBED(const std::string file);
  
  /** Reduces the GenomicRegion objects to minimal set
   */
  void mergeOverlappingIntervals();

  /** Return the number of GenomicRegions stored 
   */
  size_t size() const { return m_grv.size(); }

  /** Add a new GenomicRegion to end
   */
  void add(const GenomicRegion& g) { m_grv.push_back(g); createTreeMap(); }

  /** Is this object empty?
   */
  bool empty() const { return m_grv.size(); }

  /** Clear out all of the GenomicRegion objects
   */
  void clear() { m_grv.clear(); m_tree.clear(); idx = 0; }

  /** Retrieve a GenomicRegion at given index. 
   * 
   * Note that this does not move the idx iterator, which is 
   * used to loop through all the regions. Throws an exception
   * if the index is out of bounds.
   * @return GenomicRegion pointed to by index i
   */
  const GenomicRegion& at(size_t i) const;

  /** Find overlaps between this vector and input GenomicRegion.
   *
   * Requires that the GenomicIntervalTreeMap have been created first
   * @param gr Region to test
   * @return Number of overlapping elements in this GenomicRegionVector
   */
  size_t findOverlapping(const GenomicRegion &gr);

  /** Add two GenomicRegionVector objects together
   */
  void concat(const GenomicRegionVector& g);

  /** Output the GenomicRegionVector to a BED format
   *
   * Currently just outputs chr   pos   pos2 with no header.
   * @return BED formated string reprsentation 
   */
  string sendToBED() const;

  /** Fill the next GenomicRegion object. 
   * @return false if add end of vector
   */
  bool getNextGenomicRegion(GenomicRegion& gr);

  /** Rewind the element pointer to the first GenomicRegion
   */
  void rewind() { idx = 0; }

 private:

  std::vector<GenomicRegion> m_grv;
  
  // always construct this object any time m_grv is modifed
  GenomicIntervalTreeMap m_tree;

  // index for current GenomicRegion
  size_t idx = 0;

};


#endif
