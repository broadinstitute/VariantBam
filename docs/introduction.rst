Introduction
------------

SnowmanSV is an algorithm to identify genomic rearrangments in whole genome
sequencing (WGS) data. The input is one or two (tumor/normal) BAM files, and the 
primary output is a set of rearrangement junctions. However, you may find the local 
assemblies produced by SnowmanSV to be useful for other applications. 

Structural variations, which include both copy-number alterations and copy-neutral
translocations, are a major source of variation in both germline and somatic tissues.
SnowmanSV identifies these variations by performing a "rolling" assembly across an
aligned BAM file, assembly short contigs (typically 200-800bp) from the sequencing reads.
The default behavior is to assembly only unmapped, clipped and discordant reads. 

The principal engine of SnowmanSV is based heavily on the String Graph Assemler (SGA)
software by Jared Simpson. To build SnowmanSV, SGA was modified to accept 
reads in small windows directly from BAM files, and to keep all processing
in memory. The default assembly window is 5000bp, with 500bp overlaps. 

SnowmanSV is divided into two modules:

1. ``snowman run``: Accesses the BAM files, performs the assembly, and sends contigs to BWA-MEM for alignment.

2. ``snowman gather``: Aligns sequencing reads to the contigs, and Filters the aligned contigs to identify high confidence SVs. 
	
In principal, one could stop after snowman run and perform your own analysis on the contigs file. 
Alternatively, any BAM file produced from a BWA-MEM alignment of an assembly to the reference could be input to snowman gather.
For example, we have found it useful to use snowman gather to analyze structural variations created from whole-genome
de novo assembly with the Discovar algorithm.

A manuscript for SnowmanSV is currently being prepared.
