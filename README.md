VariantBam: One-pass extraction of sequencing reads from a BAM file using cascading rules
=========================================================================================

**License:** [GNU GPLv3][license]


Description
-----------

VariantBam is a tool to extract interesting sequencing reads from a BAM file. To save money, 
diskspace and future I/O, one may not want to store an entire BAM on disk after all the relevant VCF, 
MAFs, etc have been created. Instead it would be more efficient to store only those read-pairs 
who intersect some region around the variant locations. Ideally, these regions could be manually inspected
or reanalyzed without having to keep the entire BAM.

VariantBam was developed as a one-pass solution for the needs of different NGS tools. For instance, after
initial processing, an indel tool might be interested in storing MAPQ > 0 reads in 30k candidate regions of interest.
A separate SNP tool might find a different 20k regions, while a structural variation tool 
might be interested in only discordant or high-quality clipped reads across the BAM 
(where high-quality means they are not clipped to do low Phred quality).  
VariantBam uses a series of rules defined on distinct regions in order to handle 
all of these situations with a single pass through the BAM.

VariantBam is implemented in C++ and relies on the BamTools [API] (Derek Barnett, (c) 2009) for BAM I/O. 
Note that the capabilities of the [BamTools] command line ``bamtools filter``  
may provide a solution more suitable your needs than VariantBam. Briefly, ``bamtools filter`` allows you to 
specify a set of filters in JSON format to be applied to a BAM. See the Bamtools documentation for more detail. 

The question then is under what situations would you use ``bamtools filter``, and when would you use VariantBam? 
Below are a list of scenarios that we feel addresses the different domains of the two tools:

> 1. Extract all MAPQ 0 reads from a BAM - Either tool (prefer ``bamtools filter``)
> 2. Extract all reads in read group A - ``bamtools filter``
> 3. Extract all reads with NM tag >= 4 - Either tool (prefer ``bamtools filter``)
> 4. Extract all reads with NM tag >= 4 in exons - VariantBam.
> 5. Remove all reads with insert sizes between 100 and 600 bp - VariantBam
> 6. Extract all reads and mate within 1000bp of a variant or set of genes - VariantBam
> 7. Extract only high-quality reads with N bases beyong phred score X - VariantBam
> 8. Reduce a BAM to only high quality reads around your MAFs, VCFs and BED files - VariantBam

Thus, the main additions of VariantBam are three-fold:

> 1. The ability to filter specifically on read clipping, orientation and insert size (all important for structural variation), while taking into account the per-base phred quality.
> 2. Use of interval trees to efficiently determine if a read or read mate overlaps a region.
> 3. The ability to provide different rules for different regions, and the ability to provide these regions as common variant files (VCF, MAF, BED)

Example
-------

We ran VariantBam like this:

```bash
options=(
    --input-bam         big.bam
    --output-bam        small.bam
    --rules-file        rules.vb
    --proc-regions-file small_chr1_mask.bed
)
variant ${options[*]}
```

Syntax
------

This section will describe the syntax used by VariantBam to specify the cascades of rules and regions 
applied to a BAM. Below is an example of a valid VariantBam script:

```bash
    ### this is a comment. The line code below defines filters to be applied to each region/rule
    region@WG
    rule@!hardclip;mapped;mapped_mate;isize[0,600];!mapq[10,100]
    rule@!hardclip;mapped;mapped_mate;clip[10,101]
```

##### Region

Let's look first at the ``region`` tag. The region@ keyword marks that what follows is a genomic region, 
which is either the keyword ``WG`` for whole genome, or a VCF, MAF, Callstats or BED file. Regions are 
treated such that they will include any read who overlaps it, even partially. Optionally,
you can specify that your region of interest is a bit bigger than is actually in the file. You can do this by "padding"
the regions around the sites. For example:

``region@myvcf.vcf;pad[1000]``

You can also state that the region applies to reads who don't necessarily overlap the region, but their pair-mate does.

``region@myvcf;pad[1000];mate``

Note that the syntax is such that you must specify the file immediately after the @, following by other options
in any order. 

##### Rules

The next two lines specify a pair of rules, marked with the ``rule@`` tag. 
The default rule is to include every read, and the conditions within the rule are to be  
thought of as exclusion criteria. Note that you can take the complement of a condition 
by prefixing with a ``!``. For example:

```bash
    # do not include hardclipped reads, reads with isize > 600, or reads with mapq between 10 and 100.
    rule@!hardclip;isize[0,600];!mapq[10,100]
    
    # an equivalent specification would be
    rule@!hardclip;mapped;!isize[601,250000000];mapq[0,9]``
```

VariantBam handles multiple rules in the following way. For each read, VariantBam 
will cycle through the rules within a region until the read satisfies a rule. When it 
does, it includes the read in the output and stops checking. The logic for the entire collection of 
rules is then as follows:

On a given rule line, the read must satisfy ALL conditions (logical AND)

Across different rules, the read nead only satisfy ONE rule (logical OR)

To illustrate this, note that there is a small discrepancy in the first rule of the above. In the BAM format, 
unmapped reads and reads with unmapped mates are given an insert size of 0. However, in the same rule 
a condition is described to keep all reads with insert sizes 0-600 inclusive. Recalling the AND logic
within a rule, VariantBam will exclude the read, because it fails the ``mapped`` criteria.

Below is another example which uses the ability of VariantBam to interpret VCFs and BED files,
and apply rules separately to them.

```bash
    ### declare that region is a VCF file with pads of 1000 on either side of the variant.
    ### use the "mate" keyword to specify that pairs whose mate falls in the region belong to this rule
    region@/home/unix/jwala/myvcf.vcf;mate;pad[1000]
    #### I want to keep all the reads (this the default). Ill be explicit with the "every" keyword
    rule@every
    #### A BED file which gives a list of exons. In here, I just want to keep "variant" reads
    region@/home/unix/jwala/myexonlist.bed 
    ## keep discordant reads
    rule@!isize[0,600];
    ## keep only unmapped reads and their mates
    rule@!mapped;!mapped_mate
    ## or keep if it is hardclipped
    rule@hardclip
    ## keep reads with a mismatch to reference, but with high mapq
    rule@nm[1,101];mapq[30,100]
    
```

##### Global

To reduce redundancy, you can also type a ``global@`` rule anywhere in the stack,
and it will append that rule to everything below. For example, to exclude hardclipped, duplicate, qcfail and 
supplementary reads in every region, you would do:

```bash
    global@!hardclip;!duplicate;!qcfail;!supplementary
    region@WG
    rule@!isize[0,600]
    rule@clip[10,101];mapq[1,60]
    region@myvcf.vcf
```

which is equivalent to

```bash
    region@WG
    rule@!isize[0,600];!hardclip;!duplicate;!qcfail;!supplementary
    rule@clip[10,101];mapq[1,60];!hardclip;!duplicate;!qcfail;!supplementary
    region@myvcf.vcf
    rule@!hardclip;!duplicate;!qcfail;!supplementary
```
	
The global tag will apply through all of the regions. If you want to reset it for everything, just add ``global@every`` 
back onto the stack.

To make things run a little faster, you can set the order so that the more inclusive regions / rules are first. This only
applies if there is an overlap among regions. This is because VariantBam will move down the list of regions
that apply to this read and stop as soon as it meets an inclusion criteria. I prefer to start with a whole-genome region / rule
set, and then add more fine-mapped regions later.

##### Command Line Script

The usual method of inputing rules is with a VariantBam script as a text file (passed to
VariantBam with the ``-r`` flag). However, sometimes it is useful to not have to write an intermediate
file and just feed rules directly in. In that case, just pass a string literal to the -r flag, and VariantBam
will parse directly. Just separate lines with either a new ``-r`` flag or with a ``%``. For instance, you might run something like the following:

```bash
variant -i big.bam -o small.bam -r 'global@!hardclip' -r 'region@WG%rule@!isize[0,600];%rule@clip[10,101];mapq[1,60]' -r 'region@myvcf.vcf'
```

Note the single quotes so that it is interpreted as a string literal in BASH.

Full list of available rules
----------------------------

```bash
    #RULE           #EXAMPLE             #DESCRIPTION OF EXAMPLE / FLAG
    nm              nm[0,4]              NM tag from BAM (number of mismatches). e.g. must be 0-4 inclusive
    isize           isize[100,500]       Insert size, where all insert sizes are converted to positive.
    len             len[80,101]          Length of the read following phred trimming
    clip            clip[0,5]            Number of clipped bases following phred trimming
    phred           phred[4,100]         Range of phred scores that are "quality" 
    duplicate       duplicate            Read must be marked as optical duplicate 
    supp            !supp                Read must be primary alignment
    qcfail          !qcfail              Read must note be marked as QC Fail
    fwd_strand      fwd_strand           Read must be mapped to forward strand
    rev_strand      rev_strand           Read must be mapped to reverse strand
    mate_fwd_strand mate_fwd_strand      Mate of read must be mapped to forward strand
    mate_rev_strand mate_rev_strand      Mate of read must be mapped to reverse strand  
    mapped          !mapped              Read must be unmapped
    mapped_mate     mapped_mate          Mate must be mapped
    ff              ff                   Read pair must have forward-forward orientation
    rr              ff                   Read pair must have reverse-reverse orientation
    fr              ff                   Read pair must have forward-reverse orientation (proper)
    rf              ff                   Read pair must have reverse-forward orientation
    ic              ic                   Read pair must have inter-chromosomal mapping
    discordant      discordant[100,600]  Shortcut for !isize[100,600] || rr || ff || rf || ic (!discordant gives "proper" pairs)
```

[license]: https://github.com/broadinstitute/variant-bam/blob/master/LICENSE

[BamTools]: https://raw.githubusercontent.com/wiki/pezmaster31/bamtools/Tutorial_Toolkit_BamTools-1.0.pdf

[API]: http://pezmaster31.github.io/bamtools/annotated.html
