#!/usr/bin/perl
use Time::Local;
use Term::ANSIColor; 
$BEDTOOLS="/home/saurabha/UTILITIES/BEDTOOLS/BEDTOOLS_2/bin";
$sortBed="$BEDTOOLS/sortBed";
$intersectBed="$BEDTOOLS/intersectBed";
####################
system ("rm *.gtf *.SAF; ls -l; sleep 1; ls -l $FEATURE_COUNTS; sleep 1;");
$UCSC_ANNOTATION="/scratch/iwase-lab/GENOMES/MM10/gencode.vM21.annotation.gtf";
#########################################
system ("cat $UCSC_ANNOTATION |grep -v Rn45s|grep -v chrM |grep -v Lars2|grep -v Rpph1|grep -v Malat1|grep -v Eef1a1|grep -v Hsp90aa1|grep -v Hspa8|grep -v Rmrp|grep -v Actb|grep -v Hsp90ab1|grep -v Gapdh |grep -v Mir > REFGENE_MM10_WITHOUT_OUTLIERS.gtf");

system ("cat $UCSC_ANNOTATION |grep -v Rn45s|grep -v chrM |grep -v Mir > REFGENE_MM10_WITH_OUTLIERS.gtf");
#What exactly are the outliers and why do they need to be eliminated? (specific to one dataset?)
#highly expressed genes could skew the data/graphs in this count function
#could be worthwhile to run with and without outliers

######### CREATE SYMLINKS ###########
$FOLDER="/nfs/value/siwase2/BONEFAS/210928_P65cKO_maleandfemale/rsem_star_align";

@SAMPLE_INFO=qw(

P61D_HIPCTX:::P61D_HIPCTX.genome.sorted.bam
P61E_HIPCTX:::P61E_HIPCTX.genome.sorted.bam
P62D_HIPCTX:::P62D_HIPCTX.genome.sorted.bam
P62E_HIPCTX:::P62E_HIPCTX.genome.sorted.bam
P63D_HIPCTX:::P63D_HIPCTX.genome.sorted.bam
P63E_HIPCTX:::P63E_HIPCTX.genome.sorted.bam
P64D_HIPCTX:::P64D_HIPCTX.genome.sorted.bam
P64E_HIPCTX:::P64E_HIPCTX.genome.sorted.bam
P65D_HIPCTX:::P65D_HIPCTX.genome.sorted.bam
P65E_HIPCTX:::P65E_HIPCTX.genome.sorted.bam
P66D_HIPCTX:::P66D_HIPCTX.genome.sorted.bam
P66E_HIPCTX:::P66E_HIPCTX.genome.sorted.bam
P67D_HIPCTX:::P67D_HIPCTX.genome.sorted.bam
P68E_HIPCTX:::P68E_HIPCTX.genome.sorted.bam
P69D_HIPCTX:::P69D_HIPCTX.genome.sorted.bam
P69E_HIPCTX:::P69E_HIPCTX.genome.sorted.bam

);

############# CREATING SHORTCUTS TO BAM FILES WITH LABELS #############
for $x(0 .. $#SAMPLE_INFO)
{
($MARK,$FILE)=split(":::",$SAMPLE_INFO[$x]);
system ("ls -lh $FOLDER/$FILE"); sleep 1;
############# CREATING SHORTCUTS TO BAM FILES WITH LABELS #############
system ("ln -s  $FOLDER/$FILE $MARK.bam");

############# MAKING A LIST OF FILENAMES #############
$FILELIST=$FILELIST." $MARK.bam";
}
print "$FILELIST\n";

############# FEATURECOUNTS FOR READ ESTIMATION #############
$FEATURECOUNTS="/home/saurabha/UTILITIES/SUBREAD/subread-1.5.0-p1-source/bin/featureCounts";
$OPTIONS="-T 12 -p -S rf -a REFGENE_MM10_WITH_OUTLIERS.gtf -F GTF"; #12 threads, annotation file from above, format of annotation file
#add -s 1 for stranded feature 

$COMMAND="$FEATURECOUNTS $OPTIONS $FILELIST -o ms_P6brains_readcount_210929";
print "$COMMAND\n";
system ("$COMMAND");
sleep 3;

#removes extra columns 2-6 (chrom and strand info)
#the below line is perfect for input into DESEQ, but you can also run the other command to get the full output if you want
system("cat ms_P6brains_readcount_210929 |grep -v Program:featureCounts |cut -f 1,7-  > ms_P6brains_readcount_210929_FOR_DESEQ_withoutliers");
#system("cat SMCX_ALL_readcount |grep -v Program:featureCounts > SMCX_ALL_readcount_FOR_DESEQ");
print "Done !!!\nDone !!!\nDone !!!\n";
