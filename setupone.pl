#!/usr/bin/perl
use Getopt::Long;

my $samtools = '/opt/samtools/bin/samtools';
my $bcftools = '/opt/samtools/bin/bcftools';
my $gatk = 'module load dev/java/jdk1.7;java -Xmx4g -jar /groups/kishony/baym/gatklocal/GenomeAnalysisTK.jar';
my $sickledir = '/home/mhb15/illumina_pipeline/sickle-master'; #sickle directory
my $filtercom = "$sickledir/sickle pe -t sanger -q 20 -l 50 -x -n";
my $aligncom = '/opt/bowtie2/bowtie2 -p 2 -X 1000 --no-mixed --dovetail --very-sensitive --n-ceil 0,0.001'; #bowtiecommand
my $cutadaptcom = 'cutadapt -a CTGTCTCTTAT';

my $source = '.'; #source filepath
my $target = '/hms/scratch1/baym/experiment1/Sample001/'; #default filepath
my $refdir = '/groups/kishony/Reference_Genomes/EcoliMG1655'; #reference directory
my $protname = 'pairedendkeio';
my $runtag = 'pipeline';
my $debug;

GetOptions ('s|source=s' => \@source,
			't|target=s' => \$target, #(the directory with sample name last)
			'f|filter=s' => \$filtercom,
			'a|align=s' => \$aligncom,
			'r|refdir=s' => \$refdir,
			'p|protname=s' => \$protname,
			'runtag=s' => \$runtag,
			'copyfirst' => \$copyfirst,
			'debug' => \$debug);
			
			
sub runcommand { #this lets me switch from testing to running easily!
if ($debug) {
	print $_[0]."\n";
} else {
	print $_[0]."\n";
	return system $_[0];
}
}
#--sam-rg "ID:'.$source.'-illumina" -- currently not working
#$target =~ /\/([^\/]+)\/$/; #pulls the directory name out to find the filename for now -> replace later
#$readbase = $target.$1."_";


$reffasta = "$refdir/genome.fasta";
$fromhereflag = 1;

#make the target directory
$target =~ m/\/([^\/]+)\/$/;
$filename = $1;
runcommand "mkdir -p $target";


$numfilesin = @source;



print "\nCutting adaptors... (at ".localtime().")\n";#trim the reads

unless (-e "$target/done.cut" &&  $fromhereflag) {$fromhereflag=0;
	$success = 0;
	for $file (@source) {
		$read1in= $file."1.fastq";
		$read2in= $file."2.fastq";
	
		$success += runcommand "$cutadaptcom $read1in >> $target"."readcut1.fastq";
		$success += runcommand "$cutadaptcom $read2in >> $target"."readcut2.fastq";
		runcommand "touch $target/done.cut"  if ($success == 0);
	}
}


print "\nTrimming reads... (at ".localtime().")\n";#trim the reads
$trimtarget=$target."sickle2051";
runcommand "mkdir -p $trimtarget";
#don't copy the fastqs to the target base directory -- get them from the source!
$fr1="$trimtarget/filter_reads_1.fastq";
$fr2="$trimtarget/filter_reads_2.fastq";
unless (-e "$trimtarget/done.trim" &&  $fromhereflag) {$fromhereflag=0;

	#$success = runcommand "$sickledir/sickle pe -f $read1in -r $read2in -t sanger -o $fr1 -p $fr2 -s $trimtarget/singles.fastq -q 20 -l 50 -x -n";
	$success = runcommand "$filtercom -f $target/readcut1.fastq -r $target/readcut2.fastq -o $fr1 -p $fr2 -s $trimtarget/singles.fastq";
	runcommand "touch $trimtarget/done.trim"  if ($success == 0);
}


print "\nAligning reads... (at ".localtime().")\n";
#align the reads, using four processors
$aligndir=$trimtarget."/".$protname;
runcommand "mkdir -p $aligndir";
$alignsam = "$aligndir/aligned.sam";
unless (-e "$aligndir/done.align" &&  $fromhereflag) {$fromhereflag=0;
	$success = runcommand "$aligncom --un-conc $aligndir/unaligned.fastq -x $refdir/genome_bowtie2 -1 $fr1 -2 $fr2 -S $alignsam";
	runcommand "touch $aligndir/done.align" if ($success == 0);
}

print "\nMaking bam files... (at ".localtime().")\n";
$alignbam = "$aligndir/aligned.bam";
$alignsort = "$aligndir/aligned.sorted";
unless ((-e "$aligndir/done.sambam") &&  ($fromhereflag)) {$fromhereflag=0;
	runcommand "$samtools view -bS $alignsam > $alignbam";						#make bam file
	runcommand "$samtools sort $alignbam $alignsort";							#make sorted bam file
	$success = runcommand "$samtools index $alignsort.bam $alignsort.bai";		#index sorted bam file
	runcommand "touch $aligndir/done.sambam" if ($success == 0);
} 


#use picard to add rungroups so gatk doesn't bork
print "\nAdding read groups with picard... (at ".localtime().")\n";
$alignwrg = "$aligndir/alignedwrg.bam";
$alignwrgindex = "$aligndir/alignedwrg.bai";
$picardcommand= "java -jar /opt/picard/AddOrReplaceReadGroups.jar";
unless (-e "$aligndir/done.picard" &&  $fromhereflag) {$fromhereflag=0;
	runcommand "$picardcommand INPUT=$alignsort.bam OUTPUT=$alignwrg SORT_ORDER=coordinate RGLB=MG1655 RGPL=Illumina RGPU=2 RGSM=$filename";
	$success = runcommand "$samtools index $alignwrg $alignwrgindex";
	runcommand "touch $aligndir/done.picard" if ($success == 0);
}

#realign with GATK
print "\nRealigning with GATK... (at ".localtime().")\n";
runcommand "module load dev/java/jdk1.7";
$realigned = "$aligndir/realigned.bam";
$intervals = "$aligndir/forIndelRealigner.intervals";
unless (-e "$aligndir/done.gatk" &&  $fromhereflag) {$fromhereflag=0;
	runcommand "$gatk -T RealignerTargetCreator -R $reffasta -I $alignwrg -o $intervals";
	$success = runcommand "$gatk -T IndelRealigner -LOD 1.0 --entropyThreshold 0.10 -R $reffasta -I $alignwrg -targetIntervals $intervals -o $realigned";
	runcommand "touch $aligndir/done.gatk" if ($success == 0);
}



#use gatk to recalibrate quality scores and remake the outputbam
#use gatk recalibrater to recalibrate == note, without bootstrapping this has only been tested for Michael's data


# print "\nRecalibrating Base Score Quality with GATK... (at ".localtime().")\n";
# 
 $target =~ m/\/([^\/]+)\/$filename\/$/;
 $basedir = $1;
# 
 $recaldir = "$aligndir/recalibrated";
# $recalbam = "$recaldir/qualrecalibrated.bam";
 runcommand "mkdir -p $recaldir";
# unless (-e "$recaldir/done.bqsr" &&  $fromhereflag) {$fromhereflag=0;
# 	runcommand "$gatk -T BaseRecalibrator -R $reffasta -I $realigned -o $recaldir/recal_data.table --run_without_dbsnp_potentially_ruining_quality";
# 	$success = runcommand "$gatk -T PrintReads -R $reffasta -I $realigned -BQSR $recaldir/recal_data.table -o $recalbam";
# 	runcommand "touch $recaldir/done.bqsr" if ($success == 0);
# }
$recalbam = $realigned;

print "\nMaking the pileup in a separate process... (at ".localtime().")\n";
$pileup = "$recaldir/strain.pileup";
$pileupouts = "$filename_pileups";
runcommand "mkdir -p $basedir/runlogs/";

$runlogbase = $basedir."/runlogs/".$pileupouts;
unless (-e "$recaldir/done.pileup" && $fromhereflag) {$fromhereflag=0;
	$pileupcommand = "$samtools mpileup -q30 -s -O -d3000 -C30 -f $reffasta $recalbam > $pileup";
	$bsubcomm = "bsub -J ".$runtag."_".int(rand(5000))."_pileup -W 3:0 -o ".$runlogbase.".out -e ".$runlogbase.".err -q short ";
	$pileupdone = "touch $recaldir/done.pileup";
	runcommand $bsubcomm.'"'.$pileupcommand.";".$pileupdone.'"';
}



#### EVERYTHING BELOW HERE IS OBSOLETE


#print "\nMaking the variant vcf in a separate process... (at ".localtime().")\n";

#use GATK to make strain and variant VCFs

#$samcommand1 = "$samtools mpileup -q30 -s -d3000 -C30 -ugf $refdir/genome.fasta $realigned  > $aligndir/strain";
#$samcommand2 = "$bcftools view -g $aligndir/strain > $aligndir/strain.vcf";
#$samcommand3 = "$bcftools view -vS $aligndir/strain.vcf > $aligndir/variant.vcf";
#$vcfouts = "$filename_varvcf";
#$runlogbase = $basedir."/runlogs/".$vcfouts;
#$bsubcomm = "bsub -J ".$runtag."_".int(rand(5000))."_variant_vcf -W 1:55 -o ".$runlogbase.".out -e ".$runlogbase.".err -q short";
#$vcfcommand1 = "$bsubcomm  \'$gatk -T HaplotypeCaller -R $reffasta -I $recalbam -o $recaldir/gatk.output.raw.recal.snps.indels.vcf\'";
#$vcfcommand2 = "$gatk -T HaplotypeCaller -R $reffasta -I $recalbam -o $recaldir/gatk.output.raw.all.vcf --emitRefConfidence BP_RESOLUTION";
#$vcfcommand1 = "$bsubcomm  \'$gatk -T UnifiedGenotyper -R $reffasta -I $recalbam -o $recaldir/gatk.output.raw.recal.snps.indels.vcf  -ploidy 1 -glm BOTH -nt 4\'";
#$vcfcommand2 = "$gatk -T UnifiedGenotyper -R $reffasta -I $recalbam -o $recaldir/gatk.output.raw.all.vcf  -ploidy 1 -glm BOTH -out_mode EMIT_ALL_SITES -nt 4";



#unless (-e "$recaldir/done.vcfs" && $fromhereflag && !(-z "$recaldir/gatk.output.raw.recal.snps.indels.vcf")) {$fromhereflag=0;
#	runcommand $vcfcommand1;
#	print "\nMaking the strain vcf... (at ".localtime().")\n";
#	$success = runcommand $vcfcommand2;
#	runcommand "touch $recaldir/done.vcfs" if ($success == 0);
#}
#$vcfdone = "touch $aligndir/done.vcfs";
#	runcommand $samcommand1;
#	runcommand $samcommand2;
#	$success = runcommand $samcommand3;
#	runcommand $vcfdone if ($success == 0);

#print "\nVCFs finished at (at ".localtime().")\n";


#
#print "\nCleaning up... \n";
#if ((-e "$aligndir/done.pileup")) {
#	if (-e "$aligndir/done.vcfs") {
#		runcommand "rm $aligndir/aligned.sam; rm $aligndir/strain";
#	} else {
#		runcommand "bsub -w \"done($target_vcf)\" -W 2 -q mini \"rm $aligndir/aligned.sam; rm $aligndir/strain\"";
#	}
#}