#!/usr/bin/perl
#usage 'perl marginalia.pl -o outdir -q queue';
use Getopt::Long;

my $startdir = '.';
my $refsource = '/groups/kishony/Reference_Genomes/';
my $outdir = '.';
my $marginalia = 'marginalia';
my $queue = 'short -n 2';
my $runtag = 'pipeline'.int(rand(1000));
my $samplesheet = "";
my $gatkbase = "module load dev/java/jdk1.7;java -Xmx8g -jar /groups/kishony/baym/gatklocal/GenomeAnalysisTK.jar -T HaplotypeCaller -nct 2 ";

GetOptions ('s|startdir=s' => \$startdir,	#directory with reads split into folders
			'q|queue=s' => \$queue,			#what orchestra queue to use?
			'r|refsource=s' => \$refsource, #where is the reference genome?
			'm|marginalia=s' => \$marginalia, #where are the marginalia scripts?
			'o|outdir=s' => \$outdir,		#to what directory should analyses be outputted?
			'runtag=s' => \$runtag,			#unique tag for these processes? (to manage bsub dependencies)
			'unixcsv' => \$unixcsv,			#are samples.csv and alignment_params.csv unix or windows?
			'postprocess' => \$postprocess,			#are samples.csv and alignment_params.csv unix or windows?
			'preprocess' => \$preprocess,			#are samples.csv and alignment_params.csv unix or windows?
			'samples=s' => \$samplesheet,	#use an alternate sample sheet name
			'debug' => \$debug);			#debug mode
			
			
sub runcommand { #this lets me switch from testing to running easily!
if ($debug) {
	print $_[0]."\n";
} else {
	return system $_[0];
}
}
			
#read in the csvs
{local $/ = "\r" unless $unixcsv; #assume it's an excel csv unless the unixcsv flag is on
if ($samplesheet eq "") {
	open SAMPLES, "$startdir/samples.csv" or die "Can't find $startdir/samples.csv\n";
} else {
	open SAMPLES, "$startdir/$samplesheet" or die "Can't find $startdir/$samplesheet\n";
	
}
$sampline = <SAMPLES>; #chuck the top line
while ($sampline = <SAMPLES>) {
	$sampnum = $.-2;
	chomp $sampline;
	@sampdata = split(",",$sampline);
	$readsdir[$sampnum]=$sampdata[0];
	$readsprefix[$sampnum]=$sampdata[5];
	$sampleoutname[$sampnum]=$sampdata[3];
	$protocol[$sampnum]=$sampdata[4];
}
close SAMPLES;

$maxsamp = $sampnum;

foreach $sampnum (0..$maxsamp) {
	push @{ $memberhash{$sampleoutname[$sampnum]} }, $sampnum;
}


open ALIGNPARAMS, "$startdir/alignment_params.csv" or die "Can't find $startdir/alignment_params.csv\n";
$paramline = <ALIGNPARAMS>;
while ($paramline = <ALIGNPARAMS>) {
	@paramdata = split(",",$paramline);
	$refgenome{$paramdata[0]} = $paramdata[2];
	#print "$refgenome{$paramdata[0]}"
}
close ALIGNPARAMS;
}

#first dispatch the preanalyzer
###set up the commands
runcommand "mkdir -p $outdir";
runcommand "mkdir -p $outdir/runlogs";

unless ($postprocess) {
	#for each one we need to do a few things:
	#compile the perl command
	$commandbase = "perl $marginalia/setupone.pl";
	$perlcom = "";
	foreach $sampnum (0..$maxsamp) {
		if ($memberhash{$sampleoutname[$sampnum]}[0]==$sampnum) {
			$perlcom = "";
			foreach $localsamp (@{$memberhash{$sampleoutname[$sampnum]}}) {
			#now need to do this for each one
				$readbase = $readsdir[$localsamp].$readsprefix[$localsamp];
				chomp $readbase;
				if (-e $readbase."R1.fastq") {
					$readbase = $readbase."R";
				} elsif (-e $readbase."_R1.fastq") {
					$readbase = $readbase."_R";
				} else {
					#defaulting to this behavior
					$readbase = $readbase."_R";
					warn "Reads not found at $readbase\n" unless $debug;
				}

				$perlcom = $perlcom." -s $readbase";
			}
		
			$perlcom = $commandbase.$perlcom;
		
											#reads directory
			$perlcom = $perlcom." -t $outdir".'/'.$sampleoutname[$sampnum].'/'; 	#target directory
			$perlcom = $perlcom." -r $refsource$refgenome{$protocol[$sampnum]}";	#reference directory
			$perlcom = $perlcom." -p $protocol[$sampnum]";							#protocol name
			$perlcom = $perlcom." -runtag $runtag";									#unique tag for this run

			#compile the bsub command
			$logbase = "$outdir/runlogs/$sampleoutname[$sampnum]";
			$bsubcommand = "bsub -o $logbase.out -e $logbase.err -q $queue -W 6:0 -J ".$runtag."_".$sampleoutname[$sampnum]."_preprocess";
			runcommand "$bsubcommand $perlcom";
		
		}
	}
}

die if $preprocess;



runcommand "mkdir -p $outdir/analyses";


#now dispatch vcfcollate, but make it wait on the vcf processes (<1 second)
#and dispatch the genbank processor (5 seconds)
	$waitcom = "";
	$waitcom = '-w "done('.$runtag.'_*_preprocess)"' unless $postprocess;
	#$filecom = '"'.$filecom.')"';
	#print $filecom."\n";
	
	
$reffasta = $refsource.$refgenome{$protocol[$sampnum]}."/genome.fasta";

#### THIS ALL NEEDS TO BE UNCOMMENTED FOR ORCHESTRA
open GENOMEDICT, $refsource.$refgenome{$protocol[$sampnum]}."/genome.dict" or die "Can't find $refsource$refgenome{$protocol[$sampnum]}/genome.dict\n";
$topline = <GENOMEDICT>;
$nextline = <GENOMEDICT>;
chomp $nextline;
$nextline =~ m/SN:(\S+)\s*LN:(\d+)/;
$contigname = $1; #print $contigname."\n";
$genomesize = $2; #print $genomesize."\n";

close GENOMEDICT;
#$genomesize=50;

$logdir = "$outdir/runlogs";
$outdir2 = "$outdir/analyses/gatktranches";

$tranchsize = int($genomesize/$maxsamp/5);
$numintervals = int($genomesize/$tranchsize);

#print "$genomesize,$tranchsize,$numintervals\n";



runcommand "mkdir -p $outdir2";

	#for each one we need to do a few things:
	#compile the perl command
	$bsubbase = "bsub -q $queue -W 11:58 ";
	#$perlcom = $perlcom." -R $refsource$refgenome{$protocol[0]}/genome.fasta";	#reference directory
	$gatkbase = $gatkbase." -R $reffasta";	#reference directory
	
	foreach $sampnum (0..$maxsamp) {
		if ($memberhash{$sampleoutname[$sampnum]}[0]==$sampnum)  {
		$gatkbase = $gatkbase." -I $outdir".'/'.$sampleoutname[$sampnum]."/sickle2051/$protocol[$sampnum]/realigned.bam"; 	#target directory	
		}	
	}
	
runcommand "mkdir -p $logdir/gatktranches";
	for (0..$numintervals) {
	### this is to add the files per tranch
	#"-o multincttranch.out -e /groups/kishony/baym/maniallgatknobqsr/runlogs/multincttranch.err"
		$tranchstart = $_*$tranchsize+1;
		$tranchfinish = ($_+1)*$tranchsize;
		$tranchfinish = $genomesize if $tranchfinish>$genomesize;
		#print $contigname.":".$tranchstart."-".$tranchfinish."\n";
		next if ($tranchstart>$tranchfinish);
		next if (-e "$outdir/analyses/gatktranches/gatk.all.multiout.".$_.".vcf");
		$bsubcom = $bsubbase . $waitcom ." -o $logdir/gatktranches/tranch".$_.".out";
		$bsubcom = $bsubcom . " -J ".$runtag."_".$_."_gatktranch";
		$gatkcom = $gatkbase . " -L \\\"".$contigname.":".$tranchstart."-".$tranchfinish."\\\"";
		#-o gatk.all.multiout1.vcf"
		$gatkcom = $gatkcom . " -o $outdir/analyses/gatktranches/gatk.all.multiout.".$_.".vcf";
		runcommand "$bsubcom \"$gatkcom\""; 
	}
	
	
	
	
	
	
#now to make the next several in one shot
	
	
#perl ../marginalia/marginalia.pl -m "../marginalia/" -o staphout2 -q short --preprocess
#perl multitrancher.pl -o staphout2 -q long

#(in a bsub -- takes 20-30 seconds)
#perl multivcfmerger.pl -s staphout2/analyses/gatktranches -o staphout2/analyses
$perlcom = "perl $marginalia/multivcfmerger.pl -s $outdir/analyses/gatktranches -o $outdir/analyses;";
#perl ../marginalia/extractmutationsfromgbk.pl -i staphout2/analyses/gatk.all.multiout.merged.vcf -o staphout2/analyses/mutations.txt -r /groups/kishony/Reference_Genomes/SaureusNCTC8325
$perlcom = $perlcom." module load dev/perl/5.18.1;";
$perlcom = $perlcom." perl $marginalia/extractmutationsfromgbk.pl -i $outdir/analyses/gatk.all.multiout.merged.vcf -o $outdir/mutations.txt -r $refsource$refgenome{$protocol[$sampnum]};";
#perl multivcfparser.pl -s staphout2/analyses/gatk.all.multiout.merged.vcf -o staphout2/analyses/
$perlcom = $perlcom." perl $marginalia/multivcfparser.pl -s $outdir/analyses/gatk.all.multiout.merged.vcf -o $outdir/analyses;";
#-> copy the readfile.m to the output directory
$perlcom = $perlcom."cp $marginalia/readfiles.m $outdir/analyses/readfiles.m;";

$logbase = "$outdir/runlogs/vcfmergeandprocess";
$waitcom = '-w "done('.$runtag.'_*_gatktranch)"';
$bsubcom = "bsub -o $logbase.out -e $logbase.err -q short -W 1:0 ".$waitcom." -J ".$runtag."_makecsvs ";
runcommand "$bsubcom \"$perlcom\"";	
	
#need to run the matlab command to input everything
	
	
open NAMELIST, ">$outdir/analyses/namelist.txt";
	foreach $sampnum (0..$maxsamp) {
		if ($memberhash{$sampleoutname[$sampnum]}[0]==$sampnum)  {
		print NAMELIST $sampleoutname[$sampnum]."\n";
		}
}
close NAMELIST;


#dispatch readstomatlab.m making it wait on all previous work
#  matlab -r "cd marginalia;readfiles('/groups/kishony/baym/symcipmani2/analyses/');quit;"
		$logname = $logbase."matlabread";
	$bsubcommand = "bsub -o $logname.out -e $logname.err -q short -W 0:30 -w \"done(".$runtag."*gatktranch)\" -J ".$runtag."_matlabin";
	$matlabcommand = 'matlab -r "cd '.$outdir.'/analyses;readfiles;quit;"';
	 runcommand $bsubcommand.' '.$matlabcommand;
die;
	
	