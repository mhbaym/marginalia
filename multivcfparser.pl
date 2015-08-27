#!/usr/bin/perl
use Getopt::Long;
use List::Util qw(max min reduce);

my $source = "gatk.all.multiout.vcf";
my $outdir = ".";

GetOptions ('s|v|invcf|source=s' => \$source,
			'o|outdir=s' => \$outdir, #(the directory with sample name last)
			'debug' => \$debug);

open MULTIVCF, $source;

open OUTDP, ">$outdir/depth.csv";
open OUTAGREE, ">$outdir/agree.csv";
open OUTCALL, ">$outdir/calls.tsv";
open OUTQUAL, ">$outdir/quals.csv";
open OUTVQMQ, ">$outdir/VqMqMqrsBqrsFsQd.csv";
open OUTREFS, ">$outdir/locsandrefs.tsv";

	while ($vcfline = <MULTIVCF>) {
	 	 next if ($vcfline =~ m/$\##/);
	 	 if ($vcfline =~ m/$\#/) {
	 	 		 @linesplitarray = split(" ",$vcfline);
	 	 		 $numitems = @linesplitarray;
	 	 		 @namearray =@linesplitarray[9..$numitems-1];
	 	 		 $numsamps = @namearray;
	 	 		 #print $namearray[0]." ".$numsamps." @namearray\n";
				 open OUTNAMES, ">$outdir/names.txt";
				 print OUTNAMES "$_\n" foreach (@namearray);
				 close OUTNAMES;
	 	 		 next;
	 	 }
	 	 
		 @linesplitarray = split(" ",$vcfline);
		 #print $vcfline."\n";
		 $pos = $linesplitarray[1];
		 $ref = $linesplitarray[3];
		 $call = $linesplitarray[4];
		 print OUTREFS "$pos\t$ref\t$call\n";
		 @calls = split(",",$call);
		 $qual = $linesplitarray[5];
		 
		 $linesplitarray[7] =~ m/MQ=(\d+.?\d*);/;
		 $mq = $1;
		 $linesplitarray[7] =~ m/MQRankSum=(\d+.?\d*);/;
		 $mqrs = $1;
		 $linesplitarray[7] =~ m/BaseQRankSum=(\d+.?\d*);/;
		 $bqrs = $1;
		 $linesplitarray[7] =~ m/FS=(\d+.?\d*);/;
		 $fs = $1;
		 $linesplitarray[7] =~ m/QD=(\d+.?\d*);/;
		 $qd = $1;
		 
		 print OUTVQMQ "$qual,$mq,$mqrs,$bqrs,$fs,$qd\n";
		 
		 $numcalls = @calls;
		 $numentries = @linesplitarray;
		 
		 
		 @tags = split(":",$linesplitarray[8]); #GT AD DP GQ PL
		 #print "@tags\n";
		my $gtarray=[];
		my $adarray=[];
		my $dparray=[];
		my $gqarray=[];
		my $plarray=[];
		for $i (0..$numsamps-1) {
		 	@values = split(":",$linesplitarray[9+$i]);
		 	#print "@tags @values";
		 	@qualtags{@tags}=@values;
		 	
		 	#$gtarray[$i] = $qualtags{"GT"};
			#$gqarray[$i] = $qualtags{"GQ"};  #these are diploid-optimized, so we'll ignore
			$adarray[$i] = $qualtags{"AD"};	#agree vs disagree
			$dparray[$i] = $qualtags{"DP"}; #read depth (sum of AD)
		 	$plarray[$i] = $qualtags{"PL"}; #
		 	# "@keys";
		 	
		 	@locpl = split(",",$qualtags{"PL"});
		 	@locad = split(",",$qualtags{"AD"});
		 	
		 	#this code uses reduce to do argmin;
		 	#$plsiz=@locpl;
		 	#$foo = reduce {$locpl[$a] < $locpl[$b] ? $a : $b } (0..($plsiz-1));
		 	
		 	$thisdp = $qualtags{"DP"};
		 	$thiscall=".";
		 	$thisqual=0;
		 	$thisagree=0;
		 	
		 	if ($linesplitarray[9+$i] eq "./.") {
		 		$thiscall=".";
		 		$thisqual=0;
		 		$thisagree=0;
		 		$thisdp=0;
		 	} elsif ($numcalls==1) {
		 		#make a call based on the PL
		 		if ($locpl[0]<$locpl[2]) { #the ref call is more likely
		 			$thiscall=$ref;
		 			$thisqual=$locpl[2]-$locpl[0];
		 			$thisagree = $locad[0];
		 		} else {
		 			$thiscall=$calls[0];
		 			$thisqual=$locpl[0]-$locpl[2];
		 			$thisagree = $locad[1];
		 		}
		 	
		 	} elsif ($numcalls==2) {
				$minval = reduce {$locpl[$a] < $locpl[$b] ? $a : $b } (0,2,5);
				#print $minval."  ";
				if ($minval == 0) {
					$thiscall=$ref;
		 			$thisqual=min($locpl[2],$locpl[5])-$locpl[0];
		 			$thisagree = $locad[0];
				} elsif ($minval == 2) {
					$thiscall=$calls[0];
		 			$thisagree = $locad[1];
		 			$thisqual=min($locpl[0],$locpl[5])-$locpl[2];
				} elsif ($minval == 5) {
					$thiscall=$calls[1];
		 			$thisagree = $locad[2];
		 			$thisqual=min($locpl[0],$locpl[2])-$locpl[5];
				}	
		 	
		 	}
		 	$thisdp = 0 if $thisdp eq ".";
		 	$thisagree = 0 if ($thisagree eq "." || $thisagree eq "");
		 	
		 	print "PROBLEM: $thispl ".$linesplitarray[9+$i] ." $thisqual $i $pos $thiscall\n" if $thisagree eq "";
		 	
		 
		 print OUTDP "$thisdp,";
		 print OUTAGREE "$thisagree,";
		 print OUTCALL "$thiscall\t";
		 print OUTQUAL "$thisqual,";
		 }
	
	print OUTDP "\n";
	print OUTAGREE "\n";
	print OUTCALL  "\n";
	print OUTQUAL  "\n";
	}
