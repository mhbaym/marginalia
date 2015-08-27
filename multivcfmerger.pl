#!/usr/bin/perl
use Getopt::Long;

my $sourcedir = ".";
my $outdir = ".";

GetOptions ('s|sourcedir=s' => \$sourcedir,
			'o|outdir=s' => \$outdir, #(the directory with sample name last)
			'debug' => \$debug);

sub runcommand { #this lets me switch from testing to running easily!
if ($debug) {
	print $_[0]."\n";
} else {
	print $_[0]."\n";
	return system $_[0];
}
}

runcommand "mkdir -p $outdir";

open OUTVCF, ">$outdir/gatk.all.multiout.merged.vcf";

$filenum=0;
$maxpos=0;

while (open INVCF, $sourcedir."/gatk.all.multiout.".$filenum.".vcf") {
	print "File number $filenum opened...\n" if ($filenum % 100 == 0);
	while ($vcfline = <INVCF>) {
	 	 next if ($vcfline =~ m/$\#/ && $filenum>0);
	 	 @linesplitarray = split(" ",$vcfline);
		 $pos = $linesplitarray[1];
		 print OUTVCF $vcfline if ($pos > $maxpos || $vcfline =~ m/$\#/);
	}
	$filenum++;
}

$filenum--;
print "Last file was $filenum.\n"