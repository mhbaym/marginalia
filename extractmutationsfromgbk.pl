# first, bring in the SeqIO module
use Bio::SeqIO;
use Bio::Seq;
use List::MoreUtils qw(any);
use Getopt::Long;
use Bio::Tools::CodonTable;
 
my $startdir = '.';
my $refsource = '.';#'/groups/kishony/Reference_Genomes/';
my $outdir = '.';
my $mutfile = "newtest5variants.txt";
my $outfile = 'mutations.txt';
my $promoterlength=150;

GetOptions ('s|startdir=s' => \$startdir,	#directory with reads split into folders
			'r|g|refsource=s' => \$refsource, #where is the reference genome?
			'i|mutfile=s' => \$mutfile,		#where is the mutations vcf
			'o|outfile=s' => \$outfile,		#to what directory should analyses be outputted?
			'runtag=s' => \$runtag,			#unique tag for these processes? (to manage bsub dependencies)
			'unixcsv' => \$unixcsv,			#are samples.csv and alignment_params.csv unix or windows?
			'postprocess' => \$postprocess,			#are samples.csv and alignment_params.csv unix or windows?
			'debug' => \$debug);			#debug mode
 
 
 
# Bring in the file and format, or die with a nice
# usage statement if one or both arguments are missing.


#my $usage  = "perl testparse.pl [file] [format]\n";
#my $file   = shift or die $usage;
#my $format = shift or die $usage;
 

open VCFFILE, $mutfile || die "can't open VCF file $mutfile\n";;

	while ($vcfline = <VCFFILE>) {
	 	 next if ($vcfline =~ m/$\#/);
	 	 
		 @linesplitarray = split(" ",$vcfline);
		 #print $vcfline."\n";
		 #print $mutpos."\n";
		 $mutpos = $linesplitarray[1];
		 $refbases{$mutpos} = $linesplitarray[3];
		 
		 
		 push @mutlocs, $mutpos;
	
	}
	
close VCFFILE;

#print "<$refsource/genome.fasta";
open GENOME, "<$refsource/genome.fasta" || warn "can't open gemome";
	$genometopline = <GENOME>;
	$genometopline =~ m/\|.+\|.+\|(.+)\|/;
	$gbkfile = $refsource."/".$1.".gb";
close GENOME;
 
 
print "Reading genbank file $gbkfile  ....\n";
my $inseq = Bio::SeqIO->new(
                            -file   => "<$gbkfile",
                            -format => 'GenBank',
                            );
$myCodonTable   = Bio::Tools::CodonTable->new();
$myCodonTable->id(11);
 
 print "Parsing genbank...\n";
 
my $seq = $inseq->next_seq;


for my $feat_object ($seq->get_SeqFeatures) {         
	my $primtag =  $feat_object->primary_tag; 
    #print "primary tag: ", $primtag , "\n";  
  
	if ($feat_object->primary_tag eq "CDS") {
		$startpos = $feat_object->start;
		$endpos = $feat_object->end;
		$mutdata{$_}{"strand"}=$feat_object->strand;
	 #if (any { ($startpos < $_ && $endpos>$_) } @mutlocs ) {
	   foreach (@mutlocs) {  
		
	  	if ($startpos < $_ && $endpos>$_) {
	  	  $mutdata{$_}{"strand"}=$feat_object->strand;
	  	  $mutdata{$_}{"start"}=$feat_object->start;
	  	  $mutdata{$_}{"end"}=$feat_object->end;
	  	  $mutdata{$_}{"promotor"}=0;
		  #print $feat_object->start,"...";
		  #print $feat_object->end," ";
		  #print (($feat_object->strand eq 1)?"F":"R")."\n";
		  for my $tag ($feat_object->get_all_tags) {             
			#print "  tag: ", $tag, "\n";
			for my $value ($feat_object->get_tag_values($tag)) {                
				#print "    value: ", $value, "\n";
				$mutdata{$_}{$tag}=$value;     
			}  
		  }
		  
		  my $sequence_string = $feat_object->seq()->seq;
      	  $mutdata{$_}{"nt_sequence"}=$sequence_string;
      	  
		  $mutpos=(($feat_object->strand == 1)? ( $_ -$startpos):($endpos-$_)); #0-indexed position		  
      	  $mutdata{$_}{"nt_pos"}=$mutpos+1;
		  
		  $mutcodonstart = $mutpos-($mutpos % 3); 
      	  $mutdata{$_}{"aa_pos"}=int($mutcodonstart/3)+1;
		  #print $mutcodonstart."x y".length($sequence_string)."z a";
      	  $mutcodon = substr($sequence_string,$mutcodonstart,3);
      	  $mutdata{$_}{"codon"}=$mutcodon;
      	  #print "$mutcodon ";
      	  $altcodon = $mutcodon;
      	  $altcodings = "";
		  foreach $base ("A","T","C","G") {
		  	#print "($altcodon) ".($mutpos % 3)." ";
        	$z = substr $altcodon,($mutpos % 3), 1,$base;
      		$altcodings = $altcodings.($myCodonTable->translate($altcodon));
      		$altcodings =~ s/(.)(.)(.)(.)/$2$1$4$3/ if ($feat_object->strand == -1);
      		}
		  $mutdata{$_}{"AAs"}=$altcodings;
		
      #
    	} else { #INTERGENIC
    	
    		if (  ( (($startpos-$promoterlength) < $_)   &&   (($_ < $startpos)   &&   ($feat_object->strand==1)  ) )  ||
    		    ( (($endpos+$promoterlength) > $_)   &&  ( ($_ > $endpos)   &&   ($feat_object->strand==-1)  ))   ){

	  	 		$mutdata{$_}{"strand"}=$feat_object->strand;
	  	  		$mutdata{$_}{"start"}=$feat_object->start;
	  	  		$mutdata{$_}{"end"}=$feat_object->end;
	  	  		$mutdata{$_}{"promotor"}=1;
    			for my $tag ($feat_object->get_all_tags) {             
			#print "  tag: ", $tag, "\n";
				for my $value ($feat_object->get_tag_values($tag)) {                
				#print "    value: ", $value, "\n";
				$mutdata{$_}{$tag}=$value;     
				}  
		  		}
    		}
    		
    	} 
      }
           
	}  
   
   
      
}

open OUTFILE, ">$outfile";

foreach (@mutlocs) {
	print OUTFILE $_."\t";
	if ($mutdata{$_}{"gene"}) {
		print OUTFILE $mutdata{$_}{"gene"};
	} else {
		print OUTFILE "yyyY";
	}
	print OUTFILE "\t"; 
	print OUTFILE $mutdata{$_}{"start"}."\t";
	print OUTFILE $mutdata{$_}{"end"}."\t";
	print OUTFILE $mutdata{$_}{"strand"}."\t";
	print OUTFILE $mutdata{$_}{"promotor"}."\t";
	print OUTFILE $refbases{$_}."\t";
	print OUTFILE $mutdata{$_}{"product"}."\t";
	print OUTFILE $mutdata{$_}{"protein_id"}."\t";
	print OUTFILE $mutdata{$_}{"nt_pos"}."\t";
	print OUTFILE $mutdata{$_}{"codon"}."\t";
	print OUTFILE $mutdata{$_}{"locus_tag"}."\t";
	print OUTFILE $mutdata{$_}{"aa_pos"}."\t";
	print OUTFILE $mutdata{$_}{"AAs"}."\t";
	print OUTFILE $mutdata{$_}{"note"};
	#print OUTFILE " ".$mutdata{$_}{"function"}."\n";
	#print OUTFILE " ".$mutdata{$_}{"note"}."\n";
	#print OUTFILE " ".$mutdata{$_}{"translation"}."\n";
	print OUTFILE "\n";
	}