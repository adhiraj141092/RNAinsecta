#!/usr/bin/perl -w -I/home/stanley/ViennaRNA-1.5/Perl/blib/arch -I/home/stanley/ViennaRNA-1.5/Perl/blib/lib
############################################################################
# AUTHOR:  	Stanley NG Kwang Loong, stanley@bii.a-star.edu.sg
# DATE:		31/07/2005
# DESCRIPTION: Generates stats from fasta
############################################################################
use warnings;
use strict;
use Getopt::Long;
use RNA;

############################################################################
# Global Parameters and initialization if any.
############################################################################

my $inFile="&STDIN";
my $outFile="&STDOUT";

#Define the monomers and dimers
my %gl_monomers = ('A' => 0,'C' => 0, 'G' => 0, 'U' => 0);
my %gl_dimers = ('AA' => 0, 'AC' => 0, 'AG' => 0, 'AU' => 0,
                 'CA' => 0, 'CC' => 0, 'CG' => 0, 'CU' => 0,
                 'GA' => 0, 'GC' => 0, 'GG' => 0, 'GU' => 0,
                 'UA' => 0, 'UC' => 0, 'UG' => 0, 'UU' => 0);
my $numseqs = 0;
############################################################################
# File IO
# Parse the command line.
############################################################################
Getopt::Long::Configure ('bundling');
GetOptions (
	'i|input_file=s' => \$inFile, 
	'o|output_file=s' => \$outFile
);

if(scalar(@ARGV) == 1 || !defined($inFile) || !defined($outFile)) { 
	die ("USAGE: $0 -i <input file> -o <output file>\n");
}

open (INFILE, "<$inFile") or die( "Cannot open input file $inFile: $!" );
open (OUTFILE, ">$outFile") or die ("Cannot open output file $outFile: $!");

# ID Len A C G U G+C A+U AA AC AG AU CA CC CG CU GA GC GG GU UA UC UG UU %A %C %G %U %G+C %A+U %AA %AC %AG %AU %CA %CC %CG %CU %GA %GC %GG %GU %UA %UC %UG %UU bp %bp mfe Nmfe Q D Subopt_size  
print (OUTFILE "ID\tLen\t");
print (OUTFILE map { "$_\t" } (sort keys(%gl_monomers)));
print (OUTFILE "G\+C\tA\+U\t");
print (OUTFILE map { "$_\t" } (sort keys(%gl_dimers)));
print (OUTFILE map { "\%$_\t" } (sort keys(%gl_monomers)));
print (OUTFILE "\%G\+C\t\%A\+U\t");
print (OUTFILE map { "\%$_\t" } (sort keys(%gl_dimers)));
print (OUTFILE "pb\tNpb\t");
print (OUTFILE "mfe\tNmfe\t");
print (OUTFILE "Q\tNQ\t");
print (OUTFILE "D\tND\t");

# Read line by line.
while (my $line = uc(<INFILE>)) {

	chomp($line);
	$line =~ s/T/U/g;
	
	# Fasta First Line
    if ($line =~ m/^>/) { }
    
    # Fasta Second Line i.e. RNA sequence
    elsif ($line =~ m/^[AaCcUuGg]/) {	    

		#Absolute Values
		my %aw_monomers = %gl_monomers;
		my %aw_dimers = %gl_dimers;		

		$numseqs++;

		#remove white space etc
		$line =~ s/[^AaCcUuGg]//g;

		my $seqLen = length($line);
		print(OUTFILE "\n$numseqs\t$seqLen\t");
		
		#compute monomer and dimer distribution
		for my $i (0..$seqLen-1) {		
			my $monomer = substr($line, $i, 1);
			$aw_monomers{$monomer}++ if defined $aw_monomers{$monomer};

			my $dimer = substr($line, $i, 2);
			$aw_dimers{$dimer}++ if defined $aw_dimers{$dimer};			
		}	
		
		#Print Absolute Values
		foreach my $monomer (sort (keys(%aw_monomers))){
			print(OUTFILE "$aw_monomers{$monomer}\t");	
		}

		my $GC = $aw_monomers{'G'} + $aw_monomers{'C'};
		my $AU = $aw_monomers{'A'} + $aw_monomers{'U'};		
		print(OUTFILE "$GC\t$AU\t");
		
		foreach my $dimer (sort (keys(%aw_dimers))){
			print(OUTFILE "$aw_dimers{$dimer}\t");	
		}
		
		#Print Percentage Values
		foreach my $monomer (sort (keys(%aw_monomers))){
			printf(OUTFILE "%.2f\t", $aw_monomers{$monomer}/$seqLen*100);
		}

		printf(OUTFILE "%.2f\t%.2f\t", $GC/$seqLen*100, $AU/$seqLen*100);
		

		foreach my $dimer (sort (keys(%aw_dimers))){
			printf(OUTFILE "%.2f\t", $aw_dimers{$dimer}/($seqLen-1)*100);
		}
		
		my ($bp, $mfe, $Q, $D, $SS) = rnaAnalysis($line);
		printf(OUTFILE "%.2f\t%.4f\t", $bp, $bp/$seqLen);
		printf(OUTFILE "%.2f\t%.4f\t", $mfe, $mfe/$seqLen);
		printf(OUTFILE "%.2f\t%.4f\t", $Q, $Q/$seqLen);
		printf(OUTFILE "%.2f\t%.4f\t", $D, $D/$seqLen);

    }
	
    else { }
  
}#end of while loop

print (OUTFILE "\n");
close (INFILE) or die( "Cannot close input file $inFile: $!" );
close (OUTFILE) or die( "Cannot close output file $outFile: $!");
exit;

sub rnaAnalysis {
	my ($seq) = shift;
	my ($seqLen, $struct, $mfe) = (length($seq), RNA::fold($seq)); 
	my $bp = $struct =~ tr/(//; 
	my $Q = 0;
	my $D = 0;
		
	$RNA::pf_scale = exp((-1)*1.2*$mfe/(0.6163207755*$seqLen)) if ($seqLen > 1000); 

	# compute partition function and pair pobabilities matrix
	RNA::pf_fold($seq);   				
	# compute sum-of-entropy and bp-distance
	foreach my $j (1..$seqLen-1) {
		foreach my $k ($j+1..$seqLen) {
			my $p = RNA::get_pr($j, $k); # points to the computed pair probabilities
			if ($p > 0) {
				$Q += (-1)*$p*(log($p)/log(2));
				$D += $p*(1 - $p);
			}
		}
	}
	
	return ($bp, $mfe, $Q, $D, 0);
}