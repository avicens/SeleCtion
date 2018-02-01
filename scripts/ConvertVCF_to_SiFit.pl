use warnings;
use strict;

my $input_file = $ARGV[0];
my $output_file = $ARGV[1];

my @individuals;
my %ind_genotype;
my @snps;
my $genotype;
open(FILE,$input_file);
while(<FILE>){
  if(/^#CHROM/){
    chomp;
    my @cols = split(/\t/);
    for (my $i=9; $i<scalar(@cols); $i++) { push @individuals, $cols[$i]; }
  }
  unless(/^#/){
    chomp;
    my @cols = split(/\t/);
    my $chr = $cols[0];
    my $pos = $cols[1];
    my $ref = $cols[3];
    my $alt = $cols[4]; 
    my $headSCG = $chr.":".$pos;
    push @snps, $headSCG;
    for (my $i=9; $i<scalar(@cols); $i++) { 
      my $ind_index = $i-9; 
      $cols[$i] =~ /^([a-z0-9.]\/[a-z0-9.])\:(\d+,\d+|\.)/;
	if ($1 eq "0/0") { $genotype = "0";}
	elsif ($1 eq "0/1") { $genotype = "1";}
        elsif ($1 eq "1/0") { $genotype = "1";}
        elsif ($1 eq "1/1") { $genotype = "1";}
	else { $genotype = "3";}
      push @{$ind_genotype{$individuals[$ind_index]}}, $genotype;
        }      
    }
}
close FILE;

print "cell_id\t".join("\t",@snps)."\n";
foreach my $ind (@individuals){
  print $ind."\t".join("\t",@{$ind_genotype{$ind}})."\n";
}

