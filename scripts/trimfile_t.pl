
## create a hash from the exon file and the key is chrN_pos
## check chr file using its 1st and 3rd col as 'key'
## against the hash. if it exist in exon file, output the entire chr line
## to the chr file with "_t" appended
## f1 is chromosome file from pileup
## f2 is exon_contig chromosome file from tso_exon_contig split

#!/usr/bin/perl

 use strict;
 use warnings;

  my ($chrfile, $exofile) = @ARGV;

  open f1, $chrfile;
  open f2, $exofile;


  my $chrfile2 = $chrfile."_t"; # trimmed file name "_t" appended 
  open(my $fh, '>', $chrfile2) or die "Could not open file '$chrfile2' $!";

 my %exonhash = ();
 my (@A,@B);


while(<f2>){
  chomp;
  @A = split(/\t/,$_,7);
 my $exonkey = $A[2]."_".$A[5];
  $exonhash{$exonkey}=0;
  }

while(<f1>){
 chomp;
   @B = split(/\t/,$_,3);
   my $str = $B[0]."_".$B[1];
   matching($_, $fh, $str, \%exonhash);
}

sub matching{
   my ($line, $file, $key, $hash_ref) = @_;
     if(defined $hash_ref -> {$key}){
        print $file $line,"\n"; 
   }
}

close(f1);
close(f2);
close($fh);

