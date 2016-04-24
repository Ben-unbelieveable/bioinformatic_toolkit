
#!/usr/bin/perl -w

=head1 Description

        Description get_fasta_in_list.pl  --  get the fasta sequence in the list from all data sequnce

=head1 Version

  Author: Ben air, 614347533@qq.com
  Version: 1.0,  Date: 2016-4-8
  Note:

=head1 Usage
 
  perl get_fasta_in_list.pl -fa <all.fa> -list <list.file> -out <out_file>
 

=head1 Exmple

  perl get_fasta_in_list.pl all.fa  list  out.fa

=cut
use Getopt::Long;

die `pod2text $0` if (@ARGV == 0 || $ARGV[0] =~/help/);
GetOptions(
  "fa|in-file=s"  =>\$pesoap,
  "list|in-file1=s"  =>\$sesoap,
  "out|out-file=s"    =>\$out_file
  );

$out_unmat=$out_file."_unmatch";
open FA, $pesoap or die "ERROR: open $pesoap: $!";
open LIST, $sesoap or die "ERROR: open $sesoap: $!";
open OUT, ">$out_file" or die "ERROR: open >$out_file: $!";
open OUTUNMATCH, ">$out_unmat" or die "ERROR: open >$out_unmat: $!";
while (defined($input=<FA>)) {
chomp ($input);
      if ($input=~/gene:(\w+)/){
              $seqname=$1;
       }else{
                $seq_name{$seqname}.=$input;
      }
}

while (defined($input=<LIST>)) {
chomp ($input);
        if(defined ($seq_name{$input})){
        $output.=">".$input."\n".$seq_name{$input}."\n";
        }else{
        $out_unmatch.=$input."\n";
        }
}
print OUT$output;
print OUTUNMATCH$out_unmatch;
