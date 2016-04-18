#!/usr/bin/perl -w
#=============================================================================
#     FileName: get_fasta_in_list.pl
#         Desc: 
#       Author: Ben air	
#        Email: 614347533@qq.com
#     HomePage:https://github.com/Ben-unbelieveable/bioinformatic_script
#      Version: 0.0.1
#   LastChange: 2016-4-8
#      History:
#=============================================================================

use strict;
use Getopt::Long;


=head1 Name

get_fasta_in_list.pl  --  get the fasta sequence in the list from all data sequnce

=head1 Description

	Description

=head1 Version

  Author: Ben air, 614347533@qq.com
  Version: 1.0,  Date: 2016-4-8
  Note:

=head1 Usage
 
  perl get_fasta_in_list.pl -fa <all.fa> -list <list.file> -out <out_file>
 

=head1 Exmple

  perl get_fasta_in_list.pl all.fa  list  out.fa

=cut

GetOptions(
  "fa|in-file=s"  =>\$pesoap,
  "list|in-file1=s"  =>\$sesoap,
  "out|out-file=s"    =>\$out_file
  );
  
die `pod2text $0` if (@ARGV == 0 || $Help);
open FA, $pesoap or die "ERROR: open $pesoap: $!";
open LIST, $sesoap or die "ERROR: open $sesoap: $!";
open OUT, ">$out_file" or die "ERROR: open >$out_file: $!";

while (defined($input=<INM>)) {
chomp ($input);
      if ($input=~/^>(\w+)/){
	      $seqname=$1;
       }else{
		$seq_name{$seqname}.=$input;
      }
}

while (defined($input=<INM>)) {
chomp ($input);
	$output.=$input."\n".$seq_name{$input}."\n";
}


print OUT$output;
