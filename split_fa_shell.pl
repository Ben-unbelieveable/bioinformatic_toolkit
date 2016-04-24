#!/usr/bin/perl -w

=head1 Description

        Description split_fa_shell.pl  --  split *.fa file in to many subfile inother to improve the speed of the analyze.
		change the command line of 53 &71, so that you could get an shell file to run the next JOBs.

=head1 Version

  Author: Ben air, 614347533@qq.com
  Version: 1.0,  Date: 2016-4-8
  Note:

=head1 Usage
 
  perl split_fa_shell.pl -i <all.fa> -n <num> -o <out_file>
 

=head1 Exmple

  perl split_fa_shell.pl -i <all.fa> -n 1000 -o <out_file>

=cut
use Getopt::Long;

die `pod2text $0` if (@ARGV == 0 || $ARGV[0] =~/help/);
use Getopt::Long;
use File::Path;
use File::Basename;
use Cwd 'abs_path';
my ($samplefile,$dir);

my ($pesoap,$out_file);
GetOptions(
  "i|in-file=s"  =>\$pesoap,
  "n=i" =>\$num,
  "o|out-file=s"    =>\$outdir
  );
$fadir="$outdir/Split";
mkpath $fadir unless ( -d $fadir);
$shdir="$outdir/SH";
mkpath $shdir unless ( -d $shdir); 
$blastdir="$outdir/Blast/";
mkpath $blastdir unless ( -d $blastdir); 
open IN, $pesoap or die "ERROR: open $pesoap: $!";
$ma=1;
open $out, ">$fadir/split.$ma.fa" or die "ERROR: open $!";
open $shout, ">$shdir/blast.$ma.sh" or die "ERROR: open $!";
print $shout "#PBS -N blast.$ma\n";
print $shout "#PBS -l nodes=1:ppn=6\n";
print $shout "#PBS -q cu\n";
print $shout "#PBS -S /bin/bash\n";
print $shout "/lustre/Work/software/alignment/ncbi-blast-2.2.29+/bin/blastp -db /lustre/Work/software/alignment/ncbi-blast-2.2.29+/db/nr -query $fadir/split.$ma.fa -out $blastdir/blast.$ma.nr -num_threads 6 -outfmt 5 -evalue 1e-5 -max_target_seqs 10 \n";
open OUTM, ">$fadir/list" or die "ERROR: open $!";
while(defined($line=<IN>)) {
	chomp $line;
	if($line=~/^>/) {
		$i++;
		if($i>$num) {
			$ma++;
			$i=1;
			close $out;
			close $shout;
			open $out, ">$fadir/split.$ma.fa" or die "ERROR: open $!";
			print OUTM "$fadir/split.$ma.fa\n";
			open $shout, ">$shdir/blast.$ma.sh" or die "ERROR: open $!";
print $shout "#PBS -N blast.$ma\n";
print $shout "#PBS -l nodes=1:ppn=6\n";
print $shout "#PBS -q cu\n";
print $shout "#PBS -S /bin/bash\n";
print $shout "/lustre/Work/software/alignment/ncbi-blast-2.2.29+/bin/blastp -db /lustre/Work/software/alignment/ncbi-blast-2.2.29+/db/nr -query $fadir/split.$ma.fa -out $blastdir/blast.$ma.nr -num_threads 6 -outfmt 5 -evalue 1e-5 -max_target_seqs 10 \n";
		}
		print $out $line,"\n";
	}else{
		print $out $line,"\n";
	}
}
close IN;


		
	
