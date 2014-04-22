#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Carp;
use File::Basename;
use File::Copy;
use File::List;
use IO::Handle;
use IO::File;
use Bio::SeqIO;
use Bio::DB::Fasta;
use Text::CSV;
use Interpolation 'arg:@->$' => \&argument;
use Time::HiRes qw(gettimeofday tv_interval);

# Flush output after every write
select( (select(STDOUT), $| = 1 )[0] );

##
## This program
## 
## run with "--help" for usage information
##
## Robin Meyers

# Forward declarations
sub parse_command_line;
sub trim_adapter;
sub align_reads_to_reference;
sub perform_sam_manipulation;
sub run_homer;


# Global variables set by command line arguments
my $fastq;
my $bt2_index;
my $outdir;
my $threads = 4;
my $mem_max = "2G";
my $mapq_min = 10;
my $bt2_options = "--sensitive -p $threads -t";
my $adapter = "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC";


# More global variables
my $basename;
my $pre_dir;
my $align_dir;
my $samtools_dir;
my $homer_dir;
my $sam;
my $bam;
my $sort_filter_bam;

#
# Start of Program
#



parse_command_line;

mkdir $outdir;

trim_adapter;

align_reads_to_reference;

perform_sam_manipulation;

run_homer;

#
# End of program
#

sub run_homer {
  $homer_dir = "$outdir/homer";
  mkdir $homer_dir;

  my $maketagdirectory_out = "$homer_dir/maketagdirectory.out";  

  my $maketagdirectory_cmd = join(" ","makeTagDirectory",
                                      $homer_dir,
                                      "-genome",$homer_genome,
                                      "-checkGC -format sam",
                                      $bam,
                                      "> $maketagdirectory_out 2>&1");

  print "$maketagdirectory_cmd\n";
  my $maketagdirectory_status = system $maketagdirectory;
  croak "Error: problem executing 'makeTagDirectory'; check $maketagdirectory_out for details" unless $maketagdirectory_status == 0;


  my $makeucscfile_out = "$homer_dir/makeucscfile_out";

  my $makeucscfile_cmd = join(" ","makeUCSCfile",
                                  $homer_dir,
                                  "-o auto",
                                  "-strand separate",
                                  "> $makeucscfile_out 2>&1");

  print "$makeucscfile_cmd\n";
  my $makeucscfile_status = system $makeucscfile_cmd;
  croak "Error: problem executing 'makeUCSCfile'; check $makeucscfile_out for details" unless $makeucscfile_status == 0;


}

sub perform_sam_manipulation {

  print "Converting to bam format...\n";



  $bam = "$align_dir/$basename.bam";
  my $bam_out = "$align_dir/${basename}_bam.out";

  my $bam_cmd = join(" ","samtools view -bSh",
                         $sam,
                         "> $bam",
                         "2> $bam_out");

  print "$bam_cmd\n";
  my $bam_status = system($bam_cmd);

  croak "Error: problem executing 'samtools view'; see $bam_out for details;" unless $bam_status == 0;

  unlink $sam;

  print "Sorting and filtering on MAPQ...\n";

  $samtools_dir = "$outdir/samtools";
  mkdir $samtools_dir;

  $sort_filter_bam = "$samtools_dir/$basename.bam";
  my $sort_filter_out = "$samtools_dir/${basename}_sort_filter.out";
  my $tmpbam = "$samtools_dir/tmp.bam";

  my $sort_filter_cmd = join(" ","samtools view",
                                 "-b -q $mapq_min",
                                 $bam,
                                 "| samtools sort",
                                 "-@",$threads,
                                 "-m",$mem_max,
                                 "-o - $tmpbam > $sort_filter_bam",
                                 "2> $sort_filter_out");

  unlink $tmpbam;
  print "$sort_filter_cmd\n";
  my $sort_filter_status = system($sort_filter_cmd);
  croak "Error: problem executing 'samtools view' or 'samtools sort'; check $sort_filter_out for details" unless $sort_filter_status == 0;



}

sub align_reads_to_reference {

  print "Aligning reads to reference genome...\n";

  $align_dir = "$outdir/align";

  mkdir $align_dir;

  $sam = "$align_dir/$basename.sam";


  my $bt2_out = "$align_dir/${basename}_bt2.out";

  my $bt2_cmd = join(" ","bowtie2",
                         $bt2_options,
                         "-x",$bt2_index,
                         "-U",$fastq,
                         "-S",$sam,
                         "> $bt2_out 2>&1");

  print "$bt2_cmd\n";
  my $bt2_status = system $bt2_cmd;
  croak "Error: problem executing bowtie2-align; check $bt2_out for details" unless $bt2_status == 0;



}

sub trim_adapter {

  print "Trimming 3' Illumina adapter...\n";

  $pre_dir = "$outdir/preprocess";
  mkdir $pre_dir;

  my $trimmed_fastq = "$pre_dir/$basename.fastq";
  my $trim_out = "$pre_dir/${basename}_trim.out";

  my $trim_cmd = join(" ","homerTools trim",
                          "-3",$adapter,
                          "-min 15",
                          "-mis 1",
                          "-minMatchLength 4",
                          "-suffix trimmed",
                          "-lenSuffix lengths",
                          $fastq,
                          "> $trim_out 2>&1");
  print "$trim_cmd\n";
  my $trim_status = system $trim_cmd;
  croak "Error: problem executing 'homerTools trim'; check $trim_out for details" unless $trim_status == 0;

  $fastq =~ s/\.gz$//;

  rename "$fastq.trimmed", $trimmed_fastq;
  rename "$fastq.lengths", "$pre_dir/${basename}_lengths.txt";

  system "gzip $trimmed_fastq";

  $fastq = "$trimmed_fastq.gz";
}

sub parse_command_line {
  my $help;

  usage() if (scalar @ARGV==0);

  my $result = GetOptions (
        "fastq=s" => \$fastq,
        "bt2-index=s" => \$bt2_index,
        "bt2-opt=s" => \$bt2_options,
        "outdir=s" => \$outdir,
        "help" => \$help
      ) ;
  
  usage() if ($help);
 
  #Check options

  unless (defined $fastq && -r $fastq) {
    print "Error: must specify a valid fastq file\n";
    usage();
  }

  if ($fastq =~ /\.(fastq|fq)(\.gz)?/) {
    $basename = basename($`);
  } else {
    croak "Error: fastq must end with .fq or .fastq, or be gzipped with .fq.gz or .fastq.gz";
  }


  unless (defined $outdir) {
    print "Error: must specify output directory\n";
    usage();
  }

  if (-d $outdir) {
    my $outdir_search = File::List->new($outdir);
    croak "Error: output directory is not empty" if @{$outdir_search->find()} > 0;
  }

  unless (defined $bt2_index) {
    print "Error: must specify bowtie2 index\n";
    usage();
  }

  croak "Error: problem reading bowtie2 index" unless system("bowtie2-inspect -n $bt2_index > /dev/null") == 0;

  exit unless $result;
}

sub argument {
  my $var = shift;
  my $description = shift;
  my $default = shift;

  return sprintf("  \%\-16s - %s\n",
    (defined $default ? sprintf("%s (%s)",$var,$default) : sprintf("%s",$var)),
    $description);
}

sub usage()
{
print<<EOF;
Title, by Robin Meyers, ddmonthyyyy

This program .


Usage: $0 --fastq FILE --bt2-index IDX --outdir DIR [--opts]

Arguments (defaults in parentheses):

$arg{"--fastq","Input fastq file"}
$arg{"--bt2-index","Specify bowtie2 genome to use"}
$arg{"--bt2-opt","Specify bowtie2 parameters"}
$arg{"--outdir","Output directory"}
$arg{"--help","This helpful help screen."}


EOF

exit 1;
}