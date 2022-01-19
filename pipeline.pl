#!/usr/bin/perl 
use warnings;
use strict;
use Data::Dumper;
use Parallel::ForkManager;

my $pm = Parallel::ForkManager->new(5);

my @list = `ls temp/*L00?*/*.fastq.gz`;
chomp @list;

#print Dumper (@list);
#exit;
my %ref;
foreach my $members(@list){
   
    if ($members =~ /(UNB.+_L00\d)/){
      my $log = $1;

       $ref{$log}{$members}=0;

    }

 }




#print Dumper (%ref);
#exit;
foreach my $sample (keys %ref){
    

    #if ($sample eq "UNB-RNA-2-1_S1_L002"){
    my $fq1;
    my $fq2;
    my $merged1;
    my $merged2;
    if ($sample =~ /(UNB.+)_L001/ ){
my $fid = $1;
my $sample2 = $fid."_L002";
foreach my $file (keys %{$ref{$sample}}){
   if ($file =~ /R1/){
foreach my $file2 (keys %{$ref{$sample2}}){
   if ($file2 =~ /R1/){
$fq1 = $file;
my $fq1_1 = $file2;
$merged1 = $fid."_R1_001.fq.gz";
system ("cat $fq1 $fq1_1 > $merged1");
   }
}
   }else{
foreach my $file2 (keys %{$ref{$sample2}}){
   if ($file2 =~ /R2/){
$fq2 = $file;
my $fq2_1 = $file2;
$merged2 = $fid."_R2_001.fq.gz";
system ("cat $fq2 $fq2_1 > $merged2");
   }
}
#$fq2 = $file;
   }
}
    
my $rcorrector_command = "perl /home/acri/tools/Rcorrector-master/run_rcorrector.pl -t 12 -1 $merged1 -2 $merged2";
    #print "$rcorrector_command\n";
my $cor1 = $fid."_R1_001.cor.fq.gz";
my $cor2 = $fid."_R2_001.cor.fq.gz";

system ("$rcorrector_command") unless -e $cor1;
die " $cor1 doesnt exist" unless -e $cor1;
# exit;
    
my $filter_command = "python ~/tools/TranscriptomeAssemblyTools-master/FilterUncorrectabledPEfastq.py -1 $cor1 -2 $cor2 -s $sample";
my $unf1 = "unfixrm_".$fid."_R1_001.cor.fq";
my $unf2 = "unfixrm_".$fid."_R2_001.cor.fq";
print "$filter_command\n";
system ("$filter_command") unless -e $unf1;
die " $unf1 doesnt exist" unless -e $unf1;
# system ("rm $cor1");
# system ("rm $cor2");
 #   print $filter_command,"\n";
my $gzip_command ="gzip $unf1";
my $gzip_command2 ="gzip $unf2";
my @gzips =($gzip_command,$gzip_command2);
my $gzipunf1 = $unf1.".gz";
my $gzipunf2 = $unf2.".gz";

for (my $a=0 ; $a <= 1; $a++){
   my $pid = $pm->start and next;
   system ("$gzips[$a]") unless -e $gzipunf1;
# print "$gzips[$a]\n";
   $pm->finish
}
#die " $gzipunf1 doesnt exist" unless -e $gzipunf1;
    
my $trimgalore_command = "~/tools/TrimGalore-0.6.6/trim_galore --paired $gzipunf1 $gzipunf2";
my $trimmed1 =  "unfixrm_".$fid."_R1_001.cor_val_1.fq.gz";
my $trimmed2 =  "unfixrm_".$fid."_R2_001.cor_val_2.fq.gz";
print "$trimgalore_command\n";
sleep(60);
system ("$trimgalore_command") unless -e $trimmed1;
die " $trimmed1 doesnt exist" unless -e $trimmed1;
    #system ("rm $gzipunf1");
    #system ("rm $gzipunf2");

my $output  = "trinity_".$fid;
my $trinity_command = "~/tools/trinityrnaseq-v2.12.0/Trinity --seqType fq --output $output --max_memory 32G --CPU 40 --left $trimmed1 --right $trimmed2";
print $trinity_command,"\n";
system ("$trinity_command");
my $trinity_output = $output."/"."Trinity.fasta";
my $stats_output = $output."/"."Trinity.stats";
my $stats = "~/tools/trinityrnaseq-v2.12.0/util/TrinityStats.pl  $trinity_output  > $stats_output";
system ("$stats");
    #    print "$stats\n";
my $index_output = $output."/"."trinity";
my $index = "bowtie2-build --threads 4 $trinity_output $index_output";
system ("$index");
my $alignment_stats = $output."/"."align_stats.txt";
my $bowtie2_output = $output."/"."bowtie2.bam";
my $bowtie2 = "bowtie2 -p 10 -q --no-unal -k 20 -x $index_output -1 $trimmed1 -2 $trimmed2  2> $alignment_stats|  samtools view -@ 10 -Sb -o $bowtie2_output";
system ("$bowtie2") unless -e $bowtie2;
#print "$index\n$bowtie2\n";
my $output_busco = $output."/"."busco";
my $busco = "~/tools/busco/build/scripts-3.6/busco -c 30 -m tran -i $trinity_output -o busco2 -l viridiplantae_odb10";
system ("$busco");
system ("mv busco2 $output_busco");
    }
    #exit;
    #}
}
