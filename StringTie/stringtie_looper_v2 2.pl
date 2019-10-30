#!/usr/bin/perl

use strict;
use warnings;

#Program takes files from S3 bucket with xs tags
#Program runs stringtie using UCSC hg19 gtf file
#Program saves gtf output to gtf folder in S3
#I will run 10 of these scripts at the same time to make stringtie run as fast as possible.

my $file = shift; 


#Open in list of files
open(IN, "<", $file) or die "can't open $file $!\n";

#Save identifier (everything after output but before bam)
while(my $line = <IN>){
    chomp $line;
    
    $line =~ m/(.+)\s(.+)\s(.+)\s(.+)_xs.bam/;

    my $filename = $4;
    
    my $bam = join'',$filename,"_xs.bam";
    my $gtf = $filename.".gtf";


    print "$bam from S3, it is labeled $filename and we will use it to create $gtf\n";
    
    my $oops = system("aws s3 cp s3://allbamfilesxs/$bam .");
    die "FAILED $!" if $oops;

    print "running stringtie command stringtie -l $filename -p 16 --rf -G genes.gtf -o $gtf  $bam\n";

    my $oops2 = system("stringtie -l $filename -p 16 --rf -G genes.gtf -o $gtf  $bam");
    die "FAILED $!" if $oops2;

    print "copying file back to amazon\n";

    my $oops3 = system("aws s3 cp $gtf s3://allgtf/");
    die "FAILED!" if $oops3;

    print "removing file\n";

    my $oops4 = system("rm $filename*");
    die "FAILED!" if $oops4;


}
