#!/usr/bin/perl

use strict;
use warnings;

#Program takes gtf files from S3 bucket
#Program runs stringtie using UCSC hg19 gtf file as guide
#Program saves gtf output to merged gtf folder in S3

my @existingfiles;

open(IN2, "<", "/home/ec2-user/existinggtffiles.txt") or die "Can't open existinggtffiles.txt $!";
while(my $line = <IN2>){
    chomp $line;
    $line =~ s/(.+)\s(.+)\s(.+)\s(.+)/$4/;
    push(@existingfiles, $line);
}


my $file = shift; 
my $file2 = $file."list";
my $output = $file."merged.gtf";

#Open in list of files
open(IN, "<", $file) or die "can't open $file $!\n";
open(OUT, ">", $file2) or die "can't write to $file2 $!\n";


#Create file for stringtie merge
while(my $line = <IN>){
    chomp $line;    
    my $input = $line.".gtf";    
    
    foreach my $gtf(@existingfiles){
	if($input eq $gtf){
	    
	    print "loading $input\n";
	    
	    my $oops = system("aws s3 cp s3://allgtf/$input .");
	    die "FAILED $!" if $oops;
	    
	    print OUT "$input\n"; 
	    
	}
    }
    
}    
print "running stringtie --merge -G genes.gtf -o $output -l $file -i $file2\n";
    
my $oops2 = system("stringtie --merge -G genes.gtf -o $output -l $file -i $file2");
die "FAILED $!" if $oops2;

print "copying file back to amazon\n";

my $oops3 = system("aws s3 cp $output s3://allmergedgtf/");
die "FAILED!" if $oops3;

