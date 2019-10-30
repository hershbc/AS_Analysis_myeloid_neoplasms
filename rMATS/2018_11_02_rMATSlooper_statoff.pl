#!/usr/bin/perl

#This script will
# 1. select download bam files that match a given search term
# 2. Generate two lists with all bam files for rMATS
# 3. run rMATS on bam file lists
# 4. save output to S3 folder 

#Create list of files
my $condition = @ARGV[0];
my $control = @ARGV[1];
my $date = @ARGV[2];
my $label = $condition."_vs_".$control;
my $outputfolder = $date."_".$label;

#Get GTF files for rMATS
my $oops2 = system("aws s3 cp s3://allmergedgtf/allmergedgtfmerged.gtf .");
die "FAILED $!" if $oops2;


#Bringing in Bam files and 
#Creating input file 1
@condition;
open(IN, "<", $condition) or die "Can't open $condition\n";
while (my $line = <IN>){
    chomp $line;
    my $input = $line;
    $input =~ s/(.+)\s(.+)\s(.+)\s(.+)/$4/;
    push @condition, $input;

    my $oops3 = system("aws s3 cp s3://allbamfilesxs/$input .");
    die "FAILED $!" if $oops3;


}

my $condition_list = join( ',', @condition ); 
open(OUT, ">", "file1.txt") or die "Can't write to file1.txt\n";
print OUT $condition_list;

#Bringing in Bam files and 
#Preparing input file 2
@control;
open(IN2, "<", $control) or die "Can't open $control\n";
while (my $line2 = <IN2>){
    chomp $line2;
    my $input2 = $line2;
    $input2 =~ s/(.+)\s(.+)\s(.+)\s(.+)/$4/;
    push @control, $input2;

    my $oops4 = system("aws s3 cp s3://allbamfilesxs/$input2 .");
    die "FAILED $!" if $oops4;


}

my $control_list = join(',', @control);
open(OUT2, ">", "file2.txt") or die "Can't write to file2.txt\n";
print OUT2 $control_list;




#Running rMATS and saving output to S3
print "python rMATS.4.0.1/rMATS-turbo-Linux-UCS4/rmats.py --b1 file1.txt --b2 file2.txt --gtf allmergedgtfmerged.gtf --od $outputfolder -t paired --readLength 100 --cstat 0.1 --libType fr-firststrand --nthread 48 --statoff\n";

my $oops6 = system("python rMATS.4.0.1/rMATS-turbo-Linux-UCS4/rmats.py --b1 file1.txt --b2 file2.txt --gtf allmergedgtfmerged.gtf --od $outputfolder -t paired --readLength 100 --cstat 0.1 --libType fr-firststrand --nthread 48 --statoff");
die "FAILED!" if $oops6;

my $oops7 = system("cp {log.txt,file1.txt,file2.txt} $outputfolder");
die "FAILED!" if $oops7;

my $oops8 = system("aws s3 sync ./$outputfolder/ s3://rmats/$outputfolder/
");
die "FAILED!" if $oops8;


#Removing condition bam files
open(IN3, "<", $condition) or die "Can't open $condition\n";
while (my $line = <IN3>){
    chomp $line;
    my $input = $line;
    $input =~ s/(.+)\s(.+)\s(.+)\s(.+)/$4/;

    my $oops9 = system("rm $input");
    die "FAILED $!" if $oops9;


}

#Removing control bam files
open(IN4, "<", $control) or die "Can't open $control\n";
while (my $line2 = <IN4>){
    chomp $line2;
    my $input2 = $line2;
    $input2 =~ s/(.+)\s(.+)\s(.+)\s(.+)/$4/;

    my $oops10 = system("rm $input2");
    die "FAILED $!" if $oops10;


}



