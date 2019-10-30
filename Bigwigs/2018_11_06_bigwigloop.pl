#!/usr/bin/perl

my $files = shift;

open(IN, "<", $files) or die "Can't open $files\n";
while(my $line = <IN>){
    chomp $line;
    my $input = $line;
    $input =~ s/(.+)\s(.+)\s(.+)\s(.+)/$4/;

    my $output = $input;
    $output =~ s/(.+)_xs.bam/$1/;
    
    my $strand2 = $output.".Signal.Unique.str2.out.bg";
    my $strand1 = $output.".Signal.Unique.str1.out.bg";
    
    my $strand1bw = $strand1;
    $strand1bw =~ s/(.+).bg/$1.bw/;
    
    my $strand2bw =$strand2;
    $strand2bw =~ s/(.+).bg/$1.bw/;

print"$strand2, $strand1, $strand1bw, $strand2bw\n";    
    
    my $oops1 = system("aws s3 cp s3://allbamfilesxs/$input .");
    die "FAILED $!" if $oops1;
    
#create bedgraph
    print "/home/ec2-user/STAR-2.6.0a/source/STAR --runMode inputAlignmentsFromBAM --inputBAMfile $input --outWigType bedGraph --outWigStrand Stranded --outFileNamePrefix $output --outWigReferencesPrefix chr\n";
    
    my $oops2 = system("/home/ec2-user/STAR-2.6.0a/source/STAR --runMode inputAlignmentsFromBAM --inputBAMfile $input --outWigType bedGraph --outWigStrand Stranded --outFileNamePrefix $output --outWigReferencesPrefix chr");
    die "FAILED $!" if $oops2;
    
#sort bedgraph
    
        
    print "sort -k1,1 -k2,2n $strand2 > $strand2\n";
    print "sort -k1,1 -k2,2n $strand1 > $strand1\n";
    
    my $oops3 = system("sort -k1,1 -k2,2n $strand2 > $strand2");
    die "FAILED $!" if $oops3;

    my $oops4 = system("sort -k1,1 -k2,2n $strand1 > $strand1");
    die "FAILED $!" if $oops4;

#bedgraphtobigwig
    
    
    print "./bedGraphToBigWig $strand1 hg19.chrom.sizes $strand1bw\n";
    print "./bedGraphToBigWig $strand2 hg19.chrom.sizes $strand2bw\n";

    my $oops5 = system("./bedGraphToBigWig $strand1 hg19.chrom.sizes $strand1bw");
    die "FAILED $!" if $oops5;

    my $oops6 = system("./bedGraphToBigWig $strand2 hg19.chrom.sizes $strand2bw");
    die "FAILED $!" if $oops6;

#Move to s3

    my $oops7 = system("aws s3 cp $strand1bw s3://allbigwigs");
    die "FAILED $!" if $oops7;

    my $oops8 = system("aws s3 cp $strand2bw s3://allbigwigs");
    die "FAILED $!" if $oops8;

    my $oops9 = system("rm *.bg");
    die "FAILED $!" if $oops9;

    my $oops10 = system("rm *.bw");
    die "FAILED $!" if $oops10;

    my $oops11 = system("rm *.bam");
    die "FAILED $!" if $oops11;


}
