#!/usr/bin/env perl
#jinh 20191111
#update for getting clear notes for each contig/scaffolds or bins
################################################################
#可能存在一些bug, 有待反馈
#输出文件名问题, 待后续搭建流程时统一修改
#Future version: 提供最优注释信息
################################################################
use strict;
# use Bio::SeqIO;
use File::Basename;

sub help{
    print "Incorrect number of arguments\n Usage:  perl $0 protein_tax protein_list scaf2bin_file output_perfix \n";
}

my $lineage_path="/mnt1/database/GTDB/protein_db/gtdb_taxonomy.tsv.mo";
unless (@ARGV==4) {
    help;
    exit;
}

my $tsv = shift;
my $Protein_path = shift;
my $Scaf2bin_path = shift;
my $prefix = shift;

# my $binid=$prefix;
$prefix=$prefix."_";
my %bin;
my $obin=$prefix."_bin";

my %scaf;

open(Scaf2bin, $Scaf2bin_path) or die "cannot open scaf2bin file... die";
while(<Scaf2bin>) {
    chomp();
    my @ss = split /\t/;
    $scaf{@ss[0]} = @ss[1];
}


open(Prot, $Protein_path) or die "cannot open protein list file... die";
while(<Prot>) {
    chomp();
    my @fs = split /\t/;
    my $conid = @fs[1];
    my $binid = $scaf{$conid};
    # print $conid;
    $bin{$binid}{fap}++;
}

my %ids;
my $lastid = "none";

#parse lineage file
my %lineage;
my %h;

open(Lin, $lineage_path) or die "cannot open lineage files... die";
print "Parsing lineage database...\n";
while(<Lin>) {
    chomp();
    my @s = split /\t/;
    $h{@s[0]}=@s[1];
}

close Linn;
print "Getting results\n";

open(IN, $tsv) || die "cannot open tsv\n";
while(<IN>) {
    chomp();
    my ($qseqid,$sseqid,$stitle,$pident,$qlen,$slen,$length,$mismatch,$gapopen,$qstart,$qend,$sstart,$send,$evalue,$bitscore) = split(/\t+/);
    my $frac = $qlen / $slen;
    my(@cdata) = split(/_/, $qseqid);
    my $id = pop @cdata;
    my $conid = join("_", @cdata);
    my $binid = $scaf{$conid};

    if ($qseqid ne $lastid) {
        # we have a top hit
        my $org = undef;
        if ($stitle =~ m/(.*?.\d+)_\d+/){
            $org = $1;
        }
        my $tax = undef;
        if (exists $h{$org}){
            $tax = $h{$org};
        }
        my ($domain, $phylum, $class, $order, $family, $genus, $species) = split(/;/, $tax);

        $bin{$binid}{proteins}++;
        if ($frac > 0.7) {
            $bin{$binid}{fulllen}++;

        $bin{$binid}{phylum}{$phylum}{sump}+= $pident;
        $bin{$binid}{class}{$class}{sump}+= $pident;
        $bin{$binid}{order}{$order}{sump}+= $pident;
        $bin{$binid}{family}{$family}{sump}+= $pident;
        $bin{$binid}{genus}{$genus}{sump}+= $pident;
        $bin{$binid}{species}{$species}{sump}+= $pident;
        
        $bin{$binid}{phylum}{$phylum}{sumpn}++;
        $bin{$binid}{class}{$class}{sumpn}++;
        $bin{$binid}{order}{$order}{sumpn}++;
        $bin{$binid}{family}{$family}{sumpn}++;
        $bin{$binid}{genus}{$genus}{sumpn}++;
        $bin{$binid}{species}{$species}{sumpn}++;
    }
    # Top hist
    $lastid = $qseqid;
    }
}
close IN;

# unless (-d $outdir) {
    # mkdir $outdir;
# }

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ for bin
open(BIN, ">$prefix".'bin');
while(my $bn = each(%bin)){
    # print $bn;
    print BIN $bn."\t".$bin{$bn}{fap}."\t".$bin{$bn}{proteins}."\t".$bin{$bn}{fulllen};
    
    my @sp = sort {$bin{$bn}{phylum}{$b}{sumpn} <=> $bin{$bn}{phylum}{$a}{sumpn} or $bin{$bn}{phylum}{$b}{sump} <=> $bin{$bn}{phylum}{$a}{sump}} keys %{$bin{$bn}{phylum}};
    my @sc = sort {$bin{$bn}{class}{$b}{sumpn} <=> $bin{$bn}{class}{$a}{sumpn} or $bin{$bn}{class}{$b}{sump} <=> $bin{$bn}{class}{$a}{sump}} keys %{$bin{$bn}{class}};
    my @so = sort {$bin{$bn}{order}{$b}{sumpn} <=> $bin{$bn}{order}{$a}{sumpn} or $bin{$bn}{order}{$b}{sump} <=> $bin{$bn}{order}{$a}{sump}} keys %{$bin{$bn}{order}};
    my @sf = sort {$bin{$bn}{family}{$b}{sumpn} <=> $bin{$bn}{family}{$a}{sumpn} or $bin{$bn}{family}{$b}{sump} <=> $bin{$bn}{family}{$a}{sump}} keys %{$bin{$bn}{family}};
    my @sg = sort {$bin{$bn}{genus}{$b}{sumpn} <=> $bin{$bn}{genus}{$a}{sumpn} or $bin{$bn}{genus}{$b}{sump} <=> $bin{$bn}{genus}{$a}{sump}} keys %{$bin{$bn}{genus}};
    my @ss = sort {$bin{$bn}{species}{$b}{sumpn} <=> $bin{$bn}{species}{$a}{sumpn} or $bin{$bn}{species}{$b}{sump} <=> $bin{$bn}{species}{$a}{sump}} keys %{$bin{$bn}{species}};


    #get mean idendity of each taxonomy, and print result of them.
    if ($bin{$bn}{phylum}{$sp[0]}{sumpn} > 0) {
        my $pmean = sprintf("%0.2f", $bin{$bn}{phylum}{$sp[0]}{sump} / $bin{$bn}{phylum}{$sp[0]}{sumpn});
        print BIN "\t", $sp[0], "\t", $bin{$bn}{phylum}{$sp[0]}{sumpn}, "\t", $pmean;
    } else {
        print BIN "\t", $sp[0], "\t", $bin{$bn}{phylum}{$sp[0]}{sumpn}, "\t", "0";
    }
    
    if ($bin{$bn}{class}{$sc[0]}{sumpn} > 0) {
        my $cmean = sprintf("%0.2f", $bin{$bn}{class}{$sc[0]}{sump} / $bin{$bn}{class}{$sc[0]}{sumpn});
        print BIN "\t", $sc[0], "\t", $bin{$bn}{class}{$sc[0]}{sumpn}, "\t", $cmean;
    } else {
        print BIN "\t", $sc[0], "\t", $bin{$bn}{class}{$sc[0]}{sumpn}, "\t", "0";
    }
    if ($bin{$bn}{order}{$so[0]}{sumpn} > 0) {
        my $omean = sprintf("%0.2f", $bin{$bn}{order}{$so[0]}{sump} / $bin{$bn}{order}{$so[0]}{sumpn});
        print BIN "\t", $so[0], "\t", $bin{$bn}{order}{$so[0]}{sumpn}, "\t", $omean;
    } else {
        print BIN "\t", $so[0], "\t", $bin{$bn}{order}{$so[0]}{sumpn}, "\t", "0";
    }
    
    if ($bin{$bn}{family}{$sf[0]}{sumpn} > 0) {
        my $fmean = sprintf("%0.2f", $bin{$bn}{family}{$sf[0]}{sump} / $bin{$bn}{family}{$sf[0]}{sumpn});
        print BIN "\t", $sf[0], "\t",$bin{$bn}{family}{$sf[0]}{sumpn}, "\t", $fmean;
    } else {
        print BIN "\t", $sf[0], "\t",$bin{$bn}{family}{$sf[0]}{sumpn}, "\t", "0";
    }
    
    if ($bin{$bn}{genus}{$sg[0]}{sumpn} > 0) {
        my $gmean = sprintf("%0.2f", $bin{$bn}{genus}{$sg[0]}{sump} / $bin{$bn}{genus}{$sg[0]}{sumpn});
        print BIN "\t", $sg[0], "\t", $bin{$bn}{genus}{$sg[0]}{sumpn}, "\t", $gmean;
    } else {
        print BIN "\t", $sg[0], "\t", $bin{$bn}{order}{$sg[0]}{sumpn}, "\t", "0";
    }
    
    if ($bin{$bn}{species}{$ss[0]}{sumpn} > 0) {
        my $smean = sprintf("%0.2f", $bin{$bn}{species}{$ss[0]}{sump} / $bin{$bn}{species}{$ss[0]}{sumpn});
        print BIN "\t", $ss[0], "\t", $bin{$bn}{species}{$ss[0]}{sumpn}, "\t", $smean;
    } else {
        print BIN "\t", $ss[0], "\t", $bin{$bn}{species}{$ss[0]}{sumpn}, "\t", "0";
    }
    print BIN "\n";
}
close BIN;

#exit;