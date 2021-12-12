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
    print "Incorrect number of arguments\n Usage:  perl $0 protein_tax protein_list output_prefix \n";
}

my $lineage_path="/mnt1/database/GTDB/protein_db/gtdb_taxonomy.tsv.mo";
unless (@ARGV==3) {
    help;
    exit;
}

my $tsv = shift;
my $Protein_path = shift;
my $prefix = shift;

# $prefix=$prefix;

my %con;

open(Prot, $Protein_path) or die "cannot open protein list file... die";
while(<Prot>) {
    chomp();
    my @fs = split /\t/;
    my $conid = @fs[1];
    $con{$conid}{fap}++;
}

my %ids;
my $lastid = "none";

#parse lineage file
#/nvmessdnode3/opt/database/uniport/uniprot_trembl.name_taxid_lineage
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
#get seq id  in diamond result file by delete _ conid means raw id (not be gene predicted id )
    my(@cdata) = split(/_/, $qseqid);
    my $id = pop @cdata;
    my $conid = join("_", @cdata);

    # print "DMD: $conid\n";

    if ($qseqid ne $lastid) {
        # we have a top hit

        my $org = undef;
        if ($stitle =~ m/(.*?.\d+)_\d+/){
            $org = $1;
        } 
        # print "$org\n";
        my $tax = undef;
        if (exists $h{$org}){
            $tax = $h{$org};
        }
        my ($domain, $phylum, $class, $order, $family, $genus, $species) = split(/;/, $tax);
        # print "$org\t$tax\n";
        $con{$conid}{proteins}++;
        if ($frac > 0.7 ) {
            # print "$species\n";
            $con{$conid}{fulllen}++;        
        $con{$conid}{domain}{$domain}{sump}+= $pident;
        $con{$conid}{phylum}{$phylum}{sump}+= $pident;
        $con{$conid}{class}{$class}{sump}+= $pident;
        $con{$conid}{order}{$order}{sump}+= $pident;
        $con{$conid}{family}{$family}{sump}+= $pident;
        $con{$conid}{genus}{$genus}{sump}+= $pident;
        $con{$conid}{species}{$species}{sump}+= $pident;
        
		# # print "$domain\t$con{$conid}{domain}{$domain}{sumpn}\n";
		$con{$conid}{domain}{$domain}{sumpn}++;
        $con{$conid}{phylum}{$phylum}{sumpn}++;
        $con{$conid}{class}{$class}{sumpn}++;
        $con{$conid}{order}{$order}{sumpn}++;
        $con{$conid}{family}{$family}{sumpn}++;
        $con{$conid}{genus}{$genus}{sumpn}++;
        $con{$conid}{species}{$species}{sumpn}++;

        # print "$conid\t$species\t$con{$conid}{species}{$species}{sumpn}\n"
    }
    # Top hist
    $lastid = $qseqid;
    }
}
close IN;

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ for contigs/scaffolds
open(CON, ">$prefix".'.con.tax');
while(my($cn,$hr) = each %con) {
    # print "$cn\tdada\n"; 
    my $fl=$con{$cn}{fulllen};
    if ( $fl > 0 ) {
        print CON $cn."\t".$con{$cn}{fap}."\t".$con{$cn}{proteins}."\t".$con{$cn}{fulllen};
        my @sd = sort {$con{$cn}{domain}{$b}{sumpn} <=> $con{$cn}{domain}{$a}{sumpn} or $con{$cn}{domain}{$b}{sump} <=> $con{$cn}{domain}{$a}{sump}} keys %{$con{$cn}{domain}};
        my @sp = sort {$con{$cn}{phylum}{$b}{sumpn} <=> $con{$cn}{phylum}{$a}{sumpn} or $con{$cn}{phylum}{$b}{sump} <=> $con{$cn}{phylum}{$a}{sump}} keys %{$con{$cn}{phylum}};
        my @sc = sort {$con{$cn}{class}{$b}{sumpn} <=> $con{$cn}{class}{$a}{sumpn} or $con{$cn}{class}{$b}{sump} <=> $con{$cn}{class}{$a}{sump}} keys %{$con{$cn}{class}};
        my @so = sort {$con{$cn}{order}{$b}{sumpn} <=> $con{$cn}{order}{$a}{sumpn} or $con{$cn}{order}{$b}{sump} <=> $con{$cn}{order}{$a}{sump}} keys %{$con{$cn}{order}};
        my @sf = sort {$con{$cn}{family}{$b}{sumpn} <=> $con{$cn}{family}{$a}{sumpn} or $con{$cn}{family}{$b}{sump} <=> $con{$cn}{family}{$a}{sump}} keys %{$con{$cn}{family}};
        my @sg = sort {$con{$cn}{genus}{$b}{sumpn} <=> $con{$cn}{genus}{$a}{sumpn} or $con{$cn}{genus}{$b}{sump} <=> $con{$cn}{genus}{$a}{sump}} keys %{$con{$cn}{genus} };
        my @ss = sort {$con{$cn}{species}{$b}{sumpn} <=> $con{$cn}{species}{$a}{sumpn} or $con{$cn}{species}{$b}{sump} <=> $con{$cn}{species}{$a}{sump}} keys %{$con{$cn}{species}};
        
        # print "$cn\t$sp[0]\t$con{$cn}{phylum}{$sp[0]}{sump}\n";
        
        #get mean idendity of each taxonomy, and print result of them.
		if ($con{$cn}{domain}{$sd[0]}{sumpn} > 0) {
            my $pmean = sprintf("%0.2f", $con{$cn}{domain}{$sd[0]}{sump} / $con{$cn}{domain}{$sd[0]}{sumpn});
            print CON "\t", $sd[0], "\t", $con{$cn}{domain}{$sd[0]}{sumpn}, "\t", $pmean;
        } else {
            print CON "\t", $sd[0], "\t", $con{$cn}{domain}{$sd[0]}{sumpn}, "\t", "0";
        }
		
        if ($con{$cn}{phylum}{$sp[0]}{sumpn} > 0) {
            my $pmean = sprintf("%0.2f", $con{$cn}{phylum}{$sp[0]}{sump} / $con{$cn}{phylum}{$sp[0]}{sumpn});
            print CON "\t", $sp[0], "\t", $con{$cn}{phylum}{$sp[0]}{sumpn}, "\t", $pmean;
        } else {
            print CON "\t", $sp[0], "\t", $con{$cn}{phylum}{$sp[0]}{sumpn}, "\t", "0";
        }
        
        if ($con{$cn}{class}{$sc[0]}{sumpn} > 0) {
            my $cmean = sprintf("%0.2f", $con{$cn}{class}{$sc[0]}{sump} / $con{$cn}{class}{$sc[0]}{sumpn});
            print CON "\t", $sc[0], "\t", $con{$cn}{class}{$sc[0]}{sumpn}, "\t", $cmean;
        } else {
            print CON "\t", $sc[0], "\t", $con{$cn}{class}{$sc[0]}{sumpn}, "\t", "0";
        }
        if ($con{$cn}{order}{$so[0]}{sumpn} > 0) {
            my $omean = sprintf("%0.2f", $con{$cn}{order}{$so[0]}{sump} / $con{$cn}{order}{$so[0]}{sumpn});
            print CON "\t", $so[0], "\t", $con{$cn}{order}{$so[0]}{sumpn}, "\t", $omean;
        } else {
            print CON "\t", $so[0], "\t", $con{$cn}{order}{$so[0]}{sumpn}, "\t", "0";
        }
        
        if ($con{$cn}{family}{$sf[0]}{sumpn} > 0) {
            my $fmean = sprintf("%0.2f", $con{$cn}{family}{$sf[0]}{sump} / $con{$cn}{family}{$sf[0]}{sumpn});
            print CON "\t", $sf[0], "\t", $con{$cn}{family}{$sf[0]}{sumpn}, "\t", $fmean;
        } else {
            print CON "\t", $sf[0], "\t", $con{$cn}{family}{$sf[0]}{sumpn}, "\t", "0";
        }
        
        if ($con{$cn}{genus}{$sg[0]}{sumpn} > 0) {
            my $gmean = sprintf("%0.2f", $con{$cn}{genus}{$sg[0]}{sump} / $con{$cn}{genus}{$sg[0]}{sumpn});
            print CON "\t", $sg[0], "\t", $con{$cn}{genus}{$sg[0]}{sumpn}, "\t", $gmean;
        } else {
            print CON "\t", $sg[0], "\t", $con{$cn}{order}{$sg[0]}{sumpn}, "\t", "0";
        }
        
        if ($con{$cn}{species}{$ss[0]}{sumpn} > 0) {
            my $smean = sprintf("%0.2f", $con{$cn}{species}{$ss[0]}{sump} / $con{$cn}{species}{$ss[0]}{sumpn});
            print CON "\t", $ss[0], "\t", $con{$cn}{species}{$ss[0]}{sumpn}, "\t", $smean;
        } else {
            print CON "\t", $ss[0], "\t", $con{$cn}{species}{$ss[0]}{sumpn}, "\t", "0";
        }
        print CON "\n";
    }
}

close CON;
