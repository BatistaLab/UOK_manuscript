#!/usr/bin/perl 
# get the longest transcripts from refGene.txt
# input file should be sorted by gene name
# 2015-09-01 
# xujin937@gmail.com

use strict;
use warnings;


if(@ARGV<1)
{
	print "Usage perl $0 <refGene.txt.sortedByName>\n";
	exit;
}


open IN, "$ARGV[0]" or die "can not open $ARGV[0]\n";


my $curr=<IN>;
while(<IN>)
{
	my @curr_a=split(/\t/,$curr);
	my @next_a=split(/\t/,$_);
	if($curr_a[12] eq $next_a[12])
	{
		# cmp.
		my $next_len=0;
		my @exon_start=split(/,/,$next_a[9]);
		my @exon_end=split(/,/,$next_a[10]);
		for(my $i=0;$i<=$#exon_start;$i++)
		{
			$next_len+=$exon_end[$i]-$exon_start[$i]+1;
		}
	
		my $curr_len=0;
		  		
	  	 @exon_start=split(/,/,$curr_a[9]);
                 @exon_end=split(/,/,$curr_a[10]);
                for(my $i=0;$i<=$#exon_start;$i++)
                {
                        $curr_len+=$exon_end[$i]-$exon_start[$i]+1;
                }

		if($next_len>=$curr_len)
		{
			$curr=$_;
		}


	}
	else
	{
		print $curr;
		$curr=$_;
	}
}

print $curr;
