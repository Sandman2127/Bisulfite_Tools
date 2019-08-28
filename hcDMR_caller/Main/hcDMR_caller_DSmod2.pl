#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;

#version 1. 20171025Yu

#################### open folder contains all input files#################
my $USAGE = "\nUSAGE: hcDMR_caller.pl 
                                   -ref *.gz Reference of multiple WTs
                                   -input *.gz 100bin file of interesting library
                                   -dif 0.1 for CHH, 0.2 for CHG, 0.4 for CHH
                                   -n minimum number of supporting libraries for each bin
                                   ";
				   
my $options = {};
GetOptions($options, "-ref=s", "-input=s", "-dif=s", "-n=s"); #, "-out=s" 
die $USAGE unless defined ($options->{ref});
die $USAGE unless defined ($options->{input});
die $USAGE unless defined ($options->{dif});
die $USAGE unless defined ($options->{n});

############################# Grobal Variables #############################
my $ref = $options->{ref};
my $input = $options->{input};
my $dif = $options->{dif};
my $n = $options->{n};

	
	open MU, 'gzip -dc '.$input.'|';
	open WT, 'gzip -dc '.$ref.'|';
	
	my @names = split /\./, $input; 
	open OUT, '>'.$names[0].'.'.$names[1].'.DMR';
	
	while ( !eof(MU) ) {
		
		#-------- WT ---------
		my(@array, $mean);	
		my $wt = <WT>;
		chomp $wt;
		my @wt = split /\s+/,$wt;	
		my $chr = $wt[0];
		my $pos = $wt[1];		
		#-------- Mutant -----
	 	my $line = <MU>;
		chomp $line;
		

			
			if($line =~ /^chr/) {
				print OUT "Chr\tbegin\tend\tDMR\tSampleMethylevel\tWTMethylevel\tmethdiff_total\tWTsum\tmethdiff_required\ttotal_libs_used\n";
			} elsif ($#wt == 2) {
			#	print OUT $chr, "\t", $pos, "\tNA\tNA\n"; 
			} else {
				
				my @lines = split /\s+/,$line;
				if ($lines[5]>=4 and $lines[2]+$lines[3]>0) {
					my $hyper = 0;
					my $hypo = 0;
				  	my $meth = $lines[2]/($lines[2]+$lines[3]);
				  	my $wtval = 0; 
					for my $i(@wt[2..$#wt]) {  
						if ($meth - $i >= $dif) {
							$hyper++;
							$wtval += $i;  
						} elsif ( $i - $meth >= $dif) {
							$hypo++;
							$wtval += $i;
						}
						else {
                                                        $wtval += $i;
                                                }
				  }
					  if ($hyper >= $n) {
							print OUT $chr,"\t",$pos-99,"\t", $pos, "\thyper\t", $meth,"\t",$wtval/$hyper,"\t",$meth-($wtval/$hyper),"\t",$wtval,"\t",$dif,"\t",$hyper,"\n"; 
					} elsif ($hypo >= $n) {
							print OUT $chr,"\t",$pos-99,"\t", $pos, "\thypo\t", $meth,"\t",$wtval/$hypo,"\t",($wtval/$hypo)-$meth,"\t",$wtval,"\t",$dif,"\t",$hypo,"\n"; 
					}
					  else {
                                                        print OUT $chr,"\t",$pos-99,"\t", $pos, "\tnot_significant\t", $meth,"\t",($wtval/54),"\t",($wtval/54)-$meth,"\t",$wtval,"\t",$dif,"\t","54","\n";
                                        }					

				}
			}	
	}	
		
	close MU;
	close WT;
	close OUT;	
	
