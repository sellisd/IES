#!/usr/bin/perl
use warnings;
use strict;
use Test::More;
use lib'.';
use functions;
my %floatingTests = (
    '01. Not floating' => {
	'ies'           => 'TACCTA',
	'macDownstream' => 'TACGCGCGCGTTA',
	'macUpstream'   => 'GGCGGGTA',
        'return'         =>  '-1'},
    '02. Floating downstream' => {
	'ies'           => 'TACCTA',
	'macDownstream' => 'TACCTAGGGG',
	'macUpstream'   => 'GGGGGGGGGGTA',
	'return'        => ['4']},
    '03. Floating upstream' => {
	'ies'           => 'TACCTA',
	'macDownstream' => 'TAGGGGGGGGG',
	'macUpstream'   => 'GGGGTACCTA',
	'return'        => ['-4']},
    '04. Repeated floating' => {
	'ies'           => 'TACCTA',
	'macDownstream' => 'TACCTACCTACCTA',
	'macUpstream'   => 'GGGGTA',
	'return'        =>  ['4','8','12']},
    '05. Possible floating with too short downstream flank' => {
	'ies'           => 'TACCTAGGGTA',
	'macDownstream' => 'TAC',
	'macUpstream'   => 'GGGGTA',
        'return'        => ['']},
    '06. Possible floating with too short upstream flank' => {
	'ies'           => 'TACCTAGGGCCCCTA',
	'macDownstream' => 'TAGGGGG',
	'macUpstream'   => 'CTA',
	'return'        => ['']},
    );
foreach my $floatTest (sort keys %floatingTests){
    print $floatTest,"\n---------------\n";
    my $iesS = $floatingTests{$floatTest}->{'ies'};
    my $macDownstreamS = $floatingTests{$floatTest}->{'macDownstream'};
    my $macUpstreamS = $floatingTests{$floatTest}->{'macUpstream'};
    my $returnValue = $floatingTests{$floatTest}->{'return'};
    print "IES: ", $iesS,"\n";
    print "Mac Downstream: ", $macDownstreamS, "\n";
    print "Mac Upstream: ", $macUpstreamS, "\n";
    eval{
	is_deeply(&isFloating($iesS, $macUpstreamS, $macDownstreamS), $returnValue);
    };
    if($@){
	print "Error:", $@;
    }
    print "\n";
}
done_testing();
