#! /usr/bin/perl

use lib qw( ../../ );
use Test::Simple tests => 3;
use Feature::File::FeatureVector;

my $a = "../../../hits/trypsin_ser_og.pos.ff";
my $b = "../../../hits/trypsin_ser_og.neg.ff";
my $c = "./data/trypsin_ser_og.pos.err.ff";

my $file = new Feature::File::FeatureVector( $a );

ok( ! $file->diff( $a ) );
ok(   $file->diff( $b ) );
ok(   $file->diff( $c ) );
