#! /usr/bin/perl

use lib qw( ../../ );
use Test::Simple tests => 2;
use Feature::File::Model;

my $a = "../../../hits/trypsin_ser_og.model";
my $b = "./data/trypsin_ser_og.err.model";

my $file = new Feature::File::Model( $a );

ok( ! $file->diff( $a ) );
ok(   $file->diff( $b ) );
