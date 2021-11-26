#! /usr/bin/perl

use lib qw( ../../ );
use Test::Simple tests => 2;
use Feature::File::Score;

my $a = "../../../hits/1bqy_ser_og.hits";
my $b = "../../../hits/1ufo_ser_og.hits";

my $file = new Feature::File::Score( $a );

ok( ! $file->diff( $a ) );
ok(   $file->diff( $b ) );
