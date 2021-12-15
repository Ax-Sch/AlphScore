package Feature::GoldenRevision;

use Cwd qw( getcwd );
use Term::ReadKey;
use File::Path;
use strict;
use warnings;

=head1 NAME

  Feature::GoldenRevision

=head1 DESCRIPTION

A wrapper to the the original FEATURE v1.9 release, with fix-patches.

=cut

# ============================================================
sub new {
# ============================================================
	my ($class) = map { ref || $_ } shift;
	my $self = bless {}, $class;
	$self->init( @_ );
	return $self
}

# ============================================================
sub init {
# ============================================================
	my $self                = shift;
	my $release             = shift;
	$self->{ revision }     = 418;
	$self->{ release }      = $release;
	$self->{ feature }      = "golden-revision-$release";
	$self->{ svn_path }     = "https://simtk.org/svn/feature/tags/golden-revision-$release";
	$self->{ svn_tools }    = "https://simtk.org/svn/feature/trunk/tools";

	my $cwd = getcwd(); $cwd =~ s/\/$//;
	$self->{ input }        = "$cwd/data";
	$self->{ output }       = "$cwd/hits";
	$self->{ featurize }    = "$cwd/$self->{ feature }/featurize";
	$self->{ buildmodel }   = "$cwd/$self->{ feature }/buildmodel";
	$self->{ scoreit }      = "$cwd/$self->{ feature }/scoreit";
	$self->{ atomselector } = "$cwd/tools/bin/atomselector.py";

	$ENV{ FEATURE_DIR } = "$cwd/$self->{ feature }/";
	$ENV{ PATH }        = "$cwd/$self->{ feature }:$cwd/$self->{ feature }/tools/bin:$cwd/tools/bin:$ENV{ PATH }";
	$ENV{ PYTHONPATH }  = "$cwd/$self->{ feature }/tools/lib:$cwd/tools/lib";
	$ENV{ PDB_DIR }     = "$cwd/data/pdb";
	$ENV{ DSSP_DIR }    = "$cwd/data/dssp";

	mkpath( $self->{ output } ) unless -e $self->{ output };

	$self->svn_checkout_revision();
	$self->svn_checkout_tools();
	$self->build();
}

# ============================================================
sub build {
# ============================================================
	my $self    = shift;
	my $feature = $self->{ feature };

	# ===== SKIP IF THE GOLDEN REVISION HAS ALREADY BEEN BUILT
	return if( -e "$feature/featurize" && -e "$feature/buildmodel" && -e "$feature/scoreit" );

	my $command =
		"cd $feature && " .
		"aclocal        && " .
		"automake       && " .
		"autoconf" ;

	`$command`;
	print "Configure script created for $feature.\n";
	`cd $feature && ./configure`;
	`cd $feature && make 2>build.log`;
}

# ============================================================
sub cleanup {
# ============================================================
	my $self = shift;

	# ===== ADD NEW GOLDEN REVISIONS AS NEW APPROACHES ARE CREATED
	return unless $self->{ feature } =~ /^golden-revision-(?:1\.9\.1|metals|svm|rf)$/;
	`rm -rf $self->{ feature }`;
}

# ============================================================
sub get_user {
# ============================================================
	my $self = shift;
	ReadMode( 'normal' );
	print "Retrieving FEATURE $self->{ release } source code from SimTK\n";
	print "SimTK Username: ";
	my $username = <>;
	chomp $username;
	ReadMode( 'noecho' );
	print "SimTK Password: ";
	my $password = ReadLine( 0 );
	ReadMode( 'normal' );
	chomp $password;
	print "\n";
	$self->{ user } = { name => $username, pass => $password };
}

# ============================================================
sub set_test {
# ============================================================
	my $self = shift;
	my $test = shift;
	$self->{ test } = $test;
}

# ============================================================
sub svn_checkout_revision {
# ============================================================
	my $self    = shift;

	# ===== SKIP CHECKOUT IF THE CODE HAS ALREADY BEEN PREVIOUSLY CHECKED-OUT
	return if( -e $self->{ feature } );

	my $user    = $self->{ user } || $self->get_user();
	my $command =
		"svn checkout "               .            # Checkout command
		"--username $user->{ name } " .            # Username
		"--password $user->{ pass } " .            # Password
		"$self->{ svn_path } $self->{ feature } "; # URL

	local $_ = `$command`;
	if( -e $self->{ feature } ) {
		print "Checked out $self->{ feature }.\n";
	} else {
		die "Can't check out $self->{ feature }.\n";
	}

	print "Checkout complete.\n";
}

# ============================================================
sub svn_checkout_tools {
# ============================================================
	my $self    = shift;

	# ===== SKIP CHECKOUT IF THE CODE HAS ALREADY BEEN PREVIOUSLY CHECKED-OUT
	return if( -e 'tools' );

	my $user    = $self->{ user } || $self->get_user();
	my $command =
		"svn checkout "               . # Checkout command
		"--username $user->{ name } " . # Username
		"--password $user->{ pass } " . # Password
		"$self->{ svn_tools } tools ";  # URL

	local $_ = `$command`;
	if( -e 'tools' ) {
		print "Checked out 'tools'.\n";
		`ln -s tools $self->{ feature }/tools` unless -e "$self->{ feature }/tools";
	} else {
		die "Can't check out 'tools'.\n";
	}

	print "Checkout complete.\n";
}

sub get_input        { my $self = shift; die "Call set_test() prior to get_input()\n" unless $self->{ test };  return "$self->{ input }/$self->{ test }"; }
sub get_output       { my $self = shift; die "Call set_test() prior to get_output()\n" unless $self->{ test }; return "$self->{ output }/$self->{ test }"; }
sub get_featurize    { my $self = shift; return $self->{ featurize }; }
sub get_buildmodel   { my $self = shift; return $self->{ buildmodel }; }
sub get_scoreit      { my $self = shift; return $self->{ scoreit }; }
sub get_atomselector { my $self = shift; return $self->{ atomselector }; }

# ============================================================
sub buildmodel {
# ============================================================
	my $self       = shift;
	my $model      = shift;
	my $sites      = "$self->{ output }/$self->{ test }/$model/$model.pos.ff";
	my $nonsites   = "$self->{ output }/$self->{ test }/$model/$model.neg.ff";
	my $model_file = "$self->{ output }/$self->{ test }/$model/$model.model";

	print STDERR "    Building model $model... ";

	`$self->{ buildmodel } $sites $nonsites > $model_file 2>/dev/null`;
	if( -e $model_file ) {
		print STDERR "ok\n";
	} else {
		print STDERR "failed\n";
		exit( -1 );
	}
}

# ============================================================
sub featurize {
# ============================================================
	my $self  = shift;
	my $model = shift;
	my $file  = shift;
	my $type  = {
		'sites.ptf'    => { description => 'positive', posneg => 'pos' },
		'nonsites.ptf' => { description => 'negative', posneg => 'neg' }
	};

	my $description = $type->{ $file }{ description };
	my $posneg      = $type->{ $file }{ posneg };
	my $fv_file     = "$self->{ output }/$self->{ test }/$model/$model.$posneg.ff";
	my $point_file  = "$self->{ input }/$self->{ test }/$model/$file";
	print STDERR "    Featurizing $description training data set ($fv_file)... ";

	`$self->{ featurize } -P $point_file > $fv_file 2>/dev/null`;
	if( -e $fv_file ) {
		print STDERR "ok\n";
	} else {
		print STDERR "failed\n";
		exit( -1 );
	}
}

# ============================================================
sub train {
# ============================================================
	my $self  = shift;
	my $model = shift;
	$self->featurize( $model, "sites.ptf" );
	$self->featurize( $model, "nonsites.ptf" );
	$self->buildmodel( $model );

	print STDERR "    Packaging baseline model\n";
	my $positive_training_file = "$self->{ output }/$self->{ test }/$model/$model.pos.ff";
	my $model_file             = "$self->{ output }/$self->{ test }/$model/$model.model";

	`mv $self->{ output }/$self->{ test }/$model/$model.neg.ff $self->{ input }/$self->{ test }/$model`;
	`gzip -f $positive_training_file`;
	`gzip -f $model_file`;

}

# ============================================================
sub score {
# ============================================================
	my $self  = shift;
	my $model = shift;
	my $pdb   = shift;

	my $output       = $self->get_output();
	my $input        = $self->get_input();
	my $featurize    = $self->{ featurize };
	my $scoreit      = $self->{ scoreit };
	my $atomselector = $self->{ atomselector };
	my $point_file   = "$output/$model/$pdb.ptf";
	my $fv_file      = "$input/$model/score-$pdb.ff";
	my $model_file   = "$output/$model/$model.model";
	my $hits_file    = "$output/$model/$pdb.hits";

	print STDERR "      Featurizing $pdb... ";
	my ($pattern, $index, $residue, $atom) = split /\./, $model;
	`$atomselector -r $residue -a $atom $pdb > $point_file 2>/dev/null`;

	`$featurize -P $point_file > $fv_file 2>/dev/null`;
	if( -e $fv_file ) {
		print STDERR "ok\n";
		if( -e $point_file ) { unlink $point_file; }
	} else {
		print STDERR "failed\n";
	}
	print STDERR "      Scoring $pdb...     ";
	`$scoreit -a $model_file $fv_file > $hits_file`;
	if( -e "$hits_file" ) {
		print STDERR "ok\n";
		`gzip -f $hits_file`;
	} else {
		print STDERR "failed\n";
	}
} 

1;
