package Feature::Test;
# ============================================================
# by Marc Sosnick, San Francisco State University
# Sets up and runs current release candidate
# ============================================================

use Test::Harness;
use Cwd qw( getcwd );

# ============================================================
sub new {
# ============================================================
	my ($class) = map { ref || $_ } shift;
	my $self = bless {}, $class;
	$self->init( @_ );
	return $self;
}

# ============================================================
sub init {
# ============================================================
	my $self   = shift;
	my $source = shift || '.';

	$self->{ cwd }     = getcwd();
	$self->{ logfile } = "$self->{ cwd }/test.log";
	$self->{ apps }    = '';
	$self->{ data }    = '';
	$self->{ release } = '';
	$self->{ tools }   = '';

	$self->find_release( $source );
	$self->find_apps( $source );
	$self->find_data( $source );
	$self->find_tools( $source );

	$ENV{ PDB_DIR }     = 'data/pdb';
	$ENV{ DSSP_DIR }    = 'data/dssp';
	$ENV{ PYTHONPATH }  = $self->tools() . '/lib';
	$ENV{ FEATURE_DIR } = $self->data();
}

# ============================================================
sub find_apps {
# ============================================================
	my $self = shift;
	my $path = shift;

	my $file = "featurize";
	my @bin = split /\n/, `ls $path/$self->{ release }/$file $path/$self->{ release }/bin/$file 2>/dev/null`;
	return unless @bin;

	my $bin = shift @bin;
	$bin =~ s/\/$file//;
	$bin =~ s/^\.\///;
	
	$self->{ apps } = $bin;
}

# ============================================================
sub find_data {
# ============================================================
	my $self = shift;
	my $path = shift;

	my $file = "residue_templates.dat";
	my @data = split /\n/, `ls $path/$self->{ release }/$file $path/$self->{ release }/data/$file 2>/dev/null`;
	return unless @data;

	my $data = shift @data;
	$data =~ s/\/$file//;
	$data =~ s/^\.\///;
	
	$self->{ data } = $data;
}

# ============================================================
sub find_release {
# ============================================================
	my $self = shift;
	my $path = shift;

	my @releases = split /\n/, `ls $path/feature-*-src.tar.gz $path/feature-*-rc.tar.gz 2>/dev/null`;
	return unless @releases;

	my $release = shift @releases;
	$release =~ s/^\.\///;
	$self->{ tarball } = $release;

	$release =~ s/\-(?:src|rc)\.tar\.gz//;
	$self->{ release } = $release;
}

# ============================================================
sub find_tools {
# ============================================================
	my $self = shift;
	my $path = shift;

	my @tools = split /\n/, `ls $path/$self->{ release }/tools/bin/atomselector.py $path/tools-*-rc.tar.gz 2>/dev/null`;
	die "Feature::Test error: Can't find any FEATURE tools release candidates in '$path' $!" unless @tools;

	my $tools = shift @tools;
	$tools =~ s/^\.\///;
	$tools =~ s/\-(?:src|rc)\.tar\.gz//;
	$tools =~ s/\/bin\/atomselector.py//;
	$self->{ tools } = $tools;
}

sub apps          { my $self = shift; die "Feature::Test error: Can't find any FEATURE applications--have they been compiled?\n" unless $self->{ apps }; return "$self->{ cwd }/$self->{ apps }";    }
sub data          { my $self = shift; return "$self->{ cwd }/$self->{ data }";              }
sub tarball       { my $self = shift; return "$self->{ cwd }/$self->{ tarball }";           }
sub release       { my $self = shift; die "Feature::Test error: Can't find any FEATURE release candidates or source distributions.\n" unless $self->{ release }; return "$self->{ cwd }/$self->{ release }";           }
sub tools         { my $self = shift; return "$self->{ cwd }/$self->{ tools }";             }
sub tools_package { my $self = shift; return "$self->{ cwd }/$self->{ tools }-rc.tar.gz";   }

sub featurize     { my $self = shift; return "$self->{ cwd }/$self->{ apps }/featurize";    }
sub buildmodel    { my $self = shift; return "$self->{ cwd }/$self->{ apps }/buildmodel";   }
sub scoreit       { my $self = shift; return "$self->{ cwd }/$self->{ apps }/scoreit";      }

sub atomselector  { my $self = shift; return "$self->{ cwd }/$self->{ tools }/bin/atomselector.py"; }


# ============================================================
sub run {
# ============================================================
# - executes system command
# - logs system commands, result's std and err output
# - returns value compatible with Test::Harness
# ------------------------------------------------------------
	my $self    = shift;
	my $command = shift;

	open LOG_FILE,">>$self->{ logfile }" or die "Can't open log file '$self->{ logfile }' $!";
	print LOG_FILE "$command\n";

	# execute 
	my $output = `$command`;

	print LOG_FILE "$output\n";	
	close LOG_FILE;

	if( $command =~ /diff/ ) {
		if( $output ) { return 0; } 
		else          { return 1; }
	}

	# return 0 if NOT successful, 1 if is (for Test::Harness)
	my $success = ( $? >> 8 || $? & 127) ? 0 : 1;
	return $success;
}

1;
