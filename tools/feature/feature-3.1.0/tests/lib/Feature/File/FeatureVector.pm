package Feature::File::FeatureVector;

our $default_properties = default_properties();

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
	my $self            = shift;
	my $file            = shift;
	my $properties_file = shift;

	my @contents = ();
	if( $file =~ /\.gz$/ ) {
		@contents = `env gzip -dc $file`;
	} else {
		open FILE, $file or die
			"Cannot open feature vector file '$file' for reading.\n$!";
		@contents = <FILE>;
		close FILE;
	}

	if( defined $properties_file ) {
		$self->{ properties }      = get_properties( $properties_file );
		$self->{ properties_file } = $properties_file;
	} else {
		$self->{ properties }      = $default_properties;
	}

	while( @contents ) {
		local $_ = shift @contents;
		chomp;
		next if /^\s*#/; # Skip metadata 
		my @vector  = split /\t/;
		my $id      = shift @vector;
		my $comment = pop @vector;
		my $tag     = pop @vector;
		my $z       = pop @vector;
		my $y       = pop @vector;
		my $x       = pop @vector;
		my $hash    = pop @vector;
		warn "Error in parsing file '$file' $!" if( $hash ne '#' );

		push @{ $self->{ ids }}, $id;
		$self->{ contents }{ $id } = {
			id      => $id,
			vector  => \@vector,
			x       => $x,
			y       => $y,
			z       => $z,
			tag     => $tag,
			comment => $comment
		};
	}
}

# ============================================================
sub diff {
# ============================================================
	my $self            = shift;
	my $file            = shift;
	my $option          = shift || '';
	my $target          = new Feature::File::FeatureVector( $file, $self->{ properties_file } );
	my $properties      = $self->{ properties };

	my @differences =
		map { { id => $_, source => 'target', field => 'id' } }
		grep { ! exists $self->{ contents }{ $_ } } @{ $target->{ ids }};

	foreach my $id (@{ $self->{ ids }}) {
		if( exists $target->{ contents }{ $id } ) {
			$a = $self->{ contents }{ $id }{ vector };
			$b = $target->{ contents }{ $id }{ vector };

			my $different = 0;
			my @columns   = ();
			foreach my $i ( 0 .. $#$a ) {
				unless( within_tolerance( $a->[ $i ], $b->[ $i ] )) {
					push @columns, $i;
					$different = 1;
				}
			}
			push @differences,
				{ id => $id, source => 'self',   columns => [ @columns ] },
				{ id => $id, source => 'target', columns => [ @columns ] }
				if $different;
		} else {
			push @differences,
				{ id => $id, source => 'self',   field => 'id' };
		}
	}

	# ===== PRINT DIFFERENCES
	if( $option ne 'no-print' ) {
		foreach my $difference (@differences) {
			my $id     = $difference->{ id };
			my $source = $difference->{ source } eq 'self' ? $self : $target;
			my $arrow  = $difference->{ source } eq 'self' ? "<" : ">";
			my $entry  = $source->{ contents }{ $id };
			if( exists $difference->{ columns } ) {
				my $columns = $difference->{ columns };
				print "$arrow ", join( "\t", (
					$entry->{ id }, print_vector( $columns, $entry->{ vector }, $properties ), '#', 
					$entry->{ x }, $entry->{ y }, $entry->{ z }, 
					$entry->{ tag }, $entry->{ comment },
				)), "\n";
			} else {
				print "$arrow ", join( "\t", (
					$entry->{ id }, '...', '#', 
					$entry->{ x }, $entry->{ y }, $entry->{ z }, 
					$entry->{ tag }, $entry->{ comment },
				)), "\n";
			}
		}
	}
	return 1 if( @differences );

	return 0;
}

# ============================================================
sub print_vector {
# ============================================================
	my $columns    = shift;
	my $values     = shift;
	my $properties = shift;
	my $k          = int @$properties;
	my $n          = $#$values;
	my $result     = '';
	my @columns    = @$columns;

	while( @columns ) {
		my $i     = shift @columns;
		my $j     = defined $columns[ 0 ] ? $columns[ 0 ] - 1 : $i;
		my $l     = $i % $k;
		my $shell = int( $i / $k );
		my $field = $properties->[ $l ];
		my $value = sprintf( "%0.2f", $values->[ $i ] );
		if( $i == $n ) {
			$result .= "($field-$shell) $value";
		} else {
			$result .= "($field-$shell) $value\t";
		}
		if( $i != 0 && $i != $j && $i != $n ) {
			$result .= "...\t";
		}
	}
	return $result;
}

# ============================================================
sub default_properties {
# ============================================================
	my @properties = ();
	
	while( <DATA> ) {
		chomp;
		next if( /^$/ );
		push @properties, $_;
	}
	return \@properties;
}

# ============================================================
sub get_properties {
# ============================================================
	my $file       = shift;
	my $found      = 0;
	my @properties = ();
	my $path;
	foreach ( '.', $ENV{ FEATURE_DIR }, '/usr/local/feature/data' ) {
		if( -e "$_/$file" ) { $found = 1; $path = $_; last; }
	}
	die "Properties file '$file' not found in current working directory, $ENV{ FEATURE_DIR }, or the FEATURE default installation directory $!" if( ! $found );

	open FILE, "$path/$file" or die "Can't open '$path/$file' for reading\n";
	while( <FILE> ) {
		chomp;
		next if /^#/;
		next if /^$/;
		s/^\s*//g; s/\s*$//g;
		push @properties, $_;
	}
	close FILE;
	return \@properties;
}

# ============================================================
sub within_tolerance {
# ============================================================
	my $a = shift;
	my $b = shift;
	my $tolerance = shift || 0.1;

	return 0 unless defined $b;
	return ( abs( $a - $b ) < $tolerance );
}

1;

__DATA__
ATOM_TYPE_IS_C
ATOM_TYPE_IS_CT
ATOM_TYPE_IS_CA
ATOM_TYPE_IS_N
ATOM_TYPE_IS_N2
ATOM_TYPE_IS_N3
ATOM_TYPE_IS_NA
ATOM_TYPE_IS_O
ATOM_TYPE_IS_O2
ATOM_TYPE_IS_OH
ATOM_TYPE_IS_S
ATOM_TYPE_IS_SH
ATOM_TYPE_IS_OTHER
PARTIAL_CHARGE
ATOM_NAME_IS_ANY
ATOM_NAME_IS_C
ATOM_NAME_IS_N
ATOM_NAME_IS_O
ATOM_NAME_IS_S
ATOM_NAME_IS_OTHER
HYDROXYL
AMIDE
AMINE
CARBONYL
RING_SYSTEM
PEPTIDE
VDW_VOLUME
CHARGE
NEG_CHARGE
POS_CHARGE
CHARGE_WITH_HIS
HYDROPHOBICITY
MOBILITY
SOLVENT_ACCESSIBILITY
RESIDUE_NAME_IS_ALA
RESIDUE_NAME_IS_ARG
RESIDUE_NAME_IS_ASN
RESIDUE_NAME_IS_ASP
RESIDUE_NAME_IS_CYS
RESIDUE_NAME_IS_GLN
RESIDUE_NAME_IS_GLU
RESIDUE_NAME_IS_GLY
RESIDUE_NAME_IS_HIS
RESIDUE_NAME_IS_ILE
RESIDUE_NAME_IS_LEU
RESIDUE_NAME_IS_LYS
RESIDUE_NAME_IS_MET
RESIDUE_NAME_IS_PHE
RESIDUE_NAME_IS_PRO
RESIDUE_NAME_IS_SER
RESIDUE_NAME_IS_THR
RESIDUE_NAME_IS_TRP
RESIDUE_NAME_IS_TYR
RESIDUE_NAME_IS_VAL
RESIDUE_NAME_IS_HOH
RESIDUE_NAME_IS_OTHER
RESIDUE_CLASS1_IS_HYDROPHOBIC
RESIDUE_CLASS1_IS_CHARGED
RESIDUE_CLASS1_IS_POLAR
RESIDUE_CLASS1_IS_UNKNOWN
RESIDUE_CLASS2_IS_NONPOLAR
RESIDUE_CLASS2_IS_POLAR
RESIDUE_CLASS2_IS_ACIDIC
RESIDUE_CLASS2_IS_BASIC
RESIDUE_CLASS2_IS_UNKNOWN
SECONDARY_STRUCTURE1_IS_3HELIX
SECONDARY_STRUCTURE1_IS_4HELIX
SECONDARY_STRUCTURE1_IS_5HELIX
SECONDARY_STRUCTURE1_IS_BRIDGE
SECONDARY_STRUCTURE1_IS_STRAND
SECONDARY_STRUCTURE1_IS_TURN
SECONDARY_STRUCTURE1_IS_BEND
SECONDARY_STRUCTURE1_IS_COIL
SECONDARY_STRUCTURE1_IS_HET
SECONDARY_STRUCTURE1_IS_UNKNOWN
SECONDARY_STRUCTURE2_IS_HELIX
SECONDARY_STRUCTURE2_IS_BETA
SECONDARY_STRUCTURE2_IS_COIL
SECONDARY_STRUCTURE2_IS_HET
SECONDARY_STRUCTURE2_IS_UNKNOWN
