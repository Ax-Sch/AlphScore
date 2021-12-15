package Feature::File::Model;

our $tolerance = 0.01;

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
	my $self = shift;
	my $file = shift;

	my @contents = ();

	if( $file =~ /\.gz$/ ) {
		@contents = `env gzip -dc $file`;
	} else {
		open FILE, $file or die
			"Cannot open model file '$file' for reading.\n$!";
		@contents = <FILE>;
		close FILE;
	}

	while( @contents ) {
		local $_ = shift @contents;
		chomp;
		if( /^#\s+(\w+)\s+(.+)$/ ) {
			my $header = lc $1;
			my $value  = $2;
			$self->{ header }{ $header } = $value;
		} else {
			my ($property, $p_value, $min, $range, @bins) = split /\t/;
			push @{ $self->{ property }}, $property;
			$self->{ model }{ $property } = {
				property => $property, p => $p_value, min => $min,  
				range => $range, bins => [ @bins ]
			};
		}
	}
}

# ============================================================
sub diff {
# ============================================================
	my $self   = shift;
	my $file   = shift;
	my $target = new Feature::File::Model( $file );

	my @differences = ();

	# ===== CHECK IF MODEL FILES ARE COMPATIBLE
	push @differences, 
		map { { property => $_, source => 'target' }; } 
		grep { ! exists $self->{ model }{ $_ } } @{ $target->{ property }};
	push @differences, 
		map { { property => $_, source => 'self' }; } 
		grep { ! exists $target->{ model }{ $_ } } @{ $self->{ property }};
	if( @differences ) {
		print "Different model properties!\n";
		foreach my $difference (@differences) {
			my $arrow = $difference->{ source } eq 'self' ? '<' : '>';
			print "$arrow $difference->{ property }";
		}
		return 1;
	}

	# ===== CHECK IF MODEL VALUES ARE EQUIVALENT
	foreach my $property (@{ $self->{ property }}) {
		my $different = 0;
		$a = $self->{ model }{ $property };
		$b = $target->{ model }{ $property };

		# ===== COMPARE NUMBER OF BINS (TYPICALLY 5 EACH)
		my $a_bins = $a->{ bins };
		my $b_bins = $b->{ bins };
		unless ( $#$a_bins == $#$b_bins ) {
			$different = 1;
		}

		# ===== COMPARE BIN VALUES
		foreach my $i ( 0 .. $#$a_bins ) {
			unless( within_tolerance( $a_bins->[ $i ], $b_bins->[ $i ] )) {
				$different = 1;
			}
		}

		# ===== COMPARE P-VALUE, MIN, AND RANGE
		FIELD: foreach my $field ( qw( p min range )) {
			unless( within_tolerance( $a->{ $field }, $b->{ $field } )) {
				$different = 1;
				last FIELD;
			}
		}
		if( $different ) {
			push @differences, { property => $property, source => 'self' };
			push @differences, { property => $property, source => 'target' };
		}
	}

	# ===== PRINT DIFFERENCES
	foreach my $difference (@differences) {
		my $property = $difference->{ property };
		my $source   = $difference->{ source } eq 'self' ? $self : $target;
		my $arrow    = $difference->{ source } eq 'self' ? "<" : ">";
		my @values   = map { sprintf( "%0.02f", $source->{ model }{ $property }{ $_ } ); } qw( p min range );
		print "$arrow $property ", join( "\t", @values, @{ $source->{ model }{ $property }{ bins }} ), "\n";
	}

	return 1 if( @differences );

	return 0;
}

# ============================================================
sub within_tolerance {
# ============================================================
	my $a = shift;
	my $b = shift;

	return ( abs( $a - $b ) < $tolerance );
}

1;
