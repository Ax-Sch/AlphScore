package Feature::File::Score;

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

	$self->{ file } = $file;
	my @contents = ();
	if( $file =~ /\.gz$/ ) {
		@contents = `env gzip -dc $file`;
	} else {
		open FILE, $file or die 
			"Score File '$file' cannot be opened for reading!\n" .
			"Check to see if the file exists and the proper permissions are set.\n$!";
		@contents = <FILE>;
		close FILE;
	}
	while( @contents ) {
		local $_ = shift @contents;
		chomp;
		next if /^$/;
		next if /^\s*#/; # Skip Metadata
		my ($environment_id, $score, $x, $y, $z, $tag, $comment) = split /\t/;
		push @{ $self->{ points }}, "$x,$y,$z";
		$self->{ contents }{ "$x,$y,$z" } = { 
			id => $environment_id, score => $score, 
			x => $x, y => $y, z => $z, tag => $tag, comment => $comment 
		};
	}
}

# ============================================================
sub diff {
# ============================================================
	my $self   = shift;
	my $file   = shift;
	my $target = new Feature::File::Score( $file );

	my @differences = 
		map { { point => $_, source => 'target' } } 
		grep { ! exists $self->{ contents }{ $_ } } @{ $target->{ points }};

	foreach my $point (@{$self->{ points }}) {
		if( exists $target->{ contents }{ $point } ) {
			my $a = $self->{ contents }{ $point }{ score };
			my $b = $target->{ contents }{ $point }{ score };
			
			# ===== NOTE THE DIFFERENCES, IF THERE ARE ANY
			push @differences, 	
				{ point => $point, source => 'self' }, 
				{ point => $point, source => 'target' }
				if( abs( $a - $b ) >= 0.01 );
		} else {
			push @differences, { point => $point, source => 'self' };
		}
	}

	# ===== PRINT DIFFERENCES
	foreach my $difference (@differences) {
		my $point  = $difference->{ point };
		my $source = $difference->{ source } eq 'self' ? $self : $target;
		my $arrow  = $difference->{ source } eq 'self' ? "<" : ">";
		my @values = ();
		push @values, $source->{ contents }{ $point }{ id };
		push @values, map { sprintf( "%0.2f", $source->{ contents }{ $point }{ $_ } ); } qw( score x y z );
		push @values, map { $source->{ contents }{ $point }{ $_ }; } qw( tag comment );
		print "$arrow ", join( "\t", @values ), "\n";
	}
	return 1 if( @differences );

	return 0;
}

1;
