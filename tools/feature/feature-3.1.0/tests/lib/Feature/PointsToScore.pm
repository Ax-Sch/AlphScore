package Feature::PointsToScore;

=head1 NAME

  Feature::PointsToScore

=head1 DESCRIPTION

Utilities to sample PDBs to score using a given model.

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
	my $feature             = shift;
	my $model               = shift;
	my $database            = shift;
	$self->{ model }        = $model;
	$self->{ input }        = $feature->get_input();
	$self->{ output }       = $feature->get_output();
	$self->{ database }     = $database;
	$self->{ atomselector } = $feature->get_atomselector();
	($self->{ name }, $self->{ position }, $self->{ residue }, $self->{ atom }) = split /\./, $model;
}

# ============================================================
sub default {
# ============================================================
	my $self = shift;
	my @pdbs  = ();

	open FILE, "$self->{ input }/$self->{ model }/pdbids-to-score.txt" or die $!;
	while( <FILE> ) {
		chomp;
		my ($pdbid) = /^(\w{4})/;
		push @pdbs, $pdbid;
		`$self->{ atomselector } -r $self->{ residue } -a $self->{ atom } $pdbid > $self->{ output }/$self->{ model }/$pdbid.ptf`;
	}
	return @pdbs;
}

# ============================================================
sub resample {
# ============================================================
	my $self   = shift;

	print STDERR "      Selecting PDBs...\n";
	my @pdbs = $self->_select_pdbs();
	open FILE, ">$self->{ input }/$self->{ model }/pdbids-to-score.txt";
	print FILE map { "$_\n"; } @pdbs;
	close FILE;
	return @pdbs;
}

# ============================================================
sub _select_pdbs {
# ============================================================
	my $self       = shift;
	my @total_pdbs = $self->_get_pdbs();
	my @pdbs       = ();
	my $remaining  = {};
	my @list       = ();
	foreach my $pdb (@total_pdbs) { $remaining->{ $pdb } = 1; }

	# ===== GET SITES
	@list = ();
	open FILE, "$self->{ input }/$self->{ model }/sites.ptf" or die $!;
	while( <FILE> ) {
		my ($pdbid) = split /\t/;
		delete $remaining->{ $pdbid };
		push @list, $pdbid;
	}
	close FILE;
	push @pdbs, _select_pdbs_from_list( 5, \@list, $residue, $atom );

	# ===== GET NONSITES
	@list = ();
	open FILE, "$self->{ input }/$self->{ model }/nonsites.ptf" or die $!;
	while( <FILE> ) {
		my ($pdbid) = split /\t/;
		delete $remaining->{ $pdbid };
		push @list, $pdbid;
	}
	close FILE;
	push @pdbs, _select_pdbs_from_list( 5, \@list, $residue, $atom );
	
	# ===== GET RANDOM SITES
	@list = sort keys %$remaining;
	push @pdbs, _select_pdbs_from_list( 10, \@list, $residue, $atom );

	return @pdbs;
}

# ============================================================
sub _select_pdbs_from_list {
# ============================================================
	my $number  = shift;
	my $list    = shift;
	my $residue = shift;
	my $atom    = shift;
	my @pdbs    = ();

	for my $i ( 1 .. $number ) {
		my $selection_succeeded = 0;
		do {
			my $j = int( rand( $#$list + 1 ));
			my $pdb = splice( @$list, $j, 1 );
			print STDERR "      Trying $pdb ($i of $number)... ";
			my $points = `$atomselector -r $residue -a $atom $pdb`;
			if( $points ) {
				print STDERR "ok\n";
				$selection_succeeded = 1;
				push @pdbs, $pdb;
				open POINTFILE, ">$pdb.ptf" or die $!;
				print POINTFILE $points;
				close POINTFILE;

			} else {
				print STDERR "no good, trying a different PDB\n";
				$selection_succeeded = 0;
			}
		} while( ! $selection_succeeded );
	}
	return @pdbs;
}

# ============================================================
sub _get_pdbs {
# ============================================================
	my $self = shift;
	opendir DIR, $self->{ database } or die $!;
	my @pdbs = map { /^(\w{4})/; $1; } grep { /\.pdb\.gz$/; } readdir DIR;
	closedir DIR;
	die "Can't find any PDB files at '$self->{ database }' $!" unless @pdbs;
	return @pdbs
}

1;
