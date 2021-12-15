package Feature::InstallTest;
# ============================================================
# by Trevor Blackstone, San Francisco State University
# Adapted from Marc Sosnick's release candidate test system
# ============================================================

use Test::Harness;
use Cwd qw( getcwd );

our $cwd                    = getcwd();

# $log_filename has name of log file output.  clear log.
# $log_must be absolute filename for test() below
our $log_filename = "$cwd/installtest.log";

# ============================================================
sub run { 
# ============================================================
# - executes system command
# - logs system commands, result's std and err output
# - returns value compatible with Test::Harness
# ------------------------------------------------------------
	my $command = shift;
	# log the line we're going to execute
	open(LOG_FILE,">>$log_filename") or die("Cannot open log file");
	print LOG_FILE "$command\n\n";

	# execute 
	$output = `$command`;

	print LOG_FILE "$output\n\n";	
	close LOG_FILE;

	if( $command =~ /diff/ ) {
		if( $output ) {
			return 0;
		} else {
			return 1;
		}
	}

	# return 0 if NOT successful, 1 if is (for Test::Harness)
	my $success = ( $? >> 8 || $? & 127) ? 0 : 1;
	return $success;
}

1;
