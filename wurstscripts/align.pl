#!/usr/bin/perl

#use lib "$ENV{HOME}/pl/lib/i586-linux-thread-multi";  # Where wurst lives after installation
#use lib "/home/stud2004/tmargraf/pl/lib/i686-linux-thread-multi";

use FindBin;
use lib "$FindBin::Bin/../wurst_new/src/Wurst/blib/arch";
use lib "$FindBin::Bin/../wurst_new/src/Wurst/blib/lib";

use Wurst;
use strict;
use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);

# ----------------------- Defaults  ---------------------------------
# These are numbers you might reasonably want to change.
# They should (will) be changeable by options.
use vars qw ($modeldir $DFLT_MODELDIR $pvecdir $classfile);
*DFLT_MODELDIR = \'modeldir';
$modeldir      = $DFLT_MODELDIR;
$pvecdir       = 'pvecs';
$classfile     = 'classfile';

# ----------------------- Alignment Constants -----------------------
# These are declared globally, and set by set_params().
use vars qw (
  $align_type
  $sw1_pgap_open
  $sw1_qgap_open
  $sw1_pgap_widen
  $sw1_qgap_widen
  $m_shift
  $gauss_err
);

# These parameters will be used for extending alignments via a
# Needleman and Wunsch

use vars qw (
  $nw_pgap_open
  $nw_qgap_open
  $nw_pgap_widen
  $nw_qgap_widen
  $nw_sec_pnlty
);

# ----------------------- set_params    -----------------------------
# This gets its own function because it can be more complicated
# if, in the future, we have a version depending on various
# options like whether or not we have secondary structure
# information.
# pgap controls penalties in the sequence, qgap in the structure.

sub set_params () {
    *sw1_pgap_open  = \3.25;
    *sw1_pgap_widen = \0.8942;
    *sw1_qgap_open  = \$sw1_pgap_open;
    *sw1_qgap_widen = \$sw1_pgap_widen;
    *m_shift        = \-0.1;
    *gauss_err      = \0.4;

    *nw_pgap_open  = \$sw1_pgap_open;
    *nw_qgap_open  = \$sw1_qgap_open;
    *nw_pgap_widen = \$sw1_pgap_widen;
    *nw_qgap_widen = \$sw1_qgap_widen;
}

# ----------------------- usage   -----------------------------------
sub usage () {
    print STDERR "Usage: \n    $0 struct_file_A struct_file_B \n";
    exit(EXIT_FAILURE);
}

# ----------------------- get_scores --------------------------------
sub get_scores ($ $ $ $) {
    my ( $pair_set, $coord1, $coord2, $to_use ) = @_;

    my ( $scr_tot, $coverage, $score1, $geo_gap, $score1_gap );
    my ( $str1, $crap );
    my ( $open_cost, $widen_cost, $nseq_gap );
    ( $score1_gap, $score1 ) = pair_set_score($pair_set);
    ( $str1,       $crap )   =
      pair_set_coverage( $pair_set, coord_size($coord1), coord_size($coord2) );
    $coverage = ( $str1 =~ tr/1// );    # This is coverage as an integer
    $coverage =
      $coverage / seq_size( coord_get_seq($coord1) )
      ;                                 #and as fraction of query structure

    my ( $k_scr2, $k_gap_geo, $k_seq_gap, $k_str_gap, $k_str_wdn );

    if ( $coverage < .05 ) {
        $geo_gap    = 0;
        $nseq_gap   = 0;
        $open_cost  = 0;
        $widen_cost = 0;
    }
    else {
        ( $open_cost, $widen_cost ) = pair_set_gap( $pair_set, 1, 1 );
    }

    $k_str_gap = 1 * $sw1_pgap_open;
    $k_str_wdn = 1 * $sw1_pgap_widen;

    $scr_tot = $score1 + $k_str_gap * $open_cost + $k_str_wdn * $widen_cost;
    return ( $scr_tot, $coverage, $score1, $score1_gap, $nseq_gap, $open_cost );
}


# ----------------------- get_alt_scores ---------------------------------
# calculates scores on random paths through the scoring matrix
# parameters: number_of_paths/scores, scoring_matrix, pair_set_of_optimal_path
# return: the scores
sub get_alt_scores($ $ $) {
    my ( $num_scrs, $scr_mat, $pair_set ) = @_;
    my @scr_fin;

    for ( my $i = 0 ; $i < $num_scrs ; $i++ ) {
        $scr_fin[$i] = find_alt_path_score_simple( $scr_mat, $pair_set );
    }

    return \@scr_fin;
}

# ----------------------- normalize_alt_scores ------------------------------
#
sub normalize_alt_scores($) {
    my ($scrs_ref) = @_;
    my $mean = 0.0;

    foreach my $scr ( @{$scrs_ref} ) {
        $mean += $scr;
    }
    $mean /= @$scrs_ref;

    my $deviation = 0.0;

    foreach my $scr (@$scrs_ref) {
        my $tmp = $scr - $mean;
        $deviation += ( $tmp * $tmp );
    }
    $deviation /= ( @$scrs_ref - 1 );
    $deviation = sqrt($deviation);

    return ( $mean, $deviation );
}

# ----------------------- get_dme_thresh ----------------------------
# Given an alignment and two proteins, return the thresholded DME
sub get_dme_thresh ($ $ $) {
    my ( $pair_set, $c1, $c2 ) = @_;
    my $model = make_model( $pair_set, coord_get_seq($c1), $c2 );
    if (! $model) {
        print STDERR "Found a broken model ", coord_name($c1), ' ', coord_name($c2), "\n";
        return 0.0;
    }
    if ( coord_size($model) < 10 ) {
        return 0.0;
    }
    my $frac;
    if ( !dme_thresh( $frac, $c1, $model, 3.0 ) ) {
        print STDERR "dme_thresh broke.\n";
    }
    return ($frac);
}

# ----------------------- do_align ----------------------------------
# This does the alignment. Although it takes secondary structure
# data as an argument, we are not yet using this in our server
# script.
sub do_align ($ $ $) {
    my ( $coord1, $coord2, $formatflag ) = @_;
    my ($pvec1, $pvec2);
    my $gauss_err = 0.4;
    my $classfcn = aa_strct_clssfcn_read($classfile, $gauss_err);
	$pvec1 = strct_2_prob_vec($coord1, $classfcn, 1);
    $pvec2 = strct_2_prob_vec($coord2, $classfcn, 1);
	my $minsize;
	if(seq_size(coord_get_seq($coord1)) < seq_size(coord_get_seq($coord2)) ){
		$minsize = seq_size(coord_get_seq($coord1));
	}
	else{
		$minsize = seq_size(coord_get_seq($coord2));
	}
	my $seq_ptr1 = coord_get_seq($coord1);
    my $seq_ptr2 = coord_get_seq($coord2);

    # build score matrix
    my $matrix = score_mat_new( prob_vec_length($pvec1), prob_vec_length($pvec2) );
	# fill the score matrix
    score_pvec( $matrix, $pvec1, $pvec2 );
	# adjust the 0-level of the matrix
    $matrix = score_mat_shift( $matrix, $m_shift );

    my (
        $sw_scr_tot,    $sw_coverage, $sw_score1,
        $sw_score1_gap, $sw_seq_gap,  $sw_strct_gap
    );
    my (
        $nw_scr_tot,    $nw_coverage, $nw_score1,
        $nw_score1_gap, $nw_seq_gap,  $nw_strct_gap
    );
	# sum up the matrix and do the traceback (Smith Waterman style)
    my $sw_pair_set = score_mat_sum_smpl(
        my $crap_mat,   $matrix,         $sw1_pgap_open, $sw1_pgap_widen,
        $sw1_qgap_open, $sw1_qgap_widen, $S_AND_W
    );
	# now, our alignment is complete. The rest of this function is just
	# icing on the cake. From here on, we just apply various scoring schemes
	# to our alignment and generate readable output.
    (
        $sw_scr_tot,    $sw_coverage, $sw_score1,
        $sw_score1_gap, $sw_seq_gap,  $sw_strct_gap
      )
      = get_scores( $sw_pair_set, $coord1, $coord2, 's_and_w' );
	# now we have the alignment score and the coverage.
    my $frac_dme = get_dme_thresh( $sw_pair_set, $coord1, $coord2 );
	# this gave us a distance matrix based similarity score.
    my $num_scrs     = 1000;
	# here we generate alternative alignments from random tracebacks 
	# through our DP-matrix and calculate scores for those alignments.
    my $alt_scrs_ref = get_alt_scores( $num_scrs, $matrix, $sw_pair_set );
    my ( $mean, $deviation ) = normalize_alt_scores($alt_scrs_ref);
    undef($alt_scrs_ref);
	# from the distribution of the scores of the random alignments, we
	# estimate a pseudo-Z-score.
    my $z_scr;
    if ( $deviation != 0 ) {
        $z_scr = ( $sw_scr_tot - $mean ) / $deviation;
    }
    else {
        $z_scr = 0.0;
    }    # Should not really happen

    #   If the alignment is tiny, one can get a ridiculous z-score
    if ( $sw_coverage < 0.03 ) {    # silently wipe these guys out
        $z_scr = 0.0;
    }

	# not used here. we only calculate local alignments. if we want 
	# global alignments, its as easy as changing the last parameter 
	# of the get_scores() and score_mat_sum_simple() calls.
    my $nw_pair_set;
    
	# here it gets a bit nasty. the @r list is completely unnecessary.
	# it is a legacy from the wurst server script and has been introduced
	# here only to be able to transplant the output formatting code. Pay no 
	# attention to the man behind the curtain...
    my @r = (
        $sw_scr_tot, $nw_pair_set, $nw_scr_tot,  $sw_coverage, $nw_coverage,
        $sw_score1,  $nw_score1,   $sw_pair_set, $z_scr,       $frac_dme
    );

	# This is completely unreadable. Perl at its worst. Sorry.
	# Here we get the sequence identity.
    $r[4] =
      get_seq_id_simple( $r[7], coord_get_seq($coord2),
        coord_get_seq($coord1) );
    $r[6] = $r[3] * seq_size( coord_get_seq($coord2) );
    if ( $r[6] < 25 ) {
        $r[8] = 0;
    }
    $r[2] = $r[9] * ( ( $r[6] ) / $minsize );

    print "\n",
"sw refers to Smith and Waterman alignment, nw refers to Needleman & Wunsch\n",
      "z scr : the z-score of the alignment with 1000 alternative alignments\n",
"sw scr: the combined results of score function + gaps for sw alignment\n",
      "sw cvr: coverage (fraction of sequence) accounted for by sw alignment\n",
      "nw cvr: coverage                        accounted for by nw alignment\n",
      "asize : length of the alignment\n",
      "rmsd  : root-mean-square-distance between aligned residues\n",
"___________________________ Alignment Scores _________________________________\n";

    printf "%8s %8s %8s %6s %6s %8s %6s %6s\n", 'z scr', 'f_dme',
      'sw scr', 'sw cvr', 'seq id', 'asize', 'rmsd', 'andr_scr';
    my $pair_set = $r[7];
	# This function superimposes the two structures and calculates the RMSD.
    ( $r[1], $r[10], my $crap_b ) =
      coord_rmsd( $r[7], $coord2, $coord1, 0 );

	# here we calculate our pseudo-Q-score (formerly known as andrew score).  
    my $q =
      ( $r[6] * $r[6] ) / ( ( 1 + ( $r[1] / 3 ) * ( $r[1] / 3 ) ) *
          seq_size( coord_get_seq($coord1) ) *
          seq_size( coord_get_seq($coord2) ) );
    my $sid;
    if ( $r[6] != 0 ) {
        $sid = $r[4] / $r[6];
    }
    else {
        $sid = 0.0;
    }
    printf "%8.3g %8.3g %8.1f %6.2f %6.2g %8.4g %6.2f %4.4f\n",
      $r[8], $r[9], $r[0], $r[6]/seq_size(coord_get_seq($coord1)), $sid, $r[6],
      $r[1], $r[2];
    print "\n";
   
	# This prints a nicely formatted alignment.
    print(pair_set_pretty_string($r[7], coord_get_seq($coord1), coord_get_seq($coord2)));
	print "writing models...\n";
	# and here we write the superimposed structures to disc.
	# this should become a command line switch. also, output
	# filenames should be something more meaningful.
	# uncomment if you want to look at the structures.
    coord_2_pdb("$ARGV[0].wurst", $coord1);
    coord_2_pdb("$ARGV[1].wurst", $r[10]);
    return EXIT_SUCCESS;
}

# ------------------- bad_exit -------------------------------------
# Kill the script and print an error message
sub bad_exit ( $ ) {
	my $msg = shift;
	print STDERR "Error: \"$msg\"\n";
	exit(EXIT_FAILURE);
}

# ----------------------- catch_kill     ----------------------------
# The main thing is, if we get a KILL or TERM, to call exit and get
# out of here. This means there is a better chance of closing files
# wherever we were up to.
sub catch_kill
{
    my ($sig) = @_;
    bad_exit ("signal $sig received");

}

# ----------------------- mymain  -----------------------------------
# Arg 1 is a structure file. Arg 2 is a structure list file.
sub mymain () {
    my $fatalflag = undef;

    if(!$ARGV[1]){
        print STDERR "Must have two structure files.\n";
		usage();
    }
    set_params();
    my $structA;
# Read the pdb-files specified on the commandline
    my $structA_name = $ARGV[0];
    if ( $structA_name =~ /\.pdb/ ) {
        $structA = pdb_read( $structA_name, '', '' );
        $structA_name =~ s/\.pdb//;
    }
    
	my $structB;
	my $structB_name = $ARGV[1];
    if ( $structB_name =~ /\.pdb/ ) {
        $structB = pdb_read( $structB_name, '', '' );
        $structB_name =~ s/\.pdb//;
    }

    print "Structures read... \n ";
    my $formatflag = 0;

	# This function does the real work. The formatflag is unused.
	# I'm keeping it so I can tell the function to write structures,
	# do global alignments instead of local ones, etc.
    my $r =
      do_align( $structA, $structB, $formatflag );
        
    if ( $r == EXIT_FAILURE ) {
        bad_exit('calculation broke');
    }

	# most people don't care about this. -another legacy from the server script.
	# the times are nice for benchmarking though.
    print
"__________________________________________________________________________\n",
      "Wurst gegessen at ", scalar( localtime() ), "\n";
    my ( $user, $system, $crap, $crap2, $host );
    ( $user, $system, $crap, $crap2 ) = times();
    printf "I took %d:%d min user and %.0f:%.0f min sys time\n", $user / 60,
      $user % 60, $system / 60, $system % 60;
    use Sys::Hostname;
    $host = hostname() || { $host = 'no_host' };
    print "Run on $host\n";
    return EXIT_SUCCESS;
}

# ----------------------- main    -----------------------------------
exit( mymain() );
