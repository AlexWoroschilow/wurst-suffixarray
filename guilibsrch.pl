#!/usr/bin/perl
#use lib "$ENV{HOME}/pl/lib/i586-linux-thread-multi";  # Where wurst lives after installation
#use lib "/home/stud2004/tmargraf/pl/lib/i686-linux-thread-multi";
use FindBin;
use lib "$FindBin::Bin/../../wurst/src/Wurst/blib/arch";
use lib "$FindBin::Bin/../../wurst/src/Wurst/blib/lib";

# ----------------------- wxPerl GUI -----------------------------------
package MyFrame;
use 5.008;
use Wx;
use Wx::Grid;
use base qw(Wx::Frame);
use Wx::Event qw(EVT_BUTTON EVT_GRID_SELECT_CELL);
use Wx qw( :everything );
use strict;
use threads;
use threads::shared;
use Config;
my $filebox;
my $results;
my $progress;
my $alignment;
my $q_path;
my $query_struct; 
my @struct_list;
my @struct_dirs;
my $r : shared;
my @structures : shared;
my $libthread;

#-----------------------Wurst structure search---------------------------

use Wurst;
use vars qw ($MATRIX_DIR $PARAM_DIR
             $RS_PARAM_FILE $FX9_PARAM_FILE );

#do "$ENV{HOME}/../../torda/c/wurst/scripts/paths.inc" || die $@;

$Config{useithreads} or die "Recompile Perl with threads to run this program.";

if ($@) {
    die "broke reading paths.inc:\n$@"; }
if ( defined ($ENV{SGE_ROOT})) {
    $MATRIX_DIR = "$ENV{HOME}/../../torda/c/wurst/matrix";
    $PARAM_DIR  = "$ENV{HOME}/../../torda/c/wurst/params";
}

#use strict;

use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);


# ----------------------- Defaults  ---------------------------------
# These are numbers you might reasonably want to change.
# They should (will) be changeable by options.
use vars qw ($N_BRIEF_LIST $N_MODELS $N_ALIGNMENTS);
$N_BRIEF_LIST   = 100;
$N_MODELS       = 0;
$N_ALIGNMENTS   = 0;
use vars qw ($modeldir $DFLT_MODELDIR $pvecdir);
*DFLT_MODELDIR =\'modeldir';
$modeldir = $DFLT_MODELDIR;
$pvecdir  = '/Users/tom/pvecs';

# ----------------------- Sequence Alignment Constants -----------------------
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

use vars qw( @DFLT_STRUCT_DIRS  @PROFILE_DIRS
             $phd_suffix $bin_suffix $prof_suffix);
*DFLT_STRUCT_DIRS = ['.', '/Users/tom/pdb90.bin'];
*PROFILE_DIRS     = ["/bm/pdb90_prof"];
*bin_suffix       = \'.bin';
*prof_suffix      = \'.prof';
*phd_suffix       = \'.phd';

# ----------------------- set_params    -----------------------------
# This gets its own function because it can be more complicated
# if, in the future, we have a version depending on various
# options like whether or not we have secondary structure
# information.
# pgap controls penalties in the sequence, qgap in the structure.

sub set_params ()
{
    *sw1_pgap_open  =  \  3.182;
    *sw1_pgap_widen =  \  0.8727;
    *sw1_qgap_open  =  \  $sw1_pgap_open;
    *sw1_qgap_widen =  \  $sw1_pgap_widen;
    *m_shift        =  \  -0.4954;
    *gauss_err      =  \  0.4;

    *nw_pgap_open  =   \  $sw1_pgap_open;
    *nw_qgap_open  =   \  $sw1_qgap_open;
    *nw_pgap_widen =   \  $sw1_pgap_widen;
    *nw_qgap_widen =   \  $sw1_qgap_widen;
}



# ----------------------- get_prot_list -----------------------------
# Go to the given filename and get a list of proteins from it.
sub get_prot_list ($)
{
    my $f = shift;
    my @a;
    if ( ! open (F, "<$f")) {
        print STDERR "Open fail on $f: $!\n";
        return undef;
    }

    while (my $line = <F>) {
        chomp ($line);
        my @words = split (' ', $line);
        if (! defined $words[0]) { next;}
        $line = $words[0];
        $line =~ s/#.*//;            # Toss comments away
        $line =~ s/\..*//;           # Toss filetypes away
        $line =~ s/^ +//;            # Leading and
        $line =~ s/ +$//;            # trailing spaces.
        if ($line eq '') {
            next; }
        substr ($line, 0, 4) = lc (substr ($line, 0, 4)); # 1AGC2 to 1agc2
        if (length ($line) == 4) {  # Convert 1abc to 1abc_
            $line .= '_'; }
        push (@a, $line);
    }
    close (F);
    return (@a);
}

# ----------------------- get_path  ---------------------------------
# We have a filename and a list of directories where it could
# be. Return the path if we can find it, otherwise return undef.
sub get_path (\@ $)
{
    my ($dirs, $fname) = @_;
    foreach my $d (@$dirs) {
        my $p = "$d/$fname";
        if ( -f $p) {
            return $p; }
    }
    return undef;
}
# ----------------------- check_dirs --------------------------------
# Given an array of directory names, check if each one
# exists. Print something if it is missing, but do not give
# up. It could be that there is some crap in the command line
# args, but all the important directories are really there.
# This function is potentially destructive !
# If a directory does not seem to exist, we actually remove it
# from the array we were passed.  This saves some futile lookups
# later on.
sub check_dirs (\@)
{
    my $a = shift;
    my $last = @$a;
    for (my $i = 0; $i < $last; $i++) {
        if ( ! -d $$a[$i]) {
            print STDERR "$$a[$i] is not a valid directory. Removing\n";
            splice @$a, $i, 1;
            $last--;
            $i--;
        }
    }
}

# ----------------------- check_files -------------------------------
# We are given an array of directories and and array of protein
# names and an extension.
# Check if all the files seem to be there.
sub check_files (\@ \@ $)
{
    my ($dirs, $fnames, $ext) = @_;
    my $errors = 0;
    foreach my $f (@$fnames) {
        my $name = "$f$ext";
        if (! get_path (@$dirs, $name)) {
            $errors++;
            print STDERR "Cannot find $name\n";
        }
    }
    return $errors;
}

# ----------------------- usage   -----------------------------------
sub usage ()
{
    print STDERR "Usage: $0 seq_file struct_file \n";
    exit (EXIT_FAILURE);
}

# ----------------------- get_scores --------------------------------
sub get_scores ($ $ $ $)
{
    my ($pair_set, $coord1, $coord2, $to_use) = @_;

    my ($scr_tot, $coverage, $score1, $geo_gap, $score1_gap);
    my ($str1, $crap);
    my ($open_cost, $widen_cost, $nseq_gap);
    ($score1_gap, $score1) = pair_set_score( $pair_set );
    ($str1, $crap) =
        pair_set_coverage ($pair_set, coord_size ($coord1), coord_size ($coord2));
    $coverage = ($str1 =~ tr/1//);       # This is coverage as an integer
    $coverage = $coverage / seq_size (coord_get_seq($coord1)); 
    #and as fraction of query structure

    my ($k_scr2, $k_gap_geo, $k_seq_gap, $k_str_gap, $k_str_wdn);


    if ($coverage  < .05 ) {
        $geo_gap = 0;
        $nseq_gap = 0;
        $open_cost = 0;
        $widen_cost = 0;
    } else {
        ($open_cost, $widen_cost) = pair_set_gap($pair_set, 1, 1);
    }

    $k_str_gap = 2 * $sw1_pgap_open;
    $k_str_wdn = 2 * $sw1_pgap_widen;

    $scr_tot = $score1 +
               $k_str_gap * $open_cost +
               $k_str_wdn * $widen_cost;
    return ( $scr_tot, $coverage, $score1, $score1_gap,
             $nseq_gap, $open_cost);
}


# ----------------------- bad_exit ----------------------------------
# This will run in a server, so if something goes wrong, we
# should at least mail back an indication.  The single parameter
# should be the error message returned by the function which was
# unhappy.
# Should we print to stderr or stdout ?
# This should not matter since we have grabbed both file handles.
sub bad_exit ( $ )
{
    my $msg = shift;
    print STDERR "Error: \"$msg\"\n";
    exit (EXIT_FAILURE);
}

# ----------------------- get_alt_scores ---------------------------------
# calculates scores on random paths through the scoring matrix
# parameters: number_of_paths/scores, scoring_matrix, pair_set_of_optimal_path
# return: the scores
sub get_alt_scores($ $ $)
{
   my ($num_scrs, $scr_mat, $pair_set) = @_;
   my @scr_fin;

   for (my $i = 0; $i < $num_scrs; $i++) {
       $scr_fin[$i] = find_alt_path_score_simple ($scr_mat, $pair_set); }

   return \@scr_fin;
}

# ----------------------- normalize_alt_scores ------------------------------
#
sub normalize_alt_scores($)
{
   my ($scrs_ref) = @_;
   my $mean = 0.0;

   foreach my $scr (@{$scrs_ref}) {
       $mean += $scr }
   $mean /= @$scrs_ref;

   my $deviation = 0.0;

   foreach my $scr (@$scrs_ref) {
       my $tmp = $scr - $mean;
       $deviation += ($tmp * $tmp);
   }
   $deviation /= (@$scrs_ref - 1);
   $deviation = sqrt($deviation);

   return ($mean, $deviation);
}

# ----------------------- get_dme_thresh ----------------------------
# Given an alignment and two proteins, return the thresholded DME
sub get_dme_thresh ($ $ $)
{
    my ($pair_set, $c1, $c2) = @_;
    my $model = make_model ($pair_set, coord_get_seq ($c1), $c2);

    if ($model && coord_size($model) < 10) { #  
	return 0.0;
    }
    my $frac;
    if (!$model || ! dme_thresh ($frac, $model, $c1, 4.0)){ #  
	print STDERR "dme_thresh broke.\n"; 
	$frac = 0.0;
    }
    return ($frac);
}

# ----------------------- do_align ----------------------------------
# This does the alignment. Although it takes secondary structure
# data as an argument, we are not yet using this in our server
# script.
sub do_align ($ $ $ $) {
    my ( $coord1, $coord2, $pvec1, $pvec2 ) = @_;
    my $seq_ptr1 = coord_get_seq($coord1);
    my $seq_ptr2 = coord_get_seq($coord2);

    # build score matrix
    my $matrix = score_mat_new( seq_size($seq_ptr1), seq_size($seq_ptr2) );

    score_pvec( $matrix, $pvec1, $pvec2 );
    $matrix = score_mat_shift( $matrix, $m_shift );

    my (
        $sw_scr_tot,    $sw_coverage, $sw_score1,
        $sw_score1_gap, $sw_seq_gap,  $sw_strct_gap
    );
    my (
        $nw_scr_tot,    $nw_coverage, $nw_score1,
        $nw_score1_gap, $nw_seq_gap,  $nw_strct_gap
    );

    my $sw_pair_set = score_mat_sum_smpl(
        my $crap_mat,   $matrix,         $sw1_pgap_open, $sw1_pgap_widen,
        $sw1_qgap_open, $sw1_qgap_widen, $S_AND_W
    );

    (
        $sw_scr_tot,    $sw_coverage, $sw_score1,
        $sw_score1_gap, $sw_seq_gap,  $sw_strct_gap
      )
      = get_scores( $sw_pair_set, $coord1, $coord2, 's_and_w' );
    my $frac_dme = get_dme_thresh( $sw_pair_set, $coord1, $coord2 );

    my $num_scrs     = 1000;
    my $alt_scrs_ref = get_alt_scores( $num_scrs, $matrix, $sw_pair_set );

    my ( $mean, $deviation ) = normalize_alt_scores($alt_scrs_ref);

    undef($alt_scrs_ref);
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

    my $sec_strct_pnlty = 0.0;
    my $nw_pair_set;
    
    my @r = (
        $sw_scr_tot, $nw_pair_set, $nw_scr_tot,  $sw_coverage, $nw_coverage,
        $sw_score1,  $nw_score1,   $sw_pair_set, $z_scr,       $frac_dme
    );
    return ( \@r );
}

# ----------------------- zero_shift_mat ----------------------------
# This shifts a substitution matrix, not an alignment score matrix.
sub zero_shift_mat ($ $) {
    my ( $sub_mat, $shift ) = @_;
    for ( my $i = 0 ; $i < 20 ; $i++ ) {
        for ( my $j = $i ; $j < 20 ; $j++ ) {
            my $t = sub_mat_get_by_i( $sub_mat, $i, $j ) + $shift;
            sub_mat_set_by_i( $sub_mat, $i, $j, $t );
        }
    }
}

# ----------------------- do_lib  -----------------------------------
# Walk over a library, doing alignments and saving interesting
# scores. The definition of interesting is a bit arbitrary.
# There is one very non-obvious coding trick.  We need to be able
# to pass the score information into the sorting functions. We
# could put everything into a big, two-dimensional array, but we
# can avoid copying data. Instead, we invent a package and put
# results into @r::r. The downside is that we have to manually
# free it up at the end by calling undef().
# The $formatflag argument is used for cases like the livebench
# server for whom we have to do a bit of file re-writing.
sub do_lib (\@ \@ $ $ $ $) {
    my ( $structlist, $struct_dirs, $query_struct, $coord1, $title,
        $formatflag ) = @_;
    my (@pair_sets);
    use File::Basename;
    my ($coord2, $pvec1, $pvec2);
    my $vecname = basename($query_struct, ".pdb");
    print "Basename is $vecname \n";
    my $vecfile = "$pvecdir/$vecname.vec";
    if ( -e $vecfile){
        print "vecfile exists \n";
        $pvec1 = prob_vec_read($vecfile);
    }
    else {
        print "$vecfile not found -recalculating.\n";
        my $gauss_err = 0.4;
        my $classfile = '/Users/tom/classfile';
        my $classfcn = get_clssfcn($classfile, $gauss_err);
        #my $filename = "/bm/pdb90_bin/$query_struct.bin";
	    #my $workdir = "/bm/pdb90_bin/";
	    #my $coord = coord_read($filename);
	    $pvec1 = strct_2_prob_vec($coord1, $classfcn);
	    #$filename =~ s/\/bm\/pdb90_bin\///g;
	    #$vecname = $filename;
	    #$vecname =~ s/\.bin/\.vec/g;
	#if(!prob_vec_write($pvec, "pvecs/$vecname")){
	#    print "FAILED ";
	#}
	#else{
	#    print "SUCCESS ";
	#}
    }
    my $minsize = seq_size( coord_get_seq($coord1) );
    for ( my $i = 0 ; $i < @$structlist ; $i++ ) {
	$progress->SetValue($progress->GetValue()+1);
	$progress->Refresh();
	if ( -e "$pvecdir/$$structlist[$i].vec"){
            $pvec2 = prob_vec_read("$pvecdir/$$structlist[$i].vec");
        }
        else {
            print "$pvecdir/$$structlist[$i].vec not found\n";
        }
        $coord2 =
          coord_read(
            get_path( @$struct_dirs, $$structlist[$i] . $bin_suffix ) );
        $r::r[$i] = do_align( $coord2, $coord1, $pvec2, $pvec1 );
        
        #if ( $crap_a != 0 ) {
        #}
        #else {
        #    $r::r[$i][1] = -1;
        #    print STDERR "PERL: error in coord_rmsd";
        #}
        $r::r[$i][4] =
          get_seq_id_simple( $r::r[$i][7], coord_get_seq($coord2),
            coord_get_seq($coord1) );
        $r::r[$i][6] = $r::r[$i][3] * seq_size( coord_get_seq($coord2) );
        if ( $r::r[$i][6] < 25 ) {
            $r::r[$i][8] = 0;
        }
#        $minsize = seq_size( coord_get_seq($coord1) );
#        if ( seq_size( coord_get_seq($coord2) ) < $minsize ) {
#            $minsize = seq_size( coord_get_seq($coord2) );
#        }
        $r::r[$i][2] = $r::r[$i][9] * ( ( $r::r[$i][6] ) / $minsize );
    }

    my @indices;
    for (my $i = 0 ; $i < @$structlist ; $i++ ) {
        $indices[$i] = $i;
    }

    @indices = sort { $r::r[$b][2] <=> $r::r[$a][2]; } @indices;

    print "\n",
"sw refers to Smith and Waterman alignment, nw refers to Needleman & Wunsch\n",
      "z scr : the z-score of the alignment with 1000 alternative alignments\n",
"sw scr: the combined results of score function + gaps for sw alignment\n",
      "sw cvr: coverage (fraction of sequence) accounted for by sw alignment\n",
      "nw cvr: coverage                        accounted for by nw alignment\n",
      "asize : length of the alignment\n",
      "rmsd  : root-mean-square-distance between aligned residues\n",
"____________ Summary of best templates   _________________________________\n";
    if ( -d $modeldir ) {
        print "Directory $modeldir exists. Adding new models.\n";
    }
    else {
        if ( !mkdir( "$modeldir", 0777 ) ) {
            bad_exit("Fail create modeldir ($modeldir): $!");
        }
    }
    my $p1 = coord_name($coord1);

    my $todo = ( @$structlist > $N_BRIEF_LIST ? $N_BRIEF_LIST : @$structlist );
    printf "%8s %8s %8s %8s %6s %6s %8s %6s %6s\n", 'struct', 'z scr', 'f_dme',
      'sw scr', 'sw cvr', 'seq id', 'asize', 'rmsd', 'andr_scr';
    my $j = 0;
    for (my $i = 0 ; $i < $todo ; $i++ ) {
        my $idx      = $indices[$i];
        my $pair_set = $r::r[$idx][7];
        my $a       = $r::r[$idx];
        my $coord2   =
          coord_read(
            get_path( @$struct_dirs, $$structlist[$idx] . $bin_suffix ) );
        ( $$a[1], $$a[10], my $crap_b ) =
          coord_rmsd( $$a[7], $coord2, $coord1, 0 );
        my $q =
          ( $$a[6] * $$a[6] ) / ( ( 1 + ( $$a[1] / 3 ) * ( $$a[1] / 3 ) ) *
              seq_size( coord_get_seq($coord1) ) *
              seq_size( coord_get_seq($coord2) ) );
        my $p2 = $$structlist[$idx];
        my $sid;
        if ( $$a[6] != 0 ) {
            $sid = $$a[4] / $$a[6];
        }
        else {
            $sid = 0.0;
        }
        printf "%8s %8.3g %8.3g %8.1f %6.2f %6.2g %8.4g %6.2f %4.4f\n",
          $$structlist[$idx], $$a[8], $$a[9], $$a[0], $$a[6]/seq_size(coord_get_seq($coord1)), $sid, $$a[6],
          $$a[1], $$a[2];
    	$$a[3]=$$a[6]/seq_size(coord_get_seq($coord1));
	$$a[4]=$sid;
	$$a[5]=$$structlist[$idx];
	#$$a[10]=$coord2;
	#@structures[$idx]=$coord2;
	$r::r[$idx]=$a;
    }
    print "\n";
	for (my $i = 0 ; $i < @$structlist ; $i++ ) {
		my $bucket = $r::r[$i];
		$r::r[$i] = $r::r[$indices[$i]];
		$r::r[$indices[$i]] = $bucket;
	}
    
    for (my $i = 0 ; $i < $N_ALIGNMENTS && $i < $todo ; $i++ ) {
        my $idx = $indices[$i];
        my $a = $r::r[$idx];
        my $coord2   =  coord_read(
            get_path( @$struct_dirs, $$structlist[$idx] . $bin_suffix ) );
        print(pair_set_pretty_string($$a[7], coord_get_seq($coord2), coord_get_seq($coord1)));
    }
    if($N_MODELS > 0){
	print "writing models to $modeldir";
        coord_2_pdb("$query_struct.pdb", $coord1);
        for (my $i = 0 ; $i < $N_MODELS && $i < $todo ; $i++ ) {
            my $idx = $indices[$i];
            my $a = $r::r[$idx];
            coord_2_pdb("$modeldir/$$structlist[$idx].pdb", $$a[10]);
        }
    }
    return ($r);
}

# ----------------------- wxPerl GUI -----------------------------------

#--------------------------PickFile----------------------------------
sub PickFile {
    my( $self, $event ) = @_;
	my $filedialog = Wx::FileDialog->new( $self, 'Select query structure', '', '', '*.pdb', 0, [-1, -1]);
	if($filedialog->ShowModal() == 5100){
	$filebox->SetValue($filedialog->GetPath());
	}
}

#------------------------ShowAlign------------------------------------------
sub ShowAlign{
    my ($idx, $col) = @_;
    my $rl = $r::r[$idx];
    $alignment->Clear();
    #print @structures[$idx]." \n";
    $alignment->WriteText(pair_set_pretty_string($$rl[7], 
	coord_get_seq($$rl[10]), coord_get_seq($query_struct)));
	coord_2_pdb("/tmp/result.pdb", $$rl[10]); 
	ShowStruct($q_path, "/tmp/result.pdb");
}

#-------------------------DoSearch--------------------------------------
sub DoSearch {
	my( $self, $event ) = @_;
	print "Searchin now... \n";

    my ($structfile, $libfile);
    my $fatalflag = undef;
    @struct_dirs = @DFLT_STRUCT_DIRS;

    set_params ();
    #my $q_path = $filebox->GetLineText(0);
    #$query_struct   = pdb_read($q_path, '', '');
    $libfile   = $ARGV[0];

    check_dirs (@struct_dirs);
    if (@struct_dirs == 0) {
        die "\nNo valid structure directory. Stopping.\n"; }

    (@struct_list  = get_prot_list($libfile)) || $fatalflag++;
    check_files (@struct_dirs, @struct_list, $bin_suffix) && $fatalflag++;
    if ($fatalflag) {
        print STDERR "struct dirs were @struct_dirs\n"; }


    if ( $fatalflag) {
        print STDERR "Fatal problems\n";
        return EXIT_FAILURE;
    }

    print
        "Library from $libfile with ", $#struct_list + 1,
        " template structures\n",
        "Brief results printed for $N_BRIEF_LIST templates\n",
        "Long alignments printed for $N_ALIGNMENTS templates\n",
        "Models made for the best $N_MODELS models and put in $modeldir\n";
    
#    print "Sequence read from \'$structfile\' with length ",
#          seq_size (coord_get_seq ($query_struct)),
#          " residues\n";

    my $title = 'temp title test';
    my $formatflag = 'not set';
    $progress->SetRange($#struct_list);
    print "qpath is $q_path \n";
    $r = do_lib (@struct_list , @struct_dirs, $q_path,
                    $query_struct, $title, $formatflag);

    $results->InsertRows(0, 9); #TODO: use variables, not constants
    for (my $i = 0 ; $i < 10 ; $i++ ) {
	    my $a = $r::r[$i];
	    my $crap = sprintf "%6s", $$a[5];			
	    $results->SetCellValue($i, 0, $crap);
	    $results->SetReadOnly( $i,  0 );
	    $crap = sprintf "%8.3g", $$a[8];
	    $results->SetCellValue($i, 1, $crap);
	    $results->SetReadOnly( $i, 1 );
	    $crap = sprintf "%8.3g", $$a[9];
	    $results->SetCellValue($i, 2, $crap);
	    $results->SetReadOnly( $i, 2 );
	    $crap = sprintf "%8.1f", $$a[0];
	    $results->SetCellValue($i, 3, $crap);
	    $results->SetReadOnly( $i, 3 );
	    $crap = sprintf "%6.2f", $$a[3];
	    $results->SetCellValue($i, 4, $crap);
	    $results->SetReadOnly( $i, 4 );
	    $crap = sprintf "%6.2g", $$a[4];
	    $results->SetCellValue($i, 5, $crap);
	    $results->SetReadOnly( $i, 5 );
	    $crap = sprintf "%8.4g", $$a[6];
	    $results->SetCellValue($i, 6, $crap);
	    $results->SetReadOnly( $i, 6 );
	    $crap = sprintf "%6.2f", $$a[1];
	    $results->SetCellValue($i, 7, $crap);
	    $results->SetReadOnly( $i, 7 );
	    $crap = sprintf "%4.4f", $$a[2];
	    $results->SetCellValue($i, 8, $crap);
	    $results->SetReadOnly( $i, 8 );
	    $results->SetCellValue($i, 9, 0);
	    $results->SetCellEditor($i, 9, Wx::GridCellBoolEditor->new());
	    $results->SetCellRenderer( $i, 9, Wx::GridCellBoolRenderer->new() );
	    $results->SetReadOnly($i, 9, 0);
    }
        print
"__________________________________________________________________________\n",
    "Wurst gegessen at ", scalar (localtime()), "\n";
    my ($user, $system, $crap, $crap2, $host);
    ($user, $system, $crap, $crap2) = times();
    printf "I took %d:%d min user and %.0f:%.0f min sys time\n",
    $user / 60, $user % 60, $system / 60, $system % 60;
    use Sys::Hostname;
    $host = hostname() || {$host = 'no_host'};
    print "Run on $host\n";
    return ($r);
}

#------------------------new-------------------------------------
sub new {
    my $ref = shift;
    my $self = $ref->SUPER::new( undef,           # parent window
                                 -1,              # ID -1 means any
                                 'Salami GUI',    # title
                                 [-1, -1],        # default position
                                 [900, 750],      # size
                                 );
    # controls should not be placed directly inside
    # a frame, use a Wx::Panel instead
    my $panel = Wx::Panel->new( $self,            		# parent window
                                -1,               		# ID
                                );
	$filebox = Wx::TextCtrl->new( $panel, -1, "", [100, 50], [200, -1], 0);
	# create a button
    my $button = Wx::Button->new( $panel,         		# parent window
                                  -1,             		# ID
                                  'Browse',    			# label
                                  [310, 50],       		# position
                                  [-1, -1],       		# default size
                                  );
	my $searchbtn = Wx::Button->new( $panel,         	# parent window
	                              -1,             		# ID
				                'Search',    	# label
				                 [140, 80],     # position
				                 [-1, -1],      # default size
				                 );
	$progress = Wx::Gauge->new( $panel, -1, 14, [100, 130], [700, 30] );
	$progress->SetValue(0);
	$results = Wx::Grid->new( $panel, -1, [20, 180], [850, 270]);
	my @header = ('struct', 'z scr', 'f_dme', 'sw scr', 'sw cvr', 'seq id',
		'asize', 'rmsd', 'andr_scr', ' ');
	$results->CreateGrid(1, 10);
	foreach my $col (0 .. 9){
		$results->SetColLabelValue($col, @header[$col]);
	}

	$alignment = Wx::TextCtrl->new( $panel, -1, "", [20, 470], [850, 200], 100);

	$alignment->SetFont(Wx::Font->new(12,wxMODERN,wxNORMAL,wxNORMAL,0));

	
	
    # register the OnClick method as an handler for the
    # 'button clicked' event. The first argument is a Wx::EvtHandler
    # that receives the event
    EVT_BUTTON( $self, $button, \&PickFile );
    EVT_BUTTON( $self, $searchbtn, sub{
		    	$q_path = $filebox->GetLineText(0);
		    	$query_struct   = pdb_read($q_path, '', '');
		        DoSearch();
				#$libthread = threads->new(\&DoSearch);
			#$libthread->detach();
			#$libthread = async{ DoSearch() };
    		});
    EVT_GRID_SELECT_CELL( $self, sub{
			ShowAlign($_[1]->GetRow,  $_[1]->GetCol);
			$_[1]->Skip;
    		}
	);
    return $self;
}

use OpenGL::Simple qw(:all);
use OpenGL::Simple::GLUT qw(:all);
use OpenGL::Simple::Viewer;
use Chemistry::File::PDB;
use Chemistry::Bond::Find ':all';

#-----------------------ShowStruct-----------------------------------
# Adapted from an original script by Simon Cozens.
# This opens an openGL window and displays the query 
# structure 
# also showing the superimposed result is a bit tricky without tobi's book
# I probably have to split this function up into an InitViewer()
# function and the ShowStruct() function below. Having a seperate 
# ClearViewer function prolly makes sense too. 


use strict;
use warnings;
my $mol;
my $mass_scale = 10;
my @ballpoints = ();
my @ballsticks = ();
my %ccache;
my @colours; 
my %element_colours;
my $sphericity = 20;
my %colour;
my @displaylists;
my $v=0;
my $qflag = 0;

sub pdb2displaylist{
	my ($structurefile, $stickcolor) =@_;
	$mol = Chemistry::MacroMol->read($structurefile);
	find_bonds($mol);
	my $dl = build_displaylist($mol, $stickcolor);
	return $dl;
}

sub InitViewer{
	my ($structA, $structB) = @_;

	%colour = (
    	red     => [ 1,   0,   0,   1 ],
    	yellow  => [ 1,   1,   0,   1 ],
    	orange  => [ 1,   0.5, 0,   1 ],
    	green   => [ 0,   1,   0,   1 ],
    	cyan    => [ 0,   1,   1,   1 ],
    	blue    => [ 0,   0,   1,   1 ],
    	magenta => [ 1,   0,   1,   1 ],
    	grey    => [ 0.5, 0.5, 0.5, 1 ],
    	white   => [ 1,   1,   1,   1 ],
		black   => [ 0,   0,   0,   1 ],
	);

	@colours = values %colour;
	my $iter = 0;

	%element_colours = (
	    C => $colour{grey},
	    O => $colour{red},
	    N => $colour{blue},
	    H => $colour{white},
	    S => $colour{yellow},
	    P => $colour{orange},

	);


	$|=1;
#	if($v == 0){	
		glutInit;
		$v  = new OpenGL::Simple::Viewer(
    	    title => 'Salami Structure Viewer',
    	    draw_geometry => sub { 
#    	            if (!$displaylist) { build_displaylist(); }
 					foreach my $dl (@displaylists){
 						glCallList($dl);
					}
    	    },
			zoomscale => 0.5,
    	    screenx => 512, screeny => 512,
		);
		
		glClearColor(0,0,0,1);
		glEnable(GL_LINE_SMOOTH);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
#	}
#	else {
#		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
#		glDeleteLists(1, 1);
#		glutPostRedisplay();
#	}
	return 0;
}

sub ShowStruct{
#	my ($structA, $structB) = @_;
	$qflag = 0;
	my $initflag = 1;
	if ($v == 0) {
		$initflag = 0;
		InitViewer();
	}
	@displaylists = ();
	foreach my $struct (@_){
		if ($qflag == 0){
			push (@displaylists, pdb2displaylist($struct, "white"));
			$qflag = 1;
		}
		else{
			push (@displaylists, pdb2displaylist($struct, "orange"));
		}
	}
	
	if ($initflag == 0) {glutMainLoop();}
	else {glutPostRedisplay();}

}

sub recenter {
    my @center;
    my @atoms = $mol->atoms;
    for (@atoms) {
        my @coords = $_->coords->array;
        $center[$_] += $coords[$_] for 0..2;
    }
    $center[$_] /= @atoms for 0..2;
    my $center_v = Math::VectorReal->new(@center);
    for (@atoms) { $_->coords($_->coords - $center_v) }
}

sub make_3D_model {
    my ($molecule) = @_;
    my $i = 0;
    recenter();
	@ballpoints = ();
	@ballsticks = ();
    my @atoms = $molecule->atoms;
    for my $atom ( @atoms ) {
		if ( $atom->name eq "CA" | $atom->name eq "C" | $atom->name eq "N"){
        	my $mass   = log( 1 + $atom->mass ) / $mass_scale;
        	my $color  = $element_colours{ $atom->symbol } || $colour{cyan};

        	my @coords = $atom->coords->array;
        	push @ballpoints, [ $color, $mass, @coords ];
    	}
	}
    for my $bond ( $mol->bonds ) {
        my ( $from, $to ) = $bond->atoms;
		if ( $from->name eq "CA" | $from->name eq "C" | $from->name eq "N"){
			if ( $to->name eq "CA" | $to->name eq "C" | $to->name eq "N"){
        		my @from = $from->coords->array;
        		my @to   = $to->coords->array;
        		push @ballsticks, [ \@from, \@to ];
    		}
		}
	}

}

sub visualize {
    my ($molecule, $stickcol) = @_;
#    if ( !@ballpoints ) { make_3D_model() }
	make_3D_model($molecule);
    for (@ballpoints) {
        my ( $color, $mass, @coords ) = @$_;
        glColor(@$color);
        glPushMatrix;
        glTranslate(@coords);
        glutSolidSphere( $mass , $sphericity, $sphericity);
        glPopMatrix;
    }

    glColor( @{ $colour{$stickcol} } );
    glLineWidth(4);
    for (@ballsticks) { glBegin(GL_LINES); glVertex(@$_) for @$_; glEnd; }
}

sub build_displaylist {
    my ($molecule, $stickcol) = @_;
	my $displaylist = glGenLists(1);
    glNewList($displaylist,GL_COMPILE);
    visualize($molecule, $stickcol);
    glEndList();
	return $displaylist
}

# viewer stuff ends here.
# now we're back to normal wxWidgets stuff.



package Guisearch;
use base qw(Wx::App);

#-----------------------OnInit--------------------------------------
# WxPerl stuff.
sub OnInit {
	# create a new frame (a frame is a top level window)
		my $frame = MyFrame->new;
		$frame->Show( 1 );
}


package main;

# create the application object, this will call OnInit
my $app = Guisearch->new;
# process GUI events from the application this function will not
# return until the last frame is closed
$app->MainLoop;
