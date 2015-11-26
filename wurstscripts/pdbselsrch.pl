#!/usr/bin/perl
##$ -clear
#$ -w e 
##$ -l arch=glinux -l short=0
#$ -p -50
#$ -S /home/torda/bin/perl
#$ -cwd
#$ -j y
#$ -m e -M torda@zbh.uni-hamburg.de
#$ -q hpc.q
# 18 March 2002
# Most of the other scripts here are for testing with
# deliberately crazy values.
# This one uses the installed, working binary and serious values
# for all parameters.
# rcsid = $Id: libsrch.pl,v 1.18 2005/10/28 11:21:40 torda Exp $

=pod

=head1 NAME

libsrch.pl - Given a structure, align it to a library of templates

=head1 SYNOPSIS

libsrch.pl [options] I<struct_file> I<struct_lib_list> I< S<[ phd_file ]> >


=head1 DESCRIPTION

Given a structure, align it to every member of a library of
templates.  The sequence is given by I<struct_file>.  The library is
a list of protein names listed in <struct_lib_list>. The last
argument is optional and contains the name of a file with
secondary structure predictions. If it is not present, the script
will look in the directory containing I<struct_file> and strip off
anything that looks like a file extension. Then, it will append
B<.phd> and try to open that (so /boo/bar/1abc.seq gives
/boo/bar/1abc.phd).

If you want to run without secondary structure, it is not enough
to omit the filename. Instead, use the B<-s> option described
below.

=head2 FILE FORMAT

The list of files which make up the template library is in a
simple format. The script will try to read anything that looks

like a four-letter protein name + chain id from the first
column. Leading white space is ignored. A valid form would look
like

   1abc_
   2qrsB
   1xyz  This text after first column is ignored

=head2 Changing library and templates.

Typically, a first run will be made with whatever library we are
using. However, one will often want to add extra .bin files
for a particular sequence. To do that,

=over

=item *

Add the new file names to the list of proteins and give it a name
like F<mylist>.

=item *

Make a directory with a name like I<templates> and put the extra
F<.bin> files in there.

=item *

Run the script with the B<-t> option like:

  perl libsrch.pl -t templates blahblah.seq mylist

=back

=head2 OPTIONS

=over

=item B<-a> I<N>

Print out details of the best I<N> alignments.

=item B<-d> I<modeldir>

Save final models in I<modeldir> instead of the default
directory, B<modeldir>.

=item B<-h> I<N>

After doing alignments, a default number (maybe 50) will be
printed out. Alternatively, you can ask for the I<N> best scoring
alignments to be printed out.

=item B<-m> I<N>

After alignments, I<N> models will be built, otherwise, a small
default number will appear. Set I<N> to zero if you do not want
any models.

=item B<-s>

Do not use secondary structure predictions. This will cause a
different set of parameters to be used in the calculation.

=item B<-t> I<dir1[,dir2,...]>

Add I<dir1> to the list of places to look for template
files. This is a comma separated list, so you can add more
directories.

=back

=head1 OUTPUT

In all output

=over

=item *

B<SW> refers to the result of the second Smith and Waterman.

=item *

B<NW> refers to the result of the Needleman and Wunsch.

=item *

B<cvr> or B<cover> refers to "coverage".  This is the fraction of
sequence residues which are aligned to part of a
structure. Continuing, B<S<sw cvr>> refers to coverage from the
Smith and Waterman calculation.

=item *

The script prints out the coverage in a somewhat pictorial form
which might look like

   ----XXXXX---XXX

where the X's mean a residue was aligned.

=back

=head1 MODELS

Models will be written in PDB format for the best few
sequences. They will get written to a directory called
F<modeldir>. Maybe this should be made an option.

=head1 OPERATION

Currently the script does

=over

=item Smith and Waterman step 1

This is a locally optimal alignment.

=item Smith and Waterman step 2

This is another locally optimal alignment, but forced to pass
through the same path as the first one. It provides a small
extension to the alignment.

=item Needleman and Wunsch

This is a globally optimal alignment, but forced to pass through
the preceding Smith and Waterman.

=back

=head1 QUESTIONS and CHANGES

=item *

The selection of which scores to print out is a bit arbitrary.

=item *

the coverage picture is very ugly. It could be
beautified.

=item *

The coverage picture corresponds to the Smith and
Waterman. Perhaps it should be the Needleman and
Wunsch. Obviously, both are possible, but just a bit ugly.

=cut

#use lib "$ENV{HOME}/pl/lib/i586-linux-thread-multi";  # Where wurst lives after installation
#use lib "/home/stud2004/tmargraf/pl/lib/i686-linux-thread-multi";

use FindBin;
use lib "$FindBin::Bin/../../wurst/src/Wurst/blib/arch";
use lib "$FindBin::Bin/../../wurst/src/Wurst/blib/lib";


use Wurst;

use vars qw ($MATRIX_DIR $PARAM_DIR
             $RS_PARAM_FILE $FX9_PARAM_FILE );

do "$ENV{HOME}/../../torda/c/wurst/scripts/paths.inc" || die $@;

if ($@) {
    die "broke reading paths.inc:\n$@"; }
if ( defined ($ENV{SGE_ROOT})) {
    $MATRIX_DIR = "$ENV{HOME}/../../torda/c/wurst/matrix";
    $PARAM_DIR  = "$ENV{HOME}/../../torda/c/wurst/params";
}

use strict;

use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);


# ----------------------- Defaults  ---------------------------------
# These are numbers you might reasonably want to change.
# They should (will) be changeable by options.
use vars qw ($N_BRIEF_LIST $N_MODELS $N_ALIGNMENTS);
$N_BRIEF_LIST   = 500;
$N_MODELS       = 100;
$N_ALIGNMENTS   = 10;
use vars qw ($modeldir $DFLT_MODELDIR);
*DFLT_MODELDIR = \ 'modeldir';
$modeldir = $DFLT_MODELDIR;

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
*DFLT_STRUCT_DIRS = ['/home/stud2004/tmargraf/pdb90_bin.050506' , '/bm/pdb90_bin' , '.'];
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
    $coverage = $coverage / seq_size (coord_get_seq($coord1)); #and as fraction of query structure

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
    my $model = make_model ($pair_set, coord_get_seq ($c2), $c1);
    if ($model == 0 || coord_size($model) < 10) {
	return 0.0;}
    my $frac;
    if ( ! dme_thresh ($frac, $model, $c1, 4.0)){
	print STDERR "dme_thresh broke.\n"; }
    
    return ($frac);
}

# ----------------------- do_align ----------------------------------
# This does the alignment. Although it takes secondary structure
# data as an argument, we are not yet using this in our server
# script.
sub do_align ($ $ $ $)
{
    my($coord1, $coord2, $pvec1, $pvec2) = @_;
    my $seq_ptr1 = coord_get_seq($coord1);
    my $seq_ptr2 = coord_get_seq($coord2);
	    
    # build score matrix
    my $matrix = score_mat_new(seq_size($seq_ptr1), seq_size($seq_ptr2));

    score_pvec($matrix, $pvec1, $pvec2);
    $matrix = score_mat_shift ($matrix, $m_shift);

    my ($sw_scr_tot, $sw_coverage, $sw_score1, $sw_score1_gap,
        $sw_seq_gap, $sw_strct_gap);
    my ($nw_scr_tot, $nw_coverage, $nw_score1, $nw_score1_gap,
        $nw_seq_gap, $nw_strct_gap);
    
    my $sw_pair_set =
        score_mat_sum_smpl (my $crap_mat, $matrix,
                            $sw1_pgap_open, $sw1_pgap_widen,
                            $sw1_qgap_open, $sw1_qgap_widen,
                            $S_AND_W);


	($sw_scr_tot, $sw_coverage, $sw_score1, $sw_score1_gap, $sw_seq_gap, $sw_strct_gap) =
	    get_scores ($sw_pair_set, $coord1, $coord2, 's_and_w');
    my $frac_dme = get_dme_thresh ($sw_pair_set, $coord1, $coord2);

    my $num_scrs = 1000;
    my $alt_scrs_ref =
        get_alt_scores($num_scrs, $matrix, $sw_pair_set);

    my ($mean, $deviation) = normalize_alt_scores($alt_scrs_ref);
    
    undef ($alt_scrs_ref);
    my $z_scr;
    if ($deviation != 0) {
        $z_scr = ($sw_scr_tot - $mean) / $deviation; }
    else {
        $z_scr = 0.0; }                          # Should not really happen
#   If the alignment is tiny, one can get a ridiculous z-score
    if ($sw_coverage < 0.03) {                   # silently wipe these guys out
        $z_scr = 0.0; }


    my $sec_strct_pnlty = 0.0;
    my $nw_pair_set =
        score_mat_sum_sec (my $result_mat2, $matrix,
                           $coord2, $sec_strct_pnlty,
                           $nw_pgap_open, $nw_pgap_widen,
                           $nw_qgap_open, $nw_qgap_widen,
                           $N_AND_W, $sw_pair_set);
    ($nw_scr_tot, $nw_coverage, $nw_score1, $nw_score1_gap,
     $nw_seq_gap, $nw_strct_gap) =
    get_scores ($nw_pair_set, $coord1, $coord2, 'n_and_w');
    my @r =
        ($sw_scr_tot, $nw_pair_set, $nw_scr_tot,
         $sw_coverage, $nw_coverage, $sw_score1,
         $nw_score1, $sw_pair_set, $z_scr, $frac_dme);
    return (\@r);
}

# ----------------------- zero_shift_mat ----------------------------
# This shifts a substitution matrix, not an alignment score matrix.
sub zero_shift_mat ($ $)
{
    my ($sub_mat, $shift) = @_;
    for (my $i = 0; $i < 20; $i++) {
        for (my $j = $i; $j < 20 ; $j++) {
            my $t = sub_mat_get_by_i ($sub_mat, $i, $j) + $shift;
            sub_mat_set_by_i ($sub_mat, $i, $j, $t);
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
sub do_lib (\@ \@ $ $ $)
{
    my ($structlist, $struct_dirs, $query_struct, $title, $formatflag) = @_;
    my (@pair_sets);
    
    my ($coord1, $coord2);
    $coord1 = coord_read($query_struct);
    my $vecname = $query_struct;
    $vecname =~ s/\.bin/\.vec/g;
    $vecname =~ s/\/bm\/pdb90_bin\///g;
    $vecname = "pvecs/$vecname";
    my $pvec1 = prob_vec_read($vecname);
    for (my $i = 0; $i < @$structlist ; $i++) {
	my $pvec2 = prob_vec_read("pvecs/$$structlist[$i].vec");	
        $coord2   = coord_read(get_path (@$struct_dirs, $$structlist[$i] . $bin_suffix));
	$r::r[$i] = do_align ($coord2, $coord1, $pvec2, $pvec1);
	($r::r[$i][1], my $crap_a, my $crap_b) = coord_rmsd($r::r[$i][7], $coord2, $coord1, 0);
	if($crap_a != 0){
        }
	else{
	    $r::r[$i][1] = -1;
	    print STDERR "PERL: error in coord_rmsd";
	}
	$r::r[$i][4]=get_seq_id_simple($r::r[$i][7], coord_get_seq($coord2), coord_get_seq($coord1));
	$r::r[$i][6]=$r::r[$i][3]*seq_size(coord_get_seq ($coord2));
	if($r::r[$i][6] < 25) {
		$r::r[$i][8] = 0;
	}
    }

    my @indices;
    for ( my $i = 0; $i < @$structlist; $i++) {
        $indices[$i] = $i; }


    @indices = sort {
        $r::r[$b][8] <=> $r::r[$a][8];
    } @indices;

    print "\n",
"sw refers to Smith and Waterman alignment, nw refers to Needleman & Wunsch\n",
"z scr : the z-score of the alignment with 1000 alternative alignments\n",
"sw scr: the combined results of score function + gaps for sw alignment\n",
"sw cvr: coverage (fraction of sequence) accounted for by sw alignment\n",
"nw cvr: coverage                        accounted for by nw alignment\n",
"asize : length of the alignment\n",
"rmsd  : root-mean-square-distance between aligned residues\n",
"____________ Summary of best templates   _________________________________\n";
    if ( -d $modeldir) { 
        print "Directory $modeldir exists. Adding new models.\n"; }
    else {
        if ( ! mkdir("$modeldir",0777)) {
	    bad_exit( "Fail create modeldir ($modeldir): $!"); } }
    my $p1 = coord_name ($coord1);

    my $todo = (@$structlist > $N_BRIEF_LIST ? $N_BRIEF_LIST : @$structlist);
    printf "%8s %8s %8s %8s %6s %6s %8s %6s %6s\n",
    'struct',  'z scr', 'f_dme', 'sw scr', 'sw cvr', 'seq id', 'asize', 'rmsd', 'qval';
    my $j = 0;
    for (my $i = 0; $i < $todo; $i++) {
        my $idx = $indices [$i];
	my $pair_set = $r::r[$idx][7];
        my $a = $r::r[$idx];
        my $coord2 = coord_read (
                  get_path (@$struct_dirs, $$structlist[$idx] . $bin_suffix));
	my $q = ($$a[6]* $$a[6]) / ((1+($$a[1]/3)*($$a[1]/3))*seq_size(coord_get_seq ($coord1))*seq_size(coord_get_seq ($coord2))); 	
	my $p2 = $$structlist[$idx];
	my $sid;
	if ( $$a[6] != 0 ){
		$sid = $$a[4]/$$a[6];
	}
	else {
		$sid = 0.0;
	}
        printf "%8s %8.3g %8.3g %8.1f %6.2f %6.2g %8.4g %6.2f %4.4f\n",
               $$structlist[$idx],
               $$a[8], $$a[9], $$a[0], $$a[3], $sid, 
               $$a[6], $$a[1], $q;
    }
    print "\n";



    undef (@r::r);
    return EXIT_SUCCESS;
}



# ----------------------- mymain  -----------------------------------
# Arg 1 is a structure file. Arg 2 is a structure list file.
sub mymain ()
{
    use Getopt::Std;
    my (%opts);
    my (@struct_list, @struct_dirs);
    my ($structfile, $libfile);
    my $fatalflag = undef;
    @struct_dirs = @DFLT_STRUCT_DIRS;

    if ( !getopts ('a:d:h:m:t:', \%opts)) {
        usage(); }
    if ( defined ( $opts { a })) {$N_ALIGNMENTS = $opts { a }}
    if ( defined ( $opts { d })) {$modeldir     = $opts { d }}
    if ( defined ( $opts { h })) {$N_BRIEF_LIST = $opts { h }}
    if ( defined ( $opts { m })) {$N_MODELS     = $opts { m }}
    if ( defined ( $opts { t })) {push (@struct_dirs,split( ',', $opts { t }))}
    undef %opts;

    set_params ();

    if ( $#ARGV < 0) {
        print STDERR "Must have at least a query structure file\n";
        usage(); }
    if ( $#ARGV < 1) {
        print STDERR "Please give me a structure library / file\n";
        usage();
    }

    my $query_struct   = coord_read($ARGV[0]);
    $libfile   = $ARGV[1];

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
    
    print "Sequence read from \'$structfile\' with length ",
          seq_size (coord_get_seq ($query_struct)),
          " residues\n";

    my $title = 'temp title test';
    my $formatflag = 'not set';

    my $r = do_lib (@struct_list , @struct_dirs, $ARGV[0],
                    $title, $formatflag);
    if ($r == EXIT_FAILURE) {
        bad_exit ('calculation broke') }

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
    return EXIT_SUCCESS;
}
# ----------------------- main    -----------------------------------
exit (mymain());
