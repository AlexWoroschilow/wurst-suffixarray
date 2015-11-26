# 18 Oct 2001
# Test alignments, up and down, left and right and every reasonable variation

eval 'use warnings';
if ($@) {
    print "Warnings not present\n";}


use FindBin;
use lib "$FindBin::Bin/../src/Wurst/blib/arch";
use lib "$FindBin::Bin/../src/Wurst/blib/lib";

use strict;
use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);

use Wurst;

my $paths_place = 'paths.inc';
do $paths_place || die "Failed to read $paths_place: $!";
if ($@) { die "Compilation fail on $paths_place:\n$@" };


my $gap_open1  = 40;
my $gap_widen1 = 11;
my $gap_open2  = $gap_open1 / 10;
my $gap_widen2 = $gap_widen1 / 10;
use vars qw(@strings);
*strings = [ 'cccccc',
             'a',
             'aaaggg',
             'gggaaa',
             'cwcwcwwcwwwcwwwwc'];


# ----------------------- string_to_seq  ----------------------------
# Write the string given as an argument to a temporary file,
# read it up and return the sequence.
sub string_to_seq ( $ )
{
    my $string = shift;
    my $tmpfile = POSIX::tmpnam();
    open (T, ">$tmpfile") || die "Open fail on $tmpfile";
    print T  $string;
    close (T) || die "close fail";
    my $t = seq_read ($tmpfile) || die "Fail on $tmpfile";
    unlink ($tmpfile) || die "Failed to unlink $tmpfile";
    return $t;
}


# ----------------------- do_pair      ------------------------------
sub do_pair ($ $ $ $)
{
    my ($s1, $s2, $al_type, $subst_mat) = @_;
    print "s1 gaps $gap_open1, $gap_widen1. s1 is:\n", seq_print ($s1), "\n\n";
    print "s2 gaps $gap_open2, $gap_widen2, s2 is:\n", seq_print ($s2), "\n\n";

    if ($al_type == $N_AND_W) {
        print "Alignment Needleman and Wunsch\n" }
    else {
        print "Alignment by Smith and Waterman\n" }
    my $scr_mat = score_mat_new (seq_size ($s1), seq_size ($s2)) ||
        die "score_mat_new fail: $!";

    score_smat ($scr_mat, $s1, $s2, $subst_mat);
    $scr_mat = score_mat_shift ($scr_mat, 3);
    my ($min, $max, $av, $std_dev) = score_mat_info ($scr_mat);
    print score_mat_string ($scr_mat, $s1, $s2);
    print "min max av std_dev: $min, $max, $av, $std_dev\n";
    my $result_mat;
    my $pair_set =
        score_mat_sum_smpl ($result_mat, $scr_mat,
                            $gap_open1, $gap_widen1,
                            $gap_open2, $gap_widen2,
                            $al_type);
    print pair_set_pretty_string ($pair_set, $s1, $s2);
    print "Short extension\n";
    pair_set_extend ($pair_set, seq_size ($s1), seq_size ($s2), $EXT_SHORT);
    print pair_set_pretty_string ($pair_set, $s1, $s2);
    print "Long extension\n";
    pair_set_extend ($pair_set, seq_size ($s1), seq_size ($s2), $EXT_LONG);
    print pair_set_pretty_string ($pair_set, $s1, $s2);
}

# ----------------------- test_seq     ------------------------------
sub test_seq ($ $ $)
{
    my ($string1, $string2, $subst_mat) = @_;
    my $s1 = string_to_seq ($string1);
    my $s2 = string_to_seq ($string2);
    foreach my $al_type ($N_AND_W, $S_AND_W) {
        do_pair ($s1, $s2, $al_type, $subst_mat);
        do_pair ($s2, $s1, $al_type, $subst_mat);
    }
}

# ----------------------- mymain       ------------------------------
sub mymain ()
{
    use vars qw ($MATRIX_DIR $MATRIX_FILE);
    $SIG{TRAP} = 'IGNORE';# kill 'TRAP', $$;

    my $matname = $MATRIX_DIR . '/' . $MATRIX_FILE;
    my $subst_mat = sub_mat_read ($matname) ||
        die "Fail getting sub mat $matname";
    my $num_seq = $#strings + 1;
    kill 'TRAP', $$;
    for (my $i = 0; $i < $num_seq; $i++) {
        for (my $j = $i + 1; $j < $num_seq; $j++) {
            test_seq ($strings[$i], $strings[$j], $subst_mat);
        }
    }

    return (EXIT_SUCCESS);
}
exit mymain();
