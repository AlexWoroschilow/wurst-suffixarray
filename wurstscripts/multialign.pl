#!/usr/bin/perl
#use lib "$ENV{HOME}/pl/lib/i586-linux-thread-multi";  # Where wurst lives after installation
#use lib "/home/stud2004/tmargraf/pl/lib/i686-linux-thread-multi";

use FindBin;
use lib "$FindBin::Bin/../../wurst/src/Wurst/blib/arch";
use lib "$FindBin::Bin/../../wurst/src/Wurst/blib/lib";

use Wurst;
use GraphViz;

use vars qw ($MATRIX_DIR $PARAM_DIR
  $RS_PARAM_FILE $FX9_PARAM_FILE );

#use strict;

use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);

# ----------------------- Defaults  ---------------------------------
# These are numbers you might reasonably want to change.
# They should (will) be changeable by options.
use vars qw ($N_BRIEF_LIST $N_MODELS $N_ALIGNMENTS $DFLT_MAX_ATTACH);
$N_BRIEF_LIST = 100;
$N_MODELS     = 5;
$N_ALIGNMENTS = 0;
$DFLT_MAX_ATTACH = 5;
use vars qw ($modeldir $DFLT_MODELDIR $pvecdir);
*DFLT_MODELDIR = \'modeldir';
$modeldir      = $DFLT_MODELDIR;
$pvecdir       = '/home/stud2004/tmargraf/andrew/scripts/pvecs';

my %structs;

# This is a top level temporary directory. We can have a cron job
# run around and clean it up at night. Each job creates its own
# temporary directory under this one.
use vars qw ($log_base $top_temp);
*top_temp   = \"./wurst_delete_able_temp";
# Where we will write logs to
*log_base   = \"log";


# Define our mail program, reply-to address and from address.
use vars qw ($mail_from_addr $mail_reply_to $mail_prog);
*mail_from_addr = \'"Wurst results" <nobody@zbh.uni-hamburg.de>';
*mail_reply_to  = \'nobody@zbh.uni-hamburg.de';
*mail_prog      = \'/usr/bin/mail';

# Switches..
# During testing, we do not want to be able to switch off things
# like the calculation, mailing... These are turned on and off
# here.

use vars qw ($redirect_io $really_mail
             $do_the_calculation $fake_args $verbose);
my $nothing = undef;
*redirect_io  =       \0;
*really_mail  =       \0;
*do_the_calculation = \1;
*fake_args          = \$nothing;
*verbose            = \$nothing;
undef $nothing;

# ----------------------- Global variables  -------------------------
# Unfortunately, we need to store some things here, mainly in
# case we have to quickly die. The bad_exit() routine can mail
# back something informative if it knows the address and job
# title.
use vars qw ($email_address $tmp_dir $title);

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
*DFLT_STRUCT_DIRS = ['/home/stud2004/tmargraf/pdbsnapshot060307']; #/bm/pdb90_bin
*bin_suffix       = \'.bin';

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

# ----------------------- get_prot_list -----------------------------
# Go to the given filename and get a list of proteins from it.
sub get_prot_list ($) {
    my $f = shift;
    my @a;
    if ( !open( F, "<$f" ) ) {
        print STDERR "Open fail on $f: $!\n";
        return undef;
    }

    while ( my $line = <F> ) {
        chomp($line);
        my @words = split( ' ', $line );
        if ( !defined $words[0] ) { next; }
        $line = $words[0];
        $line =~ s/#.*//;     # Toss comments away
        $line =~ s/\..*//;    # Toss filetypes away
        $line =~ s/^ +//;     # Leading and
        $line =~ s/ +$//;     # trailing spaces.
        if ( $line eq '' ) {
            next;
        }
        substr( $line, 0, 4 ) = lc( substr( $line, 0, 4 ) );    # 1AGC2 to 1agc2
        if ( length($line) == 4 ) {    # Convert 1abc to 1abc_
            $line .= '_';
        }
        push( @a, $line );
    }
    close(F);
    return (@a);
}

# ----------------------- get_path  ---------------------------------
# We have a filename and a list of directories where it could
# be. Return the path if we can find it, otherwise return undef.
sub get_path (\@ $) {
    my ( $dirs, $fname ) = @_;
    foreach my $d (@$dirs) {
        my $p = "$d/$fname";
        if ( -f $p ) {
            return $p;
        }
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
sub check_dirs (\@) {
    my $a    = shift;
    my $last = @$a;
    for ( my $i = 0 ; $i < $last ; $i++ ) {
        if ( !-d $$a[$i] ) {
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
sub check_files (\@ \@ $) {
    my ( $dirs, $fnames, $ext ) = @_;
    my $errors = 0;
    foreach my $f (@$fnames) {
        my $name = "$f$ext";
        if ( !get_path( @$dirs, $name ) ) {
            $errors++;
            print STDERR "Cannot find $name\n";
        }
    }
    return $errors;
}

# ----------------------- usage   -----------------------------------
sub usage () {
    print STDERR "Usage: \n    $0 -q query_struct_file -l struct_library \n";
    print STDERR "Optional parameters:\n    -a # of Alignments printed\n";
    print STDERR "    -s # of superimposed structures written\n";
    print STDERR "    -h # of results in the list\n";
    print STDERR "    -d directory structures are written to\n";
    exit(EXIT_FAILURE);
}

# ----------------------- get_scores --------------------------------
sub get_scores ($ $ $ $) {
    my ( $pair_set, $size1, $size2, $to_use ) = @_;

    my ( $scr_tot, $coverage, $score1, $geo_gap, $score1_gap );
    my ( $str1, $crap );
    my ( $open_cost, $widen_cost, $nseq_gap );
    ( $score1_gap, $score1 ) = pair_set_score($pair_set);
    ( $str1,       $crap )   =
      pair_set_coverage( $pair_set, $size1, $size2 );
    $coverage = ( $str1 =~ tr/1// );    # This is coverage as an integer
    $coverage = $coverage / $size1;     # and as fraction of query structure
      

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
    restore_handlers();  # otherwise an ugly loop is possible
    print STDERR "Error: \"$msg\"\n";
    exit (EXIT_FAILURE);
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

# ----------------------- do_align ----------------------------------
# This does the alignment. 

sub do_align ($ $) {
    my ( $pvec1, $pvec2 ) = @_;

    # build score matrix
    my $matrix = score_mat_new( prob_vec_length($pvec1), prob_vec_length($pvec2) );

    score_pvec( $matrix, $pvec1, $pvec2 );

    my ($min, $max, $av, $std_dev) = score_mat_info ($matrix);
    $matrix = score_mat_shift( $matrix, $m_shift );

    my (
        $sw_scr_tot,    $sw_coverage, $sw_score1,
        $sw_score1_gap, $sw_seq_gap,  $sw_strct_gap
    );

    my $sw_pair_set  = score_mat_sum_smpl(
        my $crap_mat,   $matrix,         $sw1_pgap_open, $sw1_pgap_widen,
        $sw1_qgap_open, $sw1_qgap_widen, $S_AND_W
    );
    my $nw_pair_set  = score_mat_sum_smpl(
        my $crap_mat,   $matrix,         $sw1_pgap_open, $sw1_pgap_widen,
        $sw1_qgap_open, $sw1_qgap_widen, $N_AND_W
    );
    #my $sw_pair_set  = score_mat_sum_smpl(
    #    my $crap_mat,   $matrix,         $av, $av,
    #    $av, $av, $S_AND_W
    #);

    (
        $sw_scr_tot,    $sw_coverage, $sw_score1,
        $sw_score1_gap, $sw_seq_gap,  $sw_strct_gap
      )
      = get_scores( $sw_pair_set, prob_vec_length($pvec1), prob_vec_length($pvec2), 's_and_w' );


    my (
        $nw_scr_tot,    $nw_coverage, $nw_score1,
        $nw_score1_gap, $nw_seq_gap,  $nw_strct_gap
    )
    = get_scores( $nw_pair_set, prob_vec_length($pvec1), prob_vec_length($pvec2), 'n_and_w' );

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
	my $asize = ( $nw_coverage * prob_vec_length($pvec1) );
	my $norm_score  = $nw_scr_tot / $asize; 

    my $sec_strct_pnlty = 0.0;
    my $r  = {
        'sw_score' => $sw_scr_tot,
		'norm_score' => $nw_scr_tot,#$norm_score, 
		'asize'      => $asize,
	#	'nw_pairset' => $nw_pair_set,
	#	'nw_score' => $nw_scr_tot,
	#	'sw_coverage' => $sw_coverage,
	#	'nw_coverage' => $nw_coverage,
    #    'sw_score1' => $sw_score1,
	#	'nw_score1' => $nw_score1,
		'sw_pairset' => $nw_pair_set,
	#	'z_score' => $z_scr,
		'pvec1' => $pvec1, 
		'pvec2' => $pvec2,
		#'coord1' => $coord1,
		#'coord2' => $coord2,
		#'name1' => $name1,
		#'name2' => $name2
    };
    return ( $r );
}

# ----------------------- make_nodes -----------------------------------
# creates nodes for each structure in the list.

sub make_nodes (\@) {
	my ( $structlist ) = @_;
	my $pvec1;
	my @nodes = ();
	# create a node for each structure...
	for ( my $i = 0 ; $i < @$structlist ; $i++ ) {	
		if ( -e "$pvecdir/$$structlist[$i].vec"){
			$pvec1 = prob_vec_read("$pvecdir/$$structlist[$i].vec");
		}
		else {
			print "$pvecdir/$$structlist[$i].vec not found\n";
		}
		my @members = ();
		push @members, $$structlist[$i];
		push @nodes, {
			'name' => $$structlist[$i],
			'pvec' => $pvec1,
			'members' => \@members,
		};
	}
	return(@nodes);
}

# ----------------------- build_lib  -----------------------------------
# Walk over a list, doing all-vs.-all alignments and saving interesting
# scores. The definition of interesting is a bit arbitrary.
# There is one very non-obvious coding trick.  We need to be able
# to pass the score information into the sorting functions. We
# no longer use the r::r mess. instead we do something else which I should
# document at some point...

sub build_lib (\@ \@ \@) {
    my ( $structlist, $struct_dirs, $nodes) = @_;
    my ($coord1, $coord2, $pvec1, $pvec2) ;
    my $classfile = '/home/stud2004/tmargraf/andrew/scripts/classfile';
    my $classfcn = get_clssfcn($classfile, $gauss_err);
	my $idx = 0;
	my $results;
	my @alns ;
	my @totscores  = ();
	foreach (@$structlist){
		push (@totscores, 0.0);
	}

	for ( my $i = 0 ; $i < @$structlist ; $i++ ) {	
		for ( my $j = 0 ; $j < $i ; $j++ ) {
			$results = do_align( $$nodes[$j]{"pvec"}, $$nodes[$i]{"pvec"} );
			$results->{'node1'} = $i;
			$results->{'node2'} = $j;
			push @alns, $results;
			$idx++;
			$totscores[$i] += $results->{norm_score};
			$totscores[$j] += $results->{norm_score};
		}

	}
	my $max_i =0;
	my $curr_i =0;
	my $max = $totscores[0];
	foreach (@totscores){
		if($_ > $max){
			$max = $_;
			$max_i = $curr_i;
		}
		$curr_i++;
	}
	print "centroid is $$structlist[$max_i]\n";	

	#superimposing structures...
	my $coord1 = coord_read(get_path( @$struct_dirs, $$structlist[$max_i] . $bin_suffix ) );
	$structs{"$$structlist[$max_i]"} = $coord1;
	coord_2_pdb("modeldir/$$structlist[$max_i].pdb", $coord1);
	my $rmsd = 0.0;
	my $totalsize = 0;
	foreach(@alns){
		if($$nodes[$_->{node1}]{name} eq $$structlist[$max_i]){
			my $coord2 = coord_read(get_path(@$struct_dirs, $$nodes[$_->{node2}]{name}.$bin_suffix));
			my($t1, $t2, $t3) = coord_rmsd($_->{sw_pairset}, $coord2, $coord1, 0);
			$rmsd += ($t1 * $_->{asize});
			$totalsize += $_->{asize};
			$structs{"$$nodes[$_->{node2}]{name}"} = $t2;

			coord_2_pdb("modeldir/".$$structlist[$max_i]."_".$$nodes[$_->{node2}]{name}.".pdb", $t2);
		}
	}
	foreach(@alns){
		if($$nodes[$_->{node2}]{name} eq $$structlist[$max_i]){
			my $coord2 = coord_read(get_path(@$struct_dirs, $$nodes[$_->{node1}]{name}.$bin_suffix));
			my($t1, $t2, $t3) = coord_rmsd($_->{sw_pairset}, $coord2, $coord1, 6);
			$rmsd += ($t1 * $_->{asize});
			$totalsize += $_->{asize};
			$structs{"$$nodes[$_->{node1}]{name}"} = $t2;
			coord_2_pdb("modeldir/$$structlist[$max_i]_$$nodes[$_->{node1}]{name}.pdb", $t2);
		}
	}
	$rmsd = $rmsd/$totalsize;
	$totalsize = $totalsize/(@$structlist-1);
#	print "RMSD is $rmsd over $totalsize residues\n";
	return(@alns, $$structlist[$max_i]);
}

# ----------------------- purge_list -----------------------------------
# remove all other alignments involving the two structures in the
# chosen alignment
# ---------------------------------------------------------------------
sub purge_list (\@ \@) {
		my ($list, $chosen) = @_;
		my %delete;
		my $r_size = @$list;
		for (my $i = 0; $i < $r_size; $i++){
			my $eflag =0;
			foreach (@$chosen){
				if ($$list[$i]{node1} eq $_ || $$list[$i]{node2} eq $_ ){ 
					$delete{sprintf("%6g", $i)}=1;}
			}
		}

		my $count = 0;	
		foreach( sort keys %delete ){  # (my $crap, my $val) = each(%delete)){
			#splice(@$list, $crap, 1);
			splice(@$list, $_ - $count, 1);
			$count++;
		}
		return(@$list);
}

# -------------------- upgma_check -----------------------------------
#
#
# --------------------------------------------------------------------
sub upgma_check(\@) {
	my ($dists) = @_;
	my $pass = 0;
	my $fail = 0;
	foreach (@$dists){
		my $i = $_;
		foreach (@$dists){
			my $j = $_;
			foreach (@$dists){
				my $k = $_;
				if($k{node1} eq $j{node2} && $i{node1} eq $j{node1}){
					if($$k{norm_score} eq $$j{norm_score} && $$i{norm_score} le $$k{norm_score}){
						$pass++;
					}
					elsif($$i{norm_score} eq $$j{norm_score} && $$k{norm_score} le $$j{norm_score}){
						$pass++;
					}
					elsif($$k{norm_score} eq $$i{norm_score} && $$j{norm_score} le $$i{norm_score}){
						$pass++;
					}
					else { $fail++;}
				}
			}
		}
	}
	print "Ultrametric check: pass: $pass, fail :$fail \n";
}

# ----------------------- aggregate_cons -----------------------------------
# Aggregate alignment library with consensus-NJ. 
# ---------------------------------------------------------------------
sub aggregate_cons (\@ \@ $) {
    my ( $structlist, $struct_dirs, $formatflag ) = @_;
	my @nodes = make_nodes( @$structlist );
	my @alns = build_lib(@$structlist, @$struct_dirs, @nodes);
	upgma_check(@alns);
	my $centroid = pop(@alns);
	my $r_size = @alns;
	my @sortedalns = @alns;

	my %totals;
	foreach (@alns){
		$totals{$_{node1}} += $$_{norm_score};
		$totals{$_{node2}} += $$_{norm_score};
	}

	my $graph = GraphViz->new();
	foreach (@nodes) {
		$graph->add_node($_->{'name'}, rank=>'leaf', fontsize => 12, style => 'bold');
		$_->{'graphnode'} = $_->{'name'};
	}
	my $node_no = 1;


#do the clustering...
	my $tmp2 = 0;
	my @chosen = ();
# do the following until there are only two alignments left
# pick the best pairwise alignment
	while( $r_size >0 ){
#	while(@nodes - @chosen > 1){
		print "-------------------------- new iteration ---------------------------------\n";
		$tmp2++;
		# sort alignments according to some criterion 
		@sortedalns = sort { 
			$b->{norm_score}-($totals{$b->{node1}}/($nodes-2)-$totals{$b->{node1}}/($nodes-2)) <=> $a->{norm_score}-($totals{$a->{node1}}/($nodes-2)-$totals{$a->{node1}}/($nodes-2))
			} @sortedalns;
		#@sortedalns = @alns;
		$r_size = @sortedalns;
		my $tmp1 = shift(@sortedalns);
		$r_size--;

# remove the chosen nodes from the distance matrix		
		push(@chosen, $$tmp1{'node1'});
		push(@chosen, $$tmp1{'node2'});
		@sortedalns = purge_list(@sortedalns, @chosen);
		$r_size = @sortedalns;


# merge the two nodes
		my $cpvec = pvec_avg($tmp1->{'pvec1'}, $tmp1->{'pvec2'},
			$tmp1->{'sw_pairset'}, 2);
		my @memb1 = @{$nodes[$tmp1->{'node1'}]{'members'}};
		my @memb2 = @{$nodes[$tmp1->{'node2'}]{'members'}};
		my @members = @memb1;
		push (@members, @memb2);
		print "@members\n";
		my $ali;
		if( $nodes[$tmp1->{'node1'}]->{'alignment'} && $nodes[$tmp1->{'node2'}]->{'alignment'} ) {
			$ali = merge_alignments($nodes[$tmp1->{'node1'}]->{'alignment'}, $nodes[$tmp1->{'node2'}]->{'alignment'}, $tmp1->{'sw_pairset'});
		}
		elsif(!$nodes[$tmp1->{'node1'}]->{'alignment'} && !$nodes[$tmp1->{'node2'}]->{'alignment'}){
			$ali = $tmp1->{'sw_pairset'};
		}
		else{
			$ali = merge_alignments($nodes[$tmp1->{'node1'}]->{'alignment'}, $tmp1->{'sw_pairset'}, $tmp1->{'sw_pairset'});
			$ali = remove_seq($ali, -2);
		}

# do graph drawing shit here...
		my $node1 = "_".$node_no."_";
		my $mynode = $graph->add_node($node1);
		
		
		$graph->add_edge($node1 => $nodes[$tmp1->{'node1'}]->{'name'}, label => sprintf(" %.1f ", $tmp1->{'norm_score'}/2), style => 'bold', arrowhead => 'none', fontsize => 16);
		$graph->add_edge($node1 => $nodes[$tmp1->{'node2'}]->{'name'}, label => sprintf(" %.1f ", $tmp1->{'norm_score'}/2), style => 'bold', arrowhead => 'none', fontsize => 16);
							  
#compute alignments between the new node and all remaining nodes.
#remaining nodes are the ones in structlist - chosen. 
		for(my $itemidx =0; $itemidx < @nodes; $itemidx++){
			my $done = 0;
			foreach (@chosen){
				if($_ eq $itemidx) {
					$done++;
					#break;
				}
			}
			if (!$done){
				my ($newal) = do_align($cpvec, $nodes[$itemidx]{'pvec'});
				$newal->{'node1'} = $#nodes+1;
				$newal->{'node2'} = $itemidx;
				$newal->{'name'}  = $node1;

				$totals{$newal{node1}} += $$newal{norm_score};
				$totals{$newal{node2}} += $$newal{norm_score};
				push @sortedalns,  $newal;
				$r_size++;
			}
		}

# append the merged alignment to the list of nodes.
		my $newnode = {
				'name' => $node1,
				'pvec' => $cpvec,
				'graphnode' => $node1,
				'alignment' => $ali,
				'members' => \@members,
		};
		push(@nodes, $newnode);
		$node_no++;
		$r_size = @sortedalns;


# lather, rinse, repeat...
	}
	print "Done...\n";
	my @labels = @{$nodes[$#nodes]{'members'}};
	print "Labels: @labels \n";
	print multal_string($nodes[$#nodes]->{'alignment'}) . "\n";
	$graph->as_svg("pretty.svg");




# calculate RMSDs
	my $rmsd = 0.0;
	my $totalsize = 0;
	my $current = 0;
	my $cent_idx=0;
	$current = 0;
	my %l_i;
	foreach my $lbl (@labels){
		$l_i{"$lbl"} = $current;
		$current++;
	}
	foreach my $lbl (@labels){
		my $coord1 = $structs{"$lbl"};
		foreach my $lbl2 (@labels){
			if($lbl ne $lbl2){
				my $coord2 = $structs{"$lbl2"};
				my $ps = split_multal($nodes[$#nodes]->{'alignment'}, $l_i{"$lbl"}, $l_i{"$lbl2"});
	
			#	print multal_string($ps)."\n";
			
				my($t1, $t2) = get_rmsd($ps, $coord1, $coord2);
			#	print "RMSD: $t1\n";
				$rmsd += $t1;
	            $totalsize+= $t2;
			}
			$current++;
		}
		$cent_idx++;
	}

	$rmsd = $rmsd/$totalsize;
	$rmsd = sqrt($rmsd);
	my $asize = length multal_string($nodes[$#nodes]->{'alignment'}) ;
	$asize--;
	$asize /= $cent_idx;
	print "The alignment length is $asize\n";
	#$totalsize = $totalsize/(@$structlist-1);
	print "RMSD is $rmsd over $totalsize residues\n";

    return EXIT_SUCCESS;
}

# ----------------------- aggregate_nj -----------------------------------
# Aggregate alignment library with NJ. 
# ---------------------------------------------------------------------
sub aggregate_nj (\@ \@ $) {
    my ( $structlist, $struct_dirs, $formatflag ) = @_;
	my @nodes = make_nodes( @$structlist );
	my @alns = build_lib(@$structlist, @$struct_dirs, @nodes);
	upgma_check(@alns);
	my $centroid = pop(@alns);
	my $r_size = @alns;
	my @sortedalns = @alns;

	my %totals;
	foreach (@alns){
		$totals{$_{node1}} += $$_{norm_score};
		$totals{$_{node2}} += $$_{norm_score};
	}

	my $graph = GraphViz->new();
	foreach (@nodes) {
		$graph->add_node($_->{'name'}, rank=>'leaf', fontsize => 12, style => 'bold');
		$_->{'graphnode'} = $_->{'name'};
	}
	my $node_no = 1;


#do the clustering...
	my $tmp2 = 0;
	my @chosen = ();
# do the following until there are only two alignments left
# pick the best pairwise alignment
	while( $r_size >0 ){
#	while(@nodes - @chosen > 1){
		print "-------------------------- new iteration ---------------------------------\n";
		$tmp2++;
		# sort alignments according to some criterion 
		@sortedalns = sort { 
			$b->{norm_score}-($totals{$b->{node1}}/($nodes-2)-$totals{$b->{node1}}/($nodes-2)) <=> $a->{norm_score}-($totals{$a->{node1}}/($nodes-2)-$totals{$a->{node1}}/($nodes-2))
			} @sortedalns;
		#@sortedalns = @alns;
		$r_size = @sortedalns;
		my $tmp1 = shift(@sortedalns);
		$r_size--;

# remove the chosen nodes from the distance matrix		
		push(@chosen, $$tmp1{'node1'});
		push(@chosen, $$tmp1{'node2'});
		@sortedalns = purge_list(@sortedalns, @chosen);
		$r_size = @sortedalns;


# merge the two nodes
		my $cpvec;
		if(prob_vec_length($tmp1->{'pvec1'}) > prob_vec_length($tmp1->{'pvec1'})){
			$cpvec = $tmp1->{'pvec1'};
		}
		else {
			$cpvec = $tmp1->{'pvec2'};
		}
		my @memb1 = @{$nodes[$tmp1->{'node1'}]{'members'}};
		my @memb2 = @{$nodes[$tmp1->{'node2'}]{'members'}};
		my @members = @memb1;
		push (@members, @memb2);
		print "@members\n";
		my $ali;
		if( $nodes[$tmp1->{'node1'}]->{'alignment'} && $nodes[$tmp1->{'node2'}]->{'alignment'} ) {
			$ali = merge_alignments($nodes[$tmp1->{'node1'}]->{'alignment'}, $nodes[$tmp1->{'node2'}]->{'alignment'}, $tmp1->{'sw_pairset'});
		}
		elsif(!$nodes[$tmp1->{'node1'}]->{'alignment'} && !$nodes[$tmp1->{'node2'}]->{'alignment'}){
			$ali = $tmp1->{'sw_pairset'};
		}
		else{
			$ali = merge_alignments($nodes[$tmp1->{'node1'}]->{'alignment'}, $tmp1->{'sw_pairset'}, $tmp1->{'sw_pairset'});
			$ali = remove_seq($ali, -2);
		}

# do graph drawing shit here...
		my $node1 = "_".$node_no."_";
		my $mynode = $graph->add_node($node1);
		
		
		$graph->add_edge($node1 => $nodes[$tmp1->{'node1'}]->{'name'}, label => sprintf(" %.1f ", $tmp1->{'norm_score'}/2), style => 'bold', arrowhead => 'none', fontsize => 16);
		$graph->add_edge($node1 => $nodes[$tmp1->{'node2'}]->{'name'}, label => sprintf(" %.1f ", $tmp1->{'norm_score'}/2), style => 'bold', arrowhead => 'none', fontsize => 16);
							  
#compute alignments between the new node and all remaining nodes.
#remaining nodes are the ones in structlist - chosen. 
		for(my $itemidx =0; $itemidx < @nodes; $itemidx++){
			my $done = 0;
			foreach (@chosen){
				if($_ eq $itemidx) {
					$done++;
					#break;
				}
			}
			if (!$done){
				my ($newal) = do_align($cpvec, $nodes[$itemidx]{'pvec'});
				$newal->{'node1'} = $#nodes+1;
				$newal->{'node2'} = $itemidx;
				$newal->{'name'}  = $node1;

				$totals{$newal{node1}} += $$newal{norm_score};
				$totals{$newal{node2}} += $$newal{norm_score};
				push @sortedalns,  $newal;
				$r_size++;
			}
		}

# append the merged alignment to the list of nodes.
		my $newnode = {
				'name' => $node1,
				'pvec' => $cpvec,
				'graphnode' => $node1,
				'alignment' => $ali,
				'members' => \@members,
		};
		push(@nodes, $newnode);
		$node_no++;
		$r_size = @sortedalns;


# lather, rinse, repeat...
	}
	print "Done...\n";
	my @labels = @{$nodes[$#nodes]{'members'}};
	print "Labels: @labels \n";
	print multal_string($nodes[$#nodes]->{'alignment'}) . "\n";
	$graph->as_svg("pretty.svg");

# calculate RMSDs
	my $rmsd = 0.0;
	my $totalsize = 0;
	my $current = 0;
	my $cent_idx=0;
	$current = 0;
	my %l_i;
	foreach my $lbl (@labels){
		$l_i{"$lbl"} = $current;
		$current++;
	}
	foreach my $lbl (@labels){
		my $coord1 = $structs{"$lbl"};
		foreach my $lbl2 (@labels){
			if($lbl ne $lbl2){
				my $coord2 = $structs{"$lbl2"};
				my $ps = split_multal($nodes[$#nodes]->{'alignment'}, $l_i{"$lbl"}, $l_i{"$lbl2"});
	
			#	print multal_string($ps)."\n";
			
				my($t1, $t2) = get_rmsd($ps, $coord1, $coord2);
			#	print "RMSD: $t1\n";
				$rmsd += $t1;
	            $totalsize+= $t2;
			}
			$current++;
		}
		$cent_idx++;
	}

	$rmsd = $rmsd/$totalsize;
	$rmsd = sqrt($rmsd);
	my $asize = length multal_string($nodes[$#nodes]->{'alignment'}) ;
	$asize--;
	$asize /= $cent_idx;
	print "The alignment length is $asize\n";
	#$totalsize = $totalsize/(@$structlist-1);
	print "RMSD is $rmsd over $totalsize residues\n";
    return EXIT_SUCCESS;
}
# ----------------------- aggregate -----------------------------------
# Aggregate alignment library. 
# ---------------------------------------------------------------------
sub aggregate (\@ \@ $) {
    my ( $structlist, $struct_dirs, $formatflag ) = @_;
	my @nodes = make_nodes( @$structlist );
	my @alns = build_lib(@$structlist, @$struct_dirs, @nodes);
	upgma_check(@alns);
	my $centroid = pop(@alns);
	my $r_size = @alns;
	my @sortedalns = @alns;

	my $graph = GraphViz->new();
	foreach (@nodes) {
		$graph->add_node($_->{'name'}, rank=>'leaf', fontsize => 12, style => 'bold');
		$_->{'graphnode'} = $_->{'name'};
	}
	my $node_no = 1;


#do the clustering...
	my $tmp2 = 0;
	my @chosen = ();
# do the following until there are only two alignments left
# pick the best pairwise alignment
	while( $r_size >0 ){
#	while(@nodes - @chosen > 1){
		print "-------------------------- new iteration ---------------------------------\n";
		$tmp2++;
		# sort alignments according to some criterion 
		@sortedalns = sort { $b->{norm_score} <=> $a->{norm_score} } @sortedalns;
		#@sortedalns = @alns;
		$r_size = @sortedalns;
		my $tmp1 = shift(@sortedalns);
		$r_size--;

# remove the chosen nodes from the distance matrix		
		push(@chosen, $$tmp1{'node1'});
		push(@chosen, $$tmp1{'node2'});
		@sortedalns = purge_list(@sortedalns, @chosen);
		$r_size = @sortedalns;


# merge the two nodes
		my $cpvec;
		if(prob_vec_length($tmp1->{'pvec1'}) > prob_vec_length($tmp1->{'pvec1'})){
			$cpvec = $tmp1->{'pvec1'};
		}
		else {
			$cpvec = $tmp1->{'pvec2'};
		}
		my @memb1 = @{$nodes[$tmp1->{'node1'}]{'members'}};
		my @memb2 = @{$nodes[$tmp1->{'node2'}]{'members'}};
		my @members = @memb1;
		push (@members, @memb2);
		print "@members\n";
		my $ali;
		if( $nodes[$tmp1->{'node1'}]->{'alignment'} && $nodes[$tmp1->{'node2'}]->{'alignment'} ) {
			$ali = merge_alignments($nodes[$tmp1->{'node1'}]->{'alignment'}, $nodes[$tmp1->{'node2'}]->{'alignment'}, $tmp1->{'sw_pairset'});
		}
		elsif(!$nodes[$tmp1->{'node1'}]->{'alignment'} && !$nodes[$tmp1->{'node2'}]->{'alignment'}){
			$ali = $tmp1->{'sw_pairset'};
		}
		else{
			$ali = merge_alignments($nodes[$tmp1->{'node1'}]->{'alignment'}, $tmp1->{'sw_pairset'}, $tmp1->{'sw_pairset'});
			$ali = remove_seq($ali, -2);
		}

# do graph drawing shit here...
		my $node1 = "_".$node_no."_";
		my $mynode = $graph->add_node($node1);
		
		
		$graph->add_edge($node1 => $nodes[$tmp1->{'node1'}]->{'name'}, label => sprintf(" %.1f ", $tmp1->{'norm_score'}/2), style => 'bold', arrowhead => 'none', fontsize => 16);
		$graph->add_edge($node1 => $nodes[$tmp1->{'node2'}]->{'name'}, label => sprintf(" %.1f ", $tmp1->{'norm_score'}/2), style => 'bold', arrowhead => 'none', fontsize => 16);
							  
#compute alignments between the new node and all remaining nodes.
#remaining nodes are the ones in structlist - chosen. 
		for(my $itemidx =0; $itemidx < @nodes; $itemidx++){
			my $done = 0;
			foreach (@chosen){
				if($_ eq $itemidx) {
					$done++;
					#break;
				}
			}
			if (!$done){
				my ($newal) = do_align($cpvec, $nodes[$itemidx]{'pvec'});
				$newal->{'node1'} = $#nodes+1;
				$newal->{'node2'} = $itemidx;
				$newal->{'name'}  = $node1;
				push @sortedalns,  $newal;
				$r_size++;
			}
		}

# append the merged alignment to the list of nodes.
		my $newnode = {
				'name' => $node1,
				'pvec' => $cpvec,
				'graphnode' => $node1,
				'alignment' => $ali,
				'members' => \@members,
		};
		push(@nodes, $newnode);
		$node_no++;
		$r_size = @sortedalns;


# lather, rinse, repeat...
	}
	print "Done...\n";
	my @labels = @{$nodes[$#nodes]{'members'}};
	print "Labels: @labels \n";
	print multal_string($nodes[$#nodes]->{'alignment'}) . "\n";
	$graph->as_svg("pretty.svg");

# calculate RMSDs
	my $rmsd = 0.0;
	my $totalsize = 0;
	my $current = 0;
	my $cent_idx=0;
	$current = 0;
	my %l_i;
	foreach my $lbl (@labels){
		$l_i{"$lbl"} = $current;
		$current++;
	}
	foreach my $lbl (@labels){
		my $coord1 = $structs{"$lbl"};
		foreach my $lbl2 (@labels){
			if($lbl ne $lbl2){
				my $coord2 = $structs{"$lbl2"};
				my $ps = split_multal($nodes[$#nodes]->{'alignment'}, $l_i{"$lbl"}, $l_i{"$lbl2"});
	
			#	print multal_string($ps)."\n";
			
				my($t1, $t2) = get_rmsd($ps, $coord1, $coord2);
			#	print "RMSD: $t1\n";
				$rmsd += $t1;
	            $totalsize+= $t2;
			}
			$current++;
		}
		$cent_idx++;
	}

	$rmsd = $rmsd/$totalsize;
	$rmsd = sqrt($rmsd);
	my $asize = length multal_string($nodes[$#nodes]->{'alignment'}) ;
	$asize--;
	$asize /= $cent_idx;
	print "The alignment length is $asize\n";
	#$totalsize = $totalsize/(@$structlist-1);
	print "RMSD is $rmsd over $totalsize residues\n";


    return EXIT_SUCCESS;
}

# ----------------------- restore_handlers --------------------------
# If we are at the stage of mailing, we no longer want to trap
# interrupts. Otherwise, they will call the bad_exit routine again.
sub restore_handlers ()
{
    $SIG{INT } = 'DEFAULT';
    $SIG{QUIT} = 'DEFAULT';
    setpriority ($PRIO_PGRP, getpgrp (0), $low_priority);
}


# ----------------------- mymain  -----------------------------------
# Arg 1 is a structure list file.
sub mymain () {
    use Getopt::Std;
    @struct_dirs = @DFLT_STRUCT_DIRS; 
	my (%opts);
	if ( !getopts( 'a:d:h:s:t:q:l:', \%opts ) ) {
	   usage();
	}
    if ( defined( $opts{d} ) ) { $modeldir     = $opts{d} }
    if ( defined( $opts{l} ) ) { $libfile      = $opts{l} }
    else{
        print STDERR "Please give me a structure library / file\n";
        usage();
    }
    if ( defined( $opts{t} ) ) {
        push( @struct_dirs, split( ',', $opts{t} ) );
    }
    undef %opts;

    set_params();
    #$libfile = $ARGV[0];
    my $formatflag = 'not set';
    check_dirs(@struct_dirs);
    if ( @struct_dirs == 0 ) {
        die "\nNo valid structure directory. Stopping.\n";
    }

    ( @struct_list = get_prot_list($libfile) ) || $fatalflag++;
	if($fatalflag){
		print "FATAL: get_prot_list $libfile \n";
	}
    check_files( @struct_dirs, @struct_list, $bin_suffix ) && $fatalflag++;
    if ($fatalflag) {
        print STDERR "struct dirs were @struct_dirs\n";
    }

    if ($fatalflag) {
        print STDERR "Fatal problems\n";
        return EXIT_FAILURE;
    }
    print
        "Library from $libfile with ", $#struct_list + 1,
        " template structures.\n";

#    my $r =
#      aggregate( @struct_list, @struct_dirs, $formatflag );
#    my $r =
#      aggregate_nj( @struct_list, @struct_dirs, $formatflag );
    my $r =
      aggregate_cons( @struct_list, @struct_dirs, $formatflag );
    if ( $r == EXIT_FAILURE ) {
        bad_exit('calculation broke');
    }

    print
"__________________________________________________________________________\n",
      "Wurst gegessen at ", scalar( localtime() ), "\n";
    my ( $user, $system, $crap, $crap2, $host );
    ( $user, $system, $crap, $crap2 ) = times();
    printf "I took %d:%d min user and %.0f:%.0f min sys time\n", $user / 60,
      $user % 60, $system / 60, $system % 60;
    use Sys::Hostname;
    $host = hostname() || { $host = 'no_host' };
}
mymain();
