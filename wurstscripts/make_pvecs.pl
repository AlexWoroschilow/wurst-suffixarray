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
#	print "$line \n";
	my $vname = "~/andrew/scripts/pvecs/$line.vec";
        if (! -e $vname){
	    push (@a, $line);
	    print "$vname\n";
	}
    }
    close (F);
    return (@a);
}

# ----------------------- mymain  -----------------------------------
sub mymain ()
{
    my $fatalflag = undef;
    my $filename;
    my $workdir;
    my $coord;
    my $vecname;
    my $classfcn;
    my $pvec;
    my $list;
    my $gauss_err = 0.4;
    my $classfile = 'classfile';   
    if ( $#ARGV < 0) {
	$workdir = ".";
	
    }
    else {
	print "$ARGV[0]\n";
	$list= $ARGV[0];
	
    }
    $classfcn = aa_strct_clssfcn_read($classfile, $gauss_err);
    my @structlist  = get_prot_list("$list");
    for (my $i = 0; $i < scalar @structlist ; $i++) {
       $filename = "/home/stud2004/tmargraf/pdbsnapshot060307/$structlist[$i].bin";
	   $coord = coord_read($filename);
	   $pvec = strct_2_prob_vec($coord, $classfcn);
	   $filename =~ s/\/home\/stud2004\/tmargraf\/pdbsnapshot060307\///g;
	   $vecname = $filename;
	   $vecname =~ s/\.bin/\.vec/g;
	if(!prob_vec_write($pvec, "pvecs/$vecname")){
	    print "FAILED ";
	}
	else{
	    print "SUCCESS ";
	}
	$filename =~ s/\.bin//g;
	
	print "$filename \n";
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
    return EXIT_SUCCESS;
}
# ----------------------- main    -----------------------------------
exit (mymain());
