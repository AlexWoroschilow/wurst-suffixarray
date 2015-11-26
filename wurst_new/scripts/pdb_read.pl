# 27 March 2002
# For testing pdb reading routines.
# $Id: pdb_read.pl,v 1.4 2007/05/26 18:26:43 shoffmann Exp $

use FindBin;
use lib "../src/Wurst/blib/arch";
use lib "../src/Wurst/blib/lib";

#use lib "$FindBin::Bin/../src/Wurst/blib/arch";
#use lib "$FindBin::Bin/../src/Wurst/blib/lib";

use strict;
use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);

use Wurst;

# ----------------------- constants ---------------------------------
use vars qw ($pdbdir @pdbfiles @chains);

sub set_constants ()
{
    if (at_home()) {
        $pdbdir = undef;

        @pdbfiles = ('../struct/pdb1lea', '../struct/1craP.pdb');
        @chains   = ('', '');
    } else {
        $pdbdir = '/projects/biodata/pdb/data/structures/divided/pdb';
        @pdbfiles = ('pdb5pgm.ent.gz', '../struct/pdb1lea.ent.gz', 'pdb2pgk.ent.gz');
        @chains   = ('C',                 ' ',          '_');
    }
}

# ----------------------- name_from_pdb -----------------------------
sub name_from_pdb ($ $)
{
    my $name;     # The name, with chain, that we will return
    my ($pdbfile, $chain) = @_;
    $name = $pdbfile;
    $name =~ s/^.*\/*ent//i;     # Leading path and "pdb" removal
    $name =~ s/\..*//;           # pdb1abc.pdb.ent.gz to pdb1abc
    my $path = $pdbdir;
    my $twocode = $name;
    $twocode = substr ($name, 1, 2);
    if ( ! ($pdbfile =~ /\//)) {   # pdb file does not have leading path
        $path = "${path}/${twocode}/";
        $path = "${path}/${pdbfile}";
        $path =~ s/\/\/+/\//g;   # Get rid of multiple slashes
    } else {
        $path = $pdbfile;
    }
    if (! -f $path) {
        print STDERR "Cannot find $path\n";
        return undef;
    } else {
        return ($path, $name);
    }
}

# ----------------------- at_home    --------------------------------
# Don't know if this should be in the checked in code..
sub at_home ()
{
    my $x = `uname -n`;
    if ($?) {
        print STDERR "uname fail ", $? >>8, "\n";
        return 0;
    }

    if ($x =~ m/kartoffel/) {
        return 1;
    } else {
        return 0;
    }
}

# ----------------------- get_uncompressed --------------------------
# If the file seems to be compressed, uncompress it to a temp file.
# Return a name that can be deleted. Either a temporary file, or just
# make a link to the original file.
sub get_uncompressed ($)
{
    my $pdbfile = shift;
    my $gunzip = 'gunzip';
    my $temppdb = "delete_me.$$";
    if ( ! -f $pdbfile ) {
        die "Cannot find $pdbfile as pdbfile"; }
    if ( $pdbfile =~ /\.gz$/) {
        my $cmdline = "$gunzip -c < $pdbfile > $temppdb";
        system ($cmdline) && die "$gunzip returned ", $? >> 8, ". Stopping";
    } else {
        symlink ($pdbfile, $temppdb) || die "link failed. $!";
    }
    return $temppdb;
}

# ----------------------- one_file ----------------------------------
sub one_file ($ $)
{
    use File::Basename;
    my ($pdbfile, $chain) = @_;
    my $temppdb = get_uncompressed ($pdbfile);
    print "temppdb is $temppdb\n";
    my $mdl = pdb_read ($temppdb, '', $chain) || die "pdb_read fail: $!";
    unlink ($temppdb) || die "unlink fail: $!";
    $DB::single = 1;
    my ($file, $path, $ext) = fileparse ($pdbfile, '\..+$');
    my $out = "del_me_${file}.pdb";
    coord_2_pdb ($out, $mdl) || die "Failed writing $out";
}

# ----------------------- mymain  -----------------------------------
sub mymain ()
{
    set_constants();
    func_int; # force loading of wurst
    $SIG{TRAP} = 'IGNORE';# kill 'TRAP', $$;
    for (my $i = 0; $i < @pdbfiles; $i++) {
        my $p = $pdbfiles[$i];
        my $chain = $chains[$i];
        my ($path, $name);
        ($path, $name) = name_from_pdb ($p, $chain);
        if ( ! defined ($path)) {
            print STDERR "Giving up on $p\n";
            next;
        }
        print "Now doing $path $chain\n";
        one_file ($path, $chain);
    }
    return EXIT_SUCCESS;
}
exit mymain();
