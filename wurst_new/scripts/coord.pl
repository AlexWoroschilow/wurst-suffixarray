# 23 oct 2001
# Play with coordinate reading routines.
# rcsid = $Id: coord.pl,v 1.3 2007/05/03 08:34:06 mahmood Exp $


use warnings;

use FindBin;
use lib "$FindBin::Bin/../src/Wurst/blib/arch";
use lib "$FindBin::Bin/../src/Wurst/blib/lib";
#use lib "$ENV{HOME}";  # Where the wurst stuff lives after installation
use Wurst;

use strict;
use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);



use vars qw($COORD_DIR);
*COORD_DIR = \"/bm/pdb90_bin";
my $coord_name   = '1b4aA.bin';
my $out_pdb_base = 'del_me';


# ----------------------- mymain  -----------------------------------
sub mymain ()
{

    $SIG{TRAP} = 'IGNORE'; kill 'TRAP', $$;
    my $name = $COORD_DIR . '/' . $coord_name;
    my $coord = coord_read ( $name) || die "coord_read fail on $name";
#   $SIG{TRAP} = 'IGNORE';    kill 'TRAP', $$;
    my $out_pdb = "$out_pdb_base.pdb";
    coord_2_pdb ($out_pdb, $coord) || die "Fail on $out_pdb: $!";
    coord_play ($coord);
    unlink ($out_pdb) || die "Fail deleting $out_pdb: $!";
    $out_pdb = "${out_pdb_base}_seq.pdb";
    coord_2_pdb ($out_pdb, $coord, coord_get_seq ($coord)) || die "Failed writing $out_pdb";
#   unlink ($out_pdb) || die "Fail deleteing $out_pdb: $!";
    return EXIT_SUCCESS;
}
exit mymain();
