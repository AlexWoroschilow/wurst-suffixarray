# 6 March 2002
# $Id: paths.pl,v 1.2 2002/04/04 04:45:48 torda Exp $

*DFLT_BIN_DIR = ['/rsc/appenzeller/data1/pdblib6apr00',
                 '/rsc/appenzeller/data1/pdb_bin_pool' ];
*DFLT_PHD_DIR = ["$ENV{HOME}/phd/results"];

use FindBin;
my $where_I_am = $FindBin::Bin;
*MATRIX_DIR  = \"$where_I_am/../../matrix";
*MATRIX_FILE = \'blosum62.50';

*FX_PARAM_DIR  = \"$where_I_am/../../params";
*FX_PARAM_FILE = \'g7.param';
