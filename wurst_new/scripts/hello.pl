# A very important test program
use FindBin;
use lib "$FindBin::Bin/../src/Wurst/blib/arch";
use lib "$FindBin::Bin/../src/Wurst/blib/lib";

use strict;

use Wurst;

print "Hello. I think I am ready to try wurst.\n";
print "Calling func_int. Should return 42. It returns ", func_int(), "\n";
print "Calling func_float. Should return 3.14. Returns ", func_float(), "\n";
print "Calling func_char. Returns \"", func_char(), "\"\n";

exit;
