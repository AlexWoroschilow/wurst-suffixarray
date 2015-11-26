# 25 Dec 2001
use FindBin;
use lib "$ENV{HOME}/pl/lib";
eval ' use Wurst;';
if ($@) {
    print "I cannot find the Wurst library\n",
    "Perhaps it is somwhere like $ENV{HOME}/pl/lib\n",
    "Then please try a line like\n",
    "setenv PERL5LIB=$ENV{HOME}/pl/lib\n",
    "or\n",
    "PERL5LIB=$ENV{HOME}/pl/lib; export PERL5LIB\n";
    exit (1);
  }

print "Hello, no errors finding wurst\n";
print "Hello. I think I am ready to try wurst.\n";
print "Calling func_int. Should return 42. It returns ", func_int(), "\n";
print "Calling func_float. Should return 3.14. Returns ", func_float(), "\n";
print "Calling func_char. Returns \"", func_char(), "\"\n";

exit;
