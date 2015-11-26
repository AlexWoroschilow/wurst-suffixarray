# This is the real source for make'ing. Makefile is derived
# from this.
# To play with:
#  * set $CC
#  * set $OPTIMIZE
# Electric Fence:
# The behaviour is very version and machine specific. As of
# today, there is no combination of options here which will link
# in electric fence so that perl uses it. Instead, use LD_PRELOAD.


use ExtUtils::MakeMaker;

#my $efence;
#if ( -f '/usr/lib/libefence.a') {
#    $efence = '-lefence'; }

if (defined $ENV{CC}) {
    $CC = $ENV{CC}; }

if (defined $ENV{OPTIMIZE}) {
    $OPTIMIZE = $ENV{OPTIMIZE}; }

if (defined $ENV{MANDEST}) {
    $man_dest = $ENV{MANDEST};
} else {
    $man_dest = "$ENV{HOME}/pl/lib"; }

WriteMakefile(
      'OPTIMIZE'      => $OPTIMIZE,
      'NAME'          => 'Wurst',
      'VERSION_FROM'  => 'Wurst.pm',     # finds $VERSION
      #'LIBS'          => ['-lm'],       # e.g., '-lm'
      'MYEXTLIB'      => '../libwurst.a ../gsldir/libminigsl.a',
      'INC'           => '-I..',
      'INSTALLMAN3DIR'=> $man_dest,
      'clean'         => { 'FILES' => '*~'},
);

