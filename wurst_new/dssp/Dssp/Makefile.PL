use ExtUtils::MakeMaker;


my $efence;
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
    'OPTIMIZE'  => "$OPTIMIZE",
    'NAME'	=> 'Dssp',
    'VERSION_FROM' => 'Dssp.pm', # finds $VERSION
    'LIBS'	=> ['-lm'],
    'MYEXTLIB'  => '../libdssp.a ../../src/libwurst.a',
    'INC'	=> '-I..',
);
