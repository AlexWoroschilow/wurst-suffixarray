use 5.00004;
use ExtUtils::MakeMaker;


if (defined $ENV{CC}) {
    $CC = $ENV{CC}; }

if (defined $ENV{OPTIMIZE}) {
    $OPTIMIZE = $ENV{OPTIMIZE}; }

if (defined $ENV{MANDEST}) {
    $man_dest = $ENV{MANDEST};
} else {
    $man_dest = "$ENV{HOME}/pl/lib"; }

if ($CC eq 'gcc') {
    $OPTIMIZE="$OPTIMIZE"; }



WriteMakefile(
    'NAME'		=> 'Wurst::Pack',
    'VERSION_FROM'	=> 'Pack.pm', 
    'PREREQ_PM'		=> {'Wurst' => 0.01 }, 
    ($] >= 5.005 ?    
      (ABSTRACT_FROM => 'Pack.pm', # retrieve abstract from module
       AUTHOR     => 'James Procter <procter@suse.de>') : ()),

              'OPTIMIZE' => $OPTIMIZE,
              
    'LIBS'		=> ['-lm'],
    'DEFINE'		=> '', 
    'INC'		=> '-I./ -I../../../../../src/ -I../../',
    'MYEXTLIB'=> "../../packw_obj.o ../../../../../src/libwurst.a",

    XSPROTOARG => '-prototypes'
);
if  (eval {require ExtUtils::Constant; 1}) {
  # If you edit these definitions to change the constants used by this module,
  # you will need to use the generated const-c.inc and const-xs.inc
  # files to replace their "fallback" counterparts before distributing your
  # changes.
  my @names = (qw());
  ExtUtils::Constant::WriteConstants(
                                     NAME         => 'Wurst::Pack',
                                     NAMES        => \@names,
                                     DEFAULT_TYPE => 'IV',
                                     C_FILE       => 'const-c.inc',
                                     XS_FILE      => 'const-xs.inc',
                                  );

}
else {
  use File::Copy;
  use File::Spec;
  foreach my $file ('const-c.inc', 'const-xs.inc') {
    my $fallback = File::Spec->catfile('fallback', $file);
    copy ($fallback, $file) or die "Can't copy $fallback to $file: $!";
  }
}
