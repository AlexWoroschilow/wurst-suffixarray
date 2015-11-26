use 5.00004;
use ExtUtils::MakeMaker;

if (defined $ENV{CC}) {
    $CC = $ENV{CC}; }

if (defined $ENV{OPTIMIZE}) {
    $OPTIMIZE = $ENV{OPTIMIZE}; }

if (defined $ENV{MOD_OBJS}) {
    $MOD_OBJS = $ENV{MOD_OBJS};
    $MOD_OBJS =~ s/\s*(\S+)/ ..\/..\/\1/g;
    $MOD_OBJS = $MOD_OBJS." ../../../../../src/libwurst.a";
}


if (defined $ENV{MANDEST}) {
    $man_dest = $ENV{MANDEST};
} else {
    $man_dest = "$ENV{HOME}/pl/lib"; }

if ($CC eq 'gcc') {
    $OPTIMIZE="$OPTIMIZE"; }

WriteMakefile(
              'NAME'		=> 'Wurst::Guts',
              'VERSION_FROM'	=> 'Guts.pm', # finds $VERSION
              'PREREQ_PM'		=> {'Wurst'=>0.01}, # e.g., Module::Name => 1.1
              ($] >= 5.005 ?    ## Add these new keywords supported since 5.005
               (ABSTRACT_FROM => 'Guts.pm', # retrieve abstract from module
                AUTHOR     => 'James Procter <procter@suse.de>') : ()),
              'LIBS'		=> ['-lm'], # e.g., '-lm'
              'MYEXTLIB' => $MOD_OBJS, 
              'OPTIMIZE' => $OPTIMIZE,
              'DEFINE'		=> '', 
              'INC'		=> '-I. -I../../ -I../../../../../src/', 
              XSPROTOARG => '-prototypes'              
              );

if  (eval {require ExtUtils::Constant; 1}) {
  # If you edit these definitions to change the constants used by this module,
  # you will need to use the generated const-c.inc and const-xs.inc
  # files to replace their "fallback" counterparts before distributing your
  # changes.
  my @names = (qw());
  ExtUtils::Constant::WriteConstants(
                                     NAME         => 'Wurst::Guts',
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
