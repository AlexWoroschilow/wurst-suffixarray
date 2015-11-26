use 5.006;
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
              'NAME'		=> 'Wurst::AlDb',
              'VERSION_FROM'	=> 'AlDb.pm', # finds $VERSION
              'PREREQ_PM'		=> { 'Wurst'=>0.01,
                                             'Class::DBI' => 0.95
                                             }, # e.g., Module::Name => 1.1
              PM => { 'AlDb.pm' => '$(INST_LIBDIR)/AlDb.pm',
                      'AlDbUtils.pm', => '$(INST_LIBDIR)/AlDbUtils.pm',
                      'SeqStrCmp.pm', => '$(INST_LIBDIR)/SeqStrCmp.pm',
                      'AlDb/Alignment.pm', => '$(INST_LIBDIR)/AlDb/Alignment.pm',
                      'AlDb/StrAlgn.pm', => '$(INST_LIBDIR)/AlDb/StrAlgn.pm'},
              
#              PMLIBDIRS => ['./','AlDb'],
              ($] >= 5.005 ?    ## Add these new keywords supported since 5.005
               (ABSTRACT_FROM => 'AlDb.pm', # retrieve abstract from module
                AUTHOR     => 'Torda/Procter <procter@zbh.uni-hamburg.de>') : ()),
              );
