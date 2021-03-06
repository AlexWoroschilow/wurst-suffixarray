#!/bin/bash

# example:
# startFarmJob.sh -pdbFilePath /Users/ap3/WORK/PDB -nrAlignments 10

# send the arguments to the java app
# allows to specify a different config file
args="$*"

JDIR=$PWD;

cpath="$JDIR/core.jar:$JDIR/structure.jar:$JDIR/structure-gui.jar:$JDIR/JmolApplet.jar:$JDIR/javaws.jar:$JDIR/biojava3-core.jar:$JDIR/biojava3-alignment.jar"

#$OSG_APP/engage/jdk1.6.0_16/bin/java -Xmx1G -cp $cpath org.biojava.bio.structure.align.FarmJob $args
#java -Xmx1G -cp $cpath org.biojava.bio.structure.align.FarmJob $args


if [ -f $OSG_APP/engage/jdk1.6.0_16/bin/java ]; then
    $OSG_APP/engage/jdk1.6.0_16/bin/java -Xmx1G -cp $cpath org.biojava.bio.structure.align.FarmJob $args
else 
    if [ -f /osg/osg-app/engage/jdk1.6.0_03/bin/java ]; then
	/osg/osg-app/engage/jdk1.6.0_03/bin/java  -Xmx1G -cp $cpath org.biojava.bio.structure.align.FarmJob $args
     else
	which java
	java -version
	java  -Xmx1G -cp $cpath org.biojava.bio.structure.align.FarmJob $args
    fi
fi

exit $?