#
# script to show the steps to partition a surface mesh
# described with SMF file and a nodal 'tag file' (TG file) 
#
# 07/2012
#

# locate openftlroot
if [ ! -n "$OPENFTLROOT" ] ; then
    echo "OPENFTLROOT is not set"
    exit
else
    echo "Using OPENFTLROOT=$OPENFTLROOT"
fi

# input
#basename=${1%.*}
basename=plate
#numpart=$2
numpart=3
numoverlap=2

# tg -> tgl*
tg2tgl=$OPENFTLROOT/tools/input/partitioners/tagFile/tg2tgl
$tg2tgl $basename

# collect tgl*
tgls=""; tags=""
for tgl in `ls plate.tgl*` ; do
    tgls=$tgls${tgl##*.}" "
    if [ -n "$tags" ] ; then tags=$tags"," ; fi
    tags=$tags${tgl##*.*tgl}
done
echo "Tags: $tags"
echo "Tag node list files: $tgls"

# smf -> psmf with node lists of tags
#smf2psmf=$OPENFTLROOT/tools/input/partitioners/metis/partitioner
smf2psmf=$OPENFTLROOT/tools/input/partitioners/rcb/mesh/partitioner
$smf2psmf $numpart $numoverlap $basename $tgls

# basename.<partId>.tgl<tag> -> basename.<partId>.ptg
tgl2ptg=$OPENFTLROOT/tools/input/partitioners/tagFile/tgl2ptg
$tgl2ptg $numpart $basename $tags
