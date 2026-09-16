# environment
if [ ! -n "$OPENFTLROOT" ] ; then
    if [ ! -n "$STATFEMROOT" ] ; then
        echo "Neither \$OPENFTLROOT nor \$STATFEMROOT is set"
        exit
    fi
    OPENFTLROOT=$STATFEMROOT
fi

# make mesh
name=out/bridge
smfs=`ls ${name}*smf`
# convert SMF to VTU
for s in $smfs ; do
    $OPENFTLROOT/tools/input/smf2vtu/smf2vtu $s
done
echo "---> generated VTU files"

