#! /bin/bash
# A. Spiga 30/03/2015

#wfile="weights_UHD_to_240_equator.nc"



################################################
opt="-t 99 -z 35 -V u -V temp -V p -v ISR"
res=720
target="all720.nc"
################################################


## WEIGHTS
if [[ "$wfile" == "" ]]; then
  wfile="weights_UHD_to_$res.nc"
fi

###############
## MAIN LOOP ##
###############

var=0
list=""

for file in "$@"
do

    var=$((var + 1))
    dafile="file"$var".nc"
    echo $var $file $dafile

    # check time flag
    ncdump -h $file | tail -n 2 | head -n 1

    ## get sol number
    #num=`basename $file '.nc' | sed s/'xios_diagfi_'/''/g | sed s/'-'/' '/g | awk '{print $2}'`
    #kron=$((num / 38052))

    # REMAP REMAP REMAP
    ./remap.py -R -W $wfile $file $res $field $opt -o $dafile > /dev/null

    # populate list
    if [ "$?" -ne "0" ]; then
      echo "FAILED"
    else
      list=$list" "$dafile
    fi

done

ncrcat -O $list $target
\rm $list
