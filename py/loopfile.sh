#! /bin/bash
# A. Spiga 30/03/2015

############
field3="u"
field3="temp"
############
res=16
############

############
field3="temp"
field2="ISR"
daz=55
res=720
############


field3="p"
field2="ps"
res=720



## ALTITUDE
if [[ "$daz" != "" ]]; then
  zed="-z $daz"
else
  zed=""
fi

## FIELD 2D and 3D

field=""
nfield=""
if [[ "$field3" != "" ]]; then
  field=$field" -V $field3"
  nfield=$nfield$field3"_"
fi
if [[ "$field2" != "" ]]; then
  field=$field" -v $field2"
  nfield=$nfield$field2"_"
fi

## MAIN LOOP

var=0
list=""

for file in "$@"
do

    var=$((var + 1))
    dafile="file"$var".nc"
    echo $var $file $dafile

    ./remap.py -R -W weights_UHD_to_$res.nc $file $res $field -t 0 $zed -o $dafile

    if [ "$?" -ne "0" ]; then
      echo "FAILED"
    else
      list=$list" "$dafile
      echo $list
    fi

done

ncrcat -O $list $res"_"$nfield$daz.nc
\rm $list
