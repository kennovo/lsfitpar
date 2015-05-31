#!/bin/sh
if [ $# -ne 1 ] ; then
  echo "Usage: `basename $0` <lsfitpar binary>"
  exit 1
 fi
if [ -d "$1" -o ! -x "$1" ] ; then
  echo "ERROR: $1 not an excutable!"
  exit 1
 fi

for i in `seq 9 -1 1` ; do
  for j in 1 2 5 ; do
    exp=${j}e-$i
    sed "s/%bias%/$exp/" illcond.lsarg.template | xargs $1 >illcond$exp.ou
    echo -n "$exp "
    awk '/^a1/{printf("%s ",$5);next}
      /^a2/{if ($7 == "180.00") printf("-%s %s\n",$5,$11)
         else printf("%s %s\n",$5,$11)
        exit}' illcond$exp.prm
   done
 done > illcond.k
