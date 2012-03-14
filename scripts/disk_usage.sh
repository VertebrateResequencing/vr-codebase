#!/bin/sh
# print vertres disk usage
echo "Type       Mounted on               Size     Used    Avail     Use%"

for d in /lustre/scratch102 /lustre/scratch105 /lustre/scratch106 `ls -d /nfs/vertres*` /warehouse/g1k-01 /warehouse/g1k-02 /warehouse/g1k-03 /warehouse/g1k-04 
do
  cd $d
  df -Th . | awk 'NF==6{printf ("%-10s %-20s %8s %8s %8s %8s\n", $1, $6, $2, $3, $4, $5)}'
done
