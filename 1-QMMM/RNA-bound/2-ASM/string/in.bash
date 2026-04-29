#!/bin/bash

# If groupfile exists, remove it
if [ -f string.groupfile ];
then
  rm string.groupfile
fi

PARM=parameters.prmtop  # path to topology file
REACT=reactants.rst   # path to reactants structure
PROD=products.rst     # path to products structure

# Number of string nodes is provided as command line argument
nodes=$1
for i in $(seq 1 $nodes); do
  # generate sander input for node $i 
  sed "s/XXXXX/$RANDOM/g" in > $i.in

  # use reactants structure for the first half of the nodes and products for the rest
  if [ $i -le $((nodes/2)) ];
  then
    crd=$REACT
  else
    crd=$PROD
  fi

  # write out the sander arguments for node $i to the groupfile
  echo "-O -rem 0 -i $i.in -o $i.out -c $crd -r $i.rst -x $i.nc -inf $i.mdinfo -p $PARM" >> string.groupfile
done
echo "N REPLICAS  = $i"   
