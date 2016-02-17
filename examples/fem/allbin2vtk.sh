#!/bin/bash
shopt -s nullglob
for f in data/*-p0000-*.bin
do
  mpiexec -n 4 bin2vtk $f
done
