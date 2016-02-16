#!/bin/bash
shopt -s nullglob
for f in data/*.bin
do
  ./bin2vtk $f
done
