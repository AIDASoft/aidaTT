#!/bin/bash

echo "changing $1: rotate by 270deg, crop and replace original"

echo

pdf270 --outfile ${1}_BAK $1
pdfcrop ${1}_BAK $1
rm ${1}_BAK
