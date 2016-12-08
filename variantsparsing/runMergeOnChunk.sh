#!/bin/bash

set -e

chrm=$1
start=$2
end=$3
fn_output=$4
run_bcf=$5
fn_tcga=$6
fn_icgc=$7


$run_bcf merge -r $chrm:$start-$end $fn_tcga $fn_icgc -o $fn_output -O z
$run_bcf index $fn_output
