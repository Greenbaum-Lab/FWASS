#!/bin/sh
module load python
cmd2run=$1
echo "$cmd2run"
eval "$cmd2run"
