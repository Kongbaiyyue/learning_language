#!/bin/bash

i=0
`mkdir -p output`
for file_name in `ls zinc-10m*.train.*.pt`
    do
        echo ${file_name}
        mv ${file_name} "./output/zinc-10m.train.${i}.pt"
        ((i++))
    done

