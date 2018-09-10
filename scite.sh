#!/bin/bash

OPTIND=1

unset MUTATIONS
unset SAMPLES
SCORETYPE='m'
TREETYPE='m'
TRUETREE=false

while getopts ":n:m:t:s" opt
do
    case "$opt" in
        t)
            if [[ $OPTARG == "ranspose" ]]
            then
                TREETYPE='t'
            else
                TRUETREE=true
            fi
            ;;
        s)
            SCORETYPE='s'
            ;;
        n)
            MUTATIONS=$OPTARG
            ;;
        m)
            SAMPLES=$OPTARG
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
    esac
done

if [[ -z ${MUTATIONS+x} ]]
then
    echo "Missing parameter n."
    exit 1
fi

if [[ -z ${SAMPLES+x} ]]
then
    echo "Missing parameter m."
    exit 1
fi

if [[ $TRUETREE = true ]]
then
    TRUETREE="-DTRUE_TREE"
else
    unset TRUETREE
fi

if [[ -z ${CXX+x} ]]
then
    echo "Unable to find C++ compiler. Please set CXX via export CXX=."
    exit 1
fi

$CXX --std=c++17 -march=native -O3 src/findBestTrees.cpp -DMUTATIONS=$MUTATIONS -DSAMPLES=$SAMPLES -DTREETYPE=\'$TREETYPE\' -DSCORETYPE=\'$SCORETYPE\' $TRUETREE -o scite || exit 1
shift
./scite "$@"
