#!/bin/bash
# Usage: ./script_biopac <python script for analysis> [<CRT participant directory>] [<CRT analysis storage directory>]

# Script absolute path for python script for analysis
script_analysis=""
# CRT participant directory
data_dir=""
# CRT analysis storage directory
dest_dir=""
# Suffix for directory and for files after analysis
suffix="_analysis"
# Separater for filename and extension
separater="."
if [ "$#" -lt 1 ]
then
    echo "No python script provided for analysis; exiting..."
    exit 0
else
    # For making absolute path
    script_analysis="$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
    # echo $script_analysis
    script_analysis=$(echo $script_analysis | sed 's/ /\\ /g')
    # echo $script_analysis
    script_analysis=$(echo $script_analysis | awk '{gsub( /[(&*|`$)]/, "\\\\&"); print $0}')
    # echo $script_analysis
    # script_analysis=$1
fi
curr_dir="`pwd`"
curr_dir="$(cd "$(dirname "$curr_dir")"; pwd)/$(basename "$curr_dir")"
curr_dir=$(echo $curr_dir | sed 's/ /\\ /g')
curr_dir=$(echo $curr_dir | awk '{gsub( /[(&*|`$)]/, "\\\\&"); print $0}')
if [ "$#" -lt 2 ]
then
    echo "No CRT participant directory provided; Setting current directory as source directory"
    data_dir=$curr_dir
else
    data_dir="$(cd "$(dirname "$2")"; pwd)/$(basename "$2")"
    data_dir=$(echo $data_dir | sed 's/ /\\ /g')
    data_dir=$(echo $data_dir | awk '{gsub( /[(&*|`$)]/, "\\\\&"); print $0}')
    echo "Source directory:" $data_dir
fi
dest_dir=""
if [ "$#" -lt 3 ]
then
    echo "No destination directory provided; Setting source directory as destination directory"
    dest_dir=$data_dir
else
    dest_dir="$(cd "$(dirname "$3")"; pwd)/$(basename "$3")"
    dest_dir=$(echo $dest_dir | sed 's/ /\\ /g')
    dest_dir=$(echo $dest_dir | awk '{gsub( /[(&*|`$)]/, "\\\\&"); print $0}')
    echo "Destination directory:" $dest_dir
fi

# echo $data_dir
eval "cd $data_dir"
# eval "ls -la"

for D in `find . -mindepth 1 -type d`
do
    echo -e
    echo $D
    # eval "ls -l $D"
    eval "mkdir $dest_dir/$D$suffix"
    for F in `find ./$D/crt*.txt`
    do
        filepath="$(cd "$(dirname "$F")"; pwd)/$(basename "$F")"
        filepath=$(echo $filepath | sed 's/ /\\ /g')
        filepath=$(echo $filepath | awk '{gsub( /[(&*|`$)]/, "\\\\&"); print $0}')
        filename=$(basename "$filepath")
        extension="${filename##*$separater}"
        filename="${filename%$separater*}"
        echo $filepath
        # echo $filename $extension
        # eval "cd $curr_dir"
        # echo "python $script_analysis $filepath >> $dest_dir/$D$suffix/$filename$suffix$separater$extension"
        eval "python $script_analysis $filepath >> $dest_dir/$D$suffix/$filename$suffix$separater$extension"
        # eval "cd $data_dir"
    done
    # break
done
