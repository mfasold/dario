#!/bin/bash
# Decompress several types of (uploaded) compressed files and renames into desired filename

if [ $# != 4 ]; then
    echo "Usage: extract RUNDIR UPLOAD_DATA_FILE ORIGINAL_FILENAME RESULT_FILENAME">&2
    exit 1
fi

rundir=$1
upload_file=$2
original_filename=$3
original_filename_lower=${original_filename,,} # lowercase filename to allow .ZIP, .GZ etc.
result_filename=$4

# Check if file exists and is a supported archive
if [ ! -f $rundir/$upload_file ] ; then
		echo "$0: file $rundir/$upload_file not found. Exiting." >&2
		exit 1
fi

if [ "$(echo $original_filename_lower | egrep '(\.gz$|\.zip$|\.tar$|\.tgz$|\.tar\.gz$)')" = "" ] ; then
		echo "Skipped file ${original_filename}: it's already decompressed (or an unsupported archive)." 
		exit 0
fi


# Create a subdirectory of rundir and decompress ZIP contents into it
unzip_folder=$rundir/ziptmp
mkdir $unzip_folder

case $original_filename_lower in 
    *.tar.gz)
        tar -C $unzip_folder -xvf $rundir/$upload_file
        ;;
    *.gz)
        cp $rundir/$upload_file $unzip_folder/$upload_file.gz
        gunzip -N $unzip_folder/$upload_file.gz
        ;;
    *.zip)
        unzip -d $unzip_folder $rundir/$upload_file
        ;;
esac

# Print error if more than one file is contained
decompressed_file=`find $unzip_folder -type f` # note: contains full path
number_of_files=`find $unzip_folder -type f|wc -l`

if [ $number_of_files -eq 0 ]; then
    echo "Error: Archive could not be decompressed or contained no files.">&2
    exit 1   
fi

if [ $number_of_files -gt 1 ]; then
    echo "Error: Archive contains more than one file.">&2 
    exit 1   
fi

# TODO : handle stupid windows or mac index files 

# Move this file into rundir and rename to upload.bam or upload.bed
decompressed_file_lower=${decompressed_file,,} # lowercase filename to allow .ZIP, .GZ etc.
case "$decompressed_file_lower" in 
    *.bed)
        mv $decompressed_file $rundir/$result_filename.bed
        ;;
    *.bam)
        mv $decompressed_file $rundir/$result_filename.bam
        ;;
    *.upload)
        # it can happen that gzip files do not contain filename - then we have to use original filename
        echo "Warning: Archive does not contain original filename.">&2       
        case $original_filename_lower in 
            *.bed.*)
                mv $decompressed_file $rundir/$result_filename.bed
                ;;
            *.bam.*)
                mv $decompressed_file $rundir/$result_filename.bam
                ;;
        esac
        ;;
    *)
        echo "Warning: Archive does not file with bed or bam suffix. Assuming it is a BED file.">&2
        mv $decompressed_file $rundir/$result_filename.bed
        ;;
esac

rm -r $unzip_folder
