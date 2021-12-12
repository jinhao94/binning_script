#!/bin/bash
# Author：Jinhao 

function display_help() {
    echo ""
    echo -e "A pipeline for preparing the input files for cluster binning pipeline"
    echo -e "Author jinh"
    echo "-_-!!!"
    echo "   -a, --assembly_folder      Folder for input assembly fasta. [*Required]"
    echo "   -o, --outfolder            Output folder.  [*Required]"
    echo "   -e, --extension            extension for input fasta files (Default: fa)"
    echo "   -threads, --threads        Number of threads will be run. (Default: 30)"
    echo "   -h, --help                 Show this message"
    echo " "
    exit 1
}

# echo $#
if [ $# -eq 0  ];then
    echo "Please input parameters!"
    display_help
fi

## default settings.
assembly_folder="None"; outfolder="None"
extension="fa"
#
while [ "$1" != "" ]; do
    case $1 in
        -a | --assembly_folder ) shift
                                  assembly_folder=$1
                                  ;;
        -o | --outfolder)      shift
                                  outfolder=$1
                                  ;;
        -e | --extension)         shift
                                  extension=$1
                                  ;;
        -t | --threads)           shift
                                  threads=$1
                                  ;;
        * )                       display_help
                                  exit 1
    esac
    shift
done

##############
if [ $assembly_folder == "None" ] || [ $outfolder == "None" ]; then
    echo "Please input required parameters!"
    display_help
fi

new_sample_list=${outfolder}.snk.yaml
tree -if $assembly_folder | grep $extension | rush -j 1 echo  $'{%.} $PWD/{}' | perl -lane '$o=(join "\t", @F[0..$#F]); print "$o"' > $new_sample_list

source activate snakemake
snakemake -s /mnt1/script/script/Bin_cluster_file_prepare.py --config workdir=${outfolder} file_names_txt=$PWD/${new_sample_list} -r -p --cores $threads -j $threads
