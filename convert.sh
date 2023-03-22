#!/bin/bash

#!/bin/bash
    
## Author(s): Yinxiu Zhan
## Contact: yinxiu.zhan@ieo.it
## This software is distributed without any guarantee under the terms of the GNU General
## Public License, either Version 2, June 1991 or Version 3, June 2007.


function usage {
    echo -e "usage : $(basename $0) -i input_mat -b binsize -s start -e end [-r running_average]  [-h]"
    echo -e "Use option -h|--help for more information"
}

function help {
    usage;
    echo 
    echo "Convert HiTC matrix into Matlab suitable format."
    echo "See HiTC documentation for format https://bioconductor.org/packages/release/bioc/html/HiTC.html"
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "   -i|--input input_mat : Hi-C matrix in HiTC format"
    echo "   -b|--binsize BINSIZE : binsize of Hi-C matrix in HiTC format"
    echo "   -e|--end END : the end coordinate of the region to extract"
    echo "   -s|--start START : the start coordinate of the region to extract"
    echo "   [-r|--running_average RUNINNG_AVERAGE] : the running average of the matrix. Default 0."
    echo "   [-h|--help]: help"
    exit;
}


# Transform long options to short ones
for arg in "$@"; do
  shift
  case "$arg" in
      "--input") set -- "$@" "-it" ;;
      "--binsize") set -- "$@" "-b" ;;
      "--end") set -- "$@" "-e" ;;
      "--running_average") set -- "$@" "-r" ;;
      "--start") set -- "$@" "-s" ;;
      "--help")   set -- "$@" "-h" ;;
       *)        set -- "$@" "$arg"
  esac
done

RV=0
INPUT=""

while getopts ":i:b:e:s:r:h" OPT
do
    case $OPT in
		i) INPUT=$OPTARG;;
        b) BINSIZE=$OPTARG;;
        e) END=$OPTARG;;
        r) RV=$OPTARG;;
        s) START=$OPTARG;;
        h) help ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            usage
            exit 1
            ;;
    esac
done

if [ $# -lt 8 ]
then
    usage
    exit
fi

if ! [ -f $INPUT ]; then
	echo "$INPUT matrix does not exist" 
	exit
fi





awk 'BEGIN{getline;getline; min=9999999999999;max=-99; bin=0.+"'"$BINSIZE"'"; s=0.+"'"$START"'"; e=0.+"'"$END"'";rv="'"$RV"'"/2.}{
	split($1,a,"-"); 
	split(a[1],b,":");
	start=b[2]+rv

    split($2,c,"-"); 
    split(c[1],d,":");
    end=d[2]+rv

	matrix[start,end]=$3
	bool[start,end]++
	if(start>max) max=start
	if(start<min) min=start
	print start, end > "aaa"
		
}END{
	for(i=min;i<=max;i+=bin) if(i>s && i<e) printf "%s ", "chrX:"i
	printf "\n"
	for(i=min;i<=max;i+=bin){
		if(i>s && i<e) printf "%s ","chrX:"i
		for(j=min;j<=max;j+=bin){
			if(i>s && i<e && j>s && j<e){ 
				if(bool[i,j]>0) printf "%s ", matrix[i,j];
				else printf "%s ", "0"  
			}
		}
		if(i>s && i<e)	printf "\n"
	}

}' $INPUT > $INPUT.matrix


