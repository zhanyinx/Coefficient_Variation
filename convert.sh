#!/bin/bash

bin=6000
s=99000000
e=103350000
#running average
rv=0
awk 'BEGIN{getline;getline; min=9999999999999;max=-99; bin=0.+"'"$bin"'"; s=0.+"'"$s"'"; e=0.+"'"$e"'";rv="'"$rv"'"/2.}{
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

}' $1 > $1.matrix


