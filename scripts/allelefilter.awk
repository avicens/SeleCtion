#Usage: awk -f allelefilter.awk [file] | awk '{if ($NF>=2) {printf $0}}' > outfile


 {sum=0;
if (NR==1) {printf $0; printf "\t";}

if(NR>1) {

for(i=1; i<=9; i++) { 
	printf $i; printf "\t";
} 

 for(i=10; i<= NF; i++) {
	printf $i; printf "\t";
	 
	split($i, subfield, ":"); 
	j=subfield[2];
	split(j, subsub, ",");
	k=subsub[2];
	if (k>=2) {
	sum++;}
}
}

if (NR==1) { printf "count"; }
else {printf sum;}
	print "";

}
