{
for(i=1; i <=NF; i++) {
	split($i, a, "=");
	j=a[2];
	printf a[2]; printf " ";
}
printf "\n";
}
