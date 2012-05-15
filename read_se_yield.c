#include<stdio.h>
#include<stdlib.h>
#include<string.h>

int main(int argc, char *argv[])
{
	FILE *fp=NULL;
	int i, j, nlines=0;
	char temp[21], line[101];
	int *partid;

	if(argc != 2) {
		printf("%s <particle-file>\n", argv[0]);
		return 0;
	}

	fp = fopen(argv[1], "r");
	if(fp == NULL) {
		printf("error opening file!\n");
		return 0;
	}


	while( !feof(fp) ) {
		fgets(line, 100, fp);
		nlines++;
	}

	rewind(fp);

	partid = (int *)malloc( nlines*sizeof(int) );
	fgets(line, 100, fp);

	for(i = 0; i < nlines; i++) {
		fscanf(fp, "%s %d %*s %*s %*e", temp, &partid[i]);
		if( strncmp(temp, "nn", strlen(temp)) == 0)
			break;
	}

	printf("particles: %d\n",i);
	for( j = 0; j< i; j++)
		printf("%d\n", partid[j]);

}
