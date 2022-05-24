#include "biofilesio.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

char *append(const char *oldstring, const char c)
{
    int result;
    char *newstring;
    result = asprintf(&newstring, "%s%c", oldstring, c);
    if (result == -1) newstring = NULL;
    return newstring;
}

tagset_t *read_segemehl_file(char *filename)
/* READ FILE WITH A SET OF MAPPED READS */
{	
	// create new set of clusters
	tagset_t *set;
	set = calloc(1, sizeof(tagset_t));
	set->tags = NULL;
	set->noofallreads = 0;
	set->noofreads = 0;
	set->nooftags = 0;
	set->expression = 0;
	
	FILE *f;
	f = fopen(filename, "r");
	if(!f)
	{
		printf("cannot open %s\n", filename); exit(0);
	}
	else
	{
		char c; int e = 0; 
		while((c=getc(f)) != EOF) {
			
			char dumchar[50], info[50], strand[5], chrom[50], id[50] = "";
			int  dumint, start, end, mappingFreq, height = 0;
			
			// check for header
			int header = 0;
			
			if(c == '#')
				header = 1;
			ungetc (c,f);
			
			if(header == 0){
				// parse line
				e = fscanf(f, "%s %s %d %d %d %d %d %d %d %s %d %d %s %s %d", dumchar, info, &dumint, &dumint, &dumint, &dumint, &dumint, &dumint, &dumint, strand, &start, &end, chrom, dumchar, &mappingFreq);
				
				// split id|freq
				char *p = strtok(info, "|"); int i = 0;
				while(p != 0)
				{
					if(i == 0) strcpy(id, p);
					if(i == 1) height = atoi(p);
					i++;
					p = strtok(NULL, "|");
				}
				if(height == 0){height = 1;}
				
				//				printf("id: %s\nchrom: %s\tstart: %d\tend: %d\tstrand: %s\nreads: %d\tmappingFreq: %d\texpression: %lf\nall reads: %d\n\n",id,chrom,start,end,strand,height,mappingFreq,((double)height/(double)mappingFreq),set->noofallreads);
				
				// get memory for new tag
				set->nooftags++;
				set->tags = realloc(set->tags, (set->nooftags+1)*sizeof(tag_t));
				
				// write tag
				// id
				set->tags[set->nooftags].id = malloc(strlen(id)+1);
				strcpy(set->tags[set->nooftags].id, id);
				// chrom
				set->tags[set->nooftags].chrom = malloc(strlen(chrom)+1);
				strcpy(set->tags[set->nooftags].chrom, chrom);
				// start
				set->tags[set->nooftags].start = start;
				// end
				set->tags[set->nooftags].end = end;
				// strand
				set->tags[set->nooftags].strand = malloc(strlen(strand)+1);
				strcpy(set->tags[set->nooftags].strand, strand);
				// reads
				set->tags[set->nooftags].noofreads = height;
				// mapping frequency
				set->tags[set->nooftags].mappingFreq = mappingFreq;
				// expression
				set->tags[set->nooftags].expression = ((double)height/(double)mappingFreq);	
				
				// write set
				set->noofreads += height;
				set->expression += ((double)height/(double)mappingFreq);
			}
			else if(header == 1){
				char type[50] = "";
				char file[50] = "";
				char freqS[50] = "";
				char dum[50] = "";
				e = fscanf(f, "%s %s %s %s", type, file, freqS, dum);
				if(strcmp(type, "#query:") == 0){
					char *p = strtok(freqS, "("); int i = 0;
					while(p != 0)
					{
						set->noofallreads = atoi(p);
						i++;
						p = strtok(NULL, "(");
					}
				}
				else{
					char line[200];
					e = fscanf(f, "%[^\n]", line);
					
				}
				fgetc(f);
			}
		}
	}
	fclose(f);
	return set;
}

geneset_t *read_bed_file(char *filename)
/* READ FILE WITH A SET OF ANNOTATIONS */
{	
	// create new set of clusters
	geneset_t *set;	
	set = calloc(1, sizeof(gene_t));
	set->genes = NULL;
	set->noofgenes = 0;
	
	FILE *f;
	f = fopen(filename, "r");
	if(!f)
	{
		printf("cannot open %s\n", filename); exit(0);
	}
	else
	{
		char c; int header = 0;  
		
		char *line = "";
		
		while((c=getc(f)) != EOF) {
			if(c == '\n'){
				//				printf("line: %s\n",line);
				
				// get memory for new gene
				set->noofgenes++;
				set->genes = realloc(set->genes, (set->noofgenes+1)*sizeof(gene_t));
				set->genes[set->noofgenes].noofexons = 0;
				set->genes[set->noofgenes].exons = 0;
				set->genes[set->noofgenes].CDSstart = 0;
				set->genes[set->noofgenes].CDSend = 0;
				set->genes[set->noofgenes].exonCount = 0;
				
				char *p = strtok(line, "\t"); int i = 0;
				while(p != 0)
				{
					if(i == 0){ 
						set->genes[set->noofgenes].chrom = malloc(strlen(p)+1); 
						strcpy(set->genes[set->noofgenes].chrom, p);
					}
					if(i == 1) 
						set->genes[set->noofgenes].start = atoi(p);
					if(i == 2) 
						set->genes[set->noofgenes].end = atoi(p);
					if(i == 3){ 
						set->genes[set->noofgenes].id = malloc(strlen(p)+1); 
						strcpy(set->genes[set->noofgenes].id, p);
					}
					if(i == 5){
						set->genes[set->noofgenes].strand = malloc(strlen(p)+1); 
						strcpy(set->genes[set->noofgenes].strand, p);
					}
					if(i == 6) 
						set->genes[set->noofgenes].CDSstart = atoi(p);
					if(i == 7) 
						set->genes[set->noofgenes].CDSend = atoi(p);
					if(i == 9) 
						set->genes[set->noofgenes].exonCount = atoi(p);
					if(i == 10) 
					{
						set->genes[set->noofgenes].exonSizes = malloc(strlen(p)+1); 
						strcpy(set->genes[set->noofgenes].exonSizes, p);
					}
					if(i == 11){
						//						printf("setStart: %s\nsetSize: %s\n",p,set->genes[set->noofgenes].exonSizes);
						set->genes[set->noofgenes].exonStarts = malloc(strlen(p)+1); 
						strcpy(set->genes[set->noofgenes].exonStarts, p);
						
						// write exons
						exon_t *exons;
						exons = calloc(1, sizeof(exon_t));
						char *q = strtok(set->genes[set->noofgenes].exonStarts, ","); int j = 0;
						while(q != 0)
						{
							set->genes[set->noofgenes].exons = realloc(set->genes[set->noofgenes].exons, (j+1)*sizeof(exon_t));
							set->genes[set->noofgenes].exons[j].start = set->genes[set->noofgenes].start + atoi(q);
							j++;
							q = strtok(NULL, ",");
						}
						
						q = strtok(set->genes[set->noofgenes].exonSizes, ","); j = 0;
						while(q != 0)
						{
							set->genes[set->noofgenes].exons[j].end = set->genes[set->noofgenes].exons[j].start + atoi(q);
							//							printf("exon %d\t%d-%d\n",j,set->genes[set->noofgenes].exons[j].start,set->genes[set->noofgenes].exons[j].end);
							j++;
							q = strtok(NULL, ",");
						}
					}						
					i++;
					p = strtok(NULL, "\t");
				}
				
				free(line); line = "";
			}
			else
				line = append(line, c); 
		}
	}
	fclose(f);
	return set;
}

