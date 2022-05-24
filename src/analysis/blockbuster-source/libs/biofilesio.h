#ifndef __READBIOFILESIO_H  
#define __READBIOFILESIO_H


/*---------------------*/
//   STRUCTURES  
/*---------------------*/
typedef struct {
	int start;
	int end;
} exon_t;
/*---------------------*/
typedef struct {
	char *id;
	char *chrom;
	int start;
	int end;
	char *strand;
	int CDSstart;
	int CDSend;
	int exonCount;
	char *exonSizes;
	char *exonStarts;
	exon_t *exons;
	int noofexons;
	double expression;
} gene_t;
/*---------------------*/
typedef struct {
	gene_t *genes;
	int noofgenes;
} geneset_t;
/*---------------------*/
typedef struct {
	char *id;
	char *chrom;
	int start;
	int end;
	char *strand;
	int noofreads;
	int mappingFreq;
	double expression;
} tag_t;
/*---------------------*/
typedef struct {
	tag_t *tags;
	int noofallreads;
	int noofreads;
	int nooftags;
	double expression;
} tagset_t;
/*---------------------*/




/*---------------------*/
//   FUNCTIONS
/*---------------------*/
extern tagset_t *read_segemehl_file(char*);
extern geneset_t *read_bed_file(char*);
extern char *append(const char*, const char);

#endif