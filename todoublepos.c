/* 
   todoublepos: convert positions in SDF files from float to double;
       copy other file contents to output unchanged

   DONE: 1) get list of names in SDF file
   DONE: 2) read all scalars
   DONE: 3) read all structure members, convert float to doubles as appropriate
   DONE: 4) calculate offsets, write adjusted members into buffer
   DONE: 5) write scalars and buffer to file
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <SDF.h>

typedef enum SDF_type_enum SDF_type_t;

typedef union {
    int i;
    float f;
    double d;
} datum_t;

static void initargs(int argc, char *argv[], SDF **sdfp, FILE **fp);
static void writeinit(FILE *fp);
static void writescalars(SDF *sdfp, FILE *fp);
static void writestructs(SDF *sdfp, FILE *fp);

int main(int argc, char *argv[])
{
    SDF *sdfp = NULL;
    FILE *fp = NULL;

    initargs(argc, argv, &sdfp, &fp);

    writeinit(fp);
    writescalars(sdfp, fp);
    writestructs(sdfp, fp);

    fclose(fp);
    SDFclose(sdfp);

    return 0;
}

static void initargs(int argc, char *argv[], SDF **sdfp, FILE **fp)
{
    char input;

    if (argc != 3) {
	fprintf(stderr, "Usage: %s SDFfile outfile\n", argv[0]);
	exit(1);
    }

    *sdfp = SDFopen(NULL, argv[1]);
    if (*sdfp == NULL) {
	fprintf(stderr, "%s: %s: %s\n", argv[0], argv[1], SDFerrstring);
	exit(2);
    }

    if (access(argv[2], F_OK) == 0) {
        fprintf(stderr, "%s: %s exists; overwrite (y/n)? ", argv[0],
		argv[2]);
        input = getc(stdin);
        if ((input != 'y') && (input != 'Y'))
            exit(3);
    }
    *fp = fopen(argv[2], "w");
    if (*fp == NULL) {
	fprintf(stderr, "%s: %s\n", argv[2], strerror(errno));
	exit(errno);
    }
}

static void writeinit(FILE *fp)
{
    fprintf(fp, "# SDF\n");
    fprintf(fp, "parameter byteorder = %#x;\n", SDFcpubyteorder());
}

static void writescalars(SDF *sdfp, FILE *fp)
{
    int i, nvecs, nrecs;
    char **vecs;
    SDF_type_t type;
    datum_t datum;

    nvecs = SDFnvecs(sdfp);
    vecs = SDFvecnames(sdfp);

    for (i = 0; i < nvecs; ++i) {
	nrecs = SDFnrecs(vecs[i], sdfp);
	if (nrecs != 1) continue;

	type = SDFtype(vecs[i], sdfp);
	/* Doesn't handle string */
	/* Also, SDFgetfloat was happy to read a double scalar and
	   convert it for me; that probably works backwards too.  I
	   don't think there's an equivalent for the SDFrdvecs family
	   though, so "read; convert type; write" is the general
	   path. */
	switch (type) {
	case SDF_INT:
	    SDFgetint(sdfp, vecs[i], &(datum.i));
	    break;
	case SDF_FLOAT:
	    SDFgetfloat(sdfp, vecs[i], &(datum.f));
	    break;
	case SDF_DOUBLE:
	    SDFgetdouble(sdfp, vecs[i], &(datum.d));
	    break;
	default:
	    fprintf(stderr, "%s: type not supported\n", vecs[i]);
/* 	    exit(-1); */
	}

	switch (type) {
	case SDF_INT:
	    fprintf(fp, "int %s = %d;\n", vecs[i], datum.i);
	    break;
	case SDF_FLOAT:
	    fprintf(fp, "float %s = %.7g;\n", vecs[i], datum.f);
	    break;
	case SDF_DOUBLE:
	    fprintf(fp, "double %s = %.16g;\n", vecs[i], datum.d);
	    break;
	default:
	    fprintf(stderr, "%s: type not supported\n", vecs[i]);
/* 	    exit(-1); */
	}
    }
}

static void writestructs(SDF *sdfp, FILE *fp)
{
    int i, j, nvecs, nmembers, nrecs;
    char **vecs, **members;
    SDF_type_t *types;
    size_t stride = 0, outstride = 0;
    void *btab, *outbtab;
    void **addrs;
    int *strides, *nobjs, *starts, *todouble, *inoffsets, *outoffsets;
    double tmp;

    nvecs = SDFnvecs(sdfp);
    vecs = SDFvecnames(sdfp);

    /* Count structure members */
    for (i = 0, nmembers = 0; i < nvecs; ++i) {
	nrecs = SDFnrecs(vecs[i], sdfp);
	if (nrecs > 1) ++nmembers;
    }

    members = (char **)malloc(nmembers * sizeof(char *));
    addrs = (void **)malloc(nmembers * sizeof(void *));
    strides = (int *)malloc(nmembers * sizeof(int));
    nobjs = (int *)malloc(nmembers * sizeof(int));
    starts = (int *)malloc(nmembers * sizeof(int));
    types = (SDF_type_t *)malloc(nmembers * sizeof(SDF_type_t));
    todouble = (int *)malloc(nmembers * sizeof(int));
    inoffsets = (int *)malloc(nmembers * sizeof(int));
    outoffsets = (int *)malloc(nmembers * sizeof(int));

    for (i = 0, nmembers = 0, stride = 0; i < nvecs; ++i) {
	nrecs = SDFnrecs(vecs[i], sdfp);
	if (nrecs > 1) {
	    members[nmembers] = vecs[i];
	    nobjs[nmembers] = nrecs;
	    starts[nmembers] = 0;  /* Not correct in parallel; use
				      NobjInitial */
	    types[nmembers] = SDFtype(members[nmembers], sdfp);
	    inoffsets[nmembers] = stride;
	    stride += SDFtype_sizes[types[nmembers]];

	    if ( ( (strncmp(vecs[i], "x", strlen(vecs[i])) == 0) ||
		   (strncmp(vecs[i], "y", strlen(vecs[i])) == 0) ||
		   (strncmp(vecs[i], "z", strlen(vecs[i])) == 0) ) &&
		 (types[nmembers] == SDF_FLOAT) ) {
		todouble[nmembers] = 1;
	    }
	    else {
		todouble[nmembers] = 0;
	    }

	    ++nmembers;
	}
    }

    btab = (void *)malloc(nobjs[0] * stride);

    for (i = 0; i < nmembers; ++i) {
	addrs[i] = (char *)btab + inoffsets[i];
	strides[i] = stride;
    }

    SDFseekrdvecsarr(sdfp, nmembers, members, starts, nobjs, addrs,
		     strides);

    for (i = 0, outstride = 0; i < nmembers; ++i) {
	outoffsets[i] = outstride;
	outstride += SDFtype_sizes[ todouble[i] ? SDF_DOUBLE : types[i] ];
    }

    outbtab = (void *)malloc(nobjs[0] * outstride);

    for (j = 0; j < nobjs[0]; ++j) {
	for (i = 0; i < nmembers; ++i) {
	    if (todouble[i]) {
		tmp = (double)(*(float *)(btab + j*stride + inoffsets[i]));
		memcpy(outbtab + j*outstride + outoffsets[i],
		       &tmp, SDFtype_sizes[SDF_DOUBLE]);
	    }
	    else {
		memcpy(outbtab + j*outstride + outoffsets[i],
		       btab + j*stride + inoffsets[i],
		       SDFtype_sizes[types[i]]);
	    }
	}
    }

    fprintf(fp, "struct {\n");
    for (i = 0; i < nmembers; ++i) {
	switch (types[i]) {
	case SDF_INT:
	    fprintf(fp, "\tint %s;\n", members[i]);
	    break;
	case SDF_FLOAT:
	    if (todouble[i]) {
		fprintf(fp, "\tdouble %s;\n", members[i]);
	    }
	    else {
		fprintf(fp, "\tfloat %s;\n", members[i]);
	    }
	    break;
	case SDF_DOUBLE:
		fprintf(fp, "\tdouble %s;\n", members[i]);
	    break;
	default:
	    fprintf(stderr, "%s: type not supported\n", members[i]);
	    exit(-1);
	}
    }
    fprintf(fp, "}[%d];\n", nobjs[0]);
    fprintf(fp, "#\n");
    fprintf(fp, "# SDF-EOH\n");

    fwrite(outbtab, outstride, nobjs[0], fp);

    free(members);
    free(addrs);
    free(strides);
    free(nobjs);
    free(starts);
    free(types);
    free(todouble);
    free(inoffsets);
    free(outoffsets);

    free(btab);
    free(outbtab);
}
