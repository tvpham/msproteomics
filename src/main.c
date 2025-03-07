/* This program reads a DiaNN ouput report file
 * the contains a column containing a modified
 * sequence column and outputs a simiar file with
 * processed fields derived from this field.
 * Also FASTA files need to be specified in order
 * to expand the modifiaction locations within all
 * proteins that contain the backbone of modified
 * peptides as tryptic peptides.
 *
 * Author: Alex Henneman
 * Date    : Wed 12 Jul 2023 07:27:34 PM CEST
 * Modified: Thu Oct 31 10:42:10 AM CET 2024
 *
 * Todos
 * =====
 *
 * 0. Get rid of ceil calls, so compilation
 *    is simpler.
 * 1. If no FASTA files are loaded, there's nothing this 
 *    program can do. It should abort then.
 * 2. If it is not able find the required columns modpep
 *    and fragment.intensities it should abort.
 * 3. Add a fag to skip fragment intensity expansion.
 * 4. Modify program for multifile processing.
 *
 * --------------- */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <libgen.h>

#include <strings.h>
#include <math.h>

/* --- */
#define VERSION_STRING "0.4.0"
/* This one is printed in the help blurp */
#define STAMP_DATE "Fri Mar 7 09:53:18 PM CEST 2025"

#define ERROR_FAILED -1

#define VERBOSE_QUIET 0
#define VERBOSE_LEVEL1 1
/* This is the header for the extra columns that
 * will be written in the function below. */
#define EXTRA_OUTPUT_COLNAMES "\tProtein.Sites\tFasta.Files\tProteins.Fasta"\
    "\tPeptide.Backbone\tPeptide.Offsets\tNr.Mods\tNr.Phospho"\
    "\tPhospho.Site.Specs\tFragment.Rel.Id\tFragment.Intensity"
/* THese are the columns when no fragment intensity
 * expansion takes place. */
#define EXTRA_OUTPUT_COLNAMES_NOEXP "\tProtein.Sites\tFasta.Files\tProteins.Fasta"\
    "\tPeptide.Backbone\tPeptide.Offsets\tNr.Mods\tNr.Phospho"\
    "\tPhospho.Site.Specs\tFragment.Intensities"
#define FRAGMENT_LIST_LEN_INIT 15

#define INIT_MAX_FASTAS   10
#define INIT_DELTA_FASTAS 10
#define DEFAULT_OUTFILE "output.tsv"
#define DIGEST_PEPTIDE_LIST_LEN_INIT 100
#define SEQUENCE_ALPHABET "ACDEFGHIJKLMNPQRSTVWYUX"

#define DEFAULT_PEPTIDE_LIST_INIT_LEN 100
#define DEFAULT_PEPTIDE_LIST_DELTA_LEN 50

#define PEPTIDE_LIST_LENGTH_INIT 1000

#define COLNAME_DEFAULT_FILE_NAME                "File.Name"
#define COLNAME_DEFAULT_RUN                      "Run"
#define COLNAME_DEFAULT_PROTEIN_GROUP            "Protein.Group"
#define COLNAME_DEFAULT_PROTEIN_IDS              "Protein.Ids"
#define COLNAME_DEFAULT_GENES                    "Genes"
#define COLNAME_DEFAULT_MODIFIED_SEQUENCE        "Modified.Sequence"
#define COLNAME_DEFAULT_STRIPPED_SEQUENCE        "Stripped.Sequence"
#define COLNAME_DEFAULT_PRECURSOR_ID             "Precursor.Id"
#define COLNAME_DEFAULT_Q_VALUE                  "Q.Value"
#define COLNAME_DEFAULT_GLOBAL_Q_VALUE           "Global.Q.Value"
#define COLNAME_DEFAULT_PROTEOTYPIC              "Proteotypic"
#define COLNAME_DEFAULT_MS1_AREA                 "Ms1.Area"
#define COLNAME_DEFAULT_FRAGMENT_QUANT_CORRECTED "Fragment.Quant.Corrected"
#define COLNAME_DEFAULT_PTM_Q_VALUE              "PTM.Q.Value"
#define COLNAME_DEFAULT_PTM_SITE_CONFIDENCE      "PTM.Site.Confidence"

#define TARGETS_LEN_INIT 20 
#define TARGETS_LEN_DELTA 10

#define DEFAULT_TAB     '\t'
#define DEFAULT_NEWLINE '\n'
#define DEFAULT_NULL    '\0'

#define DEFAULT_NODE_FASTA_LIST_DELTA 5
#define DEFAULT_NODE_PROTEIN_LIST_DELTA 5

#define DEFAULT_NODE_FASTA_LIST_INIT 5
#define DEFAULT_NODE_PROTEIN_LIST_INIT 5

#define PROC_MODE_STATIC_REPEAT 1
#define PROC_MODE_MODPEP_QUERY  2
#define PROC_MODE_FRAG_EXPAND   3

#define HEADER_LIST_LEN_INIT 20
#define HEADER_LIST_LEN_DELTA 10

#define UNFOUND_INDEX -1

#define DEFAULT_MIN_PEPTIDE_LEN 7
#define DEFAULT_MAX_PEPTIDE_LEN 50
#define DEFAULT_MAX_MISSED_CLEAVAGES 2

#define PIPE_IO_FILENAME "-"

#define SITESPEC_LEN_INIT 1000
#define SITESPEC_LEN_DELTA 1000
#define PROTEINS_FOUND_BUFFER_LEN_INIT 1000
#define PROTEINS_FOUND_BUFFER_LEN_DELTA 1000
#define MOD_VEC_LEN_INIT 100
#define MOD_VEC_LEN_DELTA 100
#define MOD_OFFSET_STRING_LEN_INIT 100
#define MOD_OFFSET_STRING_LEN_DELTA 100

#define APPEND_TO_UINT_VECTOR_T_DELTA 100
/* --- */

typedef struct {

    unsigned int   len;
    unsigned int   maxlen;
    char         **string;

} string_list_t;

/* --- */

/* Here we put all program parameters
 * and settings. We also put in here 
 * obejects that are used and reused during
 * program execution. */

typedef struct Params {

    int    n_fastas;
    int    max_fastas;
    char **fasta;

    char  *infile;
    char  *outfile;

    string_list_t *targets;

    char *expand_name;
    char *modpep_name;

    char column_separator;
    char line_terminator;
    char null_char;

    int verbose_level;

    int minlen;
    int maxlen;
    int max_missed_cleavages;

    /* The ones above are specified
     * by user for the desired operation
     * of this program.
     *
     * The ones here below are more like 
     * reading implementation
     * helpers... */
    char         *buffer;
    unsigned int  buffer_len;

    /* These are derived from 
     * the targets, after processing
     * an input file. */
    int           *column_index;
    int           *proc_mode;

/* Whether fragment expansion is on/off */
    int fragment_intensity_expansion;

    int  expand_indx;
    int  modpep_indx;

    int  ncols;

} params_t;

/* ------ */

/* This is the list of peptides returned after
 * digestion together with their corresponding
 * offset in the parent protein. */

typedef struct {

    unsigned int   len;
    unsigned int   maxlen;
    char         **peptide;
    unsigned int  *offset;

} peptide_list_t;

/* ------ */

/* This object carries infor from
 * FASTA file to node, so there's
 * only a single FASTA file and
 * protein. */

typedef struct {

    int   fasta_id;
    char *protein;
    int   offset;

} message_t;

/* ------ */

/* This is where we store the data when we
 * reach a specific node in the char_tree. */

typedef struct {

    int           *fasta_id;
    int            nfastas;
    int            maxfastas;

    unsigned int   maxprots;
    unsigned int   nprots;
    char         **protein;
    int           *offset;

} node_data_t;

/* ------ */

/* In this object we store all char_tree
 * query results, all the way down to the
 * strings that are going to be printed in
 * the output for the current input line. */

typedef struct {

    /* This one is non NULL after a
     * hit in char_tree_query */
    node_data_t *node_data;

    /* These we set at a higher level after a
     * succesfull query to what we want to 
     * print for a specific modpep. */
    char *sitespec;
    char *phosphospec;

    char *fastas;
    char *proteins;
    char *backbone;

    char *modoffs;

    int nmods;
    int nphospho;

    int found_flag;

} ptm_locations_t;

/* ------ */

/* This one is to hold pointers to the strings
 * contsining the contents of the reproduced filed
 * values in the current line. */

typedef struct {

    int    len;
    int    maxlen;

    char **string;

} static_fields_t;

/* ------ */

typedef struct {

    unsigned int  len;
    char         *alphabet;

} char_tree_definition_t;

/* ------ */

typedef struct CharTreeNode {

    char_tree_definition_t  *def;
    struct CharTreeNode    **branch;
    void                    *payload;

} char_tree_node_t;

/* ------ */

typedef struct UintVector {

    unsigned int  len;
    unsigned int  maxlen;
    unsigned int *element;

} uint_vector_t;

/* ------ */

typedef struct ProteinSequence {

    char                   *sequence;
    char                   *name;

} protein_sequence_t;

/* ------ */

typedef struct ProteinList {

    unsigned int            len;
    unsigned int            maxlen;

    /* Maybe change this into 
     * protein entry or something */

    struct ProteinSequence *protein;

} protein_list_t;

/* ------ */

/* === FUNCTIONS === */

/* This function allocates a new
 * protein_list_t object.
 *
 * TODO add code for minimal size
 * requirement. Enlarge object if
 * pointer is not zero and length is
 * smaller than requested. */

int allocate_protein_list_t (

        protein_list_t **plist_pp,
        unsigned int     len

){
    int             retval=0;
    protein_list_t *plist;
    size_t          size;

    if ((plist = calloc(1,sizeof(protein_list_t))) == NULL){
        fprintf(stderr,"Failed allocating header in allocate_protein_list_t.\n");
        retval = ERROR_FAILED;
    }
    else {
        plist->protein = calloc(len,sizeof(protein_sequence_t));
#if __DEBUG
        fprintf(stderr,"The protein_list_t data appended to %p is located at %p and has "\
                "size of %u bytes.\n",plist,plist->protein,size);
#endif

        if ((plist->protein)==NULL){
            fprintf(stderr,"Failed allocating data in allocate_protein_list_t.\n");
            retval = ERROR_FAILED;
        }
        else {

            /* Everything went ok. So
             * we finish up now. */
            plist->len    = 0;
            plist->maxlen = len;
            /* Here we go */
            *plist_pp = plist;
        }
    }

    return retval;
}

/* ---------- */

/* This function enlarges the
 * protein list array by the
 * amount specified by delta. */

int expand_protein_list (

        protein_list_t *list,
        unsigned int    delta

){
    protein_sequence_t *new_array;
    unsigned int        new_len;
    size_t              new_size;
    int                 retval;

    retval = 0;

    new_len   = list->maxlen + delta;
    new_size  = new_len*sizeof(protein_sequence_t);
    new_array = realloc(list->protein,new_size);

    if (new_array != NULL){

        list->maxlen  = new_len;
        list->protein = new_array;
    }
    else retval = ERROR_FAILED;

    return retval;
}

/* ---------- */

/* The main purpose of this function is to
 * free up all memory allocated for the peptide
 * list passed as an argument. It will set its
 * return value different to 0 if any irregularities
 * are encountered. It only free up memory for
 * the proteins from index 0 to index list->len
 * and NOT maxlen!.   */

int free_protein_list_t (

        protein_list_t *plist

){
    unsigned int i;
    int retval;

    retval = 0;

    if (plist==NULL){
        fprintf(stderr,"NULL pointer to protein list in free_protein_list_t.\n");
        retval = ERROR_FAILED;
    }
    else {
        for (i=0;i<(plist->len);i++){
            if (plist->protein[i].name != NULL)
                free(plist->protein[i].name);
            else {
                fprintf(stderr,"NULL name string in free_protein_list_t.\n");
                retval = ERROR_FAILED;
            }
            if (plist->protein[i].sequence != NULL)
                free(plist->protein[i].sequence);
            else {
                fprintf(stderr,"NULL sequence string in free_protein_list_t.\n");
                retval = ERROR_FAILED;
            }

        }
    }

    return(retval);

}

/* ---------- */

/* This fucntion reads a full line expanding if
 * necessary (unbounded!) the line buffer if
 * necessary. If the linebuffer pointer is set
 * to NULL prior to calling this function a new
 * buffer is allocated with length equal to
 * the amount LINE_BUFFER_LENGTH_INIT. If the
 * buffer needs to expand it is expanded by the
 * amount LINE_BUFFER_LENGTH_DELTA.  */
#define LINE_BUFFER_LENGTH_INIT  10000000
#define LINE_BUFFER_LENGTH_DELTA 10000000

int fgets_full_line (

        FILE          *fd,
        char         **linebuffer_p,
        unsigned int  *buflen_p
){
    char         *linebuffer;
    unsigned int  buflen;
    int           retval=0;
    off_t         start_offset;

    if ((*linebuffer_p==NULL)||(*buflen_p==0)){
        buflen = LINE_BUFFER_LENGTH_INIT;
        if ((linebuffer=calloc(buflen,sizeof(char)))==NULL){
            fprintf(stderr,"Failed allocating line buffer in "\
                    "fgets_full_line.\n");
            retval = ERROR_FAILED;
        }
        else {
            *linebuffer_p = linebuffer;
            *buflen_p     = buflen;
        }
    }
    else {
        linebuffer = *linebuffer_p;
        buflen     = *buflen_p;
    }
    if(retval)return(retval);

    start_offset = ftello(fd);
    if (fgets(linebuffer,buflen,fd)==NULL){
        /* This is a quiet error because we have
         * reached the end of the file. */
        retval = ERROR_FAILED;
        return(retval);
    }
    while (strlen(linebuffer) >= (buflen-1)){
        buflen += LINE_BUFFER_LENGTH_DELTA;
        if ((linebuffer=realloc(linebuffer,buflen))==NULL){
            fprintf(stderr,"Failed expanding line buffer "\
                    "in fgets_full_line.\n");
            retval = ERROR_FAILED;
            break;
        }
        else {
            *linebuffer_p = linebuffer;
            *buflen_p   = buflen;
        }
        fseeko(fd,start_offset,SEEK_SET);
        if (fgets(linebuffer,buflen,fd)==NULL){
            fprintf(stderr,"Failed getting line in the "\
                    "expansion loop of fgets_full_line.\n");
            retval = ERROR_FAILED;
            break;
        }
    }

    return(retval);
}

/* ------- */

#define PROTEIN_LIST_LEN_INIT 5
#define PROTEIN_LIST_LENGTH_DELTA 100
#define LINE_BUFFER_LENGTH 50000
#define SEQUENCE_BUFFER_LENGTH_INIT  2000
#define SEQUENCE_BUFFER_LENGTH_DELTA 2000

/* This is the function that actually reads the
 * fasta file. */

int parse_fasta_file (

        FILE            *protein_fd,
        protein_list_t **plist_pp

){
    protein_list_t *plist_p;
    unsigned int    max_fill,sequence_buffer_fill,
                    last_char_index,line_buffer_len;
    char           *sequence_buffer,*line_buffer;
    int             retval,in_protein;

    retval = 0;

    /* This is the length of an expandable buffer
     * for gathering protein sequences. */
    max_fill = SEQUENCE_BUFFER_LENGTH_INIT;
    sequence_buffer_fill = 0;

    /* Preliminaries before start parsing the
     * fasta file. */

    if ( (sequence_buffer = calloc((max_fill+1),sizeof(char)))==NULL ){
        fprintf(stderr,"Failed allocating sequence buffer in parse_fasta_file.\n");
        return ERROR_FAILED;
    }
    sequence_buffer[0] = '\0';

    if ( allocate_protein_list_t(&plist_p, PROTEIN_LIST_LEN_INIT) ) {
        fprintf(stderr,"Failed allocating  protein list in parse_fasta_file.\n");
        return ERROR_FAILED;
    }

    /* Now we start parsing the fasta file. */

    in_protein = 0;

    line_buffer_len = 0;
    line_buffer     = NULL;

    int line_error=0;

    while (!fgets_full_line(protein_fd,&line_buffer,&line_buffer_len)){
        if(line_error){
            fprintf(stderr,"Error: Detected intermediate line without terminating newline.\n");
            retval = ERROR_FAILED;
            break;
        }
        /* Here we check whether we get the whole line or not.
         * If we don't we issue a warning and set retval.*/
        last_char_index = strlen(line_buffer)-1;
        if( line_buffer[last_char_index] != '\n' ){
            line_error = 1;
        }
        else {
            /* And we take off the newline. */
            line_buffer[last_char_index] = '\0';
        }

        /* Decide what kind of line this is. */
        if (line_buffer[0] == '>' ){
            if (in_protein){

                /* Time to close old protein. Store the sequence
                 * and update protein counter. */
                if ((plist_p->protein[plist_p->len].sequence =
                          strdup(sequence_buffer)) == NULL){
                    fprintf(stderr,"Failed duplicating sequence string in parse_fasta_file.\n");
                    retval = ERROR_FAILED;
                    break;
                }
                /* See whether protein counter is a valid index for
                 * next available protein entry. */
                if (plist_p->len >= (plist_p->maxlen-1)) {
                    if ( expand_protein_list(plist_p,
                                    PROTEIN_LIST_LENGTH_DELTA) ){
                        fprintf(stderr,"Failed expanding protein list in parse_fasta_file.\n");
                        retval = ERROR_FAILED;
                        break;
                    }
                }
                plist_p->len++;

            }
            /* Start a new one */
            if ((plist_p->protein[plist_p->len].name =
                        strdup(line_buffer))==NULL){
                 fprintf(stderr,"Failed duplicating protein name in parse_fasta_file.\n");
                 retval = ERROR_FAILED;
                 break;
            }
            in_protein           = 1;
            sequence_buffer_fill = 0;
            sequence_buffer[0]   = '\0';
        }
        else {
            if (!in_protein) {
                printf("Error: non-header line found when not in protein in parse_fasta_file.\n");
                retval = ERROR_FAILED;
                break;
            }
            /* The only thing we can do here is add the new line
             * to the sequence buffer. First see if it fits. */
            sequence_buffer_fill += strlen(line_buffer);
            /* Set equal or garter just to see if it fixes
             * anything */
            if (sequence_buffer_fill >= max_fill){
                max_fill = sequence_buffer_fill;
                if ((sequence_buffer = realloc(sequence_buffer,
                                        (max_fill+1)*sizeof(char)))==NULL){
                    fprintf(stderr,"Failed expanding sequence buffer in parse_fasta_file.\n");
                    retval = ERROR_FAILED;
                    break;
                }
            }
            strncat(sequence_buffer,line_buffer,LINE_BUFFER_LENGTH);
            sequence_buffer[sequence_buffer_fill] = '\0';
        }
    }
    /* Handle last protein is we were handling any. */
    if (in_protein){
        /* Handle the sequence. */
        if ((plist_p->protein[plist_p->len].sequence =
                  strdup(sequence_buffer)) == NULL){
            fprintf(stderr,"Failed storing protein sequence in parse_fasta_file.\n");
            retval = ERROR_FAILED;
        }
        /* and update count. */
        plist_p->len++;

    }

    *plist_pp = plist_p;

    return retval;
}

/* ---------- */

int read_fasta_file (

        char            *filename,
        protein_list_t **plist_pp
){
    int retval=0;
    FILE   *protein_fd;

    if ( (protein_fd=fopen(filename,"r"))==NULL ) {
        fprintf(stderr,"Failed opening fasta file %s.\n",filename);
        return ERROR_FAILED;
    }
    else {
        if (parse_fasta_file(protein_fd,plist_pp)){
            fprintf(stderr,"Failed to parse fasta file "\
                    "in read_fasta_file.\n");
            retval = ERROR_FAILED;
        }
        fclose(protein_fd);
    }
    return retval;
}

/* ---------- */

int allocate_uint_vector_t (
        
        unsigned int    len,
        uint_vector_t **vec_pp
        
){
    int            retval = 0;
    uint_vector_t *vec_p;

    if ((vec_p=calloc(1,sizeof(uint_vector_t)))==NULL){
        fprintf(stderr,"Failed allocating header in "\
                "allocate_uint_vector_t.\n");
        retval = ERROR_FAILED;
    }
    else {
        if ((vec_p->element=calloc(len,sizeof(unsigned int)))==NULL){
            fprintf(stderr,"Failed allocating data in "\
                    "allocate_uint_vector_t.\n");
            retval = ERROR_FAILED;
        }
        else {
            vec_p->len    = 0;
            vec_p->maxlen = len;
            *vec_pp       = vec_p;
        }
    }

    return retval;
}

/* ------------- */

int expand_uint_vector_t (

        uint_vector_t *vec,
        unsigned int   delta

){
    int retval=0;

    if (vec==NULL){
        fprintf(stderr,"NULL header pointer in expand_uint_vector_t.\n");
        retval = ERROR_FAILED;
    }
    else {
        vec->maxlen += delta;
        if ((vec->element=realloc(vec->element,
                   (vec->maxlen)*sizeof(unsigned int)))==NULL){
            fprintf(stderr,"Failed to expand uint_vector_t data.\n");
            retval = ERROR_FAILED;
        }
    }

    return retval;
}

/* ------------- */

int append_to_uint_vector_t (
        
        unsigned int value,
        uint_vector_t *clist
        
){
    int retval = 0;

    if (clist->len >= clist->maxlen){
        if (expand_uint_vector_t(clist,APPEND_TO_UINT_VECTOR_T_DELTA)){
            fprintf(stderr,"Failed to expand candidate list.\n");
            retval = ERROR_FAILED;
        }
    }
    clist->element[clist->len] = value;
    clist->len++;

    return retval;
}

/* ---------- */
#define DIGEST_INDEX_LIST_LEN_INIT  100
#define DIGEST_INDEX_LIST_LEN_DELTA 20

int list_trypsin_digestion_candidates (

        char            *sequence,
        uint_vector_t **digest_locs_p

){
    uint_vector_t *digest_locs;
    unsigned int   len,i;
    int            retval=0;

    if (allocate_uint_vector_t( DIGEST_INDEX_LIST_LEN_INIT ,
                digest_locs_p)){
        fprintf(stderr,"Failed allocating digestion index "\
                "location vector in list_trypsin_digestion_candidates.\n");
        retval = ERROR_FAILED;
    }
    else {
        digest_locs = *digest_locs_p;
    }

    if (retval) return retval;

    /* We add the first amino acid in the sequence */
    if (append_to_uint_vector_t(0,digest_locs)){
        fprintf(stderr,"Failed appending to candidate "\
                "list in list_trypsin_digestion_candidates.\n");
        retval = ERROR_FAILED;
    }

    if (sequence[0] == 'M'){
        if (append_to_uint_vector_t(1,digest_locs)){
            fprintf(stderr,"Failed appending to candidate list "\
                    "in list_trypsin_digestion_candidates.\n");
            retval = ERROR_FAILED;
        }
    }

    len = strlen(sequence) - 1;

    for (i=0;i<len;i++){
        if (sequence[i] == 'R'){
            if (append_to_uint_vector_t(i+1,digest_locs)){
                fprintf(stderr,"Failed appending to candidate list "\
                        "in list_trypsin_digestion_candidates.\n");
                retval = ERROR_FAILED;
            }
        }
        if (sequence[i] == 'K'){
            if (append_to_uint_vector_t(i+1,digest_locs)){
                fprintf(stderr,"Failed appending to candidate list "\
                        "in list_trypsin_digestion_candidates.\n");
                retval = ERROR_FAILED;
            }
        }
    }
    if (append_to_uint_vector_t(len+1,digest_locs)){
        fprintf(stderr,"Failed appending to candidate list "\
                "in list_trypsin_digestion_candidates.\n");
        retval = ERROR_FAILED;
    }

    return retval;
}

/* -------------- */

#define APPEND_TO_STRING_LIST_DELTA 500

int append_to_string_list_t (

        string_list_t *slist,
        char          *field
){
    int retval=0;

    if (slist->len >= slist->maxlen){
        slist->maxlen += APPEND_TO_STRING_LIST_DELTA;
        if ((slist->string = realloc(slist->string,
                    (slist->maxlen)*sizeof(char *)))==NULL){
            fprintf(stderr,"Failed to expand data in "\
                    "append_to_string_list_t.\n");
            retval = ERROR_FAILED;
        }
    }
    if(retval) return(retval);

    slist->string[slist->len] = field;
    slist->len++;

    return(retval);
}

/* ---------- */

int append_to_string_list_t_dup (

        string_list_t *slist,
        char          *field
){
    char *copy;
    int   retval=0;

    if((copy = strdup(field))==NULL){
        fprintf(stderr,"Failed making copy in "\
                "append_to_string_list_t_dup.\n");
        retval = ERROR_FAILED;
    }
    else {
        if(append_to_string_list_t(slist,copy)){
            fprintf(stderr,"Failed storing copy in "\
                    "append_to_string_list_t_dup.\n");
            retval = ERROR_FAILED;
        }
    }

    return(retval);
}

/* ---------- */

int generate_digested_subsequences (

        char          *sequence,
        uint_vector_t *locs,
        unsigned int   minlen,
        unsigned int   maxlen,
        unsigned int   mclvgs,
        string_list_t *peplist

){
    unsigned int      end_loc,pep_len,last,i,j;
    char             *new1,tmp;
    int               retval=0;

    last      = locs->len - 1;

    for (i=0;i<last;i++){
        for(j=i+1;j<(locs->len);j++){
            pep_len = locs->element[j] - locs->element[i];
            if (pep_len < minlen) continue;
            if (pep_len > maxlen) break;
            if ( (j-i-1) > mclvgs ) break;

            /*
            end_loc           = locs->element[j] + 1;
            */
            end_loc           = locs->element[j];
            tmp               = sequence[end_loc];
            sequence[end_loc] = '\0';

            if ((new1 = strdup((char *)(sequence+locs->element[i])))==NULL){
                fprintf(stderr,"Failed copying peptide for list in "\
                        "generate_digested_subsequences.\n");
                retval = ERROR_FAILED;
            }

            sequence[end_loc] = tmp;

            if (append_to_string_list_t(peplist,new1)){
                fprintf(stderr,"Failed storing peptide in list in "\
                        "generate_digested_subsequences.\n");
                retval = ERROR_FAILED;
            }
        }
    }

    return retval;
}

/* -------------- */

int free_uint_vector_t (

        uint_vector_t *vec

){
    int retval=0;

    if (vec == NULL){
        fprintf(stderr,"NULL header pointer in free_uint_vector_t.\n");
        retval = ERROR_FAILED;
    }
    else {
        if ((vec->element) == NULL) {
            fprintf(stderr,"Failed deleting data in free_uint_vector_t.\n");
            retval = ERROR_FAILED;
        }
        else {
            free(vec->element);
        }
        free(vec);
    }

    return retval;
}

/* ------------- */

#define PEPTIDE_LIST_LENGTH_INIT 1000

/* This is the new preferred digestion function. */

int digest_sequence_trypsin (

        char          *sequence,
        unsigned int   minlen,
        unsigned int   maxlen,
        unsigned int   mclvgs,
        string_list_t *peplist

){
    uint_vector_t * digest_locs;
    int retval=0;


    if (list_trypsin_digestion_candidates(
                sequence,&digest_locs)){
        fprintf(stderr,"Failed listing digestion candidates "\
                "in digest_sequence_trypsin.\n");
        retval = ERROR_FAILED;
    }

    if (generate_digested_subsequences( sequence,
            digest_locs,minlen,maxlen,mclvgs,peplist)){
        fprintf(stderr,"Failed generating digests in "\
                "digest_sequence_trypsin.\n");
        retval = ERROR_FAILED;
    }

    if (free_uint_vector_t(digest_locs)){
        fprintf(stderr,"Failed cleaning up digestion location "\
                "vector in digest_sequence_trypsin.\n");
        retval = ERROR_FAILED;
    }

    return retval;
}

/* -------------- */

int allocate_string_list_t (

        unsigned int    len,
        string_list_t **sl_p

){
    string_list_t *sl;
    unsigned int   i;
    int            retval=0;

    if ((sl=calloc(1,sizeof(string_list_t)))==NULL){
        fprintf(stderr,"Failed to allocate header in allocate_string_list_t.\n");
        retval = ERROR_FAILED;
    }
    else {
        if ((sl->string=calloc(len,sizeof(char *)))==NULL){
            fprintf(stderr,"Failed to allocate data in allocate_string_list_t.\n");
            retval=ERROR_FAILED;
        }
        else {
            sl->len    = 0;
            sl->maxlen = len;
            *sl_p      = sl;
            for(i=0;i<(sl->maxlen);i++){
                sl->string[i] = NULL;
            }
        }
    }
    return(retval);
}

/* -------------- */

int free_string_list_t (

        string_list_t *sl

){
    unsigned long i;
    int           retval = 0;

    if (sl==NULL){
        fprintf(stderr,"NULL header pointer in free_string_list_t.\n");
        retval = ERROR_FAILED;
    }
    else {
        if (sl->string==NULL){
            fprintf(stderr,"NULL data pointer in free_string_list_t.\n");
            retval = ERROR_FAILED;
        }
        else {
            for(i=0;i<(sl->len);i++){
                if (sl->string[i]) free(sl->string[i]);
            }
        }
        free (sl);
    }

    return retval;
}

/* ----------  */

int clear_string_list_t (

        string_list_t *list

){
    unsigned int i;
    int          retval = 0;

    if (list==NULL){
        fprintf(stderr,"NULL header pointer in clear_string_list_t.\n");
        retval = ERROR_FAILED;
    }
    else {
        if ((list->string)==NULL){
            fprintf(stderr,"NULL data pointer in clear_string_list_t.\n");
            retval= ERROR_FAILED;
        }
        else {
            for(i=0;i<list->len;i++){
                if (list->string[i]) free(list->string[i]);
            }
            list->len = 0;
        }
    }

    return retval;
}

/* ------ */


/* This function determines the index
 * of a branch in the branch node array
 * given a chracter "symbol". */

int char_tree_getindex (
        
        char                   *alphabet,
        char                    symbol,
        unsigned int           *idx_p
        
){
    char *found;
    int   retval = 0;

    if ((found=strchr(alphabet,symbol))==NULL){
        fprintf(stderr,"Index of %c not found in "\
                "char_tree_getindex.\n",symbol);
        retval = ERROR_FAILED;
    }
    else {
        *idx_p = found - alphabet;
    }

    return retval;
}

/* ------ */

int char_tree_node_create (
        
        char_tree_definition_t  *def,
        char_tree_node_t       **node_p
        
){
    char_tree_node_t *node;
    int               retval=0;

    if ((node=calloc(1,sizeof(char_tree_node_t)))==NULL){
        fprintf(stderr,"FAiled allcating char node in char_tree_node_create.\n");
        retval = ERROR_FAILED;
    }
    else {
        node->def       = def;
        /*
        node->nbranches = 0;
        */
        node->branch    = NULL;
        node->payload   = NULL;
        *node_p = node;
    }

    return retval;
}

/* ------ */

int char_tree_init (

        char              *alph,
        char_tree_node_t **root_p

){
    char_tree_node_t       *root;
    char_tree_definition_t *def;
    int                     retval=0;

    if (alph==NULL) {
        fprintf(stderr,"NULL aplhabet pointer in char_tree_init.\n");
        retval = ERROR_FAILED;
    }
    if ((def=calloc(1,
                sizeof(char_tree_definition_t)))==NULL){
        fprintf(stderr,"Failed allocating alphabet "\
                "structure in char_tree_init.\n");
        retval = ERROR_FAILED;
    }

    if (retval)return retval;

    if (char_tree_node_create(def,&root)){
        fprintf(stderr,"Failed creating root node in char_tree_init.\n");
        retval = ERROR_FAILED;
    }
    else {
        root->def           = def;
        root->def->alphabet =  strdup(alph);
        root->def->len      =  strlen(root->def->alphabet);
        *root_p = root;
    }

    return retval;
}

/* ------ */

int char_tree_enter (

        char              *path,
        char_tree_node_t  *node,
        int              (*insertfunc) (char_tree_node_t *, void *),
        void              *payload

){
    unsigned int idx;
    int          retval = 0;

    if (strlen(path)){
        if (char_tree_getindex(node->def->alphabet,path[0],&idx)){
            fprintf(stderr,"Failed determining child index in char_tree_enter.\n");
            retval = ERROR_FAILED;
        }
        else {
            if (node->branch==NULL){
                if ((node->branch=calloc(node->def->len,
                                sizeof(char_tree_node_t *)))==NULL){
                    fprintf(stderr,
                            "Failed allocating node branches in char_tree_enter.\n");
                    retval = ERROR_FAILED;
                }
            }
            if (node->branch[idx]==NULL){
                if(char_tree_node_create(node->def,&(node->branch[idx]))){
                    fprintf(stderr,"Failed allocating child node in char_tree_enter.\n");
                    retval = ERROR_FAILED;
                }
            }
            path++;
            if (char_tree_enter(path,node->branch[idx],insertfunc,payload)){
                retval = ERROR_FAILED;
            }
        }
    }
    else {
        if (insertfunc(node,payload)){
            fprintf(stderr,"Failed updating node payload in char_tree_enter.\n");
            retval = ERROR_FAILED;
        }
    }

    return retval;
}

/* ------ */

/* Command-line tool is not cleaning up, but
 * it might be a good idea to clean up after use. */

int char_tree_destroy_node (

        char_tree_node_t  *node,
        int              (*destroy_func) (void *)

){
    unsigned int i;
    int          retval=0;

    if ((node->branch)==NULL){
        if (destroy_func(node->payload)){
            fprintf(stderr,"Failed destroying payload in char_tree_destroy_node.\n");
            retval = ERROR_FAILED;
        }
    }
    else {
        for(i=0;i<(node->def->len);i++){
            if (node->branch[i]!=NULL){
                if(char_tree_destroy_node(node->branch[i],destroy_func)){
                    retval = ERROR_FAILED;
                }
            }
        }
    }
    free(node);

    return retval;
}

/* ------ */

/* This is probably an entry point function
 * to destroy un used tree. NOTE that this
 * function should only be used once on three
 * root because it destrys the branching definition
 * structure that is is referenced allover the tree.  */

int char_tree_destroy (

        char_tree_node_t  *root,
        int              (*destroy_func) (void *)

){
    char_tree_definition_t *def;
    int                     retval=0;

    if (root==NULL){
        fprintf(stderr,"NULL root poiunter in char_tree_destroy.\n");
        retval = ERROR_FAILED;
    }
    else {
        if ((root->def)==NULL){
            fprintf(stderr,"NULL definition pointer in char_tree_destroy.\n");
            retval = ERROR_FAILED;
        }
        else {
            def = root->def; 
            if (char_tree_destroy_node(root,destroy_func)){
                fprintf(stderr,"Failed destroying root tree in char_tree_destroy.\n");
                retval = ERROR_FAILED;
            }
            if (def->alphabet) free(def->alphabet);
            free(def);
        }
    }

    return retval;
}

/* ------ */

int char_tree_query (
        
        char_tree_node_t  *node,
        char              *path,
        int              (*read_function)(void*,void*),
        void              *result
){
    unsigned int idx;
    int          retval=0;

    if (strlen(path)){
        if ((node->branch)==NULL){
            retval = ERROR_FAILED;
        }
        else {
            if (char_tree_getindex(node->def->alphabet,path[0],&idx)){
                fprintf(stderr,"Failed getting index for %c in "\
                        "char_tree_query.\n",path[0]);
                retval = ERROR_FAILED;
            }
            else {
                if ((node->branch[idx])==NULL){
                    retval = ERROR_FAILED;
                }
                else {
                    if (char_tree_query(node->branch[idx],path+1,
                                read_function,result)){
                        retval = ERROR_FAILED;
                    }
                }
            }
        }
    }
    else {
        if (read_function(node->payload,result)){
            fprintf(stderr,"Failed reading payload in char_tree_query.\n");
            retval = ERROR_FAILED;
        }
    }

    return retval;
}

/* ------ */



int free_peptide_list_t(
        
        peptide_list_t *peplist
        
){
    int i,retval=0;

    if(peplist == NULL){
        fprintf(stderr,"Failed relasing main object "\
                "in free_peptide_list_t.\n");
        retval = ERROR_FAILED;
    }
    else {
        for (i=0;i<peplist->len;i++){
            if (peplist->peptide[i] != NULL){
                free(peplist->peptide[i]);
            }
        }
        if(peplist->peptide != NULL){
            free(peplist->peptide);
        }
        if(peplist->offset != NULL){
            free(peplist->offset);
        }
        free(peplist);
    }

    return retval;
}

/* -------------- */

int allocate_peptide_list_t (

        unsigned int len,
        peptide_list_t **peplist_p

){
    peptide_list_t  *peplist;
    char           **seq;
    unsigned int    *ovec;
    int              retval = 0;

    if((peplist=calloc(1,sizeof(peptide_list_t)))==NULL){
        fprintf(stderr,"Failed allocating main struct "\
                "in allocate_peptide_list_t.\n");
        retval = ERROR_FAILED;
    }
    else {
        peplist->maxlen = DEFAULT_PEPTIDE_LIST_INIT_LEN;
        peplist->len    = 0;
        if((seq=calloc((peplist->maxlen),sizeof(char*)))==NULL){
            fprintf(stderr,"Failed allocating sequence vector "\
                    "allocate_peptide_list_t.\n");
            retval = ERROR_FAILED;
        }
        else {
            peplist->peptide = seq;
        }
        if((ovec=calloc((peplist->maxlen),sizeof(int)))==NULL){
            fprintf(stderr,"Failed allocating offset vector "\
                    "in allocate_peptide_list_t.\n");
            retval = ERROR_FAILED;
        }
        else {
            peplist->offset = ovec;
        }
        *peplist_p = peplist;
    }

    return retval;
}

/* ----- */

int clear_peptide_list_t (

        peptide_list_t *peplist

){
    int i,retval=0;

    for(i=0;i<(peplist->len);i++){
        if(peplist->peptide[i]!=NULL){
            free(peplist->peptide[i]);
            peplist->peptide[i] = NULL;
        }
    }
    peplist->len = 0;

    return retval;
}

/* ----- */


int append_to_peptide_list_t (
        
        peptide_list_t *peplist,
        char           *new_copy,
        unsigned int    offset
        
){
    int retval=0;

    if( (peplist->len) >= (peplist->maxlen) ){

        peplist->maxlen += DEFAULT_PEPTIDE_LIST_DELTA_LEN;

        if ((peplist->peptide=realloc(peplist->peptide,(peplist->maxlen)*sizeof(char*)))==NULL){
            fprintf(stderr,"Failed extending sequence vector in "\
                    "append_to_peptide_list_t.\n");
            retval = ERROR_FAILED;
        }
        if ((peplist->offset=realloc(peplist->offset,(peplist->maxlen)*sizeof(int)))==NULL){
            fprintf(stderr,"Failed extending offset vector in "\
                    "append_to_peptide_list_t.\n");
            retval = ERROR_FAILED;
        }
    }
    peplist->peptide[ peplist->len ] = new_copy;
    peplist->offset[ peplist->len ]  = offset;

    peplist->len++;

    return retval;
}

/* -------------- */

int generate_digested_subsequences_w_offsets (

        char           *sequence,
        uint_vector_t  *locs,
        unsigned int    minlen,
        unsigned int    maxlen,
        unsigned int    mclvgs,
        peptide_list_t *peplist

){
    unsigned int      end_loc,pep_len,last,i,j,
                      offset;
    char             *new1,tmp;
    int               retval=0;

    last      = locs->len - 1;

    for (i=0;i<last;i++){
        for(j=i+1;j<(locs->len);j++){

            offset = locs->element[i];

            pep_len = locs->element[j] - locs->element[i];
            if (pep_len < minlen) continue;
            if (pep_len > maxlen) break;
            if ( (j-i-1) > mclvgs ) break;

            end_loc           = locs->element[j];
            tmp               = sequence[end_loc];
            sequence[end_loc] = '\0';

            if ((new1 = strdup((char *)(sequence+locs->element[i])))==NULL){
                fprintf(stderr,"Failed copying peptide for list in "\
                        "generate_digested_subsequences.\n");
                retval = ERROR_FAILED;
            }

            sequence[end_loc] = tmp;

            if (append_to_peptide_list_t(peplist,new1,offset)){
                fprintf(stderr,"Failed storing peptide in list in "\
                        "generate_digested_subsequences.\n");
                retval = ERROR_FAILED;
            }
        }
    }

    return retval;
}

/* -------------- */


/* This is the new preferred digestion function. */

int digest_w_offset (

        char           *sequence,
        unsigned int    minlen,
        unsigned int    maxlen,
        unsigned int    mclvgs,
        peptide_list_t *peplist

){
    uint_vector_t * digest_locs;
    int retval=0;


    if (list_trypsin_digestion_candidates(
                sequence,&digest_locs)){
        fprintf(stderr,"Failed listing digestion candidates "\
                "in digest_w_offset.\n");
        retval = ERROR_FAILED;
    }

    if (generate_digested_subsequences_w_offsets(sequence,
            digest_locs,minlen,maxlen,mclvgs,peplist)){
        fprintf(stderr,"Failed generating digests in "\
                "digest_w_offset.\n");
        retval = ERROR_FAILED;
    }

    if (free_uint_vector_t(digest_locs)){
        fprintf(stderr,"Failed cleaning up digestion location "\
                "vector in digest_w_offset.\n");
        retval = ERROR_FAILED;
    }

    return retval;
}

/* ---- */

int create_new_fasta_list ( 
        
        params_t     *params,
        unsigned int  len
        
){
    char **fasta_list;
    int    retval=0;

    if ((fasta_list = calloc(len,sizeof(char*))) == NULL){
        fprintf(stderr,"Failed allocating fasta file list in create_new_fasta_list.\n");
        retval = ERROR_FAILED;
    }
    else {

        params->fasta      = fasta_list;
        params->max_fastas = len;
        params->n_fastas   = 0;

    }

    return retval;
}

/* ------ */

int extend_fasta_list (
        
        params_t     *params,
        unsigned int  len
        
){
    int retval=0;
    unsigned int new_len;

    new_len = params->max_fastas + len;

    if((params->fasta=realloc(params->fasta,new_len*sizeof(char*)))==NULL){
        fprintf(stderr,"Failed expanding fasta file list "\
                "in extend_fasta_list.\n");
        retval = ERROR_FAILED;
    }
    else {
        params->max_fastas = new_len;
    }

    return retval;
}

/* ------ */

int append_to_fastalist (
        
        char     *optarg,
        params_t *params
        
){
    int    retval=0;
    int    new_len;

    /* On first call we create a new entry list
     * to store the names of the fasta files. */
    if (params->fasta == NULL){
        if(create_new_fasta_list(params,INIT_MAX_FASTAS)){
            fprintf(stderr,"Failed creating new FASTA file list "\
                    "in append_to_fastalist.\n");
            retval = ERROR_FAILED;
        }
    }

    if ( params->n_fastas >= params->max_fastas ){
        if(extend_fasta_list(params,INIT_DELTA_FASTAS)){
            fprintf(stderr,"Failed extending FASTA file list in "\
                    "append_to_fastalist.\n");
            retval = ERROR_FAILED;
        }
    }

    params->fasta[ params->n_fastas ] = strdup(optarg);
    params->n_fastas++;

    return retval;
}

/* ------ */
 
int set_default_target_names (
        
        params_t *params
        
){
    char         **target_list;
    unsigned int   n_targets;
    string_list_t *targets;
    int            retval=0;

    n_targets = 15;

    if(allocate_string_list_t(TARGETS_LEN_INIT,&targets)){
        fprintf(stderr,"Failed allocating targets list in "\
                "set_default_target_names.\n");
        retval = ERROR_FAILED;
    }
    else{
        params->targets = targets;

        if(append_to_string_list_t_dup(targets,COLNAME_DEFAULT_PROTEIN_GROUP)){
            fprintf(stderr,"Failed appending 1st target field.\n");
            retval = ERROR_FAILED;
        }
        if(append_to_string_list_t_dup(targets,COLNAME_DEFAULT_PROTEIN_IDS)){
            fprintf(stderr,"Failed appending 2nd target field.\n");
            retval = ERROR_FAILED;
        }
        if(append_to_string_list_t_dup(targets,COLNAME_DEFAULT_GENES)){
            fprintf(stderr,"Failed appending 3rd target field.\n");
            retval = ERROR_FAILED;
        }
        if(append_to_string_list_t_dup(targets,COLNAME_DEFAULT_MODIFIED_SEQUENCE)){
            fprintf(stderr,"Failed appending 4th target field.\n");
            retval = ERROR_FAILED;
        }
        if(append_to_string_list_t_dup(targets,COLNAME_DEFAULT_STRIPPED_SEQUENCE)){
            fprintf(stderr,"Failed appending 5th target field.\n");
            retval = ERROR_FAILED;
        }
        if(append_to_string_list_t_dup(targets,COLNAME_DEFAULT_PRECURSOR_ID)){
            fprintf(stderr,"Failed appending 6th target field.\n");
            retval = ERROR_FAILED;
        }
        if(append_to_string_list_t_dup(targets,COLNAME_DEFAULT_Q_VALUE)){
            fprintf(stderr,"Failed appending 7th target field.\n");
            retval = ERROR_FAILED;
        }
        if(append_to_string_list_t_dup(targets,COLNAME_DEFAULT_GLOBAL_Q_VALUE)){
            fprintf(stderr,"Failed appending 8th target field.\n");
            retval = ERROR_FAILED;
        }
        if(append_to_string_list_t_dup(targets,COLNAME_DEFAULT_PROTEOTYPIC)){
            fprintf(stderr,"Failed appending 9th target field.\n");
            retval = ERROR_FAILED;
        }
        if(append_to_string_list_t_dup(targets,COLNAME_DEFAULT_MS1_AREA)){
            fprintf(stderr,"Failed appending 10th target field.\n");
            retval = ERROR_FAILED;
        }
        if(append_to_string_list_t_dup(targets,COLNAME_DEFAULT_FRAGMENT_QUANT_CORRECTED)){
            fprintf(stderr,"Failed appending 11th target field.\n");
            retval = ERROR_FAILED;
        }
        if(append_to_string_list_t_dup(targets,COLNAME_DEFAULT_PTM_Q_VALUE)){
            fprintf(stderr,"Failed appending 12th target field.\n");
            retval = ERROR_FAILED;
        }
        if(append_to_string_list_t_dup(targets,COLNAME_DEFAULT_PTM_SITE_CONFIDENCE)){
            fprintf(stderr,"Failed appending 13th target field.\n");
            retval = ERROR_FAILED;
        }
        if(append_to_string_list_t_dup(targets,COLNAME_DEFAULT_FILE_NAME)){
            fprintf(stderr,"Failed appending 14th target field.\n");
            retval = ERROR_FAILED;
        }
        if(append_to_string_list_t_dup(targets,COLNAME_DEFAULT_RUN)){
            fprintf(stderr,"Failed appending 15th target field.\n");
            retval = ERROR_FAILED;
        }
    }

    params->expand_name = strdup(COLNAME_DEFAULT_FRAGMENT_QUANT_CORRECTED);
    params->modpep_name = strdup(COLNAME_DEFAULT_MODIFIED_SEQUENCE);
    /* The minus one value denotes that the headers for these columns
     * have not yet been detected in the input file. */
    params->expand_indx = -1;
    params->modpep_indx = -1;

    params->proc_mode = NULL;

    return retval;
}

/* ------ */

/* TODO move this earlier in program
 * these are defaults init_params or so */

int set_input_terminators(

        params_t *params,
        char      tab,
        char      nline,
        char      nullchar

){
    int retval=0;

    params->null_char        = nullchar;
    params->column_separator = tab;
    params->line_terminator  = nline;

    return retval;
}

/* ------ */
#define DEFAULT_FRAGMENT_INTENSITY_EXPANSION_FLAG 1

int init_params ( 

        params_t **params_p

){
    params_t *params;
    int       retval=0;

    if((params = calloc(1,sizeof(params_t))) == NULL){
        fprintf(stderr,"Failed allocating basic parameter structure in "\
                "init_params.\n");
        retval = ERROR_FAILED;
    }
    else {

        /* This 0 is also important as it is the
         * initialization of the counter. */
        params->n_fastas   = 0;
        params->max_fastas = 0;
        /* THis NULL is important, otherwise
         * a new list won't be allocated at command-line
         * parsing. */
        params->fasta      = NULL;

        params->outfile    = strdup(DEFAULT_OUTFILE);
        params->infile     = NULL;

        if(set_default_target_names(params)){
            fprintf(stderr,"Failed adding default target columns in "\
                    "init_params.\n");
            retval = ERROR_FAILED;
        }

        if(set_input_terminators(params,DEFAULT_TAB,DEFAULT_NEWLINE,
                    DEFAULT_NULL)){
            fprintf(stderr,"Failed setting default table terminators"\
                    "in init_params.\n");
            retval = ERROR_FAILED;
        }

        params->buffer     = NULL;
        params->buffer_len = 0;

        params->verbose_level = VERBOSE_QUIET;

        params->minlen = DEFAULT_MIN_PEPTIDE_LEN;
        params->maxlen = DEFAULT_MAX_PEPTIDE_LEN;
        params->max_missed_cleavages = DEFAULT_MAX_MISSED_CLEAVAGES;

        /* Default value whether to expand fragment intensities
         * or not */
        params->fragment_intensity_expansion = DEFAULT_FRAGMENT_INTENSITY_EXPANSION_FLAG;

        *params_p = params;
        
    }

    return retval;
}

/* ------ */

int print_usage (
        
        char *basename
        
){
    int retval = 0;

    printf("\nThis is %s version %s written by Alex Henneman.\n",
            basename, VERSION_STRING);
    printf("\n");

    printf("Usage: %s [options] [inputfile]\n",basename);
    printf("\n");
    printf("This program reads an output file from DiaNN and uses the present\n");
    printf("Modified.Sequence column to detect modifications (currently only UniMod\n");
    printf("notation supported), lookup peptide offsets in supplied FASTA file (need\n");
    printf("to be specified using the -f flag.\n");
    printf("\nGuided by presence of the tryptic digest sequence in the protein databases,\n");
    printf("this program will generate info like offset to modification, in which proteins\n");
    printf("the peptide can be found and several other fields.\n");
    printf("\nFurthermore, the Fragment.Quant.Corrected column will be expanded into single\n");
    printf("fragment intensities on each row, and adding a precursor unique id field called.\n");
    printf("Fragment.Rel.Id, which has the same meaning for each precursor according to DiaNN\n");
    printf("documentation.\n");
    printf("\n");
    printf("Options:\n\n");
    
    printf(" -E <colname>     Set fragment intensity column name to <colname> [default = %s].\n", COLNAME_DEFAULT_FRAGMENT_QUANT_CORRECTED);
    printf(" -f <fasta_file>  Add <fasta_file> to peptide search space (NOTE Need\n");
    printf("                  at least one FASTA file to work).\n");
    printf(" -h               Print this kind of stuff.\n");
    printf(" -k <n>           Set digested peptide maximum number of missed cleavages\n");
    printf("                  equal to <n> [default = %d].\n", DEFAULT_MAX_MISSED_CLEAVAGES);
    printf(" -l <len>         Set digested peptide minimum length equal to <len> [default = %d].\n", DEFAULT_MIN_PEPTIDE_LEN);
    printf(" -L <len>         Set digested peptide maximum length equal to <len> [default = %d].\n", DEFAULT_MAX_PEPTIDE_LEN);
    printf(" -M <colname>     Set modified peptide column name to <colname> [default = %s].\n", COLNAME_DEFAULT_MODIFIED_SEQUENCE);
    printf(" -o <outfile>     Write output to file <outfile> [default = %s].\n", DEFAULT_OUTFILE);
    printf(" -v               Be verbose. Print unmatched backbone sequences.\n");
    printf(" -x               Turn fragment intensity expansion off.\n");
    printf("\n");
    printf("DESCRIPTION\n");
    printf("\n");
    printf("The program reads in a tsv format table and expects to find a set of specified column names. If\n");
    printf("found, these columns are copied to output, otherwise absent. These column names are File.Name,\n");
    printf("Run, Protein.Group, Protein.Ids, Genes, Modified.Sequence, Stripped.Sequence, Precursor.Id,\n");
    printf("Q.Value, Global.Q.Value, Proteotypic, Ms1.Area, Fragment.Quant.Corrected, PTM.Q.Value and\n");
    printf("PTM.Site.Confidence.\n");
    printf("\nOf the above column names, two are mandatory, and to prevent program termination their names\n");
    printf("can be changed from the command line. The first mandatory input column is Fragment.Quant.Corrected,\n");
    printf("which contains all fragment ion intensities and is simply expanded in the output into one intensity\n");
    printf("value per row. Output columns Fragment.Rel.Id and Fragment.Intensity contain an arbitrary generated\n");
    printf("fragment id number (to distinguish fragments belonging to the same modified peptide) and the expanded\n");
    printf("fragment intensity value.\n");
    printf("\nThe other mandatory column is the Modified.Sequence column, which is the main input and\n");
    printf("are row-wise processed into an unmodified peptide backbone sequence and a list of modifications\n");
    printf("localized on this peptide. At initiation, the program digests one or more FASTA protein databases\n");
    printf("into a quick-lookup data structure, and looks up the unmodified peptide, reconstructing all full\n");
    printf("protein offsets of the input peptide modifications. In the output the columns Peptide.Backbone,\n");
    printf("Peptide.Offsets, Nr.Mods and Nr.Phospho are appended specifying the peptide sequence that is looked\n");
    printf("up, the offsets of the modifications, number of modifications and how many of these are phosphorylations,\n");
    printf("respectively. If the lookup fails, the row will be missing from the output. If the lookup succeeds, the\n");
    printf("columns Fasta.Files and Proteins.Fasta contain the FASTA file(s) and protein names in which\n");
    printf("the peptide backbone is found. The columns Protein.Sites and Phospho.Site.Specs encode for each of the\n");
    printf("Proteins.Fasta entries, the full offsets of all modifications and only phosphorylations, respectively.\n");
    printf("In these columns, for each protein, parenthesis enclose a comma-separated list\n");
    printf("of amino-acids full protein offset specifications. Multiple proteins in Proteins.Fasta\n");
    printf("filed lead to a semi-colon separated list of such offset lists. In Protein.Sites additionally,\n");
    printf("the type of modification is specified by a small case letter prepended to amino acid offsets. The\n");
    printf("symbols n, p, m, d, u, g, t, e, q, a and c denote UniMod:1, UniMod:21, UniMod:35, UniMod:36, UniMod:121,\n");
    printf("UniMod:34, UniMod:37, UniMod:27, UniMod:28, UniMod:385, and uUnimod:4 respectively. The symbol x is used for all\n");
    printf("unrecognized modifications.\n");
    printf("\n");
    printf("Version %s %s\n",VERSION_STRING,STAMP_DATE);
    printf("\n");


    return retval;
}

/* ------ */

int string_list_t_lookup (
        
        char          *pin,
        string_list_t *hay,
        int           *loc_p,
        int           *found_p
        
){
    int loc,found,i,retval=0;

    if( hay == NULL){
        retval = ERROR_FAILED;
    }

    if(retval) return retval;

    loc   = 0;
    found = 0;

    for(i=0;i<(hay->len);i++){
        if(strcmp(pin,hay->string[i])==0){
            loc   = i;
            found = 1;
            break;
        }
    }
    *found_p = found;
    *loc_p   = loc;

    return retval;
}

/* ------ */

int string_list_t_replace ( 
        
        char          *old,
        char          *new1,
        string_list_t *targets
        
){
    char *buf;
    int   indx,found,retval=0;

    if(string_list_t_lookup(old,targets,&indx,&found)){
        fprintf(stderr,"Failed looking up expand column name "\
                " in string_list_t_replace.\n");
        retval = ERROR_FAILED;
    }
    else {
        if (found == 1){
            free(targets->string[indx]);
            if((buf=strdup(new1))==NULL){
                fprintf(stderr,"Failed allocating new entry in "\
                        "string_list_t_replace.\n");
                retval = ERROR_FAILED;
            }
            else {
                targets->string[indx] = buf;
            }
        }
        else {
            fprintf(stderr,"Name not found in string_list_t_replace.\n");
            retval = ERROR_FAILED;
        }
    }

    return retval;
}

/* ------ */

int set_modified_sequence_column_name (
        
        char     *newname,
        params_t *params
        
){
    char          *oldname,*buf;
    string_list_t *targets;
    int            retval=0,indx,found;

    if(string_list_t_replace(params->modpep_name,newname,params->targets)){
        fprintf(stderr,"Failed replacing target entry in "\
                "set_modified_sequence_column_name.\n");
        retval = ERROR_FAILED;
    }
    else {
        free(params->modpep_name);
        if((params->modpep_name=strdup(newname))==NULL){
            fprintf(stderr,"Failed storing new modified sequence column name "\
                   " command-line struct.\n");
            retval = ERROR_FAILED;
        }
    }

    return retval;
}
/* ------ */

int set_expand_column_name(
        
        char     *newname,
        params_t *params
        
){ 
    char          *oldname,*buf;
    string_list_t *targets;
    int            retval=0,indx,found;

    if(string_list_t_replace(params->expand_name,newname,params->targets)){
        fprintf(stderr,"Failed replacing target entry in "\
                "set_expand_column_name.\n");
        retval = ERROR_FAILED;
    }
    else {
        free(params->expand_name);
        if((params->expand_name=strdup(newname))==NULL){
            fprintf(stderr,"Failed storing new expand column name "\
                   " command-line struct.\n");
            retval = ERROR_FAILED;
        }
    }

    return retval;
}

/* ------ */

/* ========================================
 * This is one of the functions called from
 * main.
 * ======================================== */


int parse_command_line (
        
        int        argc,
        char     **argv,
        params_t **params_p
        
){
    char      c;
    params_t *params;
    int       retval=0;

    /* First we intialize program parameters. */
    if(init_params(params_p)){
        fprintf(stderr,"Failed intializing parameter struct "\
                "in parse_command_line.\n");
        retval = ERROR_FAILED;
    }
    else {
        params = *params_p;
    }

    if(retval) return(retval);

    while((c=getopt(argc,argv,"xl:L:k:E:M:vo:f:h"))!=-1){
        switch(c){
            case 'h':
                print_usage(basename(argv[0]));
                exit(0);
            case 'f':
                if(append_to_fastalist(optarg,params)){
                    fprintf(stderr,"Failed adding fasta file.\n");
                    retval = ERROR_FAILED;
                }
                break;
            case 'E':
                if(set_expand_column_name(optarg,params)){ 
                    fprintf(stderr,"Failed setting expand column name.\n");
                    retval = ERROR_FAILED;
                }
                break;
            case 'l':
                if(sscanf(optarg,"%d",&(params->minlen))!=1){
                    fprintf(stderr,"Failed setting peptide minimum length "\
                            "equal to %s.\n",optarg);
                    retval = ERROR_FAILED;
                }
                break;
            case 'L':
                if(sscanf(optarg,"%d",&(params->maxlen))!=1){
                    fprintf(stderr,"Failed setting peptide maximum length "\
                            "equal to %s.\n",optarg);
                    retval = ERROR_FAILED;
                }
                break;
            case 'k':
                if(sscanf(optarg,"%d",&(params->max_missed_cleavages))!=1){
                    fprintf(stderr,"Failed setting peptide maximum number of "\
                            "missed cleavages equal to %s.\n",optarg);
                    retval = ERROR_FAILED;
                }
                break;
            case 'M':
                if(set_modified_sequence_column_name(optarg,params)){
                    fprintf(stderr,"Failed setting modified sequence column name.\n");
                    retval = ERROR_FAILED;
                }
                break;
            case 'o':
                if((params->outfile = strdup(optarg))==NULL){
                    fprintf(stderr,"Failed setting output file name.\n");
                    retval = ERROR_FAILED;
                }
                break;
            case 'v':
                params->verbose_level = VERBOSE_LEVEL1;
                break;
            case 'x':
                params->fragment_intensity_expansion = 0;
                break;
            default:
                break;


        }
    }

    if (argc - optind < 1){
        params->infile = strdup(PIPE_IO_FILENAME);
    }
    else {
        if (argc - optind == 1){
            params->infile = argv[optind];
        }
        else {
            fprintf(stderr,"Multi input file mode is not implemented yet!.\n");
            exit(0);
        }

    }

    return retval;
}

/* ------ */

int append_protein_loc_to_leaf (

        char         *protein,
        int           offset,
        node_data_t *leaf_data

){
    int retval=0;

    if(leaf_data->nprots >= leaf_data->maxprots){

        leaf_data->maxprots += DEFAULT_NODE_PROTEIN_LIST_DELTA;

        if((leaf_data->protein=realloc(leaf_data->protein,(leaf_data->maxprots)*sizeof(char*)))==NULL){
            fprintf(stderr,"Failed extending protein list in "\
                    "append_protein_loc_to_leaf.\n");
            retval = ERROR_FAILED;
        }
        if((leaf_data->offset=realloc(leaf_data->offset,(leaf_data->maxprots)*sizeof(int)))==NULL){
            fprintf(stderr,"Failed extending offset list in "\
                    "append_protein_loc_to_leaf.\n");
            retval = ERROR_FAILED;
        }
    }

    if(retval) return retval;

    leaf_data->protein[leaf_data->nprots] = protein;
    leaf_data->offset[leaf_data->nprots]  = offset;
    leaf_data->nprots +=1;

    return retval;
}

/* ------ */

int append_fasta_id_to_leaf (

        int          fasta_id,
        node_data_t *leaf_data

){
    int i,is_new,retval=0;

    /* First determine whether id is allready
     * in here. */
    is_new = 1;
    for(i=0;i<(leaf_data->nfastas);i++){
        if( fasta_id == leaf_data->fasta_id[i]){
            is_new = 0;
            break;
        }
    }
    /* If it is not in here, it is new and we add
     * it. */
    if(is_new){
        if( leaf_data->maxfastas <= leaf_data->nfastas){
    
            leaf_data->maxfastas += DEFAULT_NODE_FASTA_LIST_DELTA;
    
            if((leaf_data->fasta_id=realloc(leaf_data->fasta_id,(leaf_data->maxfastas)*sizeof(int)))==NULL){
                fprintf(stderr,"Failed extending node fasta list in "\
                        "append_fasta_id_to_leaf.\n");
                retval = ERROR_FAILED;
            }
        }
    
        if(retval) return retval;
    
        leaf_data->fasta_id[leaf_data->nfastas] = fasta_id;
        leaf_data->nfastas += 1;
    }

    return retval;
}

/* ------ */

int append_to_leaf_data (

        message_t        *message,
        char_tree_node_t *leaf

){
    node_data_t *leaf_data;
    int          retval=0;

    leaf_data = (node_data_t *) leaf->payload;

    if(append_protein_loc_to_leaf(message->protein,message->offset,leaf_data)){
        fprintf(stderr,"Failed adding protein to leaf data "\
                "in append_to_leaf_data.\n");
        retval = ERROR_FAILED;
    }
    if(append_fasta_id_to_leaf(message->fasta_id,leaf_data)){
        fprintf(stderr,"Failed adding fasta id to leaf data in "\
                "append_to_leaf_data.\n");
        retval = ERROR_FAILED;
    }
    return retval;
}

/* ------ */

/* This fucntion creates a new node data object
 * and fills it in with fasta id, protein and
 * offset. */

int create_new_leaf_data (

        message_t        *message,
        char_tree_node_t *leaf

){
    node_data_t *leaf_data;
    int          retval=0;

    if ((leaf_data=calloc(1,sizeof(node_data_t)))==NULL){
        fprintf(stderr,"Failed allocating space for node "\
                "data in create_new_leaf_data.\n");
        retval = ERROR_FAILED;
    }
    else {
        leaf_data->maxfastas = DEFAULT_NODE_FASTA_LIST_INIT;
        if((leaf_data->fasta_id=calloc((leaf_data->maxfastas),sizeof(int)))==NULL){
            fprintf(stderr,"Failed allocating fasta list in "\
                    "create_new_leaf_data.\n");
            retval = ERROR_FAILED;
        }
        else {
            leaf_data->nfastas = 0;
            if(append_fasta_id_to_leaf(message->fasta_id,leaf_data)){
                fprintf(stderr,"Failed appending FASTA id in "\
                        "create_new_leaf_data.\n");
                retval = ERROR_FAILED;
            }
        }
        leaf_data->maxprots = DEFAULT_NODE_PROTEIN_LIST_INIT;
        if((leaf_data->protein=calloc((leaf_data->maxprots),sizeof(char*)))==NULL){
            fprintf(stderr,"Failed allocating protein list in "\
                    "create_new_leaf_data.\n");
            retval = ERROR_FAILED;
        }
        if((leaf_data->offset=calloc((leaf_data->maxprots),sizeof(int)))==NULL){
            fprintf(stderr,"Failed allocating offset list in "\
                    "create_new_leaf_data.\n");
            retval = ERROR_FAILED;
        }
        if(retval == 0) {
            leaf_data->nprots = 0;
            if(append_protein_loc_to_leaf(message->protein,message->offset,leaf_data)){
                fprintf(stderr,"Failed appening protein to leaf "\
                        "in create_new_leaf_data.\n");
                retval = ERROR_FAILED;

            }
        }
        leaf->payload = leaf_data;
    }

    return retval;
}

/* ------ */

/* This is the char_tree_node_t
 * data insertion function in which
 * payload ia  pointer to the message_t
 * struct. This data needs to be processed
 * and stored into the node_data_t struct
 * pertaining to leaf. */

int insert_digest_peptide (
        
        char_tree_node_t *leaf_node, 
        void             *payload
        
){
    message_t *message;
    int        retval=0;

    message = (message_t*) payload;

    if(leaf_node->payload == NULL){
        /* Leaf is visited for first time, so we
         * need to create everything. */
        if(create_new_leaf_data(message,leaf_node)){
            fprintf(stderr,"Failed creating new leaf_node data "\
                    "in insert_digest_peptide.\n");
            retval = ERROR_FAILED;
        }
    }
    else {
        /* Leaf is visited for a second time. */
        if(append_to_leaf_data(message,leaf_node)){
            fprintf(stderr,"Failed appending leaf_node data "\
                    "in insert_digest_peptide.\n");
            retval = ERROR_FAILED;
        }
    }

    return retval;
}

/* ------ */

int get_protein_name_copy (
        
        char *header,char **name_p
        
){
    char *start,*end,*copy;
    int   retval=0;

    if((start = strchr(header,'|'))==NULL){
        if((copy = strdup(header+1))==NULL){
            fprintf(stderr,"Failed making copy (1) of raw header "
                    "in get_protein_name_copy.\n");
            retval = ERROR_FAILED;
        }
        else{
            *name_p = copy;
        }
    }
    else {
        start++;
        if ((end = strchr(start,'|'))==NULL){
            if((copy = strdup(header+1))==NULL){
                fprintf(stderr,"Failed making copy (2) of raw header "
                        "in get_protein_name_copy.\n");
                retval = ERROR_FAILED;
            }
        else{
            *name_p = copy;
        }
        }
        else {
            *end = '\0';
            if((copy = strdup(start))==NULL){
                fprintf(stderr,"Failed cloning protein name in "\
                        "get_protein_name_copy.\n");
                retval = ERROR_FAILED;
            }
            else{
                *name_p = copy;
            }
        }
    }

    return retval;
}

/* ------ */

int store_fasta_into_peptide_tree (

        params_t         *params,
        char             *fname,
        int               fasta_id,
        char_tree_node_t *peptree
        
){
    char *protein_name;
    message_t      *message;
    peptide_list_t  *peplist;
    protein_list_t *prots;
    void           *payload;
    unsigned int i,j;
    int             maxlen,minlen,mclvgs;
    int             retval=0;

    peplist = NULL;

    if(read_fasta_file(fname,&prots)){
        fprintf(stderr,"Failed reading FASTA file %s.\n",fname);
        retval = ERROR_FAILED;
    }
    if(allocate_peptide_list_t(DIGEST_PEPTIDE_LIST_LEN_INIT,&peplist)){
        fprintf(stderr,"Failed allocating peptide list in "\
                "store_fasta_into_peptide_tree.\n");
        retval = ERROR_FAILED;
    }

    if(retval) return retval;


    if((message=calloc(1,sizeof(message_t)))==NULL){
        fprintf(stderr,"Failed allocating payload object in "\
                "store_fasta_into_peptide_tree.\n");
        retval = ERROR_FAILED;
    }
    else {
        payload = (void *) message;
    }

    minlen = params->minlen;
    maxlen = params->maxlen;
    mclvgs = params->max_missed_cleavages;

    for (i=0;i<(prots->len);i++){

        if(digest_w_offset(prots->protein[i].sequence,minlen,maxlen,mclvgs,peplist)){
            fprintf(stderr,"Failed digesting a protein in store_fasta_into_peptide_tree.\n");
            retval = ERROR_FAILED;
            break;
        }
        else {

            if(get_protein_name_copy(prots->protein[i].name,&protein_name)){
                fprintf(stderr,"Failed extracting protein name in store_fasta_into_peptide_tree.\n");
                retval = ERROR_FAILED;
            }


            message->protein = protein_name;
            message->fasta_id = fasta_id;

            for(j=0;j<(peplist->len);j++){

                message->offset   = peplist->offset[j];

                if(char_tree_enter(peplist->peptide[j],peptree,insert_digest_peptide,payload)){
                    fprintf(stderr,"Skipping peptide %s.\n",peplist->peptide[j]);
                    /* Let's not be so panicky about this...  retval = ERROR_FAILED; */
                    break;
                }

            }
            if(clear_peptide_list_t(peplist)){
                fprintf(stderr,"Failed emptying peptide list in store_fasta_into_peptide_tree.\n");
                retval = ERROR_FAILED;
            }

        }
        
    }

    fprintf(stderr,"Inserted %d proteins into peptide tree.\n",prots->len);

    if(free_protein_list_t(prots)){
        fprintf(stderr,"Failed freeing up protein list in store_fasta_into_peptide_tree.\n");
        retval = ERROR_FAILED;
    }
    if(free_peptide_list_t(peplist)){
        fprintf(stderr,"Failed freein up peptide list in store_fasta_into_peptide_tree.\n");
        retval = ERROR_FAILED;
    }

    return retval;
}

/* ------ */

int load_fastas_into_tree (
        
        params_t          *params,
        char_tree_node_t **peptree_p

){
    int i;
    int retval=0;
    char_tree_node_t *peptree;

    if(char_tree_init(SEQUENCE_ALPHABET,&peptree)){
        fprintf(stderr,"Failed intializing peptide tree in load_fastas_into_tree.\n");
        retval = ERROR_FAILED;
    }
    else {
        *peptree_p = peptree;
    }

    if(retval) return retval;

    fprintf(stderr,"Loading FASTA files:\n");
    for(i=0;i<params->n_fastas;i++){

        fprintf(stderr,"%i) Processing %s\n",i+1,params->fasta[i]);

        if(store_fasta_into_peptide_tree(params,params->fasta[i],i,peptree)){
            fprintf(stderr,"Failed storing FASTA file in peptide tree.\n");
            retval = ERROR_FAILED;
        }
    }
    fprintf(stderr,"\n");

    return retval;
}

/* ------ */

int allocate_ptm_locations_t (

        ptm_locations_t **poffsets_p

){
    ptm_locations_t *poff;
    int              retval=0;

    if((poff=calloc(1,sizeof(ptm_locations_t)))==NULL){
        fprintf(stderr,"Failed allocating main struct "\
                "in allocate_ptm_locations_t.\n");
        retval = ERROR_FAILED;
    }
    else {
        poff->node_data = NULL;
        poff->sitespec  = NULL;
        poff->fastas    = NULL;
        poff->proteins  = NULL;
        poff->modoffs   = NULL;
        poff->backbone  = NULL;
        poff->nmods = 0;

        *poffsets_p = poff;
    }

    return retval;
}

/* ------ */

/* This is the char_tree_node_t query
 * function used. It just copies a pointer
 * to the node data. */

int get_peptide_info( 

        void *payload,
        void *request

){
    node_data_t     *nodedata;
    ptm_locations_t *poff;
    int              retval=0;

    nodedata = (node_data_t *) payload;
    poff     = (ptm_locations_t *) request;

    poff->node_data  = nodedata; 

    return retval;
}

/* ------ */


int sprint_sitespec (
        
        char **sitespec_p,
        int   *reslen_p,
        char **start_p,
        int    peptide_offset,
        int   *mod_offs,
        char  *mtypes,
        int    nmods,
        char  *backbone
        
){
    char amino_acid;
    int  prot_offset,aa_pos,totlen,
         nprint,alarmlen,j,retval = 0;

    totlen = *start_p - *sitespec_p;
    alarmlen = (int) ceil((*reslen_p) * 0.95);

    for(j=0;j<nmods;j++){

        if(mod_offs[j] == 0){
            /* N-terminal, I think. */
            aa_pos = 1;
        }
        else {
            aa_pos = mod_offs[j];
        }
        amino_acid   = backbone[ aa_pos - 1 ];
        prot_offset  = aa_pos + peptide_offset;
        if(j>0){
            nprint = sprintf(*start_p,",%c%c%d",mtypes[j],amino_acid,prot_offset);
        }
        else {
            nprint = sprintf(*start_p,"%c%c%d",mtypes[j],amino_acid,prot_offset);
        }
        *start_p += nprint;
        totlen += nprint;

        if(totlen >= alarmlen){

            *reslen_p += SITESPEC_LEN_DELTA ;
            alarmlen = (int) ceil((*reslen_p) * 0.95);
            if( (*sitespec_p = realloc(*sitespec_p, (*reslen_p) *sizeof(char))) == NULL){
                fprintf(stderr,"Failed expanding sitespec buffer in "\
                        "sprint_sitespec.\n");
                retval = ERROR_FAILED;
            }
            else {
                *start_p = *sitespec_p + totlen;
            }
        }
    }

    return retval;
}

/* ------ */

/* NOTE This is a simplistic implementation
 * generating site amino acids. When location
 * is zeroth amino acid (PTM locations start
 * with amino acid number 1!), the modofication
 * is interpreted to apply to NEXT amino acid. */

int gen_sitespec (

        char         *backbone,
        node_data_t  *ndata,
        int           nmods,
        int          *mod_offs,
        char         *mtypes,
        char        **sitespec_p

){
    char         *start,*sitespec,amino_acid;
    int           nprint,res_len,prot_offset,
                  i,j,aa_pos,retval = 0;

    res_len = SITESPEC_LEN_INIT;

    if((sitespec=calloc(res_len,sizeof(char)))==NULL){
        fprintf(stderr,"Failed allocating response buffer in "\
                "gen_sitespec.\n");
        retval = ERROR_FAILED;
    }

    if(retval) return ERROR_FAILED;

    /* To keep track of written bytes. */
    start = sitespec;

    for(i=0;i<(ndata->nprots);i++){
        if(i>0){
            if (strcmp(ndata->protein[i],ndata->protein[i-1])==0 ){

                /* This is at least second one 
                 * for a same protein. */
                nprint = sprintf(start,"(");
                start += nprint;

                if(sprint_sitespec(&sitespec,&res_len,&start,
                            ndata->offset[i],mod_offs,mtypes,nmods,backbone)){
                    fprintf(stderr,"Failed writing site spec in "\
                            "gen_sitespec.\n");
                    retval = ERROR_FAILED;
                }

                nprint = sprintf(start,")");
                start += nprint;

            }
            else {

                /* This is at least second one 
                 * for a new protein. */
                nprint = sprintf(start,";(");
                start += nprint;

                if(sprint_sitespec(&sitespec,&res_len,&start,
                            ndata->offset[i],mod_offs,mtypes,nmods,backbone)){
                    fprintf(stderr,"Failed writing site spec in "\
                            "gen_sitespec.\n");
                    retval = ERROR_FAILED;
                }

                nprint = sprintf(start,")");
                start += nprint;

            }
        }
        else {
            /* This is the first one */
            nprint = sprintf(start,"(");
            start += nprint;

            if(sprint_sitespec(&sitespec,&res_len,&start,
                        ndata->offset[i],mod_offs,mtypes,nmods,backbone)){
                fprintf(stderr,"Failed writing site spec in "\
                        "gen_sitespec.\n");
                retval = ERROR_FAILED;
            }
    
            nprint = sprintf(start,")");
            start += nprint;
        }
    }
    *start = '\0';
    
    *sitespec_p = sitespec;

    return retval;
}

/* ------ */


int gen_fastas_found (

        params_t     *params,
        node_data_t  *ndata,
        char        **fastas_p

){
    char *start,*buffer,*fastaname;
    int   nwrite,fasta_indx,buffer_len,
          tot_len,i,retval=0;

    *fastas_p = NULL;
    tot_len = 0;

    /* First calculate the amount of
     * characters we neen for all the names */
    for(i=0;i<(ndata->nfastas);i++){

        fasta_indx = ndata->fasta_id[i];
        fastaname = basename(params->fasta[ fasta_indx ]);
        tot_len += strlen(fastaname) + 1;
    }
    buffer_len = tot_len + (ndata->nfastas -1) + 5;

    /* And we create a buffer to list the names. */
    if((buffer=calloc(buffer_len,sizeof(char)))==NULL){
        fprintf(stderr,"Failed allocating fasta file buffer "\
                "in gen_fastas_found.\n");
        retval = ERROR_FAILED;
    }
    else {

        start = buffer;

        for(i=0;i<(ndata->nfastas);i++){

            fasta_indx = ndata->fasta_id[i];
            fastaname = basename(params->fasta[ fasta_indx ]);
            if(i>0){
                nwrite = sprintf(start,";%s",fastaname);
            }
            else {
                nwrite = sprintf(start,"%s",fastaname);
            }
            start += nwrite;
        }
        *fastas_p = buffer;
    }

    return retval;
}

/* ------ */

/* These functions are now safe
 * up to an additonal ")" or so.
 *
 * 1) gen_sitespec
 * 2) gen_proteins_found
 *
 * */


int gen_proteins_found (

        node_data_t  *ndata,
        char        **proteins_p

){
    char *proteins,*start;
    int   tot_len,res_len,nprint,i,retval=0,
          alarm_len;

    res_len = PROTEINS_FOUND_BUFFER_LEN_INIT;

    if((proteins=calloc(res_len,sizeof(char)))==NULL){
        fprintf(stderr,"Failed allocating protein spec string "\
                "in gen_proteins_found.\n");
        retval = ERROR_FAILED;
    }

    if(retval) return(retval);

    start = proteins;
    tot_len = 0;
    alarm_len = (int) ceil( res_len * 0.9);

    for(i=0;i<(ndata->nprots);i++){
        if(i>0){
            nprint = sprintf(start,";%s",ndata->protein[i]);
            start += nprint;
        }
        else {
            nprint = sprintf(start,"%s",ndata->protein[i]);
            start += nprint;
        }
        tot_len += nprint;

        if ( tot_len >= alarm_len ){

            res_len += PROTEINS_FOUND_BUFFER_LEN_DELTA ;

            if( (proteins=realloc(proteins,res_len*sizeof(char))) == NULL){
                fprintf(stderr,"Failed expanding found protein "\
                        "buffer length in gen_proteins_found.\n");
                retval = ERROR_FAILED;
            }
            else {
                start = proteins + tot_len;
            }
        }

    }
    *start = '\0';

    *proteins_p = proteins;

    return retval;
}

/* ------ */

/*
 * TODO
 * ====
 *
 * Perhaps it is a good idea to show 
 * this in the blurp
 *
 * # UniMod:1   acetylation          n
 * # UniMod:21  phosphorylation      p
 * # UniMod:4   carbamidomethylation c                  c
 * # UniMod:35  oxydation            m
 * # UniMod:36  dimethylation        d
 * # UniMod:121 ubiquitinilation     u
 * # UniMod:34  methylation          g
 * # UniMod:37  tri-methylation      t
 * # UniMod:27  pyro_glu from        e
 * # UniMod:28  pyro_glu from        q
 * # UniMod:385 ammonia loss         a
 *
 * Unknown UniMod is processed as x
 * */

#define UNDEFINED_UNIMOD_RETVAL 'x'

char get_mod_type (
        
        char* unimod_descr
        
){
    char retval=UNDEFINED_UNIMOD_RETVAL;

    if( strcmp(unimod_descr,"UniMod:1") == 0){
        retval = 'n';
    }
    if( strcmp(unimod_descr,"UniMod:21") == 0){
        retval = 'p';
    }
    if( strcmp(unimod_descr,"UniMod:4") == 0){
        retval = 'c';
    }
    if( strcmp(unimod_descr,"UniMod:35") == 0){
        retval = 'm';
    }
    /* --- These are the new ones ---- */
    if( strcmp(unimod_descr,"UniMod:36") == 0){
        retval = 'd';
    }
    if( strcmp(unimod_descr,"UniMod:121") == 0){
        retval = 'u';
    }
    if( strcmp(unimod_descr,"UniMod:34") == 0){
        retval = 'g';
    }
    if( strcmp(unimod_descr,"UniMod:37") == 0){
        retval = 't';
    }
    if( strcmp(unimod_descr,"UniMod:27") == 0){
        retval = 'e';
    }
    if( strcmp(unimod_descr,"UniMod:28") == 0){
        retval = 'q';
    }
    if( strcmp(unimod_descr,"UniMod:385") == 0){
        retval = 'a';
    }

    if(retval == UNDEFINED_UNIMOD_RETVAL){
        fprintf(stderr,"Warning: unrecognized Unimod code:%s\n",
                unimod_descr);
    }

    return retval;
}

/* ------ */

int extract_backbone (

        char  *modpep,
        char **backbone_p,
        int   *nmods_p,
        int  **mod_offs_p,
        char **mtype_p

){
    char *mtypes,*src,*dest,*next_start,*next_end,
         *last_null,*start,*backbone;
    int chars_to_go,nmods,maxmods,*moffs,retval=0;

    maxmods = MOD_VEC_LEN_INIT ;

    /* Here we stoe the peptide offsets of
     * the fond modifications. */
    if((moffs=calloc(maxmods,sizeof(int)))==NULL){
        fprintf(stderr,"Failed allocating mod array in "\
                "extract_backbone.\n");
        retval = ERROR_FAILED;
    }
    if((mtypes=calloc(maxmods,sizeof(char)))==NULL){
        fprintf(stderr,"Failed allocating mod type array in "\
                "extract_backbone.\n");
        retval = ERROR_FAILED;
    }
    if((backbone=strdup(modpep))==NULL){
        fprintf(stderr,"Failed duplicating modpep in extract_backbone.\n");
        retval = ERROR_FAILED;
    }

    if(retval) return retval;

    nmods = 0;
    start = backbone;

    /* We first determine the end of the sequence for the shifting. */
    if((last_null = strchr(start,'\0'))==NULL){
        fprintf(stderr,"Failed finding end of mod peptide string "\
                "in extract_backbone.\n");
        retval = ERROR_FAILED;
    }

    while(1){
        if ((next_start = strchr(start,'(')) == NULL){
            /* No more mods */
            break;
        }
        else {
            /* We found one */
            if ((next_end = strchr(next_start,')')) == NULL){
                fprintf(stderr,"Unable to find Unimod closing character "\
                        "in extract_backbone.\n");
                retval = ERROR_FAILED;
                break;
            }
            else {
                *next_end = '\0';
                if(nmods >= maxmods){
                    maxmods += MOD_VEC_LEN_DELTA ;
                    if ( (moffs=realloc(moffs,maxmods*sizeof(int)))==NULL ){
                        fprintf(stderr,"Failed axpanding mod vector in "\
                                "extract_backbone.\n");
                        retval = ERROR_FAILED;
                    }
                    if ( (mtypes=realloc(mtypes,maxmods*sizeof(char)))==NULL ){
                        fprintf(stderr,"Failed axpanding mod type vector in "\
                                "extract_backbone.\n");
                        retval = ERROR_FAILED;
                    }
                }
                /* Now we know for sure there's room */
                moffs[nmods] = next_start - backbone;
                mtypes[nmods] = get_mod_type(next_start+1);
                nmods++;

                /* Now we shft next_end+1 - last
                 * to next_start */
                chars_to_go = last_null - next_end;
                src  = next_end + 1;
                dest = next_start;
                while(chars_to_go >0){
                    *dest = *src;
                    src++;
                    dest++;
                    chars_to_go--;
                }

                start = next_start;
            }
        }
    }
    
    /* Passing  all results upwards */
    *nmods_p    = nmods;
    *mod_offs_p = moffs;
    *backbone_p = backbone;
    *mtype_p    = mtypes;

    return retval;
}

/* ------ */

/* NOTE I don't know whether thhese
 * results are used at all downstream.... 
 * FIXME. */

int gen_local_modoffsets (

        int    nmods,
        int   *mod_offs,
        char **modoffs_p

){
    char *start,*mspec;
    int   alarmlen,totlen,res_len,i,nprint,retval=0;

    res_len = MOD_OFFSET_STRING_LEN_INIT ;

    if((mspec=calloc(res_len,sizeof(char)))==NULL){
        fprintf(stderr,"Failed allocating local mod offset "\
                "output table value in gen_local_modoffsets.\n");
        retval = ERROR_FAILED;
    }
    else {
        totlen = 0;
        alarmlen = (int) ceil(0.9*res_len);

        start = mspec;

        for(i=0;i<nmods;i++){
            if (i>0){
                nprint = sprintf(start,";%d",mod_offs[i]);
                start += nprint;
            }
            else {
                nprint = sprintf(start,"%d",mod_offs[i]);
                start += nprint;
            }
            totlen += nprint;

            if( totlen >= alarmlen ){

                res_len += MOD_OFFSET_STRING_LEN_DELTA ;
                if((mspec=realloc(mspec,res_len*sizeof(char)))==NULL){
                    fprintf(stderr,"Failed expanding mod offset buffer "\
                            "in gen_local_modoffsets.\n");
                    retval = ERROR_FAILED;
                }
                else {
                    start = mspec + totlen;
                }
            }
        }
        *start = '\0';

        *modoffs_p = mspec;
    }


    return retval;
}

/* ------ */

int sprintf_phospho_spec (
        
        char **phspec_p,
        int   *reslen_p,
        char **start_p,
        int    pep_offset,
        int   *ph_offs,
        int    nphosphos,
        char  *backbone
        
){
    char amino_acid;
    int prot_offset,i,nprint,totlen,alarmlen,retval=0;

    alarmlen = (int) ceil(0.9*(*reslen_p));
    totlen = *start_p - *phspec_p;

    for(i=0;i<nphosphos;i++){

        amino_acid = backbone[ ph_offs[i] - 1 ];
        prot_offset = ph_offs[i] + pep_offset;

        if(i>0){
            nprint = sprintf(*start_p,",%c%d",amino_acid,prot_offset);
        }
        else {
            nprint = sprintf(*start_p,"%c%d",amino_acid,prot_offset);
        }
        *start_p += nprint;
        totlen   += nprint;

        if(totlen >= alarmlen){
            *reslen_p += SITESPEC_LEN_DELTA ;
            if ((*phspec_p = realloc( *phspec_p ,(*reslen_p)*sizeof(char)))==NULL){
                fprintf(stderr,"Failed expanding phospho spec buffer in "\
                        "sprintf_phospho_spec.\n");
                retval = ERROR_FAILED;
            }
            else {
                *start_p = *phspec_p + totlen;
            }
        }
    }

    return retval;
}

/* ------ */

int gen_phospho_spec (
        
        int              nmods,
        int             *mod_offs,
        char            *mtypes,
        ptm_locations_t *query_result
        
){
    node_data_t *ndata;
    char        *backbone,*phspec,*start;
    int         *ph_offs,pcount,reslen,i,nphospho,retval=0;

    ndata    = query_result->node_data;
    backbone = query_result->backbone;

    nphospho=0;
    for(i=0;i<nmods;i++){
        if (mtypes[i] == 'p'){
            nphospho++;
        }
    }
    query_result->nphospho    = nphospho;
    query_result->phosphospec = NULL;

    if(nphospho > 0){
        reslen = SITESPEC_LEN_INIT;
        if((phspec=calloc(reslen,sizeof(char)))==NULL){
            fprintf(stderr,"Failed allocating phosho spec in "\
                    "gen_phospho_spec.\n");
            retval = ERROR_FAILED;
        }
        else{
            start = phspec;

            if((ph_offs=calloc(nphospho,sizeof(int)))==NULL){
                fprintf(stderr,"Failed allocating phospho offsets "\
                        "in gen_phospho_spec.\n");
                retval = ERROR_FAILED;
            }
            else {
                pcount = 0;
                for(i=0;i<=nmods;i++){
                    if (mtypes[i] == 'p'){
                        ph_offs[pcount] = mod_offs[i];
                        pcount++;
                    }
                }
            }
            for(i=0;i<(ndata->nprots);i++){
                if(i>0){
                    if (strcmp(ndata->protein[i],ndata->protein[i-1])==0 ){
                        start += sprintf(start,"(");
                        if(sprintf_phospho_spec(&phspec,&reslen,&start,
                                    ndata->offset[i],ph_offs,pcount,backbone)){

                            fprintf(stderr,"Failed adding to phosho spec from "\
                                    "same protein in gen_phospho_spec.\n");
                            retval = ERROR_FAILED;
                        }
                        start += sprintf(start,")");
                    }
                    else {
                        start += sprintf(start,";(");
                        if(sprintf_phospho_spec(&phspec,&reslen,&start,
                                    ndata->offset[i],ph_offs,pcount,backbone)){

                            fprintf(stderr,"Failed adding to phosho spec from "\
                                    "different protein in gen_phospho_spec.\n");
                            retval = ERROR_FAILED;
                        }
                        start += sprintf(start,")");
                    }
                }
                else {
                    start += sprintf(start,"(");
                    if(sprintf_phospho_spec(&phspec,&reslen,&start,ndata->offset[i],ph_offs,pcount,backbone)){
                        fprintf(stderr,"Failed adding to phosho spec from "\
                                "first protein in gen_phospho_spec.\n");
                        retval = ERROR_FAILED;
                    }
                    start += sprintf(start,")");
                }
            }
            query_result->phosphospec = phspec;
        }
    }

    return retval;
}

/* ------ */

/* NOTE That in this function we introduce a flag
 * in query_result->found_flag to indicate upwards
 * whether peptide sequence was foun in trie or not. */

int lookup_mod_peptide (

        params_t         *params,
        char             *modpep,
        char_tree_node_t *peptree,
        ptm_locations_t  *query_result

){
    char *mtypes,*modoffs,*backbone,*sitespec,
         *fastas,*proteins;
    int   nmods,*mod_offs;
    int   retval=0;

    query_result->backbone = NULL;

    /* NOTE That this function allocates a modification
     * offset vector of some specified length. This int 
     * vector need to be freed. */
    if(extract_backbone(modpep,&backbone,&nmods,&mod_offs,&mtypes)){
        fprintf(stderr,"Failed getting peptide backbone "\
                "in lookup_mod_peptide.\n");
        retval = ERROR_FAILED;
    }
    else { 
        query_result->nmods    = nmods;
        query_result->backbone = backbone;

        if(char_tree_query(peptree,backbone,get_peptide_info,(void*)query_result)){
            query_result->found_flag = 0;
        }
        else {
            query_result->found_flag = 1;

            if ( nmods <= 0 ){
                /* NO mods! */
                query_result->modoffs  = NULL;
                query_result->sitespec = NULL;

                query_result->nmods   = 0;
                query_result->nphospho = 0;
                query_result->phosphospec = NULL;
            } else {
                if(gen_phospho_spec(nmods,mod_offs,mtypes,query_result)){
                    fprintf(stderr,"Failed generatng phospho spec in "\
                            "lookup_mod_peptide.\n");
                    retval = ERROR_FAILED;
                }
                /* This peptide contains modifications
                 * so we can generate these structures */

                if(gen_local_modoffsets(nmods,mod_offs,&modoffs)){
                    fprintf(stderr,"Failed generating local modification offsets string "
                            "in lookup_mod_peptide.\n");
                    retval = ERROR_FAILED;
                }
                else {
                    query_result->modoffs = modoffs;
                }
                if(gen_sitespec(backbone,query_result->node_data,nmods,mod_offs,mtypes,&sitespec)){
                    fprintf(stderr,"Failed creating site spec in lookup_mod_peptide.\n");
                    retval = ERROR_FAILED;
                }
                else {
                    query_result->sitespec = sitespec;
                }
            }

            if(gen_fastas_found(params,query_result->node_data,&fastas)){
                fprintf(stderr,"Failed making fasta field in lookup_mod_peptide.\n");
                retval = ERROR_FAILED;
            }
            else {
                query_result->fastas = fastas;
            }
            if(gen_proteins_found(query_result->node_data,&proteins)){
                fprintf(stderr,"Failed generating proteins spec in lookup_mod_peptide.\n");
                retval = ERROR_FAILED;
            }
            else {
                 query_result->proteins = proteins;
            }
        }

        free(mod_offs);
        free(mtypes);
    }

    return retval;
}

/* ------ */

/* NOTE That this function counts on the fragment
 * vectors to be separated with ";" characters and
 * also provided of a terminating ";" character. */

int expand_intensities_field (

        char          *frag_intens,
        string_list_t *exp_intens

){
    char *start,field_separator,*field,
         *next,tmp;
    int   retval=0;

    field_separator = ';';

    if(clear_string_list_t(exp_intens)){
        fprintf(stderr,"Failed clearing fragment list "\
                "in expand_intensities_field.\n");
        retval = ERROR_FAILED;
    }

    /* A string_list_t carries copies
     * of the strings not just pointers. */
    start = frag_intens;
    while((next=strchr(start,field_separator))!=NULL){

        tmp = *next;
        *next = '\0';

        if( ( field = strdup(start) )==NULL ){
            fprintf(stderr,"Failed duplicating field in "\
                    "expand_intensities_field.\n");
            retval = ERROR_FAILED;
        }
        else {
            if(append_to_string_list_t(exp_intens,field)){
                fprintf(stderr,"FAiled appending fragment "\
                        "intensity in expand_intensities_field.\n");
                retval = ERROR_FAILED;
            }
        }
        *next = tmp;
        start = next+1;
    }

    return retval;
}

/* ------ */

/* NOTE That clearing this object leads to setting
 * all allocated data pointer back to NULL. This is
 * the signal that there's no data (mem) hanging around. */

int clear_ptm_locations_t (

        ptm_locations_t *prot_offsets

){
    int retval=0;

    if(prot_offsets->sitespec != NULL){
        free(prot_offsets->sitespec);
        prot_offsets->sitespec = NULL;
    }
    if(prot_offsets->fastas != NULL){
        free(prot_offsets->fastas);
        prot_offsets->fastas = NULL;
    }
    if(prot_offsets->proteins != NULL){
        free(prot_offsets->proteins);
        prot_offsets->proteins = NULL;
    }

    prot_offsets->node_data = NULL;

    return retval;
}

/* ------ */


int determine_proc_mode (

        params_t *params

){
    int *proc_mods,i,n_targets;
    int  retval=0;

    n_targets = params->targets->len;

    if((proc_mods=calloc(n_targets,sizeof(int)))==NULL){
        fprintf(stderr,"Failed allocating processing mode vector in "\
                "determine_proc_mode.\n");
        retval = ERROR_FAILED;
    }
    else {

        params->proc_mode = proc_mods;

        for(i=0;i<n_targets;i++){
            /* This is the default */
            params->proc_mode[i] = PROC_MODE_STATIC_REPEAT;
            if( i == params->expand_indx ){
                params->proc_mode[i] = PROC_MODE_FRAG_EXPAND ;
            }
            if( i == params->modpep_indx ){
                params->proc_mode[i] = PROC_MODE_MODPEP_QUERY ;
            }
        }
    }

    return retval;
}

/* ------ */

int find_next_target (

        char      **start_p,
        int         target_indx,
        params_t   *params,
        char      **next_start_p

){
    char *start,*next,flex_separator;
    int   retval=0,skips;

    start = *start_p;
    if(target_indx == 0){
        skips = params->column_index[0] - 1;
    }
    else {
        skips = params->column_index[target_indx] - params->column_index[target_indx-1] - 1;
    }

    /* Here we fast forward to next target */
    while(skips>0){
        if((next=strchr(start,params->column_separator))==NULL){
            fprintf(stderr,"Failed finding next tab in "\
                    "find_next_target.\n");
            retval = ERROR_FAILED;
            break;
        }
        else {
            start = next+1;
        }
        skips--;
    }
    /* Here we extract the field value */
    flex_separator = params->column_separator;
    if( target_indx == (params->targets->len - 1) ){
        if(params->column_index[params->targets->len - 1] == params->ncols){
            flex_separator = params->line_terminator;
        }
    }

    if((next=strchr(start,flex_separator))==NULL){
        fprintf(stderr,"Failed finding end tab in "\
                "find_next_target.\n");
        retval = ERROR_FAILED;
    }
    else {
        *next = params->null_char;
    }

    *next_start_p = next+1;
    *start_p = start;

    return retval;
}

/* ------ */

int append_to_static_fields_t (

        char *field,
        static_fields_t *stat

){
    int retval=0;

    if( stat->len >= stat->maxlen ){
        stat->maxlen += HEADER_LIST_LEN_DELTA;
        if((stat->string=realloc(stat->string,(stat->maxlen)*sizeof(char *)))==NULL){
            fprintf(stderr,"Failed expanding pointer "\
                    "vector in append_to_static_fields_t.\n");
            retval = ERROR_FAILED;
        }
    }

    stat->string[ stat->len ] = field;
    stat->len++;

    return retval;
}

/* ------ */

int allocate_static_fields_t (

        int               len,
        static_fields_t **stat_fields_p

){
    static_fields_t  *sf;
    char            **pvec;
    int               retval=0;

    if((sf=calloc(1,sizeof(static_fields_t)))==NULL){
        fprintf(stderr,"Failed allocating main struct in "\
                "allocate_static_fields_t.\n");
        retval = ERROR_FAILED;
    }
    else {
        if((pvec=calloc(len,sizeof(char*)))==NULL){
            fprintf(stderr,"Failed allocating pointer array in "\
                    "allocate_static_fields_t.\n");
            retval = ERROR_FAILED;
        }
        else {
            sf->string = pvec;
            sf->len    = 0;
            sf->maxlen = len;

            *stat_fields_p = sf;
        }
    }

    return retval;
}

/* ------ */

int free_static_fields_t (

        static_fields_t *sf

){
    int retval=0;

    if(sf==NULL){
        fprintf(stderr,"Failed releasing static_fields_t "\
                "base structure in free_static_fields_t.\n");
        retval = ERROR_FAILED;
    }
    else{
        if (sf->string == NULL){
            fprintf(stderr,"Failed relasing strng array in "\
                    "free_static_fields_t.\n");
            retval = ERROR_FAILED;
        }
    }
    return retval;
}

/* ------ */

int process_input_line (

        char             *line,
        char_tree_node_t *peptree,
        params_t         *params,
        ptm_locations_t  *ptm_loc_info,
        string_list_t    *exp_intens,
        static_fields_t  *stat_fields

){
    char *start,*next,flex_separator;
    int   i,retval=0;

    start = line;

    /* We loop over ordered targets in this line
     * and should have encountered them all in 
     * order by the end */
    for(i=0;i<(params->targets->len);i++){

        if(find_next_target(&start,i,params,&next)){
            fprintf(stderr,"Failed finding next target in "\
                    "process_input_line.\n");
            retval = ERROR_FAILED;
            break;
        }
        else {
            /* Depending on the TARGET column index (in target list) we
             * handle columns accordingly. */
            switch(params->proc_mode[i]){
                case PROC_MODE_FRAG_EXPAND :
                    if(params->fragment_intensity_expansion){
                        if(expand_intensities_field(start,exp_intens)){
                            fprintf(stderr,"Failed expanding fragment "\
                                    "intensities in process_input_line.\n");
                            retval = ERROR_FAILED;
                        }
                    }
                    else {
                        if(clear_string_list_t(exp_intens)){
                            fprintf(stderr,"Failed clearing intensity list "\
                                    "in process_input_line\n");
                            retval = ERROR_FAILED;
                        }
                        else {
                            if(append_to_string_list_t_dup(exp_intens,start)){
                                fprintf(stderr,"Failed adding full intensity "\
                                        "field in process_input_line.\n");
                                retval= ERROR_FAILED;
                            }
                        }
                    }
                    break;
                case PROC_MODE_MODPEP_QUERY :
                    if(lookup_mod_peptide(params,start,peptree,ptm_loc_info)){
                        fprintf(stderr,"Failed processing modified peptide "\
                                "in process_input_line.\n");
                        retval = ERROR_FAILED;
                    }
                    break;
                case PROC_MODE_STATIC_REPEAT:
                    break;
                default:
                    fprintf(stderr,"Unmatched processing mode in "\
                            "process_input_line\n");
                    retval = ERROR_FAILED;
                    break;
            }
            if(append_to_static_fields_t(start,stat_fields)){
                fprintf(stderr,"Failed adding current field to static "\
                        "fields list.\n");
                retval = ERROR_FAILED;
            }
            /* This is for the next round */
            start = next;
        }

    }

    return retval;
}

/* ------ */

/* NOTE THat we read,process and write each
 * line in the input, line by line. */

int write_expanded_output_line (

        FILE            *out_fd,
        int              expand_flag,
        ptm_locations_t *ptm_locs,
        string_list_t   *exp_intens,
        static_fields_t *stat_fields

){
    int i,j,retval=0;

    /* The outter loop generates a dedicated line for each
     * fragment. */
    if(expand_flag){
        for(j=0;j<(exp_intens->len);j++){
            /* Here we need to write out all desired output
             * lines. 
             * First the static stuff */
            for(i=0;i<(stat_fields->len);i++){
                if(i>0){
                    fprintf(out_fd,"\t%s",stat_fields->string[i]);
                }
                else {
                    fprintf(out_fd,"%s",stat_fields->string[i]);
                }
            }
            /* And here the dynamic stuff. */
    
            if( ptm_locs->nmods > 0){
                fprintf(out_fd,"\t%s",ptm_locs->sitespec);
            } 
            else {
                fprintf(out_fd,"\t");
            }
    
            fprintf(out_fd,"\t%s",ptm_locs->fastas);
            fprintf(out_fd,"\t%s",ptm_locs->proteins);
            fprintf(out_fd,"\t%s",ptm_locs->backbone);
    
            if( ptm_locs->nmods > 0){
                fprintf(out_fd,"\t%s",ptm_locs->modoffs);
            }
            else {
                fprintf(out_fd,"\t");
            }
    
            fprintf(out_fd,"\t%d",ptm_locs->nmods);
    
            fprintf(out_fd,"\t%d",ptm_locs->nphospho);
            if(ptm_locs->nphospho > 0){
                fprintf(out_fd,"\t%s",ptm_locs->phosphospec);
            }
            else {
                fprintf(out_fd,"\t");
            }
    
            /* Note for fragment intensities we need
             * to introduce an arbitraty fragment id. */
            fprintf(out_fd,"\t%d",j+1);
            fprintf(out_fd,"\t%s",exp_intens->string[j]);
            /* The end. */
            fprintf(out_fd,"\n");
        }
    }
    else {
        /* This is the case of no fragment expansion */
        /* Static crap */
        for(i=0;i<(stat_fields->len);i++){
            if(i>0){
                fprintf(out_fd,"\t%s",stat_fields->string[i]);
            }
            else {
                fprintf(out_fd,"%s",stat_fields->string[i]);
            }
        }
        /* And here the dynamic stuff. */

        if( ptm_locs->nmods > 0){
            fprintf(out_fd,"\t%s",ptm_locs->sitespec);
        } 
        else {
            fprintf(out_fd,"\t");
        }

        fprintf(out_fd,"\t%s",ptm_locs->fastas);
        fprintf(out_fd,"\t%s",ptm_locs->proteins);
        fprintf(out_fd,"\t%s",ptm_locs->backbone);

        if( ptm_locs->nmods > 0){
            fprintf(out_fd,"\t%s",ptm_locs->modoffs);
        }
        else {
            fprintf(out_fd,"\t");
        }

        fprintf(out_fd,"\t%d",ptm_locs->nmods);

        fprintf(out_fd,"\t%d",ptm_locs->nphospho);
        if(ptm_locs->nphospho > 0){
            fprintf(out_fd,"\t%s",ptm_locs->phosphospec);
        }
        else {
            fprintf(out_fd,"\t");
        }

        /* In the non-expansion case the first entry
         * is the original field value */
        fprintf(out_fd,"\t%s",exp_intens->string[0]);
        /* The end. */
        fprintf(out_fd,"\n");
    }

    return retval;
}

/* ------ */

int destroy_ptm_locations_t(
        
        ptm_locations_t *ptm_locs
        
){
    int retval=0;

    if(ptm_locs==NULL){
        fprintf(stderr,"Failed releasing main structure in "\
                "destroy_ptm_locations_t.\n");
        retval = ERROR_FAILED;
    }
    else {
        if(clear_ptm_locations_t(ptm_locs)){
            fprintf(stderr,"Failed clearing contents in "\
                    "destroy_ptm_locations_t.\n");
            retval = ERROR_FAILED;
        }
        else {
            free(ptm_locs);
        }
    }

    return retval;
}

/* ------ */

int process_input_bulk (
        
        FILE             *in_fd,
        params_t         *params,
        char_tree_node_t *peptree,
        FILE             *out_fd
        
){
    int              last_col_index,retval=0,skips;
    int              linenr=0,unfound_backbones = 0;
    string_list_t   *exp_intens;
    ptm_locations_t *ptm_locs;
    static_fields_t *stat_fields;

    /* This we need to do until file
     * is processed. */
    if(allocate_string_list_t(FRAGMENT_LIST_LEN_INIT,&exp_intens)){
        fprintf(stderr,"Failed allocating new fragment list "\
                "in process_input_bulk.\n");
        retval = ERROR_FAILED;

    }
    if(allocate_ptm_locations_t(&ptm_locs)){
        fprintf(stderr,"Failed allocating protein offstes "\
                "list in process_input_bulk.\n");
        retval = ERROR_FAILED;
    }
    if(determine_proc_mode(params)){
        fprintf(stderr,"Failed determining processing mode "\
                "in process_input_bulk.\n");
        retval = ERROR_FAILED;
    }

    if(allocate_static_fields_t(HEADER_LIST_LEN_INIT,&stat_fields)){
        fprintf(stderr,"Failed alocating static fields object in "\
                "process_input_bulk.\n");
        retval = ERROR_FAILED;
    }
    /* NOTE we hardcode that there are two variable fields */
    if(retval) return retval;

    while(!fgets_full_line(in_fd,&(params->buffer),&(params->buffer_len))){
        linenr++;

        if(process_input_line(params->buffer,peptree,params,ptm_locs,
                    exp_intens,stat_fields)){

            fprintf(stderr,"Failed parsing line in "\
                    "process_input_bulk.\n");
            retval = ERROR_FAILED;
        }
        else {
            if(ptm_locs->found_flag == 1){
                if(write_expanded_output_line(out_fd,
                            params->fragment_intensity_expansion,
                            ptm_locs,exp_intens,stat_fields)){
                    fprintf(stderr,"Failed writing output in "\
                            "process_input_bulk.\n");
                    retval = ERROR_FAILED;
                }
            }
            else {
                unfound_backbones++;
                if(params->verbose_level > VERBOSE_QUIET){
                    fprintf(stderr,"Skipping input line nr %d. "\
                            "Sequence look-up error nr %d :%s\n",linenr,
                        unfound_backbones,
                        ptm_locs->backbone);
                }
            }
        }

        /* Clean this up for the next line */
        stat_fields->len = 0;
        if(clear_ptm_locations_t(ptm_locs)){
            fprintf(stderr,"Failed clearing ptm_locations object in "\
                    "process_input_bulk.\n");
            retval = ERROR_FAILED;
        }

    }
    if(unfound_backbones){
        fprintf(stderr,"Skipped %d lines with peptide backbone sequences"\
                " not found in FASTA files.\n",unfound_backbones);
    }

    /* Cleaning up here 
     *
     * Todo
     * ====
     * Cleanup:
     * 1) stat_fields
     * 2) ptm_locs
     *
     * */
    if(free_string_list_t(exp_intens)){
        fprintf(stderr,"Failed cleaning up expanded intensities list "\
                "in process_input_bulk.\n");
        retval = ERROR_FAILED;
    }
    if(free_static_fields_t(stat_fields)){
        fprintf(stderr,"Failed releasing static fields value "\
                "list in in process_input_bulk.\n");
        retval = ERROR_FAILED;
    }
    if(destroy_ptm_locations_t(ptm_locs)){
        fprintf(stderr,"Failed releasing ptm query search structure "\
                "in process_input_bulk.\n");
        retval = ERROR_FAILED;
    }

/* ------ */


    return retval;
}

/* ------ */

int header_is_target (
        
        char          *header,
        string_list_t *targets,
        int           *found_indx_p
        
){
    int i,found=0;

    for(i=0;i<(targets->len);i++){
        if(strcmp(header,targets->string[i])==0){
            found = 1;
            *found_indx_p = i;
            break;
        }
    }

    return found;
}

/* ------ */

int cleanup_targets (

        params_t *params

){
    int i,deleted,retval=0;

    deleted = 0;

    for(i=0;i<(params->targets->len);i++){

        if (params->column_index[i] == -1){
            fprintf(stderr,"Column %s is absent from input. Deleting from targets list.\n",
                    params->targets->string[i]);
            deleted++;
        }
        else {
            params->targets->string[i-deleted] = params->targets->string[i];
            params->column_index[i-deleted]    = params->column_index[i];
        }
    }
    if(deleted>0){
        fprintf(stderr,"Deleted %d targets.\n",deleted);
    }
    params->targets->len -= deleted;

    return retval;
}

/* ------ */

int init_column_index ( 
        
        params_t *params 
        
){
    int *column_index,i,retval = 0;

    if((column_index=calloc((params->targets->len),sizeof(int)))==NULL){
        fprintf(stderr,"Failed allocating column index vector "\
                "in init_column_index.\n");
        retval = ERROR_FAILED;
    }
    else {
        params->column_index = column_index;

        for(i=0;i<(params->targets->len);i++){
            column_index[i] = UNFOUND_INDEX;
        }
    }

    return retval;
}

/* ------ */

int find_all_target_column_names(

        params_t      *params
        
){
    char         *start,*next,tmp;
    int          *column_index,ready,col_counter,
                 found_indx,retval=0;

    if(init_column_index(params)){
        fprintf(stderr,"Failed initializing column index in "\
                "find_all_target_column_names.\n");
        retval = ERROR_FAILED;
    }

    start       = params->buffer;
    col_counter = 0;
    ready       = 0;

    while (1){

        if((next=strchr(start,params->column_separator))==NULL){

            if((next=strchr(start,params->line_terminator))==NULL){

                fprintf(stderr,"Failed finding line terminator in "\
                        "find_all_target_column_names.\n");
            }
            else {
                ready = 1;
            }
        }
        col_counter++;

        tmp = *next;
        *next = params->null_char;

        /* Here we check whether this is a required field. */
        if(header_is_target(start,params->targets,&found_indx)){
            params->column_index[found_indx] = col_counter;
        }

        *next = tmp;
        start = next+1;

        if(ready){
            break;
        }
    }
    /* Note that we get rid of the targets we coud not find */
    if(cleanup_targets(params)){
        fprintf(stderr,"Failed cleaning up targets in find_all_target_column_name.\n");
        retval = ERROR_FAILED;
    }

    params->ncols = col_counter;

    return retval;
}

/* ------ */

/* This function swaps header clums in order
 * for the column indices to be ordered. */

int sort_target_headers(
        
        params_t *parms
        
){
    char *tmp;
    int   i,j,tmp2;
    int   retval=0;

    for(i=0;i<(parms->targets->len-1);i++){
        for(j=i+1;j<(parms->targets->len);j++){

            if ( parms->column_index[i] > parms->column_index[j] ){

                /* This is the corresponding column header being swapped */
                tmp                          = parms->targets->string[i];
                parms->targets->string[i]    = parms->targets->string[j];
                parms->targets->string[j]    = tmp;

                /* This is the column index being swapped */
                tmp2                   = parms->column_index[i];
                parms->column_index[i] = parms->column_index[j];
                parms->column_index[j] = tmp2;

            }
        }
    }

    return retval;
}

/* ------ */

int check_necessary_targets (
        
        params_t *parms

){
    int   i,retval=ERROR_FAILED;

    /* Thes are the necessary targets. If any of the two are
     * unchaged after this, an error is issued. */

    parms->expand_indx = -1;
    parms->modpep_indx = -1;

    for(i=0;i<(parms->targets->len);i++){

        if( strcmp(parms->expand_name,parms->targets->string[i])==0){
            parms->expand_indx = i;
        }
        if( strcmp(parms->modpep_name,parms->targets->string[i])==0){
            parms->modpep_indx = i;
        }
    }
    if( parms->expand_indx != -1 & parms->modpep_indx != -1){
        retval = 0;
    }
    else {
        fprintf(stderr,"Unable to find required targets in list in "\
                "sort_target_headers.\n");
        retval = ERROR_FAILED;
    }

    return retval;
}

/* ------ */

int parse_header (
        
        FILE      *in_fd,
        params_t  *params
        
){
    unsigned int  buffer_len;
    char         *buffer;
    unsigned int *column_index;
    int           ncols;
    int           retval=0;

    buffer = NULL;

    if(fgets_full_line(in_fd,&buffer,&buffer_len)){
        fprintf(stderr,"Failed reading first line in parse_header.\n");
        retval = ERROR_FAILED;
    }
    else {
        params->buffer     = buffer;
        params->buffer_len = buffer_len;

        if(find_all_target_column_names(params)){
            fprintf(stderr,"Failed locating all target columns in "\
                    "parse_header.\n");
            retval = ERROR_FAILED;
        }
        else {
            if(sort_target_headers(params)){
                fprintf(stderr,"Failed sorting target headers in "\
                        "parse_header.\n");
                retval = ERROR_FAILED;
            }
            else {
                if(check_necessary_targets(params)){
                    fprintf(stderr,"Necessary target headers missing.\n");
                    retval = ERROR_FAILED;
                }
            }
        }
    }

    return retval;
}

/* ------ */

/* NOTE That here we write the columns in the 
 * output file, specifying their order. This order
 * should be reflected elsewhere in the code. */

int write_header (
        
        params_t *params,
        FILE     *out_fd
        
){
    unsigned int i;
    int retval=0;

    for(i=0;i<(params->targets->len);i++){
        if(i == 0){
            fprintf(out_fd,"%s",params->targets->string[i]);
        }
        else {
            fprintf(out_fd,"\t%s",params->targets->string[i]);
        }
    }
    /* This is the header corresponding to
     * the new additional output columns. */
    if (params->fragment_intensity_expansion){
        fprintf(out_fd,EXTRA_OUTPUT_COLNAMES);
    }
    else {
        fprintf(out_fd,EXTRA_OUTPUT_COLNAMES_NOEXP);
    }

    fprintf(out_fd,"\n");

    return retval;
}

/* ------ */

/* In this funcition we check wether input and
 * output files are real files of stdin,stdout,
 * and open or pass the corresponding files upwards. */

int open_input_files (

        params_t  *params,
        FILE     **in_p,
        FILE     **out_p

){

    FILE *in_fd,*out_fd;
    int retval=0;

    if (strcmp(params->infile,PIPE_IO_FILENAME)==0){
        *in_p = stdin;
        fprintf(stderr,"Using stdin as input file.\n");
    }
    else {
        if((in_fd=fopen(params->infile,"r"))==NULL){
            fprintf(stderr,"Failed opeing input file in open_input_files.\n");
            retval = ERROR_FAILED;
        }
        else {
            *in_p =in_fd;
        }
    }

    if(strcmp(params->outfile,PIPE_IO_FILENAME)==0){
        *out_p = stdout;
        fprintf(stderr,"Writing to stdout as output file.\n");
    }
    else {
        if((out_fd=fopen(params->outfile,"w+"))==NULL){
            fprintf(stderr,"Failed opening output file in open_input_files.\n");
            retval = ERROR_FAILED;
        }
        else{
            *out_p = out_fd;
        }
    }

    return retval;
}

/* ------ */

int close_files(
        
        FILE *in_fd,
        FILE *out_fd
        
){
    int retval=0;

    if (fclose(in_fd)){
        fprintf(stderr,"Failed closinf input file in close_files.\n");
        retval  = ERROR_FAILED;
    }
    if (fclose(out_fd)){
        fprintf(stderr,"Failed closing output file in close_files.\n");
        retval  = ERROR_FAILED;
    }

    return retval;
}

/* ------ */

int process_input_file (
        
        params_t         *params,
        char_tree_node_t *peptree
        
){
    FILE     *in_fd,*out_fd;
    int       retval=0;

    if(open_input_files(params,&in_fd,&out_fd)){
        fprintf(stderr,"Failed opening files in process_input_file.\n");
        retval = ERROR_FAILED;
    }

    if(retval) return retval;

    fprintf(stderr,"Processing input file %s.\n",params->infile);

    /* In this step we determine the location of
     * the columns we are intersted in, and materialize all
     * in the params object. */
    if(parse_header(in_fd,params)){
        fprintf(stderr,"failed parsing header in process_input_file.\n");
        retval = ERROR_FAILED;
    }
    else {
        /* Here we just write out what we are going to do, that is
         * allready encoded in the params object. */
        if(write_header(params,out_fd)){
            fprintf(stderr,"Failed writing header in process_input_file.\n");
            retval = ERROR_FAILED;
        }
        /* Here we  process line-wise following our params encoding. */
        if(process_input_bulk(in_fd,params,peptree,out_fd)){
            fprintf(stderr,"Failed proceesing input file bulk.\n");
            retval = ERROR_FAILED;

        }
    }
    /* And clean up alle shit we left over. */
    if(close_files(in_fd,out_fd)){
        fprintf(stderr,"Failed closing files in process_input_file.\n");
        retval = ERROR_FAILED;
    }
    fprintf(stderr,"Wrote output to file %s.\n",params->outfile);

    return retval;
}

/* ------ */

/* === MAIN PROGRAM === */

int main(
        
        int    argc,
        char **argv
        
){
    char_tree_node_t *peptree;
    params_t         *params;
    int               retval=0;

    if(parse_command_line(argc,argv,&params)){
        fprintf(stderr,"FAiled parsing command-line.\n");
        retval = ERROR_FAILED;
    }
    if (load_fastas_into_tree(params,&peptree)){
        fprintf(stderr,"Failed loading FASTA files.\n");
        retval = ERROR_FAILED;
    }
    if (process_input_file(params,peptree)){
        fprintf(stderr,"Failed processing input.\n");
        retval = ERROR_FAILED;
    }

    return retval;
}
