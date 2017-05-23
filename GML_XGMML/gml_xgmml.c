/*
 * Demo for the GML_parser. Usage: gml_xgmml.exe <gml_file>
 * Translates gml to xgmml (xml application)
 * John Punin (puninj@cs.rpi.edu) 
 * Original code by Marcus Raitner (raitner@fmi.uni-passau.de)
 * 
 * Date: March 05 2000
 */

#include "gml_parser.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define XML_HEADER "<?xml version=\"1.0\"?>\n<!DOCTYPE graph SYSTEM \"xgmml.dtd\">"

char * to_xmlstring(char *str)
{
    char *xstr = NULL, *p = NULL, *px = NULL;
    int len = strlen(str);
    int camp = 0, cgt = 0, clt = 0, cquot = 0;
    int ctot =0;
    p = str;
    while(*p != '\0') {
	switch(*p) {
	case '&':
	    camp++;
	    break;
	case '>':
	    cgt++;
	    break;
	case '<':
	    clt++;
	    break;
	case '\"': /* " */
	    cquot++;
	    break;
	default:
	    break;
	}
	p++;
    }

    ctot = camp+cgt+clt+cquot;
    if(ctot == 0)
	return (char *)NULL;
    xstr = malloc(len+ctot+1);
    if(!xstr) {
	fprintf(stderr,"No memory\n");
	exit(0);
    }
    p = str; px = xstr;

    while(*p != '\0') {
	switch(*p) {
	case '&':
	    strcpy(px,"&amp;");
	    px+=5;
	    break;
	case '>':
	    strcpy(px,"&gt;");
	    px+=4;
	    break;
	case '<':
	    strcpy(px,"&lt;");
	    px+=4;
	    break;
	case '\"': /* " */
	    strcpy(px,"&quot;");
	    px+=6;
	    break;
	default:
	    *px = *p;
	    px++;
	    break;
	}
	p++; 
    }
    *px = '\0';
    return xstr;
}

void GML_print_xgmml (struct GML_pair* list, int level) {
    
    struct GML_pair* tmp = list;
    int i;
    int flag = 0;

    if(level == 0)
	printf(XML_HEADER);

    while (tmp) {
	char *tstr = NULL;
	if(tmp->kind != GML_LIST)
	    printf (" %s", tmp->key);

	switch (tmp->kind) {
	case GML_INT:
	    printf ("=\"%ld\"",tmp->value.integer);
	    break;

	case GML_DOUBLE:
	    printf ("=\"%f\"", tmp->value.floating);
	    break;

	case GML_STRING:
	    tstr = to_xmlstring(tmp->value.string);
	    if(tstr) {
		printf ("=\"%s\"", tstr);
		free(tstr); tstr = NULL;
	    } else {
		printf ("=\"%s\"", tmp->value.string);
	    }
	    break;
	    
	case GML_LIST:
	    if(flag == 0 && level!=0)
		printf(">");
	    flag = 1;
	    printf ("\n<%s", tmp->key);
	    GML_print_xgmml (tmp->value.list, level+1);
	    printf ("\n</%s>", tmp->key);
	    break;

	default:
	    break;
	}
	
	tmp = tmp->next;
    }
    if(flag == 0)
	printf(">");
    if(level == 0)
	printf("\n");
    return;
}
	

int main (int argc, char* argv[]) {
  
    struct GML_pair* list;
    struct GML_stat* stat=(struct GML_stat*)malloc(sizeof(struct GML_stat));
    stat->key_list = NULL;
    
    if (argc != 2) fprintf (stderr,"Usage: gml_xgmml.exe <gml_file> \n");
    else {
	FILE* file = fopen (argv[1], "r");
	if (file == 0) fprintf (stderr,"\n No such file: %s", argv[1]);
	else {
	    GML_init ();
	    list = GML_parser (file, stat, 0);

	    if (stat->err.err_num != GML_OK) {
		fprintf (stderr,"An error occured while reading line %d column %d of %s:\n", stat->err.line, stat->err.column, argv[1]);
		
		switch (stat->err.err_num) {
		case GML_UNEXPECTED:
		    fprintf (stderr,"UNEXPECTED CHARACTER");
		    break;
		    
		case GML_SYNTAX:
		    fprintf (stderr,"SYNTAX ERROR"); 
		    break;
		    
		case GML_PREMATURE_EOF:
		    fprintf (stderr,"PREMATURE EOF IN STRING");
		    break;
		    
		case GML_TOO_MANY_DIGITS:
		    fprintf (stderr,"NUMBER WITH TOO MANY DIGITS");
		    break;
		    
		case GML_OPEN_BRACKET:
		    fprintf (stderr,"OPEN BRACKETS LEFT AT EOF");
		    break;
		    
		case GML_TOO_MANY_BRACKETS:
		    fprintf (stderr,"TOO MANY CLOSING BRACKETS");
		    break;
		
		default:
		    break;
		}
		
		fprintf (stderr,"\n");
	    }      
	    GML_print_xgmml(list,0);
	    GML_free_list (list, stat->key_list);
	}
    }
    return 0;
}

