#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "readini.h"

#define BUFLEN 200

/*
 * Co nam brakuje w programie pod Linucha
 */
int
getch() {
    return getchar();
}

int
stricmp(const char *s1, const char *s2) {
    return strcmp(s1, s2);
}

/**
 * funkcja usuwa z poczatku i konca ciagu znaki nowej linii i tabulatora
 * @param s ciag znakow
 * @return bufor wejsciowy
 */
char *strbtrim(char *s) {
    register char  *p;
    register char  *q;

    p = q = s;
    while((' ' == *p || *p == '\x0A' || *p == '\x09' || *p == '\r' || *p == '\t') && (*p != '\0'))
        p++;

    while(*p != '\0')
        *q++ = *p++;

    *q = '\0';

    q--;

    while((q >= s) && (*q==' ' || *q == '\x0A' || *q == '\x09' || *q == '\r'  || *p == '\t')) {
        q--;
    }

    q++;
    *q = '\0';

    return(s);
}


#define FREE(a,b,c)             free(a);
#define MALLOC(a,b,c)           malloc(a);
#define REALLOC(a,b,c,d)        realloc(a,b);
#define STRDUP(a,b,c)           strdup(a);

//-----------------------------------------------------------------
ReadINI::ReadINI (const char *nazwa) {
    name=nazwa;
    il=0;
}

//-----------------------------------------------------------------
ReadINI::~ReadINI() {
    for ( int i=0; i<il; i++ ) {
        delete nazwa[i];
        delete wart[i];
    }
}

//-----------------------------------------------------------------
char * ReadINI::getline(FILE *cfg) {
    char *ptr;
    static char buf[BUFLEN];
    char *trimmed_line;

    fgets(buf, BUFLEN-1, cfg);

    ptr = strchr(buf, ';');
    if (ptr)
        *ptr = 0;
    ptr = strchr(buf, '\n');
    if (ptr)
        *ptr = 0;
    trimmed_line = strbtrim(buf);
    if ( trimmed_line[0]=='[' && trimmed_line[strlen(buf)-1]!=']' ) {
        printf("\nBlad w formacie pliku \"%s\"\n", name);
        return 0;
    }
    // return STRDUP(buf, __FILE__, __LINE__);
    return buf;
}

//-----------------------------------------------------------------
char * ReadINI::GetPar(const char * group, const char * par) {
    FILE * cfg;


    cfg = fopen(name, "r");
    if (!cfg) {
        printf("\nBlad otwierania pliku \"%s\"\n", name);
        return 0;
    }


    char * group0 = new char[strlen(group)+3];
    if ( !group0 ) {
        printf("\nZbyt malo pami)ci - ReadINI 1.\n");
        getch();
        exit(0);
    }


    strcpy(group0,"[");
    strcat(group0,group);
    strcat(group0,"]");
    char *ptr = NULL;
    char *ptr1;
    char *ptr2;
    char *trimmed_line;


    while ( !feof(cfg) ) {


        ptr = getline(cfg);

        if(ptr == NULL)
            break;

        trimmed_line = strbtrim(ptr);


        if ( stricmp(group0,trimmed_line) == 0 ) {

            ptr = 0;
            delete [] group0;
            group0=0; // zerowanie wskaznika, gdy? delete wystjpuje dwukrotnie


            //IJ By3o:     while ( !feof(cfg) && strncmp(ptr,"[",1) != 0 )
            while ( !feof(cfg) && (ptr == 0 || strncmp(ptr,"[",1) != 0)) {


                ptr=getline(cfg);

                ptr1=STRDUP(ptr, __FILE__, __LINE__);

                if ( !ptr1 ) {
                    printf("\nZbyt malo pami)ci - ReadINI. 2\n");
                    getch();
                    exit(0);
                }

                ptr2 = strchr(ptr1, '=');
                if ( ptr2 )
                    *ptr2 = 0;

                strbtrim(ptr1);

                if ( stricmp(par,ptr1) == 0 ) {
                    FREE(ptr1, __FILE__, __LINE__);
                    break;
                }

                //IJ
                FREE(ptr1, __FILE__, __LINE__);

                ptr1=0;
                ptr=0;


            }

            break;
        }


        ptr=0;
    }


    fclose(cfg);
    delete [] group0;
    if ( !ptr )
        return 0;


    ptr1 = strchr(ptr, '=');
    if (!ptr1) {
        printf("\nBlad struktury pliku \"%s\", brak znaku \"=\".", name);
        return 0;
    }


    ptr1 = STRDUP(ptr1+1, __FILE__, __LINE__);

    strbtrim(ptr1);


    return ptr1;
}

//-----------------------------------------------------------------
int ReadINI::GetInt(const char * group, const char * par) {
    char * ptr;
    int i;
    ptr = GetPar(group,par);
    if ( !ptr )
        return 0;

    i = atoi(ptr);

    FREE(ptr, __FILE__, __LINE__);

    return i;
}

//-----------------------------------------------------------------
double ReadINI::GetFloat(const char * group, const char * par) {
    char * ptr;
    double i;
    ptr = GetPar(group,par);
    if ( !ptr )
        return 0;
    i = strtod(ptr,NULL);
    FREE(ptr, __FILE__, __LINE__);
    return i;
}

//-----------------------------------------------------------------
void ReadINI::Spacje (char *str) {
    int i=(strlen(str)-1);
    while ( str[i] == ' ' ) {
        str[i]=0;
        i--;
    }
    while ( str[0] == ' ' )
        for ( i=0; i<=(int)strlen(str); i++)
            str[i]=str[i+1];
}

//-----------------------------------------------------------------
char * ReadINI::GetPar_(const char * par) {
    for ( int i=0; i<il; i++ )
        if ( stricmp(par,nazwa[i]) == 0 )
            return wart[i];
    return 0;
}

//-----------------------------------------------------------------
char * ReadINI::GetPar_ (const char * group, const char * par) {
    for ( int i=0; i<il; i++ )
        if ( stricmp(par,nazwa[i]) == 0 )
            return wart[i];

    // char * p=GetPar( group, par);
    // if ( p )
    //       return p;
    // return 0;
    return GetPar(group, par);
}

//-----------------------------------------------------------------
int  ReadINI::Get2ram  (const char * group) {
    FILE * cfg;
    cfg = fopen(name, "r");
    if (!cfg) {
        printf("\nBlad otwierania pliku \"%s\"\n", name);
        return 0;
    }
    char *group0;
    group0 = (char *) new char[strlen(group)+3];
    if ( !group0 ) {
        printf("\nZbyt malo pami)ci - ReadINI. 3\n");
        getch();
        exit(0);
    }
    strcpy(group0,"[");
    strcat(group0,group);
    strcat(group0,"]");
    char *ptr;
    char *ptr1;
    il=0;
    while ( !feof(cfg) ) {
        ptr = getline(cfg);
        if ( stricmp(group0,ptr) == 0 ) {
            wart = (char **) MALLOC(sizeof(char*), __FILE__, __LINE__);
            //IJ     if (**wart==NULL )   {
            if (wart==NULL ) {
                printf("\nZbyt malo pami)ci - ReadINI. 4\n");
                getch();
                exit(0);
            }
            nazwa = (char **) MALLOC(sizeof(char*), __FILE__, __LINE__);
            //IJ     if (**nazwa==NULL)   {
            if (nazwa==NULL) {
                printf("\nZbyt malo pami)ci - ReadINI.5 \n");
                getch();
                exit(0);
            }
            ptr=0;
            ptr=getline(cfg);
            while ( !feof(cfg) && strncmp(ptr,"[",1)!=0 ) {
                ptr1 = strchr(ptr, '=');
                if ( ptr1 ) {
                    nazwa = (char **) REALLOC(nazwa,(il+1)*sizeof(char*), __FILE__, __LINE__);
                    //IJ         if (**nazwa==NULL)   {
                    if (nazwa==NULL) {
                        printf("\nZbyt malo pami)ci - ReadINI. 6\n");
                        getch();
                        exit(0);
                    }
                    wart = (char **) REALLOC(wart,(il+1)*sizeof(char*), __FILE__, __LINE__);
                    //IJ         if (**wart==NULL )   {
                    if (wart==NULL ) {
                        printf("\nZbyt malo pami)ci - ReadINI. 7\n");
                        getch();
                        exit(0);
                    }

                    //MG
                    wart[il]=STRDUP(ptr1+1, __FILE__, __LINE__);

                    if (wart[il] == NULL) {
                        printf("\nZbyt malo pami)ci - ReadINI. 8\n");
                        getch();
                        exit(0);
                    }
                    *ptr1 = 0;

                    //MG
                    nazwa[il]=STRDUP(ptr, __FILE__, __LINE__);

                    if (nazwa[il] == NULL ) {
                        printf("\nZbyt malo pami)ci - ReadINI. 9\n");
                        getch();
                        exit(0);
                    }
                    il++;
                }
                ptr1=0;
                ptr=0;
                ptr=getline(cfg);
            }
            break;
        }
        ptr=0;
    }
    fclose(cfg);
    delete [] group0;
    return 1;
}
//-----------------------------------------------------------------
//-----------------------------------------------------------------

extern char *PathIni(const char *ini_file, const char *group, const char *param) {

    char         *wsk;
    static  char zwrot[256];


    ReadINI *readini = new ReadINI(ini_file);

    if(!readini)
        return NULL;


    wsk = readini->GetPar(group, param);


    if (wsk != NULL)
        wsk = strcpy(zwrot, wsk);


    delete readini;

    return wsk;
}
