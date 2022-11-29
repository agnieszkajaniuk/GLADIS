#ifndef __READINI__
#define __READINI__


extern char *PathIni(const char *ini_file, const char *group, const char *param);

/*-----------------------------------------------------------------*/

class ReadINI {

private:

    int il;
    char ** nazwa;
    char ** wart;

    const char * name;
    char * getline(FILE *cfg);
    void Spacje (char *str);

public:

    ReadINI(const char * nazwa);
    /* ustala nazwj pliku INI  */
    ~ReadINI();
    char*  GetPar  (const char * group, const char * par);
    /* pobiera wartof parametru "par" z grupy "group", przeszukuje plik,
       zwraca: wskazanie na wartof "par" - sukces, 0 - error   */
    int  GetInt  (const char * group, const char * par);
    /* j.w., konwertuje zwrscon9 wartof na int
       zwraca: wartof "par" - sukces, 0 - error    */
    double GetFloat(const char * group, const char * par);
    /* j.w., konwertuje zwrscon9 wartof na double
       zwraca: wartof "par" - sukces, 0 - error    */
    int Get2ram (const char * group);
    /* wczytuje do pamijci nazwy i wartoci wszystkich parametrsw z grupy
       "group", grupa nie powinna mief wijcej ni? kilkadziesi9t parametrsw,
       zwraca: 1 - sukces, 0 - error */
    char*  GetPar_ (const char * par);
    /* pobiera wartof parametru "par" z grupy wczytanej do pamijci
       za pomoc9 funkcji "Get2ram(..)", je?eli "par" nie zosta3
       znaleziony (lub nic nie zosta3o wczytane) funkcja zwraca 0,
       zwraca: wskazanie na wartof "par" - sukces, 0 - error   */
    char*  GetPar_ (const char * group, const char * par);
    /* pobiera wartof parametru "par" z grupy wczytanej do pamijci
       za pomoc9 funkcji "Get2ram(..)", je?eli "par" nie zosta3
       znaleziony (lub nic nie zosta3o wczytane) plik INI jest
       przeszukiwany funkcj9 GetPar(..),
       zwraca: wskazanie na wartof "par" - sukces, 0 - error   */

};

#endif
