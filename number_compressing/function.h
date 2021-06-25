#ifndef DEF_FUNCTION
#define DEF_FUNCTION

#define len_string 100 
#define max_char_read 100
#include <stdbool.h>
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<ctype.h>

char *trim_space(char *str);
char** str_split(char* str, char* car, char*reslt[]);
void print_string(char* elt);
void print_int(int*elt);
#endif
