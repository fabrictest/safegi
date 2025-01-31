
/* A Bison parser, made by GNU Bison 2.4.1.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C
   
      Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.
   
   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.4.1"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1

/* Using locations.  */
#define YYLSP_NEEDED 0



/* Copy the first part of user declarations.  */

/* Line 189 of yacc.c  */
#line 24 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"

#include "api.h"
#include "pbrt.h"
#include "paramset.h"
#include <stdarg.h>

#ifdef WIN32
#pragma warning(disable:4065)
#pragma warning(disable:4996)
#pragma warning(disable:4018)
#endif // WIN32

extern int yylex();
int line_num = 0;
string current_file;

#define YYMAXDEPTH 100000000

void yyerror(const char *str) {
    Severe("Parsing error: %s", str);
}



struct ParamArray {
    ParamArray() {
        isString = false;
        element_size = allocated = nelems = 0;
        array = NULL;
    }
    bool isString;
    int element_size;
    int allocated;
    int nelems;
    void *array;
};



struct ParamListItem {
    ParamListItem(const char *t, ParamArray *array) {
        arg = array->array;
        name = t;
        size = array->nelems;
        isString = array->isString;
        array->allocated = 0;
        array->nelems = 0;
        array->array = NULL;
    }
    const char *name;
    void *arg;
    int size;
    bool isString;
};



static std::vector<ParamListItem> cur_paramlist;

static ParamArray *cur_array = NULL;

static void AddArrayElement(void *elem) {
    if (cur_array->nelems >= cur_array->allocated) {
        cur_array->allocated = 2*cur_array->allocated + 1;
        cur_array->array = realloc(cur_array->array,
            cur_array->allocated*cur_array->element_size);
    }
    char *next = ((char *)cur_array->array) + cur_array->nelems * cur_array->element_size;
    Assert(cur_array->element_size == 4 || cur_array->element_size == 8);
    if (cur_array->element_size == 4)
        *((uint32_t *)next) = *((uint32_t *)elem);
    else
        *((uint64_t *)next) = *((uint64_t *)elem);
    cur_array->nelems++;
}



static void ArrayFree(ParamArray *ra) {
    if (ra->isString && ra->array)
        for (int i = 0; i < ra->nelems; ++i) free(((char **)ra->array)[i]);
    free(ra->array);
    delete ra;
}



static void FreeArgs() {
    for (uint32_t i = 0; i < cur_paramlist.size(); ++i)
        free((char *)cur_paramlist[i].arg);
    cur_paramlist.erase(cur_paramlist.begin(), cur_paramlist.end());
}



static bool VerifyArrayLength(ParamArray *arr, int required,
    const char *command) {
    if (arr->nelems != required) {
        Error("\"%s\" requires a %d element array! (%d found)",
                    command, required, arr->nelems);
        return false;
    }
    return true;
}


enum { PARAM_TYPE_INT, PARAM_TYPE_BOOL, PARAM_TYPE_FLOAT, PARAM_TYPE_POINT,
    PARAM_TYPE_VECTOR, PARAM_TYPE_NORMAL, PARAM_TYPE_RGB, PARAM_TYPE_XYZ,
    PARAM_TYPE_BLACKBODY, PARAM_TYPE_SPECTRUM,
    PARAM_TYPE_STRING, PARAM_TYPE_TEXTURE };
static const char *paramTypeToName(int type);
static void InitParamSet(ParamSet &ps, SpectrumType);
static bool lookupType(const char *name, int *type, string &sname);
#define YYPRINT(file, type, value)  { \
    if ((type) == ID || (type) == STRING) \
        fprintf ((file), " %s", (value).string); \
    else if ((type) == NUM) \
        fprintf ((file), " %f", (value).num); \
}




/* Line 189 of yacc.c  */
#line 197 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.cpp"

/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 1
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif


/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     STRING = 258,
     ID = 259,
     NUM = 260,
     LBRACK = 261,
     RBRACK = 262,
     ACCELERATOR = 263,
     ACTIVETRANSFORM = 264,
     ALL = 265,
     AREALIGHTSOURCE = 266,
     ATTRIBUTEBEGIN = 267,
     ATTRIBUTEEND = 268,
     CAMERA = 269,
     CONCATTRANSFORM = 270,
     COORDINATESYSTEM = 271,
     COORDSYSTRANSFORM = 272,
     ENDTIME = 273,
     FILM = 274,
     IDENTITY = 275,
     LIGHTSOURCE = 276,
     LOOKAT = 277,
     MAKENAMEDMATERIAL = 278,
     MATERIAL = 279,
     NAMEDMATERIAL = 280,
     OBJECTBEGIN = 281,
     OBJECTEND = 282,
     OBJECTINSTANCE = 283,
     PIXELFILTER = 284,
     RENDERER = 285,
     REVERSEORIENTATION = 286,
     ROTATE = 287,
     SAMPLER = 288,
     SCALE = 289,
     SHAPE = 290,
     STARTTIME = 291,
     SURFACEINTEGRATOR = 292,
     TEXTURE = 293,
     TRANSFORMBEGIN = 294,
     TRANSFORMEND = 295,
     TRANSFORMTIMES = 296,
     TRANSFORM = 297,
     TRANSLATE = 298,
     VOLUME = 299,
     VOLUMEINTEGRATOR = 300,
     WORLDBEGIN = 301,
     WORLDEND = 302,
     HIGH_PRECEDENCE = 303
   };
#endif



#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{

/* Line 214 of yacc.c  */
#line 147 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"

char string[1024];
float num;
ParamArray *ribarray;



/* Line 214 of yacc.c  */
#line 289 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.cpp"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif


/* Copy the second part of user declarations.  */


/* Line 264 of yacc.c  */
#line 301 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.cpp"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int yyi)
#else
static int
YYID (yyi)
    int yyi;
#endif
{
  return yyi;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef _STDLIB_H
#      define _STDLIB_H 1
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined _STDLIB_H \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef _STDLIB_H
#    define _STDLIB_H 1
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T yyi;				\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (YYID (0))
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)				\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack_alloc, Stack, yysize);			\
	Stack = &yyptr->Stack_alloc;					\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  73
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   117

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  49
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  20
/* YYNRULES -- Number of rules.  */
#define YYNRULES  65
/* YYNRULES -- Number of states.  */
#define YYNSTATES  134

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   303

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint8 yyprhs[] =
{
       0,     0,     3,     5,     6,     7,     8,    10,    12,    17,
      19,    22,    25,    27,    30,    35,    37,    40,    43,    45,
      48,    51,    52,    55,    56,    59,    62,    64,    68,    71,
      74,    77,    81,    83,    85,    89,    92,    95,    98,   102,
     104,   108,   119,   123,   127,   130,   133,   135,   138,   142,
     146,   148,   154,   158,   163,   167,   171,   177,   179,   181,
     185,   188,   193,   197,   201,   203
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      50,     0,    -1,    67,    -1,    -1,    -1,    -1,    55,    -1,
      59,    -1,    51,     6,    57,     7,    -1,    56,    -1,    51,
      58,    -1,    57,    58,    -1,    58,    -1,    52,     3,    -1,
      51,     6,    61,     7,    -1,    60,    -1,    51,    62,    -1,
      61,    62,    -1,    62,    -1,    53,     5,    -1,    64,    65,
      -1,    -1,    66,    65,    -1,    -1,     3,    54,    -1,    67,
      68,    -1,    68,    -1,     8,     3,    63,    -1,     9,    10,
      -1,     9,    18,    -1,     9,    36,    -1,    11,     3,    63,
      -1,    12,    -1,    13,    -1,    14,     3,    63,    -1,    15,
      59,    -1,    16,     3,    -1,    17,     3,    -1,    19,     3,
      63,    -1,    20,    -1,    21,     3,    63,    -1,    22,     5,
       5,     5,     5,     5,     5,     5,     5,     5,    -1,    23,
       3,    63,    -1,    24,     3,    63,    -1,    25,     3,    -1,
      26,     3,    -1,    27,    -1,    28,     3,    -1,    29,     3,
      63,    -1,    30,     3,    63,    -1,    31,    -1,    32,     5,
       5,     5,     5,    -1,    33,     3,    63,    -1,    34,     5,
       5,     5,    -1,    35,     3,    63,    -1,    37,     3,    63,
      -1,    38,     3,     3,     3,    63,    -1,    39,    -1,    40,
      -1,    41,     5,     5,    -1,    42,    59,    -1,    43,     5,
       5,     5,    -1,    45,     3,    63,    -1,    44,     3,    63,
      -1,    46,    -1,    47,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   170,   170,   176,   184,   192,   200,   206,   213,   220,
     228,   234,   239,   245,   253,   260,   268,   274,   279,   285,
     293,   299,   312,   318,   323,   331,   336,   342,   351,   357,
     363,   369,   378,   384,   390,   399,   407,   413,   419,   428,
     434,   443,   449,   458,   467,   473,   479,   485,   491,   500,
     509,   515,   521,   530,   536,   545,   554,   563,   569,   575,
     581,   589,   595,   604,   613,   619
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "STRING", "ID", "NUM", "LBRACK",
  "RBRACK", "ACCELERATOR", "ACTIVETRANSFORM", "ALL", "AREALIGHTSOURCE",
  "ATTRIBUTEBEGIN", "ATTRIBUTEEND", "CAMERA", "CONCATTRANSFORM",
  "COORDINATESYSTEM", "COORDSYSTRANSFORM", "ENDTIME", "FILM", "IDENTITY",
  "LIGHTSOURCE", "LOOKAT", "MAKENAMEDMATERIAL", "MATERIAL",
  "NAMEDMATERIAL", "OBJECTBEGIN", "OBJECTEND", "OBJECTINSTANCE",
  "PIXELFILTER", "RENDERER", "REVERSEORIENTATION", "ROTATE", "SAMPLER",
  "SCALE", "SHAPE", "STARTTIME", "SURFACEINTEGRATOR", "TEXTURE",
  "TRANSFORMBEGIN", "TRANSFORMEND", "TRANSFORMTIMES", "TRANSFORM",
  "TRANSLATE", "VOLUME", "VOLUMEINTEGRATOR", "WORLDBEGIN", "WORLDEND",
  "HIGH_PRECEDENCE", "$accept", "start", "array_init", "string_array_init",
  "num_array_init", "array", "string_array", "single_element_string_array",
  "string_list", "string_list_entry", "num_array",
  "single_element_num_array", "num_list", "num_list_entry", "paramlist",
  "paramlist_init", "paramlist_contents", "paramlist_entry",
  "pbrt_stmt_list", "pbrt_stmt", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    49,    50,    51,    52,    53,    54,    54,    55,    55,
      56,    57,    57,    58,    59,    59,    60,    61,    61,    62,
      63,    64,    65,    65,    66,    67,    67,    68,    68,    68,
      68,    68,    68,    68,    68,    68,    68,    68,    68,    68,
      68,    68,    68,    68,    68,    68,    68,    68,    68,    68,
      68,    68,    68,    68,    68,    68,    68,    68,    68,    68,
      68,    68,    68,    68,    68,    68
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     0,     0,     0,     1,     1,     4,     1,
       2,     2,     1,     2,     4,     1,     2,     2,     1,     2,
       2,     0,     2,     0,     2,     2,     1,     3,     2,     2,
       2,     3,     1,     1,     3,     2,     2,     2,     3,     1,
       3,    10,     3,     3,     2,     2,     1,     2,     3,     3,
       1,     5,     3,     4,     3,     3,     5,     1,     1,     3,
       2,     4,     3,     3,     1,     1
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     0,     0,     0,    32,    33,     0,     3,     0,     0,
       0,    39,     0,     0,     0,     0,     0,     0,    46,     0,
       0,     0,    50,     0,     0,     0,     0,     0,     0,    57,
      58,     0,     3,     0,     0,     0,    64,    65,     0,     2,
      26,    21,    28,    29,    30,    21,    21,     5,    35,    15,
      36,    37,    21,    21,     0,    21,    21,    44,    45,    47,
      21,    21,     0,    21,     0,    21,    21,     0,     0,    60,
       0,    21,    21,     1,    25,    27,    23,    31,    34,     5,
       0,    16,    38,    40,     0,    42,    43,    48,    49,     0,
      52,     0,    54,    55,     0,    59,     0,    63,    62,     3,
      20,    23,     5,    18,    19,     0,     0,    53,    21,    61,
       4,    24,     6,     9,     7,    22,    14,    17,     0,    51,
      56,     4,     0,    10,     0,     4,    12,    13,     0,     8,
      11,     0,     0,    41
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,    38,    47,   122,    80,   111,   112,   113,   125,   123,
      48,    49,   102,    81,    75,    76,   100,   101,    39,    40
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -116
static const yytype_int8 yypact[] =
{
      57,     1,    -5,     4,  -116,  -116,    15,  -116,    17,    20,
      21,  -116,    22,    24,    27,    29,    30,    31,  -116,    32,
      33,    34,  -116,    35,    36,    37,    38,    40,    41,  -116,
    -116,    42,  -116,    43,    46,    47,  -116,  -116,    51,    57,
    -116,  -116,  -116,  -116,  -116,  -116,  -116,    48,  -116,  -116,
    -116,  -116,  -116,  -116,    50,  -116,  -116,  -116,  -116,  -116,
    -116,  -116,    52,  -116,    53,  -116,  -116,    49,    54,  -116,
      55,  -116,  -116,  -116,  -116,  -116,    58,  -116,  -116,  -116,
      70,  -116,  -116,  -116,    88,  -116,  -116,  -116,  -116,   100,
    -116,   101,  -116,  -116,    59,  -116,   102,  -116,  -116,  -116,
    -116,    58,    39,  -116,  -116,   103,   104,  -116,  -116,  -116,
       9,  -116,  -116,  -116,  -116,  -116,  -116,  -116,   105,  -116,
    -116,   106,    60,  -116,   107,   108,  -116,  -116,   109,  -116,
    -116,   111,   112,  -116
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
    -116,  -116,   -61,  -116,  -116,  -116,  -116,  -116,  -116,  -115,
     -32,  -116,  -116,   -76,   -44,  -116,   -48,  -116,  -116,    74
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -6
static const yytype_int16 yytable[] =
{
      69,    77,    78,   103,    41,    42,   126,    45,    82,    83,
     130,    85,    86,    43,    -5,   121,    87,    88,    46,    90,
      50,    92,    93,    51,    52,    53,   117,    97,    98,    54,
      55,    44,    56,    57,    58,    59,    60,    61,   110,    63,
      62,    65,    64,    66,    67,   103,   116,    68,    70,    71,
      72,    73,    94,   115,    79,    84,     0,    89,    91,    95,
      96,    99,   108,   127,   120,     1,     2,   114,     3,     4,
       5,     6,     7,     8,     9,   104,    10,    11,    12,    13,
      14,    15,    16,    17,    18,    19,    20,    21,    22,    23,
      24,    25,    26,   105,    27,    28,    29,    30,    31,    32,
      33,    34,    35,    36,    37,   106,   107,   109,   118,   119,
     124,    -5,   128,    74,   131,   129,   132,   133
};

static const yytype_int8 yycheck[] =
{
      32,    45,    46,    79,     3,    10,   121,     3,    52,    53,
     125,    55,    56,    18,     5,     6,    60,    61,     3,    63,
       3,    65,    66,     3,     3,     3,   102,    71,    72,     5,
       3,    36,     3,     3,     3,     3,     3,     3,    99,     3,
       5,     3,     5,     3,     3,   121,     7,     5,     5,     3,
       3,     0,     3,   101,     6,     5,    -1,     5,     5,     5,
       5,     3,     3,     3,   108,     8,     9,    99,    11,    12,
      13,    14,    15,    16,    17,     5,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    29,    30,    31,    32,
      33,    34,    35,     5,    37,    38,    39,    40,    41,    42,
      43,    44,    45,    46,    47,     5,     5,     5,     5,     5,
       5,     5,     5,    39,     5,     7,     5,     5
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     8,     9,    11,    12,    13,    14,    15,    16,    17,
      19,    20,    21,    22,    23,    24,    25,    26,    27,    28,
      29,    30,    31,    32,    33,    34,    35,    37,    38,    39,
      40,    41,    42,    43,    44,    45,    46,    47,    50,    67,
      68,     3,    10,    18,    36,     3,     3,    51,    59,    60,
       3,     3,     3,     3,     5,     3,     3,     3,     3,     3,
       3,     3,     5,     3,     5,     3,     3,     3,     5,    59,
       5,     3,     3,     0,    68,    63,    64,    63,    63,     6,
      53,    62,    63,    63,     5,    63,    63,    63,    63,     5,
      63,     5,    63,    63,     3,     5,     5,    63,    63,     3,
      65,    66,    61,    62,     5,     5,     5,     5,     3,     5,
      51,    54,    55,    56,    59,    65,     7,    62,     5,     5,
      63,     6,    52,    58,     5,    57,    58,     3,     5,     7,
      58,     5,     5,     5
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK (1);						\
      goto yybackup;						\
    }								\
  else								\
    {								\
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (YYID (0))
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
	      (Loc).first_line, (Loc).first_column,	\
	      (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
	break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
#else
static void
yy_stack_print (yybottom, yytop)
    yytype_int16 *yybottom;
    yytype_int16 *yytop;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule)
#else
static void
yy_reduce_print (yyvsp, yyrule)
    YYSTYPE *yyvsp;
    int yyrule;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into YYRESULT an error message about the unexpected token
   YYCHAR while in state YYSTATE.  Return the number of bytes copied,
   including the terminating null byte.  If YYRESULT is null, do not
   copy anything; just return the number of bytes that would be
   copied.  As a special case, return 0 if an ordinary "syntax error"
   message will do.  Return YYSIZE_MAXIMUM if overflow occurs during
   size calculation.  */
static YYSIZE_T
yysyntax_error (char *yyresult, int yystate, int yychar)
{
  int yyn = yypact[yystate];

  if (! (YYPACT_NINF < yyn && yyn <= YYLAST))
    return 0;
  else
    {
      int yytype = YYTRANSLATE (yychar);
      YYSIZE_T yysize0 = yytnamerr (0, yytname[yytype]);
      YYSIZE_T yysize = yysize0;
      YYSIZE_T yysize1;
      int yysize_overflow = 0;
      enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
      char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
      int yyx;

# if 0
      /* This is so xgettext sees the translatable formats that are
	 constructed on the fly.  */
      YY_("syntax error, unexpected %s");
      YY_("syntax error, unexpected %s, expecting %s");
      YY_("syntax error, unexpected %s, expecting %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
# endif
      char *yyfmt;
      char const *yyf;
      static char const yyunexpected[] = "syntax error, unexpected %s";
      static char const yyexpecting[] = ", expecting %s";
      static char const yyor[] = " or %s";
      char yyformat[sizeof yyunexpected
		    + sizeof yyexpecting - 1
		    + ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
		       * (sizeof yyor - 1))];
      char const *yyprefix = yyexpecting;

      /* Start YYX at -YYN if negative to avoid negative indexes in
	 YYCHECK.  */
      int yyxbegin = yyn < 0 ? -yyn : 0;

      /* Stay within bounds of both yycheck and yytname.  */
      int yychecklim = YYLAST - yyn + 1;
      int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
      int yycount = 1;

      yyarg[0] = yytname[yytype];
      yyfmt = yystpcpy (yyformat, yyunexpected);

      for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	  {
	    if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
	      {
		yycount = 1;
		yysize = yysize0;
		yyformat[sizeof yyunexpected - 1] = '\0';
		break;
	      }
	    yyarg[yycount++] = yytname[yyx];
	    yysize1 = yysize + yytnamerr (0, yytname[yyx]);
	    yysize_overflow |= (yysize1 < yysize);
	    yysize = yysize1;
	    yyfmt = yystpcpy (yyfmt, yyprefix);
	    yyprefix = yyor;
	  }

      yyf = YY_(yyformat);
      yysize1 = yysize + yystrlen (yyf);
      yysize_overflow |= (yysize1 < yysize);
      yysize = yysize1;

      if (yysize_overflow)
	return YYSIZE_MAXIMUM;

      if (yyresult)
	{
	  /* Avoid sprintf, as that infringes on the user's name space.
	     Don't have undefined behavior even if the translation
	     produced a string with the wrong number of "%s"s.  */
	  char *yyp = yyresult;
	  int yyi = 0;
	  while ((*yyp = *yyf) != '\0')
	    {
	      if (*yyp == '%' && yyf[1] == 's' && yyi < yycount)
		{
		  yyp += yytnamerr (yyp, yyarg[yyi++]);
		  yyf += 2;
		}
	      else
		{
		  yyp++;
		  yyf++;
		}
	    }
	}
      return yysize;
    }
}
#endif /* YYERROR_VERBOSE */


/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  YYUSE (yyvaluep);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {

      default:
	break;
    }
}

/* Prevent warnings from -Wmissing-prototypes.  */
#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */


/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;



/*-------------------------.
| yyparse or yypush_parse.  |
`-------------------------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{


    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       `yyss': related to states.
       `yyvs': related to semantic values.

       Refer to the stacks thru separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yytoken = 0;
  yyss = yyssa;
  yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */
  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;

	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),
		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss_alloc, yyss);
	YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  *++yyvsp = yylval;

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:

/* Line 1455 of yacc.c  */
#line 171 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
;}
    break;

  case 3:

/* Line 1455 of yacc.c  */
#line 177 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    if (cur_array) Severe("MUH");
    cur_array = new ParamArray;
;}
    break;

  case 4:

/* Line 1455 of yacc.c  */
#line 185 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    cur_array->element_size = sizeof(const char *);
    cur_array->isString = true;
;}
    break;

  case 5:

/* Line 1455 of yacc.c  */
#line 193 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    cur_array->element_size = sizeof(float);
    cur_array->isString = false;
;}
    break;

  case 6:

/* Line 1455 of yacc.c  */
#line 201 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    (yyval.ribarray) = (yyvsp[(1) - (1)].ribarray);
;}
    break;

  case 7:

/* Line 1455 of yacc.c  */
#line 207 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    (yyval.ribarray) = (yyvsp[(1) - (1)].ribarray);
;}
    break;

  case 8:

/* Line 1455 of yacc.c  */
#line 214 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    (yyval.ribarray) = cur_array;
    cur_array = NULL;
;}
    break;

  case 9:

/* Line 1455 of yacc.c  */
#line 221 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    (yyval.ribarray) = cur_array;
    cur_array = NULL;
;}
    break;

  case 10:

/* Line 1455 of yacc.c  */
#line 229 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
;}
    break;

  case 11:

/* Line 1455 of yacc.c  */
#line 235 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
;}
    break;

  case 12:

/* Line 1455 of yacc.c  */
#line 240 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
;}
    break;

  case 13:

/* Line 1455 of yacc.c  */
#line 246 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    char *to_add = strdup((yyvsp[(2) - (2)].string));
    AddArrayElement(&to_add);
;}
    break;

  case 14:

/* Line 1455 of yacc.c  */
#line 254 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    (yyval.ribarray) = cur_array;
    cur_array = NULL;
;}
    break;

  case 15:

/* Line 1455 of yacc.c  */
#line 261 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    (yyval.ribarray) = cur_array;
    cur_array = NULL;
;}
    break;

  case 16:

/* Line 1455 of yacc.c  */
#line 269 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
;}
    break;

  case 17:

/* Line 1455 of yacc.c  */
#line 275 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
;}
    break;

  case 18:

/* Line 1455 of yacc.c  */
#line 280 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
;}
    break;

  case 19:

/* Line 1455 of yacc.c  */
#line 286 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    float to_add = (yyvsp[(2) - (2)].num);
    AddArrayElement(&to_add);
;}
    break;

  case 20:

/* Line 1455 of yacc.c  */
#line 294 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
;}
    break;

  case 21:

/* Line 1455 of yacc.c  */
#line 300 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    for (uint32_t i = 0; i < cur_paramlist.size(); ++i) {
        if (cur_paramlist[i].isString) {
            for (uint32_t j = 0; j < cur_paramlist[i].size; ++j)
                free(((char **)cur_paramlist[i].arg)[j]);
        }
    }
    cur_paramlist.erase(cur_paramlist.begin(), cur_paramlist.end());
;}
    break;

  case 22:

/* Line 1455 of yacc.c  */
#line 313 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
;}
    break;

  case 23:

/* Line 1455 of yacc.c  */
#line 318 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
;}
    break;

  case 24:

/* Line 1455 of yacc.c  */
#line 324 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    cur_paramlist.push_back(ParamListItem((yyvsp[(1) - (2)].string), (yyvsp[(2) - (2)].ribarray)));
    ArrayFree((yyvsp[(2) - (2)].ribarray));
;}
    break;

  case 25:

/* Line 1455 of yacc.c  */
#line 332 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
;}
    break;

  case 26:

/* Line 1455 of yacc.c  */
#line 337 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
;}
    break;

  case 27:

/* Line 1455 of yacc.c  */
#line 343 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    ParamSet params;
    InitParamSet(params, SPECTRUM_REFLECTANCE);
    pbrtAccelerator((yyvsp[(2) - (3)].string), params);
    FreeArgs();
;}
    break;

  case 28:

/* Line 1455 of yacc.c  */
#line 352 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    pbrtActiveTransformAll();
;}
    break;

  case 29:

/* Line 1455 of yacc.c  */
#line 358 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    pbrtActiveTransformEndTime();
;}
    break;

  case 30:

/* Line 1455 of yacc.c  */
#line 364 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    pbrtActiveTransformStartTime();
;}
    break;

  case 31:

/* Line 1455 of yacc.c  */
#line 370 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    ParamSet params;
    InitParamSet(params, SPECTRUM_ILLUMINANT);
    pbrtAreaLightSource((yyvsp[(2) - (3)].string), params);
    FreeArgs();
;}
    break;

  case 32:

/* Line 1455 of yacc.c  */
#line 379 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    pbrtAttributeBegin();
;}
    break;

  case 33:

/* Line 1455 of yacc.c  */
#line 385 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    pbrtAttributeEnd();
;}
    break;

  case 34:

/* Line 1455 of yacc.c  */
#line 391 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    ParamSet params;
    InitParamSet(params, SPECTRUM_REFLECTANCE);
    pbrtCamera((yyvsp[(2) - (3)].string), params);
    FreeArgs();
;}
    break;

  case 35:

/* Line 1455 of yacc.c  */
#line 400 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    if (VerifyArrayLength((yyvsp[(2) - (2)].ribarray), 16, "ConcatTransform"))
        pbrtConcatTransform((float *) (yyvsp[(2) - (2)].ribarray)->array);
    ArrayFree((yyvsp[(2) - (2)].ribarray));
;}
    break;

  case 36:

/* Line 1455 of yacc.c  */
#line 408 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    pbrtCoordinateSystem((yyvsp[(2) - (2)].string));
;}
    break;

  case 37:

/* Line 1455 of yacc.c  */
#line 414 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    pbrtCoordSysTransform((yyvsp[(2) - (2)].string));
;}
    break;

  case 38:

/* Line 1455 of yacc.c  */
#line 420 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    ParamSet params;
    InitParamSet(params, SPECTRUM_REFLECTANCE);
    pbrtFilm((yyvsp[(2) - (3)].string), params);
    FreeArgs();
;}
    break;

  case 39:

/* Line 1455 of yacc.c  */
#line 429 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    pbrtIdentity();
;}
    break;

  case 40:

/* Line 1455 of yacc.c  */
#line 435 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    ParamSet params;
    InitParamSet(params, SPECTRUM_ILLUMINANT);
    pbrtLightSource((yyvsp[(2) - (3)].string), params);
    FreeArgs();
;}
    break;

  case 41:

/* Line 1455 of yacc.c  */
#line 444 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    pbrtLookAt((yyvsp[(2) - (10)].num), (yyvsp[(3) - (10)].num), (yyvsp[(4) - (10)].num), (yyvsp[(5) - (10)].num), (yyvsp[(6) - (10)].num), (yyvsp[(7) - (10)].num), (yyvsp[(8) - (10)].num), (yyvsp[(9) - (10)].num), (yyvsp[(10) - (10)].num));
;}
    break;

  case 42:

/* Line 1455 of yacc.c  */
#line 450 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    ParamSet params;
    InitParamSet(params, SPECTRUM_REFLECTANCE);
    pbrtMakeNamedMaterial((yyvsp[(2) - (3)].string), params);
    FreeArgs();
;}
    break;

  case 43:

/* Line 1455 of yacc.c  */
#line 459 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    ParamSet params;
    InitParamSet(params, SPECTRUM_REFLECTANCE);
    pbrtMaterial((yyvsp[(2) - (3)].string), params);
    FreeArgs();
;}
    break;

  case 44:

/* Line 1455 of yacc.c  */
#line 468 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    pbrtNamedMaterial((yyvsp[(2) - (2)].string));
;}
    break;

  case 45:

/* Line 1455 of yacc.c  */
#line 474 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    pbrtObjectBegin((yyvsp[(2) - (2)].string));
;}
    break;

  case 46:

/* Line 1455 of yacc.c  */
#line 480 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    pbrtObjectEnd();
;}
    break;

  case 47:

/* Line 1455 of yacc.c  */
#line 486 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    pbrtObjectInstance((yyvsp[(2) - (2)].string));
;}
    break;

  case 48:

/* Line 1455 of yacc.c  */
#line 492 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    ParamSet params;
    InitParamSet(params, SPECTRUM_REFLECTANCE);
    pbrtPixelFilter((yyvsp[(2) - (3)].string), params);
    FreeArgs();
;}
    break;

  case 49:

/* Line 1455 of yacc.c  */
#line 501 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    ParamSet params;
    InitParamSet(params, SPECTRUM_REFLECTANCE);
    pbrtRenderer((yyvsp[(2) - (3)].string), params);
    FreeArgs();
;}
    break;

  case 50:

/* Line 1455 of yacc.c  */
#line 510 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    pbrtReverseOrientation();
;}
    break;

  case 51:

/* Line 1455 of yacc.c  */
#line 516 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    pbrtRotate((yyvsp[(2) - (5)].num), (yyvsp[(3) - (5)].num), (yyvsp[(4) - (5)].num), (yyvsp[(5) - (5)].num));
;}
    break;

  case 52:

/* Line 1455 of yacc.c  */
#line 522 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    ParamSet params;
    InitParamSet(params, SPECTRUM_REFLECTANCE);
    pbrtSampler((yyvsp[(2) - (3)].string), params);
    FreeArgs();
;}
    break;

  case 53:

/* Line 1455 of yacc.c  */
#line 531 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    pbrtScale((yyvsp[(2) - (4)].num), (yyvsp[(3) - (4)].num), (yyvsp[(4) - (4)].num));
;}
    break;

  case 54:

/* Line 1455 of yacc.c  */
#line 537 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    ParamSet params;
    InitParamSet(params, SPECTRUM_REFLECTANCE);
    pbrtShape((yyvsp[(2) - (3)].string), params);
    FreeArgs();
;}
    break;

  case 55:

/* Line 1455 of yacc.c  */
#line 546 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    ParamSet params;
    InitParamSet(params, SPECTRUM_REFLECTANCE);
    pbrtSurfaceIntegrator((yyvsp[(2) - (3)].string), params);
    FreeArgs();
;}
    break;

  case 56:

/* Line 1455 of yacc.c  */
#line 555 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    ParamSet params;
    InitParamSet(params, SPECTRUM_REFLECTANCE);
    pbrtTexture((yyvsp[(2) - (5)].string), (yyvsp[(3) - (5)].string), (yyvsp[(4) - (5)].string), params);
    FreeArgs();
;}
    break;

  case 57:

/* Line 1455 of yacc.c  */
#line 564 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    pbrtTransformBegin();
;}
    break;

  case 58:

/* Line 1455 of yacc.c  */
#line 570 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    pbrtTransformEnd();
;}
    break;

  case 59:

/* Line 1455 of yacc.c  */
#line 576 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    pbrtTransformTimes((yyvsp[(2) - (3)].num), (yyvsp[(3) - (3)].num));
;}
    break;

  case 60:

/* Line 1455 of yacc.c  */
#line 582 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    if (VerifyArrayLength( (yyvsp[(2) - (2)].ribarray), 16, "Transform" ))
        pbrtTransform( (float *) (yyvsp[(2) - (2)].ribarray)->array );
    ArrayFree((yyvsp[(2) - (2)].ribarray));
;}
    break;

  case 61:

/* Line 1455 of yacc.c  */
#line 590 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    pbrtTranslate((yyvsp[(2) - (4)].num), (yyvsp[(3) - (4)].num), (yyvsp[(4) - (4)].num));
;}
    break;

  case 62:

/* Line 1455 of yacc.c  */
#line 596 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    ParamSet params;
    InitParamSet(params, SPECTRUM_REFLECTANCE);
    pbrtVolumeIntegrator((yyvsp[(2) - (3)].string), params);
    FreeArgs();
;}
    break;

  case 63:

/* Line 1455 of yacc.c  */
#line 605 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    ParamSet params;
    InitParamSet(params, SPECTRUM_REFLECTANCE);
    pbrtVolume((yyvsp[(2) - (3)].string), params);
    FreeArgs();
;}
    break;

  case 64:

/* Line 1455 of yacc.c  */
#line 614 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    pbrtWorldBegin();
;}
    break;

  case 65:

/* Line 1455 of yacc.c  */
#line 620 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"
    {
    pbrtWorldEnd();
;}
    break;



/* Line 1455 of yacc.c  */
#line 2237 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.cpp"
      default: break;
    }
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
      {
	YYSIZE_T yysize = yysyntax_error (0, yystate, yychar);
	if (yymsg_alloc < yysize && yymsg_alloc < YYSTACK_ALLOC_MAXIMUM)
	  {
	    YYSIZE_T yyalloc = 2 * yysize;
	    if (! (yysize <= yyalloc && yyalloc <= YYSTACK_ALLOC_MAXIMUM))
	      yyalloc = YYSTACK_ALLOC_MAXIMUM;
	    if (yymsg != yymsgbuf)
	      YYSTACK_FREE (yymsg);
	    yymsg = (char *) YYSTACK_ALLOC (yyalloc);
	    if (yymsg)
	      yymsg_alloc = yyalloc;
	    else
	      {
		yymsg = yymsgbuf;
		yymsg_alloc = sizeof yymsgbuf;
	      }
	  }

	if (0 < yysize && yysize <= yymsg_alloc)
	  {
	    (void) yysyntax_error (yymsg, yystate, yychar);
	    yyerror (yymsg);
	  }
	else
	  {
	    yyerror (YY_("syntax error"));
	    if (yysize != 0)
	      goto yyexhaustedlab;
	  }
      }
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (yychar == YYEOF)
	    YYABORT;
	}
      else
	{
	  yydestruct ("Error: discarding",
		      yytoken, &yylval);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;


      yydestruct ("Error: popping",
		  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  *++yyvsp = yylval;


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#if !defined(yyoverflow) || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
     yydestruct ("Cleanup: discarding lookahead",
		 yytoken, &yylval);
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}



/* Line 1675 of yacc.c  */
#line 625 "D:/graphics/nprojects/safegi/final_pack/SafePBRT/typed/src/pbrt/core/pbrtparse.yy"

static const char *paramTypeToName(int type) {
    switch (type) {
    case PARAM_TYPE_INT: return "int";
    case PARAM_TYPE_BOOL: return "bool";
    case PARAM_TYPE_FLOAT: return "float";
    case PARAM_TYPE_POINT: return "point";
    case PARAM_TYPE_VECTOR: return "vector";
    case PARAM_TYPE_NORMAL: return "normal";
    case PARAM_TYPE_RGB: return "rgb/color";
    case PARAM_TYPE_XYZ: return "xyz";
    case PARAM_TYPE_BLACKBODY: return "blackbody";
    case PARAM_TYPE_SPECTRUM: return "spectrum";
    case PARAM_TYPE_STRING: return "string";
    case PARAM_TYPE_TEXTURE: return "texture";
    default: Severe("Error in paramTypeToName"); return NULL;
    }
}


static void InitParamSet(ParamSet &ps, SpectrumType type) {
    ps.Clear();
    for (uint32_t i = 0; i < cur_paramlist.size(); ++i) {
        int type;
        string name;
        if (lookupType(cur_paramlist[i].name, &type, name)) {
            if (type == PARAM_TYPE_TEXTURE || type == PARAM_TYPE_STRING ||
                type == PARAM_TYPE_BOOL) {
                if (!cur_paramlist[i].isString) {
                    Error("Expected string parameter value for parameter \"%s\" with type \"%s\". Ignoring.",
                          name.c_str(), paramTypeToName(type));
                    continue;
                }
            }
            else if (type != PARAM_TYPE_SPECTRUM) { /* spectrum can be either... */
                if (cur_paramlist[i].isString) {
                    Error("Expected numeric parameter value for parameter \"%s\" with type \"%s\".  Ignoring.",
                          name.c_str(), paramTypeToName(type));
                    continue;
                }
            }
            void *data = cur_paramlist[i].arg;
            int nItems = cur_paramlist[i].size;
            if (type == PARAM_TYPE_INT) {
                // parser doesn't handle ints, so convert from floats here....
                int nAlloc = nItems;
                int *idata = new int[nAlloc];
                float *fdata = (float *)cur_paramlist[i].arg;
                for (int j = 0; j < nAlloc; ++j)
                    idata[j] = int(fdata[j]);
                ps.AddInt(name, idata, nItems);
                delete[] idata;
            }
            else if (type == PARAM_TYPE_BOOL) {
                // strings -> bools
                int nAlloc = cur_paramlist[i].size;
                bool *bdata = new bool[nAlloc];
                for (int j = 0; j < nAlloc; ++j) {
                    string s(((const char **)data)[j]);
                    if (s == "true") bdata[j] = true;
                    else if (s == "false") bdata[j] = false;
                    else {
                        Warning("Value \"%s\" unknown for boolean parameter \"%s\"."
                            "Using \"false\".", s.c_str(), cur_paramlist[i].name);
                        bdata[j] = false;
                    }
                }
                ps.AddBool(name, bdata, nItems);
                delete[] bdata;
            }
            else if (type == PARAM_TYPE_FLOAT) {
                ps.AddFloat(name, (float *)data, nItems);
            } else if (type == PARAM_TYPE_POINT) {
                if ((nItems % 3) != 0)
                    Warning("Excess values given with point parameter \"%s\". "
                            "Ignoring last %d of them", cur_paramlist[i].name, nItems % 3);
                ps.AddPoint(name, (tuple3 *)data, nItems / 3);
            } else if (type == PARAM_TYPE_VECTOR) {
                if ((nItems % 3) != 0)
                    Warning("Excess values given with vector parameter \"%s\". "
                            "Ignoring last %d of them", cur_paramlist[i].name, nItems % 3);
                ps.AddVector(name, (tuple3 *)data, nItems / 3);
            } else if (type == PARAM_TYPE_NORMAL) {
                if ((nItems % 3) != 0)
                    Warning("Excess values given with normal parameter \"%s\". "
                            "Ignoring last %d of them", cur_paramlist[i].name, nItems % 3);
                ps.AddNormal(name, (tuple3 *)data, nItems / 3);
            } else if (type == PARAM_TYPE_RGB) {
                if ((nItems % 3) != 0)
                    Warning("Excess RGB values given with parameter \"%s\". "
                            "Ignoring last %d of them", cur_paramlist[i].name, nItems % 3);
                ps.AddRGBSpectrum(name, (float *)data, nItems);
            } else if (type == PARAM_TYPE_XYZ) {
                if ((nItems % 3) != 0)
                    Warning("Excess XYZ values given with parameter \"%s\". "
                            "Ignoring last %d of them", cur_paramlist[i].name, nItems % 3);
                ps.AddXYZSpectrum(name, (float *)data, nItems);
            } else if (type == PARAM_TYPE_BLACKBODY) {
                if ((nItems % 2) != 0)
                    Warning("Excess value given with blackbody parameter \"%s\". "
                            "Ignoring extra one.", cur_paramlist[i].name);
                ps.AddBlackbodySpectrum(name, (float *)data, nItems);
            } else if (type == PARAM_TYPE_SPECTRUM) {
                if (cur_paramlist[i].isString) {
                    ps.AddSampledSpectrumFiles(name, (const char **)data, nItems);
                }
                else {
                    if ((nItems % 2) != 0)
                        Warning("Non-even number of values given with sampled spectrum "
                                "parameter \"%s\". Ignoring extra.", cur_paramlist[i].name);
                    ps.AddSampledSpectrum(name, (float *)data, nItems);
                }
            } else if (type == PARAM_TYPE_STRING) {
                string *strings = new string[nItems];
                for (int j = 0; j < nItems; ++j)
                    strings[j] = string(((const char **)data)[j]);
                ps.AddString(name, strings, nItems);
                delete[] strings;
            }
            else if (type == PARAM_TYPE_TEXTURE) {
                if (nItems == 1) {
                    string val(*((const char **)data));
                    ps.AddTexture(name, val);
                }
                else
                    Error("Only one string allowed for \"texture\" parameter \"%s\"",
                        name.c_str());
            }
        }
        else
            Warning("Type of parameter \"%s\" is unknown",
                cur_paramlist[i].name);
    }
}


static bool lookupType(const char *name, int *type, string &sname) {
    Assert(name != NULL);
    *type = 0;
    const char *strp = name;
    while (*strp && isspace(*strp))
        ++strp;
    if (!*strp) {
        Error("Parameter \"%s\" doesn't have a type declaration?!", name);
        return false;
    }
#define TRY_DECODING_TYPE(name, mask) \
        if (strncmp(name, strp, strlen(name)) == 0) { \
            *type = mask; strp += strlen(name); \
        }
         TRY_DECODING_TYPE("float",     PARAM_TYPE_FLOAT)
    else TRY_DECODING_TYPE("integer",   PARAM_TYPE_INT)
    else TRY_DECODING_TYPE("bool",      PARAM_TYPE_BOOL)
    else TRY_DECODING_TYPE("point",     PARAM_TYPE_POINT)
    else TRY_DECODING_TYPE("vector",    PARAM_TYPE_VECTOR)
    else TRY_DECODING_TYPE("normal",    PARAM_TYPE_NORMAL)
    else TRY_DECODING_TYPE("string",    PARAM_TYPE_STRING)
    else TRY_DECODING_TYPE("texture",   PARAM_TYPE_TEXTURE)
    else TRY_DECODING_TYPE("color",     PARAM_TYPE_RGB)
    else TRY_DECODING_TYPE("rgb",       PARAM_TYPE_RGB)
    else TRY_DECODING_TYPE("xyz",       PARAM_TYPE_XYZ)
    else TRY_DECODING_TYPE("blackbody", PARAM_TYPE_BLACKBODY)
    else TRY_DECODING_TYPE("spectrum",  PARAM_TYPE_SPECTRUM)
    else {
        Error("Unable to decode type for name \"%s\"", name);
        return false;
    }
    while (*strp && isspace(*strp))
        ++strp;
    sname = string(strp);
    return true;
}



