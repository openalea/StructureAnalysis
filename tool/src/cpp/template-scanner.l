/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2000 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): Ch. Godin (christophe.godin@cirad.fr)
 *
 *       $Source$
 *       $Id: template-scanner.l 3254 2007-06-06 08:27:11Z pradal $
 *
 *       Forum for AMAPmod developers    : amldevlp@cirad.fr
 *
 *  ----------------------------------------------------------------------------
 *
 *                      GNU General Public Licence
 *
 *       This program is free software; you can redistribute it and/or
 *       modify it under the terms of the GNU General Public License as
 *       published by the Free Software Foundation; either version 2 of
 *       the License, or (at your option) any later version.
 *
 *       This program is distributed in the hope that it will be useful,
 *       but WITHOUT ANY WARRANTY; without even the implied warranty of
 *       MERCHANTABILITY or FITNESS For A PARTICULAR PURPOSE. See the
 *       GNU General Public License for more details.
 *
 *       You should have received a copy of the GNU General Public
 *       License along with this program; see the file COPYING. If not,
 *       write to the Free Software Foundation, Inc., 59
 *       Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 *  ----------------------------------------------------------------------------
 */


%{ //-*-c++-*-


  // flags that control the debugging

#if 0
#define LEXDEBUG
#if 1
#define PRINT_LINE_AND_COLUMN
#endif
#endif

#include <fstream.h>
#include "gscanner.h"           // to access to the generic scanner
#ifdef BISON_HPP
#include "template-parser.hpp"          // definition of the parser tokens
#else
#include "template-parser.cpp.h"
#endif


%}

/*
   %option c++ is equivalent to compiling with flex++.
   %option yyclass only applies when generating a C++ scanner (-+
   option).  It informs flex that you have derived foo as a subclass
   of yyFlexLexer, so flex will place your actions in the member func­
   tion foo::yylex() instead of yyFlexLexer::yylex().  It also
   generates a yyFlexLexer::yylex() member function that emits a
   run-time error (by invoking yyFlexLexer::LexerError()) if called.

   %pointer yytext grows dynamically to accommodate large tokens.
    While this means your %pointer scanner can accommodate very
    large tokens (such as matching entire blocks of com­ ments),
    bear in mind that each time the scanner must resize yytext it
    also must rescan the entire token from the beginning, so
    matching such tokens can prove slow.


*/

%option c++
%option yyclass="GENERIC_LEXER"
%option yylineno
%pointer

/* lex states */

%x include
%x comment

S       [a-zA-Z_]([0-9]|[a-zA-Z_])*
B       [ ]+
T       [\t]+
D       [0-9]+
R       ({D}"."{D}?)|({D}?"."{D})|({D}("."{D}?)?[Ee][+-]?{D})|("."{D}[Ee][+-]?{D})
wsstr   ([^=].*)?
quoted_wsstr  \"[^\"]*\"
quoted_str  \"[^\"\t\n]*\"

%%

 /* yyterminate() could be used directly within the lexer if
    used in stand alone mode instead of return TokQUIT */

:[qQ]uit        {TRACE;return TokQUIT;}
:[eE]xit        {TRACE;return TokQUIT;}
:[bB]ye         {TRACE;return TokQUIT;}

:[cC]lear       {TRACE; return TokCLEAR;}
:[pP]rint       {TRACE; return TokPRINT;}
:[iI]nclude     {TRACE; BEGIN(include);}

!{wsstr}$       {TRACE;
                  int result = ::system(YYText()+1);
                  if (result == -1) lostream << "scanner: command `system' Failed" << endl;
                }
\?              {TRACE; return TokINFO;}

{B}             {TRACE;}
{T}             {_columno += strlen(YYText());TRACE;}

\(              {TRACE; return '(';}
\)              {TRACE; return ')';}
\,              {TRACE; return ',';}
\;              {TRACE; return ';';}
{               {TRACE; return '{';}
}               {TRACE; return '}';}
=               {TRACE; return '=';}
\-              {TRACE; return TokMINUS;}
\+              {TRACE; return TokPLUS;}
\*              {TRACE; return TokTIMES;}
"/"             {TRACE; return TokSLASH;}

\'.\'           {TRACE; VAL->character = YYText()[1]; return TokCharacter;}
{S}             {TRACE;
                 // VAL is a pointer to th union YYSTYPE defined by the parser
                 // The user must make sure that the allocated string will be released later
                 // outside the parser
                 VAL->string = new char[strlen(YYText())+1];
                 strcpy(VAL->string, YYText());
                 return TokName;}
{quoted_wsstr}  {TRACE;
                 VAL->string = new char[strlen(YYText())-1];
                 strncpy(VAL->string, YYText()+1,strlen(YYText())-2);
                 VAL->string[strlen(YYText())-2] = '\0';
                 return TokWSString;}
{quoted_str}    {TRACE;
                 VAL->string = new char[strlen(YYText())-11];
                 strncpy(VAL->string, YYText()+1,strlen(YYText())-2);
                 VAL->string[strlen(YYText())-2] = '\0';
                 return TokString;}
{D}             {TRACE; VAL->integer = atoi(YYText()); return TokInt;}
{R}             {TRACE; VAL->real = atof(YYText()); return TokReal;}

<include>{B}+
<include>{T}+   {_columno += strlen(YYText());TRACE;}
<include>[^ \t\n]+ {TRACE; // reading the filename
                    istream* fs = new ifstream(YYText());
                    if (!*fs) {
                      delete fs;
                      lostream << "scanner: error: cannot include file: " << YYText() << endl;
                      lostream << "scanner: error: skipping the include command. " << endl;
                    }
                    else
                    {
                      if (!pushStream(*fs,YYText())) {
                        lostream << "error: Includes nested too deeply (max depth = "
                                 << MAX_INCLUDE_DEPTH << ")" << endl;
                      }
                    }
                    BEGIN(INITIAL);
                   }

"(#"              {Trace("comment");BEGIN(comment); /* (# commentaire a la C++ #) */}
"#".*$            {Trace("comment");/* # commentaire de ligne a la C++ */}
<comment>[^#\n]*        /* eat anything that's not a '#' */
<comment>"#"+[^#)\n]*   /* eat up '#'s not followed by ')'s */
<comment>\n
<comment>"#"+")"  { BEGIN(INITIAL); }
<comment><<EOF>>  { lostream << "scanner: unexpected end in comment" << endl; if (!popStream()) yyterminate();}

\n              {_columno=1;TRACE;}
<<EOF>>         {Trace("EOF");if (!popStream()) yyterminate();}
.               {Trace("unmatched");
#ifdef LEXDEBUG
                 printChar(*YYText());
#endif
                 return TokError;}

%%

