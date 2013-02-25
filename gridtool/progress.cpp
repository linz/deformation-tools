#include <iostream>
#include <string.h>

using namespace std;

#include "progress.hpp"


#define HEADER_STRING  ""
#define TAIL_STRING    ""

#if defined(__MSDOS__) || defined(__WATCOMC__)
#define SPACE_CHAR ((unsigned char) 0xB0)
#define FILL_CHAR  ((unsigned char) 0xDB)
#else
#define SPACE_CHAR '-'
#define FILL_CHAR '*'
#endif



void ProgressMeter::Start( char *status, long maxValue, long value )
{
    nesting++;
    if( nesting > 1 ) return;
    size = maxValue;
    curvalue = 0;
    display = 0;
    CalcDisplayValue( value );
    Show( status );
}

short ProgressMeter::Update( long newValue )
{
    if( checkAbort && Abort() ) return 0;
    if( nesting != 1 ) return 1;
    if( CalcDisplayValue(newValue)) Redisplay();
    return 1;
}

void ProgressMeter::Finish()
{
    nesting--;
    if( nesting != 0 ) return;
    Hide();
}

short ProgressMeter::CalcDisplayValue( long newValue )
{
    short newDisplay;
    newDisplay = (short) (newValue * resolution / size);
    if( newDisplay < 0 ) newDisplay = 0;
    if( newDisplay > resolution ) newDisplay = resolution;
    if( newDisplay )
    {
        display = newDisplay;
        return 1;
    }
    return 0;
}

AsciiBarMeter::AsciiBarMeter( char *Prefix, short size )
{
    if( ! Prefix ) Prefix = "";
    prefix = new char[strlen(Prefix)+1];
    strcpy( prefix, Prefix );
    if( size > 0 ) resolution = size;
}

AsciiBarMeter::~AsciiBarMeter()
{
    if( prefix ) delete [] prefix;
}

void AsciiBarMeter::Show( char *status )
{
    if( status ) cout << prefix << status << "\n";
    cout << prefix;
    cout << HEADER_STRING;
    for( short i = 0; i < resolution; i++ ) cout << SPACE_CHAR;
    cout << TAIL_STRING;
    cout << '\r';
    cout << prefix;
    cout << HEADER_STRING;
    cout.flush();
    curPos = 0;
    Redisplay();
}


void AsciiBarMeter::Redisplay()
{
    while( curPos < display )
    {
        cout << FILL_CHAR;
        curPos++;
    }
    cout.flush();
}

void AsciiBarMeter::Hide()
{
    cout << '\r';
    for( short i = strlen(prefix) + strlen( HEADER_STRING ) + strlen( TAIL_STRING) + resolution;
            i--;
       ) cout << ' ';
    cout << '\r';
    cout.flush();
}

