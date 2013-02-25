#ifndef _PROGRESS_HPP
#define _PROGRESS_HPP

class ProgressMeter
{
private:
    short CalcDisplayValue( long newValue );
protected:
    long size;
    long curvalue;
    short resolution;
    short display;
    short nesting;
    short showing;
    short checkAbort;
    virtual void Show( char *status ) = 0;
    virtual void Redisplay() = 0;
    virtual void Hide() = 0;
    virtual short Abort()
    {
        return 0;
    }
public:
    ProgressMeter() :
        size(0),
        nesting(0),
        resolution(100),
        showing(0),
        checkAbort(0) {}
    virtual ~ProgressMeter() {}
    void Start(char *status, long range, long value = 0);
    short Update(long value);
    void Finish();
};


class AsciiBarMeter : public ProgressMeter
{
    short curPos;
    char *prefix;
protected:
    virtual void Show( char *status );
    virtual void Redisplay();
    virtual void Hide();
public:
    AsciiBarMeter( char *prefix = 0, short size = 50 );
    ~AsciiBarMeter();
};

#endif
