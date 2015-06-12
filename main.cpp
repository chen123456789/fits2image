/*
 * This program converts FITS images to JPG and/or PNGs.
 *
 * Copyright (c) 2010-2013, www.cyberska.org. All rights reserved.
 * Written by: Pavol Federl, federl@gmail.com
 *
 */

#include <iostream>
#include <limits>
#include <vector>
#include <utility>
#include <cmath>
#include <QApplication>
#include <QString>
#include <QStringList>
#include <QFile>
#include <QRect>
#include <QImage>
#include <QColor>
#include <QVariant>
#include <QFileInfo>

namespace fits2image {

static const double NaNd = std::numeric_limits<double>::quiet_NaN();

// provide convenience overload for std::ostream << QString
std::basic_ostream<char> & operator << ( std::basic_ostream<char> & out, const QString & s)
{
    out << s.toStdString();
    return out;
}


template <typename T>
inline
T clamp( const T & v, const T & v1, const T & v2)
{
    if( v < v1) return v1;
    if( v > v2) return v2;
    return v;
}

template <typename T>
inline void
swap_ordered( T & v1, T & v2)
{
    if( v1 > v2) std::swap( v1, v2);
}

// simple Rgb class (8 bit per color, 0..255)
struct Rgb {
    int r, g, b;
    Rgb() { r = g = b = 0; }
    Rgb( int red, int green, int blue) { r = red; g = green; b = blue; }

    QString toString() { return QString("(%1,%2,%3)").arg (r).arg (g).arg (b); }
};


class PWLinear {
public:
    PWLinear & add( double x, double y);
    double operator()(double x);
private:
    QList<QPointF> points_;
};

class ColormapFunction {
public:
    virtual Rgb operator()( double x);
    ColormapFunction( PWLinear & red, PWLinear & green, PWLinear & blue);
    void setInvert( bool);
    void setReverse( bool);
    void setRgbScales( double r, double g, double b);
    void setScales( double s1, double s2);
    ColormapFunction();

    // predefined colormaps
    static ColormapFunction gray();
    static ColormapFunction heat();
    static ColormapFunction fire();
    static ColormapFunction spring();
    static ColormapFunction sea();
    static ColormapFunction sunbow();
    static ColormapFunction mutant();
    static ColormapFunction aberration();
    static ColormapFunction rgb();
    static ColormapFunction velocity();
    static ColormapFunction cubeHelix( double start, double rots, double hue, double gamma);


private:
    PWLinear red_, green_, blue_;
    bool invert_, reverse_;
    double redScale_, greenScale_, blueScale_;
    double scale1_, scale2_;
};

class HistFunctor
{
public:
    HistFunctor( double min, double max);
    Rgb operator() ( double v) const ;
private:
    double min_, max_;
    double la_, lb_;
};

// combines histogram & colormap functions into one
class HistogramColormapFunctor
{
public:
    HistogramColormapFunctor( double min, double max, ColormapFunction cmap);
    HistogramColormapFunctor();
    Rgb operator() (double x);

private:
    PWLinear histogramFunction_;
    ColormapFunction colormapFunction_;
};

// caches a double -> Rgb function by precomputing values
/**
 * precalculates values for a function to speed up computation
 *
 * function is (x) --> y
 *
 * Parameters: {
 *   fn: function(x) --> y
 *   min: minimum value for x
 *   max: maximum value for y
 *   n: number of intervals to cache
 *   nanColor: what to return if x == NaN
 * }
 */
class CachedRgbFunction {
public:
    CachedRgbFunction();
    CachedRgbFunction( HistogramColormapFunctor fn, double min, double max, int n, Rgb nanColor);
    Rgb operator()( double x);

private:
    QVector<Rgb> cache;
    double min_, max_;
    Rgb nanColor_;
    double diff, dInvN1;
    Rgb minRes, maxRes;
    int n1;

    void init( HistogramColormapFunctor fn, double min, double max, int n, Rgb nanColor);

};

inline Rgb
CachedRgbFunction::operator () (double x)
{
    if( ! std::isfinite(x))
        return nanColor_;

    if( diff <= 0) {
        if( x < min_) return minRes; else return maxRes;
    }
    int ind = round( (x-min_)*dInvN1);
    if( ind < 0) return minRes;
    if( ind > n1) return maxRes;
    return cache[ind];

}


PWLinear & PWLinear::add( double x, double y) {
    points_.push_back ( QPointF(x,y));
    return * this; // allow chaining
}

double PWLinear::operator()(double x) {
    if( ! std::isfinite(x) || points_.empty ())
        return std::numeric_limits<double>::quiet_NaN();
    // test boundary conditions
    if( x <= points_.first ().x ())
        return points_.first ().y ();
    if( x >= points_.last ().x ())
        return points_.last ().y ();
    // find the segment and interpolate within it
    for( int i = 1 ; i < points_.size() ; i ++ ) {
        if( x <= points_[i].x ()) {
            double a = (points_[i-1].y() - points_[i].y()) / (points_[i-1].x() - points_[i].x());
            double b = points_[i].y() - a * points_[i].x();
            return a * x + b;
        }
    }
    //    dbg(0) << ConsoleColors::error ()
    //           << "Weird. PWLinear could not find segment for x = " << x
    //           << ConsoleColors::resetln ();
    return std::numeric_limits<double>::quiet_NaN();
}

ColormapFunction::ColormapFunction()
{
    invert_ = false;
    reverse_ = false;
    redScale_ = greenScale_ = blueScale_ = 1.0;
    scale1_ = scale2_ = 0.0;
}
ColormapFunction::ColormapFunction(
        PWLinear & red, PWLinear & green, PWLinear & blue)
{
    red_ = red; green_ = green; blue_ = blue;
    invert_ = false;
    reverse_ = false;
    redScale_ = greenScale_ = blueScale_ = 1.0;
}
Rgb ColormapFunction::operator()( double x) {
    if( ! std::isfinite( x))
        return Rgb(0,0,0);

    x = x - 0.5 * scale1_;
    x = (x-0.5) * exp( - 3 * scale2_) + 0.5;

    if( reverse_) x = 1 - x;
    Rgb rgb( 255 * red_(x), 255 * green_(x), 255 * blue_(x));
    if( invert_) {
        rgb.r = 255 - rgb.r;
        rgb.g = 255 - rgb.g;
        rgb.b = 255 - rgb.b;
    }
    rgb.r = round ( rgb.r * redScale_);
    rgb.g = round ( rgb.g * greenScale_);
    rgb.b = round ( rgb.b * blueScale_);
    return rgb;
}
void ColormapFunction::setInvert (bool on) {
    invert_ = on;
}
void ColormapFunction::setReverse (bool on) {
    reverse_ = on;
}
void ColormapFunction::setRgbScales (double r, double g, double b) {
    redScale_ = r; greenScale_ = g; blueScale_ = b;
}

void ColormapFunction::setScales(double s1, double s2)
{
    scale1_ = clamp<double>( s1, -1.0, 1.0);
    scale2_ = clamp<double>( s2, -1.0, 1.0);
}

HistFunctor::HistFunctor( double min, double max)
{
    min_ = min; max_ = max;
    la_ = 255 / (max_ - min_);
    lb_ = - 255 * min_ / (max_ - min_) + 0.5;
}
Rgb HistFunctor::operator() ( double v) const {
    if( ! std::isfinite(v))
        return Rgb(0,0,0);
    if( v < min_)
        return Rgb(255,0,0);
    if( v > max_)
        return Rgb(0,255,0);
    int gray = round(v * la_ + lb_);
    if( gray < 0) gray = 0; if( gray > 255) gray = 255;

    return Rgb(gray, gray, gray);
}

HistogramColormapFunctor::HistogramColormapFunctor()
{
    histogramFunction_ = PWLinear().add (0,0.0).add (1,1.0);
    colormapFunction_ = ColormapFunction::gray ();
}

// combines histogram & colormap functions into one
HistogramColormapFunctor::HistogramColormapFunctor(
        double histMin, double histMax, ColormapFunction cmap)
{
    histogramFunction_ = PWLinear().add (histMin,0.0).add (histMax,1.0);
    colormapFunction_ = cmap;
}
Rgb HistogramColormapFunctor::operator() (double x)
{
    return colormapFunction_( histogramFunction_( x));
}

// PREDEFINED COLORMAPS
// --------------------

ColormapFunction
ColormapFunction::gray ()
{
    static ColormapFunction map;
    static bool initialized = false;
    if( ! initialized) {
        initialized = true;
        map = ColormapFunction(
                    PWLinear().add(0,0).add(1,1),
                    PWLinear().add(0,0).add(1,1),
                    PWLinear().add(0,0).add(1,1)
                    );
    }
    return map;
}

ColormapFunction
ColormapFunction::heat ()
{
    static ColormapFunction map;
    static bool initialized = false;
    if( ! initialized) {
        initialized = true;
        map = ColormapFunction(
                    PWLinear().add(0,0).add(1.0/3,1),
                    PWLinear().add(0,0).add(1,1),
                    PWLinear().add(2.0/3,0).add(1,1)
                    );
    }
    return map;
}

ColormapFunction
ColormapFunction::fire ()
{
    static ColormapFunction map;
    static bool initialized = false;
    if( ! initialized) {
        initialized = true;
        map = ColormapFunction(
                    PWLinear().add(0,0).add(1.0/2,1),
                    PWLinear().add(1.0/4,0).add(3.0/4,1),
                    PWLinear().add(1.0/2,0).add(1,1)
                    );
    }
    return map;
}

ColormapFunction
ColormapFunction::spring ()
{
    static ColormapFunction map;
    static bool initialized = false;
    if( ! initialized) {
        initialized = true;
        map = ColormapFunction(
                    PWLinear().add(0.25,0).add(0.75,1),
                    PWLinear().add(0,0).add(0.5,1),
                    PWLinear().add(0.5,0).add(1,1)
                    );
    }
    return map;
}

ColormapFunction
ColormapFunction::sea ()
{
    static ColormapFunction map;
    static bool initialized = false;
    if( ! initialized) {
        initialized = true;
        map = ColormapFunction(
                    PWLinear().add(0.5,0).add(1,1),
                    PWLinear().add(0,0).add(0.5,1),
                    PWLinear().add(0.25,0).add(0.75,1)
                    );
    }
    return map;
}

ColormapFunction
ColormapFunction::sunbow ()
{
    static ColormapFunction map;
    static bool initialized = false;
    if( ! initialized) {
        initialized = true;
        map = ColormapFunction(
                    PWLinear().add(0.5,1).add(0.75,0),
                    PWLinear().add(0,0).add(0.25,1).add (0.75,1).add (1,0),
                    PWLinear().add(0.25,0).add(0.5,1)
                    );
    }
    return map;
}

ColormapFunction
ColormapFunction::mutant ()
{
    static ColormapFunction map;
    static bool initialized = false;
    if( ! initialized) {
        initialized = true;
        map = ColormapFunction(
                    PWLinear().add(0,0).add (0.25,1).add (0.75,1).add (1,0),
                    PWLinear().add(0.5,1).add (0.75,0),
                    PWLinear().add(0.25,0).add(0.5,1)
                    );
    }
    return map;
}

ColormapFunction
ColormapFunction::aberration ()
{
    static ColormapFunction map;
    static bool initialized = false;
    if( ! initialized) {
        initialized = true;
        map = ColormapFunction(
                    PWLinear().add(0.5,0).add (0.75,1),
                    PWLinear().add(0.25,0).add (0.5,1),
                    PWLinear().add(0,0).add(0.25,1).add (0.75,1).add (1,0)
                    );
    }
    return map;
}

ColormapFunction
ColormapFunction::rgb ()
{
    static ColormapFunction map;
    static bool initialized = false;
    if( ! initialized) {
        initialized = true;
        map = ColormapFunction(
                    PWLinear().add(0,0).add (0.25,1).add (0.5,0),
                    PWLinear().add(0.25,0).add (0.5,1).add (0.75,0),
                    PWLinear().add(0.5,0).add (0.75,1).add (1,0)
                    );
    }
    return map;
}

ColormapFunction
ColormapFunction::velocity ()
{
    static ColormapFunction map;
    static bool initialized = false;
    if( ! initialized) {
        initialized = true;
        map = ColormapFunction(
                    PWLinear().add(0.5,0).add (1,1),
                    PWLinear().add(0,0).add (0.5,1).add (1,0),
                    PWLinear().add(0,1).add (0.5,0)
                    );
    }
    return map;
}

CachedRgbFunction::CachedRgbFunction()
{
    init( HistogramColormapFunctor(0,1,ColormapFunction::heat()), 0, 1, 10000, Rgb(0,0,0));
}

void
CachedRgbFunction::init (
        HistogramColormapFunctor fn,
        double min,
        double max,
        int n,
        Rgb nanColor)
{
    min_ = min;
    max_ = max;
    nanColor_ = nanColor;

    cache.resize ( n);
    diff = max - min;
    n1 = n - 1;
    dInvN1 = n1 / diff ;
    double x; int i;
    double delta = diff / n;
    for( i = 0 , x = min ; i < n ; i ++ , x+= delta ) {
        cache[i] = fn( x);
    }
    if( diff <= 0) {
        minRes = fn( min - 1e-9);
        maxRes = fn( max + 1e-9);
    } else {
        minRes = cache.first ();
        maxRes = cache.last ();
    }
}

CachedRgbFunction::CachedRgbFunction(
        HistogramColormapFunctor fn,
        double min,
        double max,
        int n,
        Rgb nanColor)
{
    init( fn, min, max, n, nanColor);
}

// cube helix function
// as described here: http://www.mrao.cam.ac.uk/~dag/CUBEHELIX/
// and in
//  Green, D. A., 2011, Bulletin of the Astronomical Society of India, Vol.39, p.289
ColormapFunction
ColormapFunction::cubeHelix (double start, double rots, double hue, double gamma)
{
    PWLinear red, green, blue;

    int nlev = 1000;
    for( int i = 0 ; i < nlev ; i ++ ) {

        double fract = double(i) / nlev;
        double angle = 2 * M_PI * ( start/3 + 1 + rots + fract);
        fract = pow( fract, gamma);

        double amp = hue * fract * (1-fract)/2.0;

        double r = fract + amp*(-0.14861*cos(angle)+1.78277*sin(angle));
        double g = fract + amp*(-0.29227*cos(angle)-0.90649*sin(angle));
        double b = fract + amp*(+1.97294*cos(angle));

        if( r < 0) r = 0;
        if( r > 1) r = 1;
        if( g < 0) g = 0;
        if( g > 1) g = 1;
        if( b < 0) b = 0;
        if( b > 1) b = 1;

        red.add( fract, r);
        green.add( fract, g);
        blue.add( fract, b);
    }
    return ColormapFunction( red, green, blue);
}


namespace RaiLib {

/// HistogramInfo holds histogram related information about a frame or file
/// i.e. min/max for 100%, 99.99%, etc...
struct HistogramInfo {
    double min, min95, min98, min99, min995, min999, min9999,
    max, max95, max98, max99, max995, max999, max9999;

    /// has this been retrieved from the db
    bool valid;
    /// default constructor - initialize invalid info
    HistogramInfo() { valid = false; }
    /// debugging
    std::string str() const;

    static const quint64 BinSize = sizeof(min) * 14 + sizeof(valid);

    bool operator == ( const HistogramInfo & h) const;
};
};



/// This class is used to do some basic parsing of FITS files. Only the primary HDU
/// is parsed. The primary purpose is to extract some 'interesting' values from
/// the header, and then to access the individual elements of the image array.
class FitsParser
{
public:

    /// @brief Relevant values extracted from the fits header.
    /// Only 'interesting values' are extracted. These are basically the ones needed to
    /// complete the extraction, and ones that the fits viewer can use
    struct HeaderInfo {
        // raw values
        int bitpix;
        int naxis;
        int naxis1, naxis2, naxis3;
        double bscale, bzero;
        double crpix1, crpix2, crpix3, crval1, crval2, crval3, cdelt1, cdelt2, cdelt3;
        QString ctype1, ctype2, ctype3, cunit3, bunit;
        double equinox;
        int blank; bool hasBlank;
        QString bmaj, bmin, bpa;
        int iquset;

        // derived values
        int bitpixSize; // bitpix converted to bytes
        bool scalingRequired; // if BLANK, BZERO or BITPIX are specified
        int totalFrames; // total number of frames = product(naxis3,...,naxisN)
        qint64 dataOffset; // calculated offset to the first byte of the data

        // all header lines stored raw
        QStringList headerLines;

        HeaderInfo();
    };

    /// constructor - creates a non-valid parser
    FitsParser();
    /// destructor - frees up resources (if any)
    ~FitsParser();
    /// Sets the file location, parses the header and returns true on success
    bool loadFile( const QString & fname);
    /// Returns a parsed header info (only interesting information though)
    inline const HeaderInfo & getHeaderInfo() const { return headerInfo_; }
    /// Extracts a value from a pixel (using image coordinates)
    /// BLANK, BZERO and BSCALE are applied to the extracted raw values
    /// and returns the result as double, regardless of bitpix
    double src( int x, int y, int z);
    /// extracts a value from a pixel and returns a pointer to the raw value, adjusted
    /// to the endianness of the local machine, but without applying BLANK,BZERO or BSCALE
    const char * srcRaw( int x, int y, int z);
    /// applies BLANK, BZERO and BSCALE to a raw value (passed in as a pointer)
    double raw2double( const char * ptr);
    /// calculates frame information
    std::pair<double,double> getFrameInfo( int frame, double preset);

    // pre-read an entire frame into internal cache
    void cacheFrame( int z);

    // format a value
    QString formatValue( double val);

    //    // return the cache for the current file
    //    FitsInfoCache & cache() { return fitsCache_; }

private:
    QFile * file_;
    //    FitsInfoCache fitsCache_;
    //    FitsInfoCache::FileInfo fileInfo_;
    //    QVector< RaiLib::HistogramInfo > frameInfos_;
    QString fLocation_;
    HeaderInfo headerInfo_;
    struct M3DBitpixFile;
    M3DBitpixFile * bitPixFile_;
    // storing these as private members helps to avoid memory allocations whenever possible
    //    QList<double> profileX_, profileY_, profileZ_;

    //    std::vector<double> frameCache_;
    std::vector<char> frameCacheRaw_;
    int cachedFrame_;

    // disable copy
    FitsParser( const FitsParser &);
    FitsParser & operator = ( const FitsParser &);

    //    // mutex to make this little more thread safe
    //    static QMutex mutex_;

    //    // casa stuff
    //    class Casa;
    //    Casa * casa_;
};

FitsParser::FitsParser()
{
    file_ = 0;
    bitPixFile_ = 0;
    cachedFrame_ = -1;
}

FitsParser::~FitsParser()
{
    if( file_ != 0) {
        if( file_ -> isOpen ())
            file_ -> close ();
        delete file_;
    }
}

FitsParser::HeaderInfo::HeaderInfo()
{
    bitpix = naxis = naxis1 = naxis2 = naxis3 = 0;
    bscale = bzero = crpix1 = crpix2 = crpix3 = crval1 = crval2 = crval3 = 0;
    cdelt1 = cdelt2 = cdelt3 = 0;
    equinox = 0;
    blank = 0; hasBlank = false;
    totalFrames = 0;
    dataOffset = 0;
    scalingRequired = false;
    iquset = -1;
}

// buffered 3D FITS cube accessor (by index x,y,z)
// almost acts as a 3D matrix but the data is stored on the disk
struct FitsParser::M3DBitpixFile {
    M3DBitpixFile( const FitsParser::HeaderInfo & fits, QFile * fp, qint64 offset, int dx, int dy, int dz) {
        _fp = fp; _offset = offset; _dx = dx; _dy = dy; _dz = dz; _fits = fits;
        _dataSize = _fits.bitpixSize;
        _buffStart = _buffEnd = -1; // indicate invalid buffer
        _buffSize = 4096;
        _buff = new char[ _buffSize];
        if( ! _buff) throw QString( "Could not allocate buffer");
    }
    ~M3DBitpixFile() { delete[] _buff; }

    inline double getDouble(int x, int y, int z) {

        const char * tmp = getRaw ( x, y, z);

        // convert the raw bytes to an initial double
        double original;
        switch( _fits.bitpix) {
        case   8: original = double( * (quint8 *) tmp); break;
        case  16: original = double( * (qint16 *) tmp); break;
        case  32: original = double( * (qint32 *) tmp); break;
        case  64: original = double( * (qint64 *) tmp); break;
        case -32: original = double( * ( float *) tmp); break;
        case -64: original = double( * (double *) tmp); break;
        default: throw QString( "Illegal value BITPIX = %1").arg( _fits.bitpix);
        }

        // check to see if we should apply BLANK
        //        if( _fits.hasBlank && original == _fits.bzero + _fits.blank * _fits.bscale)
        if( _fits.hasBlank && original == _fits.blank )
            return std::numeric_limits<double>::quiet_NaN();
        else
            return _fits.bzero + _fits.bscale * original;
    }


    // place to store results from getRaw
    char resBuffer[8];

    // returns a pointer where the raw bytes can be extracted from, already
    // converted to the native endian
    // on error an exception is thrown
    inline const char * getRaw(int x, int y, int z) {
        qint64 ind = z; ind = ind * _dy + y; ind = ind * _dx + x; ind *= _dataSize; ind += _offset;
        // if buffer does not yet point to what we need, read in stuff into buffer
        if( ind < _buffStart || ind + _dataSize -1 > _buffEnd ) {
            qint64 req = _buffSize;
            if( ! _fp-> seek( ind)) throw QString( "Failed to seek to extract data");
            qint64 len = _fp-> read( _buff, req);
            if( len < _dataSize) {
                throw QString( "failed to read the minimum (%1").arg( _dataSize);
            }
            _buffStart = ind;
            _buffEnd = _buffStart + len - 1;
        }
        char * p = _buff + (ind - _buffStart);
        // endian swap
#if Q_BYTE_ORDER == Q_BIG_ENDIAN
        // no wap, eg. on powerpc
        return p;
#else
        // big to small endian copy. e.g. on intel
        std::reverse_copy( p, p + _dataSize, resBuffer);
        return resBuffer;
#endif
    }


    QFile * _fp; qint64 _offset; int _dx, _dy, _dz, _dataSize;
    char * _buff;
    qint64 _buffStart, _buffEnd, _buffSize;
    FitsParser::HeaderInfo _fits;

};


// given a raw bitpix value (as found in FITS files), returns the number of bytes
// which is basically abs(bitpix)/8
inline int bitpixToSize( int bitpix )
{
    if( bitpix == 8 ) return 1;
    if( bitpix == 16 ) return 2;
    if( bitpix == 32 ) return 4;
    if( bitpix == -32 ) return 4;
    if( bitpix == -64 ) return 8;

    throw QString( "Illegal value BITPIX = %1").arg( bitpix);
}

// FitsLine represents a single entry in the Fits header (I think it's called a card... :)
struct FitsLine {
    FitsLine( const QString & rawLine ) {
        _raw = rawLine;
    }
    QString raw() { return _raw; }
    QString key() { QString k, v, c; parse( k, v, c); return k; }
    QString value() { QString k, v, c; parse( k, v, c); return v; }
    QString comment() { QString k, v, c; parse( k, v, c); return c; }
    // parse the line into key/value/comment
    void parse( QString & key, QString & value, QString & comment ) {
        // key is the first 8 characters (trimmed)
        key = _raw.left(8).trimmed();
        // by default, value & comment are empty
        value = comment = QString();
        // if there is no equal sign present, return the default values for value/comment, which is empty
        if( _raw.mid( 8, 2).trimmed() != "=") return;
        // find the start/end of the value
        //   start = first non-white character
        //   end   = last character of the value (if string, it's the closing quote, otherwise it's the last non-space
        int vStart = 10, vEnd = -1;
        while( _raw[vStart].isSpace()) { vStart ++; if( vStart >= 80) { vStart = -1; break; }}
        if( vStart == -1) // entire line is empty after the '='
            return;
        if( _raw[vStart] != '\'') { // it's an unquoted value
            // non-string value, find the end
            vEnd = _raw.indexOf( '/', vStart + 1); if( vEnd != -1) vEnd --; else vEnd = 79;
            //            vEnd = vStart + 1;
            //            while( ! _raw[vEnd].isSpace()) { if( vEnd >= 80) break; else vEnd ++; }
            //            vEnd --;
        } else { // it's s quoted string
            // temporarily remove all occurrences of double single-quotes and then find the next single quote
            QString tmp = _raw; for(int i=0;i<=vStart;i++){tmp[i]=' ';} tmp.replace( "''", "..");
            vEnd = tmp.indexOf( '\'', vStart + 1);
            if( vEnd == -1) // we have an unterminated string here
                throw QString( "Unterminated string in header for %1").arg(key);
        }
        // now that we know start/end, get the value
        value = _raw.mid( vStart, vEnd - vStart + 1).trimmed();

        // if this was a string value, get rid of the double single-quotes permanently, and remove the surrounding quotes too
        //if( value[0] == '\'') value = value.mid( 1, value.length()-2).replace( "''", "'");

        // is there a comment?
        comment = _raw.mid( vEnd + 1).trimmed();
        if( ! comment.isEmpty()) {
            if( comment[0] != '/')
                throw ("Syntax error in header: " + _raw.trimmed());
            else
                comment.remove(0,1);
        }
    }


protected:
    QString _raw;
};

// represents a FITS header
struct FitsHeader
{
    // do the parse of the fits file
    static FitsHeader parse( QFile & f );
    // write the header to a file
    bool write( QFile & f);
    // was the parse successful?
    bool isValid() const { return _valid; }
    // find a line with a given key
    int findLine( const QString & key ) {
        for( size_t i = 0 ; i < _lines.size() ; i ++ )
            if( _lines[i].key() == key )
                return i;
        return -1;
    }
    qint64 dataOffset() const { return _dataOffset; }
    std::vector< FitsLine > & lines() { return _lines; }

    // add a raw line to the header
    void addRaw( const QString & line );

    // sets a value in the header
    void setIntValue( const QString & key, int value, const QString & comment = QString());
    void setDoubleValue(const QString & pkey, double value, const QString & pcomment = QString());

    // general access function to key/values, does not throw exceptions but can return
    // variant with isValid() = false
    QVariant getValue( const QString & key, QVariant defaultValue = QVariant());

    // convenience functions that lookup key and convert it to requested type
    // all these throw exceptions if (a) key is not defined (b) key does not have a value
    // that can be converted to the requested type:

    // find a line with 'key' and conver it's 'value' to integer
    int intValue( const QString & key );
    int intValue( const QString & key, int defaultValue);
    QString stringValue( const QString & key );
    QString stringValue( const QString & key, const QString & defaultValue);
    double doubleValue( const QString & key );
    double doubleValue( const QString & key, double defaultValue);


protected:
    // where does the data start? This is set only in parse()! Bad design, I know.
    qint64 _dataOffset;
    // is this header valid? This is also only set in parse();
    bool _valid;
    // the lines
    std::vector< FitsLine > _lines;
    // protected constructor
    FitsHeader() { _valid = false; _dataOffset = 0; }
    // convenienty 80 spaces string
    static QString space80;
};

QString FitsHeader::space80 = "                                                                                ";

// wrapper around regular QFile::read() - it makes sure to read in requested size 's' if possible
static
bool blockRead( QFile & f, char * ptr, qint64 s)
{
    qint64 remaining = s;
    while( remaining > 0 ) {
        qint64 d = f.read( (char *) ptr, remaining);
        if( d <= 0 ) {
            std::cerr << "Error: blockRead(): could not read another block.\n";
            return false;
        }
        // update remaining & ptr
        ptr += d;
        remaining -= d;
    }
    return true;
}

// fits header parser
FitsHeader FitsHeader::parse( QFile & f)
{
    FitsHeader hdr;

    // read in header one 'line' (or card) at a time, which is 80 bytes
    // until we find one that conains END
    while( 1)
    {
        // read in another header block
        char block[80];
        if( ! blockRead( f, block, 80)) {
            std::cerr << "Error: FitsHeader::parse() could not read card.\n";
            return hdr;
        }

        // clean up the block by converting anything outside of ASCII [32..126]
        // to spaces
        for( size_t i = 0 ; i < sizeof( block) ; i ++ )
            if( block[i] < 32 || block[i] > 126)
                block[i] = ' ';

        // data offset moves
        hdr._dataOffset += sizeof( block);

        // parse the block: one line at a time (there are 36 lines of 80 chars each,
        // but they are not \0 terminated!!!)
        QString rawLine = QByteArray( (char *) block, sizeof( block) );
        // add this line to the header
        hdr._lines.push_back( rawLine);
        // if this is the 'END' line, terminate the parse
        if( rawLine.startsWith( "END     " ))
            break;
    }
    // adjust offset to be a multiple of 2880
    hdr._dataOffset = ((hdr._dataOffset -1)/ 2880 + 1) * 2880;
    // return this header
    hdr._valid = true;
    return hdr;
}



// get a value from the header as int - throwing an exception if this fails!
int FitsHeader::intValue( const QString & key)
{
    QVariant value = getValue( key);
    if( ! value.isValid())
        throw QString("Could not find key %1 in fits file.").arg(key);
    bool ok;
    int result = value.toInt( & ok);
    if( ! ok )
        throw QString("Found %1=%2 in fits file but expected an integer.").arg(key).arg(value.toString());

    // value converted, return it
    return result;
}

// get a value from the header as int - throwing an exception if this fails!
int FitsHeader::intValue( const QString & key, int defaultValue)
{
    QVariant value = getValue( key);
    if( ! value.isValid())
        return defaultValue;
    bool ok;
    int result = value.toInt( & ok);
    if( ! ok )
        throw QString("Found %1=%2 in fits file but expected an integer.").arg(key).arg(value.toString());

    // value converted, return it
    return result;
}


// get a value from the header as double - throwing an exception if this fails!
double FitsHeader::doubleValue( const QString & key)
{
    QVariant value = getValue( key);
    if( ! value.isValid())
        throw QString("Could not find key %1 in fits file.").arg(key);
    bool ok;
    double result = value.toDouble( & ok);
    if( ! ok )
        throw QString("Found %1=%2 in fits file but expected a double.").arg(key).arg(value.toString());

    // value converted, return it
    return result;
}

// get a value from the header as double - substituting default value if needed!
double FitsHeader::doubleValue( const QString & key, double defaultValue)
{
    QVariant value = getValue( key);
    if( ! value.isValid())
        return defaultValue;
    bool ok;
    double result = value.toDouble( & ok);
    if( ! ok )
        throw QString("Found %1=%2 in fits file but expected a double.").arg(key).arg(value.toString());

    // value converted, return it
    return result;
}


// get a value from the header as string - throwing an exception if this fails!
QString FitsHeader::stringValue( const QString & key)
{
    QVariant value = getValue( key);
    if( ! value.isValid())
        throw QString("Could not find key %1 in fits file.").arg(key);
    return value.toString();
}

// get a value from the header as string - throwing an exception if this fails!
QString FitsHeader::stringValue( const QString & key, const QString & defaultValue)
{
    QVariant value = getValue( key);
    if( ! value.isValid())
        return defaultValue;
    return value.toString();
}


// get a value from the header as int
QVariant FitsHeader::getValue( const QString & key, QVariant defaultValue)
{
    // find the line with this key
    int ind = findLine( key);

    // if there is no such line, report error
    if( ind < 0 )
        return defaultValue;

    // return the value as qvariant
    return QVariant( _lines[ind].value());
}

// set an integer value
void FitsHeader::setIntValue(const QString & pkey, int value, const QString & pcomment)
{
    QString key = (pkey + space80).left(8);
    QString comment = (pcomment + space80).left( 47);
    // construct a line based on the parameters
    QString rawLine = QString( "%1= %2 /  %3").arg( key, -8).arg( value, 20).arg( comment);
    rawLine = (rawLine + space80).left(80); // just in case :)
    // find a line with this key so that we can decide if we are adding a new line or
    // replacing an existing one
    int ind = findLine( pkey);
    if( ind < 0 )
        _lines.push_back( rawLine);
    else
        _lines[ind] = rawLine;
}

// set a double value
void FitsHeader::setDoubleValue(const QString & pkey, double value, const QString & pcomment)
{
    QString space80 = "                                                                                ";
    QString key = (pkey + space80).left(8);
    QString comment = (pcomment + space80).left( 47);
    // construct a line based on the parameters
    QString rawLine = QString( "%1= %2 /  %3").arg( key, -8).arg( value, 20, 'G', 10).arg( comment);
    rawLine = (rawLine + space80).left(80); // just in case :)
    // find a line with this key so that we can decide if we are adding a new line or
    // replacing an existing one
    int ind = findLine( pkey);
    if( ind < 0 )
        _lines.push_back( rawLine);
    else
        _lines[ind] = rawLine;
}

// insert a raw line into fits - no syntax checking is done, except making sure it's padded to 80 chars
void FitsHeader::addRaw(const QString & line)
{
    _lines.push_back( (line + space80).left(80));
}

// convenience function to convert fits string to a raw string (removing intial quotes & replacing all double quotes with single ones)
static QString fitsString2raw( const QString & s)
{
    if( s.length() < 2) throw "fitsString2raw - string less than 2 characters.";
    QString res = s;
    // remove the leading and ending quotes
    res[0] = res[ res.length()-1] = ' ';
    // replace all double single-quotes with a single quote
    res.replace( "''", "'");

    return res;
}

// remove leading/trailing spaces from a fits string
static QString fitsStringTrimmed( const QString & s)
{
    return QString( "'%1'").arg(fitsString2raw(s).trimmed());
}


// loads a fits file
// - opens it up
// - parses it's header & stores this
bool
FitsParser::loadFile(const QString & fname)
{
    // if we had a file open previously, get rid of it
    if( file_ != 0) delete file_;
    file_ = 0;
    if( bitPixFile_ != 0) delete bitPixFile_;
    bitPixFile_ = 0;
    cachedFrame_ = -1;
    frameCacheRaw_.clear ();

    // open the new file
    fLocation_ = fname;

    try {
        file_ = new QFile();
        file_-> setFileName( fname);

        if( ! file_-> open( QFile::ReadOnly)) return false;

        // parse the header
        FitsHeader hdr = FitsHeader::parse( * file_);
        if( ! hdr.isValid())
            throw QString( "Could not parse FITS header");

        // extract some parameters from the fits file and also validate it a bit
        if( hdr.stringValue("SIMPLE") != "T" )
            throw "Input FITS file does not have 'SIMPLE = T'";
        headerInfo_ = HeaderInfo();
        headerInfo_.bitpix = hdr.intValue( "BITPIX");
        headerInfo_.bitpixSize = bitpixToSize( headerInfo_.bitpix);
        headerInfo_.naxis = hdr.intValue( "NAXIS");
        if( headerInfo_.naxis < 2)
            throw QString( "Cannot deal with files that have NAXIS = %1").arg(headerInfo_.naxis);
        headerInfo_.naxis1 = hdr.intValue( "NAXIS1");
        headerInfo_.naxis2 = hdr.intValue( "NAXIS2");
        headerInfo_.naxis3 = hdr.intValue( "NAXIS3", 1);
        if( headerInfo_.naxis == 2) {
            headerInfo_.naxis = 3; headerInfo_.naxis3 = 1; // simulate 3d cube with 1 frame
        }

        QVariant blank = hdr.getValue( "BLANK");
        if( blank.isValid()) {
            headerInfo_.blank = hdr.intValue( "BLANK");
            headerInfo_.hasBlank = true;
            // blank is only supported for BITPIX > 0
            if( headerInfo_.bitpix < 0)
                throw QString( "Invalid use of BLANK = %1 keyword with BITPIX = %2.").arg( headerInfo_.blank).arg( headerInfo_.bitpix);
        } else {
            headerInfo_.hasBlank = false;
        }

        // calculate the number of frames in this fits file (this is basically a product
        // of NAXIS3 * NAXIS4 * ... * NAXISn
        int nFrames = 1;
        for( int i = 3 ; i <= hdr.intValue( "NAXIS") ; i ++ )
            nFrames *= hdr.intValue( QString("NAXIS%1").arg(i));
        headerInfo_.totalFrames = nFrames;

        headerInfo_.bzero = hdr.doubleValue( "BZERO", 0);
        headerInfo_.bscale = hdr.doubleValue( "BSCALE", 1);
        headerInfo_.crval1 = hdr.doubleValue( "CRVAL1", 0);
        headerInfo_.crval2 = hdr.doubleValue( "CRVAL2", 0);
        headerInfo_.crval3 = hdr.doubleValue( "CRVAL3", 0);
        headerInfo_.cdelt1 = hdr.doubleValue( "CDELT1", 1);
        headerInfo_.cdelt2 = hdr.doubleValue( "CDELT2", 1);
        headerInfo_.cdelt3 = hdr.doubleValue( "CDELT3", 1);
        headerInfo_.crpix1 = hdr.doubleValue( "CRPIX1", 0);
        headerInfo_.crpix2 = hdr.doubleValue( "CRPIX2", 0);
        headerInfo_.crpix3 = hdr.doubleValue( "CRPIX3", 0);
        headerInfo_.ctype1 = hdr.stringValue( "CTYPE1", "''");
        headerInfo_.ctype2 = hdr.stringValue( "CTYPE2", "''");
        headerInfo_.ctype3 = hdr.stringValue( "CTYPE3", "''");
        headerInfo_.cunit3 = hdr.stringValue( "CUNIT3", "''");
        headerInfo_.bunit = hdr.stringValue( "BUNIT", "''"); headerInfo_.bunit = fitsStringTrimmed( headerInfo_.bunit);
        headerInfo_.equinox = hdr.doubleValue( "EQUINOX", 2000.0);
        headerInfo_.bmaj = hdr.stringValue( "BMAJ", "undefined" );
        headerInfo_.bmin = hdr.stringValue( "BMIN", "undefined" );
        headerInfo_.bpa = hdr.stringValue( "BPA", "undefined" );

        {
            QString iquset = hdr.stringValue( "IQUSET", "File -1");
            QStringList iqa = iquset.split( ' ', QString::SkipEmptyParts);
            if( iqa.length() > 1) {
                QString s = iqa[1];
                bool ok;
                int x = s.toInt( & ok);
                if( ok) {
                    headerInfo_.iquset = x;
                }
            }
        }

        // make sure the data segment following header is big enough for the data
        qint64 inputSize = file_-> size();
        qint64 s = headerInfo_.naxis1 * headerInfo_.naxis2 * headerInfo_.naxis3 * headerInfo_.bitpixSize;
        headerInfo_.dataOffset = hdr.dataOffset();
        if( headerInfo_.dataOffset + s > inputSize)
            throw QString( "Invalid fits file size. Maybe accidentally truncated?");

        // position the input to the offset
        if( ! file_-> seek( headerInfo_.dataOffset))
            throw QString( "Could not read the data (seek failed)");

        // set the header lines
        for( size_t i = 0 ; i < hdr.lines ().size () ; i ++ )
            headerInfo_.headerLines.push_back ( hdr.lines ()[i].raw());

        // figure out whether scaling is actually required
        if( headerInfo_.hasBlank || headerInfo_.bzero != 0 || headerInfo_.bscale != 1.0)
            headerInfo_.scalingRequired = true;
        else
            headerInfo_.scalingRequired = false;
    }
    catch( const QString & err) {
        //        dbg(0) << ConsoleColors::warning () << err << ConsoleColors::resetln ();
        std::cerr << err << "\n";
        return false;
    }
    catch( ...) {
        return false;
    }

    bitPixFile_ = new M3DBitpixFile( headerInfo_, file_, headerInfo_.dataOffset, headerInfo_.naxis1, headerInfo_.naxis2, headerInfo_.totalFrames);

    //    fileInfo_.totalFrames = headerInfo_.totalFrames;
    //    if( ! fitsCache_.setFileInfo( fLocation, fileInfo_)) {
    //        dbg(1) << ConsoleColors::error() << "Failed to update cache with fits info"
    //               << ConsoleColors::resetln()
    //               << "  -> " << fitsCache_.getErrorMessage () << "\n";
    //    }
    //    frameInfos_ = fitsCache_.getAllFrames( fileInfo_);

    return true;
}

std::pair< double, double>
FitsParser::getFrameInfo(
        int frame, double preset)
{

    // put all values from the frame into an array
    QVector<double> values;
    try {
        values.reserve(headerInfo_.naxis1 * headerInfo_.naxis2); // make push_back a little faster
    } catch ( ... ) {
        // TODO: if the image is really big, we need to use offline algorithm
        // instead of nth_element()
        std::cerr << "Out of memory??? Buy more ram.\n";
        std::cerr << "Or ontact developer to write an offline selection algorithm.\n";
        exit(-1);
    }

    for( int sy = 0 ; sy < headerInfo_.naxis2 ; sy ++ ) {
        for( int sx = 0 ; sx < headerInfo_.naxis1 ; sx ++ ) {
            double v = src(sx, sy, frame);
            if( std::isfinite(v)) {
                values.push_back( v);
            }
        }
    }

    // in g++ we cannot have local functions, but we can have local structures
    // with static methods... go figure
    struct local {
        // compute auto-histogram values for the given percentage
        static void getMinMax( QVector<double> & values,
                               double percent,
                               double & outMin, double & outMax )
        {
            if( values.size() == 0) {
                // special case when we have no values... return NaNs
                outMin = outMax = std::numeric_limits<double>::quiet_NaN();
                return;
            }
            // otherwise compute the histogram values using selection algorithm
            // STL is nice enough to implement this for us
            int ind1 = round(((1-percent)/2) * values.size());
            int ind2 = values.size() - ind1 - 1;
            std::nth_element( values.begin(), values.begin() + ind1, values.end());
            outMin = values[ind1];
            std::nth_element( values.begin(), values.begin() + ind2, values.end());
            outMax = values[ind2];
        }
    };

    std::pair< double, double > res;

    //    RaiLib::HistogramInfo hist;
    local::getMinMax( values, preset, res.first, res.second);
    //    local::getMinMax( values, 0.9800, hist.min98, hist.max98);
    //    local::getMinMax( values, 0.9900, hist.min99, hist.max99);
    //    local::getMinMax( values, 0.9950, hist.min995, hist.max995);
    //    local::getMinMax( values, 0.9990, hist.min999, hist.max999);
    //    local::getMinMax( values, 0.9999, hist.min9999, hist.max9999);
    //    local::getMinMax( values, 1.0000, hist.min, hist.max);
    //    hist.valid = true;

    return res;
}

inline double
FitsParser::src(int x, int y, int z)
{
    if( bitPixFile_ == 0) throw "FitsParser::src() bitPixFile is null.";
    const char * ptr = 0;
    if( cachedFrame_ == z) {
        qint64 ind = (y * headerInfo_.naxis1 + x) * headerInfo_.bitpixSize;
        if( ind < 0 || ind >= qint64(frameCacheRaw_.size ()))
            throw "FitsParser::src() index out of bounds";
        ptr = & frameCacheRaw_[ind];
        //        return raw2double ( & frameCacheRaw_[ind]);
    } else {
        ptr = bitPixFile_->getRaw ( x, y, z);
    }
    //    return bitPixFile_->operator ()( x, y, z);
    return raw2double ( ptr);
}

// applies BLANK, BZERO and BSCALE to a raw value (passed in as a pointer)
inline double FitsParser::raw2double( const char * ptr)
{
    // convert the raw bytes to an initial double
    double original;
    switch( headerInfo_.bitpix) {
    case   8: original = double( * (quint8 *) ptr); break;
    case  16: original = double( * (qint16 *) ptr); break;
    case  32: original = double( * (qint32 *) ptr); break;
    case  64: original = double( * (qint64 *) ptr); break;
    case -32: original = double( * ( float *) ptr); break;
    case -64: original = double( * (double *) ptr); break;
    default: throw QString( "Illegal value BITPIX = %1").arg( headerInfo_.bitpix);
    }

    // if blank, bzero or bscale are not specified (& nontrivial), apply scaling,
    // otherwise return the results directly
    if( ! headerInfo_.scalingRequired) {
        return original;
    }

    // check to see if we should apply BLANK
    if( headerInfo_.hasBlank && original == headerInfo_.blank )
        return std::numeric_limits<double>::quiet_NaN();
    else
        return headerInfo_.bzero + headerInfo_.bscale * original;

}

// load an entire frame into cache
void
FitsParser::cacheFrame (int z)
{

    //    dbg(2) << "FitsParser::cacheFrame " << z << ":"
    //           << headerInfo_.naxis1 << " x " << headerInfo_.naxis2 << " pixels ...\n";
    if( cachedFrame_ == z) return;
    cachedFrame_ = z;
    frameCacheRaw_.resize( headerInfo_.naxis1 * headerInfo_.naxis2 * headerInfo_.bitpixSize);
    std::vector<char>::iterator it = frameCacheRaw_.begin ();
    const char * rawPtr = 0;
    for( int y = 0 ; y < headerInfo_.naxis2 ; y ++ ) {
        for( int x = 0 ; x < headerInfo_.naxis1 ; x ++ ) {
            rawPtr = bitPixFile_->getRaw ( x, y, z);
            it = std::copy( rawPtr, rawPtr + headerInfo_.bitpixSize, it);
        }
    }
    //    dbg(2) << "finished in " << t.elapsed () / 1000.0 << "s.\n";
}


// converts a fits frame to qimage, but only the 'rect' part, using RgbFunctor to
// do the conversion. RgbFunctor has to supply operator() which takes a double
// and return Rgb. It must handle all possible values, including NaNs.
template <typename RgbFunctor>
static
QImage fits2image( FitsParser & parser, int frame, QRect rect, RgbFunctor & rgbFunctor)
{
    const FitsParser::HeaderInfo & headerInfo = parser.getHeaderInfo ();
    if( frame < 0 || frame >= headerInfo.totalFrames )
        return QImage();

    if( rect.isNull())
        rect = QRect( 0, 0, headerInfo.naxis1, headerInfo.naxis2);

    // make an image out of the raw data
    //    QImage img( headerInfo.naxis1, headerInfo.naxis2, QImage::Format_RGB888);
    QImage img( headerInfo.naxis1, headerInfo.naxis2, QImage::Format_ARGB32);
    img.fill( QColor("black").rgb());

    // make sure this frame is cached, to make src() faster
    parser.cacheFrame ( frame);

    for( int y = rect.top() ; y <= rect.bottom() ; y ++ ) {
        QRgb * scanLine = (QRgb *) (img.scanLine (headerInfo.naxis2-y-1));
        for( int x = rect.left() ; x <= rect.right() ; x ++ ) {
            double val = parser.src(x,y,frame);
            Rgb rgb = rgbFunctor( val);
            // img.setPixel(x,headerInfo.naxis2-y-1,qRgb(rgb.r,rgb.g,rgb.b));
            //            scanLine[x - rect.left ()] = qRgb( rgb.r, rgb.g, rgb.b);
            scanLine[x] = qRgb( rgb.r, rgb.g, rgb.b);
        }
    }

    return img;
}

struct Args {
    QString input, output;
    enum { JPEG, PNG } outputType;
    int jpegQuality;
    int startFrame, endFrame;
    double clipMin, clipMax, clipPreset;
    bool ignoreOverwrite;
    ColormapFunction cmap;
    Args() {
        output = "%F-%N.%T";
        outputType = JPEG;
        jpegQuality = 90;
        startFrame = endFrame = -1;
        clipMin = clipMax = NaNd;
        clipPreset = 99.0;
        ignoreOverwrite = false;
        cmap = ColormapFunction::heat();
    }
};

static
void usage()
{
    std::cerr
            << "Usage: fits2image {option=value}*\n"
            << "\n"
               //            << "Options are:\n"
               //            << "\n"
            << "To specify input file 'x.fits', use:\n"
            << "  input=x.fits\n"
               //            << "\n"
            << "The output files are JPG by default, to change it to png:\n"
            << "  outputType=png\n"
               //            << "\n"
            << "To change JPG quality (default is 95%):\n"
            << "  jpegQuality=80\n"
               //            << "\n"
            << "Change which frames to convert:\n"
            << "  frames=frame-number for a single frame, or\n"
            << "  frames=startFrame,endFrame for a range\n"
            << "  frames=all for all frames (default)\n"
               //            << "\n"
            << "To change clipping values:\n"
            << "  clip=95% does a 95% preset calculation for each frame (default), or\n"
            << "  clip=-0.3,3.0 does a manual clip setting, same for all frames\n"
               //            << "\n"
            << "Change output filenames:\n"
            << "  output=%F-%N.%T (default) will generate output:\n"
            << "       x-00001.png or x-00001.jpg\n"
            << "  %F is the name of the input file, without the path\n"
            << "  %N is the frame number padded with '0'\n"
            << "  %n is the frame number without padding\n"
            << "  %T is the output type, 'png' or 'jpg'\n"
               //            << "\n"
            << "Change colormap:\n"
            << "  cmap=heat (default)\n"
            << "  other colormaps are gray, fire, spring, sea, sunbow, mutant\n"
            << "  aberration, rgb, velocity, helix1, helix2\n"
            << "By default files will not be overwritten. To enable overwriting:\n"
            << "  overwrite=true\n"
               ;
    exit( -1);
}

static
Args parseArgs( int argc, char ** argv)
{
    Args res;

    if( argc < 2) {
        usage();
    }

    bool argError = false;
    int errorCount = 0;
    for( int i = 1 ; i <= argc ; i ++ ) {
        if( argError) {
            std::cerr << "Error: I do not understand '" << argv[i-1] << "'\n";
            errorCount ++;
        }
        if( i == argc) break;
        argError = true; // any continue after this will be reported as error
        QString arg( argv[i]);
        int eqPos = arg.indexOf( '=');
        if( eqPos < 1) continue;
        QString par = arg.left( eqPos);
        QString val = arg.mid( eqPos + 1);
        //        std::cout << "par=" << par << " val=" << val << "\n";
        if( par == "input") {
            res.input = val;
        }
        else if( par == "output") {
            res.output = val;
        }
        else if( par == "overwrite") {
            res.ignoreOverwrite = val.toLower() == "true";
        }
        else if( par == "cmap") {
            if( val.toLower() == "heat")
                res.cmap = ColormapFunction::heat();
            else if( val.toLower() == "gray")
                res.cmap = ColormapFunction::gray();
            else if( val.toLower() == "fire")
                res.cmap = ColormapFunction::fire();
            else if( val.toLower() == "spring")
                res.cmap = ColormapFunction::spring();
            else if( val.toLower() == "sea")
                res.cmap = ColormapFunction::sea();
            else if( val.toLower() == "sunbow")
                res.cmap = ColormapFunction::sunbow();
            else if( val.toLower() == "mutant")
                res.cmap = ColormapFunction::mutant();
            else if( val.toLower() == "aberration")
                res.cmap = ColormapFunction::aberration();
            else if( val.toLower() == "rgb")
                res.cmap = ColormapFunction::rgb();
            else if( val.toLower() == "velocity")
                res.cmap = ColormapFunction::velocity();
            else if( val.toLower() == "helix1")
                res.cmap = ColormapFunction::cubeHelix( 0.5, -1.5, 1.5, 1.0);
            else if( val.toLower() == "helix2")
                res.cmap = ColormapFunction::cubeHelix( 0.5, -1.5, 1.0, 0.8);
            else
                continue;
        }
        else if( par == "jpegQuality") {
            bool ok;
            res.jpegQuality = val.toDouble( & ok);
            if( ! ok) continue;
            res.jpegQuality = clamp( res.jpegQuality, 0, 100);
        }
        else if( par == "frames") {
            if( val.toLower() == "all") {
                res.startFrame = res.endFrame = -1;
            }
            else {
                QStringList lst = val.split(",");
                bool ok;
                if( lst.size() == 1) {
                    res.startFrame = res.endFrame = lst[0].toInt( & ok);
                    if( ! ok) continue;
                }
                else if( lst.size() == 2) {
                    res.startFrame = lst[0].toInt( & ok);
                    if( ! ok) continue;
                    res.endFrame = lst[1].toInt( & ok);
                    if( ! ok) continue;

                }
                else {
                    continue;
                }
            }
        }
        else if( par == "outputType") {
            if( val.toLower() == "jpeg" || val.toLower() == "jpg") {
                res.outputType = res.JPEG;
            }
            else if( val.toLower() == "png" || val.toLower() == "png") {
                res.outputType = res.PNG;
            }
            else {
                continue;
            }
        }
        else if( par == "clip") {
            if( val.endsWith('%')) {
                res.clipMin = res.clipMax = NaNd;
                val = val.left( val.size()-1);
                bool ok;
                res.clipPreset = val.toDouble( & ok);
                if( ! ok) continue;
            } else if( val.count(',') == 1) {
                res.clipPreset = NaNd;
                QStringList lst = val.split(",");
                if( lst.size() != 2) continue;
                bool ok;
                res.clipMin = lst[0].toDouble( & ok);
                if( ! ok) continue;
                res.clipMax = lst[1].toDouble( & ok);
                if( ! ok) continue;
            } else {
                continue;
            }
        }
        else {
            continue;
        }
        argError = false;
    }

    if( errorCount > 0) {
        std::cerr << "\n";
        usage();
    }

    return res;
}

static
int cppmain(int argc, char *argv[]){
    // construct a non-gui app (initialize Qt)
    QApplication app(argc, argv, false);

    // parse arguments
    Args args = parseArgs( argc, argv);

    // try to open up the fits file
    FitsParser parser;
    if( ! parser.loadFile( args.input)) {
        std::cerr << "Could not open/parse file " << args.input << "\n";
        exit( -1);
    }

    // fix up frame range
    int nFrames = parser.getHeaderInfo().totalFrames;
    if( args.startFrame == args.endFrame && args.startFrame == -1) {
        args.startFrame = 0;
        args.endFrame = nFrames - 1;
    }
    swap_ordered( args.startFrame, args.endFrame);
    args.startFrame = clamp( args.startFrame, 0, nFrames - 1);
    args.endFrame = clamp( args.endFrame, 0, nFrames - 1);

    if( ! isnan( args.clipPreset)) {
        args.clipPreset = clamp( args.clipPreset, 0.0, 100.0);
//        std::cout << "Doing auto preset at " << args.clipPreset << "%\n";
        args.clipPreset /= 100.0;
    }
//    else {
//        std::cout << "Doing manual clip range " << args.clipMin
//                  << " .. " << args.clipMax << "\n";
//    }

    ColormapFunction cmap_ = args.cmap;
    HistogramColormapFunctor hcMap_;
    CachedRgbFunction cachedHcMap_;
    if( isnan( args.clipPreset)) {
        hcMap_ = HistogramColormapFunctor( args.clipMin, args.clipMax, cmap_);
        cachedHcMap_ = CachedRgbFunction(
                    hcMap_, args.clipMin, args.clipMax, 10000, Rgb(255,0,0));
    }

    //    ColormapFunction cmap_ = ColormapFunction::heat();
    //    HistogramColormapFunctor hcMap_
    //            = HistogramColormapFunctor( args.clipMin, args.clipMax, cmap_);
    //    CachedRgbFunction cachedHcMap_ = CachedRgbFunction(
    //                hcMap_, args.clipMin, args.clipMax, 10000, Rgb(255,0,0));

    for( int frame = args.startFrame ; frame <= args.endFrame ; frame ++ ) {
        parser.cacheFrame( frame);
        std::cout << "Doing frame " << frame << "\n";
        if( ! isnan( args.clipPreset)) {
            // computing clip values
//            std::cout << "  computing preset\n";
            std::pair<double,double> clp = parser.getFrameInfo( frame, args.clipPreset);
            args.clipMin = clp.first;
            args.clipMax = clp.second;
            std::cout << "  computed clips are " << args.clipMin << " .. "
                      << args.clipMax << "\n";

            // update functor map
            hcMap_= HistogramColormapFunctor( args.clipMin, args.clipMax, cmap_);
            cachedHcMap_ = CachedRgbFunction(
                        hcMap_, args.clipMin, args.clipMax, 10000, Rgb(255,0,0));
        }


        // convert frame to image
//        std::cout << "  converting to image\n";
        QImage img = fits2image( parser, frame, QRect(), cachedHcMap_ );

        // write out
        QFileInfo finfo( args.input);


        // QString QString::arg ( int a, int fieldWidth = 0, int base = 10, const QChar & fillChar = QLatin1Char( ' ' ) ) const

        QString output = args.output;
        output.replace( "%N", QString( "%1").arg( frame, 5, 10, QChar('0')));
        output.replace( "%n", QString::number(frame));
        output.replace( "%T", args.outputType == args.JPEG ? "jpg" : "png");
        output.replace( "%F", finfo.fileName());
        //        QString pF = finfo.fileName();
        //        QString pN = QString( "%1").arg( frame, 5, 10, QChar('0'))
        //        QString output = QString("%1-frame%2.%3")
        //                .arg( finfo.fileName())
        //                .arg( frame, 5, 10, QChar('0'))
        //                .arg( args.outputType == args.JPEG ? "jpg" : "png");
        std::cout << "  writing to file " << output << "\n";

        if( ! args.ignoreOverwrite) {
            if( QFileInfo(output).exists()) {
                std::cerr << "File " << output << " already exists. Aborting.\n";
                exit(-1);
            }
        }
        bool writeOk;
        if( args.outputType == args.JPEG) {
            writeOk = img.save( output, "JPG", args.jpegQuality);
        } else {
            writeOk = img.save( output, "PNG", -1);
        }
        if( ! writeOk) {
            std::cerr << "Could not write output. Aborting.\n";
            exit(-1);
        }

    }

    return 0;

}

};

int main(int argc, char *argv[])
{
    try {
        return fits2image::cppmain( argc, argv);
    } catch( ... ) {
        std::cerr << "Uncaught exception. Bailing out.\n";
    }
    exit(-1);
}
