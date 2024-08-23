/**************************************************************/
/* Static Utility class for doing a variety of pdfconversions */
/**************************************************************/
#ifndef __OXSX_DIST_TOOLS__
#define __OXSX_DIST_TOOLS__

class BinnedED;
class PDF;
class AxisCollection;
class Histogram;
class TH2D;
class TH1D;
class TH3D;

class DistTools
{
public:
    static Histogram ToHist(const PDF &, const AxisCollection &axes_);
    static Histogram ToHist(const TH1D &);
    static Histogram ToHist(const TH2D &);
    static Histogram ToHist(const TH3D &);
    static TH1D ToTH1D(const BinnedED &, const bool widthCorrect_ = false);
    static TH1D ToTH1D(const Histogram &, const bool widthCorrect_ = false);
    static TH2D ToTH2D(const Histogram &);
    static TH2D ToTH2D(const BinnedED &);
    static TH3D ToTH3D(const Histogram &);
    static TH3D ToTH3D(const BinnedED &);
};
#endif
