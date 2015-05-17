#ifndef __COLORSPACE_HH__
#define __COLORSPACE_HH__

enum eRenderingIntent {
  eRenderingIntentPerceptual = 0,
  eRenderingIntentRelativeColorimetric = 1,
  eRenderingIntentSaturation = 2,
  eRenderingIntentAbsoluteColorimetric = 3
};

#include <qpdf/QPDF.hh>
#include <qpdf/QPDFObjectHandle.hh>
#include <qpdf/QUtil.hh>

class ColorSpace : QPDFObjectHandle
{
public:
  QPDFObjectHandle createObject          (const PNGImage& src, QPDF& qpdf);

  QPDFObjectHandle createCalRGBObject    (const PNGImage& src);
  QPDFObjectHandle createCalGrayObject   (const PNGImage& src);
  QPDFObjectHandle createICCBasedObject  (const PNGImage& src, QPDF& qpdf);
  QPDFObjectHandle createsRGBObject      (const PNGImage& src);

private:
  static QPDFObjectHandle make_param_Cal (bool isRGB,
                                          double G, /* Gamma */
                                          double xw, double yw,
                                          double xr, double yr,
                                          double xg, double yg,
                                          double xb, double yb);
}

#endif // __COLORSPACE_HH__
