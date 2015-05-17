

// With USE_PHOTOSHOP_GAMMA pngtopdf just assumes gamma value of 2.2
// when cHRM chunk exists but gAMA does not exist.
#define USE_DEFAULT_GAMMA_22


#include <math.h>

#include <cassert>
#include <fstream>

#include "PNGImage.hh"

#include "ColorSpace.hh"

// Approximated sRGB
QPDFObjectHandle
ColorSpace::createsRGBObject (const PNGImage& src)
{
  QPDFObjectHandle colorspace;
  QPDFObjectHandle cal_param;
  bool             isRGB = false;

  assert( src.getCalibrationType() == calibration_srgb );

  isRGB = (src.getNComps() >= 3 || src.hasPalette()) ? true : false;

  // Parameters taken from PNG spec. section 4.2.2.3.
  cal_param = make_param_Cal(isRGB,
                             2.2,
                             0.3127, 0.329,
                             0.64, 0.33, 0.3, 0.6, 0.15, 0.06);
  if (cal_param.isNull())
    return QPDFObjectHandle::newNull();

  colorspace = QPDFObjectHandle::newArray();
  if (isRGB) {
    colorspace.appendItem(QPDFObjectHandle::newName("/CalRGB"));
  } else {
    colorspace.appendItem(QPDFObjectHandle::newName("/CalGray"));
  }
  colorspace.appendItem(cal_param);

  return  colorspace;
}

QPDFObjectHandle
ColorSpace::createICCBasedObject (const PNGImage& src, QPDF& qpdf)
{
  assert( src.getCalibrationType() == calibration_profile );

  int num_comp = src.hasPalette() ? 3 : src.getNComps();
  const std::string profile = std::string(src.getICCProfile().begin(),
                                          src.getICCProfile().end());
  QPDFObjectHandle  colorspace;

  if (profile.size() == 0)
    colorspace = QPDFObjectHandle::newNull();
  else {
    QPDFObjectHandle iccp = QPDFObjectHandle::newStream(&qpdf, profile);
    QPDFObjectHandle dict = iccp.getDict();
    dict.replaceKey("/N", QPDFObjectHandle::newInteger(num_comp));
    colorspace = QPDFObjectHandle::newArray();
    colorspace.appendItem(QPDFObjectHandle::newName("/ICCBased"));
    colorspace.appendItem(qpdf.makeIndirectObject(iccp));
  }

  return  colorspace;
}

// gAMA, cHRM:
//
//   If cHRM is present, we use CIE-Based color space. gAMA is also used here
//   if available.
#define INVALID_CHRM_VALUE(xw,yw,xr,yr,xg,yg,xb,yb) (\
  (xw) <= 0.0 || (yw) < 1.0e-10 || \
  (xr) < 0.0  || (yr) < 0.0 || (xg) < 0.0 || (yg) < 0.0 || \
  (xb) < 0.0  || (yb) < 0.0)

QPDFObjectHandle
ColorSpace::createCalRGBObject (const PNGImage& src)
{
  QPDFObjectHandle colorspace;
  QPDFObjectHandle cal_param;
  double  xw, yw, xr, yr, xg, yg, xb, yb;
  double  G;

  assert( src.getCalibrationType() == calibration_matrix );

  std::vector<float> v = src.getChromaticity();
  xw = v[0]; yw = v[1]; xr = v[2]; yr = v[3];
  xg = v[4]; yg = v[5]; xb = v[6]; yb = v[7];

  if (xw <= 0.0 || yw < 1.0e-10 ||
      xr < 0.0  || yr < 0.0 || xg < 0.0 || yg < 0.0 || xb < 0.0 || yb < 0.0) {
    std::cerr << "Invalid cHRM chunk parameters found." << std::endl;
    return  QPDFObjectHandle::newNull();
  }

  if (src.hasGamma()) {
    G = src.getGamma();
    if (G < 1.0e-2) {
      std::cerr << "Unusual Gamma value: 1.0 / " << G << std::endl;
      return QPDFObjectHandle::newNull();
    }
  } else {
  // Adobe PhotoShop CC assumes gAMA value of 2.2 to be used if
  // gAMA chunk does not exist?
#ifdef USE_DEFAULT_GAMMA_22
    G = 2.2;
#else
    G = 1.0;
#endif
  }

  cal_param = make_param_Cal(true, // isRGB
                             G, xw, yw, xr, yr, xg, yg, xb, yb);

  if (cal_param.isNull())
    return QPDFObjectHandle::newNull();

  colorspace = QPDFObjectHandle::newArray();
  colorspace.appendItem(QPDFObjectHandle::newName("/CalRGB")),
  colorspace.appendItem(cal_param);

  return  colorspace;
}

QPDFObjectHandle
ColorSpace::createCalGrayObject (const PNGImage& src)
{
  QPDFObjectHandle colorspace;
  QPDFObjectHandle cal_param;
  double  xw, yw, xr, yr, xg, yg, xb, yb;
  double  G;

  std::vector<float> v = src.getChromaticity();
  xw = v[0]; yw = v[1]; xr = v[2]; yr = v[3];
  xg = v[4]; yg = v[5]; xb = v[6]; yb = v[7];

  if (xw <= 0.0 || yw < 1.0e-10 ||
      xr < 0.0  || yr < 0.0 || xg < 0.0 || yg < 0.0 || xb < 0.0 || yb < 0.0) {
    std::cerr << "Invalid cHRM chunk parameters found." << std::endl;
    return  QPDFObjectHandle::newNull();
  }

  if (src.hasGamma()) {
    G = src.getGamma();
    if (G < 1.0e-2) {
      std::cerr << "Unusual Gamma value: 1.0 / " << G << std::endl;
      return  QPDFObjectHandle::newNull();
    }
  } else {
#ifdef USE_DEFAULT_GAMMA_22
    G = 2.2;
#else
    G = 1.0;
#endif
  }

  cal_param = make_param_Cal(false, // isRGB
                             G, xw, yw, xr, yr, xg, yg, xb, yb);

  if (cal_param.isNull())
    return  QPDFObjectHandle::newNull();

  colorspace = QPDFObjectHandle::newArray();
  colorspace.appendItem(QPDFObjectHandle::newName("/CalGray"));
  colorspace.appendItem(cal_param);

  return  colorspace;
}

QPDFObjectHandle
ColorSpace::make_param_Cal (bool isRGB,
                            double G, /* Gamma */
                            double xw, double yw,
                            double xr, double yr,
                            double xg, double yg,
                            double xb, double yb)
{
  QPDFObjectHandle cal_param;
  QPDFObjectHandle white_point, matrix, dev_gamma;
  double Xw, Yw, Zw; /* Yw = 1.0 */
  double Xr, Xg, Xb, Yr, Yb, Yg, Zr, Zg, Zb;

#if 1
   // Conversion found in
   //
   //  com.sixlegs.image.png - Java package to read and display PNG images
   //  Copyright (C) 1998, 1999, 2001 Chris Nokleberg
   //
   //  http://www.sixlegs.com/software/png/
   //
  {
    double zw, zr, zg, zb;
    double fr, fg, fb;
    double det;

    // WhitePoint
    zw = 1 - (xw + yw);
    zr = 1 - (xr + yr); zg = 1 - (xg + yg); zb = 1 - (xb + yb);
    Xw = xw / yw; Yw = 1.0; Zw = zw / yw;

    // Matrix
    det = xr * (yg * zb - zg * yb) -
          xg * (yr * zb - zr * yb) + xb * (yr * zg - zr * yg);
    if (fabs(det) < 1.0e-10) {
      std::cerr << "Non invertible matrix: "
                   "Maybe invalid value(s) specified in cHRM chunk." <<
                   std::endl;
      return  QPDFObjectHandle::newNull();
    }
    fr  = (Xw * (yg * zb - zg * yb) -
           xg * (zb - Zw * yb) + xb * (zg - Zw * yg)) / det;
    fg  = (xr * (zb - Zw * yb) -
           Xw * (yr * zb - zr * yb) + xb * (yr * Zw - zr)) / det;
    fb  = (xr * (yg * Zw - zg) -
           xg * (yr * Zw - zr) + Xw * (yr * zg - zr * yg)) / det;
    Xr = fr * xr; Yr = fr * yr; Zr = fr * zr;
    Xg = fg * xg; Yg = fg * yg; Zg = fg * zg;
    Xb = fb * xb; Yb = fb * yb; Zb = fb * zb;
  }
#else
  {
    double z = yw * ((xg - xb) * yr - (xr - xb) * yg + (xr - xg) * yb);
    Yr = yr * ((xg - xb) * yw - (xw - xb) * yg + (xw - xg) * yb) / z;
    Xr = Yr * xr / yr;
    Zr = Yr * ((1. - xr) / yr - 1.);
    Yg = -yg * ((xr - xb) * yw - (xw - xb) * yr + (xw - xr) * yb) / z;
    Xg = Yg * xg / yg;
    Zg = Yg * ((1. - xg) / yg - 1.);
    Yb = yb * ((xr - xg) * yw - (xw - xg) * yw + (xw - xr) * yg) / z;
    Xb = Yb * xb / yb;
    Zb = Yb * ((1. - xb) / yb - 1.);
    Xw = Xr + Xg + Xb;
    Yw = 1.0;
    Zw = Zr + Zg + Zb;
  }
#endif

  if (G < 1.0e-2) {
    std::cerr << "Unusual Gamma specified: 1.0 / " << G << std::endl;
    return  QPDFObjectHandle::newNull();
  }

  cal_param = QPDFObjectHandle::newDictionary();

  // White point is always required.
  white_point = QPDFObjectHandle::newArray();
  white_point.appendItem(QPDFObjectHandle::newReal(Xw));
  white_point.appendItem(QPDFObjectHandle::newReal(Yw));
  white_point.appendItem(QPDFObjectHandle::newReal(Zw));
  cal_param.replaceKey("/WhitePoint", white_point);

  // Matrix - default: Identity
  if (isRGB) {
    if (G != 1.0) {
      dev_gamma = QPDFObjectHandle::newArray();
      dev_gamma.appendItem(QPDFObjectHandle::newReal(G));
      dev_gamma.appendItem(QPDFObjectHandle::newReal(G));
      dev_gamma.appendItem(QPDFObjectHandle::newReal(G));
      cal_param.replaceKey("/Gamma", dev_gamma);
    }

    matrix = QPDFObjectHandle::newArray();
    matrix.appendItem(QPDFObjectHandle::newReal(Xr));
    matrix.appendItem(QPDFObjectHandle::newReal(Yr));
    matrix.appendItem(QPDFObjectHandle::newReal(Zr));
    matrix.appendItem(QPDFObjectHandle::newReal(Xg));
    matrix.appendItem(QPDFObjectHandle::newReal(Yg));
    matrix.appendItem(QPDFObjectHandle::newReal(Zg));
    matrix.appendItem(QPDFObjectHandle::newReal(Xb));
    matrix.appendItem(QPDFObjectHandle::newReal(Yb));
    matrix.appendItem(QPDFObjectHandle::newReal(Zb));
    cal_param.replaceKey("/Matrix", matrix);
  } else { // Gray
    if (G != 1.0)
      cal_param.replaceKey("/Gamma",
		                       QPDFObjectHandle::newReal(G));
  }

  return  cal_param;
}

QPDFObjectHandle
ColorSpace::createObject (const PNGImage& src, QPDF &qpdf)
{
  QPDFObjectHandle colorspace;

  switch (src.getCalibrationType()) {
  case calibration_profile:
    colorspace = createICCBasedObject(src, qpdf);
    break;
  case calibration_srgb:
    colorspace = createsRGBObject(src);
    break;
  case calibration_matrix:
    if (src.getNComps() >= 3)
      colorspace = createCalRGBObject(src);
    else
      colorspace = createCalGrayObject(src);
    break;
  case calibration_none:
  default:
    switch (src.getNComps()) {
    case 3: case 4:
      colorspace = QPDFObjectHandle::newName("/DeviceRGB");
      break;
    case 1: case 2:
      colorspace = QPDFObjectHandle::newName("/DeviceGray");
      break;
    }
  break;
  }

  return  colorspace;
}
