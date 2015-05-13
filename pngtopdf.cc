//  pngtopdf.cc
//
//  This is an adaptation of dvipdfmx PNG support code written by myself.
//
//  TODO: XMP metadata
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef _WIN32
#  define FOPEN_RBIN_MODE "rb"
#else
#  define FOPEN_RBIN_MODE "r"
#endif // _WIN32

// With USE_PHOTOSHOP_GAMMA pngtopdf just assumes gamma value of 2.2
// when cHRM chunk exists but gAMA does not exist.
#define USE_DEFAULT_GAMMA_22

#include <unistd.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <sstream>
#include <iostream>
#include <fstream>
#include <cassert>

#include <lcms2.h>

#include "Image.hh"
#include "PNGImage.hh"

// TODO: Need ColorSpace class
enum eRenderingIntent {
  eRenderingIntentPerceptual = 0,
  eRenderingIntentRelativeColorimetric = 1,
  eRenderingIntentSaturation = 2,
  eRenderingIntentAbsoluteColorimetric = 3
};

// TODO: Using struct may increase maintainace cost...
static struct
{
  struct {
    std::string  version;
    bool         autoIncrementVersion;
    bool         useObjectStream;
  } PDF;

  struct {
    std::string  color; // location of color profiles
  } resourceDirectory;

  struct {
    bool         convertToCMYK;

    std::string  defaultRGBProfilePath;
    std::string  sRGBProfilePath;
    std::string  defaultCMYKProfilePath;

    enum eRenderingIntent renderingIntent;
    bool useBlackPointCompensation;
    bool gammaCorrect; // do pre-compensation when gAMA chunk exists but
                       // cHRM does *not* exist.
  } colorManagement;

  struct {
    bool useFlatePredictorTIFF2;
  } options;

} config = {
  { "1.7", false, true },

  {
#if defined(_WIN32) || defined(_WIN64)
    "C:\\Windows\\System32\\Spool\\Drivers\\Color\\"
#else
    "/usr/share/color/icc/Adobe ICC Profiles/"
#endif
  },

  {
     false,

     "sRGB Color Space Profile.icm",
     "sRGB Color Space Profile.icm",
     "JapanColor2001Coated.icc", // Please change as appropriate for you.

     eRenderingIntentPerceptual,
     false,
     true,
  },

  {
    false  // Current version of QPDF disallows the use of predictors...
  }
};

class Margins
{
public:
  Margins( float left = 0.0, float right  = 0.0,
           float top  = 0.0, float bottom = 0.0 );
public:
  float left, right, top, bottom;
};

Margins::Margins(float top, float right, float bottom, float left) :
    left(left), right(right), top(top), bottom(bottom)
{
}

static void message_unavailable (const std::string& feature,
                                 const std::string& version,
                                 const std::string& required,
                                 const std::string& treatment)
{
  std::cerr << "PDF version \"" << version
            << "\" does not support " << feature << " (version "
            << required << " required). " << feature << " "
            << treatment << "." << std::endl;
}

#include <qpdf/QPDF.hh>
#include <qpdf/QPDFExc.hh>
#include <qpdf/QPDFObjectHandle.hh>
#include <qpdf/QPDFWriter.hh>
#include <qpdf/QUtil.hh>

// FIXME
static struct
{
  QPDFObjectHandle CMYKProfile; // Indirect object
} docResources;

static QPDFObjectHandle make_param_Cal (bool isRGB,
                                        double G, /* Gamma */
                                        double xw, double yw,
                                        double xr, double yr,
                                        double xg, double yg,
                                        double xb, double yb);
static QPDFObjectHandle create_colorspace_CalRGB  (const PNGImage& src);
static QPDFObjectHandle create_colorspace_CalGray (const PNGImage& src);
static QPDFObjectHandle create_colorspace_ICCBased(const PNGImage& src,
                                                   QPDF& qpdf);
static QPDFObjectHandle create_colorspace_sRGB    (const PNGImage& src);
static QPDFObjectHandle create_colorspace_Indexed (const PNGImage& src,
                                                   QPDF& qpdf);
static QPDFObjectHandle create_colorspace         (const PNGImage& src,
                                                   QPDF& qpdf);

static std::pair<std::string, std::string> strip_soft_mask(const PNGImage& src);
static std::string create_soft_mask(const PNGImage& src);
static bool        need_soft_mask(const PNGImage& src);

// Create an stream object with stream data source from file
static QPDFObjectHandle
newStreamFromFile (QPDF& qpdf, const std::string filename)
{
  std::ifstream ifs(filename.c_str(), std::ios::binary | std::ios::in);
  std::string   data((std::istreambuf_iterator<char>(ifs)),
                      std::istreambuf_iterator<char>());
  return QPDFObjectHandle::newStream(&qpdf, data);
}

// TODO: Cache profile
static cmsHTRANSFORM
setup_transform (const PNGImage& src)
{
  cmsHTRANSFORM hTransform = NULL;

  cmsHPROFILE hInProfile = NULL;
  switch (src.getCalibrationType()) {
  case calibration_profile:
    {
      std::vector<unsigned char> profile = src.getICCProfile();
      hInProfile = cmsOpenProfileFromMem(profile.data(), profile.size());
    }
    break;
  case calibration_srgb:
    {
      // first try loading sRGB ICC Profile from file.
      if (!config.colorManagement.sRGBProfilePath.empty()) {
        hInProfile =
            cmsOpenProfileFromFile(
              (config.resourceDirectory.color +
                config.colorManagement.sRGBProfilePath).c_str(), "r");
      }
      if (!hInProfile) {
        hInProfile = cmsCreate_sRGBProfile();
      }
    }
    break;
  case calibration_matrix:
    {
      // Photoshop assumes implicit gamma 2.2?
      double           gamma = src.hasGamma() ? src.getGamma() : 2.2;
      cmsCIExyY        white_point;
      cmsCIExyYTRIPLE  primaries;
      cmsToneCurve    *gamma_table[3];

      std::vector<float> v = src.getChromaticity();
      white_point.x     = v[0]; white_point.y     = v[1];
      primaries.Red.x   = v[2]; primaries.Red.y   = v[3];
      primaries.Green.x = v[4]; primaries.Green.y = v[5];
      primaries.Blue.x  = v[6]; primaries.Blue.y  = v[7];
      gamma_table[0] = gamma_table[1] = gamma_table[2]
          = cmsBuildGamma(NULL, gamma);
      hInProfile = cmsCreateRGBProfile(&white_point, &primaries, gamma_table);
      cmsFreeToneCurve(gamma_table[0]);
    }
    break;
  case calibration_gamma_only:
    {
      // Don't know how to treat this.
      cmsToneCurve    *gamma_table[3];
      cmsCIExyY        white_point = { .3127, .3290, 1.00 };
      cmsCIExyYTRIPLE  primaries   = {
                                       { .64, .33, 1.0 },
                                       { .30, .60, 1.0 },
                                       { .15, .06, 1.0 }
                                     };
     gamma_table[0] = gamma_table[1] = gamma_table[2]
         = cmsBuildGamma(NULL, src.getGamma());
     hInProfile = cmsCreateRGBProfile(&white_point, &primaries, gamma_table);
     cmsFreeToneCurve(gamma_table[0]);
   }
   break;
  default:
    // png_colorspace_device: Device dependent color
    if (!config.colorManagement.defaultRGBProfilePath.empty()) {
      hInProfile =
          cmsOpenProfileFromFile(
            (config.resourceDirectory.color +
              config.colorManagement.defaultRGBProfilePath).c_str(), "r");
    } else {
      hInProfile = cmsCreate_sRGBProfile();
    }
    break;
  }
  // Fallback to built-in sRGB
  if (hInProfile == NULL)
    hInProfile = cmsCreate_sRGBProfile();

  cmsHPROFILE hOutProfile =
    cmsOpenProfileFromFile(
      (config.resourceDirectory.color +
        config.colorManagement.defaultCMYKProfilePath).c_str(), "r");
  if (cmsGetColorSpace(hOutProfile) != cmsSigCmykData) {
    std::cerr << "ICC profile \""
              << config.colorManagement.defaultCMYKProfilePath
              << "\" not for CMYK." << std::endl;
  } else {
    hTransform =
      cmsCreateTransform(
          hInProfile,  src.getBPC() == 8 ? TYPE_RGB_8  : TYPE_RGB_16_SE,
          hOutProfile, src.getBPC() == 8 ? TYPE_CMYK_8 : TYPE_CMYK_16_SE,
          config.colorManagement.renderingIntent,
          config.colorManagement.useBlackPointCompensation ?
              cmsFLAGS_BLACKPOINTCOMPENSATION : 0);
  }
  cmsCloseProfile(hInProfile);
  cmsCloseProfile(hOutProfile);

  return  hTransform;
}

// Approximated sRGB
static QPDFObjectHandle
create_colorspace_sRGB (const PNGImage& src)
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

static QPDFObjectHandle
create_colorspace_ICCBased (const PNGImage& src, QPDF& qpdf)
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

static QPDFObjectHandle
create_colorspace_CalRGB (const PNGImage& src)
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

static QPDFObjectHandle
create_colorspace_CalGray (const PNGImage& src)
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

static QPDFObjectHandle
make_param_Cal (bool isRGB,
                double G, /* Gamma */
                double xw, double yw,
                double xr, double yr, double xg, double yg, double xb, double yb)
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

static QPDFObjectHandle
create_colorspace_Indexed (const PNGImage& src, QPDF& qpdf)
{
  std::vector<Color> palette = src.getPalette();
  std::string lookup(3*palette.size(), 0);
  for (size_t i = 0; i < palette.size(); i++) {
    lookup[3*i  ] = palette[i].v[0];
    lookup[3*i+1] = palette[i].v[1];
    lookup[3*i+2] = palette[i].v[2];
  }
  QPDFObjectHandle base;
  QPDFObjectHandle colorspace = QPDFObjectHandle::newArray();
  colorspace.appendItem(QPDFObjectHandle::newName("/Indexed"));
  if (config.colorManagement.convertToCMYK) {
    cmsHTRANSFORM hTransform = setup_transform(src);
    if (!hTransform)
      base = QPDFObjectHandle::newName("/DeviceRGB");
    else {
      size_t data_size = 4 * palette.size();
      char  *buffer    = new char[data_size];
      cmsDoTransform(hTransform, lookup.data(), buffer, palette.size());
      lookup = std::string(buffer, data_size);
      delete buffer;
      cmsDeleteTransform(hTransform);
      base = QPDFObjectHandle::newName("/DeviceCMYK");
    }
  } else {
    base = QPDFObjectHandle::newName("/DeviceRGB");
  }
  colorspace.appendItem(base);
  colorspace.appendItem(QPDFObjectHandle::newInteger(palette.size()-1));
  colorspace.appendItem(QPDFObjectHandle::newString(lookup));

  return colorspace;
}

static QPDFObjectHandle
create_colorspace (const PNGImage& src, QPDF &qpdf)
{
  QPDFObjectHandle colorspace;

  switch (src.getCalibrationType()) {
  case calibration_profile:
    colorspace = create_colorspace_ICCBased(src, qpdf);
    break;
  case calibration_srgb:
    colorspace = create_colorspace_sRGB(src);
    {
      std::string intent = src.getsRGBIntent();
    }
    break;
  case calibration_matrix:
    if (src.getNComps() >= 3)
      colorspace = create_colorspace_CalRGB(src);
    else
      colorspace = create_colorspace_CalGray(src);
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

static bool
need_soft_mask (const PNGImage& src)
{
  assert( src.hasPalette() );

  std::vector<Color> palette = src.getPalette();

  return (palette.size() > 0 && palette[0].n == 4) ? true : false;
}

// This is only for palette color images.
// PNG associates transparency with color instead of pixels for palette
// color.
static std::string
create_soft_mask (const PNGImage& src)
{
  assert( src.hasPalette() );

  int32_t     num_pixel = src.getWidth() * src.getHeight();
  std::string smask(num_pixel, 0xff);
  std::vector<Color> palette = src.getPalette();

  for (int32_t j = 0; j < src.getHeight(); j++) {
    for (int32_t i = 0; i < src.getWidth(); i++) {
      uint8_t n = round(255 * (src.getPixel(i, j).v[0] / 65535.)); // FIXME
      smask[src.getWidth() * j + i] =
        (n < palette.size()) ? round(255 * (palette[n].v[4] / 65535.)) : 0xffu;
    }
  }

  return  smask;
}

// Split the alpha channel from raster image data.
// Returns <image data without alpha channel, soft mask data>
static std::pair<std::string, std::string>
strip_soft_mask (const PNGImage& src)
{
  // We must be sure that image has alpha channel.
  // Bit depth must be either 8 or 16 here.
  assert( src.getNComps() == 2 || src.getNComps() == 4 );
  assert( src.getBPC() == 8 || src.getBPC() == 16 );

  size_t size = (src.getBPC() / 8) * src.getWidth() * src.getHeight();
  char  *smask_data_ptr = new char[size];
  char  *image_data_ptr = new char[size*(src.getNComps()-1)];

  switch (src.getNComps()) {
  case 4: // RGB
    if (src.getBPC() == 8) {
      for (int32_t j = 0; j < src.getHeight(); j++) {
        for (int32_t i = 0; i < src.getWidth(); i++) {
          Color   pixel = src.getPixel(i, j);
          int32_t pos   = src.getWidth() * j + i;
          image_data_ptr[3*pos  ] = round(255 * (pixel.v[0] / 65535.));
          image_data_ptr[3*pos+1] = round(255 * (pixel.v[1] / 65535.));
          image_data_ptr[3*pos+2] = round(255 * (pixel.v[2] / 65535.));;
          smask_data_ptr[  pos  ] = round(255 * (pixel.v[3] / 65535.));
        }
      }
    } else if (src.getBPC() == 16) {
      for (int32_t j = 0; j < src.getHeight(); j++) {
        for (int32_t i = 0; i < src.getWidth(); i++) {
          Color   pixel = src.getPixel(i, j);
          int32_t pos   = src.getWidth() * j + i;
          image_data_ptr[6*pos  ] = (pixel.v[0] >> 8) & 0xff;
          image_data_ptr[6*pos+1] =  pixel.v[0] & 0xff;
          image_data_ptr[6*pos+2] = (pixel.v[1] >> 8) & 0xff;
          image_data_ptr[6*pos+3] =  pixel.v[1] & 0xff;
          image_data_ptr[6*pos+4] = (pixel.v[2] >> 8) & 0xff;
          image_data_ptr[6*pos+5] =  pixel.v[2] & 0xff;
          smask_data_ptr[2*pos  ] = (pixel.v[3] >> 8) & 0xff;
          smask_data_ptr[2*pos+1] =  pixel.v[3] & 0xff;
        }
      }
    }
    break;
  case 2: // Gray
    if (src.getBPC() == 8) {
      for (int32_t j = 0; j < src.getHeight(); j++) {
        for (int32_t i = 0; i < src.getWidth(); i++) {
          Color   pixel = src.getPixel(i, j);
          int32_t pos   = src.getWidth() * j + i;
          image_data_ptr[pos] = round(255 * (pixel.v[0] / 65535.));
          smask_data_ptr[pos] = round(255 * (pixel.v[1] / 65535.));
        }
      }
    } else if (src.getBPC() == 16) {
      for (int32_t j = 0; j < src.getHeight(); j++) {
        for (int32_t i = 0; i < src.getWidth(); i++) {
          Color   pixel = src.getPixel(i, j);
          int32_t pos   = src.getWidth() * j + i;
          image_data_ptr[2*pos  ] = (pixel.v[0] >> 8) & 0xff;
          image_data_ptr[2*pos+1] =  pixel.v[0] & 0xff;
          smask_data_ptr[2*pos  ] = (pixel.v[1] >> 8) & 0xff;
          smask_data_ptr[2*pos+1] =  pixel.v[1] & 0xff;
        }
      }
    }
    break;
  }

  std::string raster =
      std::string(image_data_ptr, size * (src.getNComps() - 1));
  std::string smask = std::string(smask_data_ptr, size);

  delete smask_data_ptr;
  delete image_data_ptr;

  return  std::make_pair(raster, smask);
}

// Returns DecodeParms dictionary
// NOTICE: "raster" modified.
QPDFObjectHandle
apply_tiff2_filter (std::string& raster,
                    int32_t width, int32_t height, int8_t bpc, int8_t num_comp)
{
  if (bpc < 8 || num_comp > 4)
    return QPDFObjectHandle::newNull(); // Not supported yet

  for (int32_t j = 0; j < height; j++) {
    uint16_t prev[4] = {0, 0, 0, 0};
    for (int32_t i = 0; i < width; i++) {
      int32_t pos = (bpc / 8) * num_comp * (width * j + i);
      for (int c = 0; c < num_comp; c++) {
        switch (bpc) {
        case 8:
          {
            uint8_t cur   = raster[pos+c];
            int32_t sub   = cur - prev[c];
            prev[c]       = cur;
            raster[pos+c] = sub;
          }
          break;


        case 16:
          {
            uint16_t cur = ((uint8_t)raster[pos+2*c])*256 +
                             (uint8_t)raster[pos+2*c+1];
            uint16_t sub  = cur - prev[c];
            prev[c]       = cur;
            raster[pos+2*c  ] = (sub >> 8) & 0xff;
            raster[pos+2*c+1] = sub & 0xff;
          }
          break;
        }
      }
    }
  }
  QPDFObjectHandle parms = QPDFObjectHandle::newDictionary();
  parms.replaceKey("/BitsPerComponent",
                   QPDFObjectHandle::newInteger(bpc));
  parms.replaceKey("/Colors",
                   QPDFObjectHandle::newInteger(num_comp));
  parms.replaceKey("/Columns",
                   QPDFObjectHandle::newInteger(width));
  parms.replaceKey("/Predictor",
                    QPDFObjectHandle::newInteger(2));

  return  parms;
}

// TODO: Rename it
int
png_include_image (QPDF& qpdf, const std::string filename, const Margins margin)
{
  int32_t     page_width, page_height;
  std::string raster, raster_mask;
  bool        has_smask = false, use_cmyk = false;
  int         png_load_options = 0;

  // Check for PDF version setting and choose appropriate PNG load options.
  // PDF version < 1.3 supports no transparency.
  // PDF version < 1.4 does not support Soft-Mask (alpha channel)
  // PDF version < 1.5 does not suuport images with bit-depth 16.
  if (config.colorManagement.gammaCorrect)
    png_load_options |= PNGImage::eLoadOptionGammaCorrect;
  if ( config.PDF.version < "1.5" &&
      !config.PDF.autoIncrementVersion )
    png_load_options |= PNGImage::eLoadOptionStrip16;
  if ( config.PDF.version < "1.3" &&
      !config.PDF.autoIncrementVersion )
    png_load_options |= PNGImage::eLoadOptionPrecomposeAlpha;

  // Load PNG image. Object is invalid when opening file failed.
  PNGImage src(filename, png_load_options);
  if (!src.valid())
    return -1;

  // Pagesize settings: PNGImage honors information from pHYs chunk.
  // If there are no pHYs chunk, 72dpi will be assumed.
  page_width  = 72. * src.getWidth()  / src.getResolutionX()
                + margin.left + margin.right  + 0.5;
  page_height = 72. * src.getHeight() / src.getResolutionY()
                + margin.top  + margin.bottom + 0.5;

  // PDF wants alpha channel to be separated.
  // Be careful! NComps may not be valied after here...
  if (src.getNComps() == 2 || src.getNComps() == 4) {
    std::pair<std::string, std::string> data = strip_soft_mask(src);
    raster      = data.first;
    raster_mask = data.second;
    has_smask   = true;
  } else if (src.hasPalette()) { // Indexed color
    raster     = src.getPixelBytes(); // MARK
    has_smask  = false;
  } else {
    raster     = src.getPixelBytes();
    has_smask  = false;
  }

  // From here we do not use raster image data in PDFImage object "src".
  // We use "src" object only for obtaining color and image information.
  if (config.colorManagement.convertToCMYK && src.getNComps() == 3) {
    cmsHTRANSFORM hTransform = setup_transform(src);
    if (hTransform) {
      size_t num_pixel = src.getWidth() * src.getHeight();
      size_t data_size = (src.getBPC() / 8) * 4 * num_pixel;
      char  *buff      = new char[data_size];
      cmsDoTransform(hTransform, raster.data(), buff, num_pixel);
      raster = std::string(buff, data_size);
      delete buff;
      cmsDeleteTransform(hTransform);
    }
    use_cmyk = true;
  }

  // If BPC = 16, config.PDF.autoIncrementVersion must be true here since load
  // option strip16 is set when PDF version < 1.5 and auto increment feature
  // is not enabled. So simply increment version here.
  if (src.getBPC() == 16 && config.PDF.version < "1.5")
          config.PDF.version = "1.5";

  // Creating an Image XObject.
  QPDFObjectHandle image = QPDFObjectHandle::newStream(&qpdf);
  QPDFObjectHandle image_dict = image.getDict();
  image_dict.replaceKey("/Type",    QPDFObjectHandle::newName("/XObject"));
  image_dict.replaceKey("/Subtype", QPDFObjectHandle::newName("/Image"));
  image_dict.replaceKey("/Width",
                        QPDFObjectHandle::newInteger(src.getWidth()));
  image_dict.replaceKey("/Height",
                        QPDFObjectHandle::newInteger(src.getHeight()));
  image_dict.replaceKey("/BitsPerComponent",
                        QPDFObjectHandle::newInteger(src.getBPC()));

  // Metadata to be inserted here
#if 0
  QPDFObjectHandle metadata = QPDFObjectHandle::newStream(qpdf);
  // We must inspect iTXt chunk of PNG image to get XMP metadata.
  // Currently no interface for doing that is available.
  image_dict.replaceKey("/Metadata",
                        qpdf.mekeIndeirectObject(metadata));
#endif

  // Handle ColorSpace: ColorSpace can be DeviceRGB, DeviceGray, CalRGB,
  // CalGray, ICCBased according as PNG chunk information, or DeviceCMYK
  // if conversion to CMYK ColorSpace is enabled. The sRGB ColorSpace is
  // currently converted to approximate CalRGB ColorSpace, not ICC profile.
  // DeviceN might be supported for supporting output devices with more
  // colorants than just CMYK four colorants.
  QPDFObjectHandle colorspace;
  if (use_cmyk) {
    colorspace = QPDFObjectHandle::newName("/DeviceCMYK");
  } else if (src.hasPalette()) { // Indexed ColorSpace
    // For creation of Indexed ColorSpace, a stream object may be created
    // for base ColorSpace ICCBased. So, we need to pass qpdf object.
    colorspace = create_colorspace_Indexed(src, qpdf);
    // tRNs chunk may contain intermediate values of opacity, which can
    // not be supported via color-key masking. We need to create Soft-Mask
    // for them.
    if (need_soft_mask(src)) {
      raster_mask = create_soft_mask(src);
      has_smask = true;
    }
  } else {
    colorspace = create_colorspace(src, qpdf);
  }
  image_dict.replaceKey("/ColorSpace", colorspace);
  // I do not understand it well though...
  // Rendering intent is used for converting RGB color to CMYK.
  // We reuse it for rendering of resulting image.
  if (config.colorManagement.convertToCMYK) {
    std::string intent;
    switch (config.colorManagement.renderingIntent) {
    case eRenderingIntentPerceptual:
      intent = "/Perceptual";
      break;
    case eRenderingIntentSaturation:
      intent = "/Saturation";
      break;
    case eRenderingIntentRelativeColorimetric:
      intent = "/RelativeColorimetric";
      break;
    case eRenderingIntentAbsoluteColorimetric:
      intent = "/AbsoluteColorimetric";
      break;
    }
    if (!intent.empty()) {
      image_dict.replaceKey("/Intent",
                            QPDFObjectHandle::newName(intent));
    }
  } else if (src.getCalibrationType() == calibration_srgb) {
    std::string intent = src.getsRGBIntent();
    if (!intent.empty()) {
      image_dict.replaceKey("/Intent",
                            QPDFObjectHandle::newName(intent));
    }
  }

  // Handle transparency: Handling transparency is a bit complicated.
  // Transparency is a latecomer for PDF. Don't forget forbid them
  // for PDF/X-1 output.
  if (src.hasColorKeyMask()) {
    if (config.PDF.version < "1.3" && config.PDF.autoIncrementVersion)
      config.PDF.version = "1.3";
    if (config.PDF.version < "1.3") {
      // Shouldn't reach here -- pre-composition is requested to PNG loader
      // when PDF version setting is less than 1.3 and auto-increment is not
      // enabled.
      message_unavailable("Color-Key Mask",
                           config.PDF.version, "1.3", "ignored");
    } else {
      // Color-key masking.
      QPDFObjectHandle colorkey = QPDFObjectHandle::newArray();
      Color color = src.getMaskColor();
      // Color-key must be transformed to CMYK too.
      if (config.colorManagement.convertToCMYK && src.getNComps() == 3) {
        cmsHTRANSFORM hTransform = setup_transform(src);
        if (hTransform) {
          char inbuf[3], outbuf[4];
          inbuf[0] = color.v[0]; // TODO: implement operator=
          inbuf[1] = color.v[1];
          inbuf[2] = color.v[2];
          cmsDoTransform(hTransform, inbuf, outbuf, 1);
          color.v[0] = outbuf[0];
          color.v[1] = outbuf[1];
          color.v[2] = outbuf[2];
          cmsDeleteTransform(hTransform);
        }
      }
      // TODO: implemet an easy way to append items to an array.
      for (int i = 0; i < color.n; i++) {
        colorkey.appendItem(QPDFObjectHandle::newInteger(color.v[i]));
        colorkey.appendItem(QPDFObjectHandle::newInteger(color.v[i]));
      }
      image_dict.replaceKey("/Mask", colorkey);
    }
  } else if (has_smask && raster_mask.size() > 0) {
    // Soft-Mask required.
    if (config.PDF.version < "1.4") {
      if (config.PDF.autoIncrementVersion)
        config.PDF.version = "1.4";
      else {
        message_unavailable("Alpha transparency",
                            config.PDF.version, "1.4", "ignored");
        has_smask = false; // Sorry for this...
      }
    }
    // Create an Image XObject representing Soft-Mask.
    if (has_smask) {
      QPDFObjectHandle smask = QPDFObjectHandle::newStream(&qpdf);
      QPDFObjectHandle dict  = smask.getDict();
      dict.replaceKey("/Type",    QPDFObjectHandle::newName("/XObjcect"));
      dict.replaceKey("/Subtype", QPDFObjectHandle::newName("/Image"));
      dict.replaceKey("/Width"  ,
                      QPDFObjectHandle::newInteger(src.getWidth()));
      dict.replaceKey("/Height" ,
                      QPDFObjectHandle::newInteger(src.getHeight()));
      dict.replaceKey("/ColorSpace",
                      QPDFObjectHandle::newName("/DeviceGray"));
      dict.replaceKey("/BitsPerComponent",
                      QPDFObjectHandle::newInteger(src.getBPC()));
      // Not working: libqpdf removes DecodeParms!
      QPDFObjectHandle parms = QPDFObjectHandle::newNull();
      // Maybe apprication of predictor filter is unnecessary.
      if (config.options.useFlatePredictorTIFF2) {
        if (src.getBPC() >= 8) {
          parms = apply_tiff2_filter(raster_mask,
                                     src.getWidth(), src.getHeight(),
                                     src.getBPC(), 1);
          dict.replaceKey("/DecodeParms", parms);
        }
      }
      smask.replaceStreamData(raster_mask,
                              QPDFObjectHandle::newNull(), parms);
      image_dict.replaceKey("/SMask", smask);
    }
  }

  // TIFF predictor 2 -- horizontal differencing
  // Not working: libqpdf removes DecodeParms!
  // With strip_soft_mask() getNComps() no longer represents correct value
  // for raster image data. Don't forget that CMYK conversion also modifies
  // actual NComps.
  QPDFObjectHandle parms = QPDFObjectHandle::newNull();
  if (config.options.useFlatePredictorTIFF2) {
    int NComps =
     use_cmyk ? 4 :
                 (has_smask ? src.getNComps() - 1 : src.getNComps()); // Ugh
    if (src.getBPC() >= 8 && src.getNComps() <= 4) {
      parms = apply_tiff2_filter(raster,
                                 src.getWidth(), src.getHeight(),
                                 src.getBPC(), NComps);
      image_dict.replaceKey("/DecodeParme", parms);
    }
  }
  image.replaceStreamData(raster,
                          QPDFObjectHandle::newNull(), parms);

  // Page resources
  // Put all procset. We don't need all though.
  QPDFObjectHandle procset =
      QPDFObjectHandle::parse("[/PDF/ImageC/ImageB/ImageI]");

  // Image XObject
  QPDFObjectHandle xobject = QPDFObjectHandle::newDictionary();
  xobject.replaceKey("/Im1", image);
  QPDFObjectHandle resources = QPDFObjectHandle::newDictionary();
  resources.replaceKey("/ProcSet", procset);
  resources.replaceKey("/XObject", xobject);

  // MediaBox and ArtBox
  QPDFObjectHandle mediabox = QPDFObjectHandle::newArray();
  mediabox.appendItem(QPDFObjectHandle::newInteger(0));
  mediabox.appendItem(QPDFObjectHandle::newInteger(0));
  mediabox.appendItem(QPDFObjectHandle::newInteger(page_width));
  mediabox.appendItem(QPDFObjectHandle::newInteger(page_height));

  QPDFObjectHandle artbox = QPDFObjectHandle::newArray();
  artbox.appendItem(QPDFObjectHandle::newInteger(0));
  artbox.appendItem(QPDFObjectHandle::newInteger(0));
  artbox.appendItem(QPDFObjectHandle::newInteger(page_width));
  artbox.appendItem(QPDFObjectHandle::newInteger(page_height));

  // Create the page content stream.
  // Image should be scaled to the physical size of image, represented
  // in PostScript point -- 1inch = 72pt
  std::string content_stream =
      "q "   +
      QUtil::double_to_string(72 * src.getWidth()  / src.getResolutionX()) +
      " 0 0 " +
      QUtil::double_to_string(72 * src.getHeight() / src.getResolutionY()) +
      " " +
      QUtil::double_to_string(margin.left)   + " " +
      QUtil::double_to_string(margin.bottom) +
       " cm /Im1 Do Q\n";
  QPDFObjectHandle contents =
      QPDFObjectHandle::newStream(&qpdf, content_stream);

  // Create Page dictionary
  QPDFObjectHandle page =
      qpdf.makeIndirectObject(QPDFObjectHandle::newDictionary());
  page.replaceKey("/Type",      QPDFObjectHandle::newName("/Page"));
  page.replaceKey("/MediaBox",  mediabox);
  page.replaceKey("/ArtBox",    artbox);
  page.replaceKey("/Contents",  contents);
  // TODO: DefaultCMYK may be added to Page Resources.
  page.replaceKey("/Resources", resources);

  // Add page to PDF file. To append a page use "false" as second argument.
  qpdf.addPage(page, false);

  return  0;
}


// Main
static const int EXIT_ERROR = 2;
static const int EXIT_OK    = 0;

static const char* const usage = "bals!";

static void error_exit (const std::string& mesg)
{
  if (!mesg.empty())
    std::cerr << mesg << std::endl;
  exit(EXIT_ERROR);
}

static int
optarg_parse_margins (const char *arg, Margins& margin)
{
  const char *p = arg;
  float v[4] = {0.0, 0.0, 0.0, 0.0};
  int   count = 0;

  while (count < 4 && p[0] != 0) {
    float unit = 1.0;
    char *nextp;
    v[count] = strtof(p, &nextp);
    if (p == nextp)
      break;
    switch (nextp[0]) {
    case ' ': case '\t': case 0:
      // no unit ... in point
      while (isspace(*p))
        p++;
      unit = 1.0;
      break;
    default:
      if (strlen(nextp) < 2) {
        error_exit(usage);
      } else if (!strncmp(nextp, "pt", 2)) {
        p = nextp + 2;
      } else if (!strncmp(nextp, "in", 2)) {
        p = nextp + 2;
        unit = 72.0;
      } else if (!strncmp(nextp, "cm", 2)) {
        p = nextp + 2;
        unit = 72.0 / 2.54;
      } else if (!strncmp(nextp, "mm", 2)) {
        p = nextp + 2;
        unit = 72.0 / 25.4;
      } else {
        error_exit(usage);
      }
      while (isspace(*p))
        p++;
      break;
    }
    v[count] *= unit;
    count++;
  }
  if (p[0] != 0) {
    error_exit(usage);
  }
  // CSS style margin specification:
  //   4 values: top right bottom left
  //   3 values: top right-left bottom
  //   2 values: top-bottom right-left
  //   1 values: top-bottom-right-left
  if (count == 0) {
    error_exit(usage);
  } else if (count == 1) {
    v[1] = v[2] = v[3] = v[0];
  } else if (count == 2) {
    v[2] = v[0];
    v[3] = v[1];
  } else if (count == 3) {
    v[3] = v[1];
  }

  margin.top    = v[0];
  margin.right  = v[1];
  margin.bottom = v[2];
  margin.left   = v[3];
  return 0;
}

struct permission
{
  enum qpdf_r3_print_e  print;
  enum qpdf_r3_modify_e modify;
  bool extract;
  bool accessibility;
};

static std::vector<std::string>
split (const std::string &str, char sep)
{
  std::vector<std::string> v;
  std::stringstream ss(str);
  std::string buffer;

  while( std::getline(ss, buffer, sep) ) {
    v.push_back(buffer);
  }

  return v;
}

static int
optarg_parse_permission (const char *arg, struct permission *perm)
{
  std::vector<std::string> list = split(arg, ' ');
  std::vector<std::string>::iterator it;

  for (it = list.begin(); it != list.end(); it++) {
    if (!it->empty()) {
      std::vector<std::string> kv = split(*it, ':');
      if (kv[0].empty()) {
        error_exit(usage);
      }
      switch (kv[0][0]) {
      case 'a': // accessibility
        if (kv[1] == "allow" || kv[1] == "yes" || kv[1] == "true") {
          perm->accessibility = true;
        } else {
          perm->accessibility = false;
        }
        break;
      case 'e': // extract
        if (kv[1] == "allow" || kv[1] == "yes" || kv[1] == "true") {
          perm->extract = true;
        } else {
          perm->extract = false;
        }
        break;
      case 'm': // modify
        if (kv[1] == "all")
          perm->modify = qpdf_r3m_all;
        else if (kv[1] == "annotate" || kv[1] == "annot")
          perm->modify = qpdf_r3m_annotate;
        else if (kv[1] == "form")
          perm->modify = qpdf_r3m_form;
        else if (kv[1] == "assemble" || kv[1] == "assem")
          perm->modify = qpdf_r3m_assembly;
        else if (kv[1] == "none") {
          perm->modify = qpdf_r3m_none;
        }
        break;
      case 'p': // print
        if (kv[1] == "full")
          perm->print = qpdf_r3p_full;
        else if (kv[1] == "low")
          perm->print = qpdf_r3p_low;
        else if (kv[1] == "none") {
          perm->print = qpdf_r3p_none;
        }
        break;
      default:
        error_exit(usage);
        break;
      }
    }
  }

  return 0;
}

int main (int argc, char* argv[])
{
  int               opt, error = 0;
  extern char      *optarg;
  extern int        optind;
  Margins           margin;
  struct permission perm = {qpdf_r3p_full, qpdf_r3m_all, true, true};
  std::string       outfile, xmpfile, upasswd, opasswd;
  bool              nocompress = false, linearize = false, encrypt = false,
                    use_RC4 = false;
  int               keysize = 40;

  while ((opt = getopt(argc,
                       argv, "bBcCfFgGsSi:o:v:m:x:zZlLeEK:U:O:P:R")) != -1)
  {
    switch (opt) {
    // Some options actually are toggle switches. Use lowercase for enable
    // and uppercase for disable.

    // Color management section.
    // User Black Point compensation for conversion to CMYK
    case 'b':
      config.colorManagement.useBlackPointCompensation = true;
      break;
    case 'B':
      config.colorManagement.useBlackPointCompensation = false;
      break;
    // Enable conversion to CMYK.
    case 'c':
      config.colorManagement.convertToCMYK = true;
      break;
    case 'C':
      config.colorManagement.convertToCMYK = false;
      break;
    // Do gamma-precompensation for some situation.
    // PDF do not support the case where simply only gamma value is supplied.
    case 'g':
      config.colorManagement.gammaCorrect = true;
      break;
    case 'G':
      config.colorManagement.gammaCorrect = false;
      break;
    // Rendering intent: I don't understand it well...
    // It can be used when conversion to CMYK happens.
    case 'i':
      if (!strcmp(optarg, "relative")) {
        config.colorManagement.renderingIntent =
            eRenderingIntentRelativeColorimetric;
      } else if (!strcmp(optarg, "absolute")) {
        config.colorManagement.renderingIntent =
            eRenderingIntentAbsoluteColorimetric;
      } else if (!strcmp(optarg, "perceptual")) {
        config.colorManagement.renderingIntent = eRenderingIntentPerceptual;
      } else if (!strcmp(optarg, "saturation")) {
        config.colorManagement.renderingIntent = eRenderingIntentSaturation;
      }
      break;

    // Not working... QPDF doesn't like DecodeParms supplied.
    case 'f':
      config.options.useFlatePredictorTIFF2 = true;
      break;
    case 'F':
      config.options.useFlatePredictorTIFF2 = false;
      break;

    // Use object stream -- requires PDF version 1.5 or higher
    case 's':
      config.PDF.useObjectStream = true;
      break;
    case 'S':
      config.PDF.useObjectStream = false;
      break;

    // Output filename.
    case 'o':
      outfile = std::string(optarg);
      break;

    // PDF version -- use appended '+' letter to suggest PDF version can
    // be automatically increased necessary for supporting various feature
    // such as transparency and encryption; e.g., "1.5+".
    // WARNING: QPDF may abort() if invalid string is supplied...
    case 'v':
      if (strlen(optarg) > 1 &&
          optarg[strlen(optarg) - 1] == '+') {
        config.PDF.version = std::string(optarg, optarg + strlen(optarg) - 1);
        config.PDF.autoIncrementVersion = true;
      } else {
        config.PDF.version = std::string(optarg);
        config.PDF.autoIncrementVersion = false;
      }
      break;

    // Page setup
    // Margins -- specified in the form "1in 30pt 2cm 20mm", CSS-style margin
    // specification.
    case 'm':
      error   = optarg_parse_margins(optarg, margin);
      break;

    // XMP metadata
    case 'x':
      xmpfile = std::string(optarg);
      break;

    // Disable compression -- useful for debugging purpose.
    case 'z': nocompress = false; break;
    case 'Z': nocompress = true;  break;

    // Linearization
    case 'l': linearize  = true;  break;
    case 'L': linearize  = false; break;

    // Encryption section
    case 'e': encrypt    = true;  break;
    case 'E': encrypt    = false; break;
    case 'K':
      keysize = atoi(optarg);
      encrypt = true;
      break;
    case 'U':
      upasswd = std::string(optarg);
      encrypt = true;
      break;
    case 'O':
      opasswd = std::string(optarg);
      encrypt = true;
      break;
    case 'P':
      error   = optarg_parse_permission(optarg, &perm);
      encrypt = true;
      break;
    // 'R' option is exception. Should be removed?
    // It is not lower/upper pair to enable/disable this feature.
    case 'R':
      use_RC4 = true;
      break;

    case ':': case '?':
      error_exit(usage);
      break;
    }
  }

  // Create an instance of QPDF object and create empty page tree node.
  QPDF qpdf;
  qpdf.emptyPDF();

  // Start processing input PNG file
  // Forcing 'o' option for multiple PNG file conversion.
  if (argc - optind > 2 && outfile.empty()) {
    // require -o option
    error_exit(usage);
  } else if (outfile.empty() && argc > 1) { // last filename is for output
    argc--;
    outfile = std::string(argv[argc]);
  } else if (argc - optind == 0) {
    error_exit(usage);
  }

  // For object reuse -- DefaultCMYK profile embedded only once. (currently
  // not embedded though)
  if (config.colorManagement.convertToCMYK) {
    QPDFObjectHandle profile =
        newStreamFromFile(qpdf, config.colorManagement.defaultCMYKProfilePath);
    profile.getDict().replaceKey("/N", QPDFObjectHandle::newInteger(4));
    docResources.CMYKProfile = qpdf.makeIndirectObject(profile);
  } else {
    docResources.CMYKProfile = QPDFObjectHandle::newNull();
  }

  // Predictor filter must be disabled if not compressed.
  if (nocompress)
      config.options.useFlatePredictorTIFF2 = false;

  // Main part
  for ( ;!error && optind < argc; optind++) {
    try
    {
      error = png_include_image(qpdf, argv[optind], margin);
    }
    catch (std::exception& e)
    {
      std::cerr << e.what() << std::endl;
      error_exit("Error occured while processing file(s). No output written.");
    }
  }

  // Document Metadata
  if (!xmpfile.empty()) {
    if (config.PDF.version < "1.4") {
      message_unavailable("Metadata", config.PDF.version, "1.4", "ignored");
    } else {
      QPDFObjectHandle catalog  = qpdf.getRoot();
      QPDFObjectHandle metadata = newStreamFromFile(qpdf, xmpfile);
      if (metadata.isStream()) {
        QPDFObjectHandle dict = metadata.getDict();
        dict.replaceKey("/Type", QPDFObjectHandle::newName("/Metadata"));
        dict.replaceKey("/Subtype", QPDFObjectHandle::newName("/XML"));
        catalog.replaceKey("/Metadata", qpdf.makeIndirectObject(metadata));
      }
    }
  }

  // Write PDF output
#if 0
  // I do not yet understand Tagged PDF and Marked Content.
  // It may be required for PDF/A-1a support.
  {
    QPDFObjectHandle catalog  = qpdf.getRoot();
    QPDFObjectHandle markinfo = QPDFObjectHandle::newDictionary();
    markinfo.replaceKey("/Marked", QPDFObjectHandle::newBool(true));
    catalog.replaceKey("/MarkInfo", markinfo);
  }
#endif
  // There should be a way to attach document-wide XMP metadata.

  // Write output PDF file.
  //
  // BE CAREFUL QPDF abort if outfile could not be opened.
  // It is often the case when output file is opened by Acrobat...
  // Do not waste your time for this!!!
  QPDFWriter w(qpdf);
  w.setOutputFilename(outfile.c_str());

  // Setup for encryption
  if (encrypt) {
    if (use_RC4) {
      if (keysize == 40) {
        w.setR2EncryptionParameters(upasswd.c_str(), opasswd.c_str(),
                                    perm.print  != qpdf_r3p_none ? true : false,
                                    perm.modify != qpdf_r3m_none ? true : false,
                                    perm.extract,
                                    perm.modify <= qpdf_r3m_annotate
                                                               ? true : false);
      } else if (keysize <= 128) {
        if ( config.PDF.version < "1.4" &&
            !config.PDF.autoIncrementVersion) {
          message_unavailable("Encryption (Key size > 40)",
                               config.PDF.version, "1.4", "disabled");
          encrypt = false;
        } else {
          w.setR3EncryptionParameters(upasswd.c_str(), opasswd.c_str(),
                                      perm.accessibility, perm.extract,
                                      perm.print, perm.modify);
        }
      } else {
        error_exit("Key length > 128 unsupported for RC4 encryption.");
      }
    } else { // AESV2 or AESV3
      if (keysize <= 128) {
        // Unencrypt Metadata unsupported yet
        if ( config.PDF.version < "1.5" &&
            !config.PDF.autoIncrementVersion) {
          message_unavailable("AESV2 encryption",
                               config.PDF.version, "1.5", "disabled");
          encrypt = false;
        } else {
          w.setR4EncryptionParameters(upasswd.c_str(), opasswd.c_str(),
                                      perm.accessibility, perm.extract,
                                      perm.print, perm.modify, true, true);
        }
      } else if (keysize == 256) {
        // Unencrypt Metadata unsupported yet
        if ( config.PDF.version <= "1.7" &&
            !config.PDF.autoIncrementVersion) {
          message_unavailable("AESV3 encryption",
                              config.PDF.version, "1.7+", "disabled");
          encrypt = false;
        } else {
          w.setR6EncryptionParameters(upasswd.c_str(), opasswd.c_str(),
                                      perm.accessibility, perm.extract,
                                      perm.print, perm.modify, true);
        }
      } else {
        error_exit("Key length > 256 unsupported.");
      }
    }
  }

  // Set PDF version
  if (config.PDF.autoIncrementVersion)
    w.setMinimumPDFVersion(config.PDF.version);
  else {
    w.forcePDFVersion(config.PDF.version);
  }

  // Other nice features of QPDF -- Object Stream and Linearization.
  if (config.PDF.useObjectStream) {
    if (config.PDF.version < "1.5") {
      if (config.PDF.autoIncrementVersion)
        config.PDF.version = "1.5";
      else {
        message_unavailable("Object Stream",
                             config.PDF.version, "1.5", "disabled");
        config.PDF.useObjectStream = false;
      }
    }
    if (config.PDF.useObjectStream)
      w.setObjectStreamMode(qpdf_o_generate);
  }
  w.setStreamDataMode(nocompress ? qpdf_s_uncompress : qpdf_s_compress);
  w.setLinearization(linearize);

  std::vector<QPDFExc> warns = qpdf.getWarnings();
  // Write
  w.write();
  for (std::vector<QPDFExc>::iterator it = warns.begin();
       it != warns.end(); it++) {
    std::cerr << "From QPDF: " << it->getMessageDetail() << std::endl;
  }

  return  EXIT_OK;
}
