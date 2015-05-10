// TODO: Preserving iTXt, 2 and 4 bit gray scale image

// strlen
#include <string.h>
#include <stdio.h>
#include <math.h>

#include <iostream>
#include <string>
#include <vector>
#include <cassert>

#ifdef _WIN32
#  define FOPEN_RBIN_MODE "rb"
#  define FOPEN_WBIN_MODE "wb"
#else
#  define FOPEN_RBIN_MODE "r"
#  define FOPEN_WBIN_MODE "w"
#endif

#define PNG_NO_MNG_FEATURES
#define PNG_NO_PROGRESSIVE_READ

#include <png.h>

#include "PNGImage.hh"

static void warn(png_structp png_ptr, png_const_charp mesg)
{
  (void)png_ptr; (void)mesg;
}

int
PNGImage::check_for_png (FILE *fp)
{
  unsigned char sigbytes[4];

  rewind(fp);
  if (fread(sigbytes, 1, sizeof(sigbytes), fp) != sizeof(sigbytes) ||
      (png_sig_cmp(sigbytes, 0, sizeof(sigbytes))))
    return 0;
  else
    return 1;
}

// Creating an instance with data read from file
PNGImage::PNGImage (const std::string filename, int options)
  : Image(0, 0, 0, 0) // dummy
{
  FILE       *fp;
  png_structp png_ptr;
  png_infop   png_info_ptr;
  png_byte    color_type, bpc, nComps;
  png_uint_32 width, height, rowbytes;

  isValid = true;
  isIndexed = hasMaskColor = false;
  calibrationType = calibration_none;
  dpi_x = dpi_y =72;

  // Shouldn't pass files which is not in PNG format to libpng.
  fp = fopen(filename.c_str(), FOPEN_RBIN_MODE);
  if(!fp || !check_for_png(fp)) {
    isValid = false;
    return;
  }
  rewind(fp);

  png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, warn);
  if (png_ptr == NULL ||
      (png_info_ptr = png_create_info_struct(png_ptr)) == NULL) {
    if (png_ptr)
      png_destroy_read_struct(&png_ptr, NULL, NULL);
    isValid = false;
    return;
  }

#if PNG_LIBPNG_VER >= 10603
  // ignore possibly incorrect CMF bytes
  png_set_option(png_ptr, PNG_MAXIMUM_INFLATE_WINDOW, PNG_OPTION_ON);
#endif

  // Inititializing file IO.
  png_init_io(png_ptr, fp);

  // Read PNG info-header and get some info.
  png_read_info(png_ptr, png_info_ptr);
  color_type = png_get_color_type  (png_ptr, png_info_ptr);
  width      = png_get_image_width (png_ptr, png_info_ptr);
  height     = png_get_image_height(png_ptr, png_info_ptr);
  bpc        = png_get_bit_depth   (png_ptr, png_info_ptr);
  nComps     = png_get_channels    (png_ptr, png_info_ptr);

  png_uint_32 xppm = png_get_x_pixels_per_meter(png_ptr, png_info_ptr);
  png_uint_32 yppm = png_get_y_pixels_per_meter(png_ptr, png_info_ptr);

  dpi_x = xppm > 0 ? xppm * 0.0254 : 72;
  dpi_y = yppm > 0 ? yppm * 0.0254 : 72;

  // Ask libpng to convert from 16 bit-depth to 8 bit-depth.
  if ((options & eLoadOptionStrip16) && bpc > 16) {
    png_set_strip_16(png_ptr);
    bpc = 8;
  }
  // Ask libpng to gamma-correct.
  // It is wrong to assume screen gamma value 2.2 but...
  // We do gamma correction here only when uncalibrated color space is used.
  if ((options & eLoadOptionGammaCorrect) &&
      !png_get_valid(png_ptr, png_info_ptr, PNG_INFO_iCCP) &&
      !png_get_valid(png_ptr, png_info_ptr, PNG_INFO_sRGB) &&
      !png_get_valid(png_ptr, png_info_ptr, PNG_INFO_cHRM) &&
       png_get_valid(png_ptr, png_info_ptr, PNG_INFO_gAMA)) {
    double G = 1.0;
    png_get_gAMA (png_ptr, png_info_ptr, &G);
    png_set_gamma(png_ptr, 2.2, G);
  }
  // Does pre-compositio of pixel into background when transparency
  // isn't supported.
  if ((options & eLoadOptionPrecomposeAlpha) &&
      (color_type & PNG_COLOR_MASK_ALPHA)) {
    png_color_16 bg;
    bg.red = 255; bg.green = 255; bg.blue  = 255; bg.gray = 255; bg.index = 0;
    png_set_background(png_ptr, &bg, PNG_BACKGROUND_GAMMA_SCREEN, 0, 1.0);
  }

  // Read colorspace information.
  bool hassRGB, hasiCCP, hasgAMA, hascHRM;
  hassRGB = hasiCCP = hasgAMA = hascHRM = false;
  // Precedence in this order -- ICC profile and sRGB chunk take precedence.
  if (png_get_valid(png_ptr, png_info_ptr, PNG_INFO_iCCP)) {
    png_charpp  name    = NULL;
    int         compression_type = 0; // Don't know how to treat it...
    png_bytepp  profile = NULL;
    png_uint_32 proflen = 0;
    if (png_get_iCCP(png_ptr, png_info_ptr,
                     name, &compression_type, profile, &proflen)) {
      colorSpace.ICCP.name = std::string(*name, strlen(*name));
      colorSpace.ICCP.profile.resize(proflen);
      for (size_t i = 0; i < proflen; i++)
        colorSpace.ICCP.profile[i] = (*profile)[i];
      calibrationType = calibration_profile;
      hasiCCP = true;
    }
  } else if (png_get_valid(png_ptr, png_info_ptr, PNG_INFO_sRGB)) {
    int  intent;
    if (png_get_sRGB(png_ptr, png_info_ptr, &intent)) {
      // FIXME: use enum instead of strings.
      switch (intent) {
      case PNG_sRGB_INTENT_PERCEPTUAL:
        colorSpace.sRGB = "Perceptual";
        break;
      case PNG_sRGB_INTENT_RELATIVE:
        colorSpace.sRGB = "RelativeColorimetric";
        break;
      case PNG_sRGB_INTENT_ABSOLUTE:
        colorSpace.sRGB = "AbsoluteColorimetric";
        break;
      case PNG_sRGB_INTENT_SATURATION:
        colorSpace.sRGB  = "Saturation";
        break;
      }
      hassRGB = true;
    }
  } else if (png_get_valid(png_ptr, png_info_ptr, PNG_INFO_cHRM)) {
    double xw, yw, xr, yr, xg, yg, xb, yb;
    if (png_get_cHRM(png_ptr, png_info_ptr, &xw, &yw,
                     &xr, &yr, &xg, &yg, &xb, &yb)) {
      colorSpace.calRGB.xw = xw;
      colorSpace.calRGB.yw = yw;
      colorSpace.calRGB.xr = xr;
      colorSpace.calRGB.yr = yr;
      colorSpace.calRGB.xg = xg;
      colorSpace.calRGB.yg = yg;
      colorSpace.calRGB.xb = xb;
      colorSpace.calRGB.yb = yb;
      colorSpace.gamma = 2.2;    // maybe wrong but PhotoShop assume this?
      calibrationType = calibration_matrix;
      hascHRM = true;
    }
  }
  // Specifing gamma when iCCP or sRGB exists is invalid.
  if (!hasiCCP && !hassRGB &&
      png_get_valid(png_ptr, png_info_ptr, PNG_INFO_gAMA)) {
    double G;
    if (png_get_gAMA(png_ptr, png_info_ptr, &G)) {
      colorSpace.gamma = 1.0 / G;
      if ((options & eLoadOptionGammaCorrect) && !hascHRM) {
        // Ask libpng to gamma-correct images (only) when gAMA chunk exists
        // but not cHRM. We should use calibrated color space supported by
        // output format whenever possible instead of pre-compensating gamma.
        png_set_gamma(png_ptr, 2.2, G);
        colorSpace.hasGamma = hasgAMA = false;
        calibrationType = calibration_none;
      } else if (!hascHRM) {
        calibrationType = calibration_gamma_only;
        colorSpace.hasGamma = hasgAMA = true;
      } else {
        colorSpace.hasGamma = hasgAMA = true;
      }
    }
  }
  png_read_update_info(png_ptr, png_info_ptr);

  // Reading palette -- must come *before* tRNS
  isIndexed = false;
  if (png_get_valid(png_ptr, png_info_ptr, PNG_INFO_PLTE)) {
    png_colorp plte = NULL;
    int        num_plte = 0;
    png_get_PLTE(png_ptr, png_info_ptr, &plte, &num_plte);
    if (plte && num_plte > 0) {
      palette.resize(num_plte);
      for (int i = 0; i < num_plte; i++) {
        palette[i].v[0] = plte[i].red;
        palette[i].v[1] = plte[i].green;
        palette[i].v[2] = plte[i].blue;
        palette[i].n    = 3;
      }
      isIndexed = true;
    }
  }

  // Transparency information tRNS
  if (!png_get_valid(png_ptr, png_info_ptr, PNG_INFO_tRNS))
    hasMaskColor = false;
  else {
    png_bytep     trans = NULL;
    int           num_trans = 0;
    png_color_16 *colors = NULL;

    if (!png_get_tRNS(png_ptr, png_info_ptr, &trans, &num_trans, &colors))
      hasMaskColor = false;
    else {
      switch (color_type) {
      case PNG_COLOR_TYPE_PALETTE:
        // Opacity is associated with each color instead of pixel for
        // indexed color images...
        if (trans) {
          if (num_trans == 1 && trans[0] == 0) {
            // Actually color key mask
            maskColor = Color(0); // NOTE: index whithin color palette
            hasMaskColor = true;
          } else {
            // TODO:
            // There is a possibility that multiple color has opacity value
            // but their opacity values are either 0 or 1. In that case, we
            // can use more simpler color-key masking instead.
            if (num_trans <= (int) palette.size()) {
              for (int8_t i = 0; i < num_trans; i++) {
                // Add alpha value to the existing color palette.
                // We must read PLTE before tRNS!
                palette[i].v[3] = round(65535 * (trans[i] / 255.0)); // alpha
                palette[i].n    = 4;
              }
              for (int8_t i = 0; i < (int8_t) palette.size(); i ++) {
                palette[i].v[3] = 65535u;
                palette[i].n    = 4;
              }
            }
            hasMaskColor = false; // Not simple color key masking
          }
        }
        break;
      case PNG_COLOR_TYPE_GRAY:
        maskColor.v[0] = colors->gray;
        maskColor.n    = 1;
        hasMaskColor = true;
        break;
      case PNG_COLOR_TYPE_RGB:
       maskColor.v[0] = colors->red;
       maskColor.v[1] = colors->green;
       maskColor.v[2] = colors->blue;
       maskColor.n    = 4;
       hasMaskColor = true;
       break;
      }
    }
  }

  // Read raster image data.
  rowbytes = png_get_rowbytes(png_ptr, png_info_ptr);
  png_bytep  data_ptr = new png_byte[rowbytes * height];
  png_bytepp rows_p   = new png_bytep[height];
  for (uint32_t i = 0; i < height; i++)
    rows_p[i] = (unsigned char *)(data_ptr + rowbytes * i);
  png_read_image(png_ptr, rows_p);
  delete rows_p;

  // Reading file finished.
  png_read_end(png_ptr, NULL);

  // Cleanup.
  if (png_info_ptr)
    png_destroy_info_struct(png_ptr, &png_info_ptr);
  if (png_ptr)
    png_destroy_read_struct(&png_ptr, NULL, NULL);
  fclose(fp);

  std::string raster(reinterpret_cast<const char *>(data_ptr),
                     rowbytes * height);
  delete data_ptr;

  // We now have complete image data.
  // Actual timing of construction of Image base class is here.
  Image::reset(width, height, nComps, bpc, raster);
}

// TODO: Indexed color
int
PNGImage::save (const std::string filename) const
{
  FILE       *fp;
  png_structp png_ptr;
  png_infop   png_info_ptr;
  png_byte    color_type = PNG_COLOR_TYPE_RGB;

  fp = fopen(filename.c_str(), FOPEN_WBIN_MODE);
  if(!fp)
    return -1;

  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (png_ptr == NULL ||
      (png_info_ptr = png_create_info_struct(png_ptr)) == NULL) {
    if (png_ptr)
      png_destroy_read_struct(&png_ptr, NULL, NULL);
    return -1;
  }

#if PNG_LIBPNG_VER >= 10603
  // ignore possibly incorrect CMF bytes
  png_set_option(png_ptr, PNG_MAXIMUM_INFLATE_WINDOW, PNG_OPTION_ON);
#endif

  switch (getNComps()) {
  case 1:
    color_type = PNG_COLOR_TYPE_GRAY;
    break;
  case 2:
    color_type = PNG_COLOR_TYPE_GRAY_ALPHA;
    break;
  case 3:
    color_type = PNG_COLOR_TYPE_RGB;
    break;
  case 4:
    color_type = PNG_COLOR_TYPE_RGB_ALPHA;
    break;
  }
  // Inititializing file IO.
  png_init_io(png_ptr, fp);
  // Write header (8 bit colour depth)
  png_set_IHDR(png_ptr, png_info_ptr, getWidth(), getHeight(),
               getBPC(), // BPC of original image data
               color_type,
               PNG_INTERLACE_NONE,
               PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

  // Write colorspace related information.
  switch (calibrationType) {
  case calibration_gamma_only:
   png_set_gAMA(png_ptr, png_info_ptr, 1.0 / colorSpace.gamma);
   break;
  case calibration_matrix:
    if (colorSpace.hasGamma)
      png_set_gAMA(png_ptr, png_info_ptr, 1.0 / colorSpace.gamma);
    png_set_cHRM(png_ptr, png_info_ptr,
                 colorSpace.calRGB.xw, colorSpace.calRGB.yw,
                 colorSpace.calRGB.xr, colorSpace.calRGB.yr,
                 colorSpace.calRGB.xg, colorSpace.calRGB.yg,
                 colorSpace.calRGB.xb, colorSpace.calRGB.yb);
    break;
  case calibration_srgb:
    {
      int intent;
      if (colorSpace.sRGB == "Perceptual") {
        intent = PNG_sRGB_INTENT_PERCEPTUAL;
      } else if (colorSpace.sRGB == "RelativeColorimetric") {
        intent = PNG_sRGB_INTENT_RELATIVE;
      } else if (colorSpace.sRGB == "AbsoluteColorimetric") {
        intent = PNG_sRGB_INTENT_ABSOLUTE;
      } else if (colorSpace.sRGB == "Saturation") {
        intent = PNG_sRGB_INTENT_SATURATION;
      } else { // unknown
        intent = PNG_sRGB_INTENT_PERCEPTUAL;
      }
      png_set_sRGB(png_ptr, png_info_ptr, intent);
    }
    break;
  case calibration_profile:
    if (colorSpace.ICCP.profile.size() > 0) {
      png_set_iCCP(png_ptr, png_info_ptr,
                   colorSpace.ICCP.name.c_str(),
                   PNG_COMPRESSION_TYPE_DEFAULT,
                   colorSpace.ICCP.profile.data(),
                   colorSpace.ICCP.profile.size());
    }
    break;
  default:
    assert(calibrationType == calibration_none);
    break;
  }
  // Resolution convert ppi to pix-per-meter
  png_set_pHYs(png_ptr, png_info_ptr,
               round((dpi_x * 10000.) / 254.),
               round((dpi_y * 10000.) / 254.),
               PNG_RESOLUTION_METER);

  png_write_info(png_ptr, png_info_ptr);

  // Write raster image body.
  std::string raster = getPixelBytes();
  char *raster_data  = const_cast<char *>(raster.data());
  png_write_image(png_ptr, (reinterpret_cast<png_bytepp>(&raster_data)));

  png_write_end(png_ptr, NULL);
  if (png_info_ptr)
    png_free_data(png_ptr, png_info_ptr, PNG_FREE_ALL, -1);
  if (png_ptr)
    png_destroy_write_struct(&png_ptr, (png_infopp) NULL);
  fclose(fp);

  return  0;
}
