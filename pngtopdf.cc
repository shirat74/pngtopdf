//  pngtopdf.cc
//
//  This is an adaptation of dvipdfmx PNG support code written by myself.
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
#define USE_PHOTOSHOP_GAMMA

#include <unistd.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ROUND(v,acc) (round(((double)(v))/(acc))*(acc))

#include <sstream>
#include <iostream>
#include <cassert>


#include <qpdf/QPDF.hh>
#include <qpdf/QPDFExc.hh>
#include <qpdf/QPDFObjectHandle.hh>
#include <qpdf/QPDFWriter.hh>
#include <qpdf/QUtil.hh>

#define PNG_NO_WRITE_SUPPORTED
#define PNG_NO_MNG_FEATURES
#define PNG_NO_PROGRESSIVE_READ

#include <png.h>

#define PDF_TRANS_TYPE_NONE   0
#define PDF_TRANS_TYPE_BINARY 1
#define PDF_TRANS_TYPE_ALPHA  2

static std::string version = "1.7";

/* ColorSpace */
static QPDFObjectHandle create_cspace_Indexed(png_structp png_ptr,
                                              png_infop info_ptr, QPDF& qpdf);
/* CIE-Based: CalRGB/CalGray */
static QPDFObjectHandle create_cspace_CalRGB (png_structp png_ptr,
                                              png_infop info_ptr);
static QPDFObjectHandle create_cspace_CalGray(png_structp png_ptr,
                                              png_infop info_ptr);
static QPDFObjectHandle make_param_Cal       (png_byte color_type,
                                              double G,
                                              double xw, double yw,
                                              double xr, double yr,
                                              double xg, double yg,
                                              double xb, double yb);
// sRGB:
//   We (and PDF) do not have direct sRGB support. The sRGB color space can be
//   precisely represented by ICC profile, but we use approximate CalRGB color
//   space.
static QPDFObjectHandle create_cspace_sRGB (png_structp png_ptr,
                                            png_infop info_ptr);
static std::string get_rendering_intent (png_structp png_ptr,
                                         png_infop info_ptr);
// ICCBased
static QPDFObjectHandle create_cspace_ICCBased (png_structp png_ptr,
                                                png_infop info_ptr,
                                                QPDF& qpdf);

static QPDFObjectHandle create_colorspace (png_structp png_ptr,
                                           png_infop   info_ptr,
                                           QPDF& qpdf);
// Transparency
static int check_transparency (png_structp png_ptr, png_infop info_ptr,
                               QPDF& qpdf);

// Color-Key Mask
static QPDFObjectHandle create_ckey_mask (png_structp png_ptr,
                                          png_infop info_ptr);
// Soft Mask:
//  create_soft_mask() is for PNG_COLOR_TYPE_PALLETE.
//  Images with alpha chunnel use strip_soft_mask().
//  An object representing mask itself is returned.
static QPDFObjectHandle create_soft_mask (png_structp png_ptr,
                                          png_infop info_ptr,
                                          png_bytep image_data_ptr,
                                          png_uint_32 width, png_uint_32 height,
                                          QPDF& qpdf);
static QPDFObjectHandle strip_soft_mask (png_structp png_ptr,
                                         png_infop info_ptr,
                                         png_bytep image_data_ptr,
                                         png_uint_32p rowbytes_ptr,
                                         png_uint_32 width, png_uint_32 height,
                                         QPDF& qpdf);

/* Read image body */
static void read_image_data (png_structp png_ptr,
                             png_bytep dest_ptr,
                             png_uint_32 height, png_uint_32 rowbytes);

int
check_for_png (FILE *fp)
{
  unsigned char sigbytes[4];

  rewind(fp);
  if (fread (sigbytes, 1, sizeof(sigbytes), fp) !=
      sizeof(sigbytes) ||
      (png_sig_cmp (sigbytes, 0, sizeof(sigbytes))))
    return 0;
  else
    return 1;
}

static void png_warn (png_structp png_ptr, png_const_charp msg)
{
  (void)png_ptr; /* Make compiler happy */
  std::cerr << msg << std::endl;
}

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

int
png_include_image (QPDF& qpdf, std::string filename, Margins margin)
{
  png_structp png_ptr;
  png_infop   png_info_ptr;
  png_byte    bpc, color_type;
  png_uint_32 width, height, rowbytes;
  double      xdensity, ydensity, page_width, page_height;
  int         trans_type;
  FILE       *fp;

  fp = fopen(filename.c_str(), FOPEN_RBIN_MODE);
  if (!fp) {
    std::cerr << "Could not open file: " << filename << std::endl;
    return -1;
  }
  if (!check_for_png(fp)) {
    std::cerr << "Not a PNG file: " << filename << std::endl;
    return -1;
  }

  rewind(fp);
  png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, png_warn);
  if (png_ptr == NULL ||
      (png_info_ptr = png_create_info_struct(png_ptr)) == NULL) {
    std::cerr << "Creating Libpng read/info struct failed." << std::endl;
    if (png_ptr)
      png_destroy_read_struct(&png_ptr, NULL, NULL);
    return -1;
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

  // Ask libpng to convert down to 8-bpc.
  if (bpc > 8) {
    if (version < "1.5") {
      std::cerr << "16-bpc PNG requires PDF version 1.5." << std::endl;
      png_set_strip_16(png_ptr);
      bpc = 8;
    }
  }

  // Ask libpng to gamma-correct.
  // It is wrong to assume screen gamma value 2.2 but...
  // We do gamma correction here only when uncalibrated color space is used.
  if (!png_get_valid(png_ptr, png_info_ptr, PNG_INFO_iCCP) &&
      !png_get_valid(png_ptr, png_info_ptr, PNG_INFO_sRGB) &&
      !png_get_valid(png_ptr, png_info_ptr, PNG_INFO_cHRM) &&
       png_get_valid(png_ptr, png_info_ptr, PNG_INFO_gAMA)) {
    double G = 1.0;
    png_get_gAMA (png_ptr, png_info_ptr, &G);
    png_set_gamma(png_ptr, 2.2, G);
  }

  // png_set_background() called whithin check_transparency
  // will modify raster image data.
  trans_type = check_transparency(png_ptr, png_info_ptr, qpdf);
  // check_transparency() does not do updata_info()
  png_read_update_info(png_ptr, png_info_ptr);
  rowbytes = png_get_rowbytes(png_ptr, png_info_ptr);

  // Determine physical size.
  {
    png_uint_32 xppm = png_get_x_pixels_per_meter(png_ptr, png_info_ptr);
    png_uint_32 yppm = png_get_y_pixels_per_meter(png_ptr, png_info_ptr);

    xdensity = xppm > 0 ? 72.0 / 0.0254 / xppm : 1.0;
    ydensity = yppm > 0 ? 72.0 / 0.0254 / yppm : 1.0;
    page_width  = xdensity * width  + margin.left + margin.right ;
    page_height = ydensity * height + margin.top  + margin.bottom;
  }

  // Creating an image XObject.
  QPDFObjectHandle image = QPDFObjectHandle::newStream(&qpdf, " ");
  QPDFObjectHandle image_dict = image.getDict();
  image_dict.replaceKey("/Type",    QPDFObjectHandle::newName("/XObject"));
  image_dict.replaceKey("/Subtype", QPDFObjectHandle::newName("/Image"));
  image_dict.replaceKey("/Width",   QPDFObjectHandle::newInteger(width));
  image_dict.replaceKey("/Height",  QPDFObjectHandle::newInteger(height));
  image_dict.replaceKey("/BitsPerComponent",
                         QPDFObjectHandle::newInteger(bpc));

  // ColorSpace
  if (png_get_valid(png_ptr, png_info_ptr, PNG_INFO_sRGB)) {
    std::string intent = get_rendering_intent(png_ptr, png_info_ptr);
    if (!intent.empty()) {
      image_dict.replaceKey("/Intent",
                             QPDFObjectHandle::newName(intent));
    }
  }
  QPDFObjectHandle colorspace = create_colorspace(png_ptr, png_info_ptr, qpdf);
  if (colorspace.isNull()) {
    std::cerr << "Unknown/unsupported colorspace???" << std::endl;
    return -1;
  }
  image_dict.replaceKey("/ColorSpace", colorspace);

  // Read image body
  png_bytep stream_data_ptr = new png_byte[rowbytes * height];
  read_image_data(png_ptr, stream_data_ptr, height, rowbytes);
  // Handle alpha channel
  if (color_type == PNG_COLOR_TYPE_RGB ||
      color_type == PNG_COLOR_TYPE_RGB_ALPHA ||
      color_type == PNG_COLOR_TYPE_GRAY ||
      color_type == PNG_COLOR_TYPE_GRAY_ALPHA) {
    switch (trans_type) {
    case PDF_TRANS_TYPE_ALPHA:
      {
        QPDFObjectHandle smask =
            strip_soft_mask(png_ptr, png_info_ptr,
                            stream_data_ptr, &rowbytes, width, height, qpdf);
        if (!smask.isNull())
          image_dict.replaceKey("/SMask", qpdf.makeIndirectObject(smask));
      }
      break;
    case PDF_TRANS_TYPE_BINARY:
      {
        QPDFObjectHandle mask = create_ckey_mask(png_ptr, png_info_ptr);
        if (!mask.isNull())
          image_dict.replaceKey("/Mask", qpdf.makeIndirectObject(mask));
      }
      break;
    }
  } else if (color_type == PNG_COLOR_TYPE_PALETTE) {
    switch (trans_type) {
    case PDF_TRANS_TYPE_ALPHA:
      {
        QPDFObjectHandle smask =
            create_soft_mask(png_ptr, png_info_ptr,
                             stream_data_ptr, width, height, qpdf);
        if (!smask.isNull())
          image_dict.replaceKey("/SMask", qpdf.makeIndirectObject(smask));
      }
      break;
    case PDF_TRANS_TYPE_BINARY:
      {
        QPDFObjectHandle mask = create_ckey_mask(png_ptr, png_info_ptr);
        if (!mask.isNull())
          image_dict.replaceKey("/Mask", qpdf.makeIndirectObject(mask));
      }
      break;
    }
  }
  std::string raster =
      std::string(reinterpret_cast<const char *>(stream_data_ptr),
                  rowbytes * height);
  image.replaceStreamData(raster,
                          QPDFObjectHandle::newNull(),
                          QPDFObjectHandle::newNull());
  delete stream_data_ptr;

  // Finally read XMP Metadata
  // See, XMP Specification Part 3, Storage in Files
  // http://www.adobe.com/jp/devnet/xmp.html
  //
  // We require libpng version >= 1.6.14 since prior versions
  // of libpng had a bug that incorrectly treat the compression
  // flag of iTxt chunks.
#if PNG_LIBPNG_VER >= 10614
  if (version >= "1.4") {
    png_textp text_ptr;
    int       num_text;
    int       have_XMP = 0;

    num_text = png_get_text(png_ptr, png_info_ptr, &text_ptr, NULL);
    for (int i = 0; i < num_text; i++) {
      if (!memcmp(text_ptr[i].key, "XML:com.adobe.xmp", 17)) {
        /* XMP found */
        if (text_ptr[i].compression != PNG_ITXT_COMPRESSION_NONE ||
            text_ptr[i].itxt_length == 0)
          std::cerr << "Invalid value(s) in iTXt chunk for XMP Metadata.";
        else if (have_XMP)
          std::cerr << "Multiple XMP Metadata. Don't know how to treat it.";
        else {
          std::string metadata(text_ptr[i].text, text_ptr[i].itxt_length);
          QPDFObjectHandle XMP_stream =
              QPDFObjectHandle::newStream(&qpdf, metadata);
          QPDFObjectHandle XMP_stream_dict = XMP_stream.getDict();
          XMP_stream_dict.replaceKey("/Type",
                                      QPDFObjectHandle::newName("/Metadata"));
          XMP_stream_dict.replaceKey("/Subtype",
                                      QPDFObjectHandle::newName("/XML"));

          image_dict.replaceKey("/Metadata",
                                 qpdf.makeIndirectObject(XMP_stream));
          have_XMP = 1;
        }
      }
    }
  }
#endif // PNG_LIBPNG_VER

  png_read_end(png_ptr, NULL);

  // Cleanup
  if (png_info_ptr)
    png_destroy_info_struct(png_ptr, &png_info_ptr);
  if (png_ptr)
    png_destroy_read_struct(&png_ptr, NULL, NULL);

  // Page resources
  QPDFObjectHandle procset =
      QPDFObjectHandle::parse("[/PDF/ImageC/ImageB/ImageI]");
  QPDFObjectHandle xobject = QPDFObjectHandle::newDictionary();
  xobject.replaceKey("/Im1", image);
  QPDFObjectHandle resources = QPDFObjectHandle::newDictionary();
  resources.replaceKey("/ProcSet", procset);
  resources.replaceKey("/XObject", xobject);

  // Media box
  QPDFObjectHandle mediabox = QPDFObjectHandle::newArray();
  mediabox.appendItem(QPDFObjectHandle::newInteger(0));
  mediabox.appendItem(QPDFObjectHandle::newInteger(0));
  mediabox.appendItem(QPDFObjectHandle::newInteger(page_width));
  mediabox.appendItem(QPDFObjectHandle::newInteger(page_height));

  // Create the page content stream
  std::string content_stream =
      "q "    + QUtil::double_to_string(width * xdensity) + " " +
      "0 0 " + QUtil::double_to_string(height * ydensity) + " " +
       QUtil::double_to_string(margin.left)   + " " +
       QUtil::double_to_string(margin.bottom) + " cm /Im1 Do Q\n";
  QPDFObjectHandle contents =
      QPDFObjectHandle::newStream(&qpdf, content_stream);
  // Create the page dictionary
  QPDFObjectHandle page =
      qpdf.makeIndirectObject(QPDFObjectHandle::newDictionary());
  page.replaceKey("/Type", QPDFObjectHandle::newName("/Page"));
  page.replaceKey("/MediaBox", mediabox);
  page.replaceKey("/Contents", contents);
  page.replaceKey("/Resources", resources);

  // Add the page to the PDF file
  qpdf.addPage(page, true);

  return 0;
}

//
// The returned value trans_type is the type of transparency to be used for
// this image. Possible values are:
//
//   PDF_TRANS_TYPE_NONE    No Masking will be used/required.
//   PDF_TRANS_TYPE_BINARY  Pixels are either fully opaque/fully transparent.
//   PDF_TRANS_TYPE_ALPHA   Uses alpha channel, requies SMask.(PDF-1.4)
//
// check_transparency() must check the current setting of output PDF version
// and must choose appropriate trans_type value according to PDF version of
// current output PDF document.
//
// If the PDF version is less than 1.3, no transparency is supported for this
// version of PDF, hence PDF_TRANS_TYPE_NONE must be returned. And when the PDF
// version is equal to 1.3, possible retrun values are PDF_TRANS_TYPE_BINARY or
// PDF_TRANS_TYPE_NONE. The latter case arises when PNG file uses alpha channel
// explicitly (color type PNG_COLOR_TYPE_XXX_ALPHA), or the tRNS chunk for the
// PNG_COLOR_TYPE_PALETTE image contains intermediate values of opacity.
//
// Finally, in the case of PDF version 1.4, all kind of translucent pixels can
// be represented with Soft-Mask.

static int
check_transparency (png_structp png_ptr, png_infop info_ptr, QPDF& qpdf)
{
  int           trans_type;
  png_byte      color_type;
  png_color_16p trans_values;
  png_bytep     trans;
  int           num_trans;

  color_type = png_get_color_type(png_ptr, info_ptr);
  // First we set trans_type to appropriate value for PNG image.
  if (color_type == PNG_COLOR_TYPE_RGB_ALPHA ||
      color_type == PNG_COLOR_TYPE_GRAY_ALPHA) {
    trans_type = PDF_TRANS_TYPE_ALPHA;
  } else if (png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS) &&
	     png_get_tRNS(png_ptr, info_ptr, &trans, &num_trans, &trans_values)) {
    // Have valid tRNS chunk.
    switch (color_type) {
    case PNG_COLOR_TYPE_PALETTE:
      // Use color-key mask if possible.
      trans_type = PDF_TRANS_TYPE_BINARY;
      while (num_trans-- > 0) {
        if (trans[num_trans] != 0x00 && trans[num_trans] != 0xff) {
          // This seems not to be binary transparency
          trans_type = PDF_TRANS_TYPE_ALPHA;
          break;
        }
      }
      break;
    case PNG_COLOR_TYPE_GRAY:
    case PNG_COLOR_TYPE_RGB:
      // RGB or GRAY, single color specified by trans_values is transparent.
      trans_type = PDF_TRANS_TYPE_BINARY;
      break;
    default:
      // Else tRNS silently ignored.
      trans_type = PDF_TRANS_TYPE_NONE;
    }
  } else { // no transparency
    trans_type = PDF_TRANS_TYPE_NONE;
  }

  // Now we check PDF version.
  // We can convert alpha cahnnels to explicit mask via user supplied alpha-
  // threshold value. But I will not do that.
  if (( version < "1.3" && trans_type != PDF_TRANS_TYPE_NONE   ) ||
      ( version < "1.4" && trans_type == PDF_TRANS_TYPE_ALPHA )) {
    // No transparency supported but PNG uses transparency, or Soft-Mask
    // required but no support for it is available in this version of PDF.
    // We must do pre-composition of image with the background image here. But,
    // we cannot do that in general since pngtopdf is not a rasterizer. What we
    // can do here is to composite image with a rectangle filled with the
    // background color. However, images are stored as an Image XObject which
    // can be referenced anywhere in the PDF document content. Hence, we cannot
    // know the correct background color at this time. So we will choose white
    // as background color, which is most probable color in our cases.
    // We ignore bKGD chunk.
    png_color_16 bg;
    bg.red = 255; bg.green = 255; bg.blue  = 255; bg.gray = 255; bg.index = 0;
    png_set_background(png_ptr, &bg, PNG_BACKGROUND_GAMMA_SCREEN, 0, 1.0);
    std::cerr <<
      "Transparency will be ignored." <<
      std::endl;
    if (version < "1.3")
      std::cerr <<
        "Please use -V 3 option to enable binary transparency support." <<
        std::endl;
    if (version < "1.4")
      std::cerr <<
        "Please use -V 4 option to enable full alpha channel support." <<
        std::endl;
    trans_type = PDF_TRANS_TYPE_NONE;
  }

  return trans_type;
}

 // sRGB:
 //
 //  If sRGB chunk is present, cHRM and gAMA chunk must be ignored.
 //
static std::string
get_rendering_intent (png_structp png_ptr, png_infop info_ptr)
{
  int         srgb_intent;
  std::string intent;

  if (png_get_valid(png_ptr, info_ptr, PNG_INFO_sRGB) &&
      png_get_sRGB (png_ptr, info_ptr, &srgb_intent)) {
    switch (srgb_intent) {
    case PNG_sRGB_INTENT_SATURATION:
      intent = "/Saturation";
      break;
    case PNG_sRGB_INTENT_PERCEPTUAL:
      intent = "/Perceptual";
      break;
    case PNG_sRGB_INTENT_ABSOLUTE:
      intent = "/AbsoluteColorimetric";
      break;
    case PNG_sRGB_INTENT_RELATIVE:
      intent = "/RelativeColorimetric";
      break;
    }
  }

  return intent;
}

// Approximated sRGB
static QPDFObjectHandle
create_cspace_sRGB (png_structp png_ptr, png_infop info_ptr)
{
  QPDFObjectHandle colorspace;
  QPDFObjectHandle cal_param;
  png_byte         color_type;

  color_type = png_get_color_type(png_ptr, info_ptr);

  // Parameters taken from PNG spec. section 4.2.2.3.
  cal_param = make_param_Cal(color_type,
                             2.2,
                             0.3127, 0.329,
                             0.64, 0.33, 0.3, 0.6, 0.15, 0.06);
  if (cal_param.isNull())
    return QPDFObjectHandle::newNull();

  colorspace = QPDFObjectHandle::newArray();

  switch (color_type) {
  case PNG_COLOR_TYPE_RGB:
  case PNG_COLOR_TYPE_RGB_ALPHA:
  case PNG_COLOR_TYPE_PALETTE:
    colorspace.appendItem(QPDFObjectHandle::newName("/CalRGB"));
    break;
  case PNG_COLOR_TYPE_GRAY:
  case PNG_COLOR_TYPE_GRAY_ALPHA:
    colorspace.appendItem(QPDFObjectHandle::newName("/CalGray"));
    break;
  }
  colorspace.appendItem(cal_param);

  return  colorspace;
}

static QPDFObjectHandle
create_cspace_ICCBased (png_structp png_ptr, png_infop info_ptr, QPDF& qpdf)
{
  png_charp  name;
  int        color_type, num_comps;
  int        compression_type;  // Manual page for libpng does not
				                        // clarify whether profile data is
                                // inflated by libpng.
#if PNG_LIBPNG_VER_MINOR < 5
  png_charp   profile;
#else
  png_bytep   profile;
#endif
  png_uint_32 proflen;

  if (!png_get_valid(png_ptr, info_ptr, PNG_INFO_iCCP) ||
      !png_get_iCCP(png_ptr, info_ptr,
                    &name, &compression_type, &profile, &proflen))
    return QPDFObjectHandle::newNull();

  color_type = png_get_color_type(png_ptr, info_ptr);
  num_comps  = (color_type == PNG_COLOR_TYPE_GRAY ||
                color_type == PNG_COLOR_TYPE_GRAY_ALPHA) ? 1 : 3; // No CMYK

  QPDFObjectHandle colorspace;
  if (proflen == 0)
    colorspace = QPDFObjectHandle::newNull();
  else {
    std::string iccp_data =
        std::string(reinterpret_cast<const char *>(profile), proflen);
    QPDFObjectHandle iccp = QPDFObjectHandle::newStream(&qpdf, iccp_data);
    QPDFObjectHandle dict = iccp.getDict();
    dict.replaceKey("/N", QPDFObjectHandle::newInteger(num_comps));
    colorspace = QPDFObjectHandle::newArray();
    colorspace.appendItem(QPDFObjectHandle::newName("/ICCBased"));
    colorspace.appendItem(qpdf.makeIndirectObject(iccp));
  }

  return colorspace;
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
create_cspace_CalRGB (png_structp png_ptr, png_infop info_ptr)
{
  QPDFObjectHandle colorspace;
  QPDFObjectHandle cal_param;
  double  xw, yw, xr, yr, xg, yg, xb, yb;
  double  G;

  if (!png_get_valid(png_ptr, info_ptr, PNG_INFO_cHRM) ||
      !png_get_cHRM(png_ptr, info_ptr, &xw, &yw, &xr, &yr, &xg, &yg, &xb, &yb))
    return  QPDFObjectHandle::newNull();

  if (xw <= 0.0 || yw < 1.0e-10 ||
      xr < 0.0  || yr < 0.0 || xg < 0.0 || yg < 0.0 || xb < 0.0 || yb < 0.0) {
    std::cerr << "Invalid cHRM chunk parameters found." << std::endl;
    return  QPDFObjectHandle::newNull();
  }

  if (png_get_valid(png_ptr, info_ptr, PNG_INFO_gAMA) &&
      png_get_gAMA (png_ptr, info_ptr, &G)) {
    if (G < 1.0e-2) {
      std::cerr << "Unusual Gamma value: 1.0 / " << G << std::endl;
      return QPDFObjectHandle::newNull();
    }
    G = 1.0 / G; /* Gamma is inverted. */
  } else {
  // Adobe PhotoShop CC assumes gAMA value of 2.2 to be used if
  // gAMA chunk does not exist?
#ifdef USE_PHOTOSHOP_GAMMA
    G = 2.2;
#else
    G = 1.0;
#endif
  }

  cal_param = make_param_Cal(PNG_COLOR_TYPE_RGB,
                             G, xw, yw, xr, yr, xg, yg, xb, yb);

  if (cal_param.isNull())
    return QPDFObjectHandle::newNull();

  colorspace = QPDFObjectHandle::newArray();
  colorspace.appendItem(QPDFObjectHandle::newName("/CalRGB")),
  colorspace.appendItem(cal_param);

  return  colorspace;
}

static QPDFObjectHandle
create_cspace_CalGray (png_structp png_ptr, png_infop info_ptr)
{
  QPDFObjectHandle colorspace;
  QPDFObjectHandle cal_param;
  double  xw, yw, xr, yr, xg, yg, xb, yb;
  double  G;

  if (!png_get_valid(png_ptr, info_ptr, PNG_INFO_cHRM) ||
      !png_get_cHRM(png_ptr, info_ptr, &xw, &yw, &xr, &yr, &xg, &yg, &xb, &yb))
    return  QPDFObjectHandle::newNull();

  if (xw <= 0.0 || yw < 1.0e-10 ||
      xr < 0.0  || yr < 0.0 || xg < 0.0 || yg < 0.0 || xb < 0.0 || yb < 0.0) {
    std::cerr << "Invalid cHRM chunk parameters found." << std::endl;
    return  QPDFObjectHandle::newNull();
  }

  if (png_get_valid(png_ptr, info_ptr, PNG_INFO_gAMA) &&
      png_get_gAMA (png_ptr, info_ptr, &G)) {
    if (G < 1.0e-2) {
      std::cerr << "Unusual Gamma value: 1.0 / " << G << std::endl;
      return  QPDFObjectHandle::newNull();
    }
    G = 1.0 / G; // Gamma is inverted.
  } else {
#ifdef USE_PHOTOSHOP_GAMMA
    G = 2.2;
#else
    G = 1.0;
#endif
  }

  cal_param = make_param_Cal(PNG_COLOR_TYPE_GRAY,
                             G, xw, yw, xr, yr, xg, yg, xb, yb);

  if (cal_param.isNull())
    return  QPDFObjectHandle::newNull();

  colorspace = QPDFObjectHandle::newArray();
  colorspace.appendItem(QPDFObjectHandle::newName("/CalGray"));
  colorspace.appendItem(cal_param);

  return  colorspace;
}

static QPDFObjectHandle
make_param_Cal (png_byte color_type,
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
  white_point.appendItem(QPDFObjectHandle::newReal(ROUND(Xw, 0.00001)));
  white_point.appendItem(QPDFObjectHandle::newReal(ROUND(Yw, 0.00001)));
  white_point.appendItem(QPDFObjectHandle::newReal(ROUND(Zw, 0.00001)));
  cal_param.replaceKey("/WhitePoint", white_point);

  // Matrix - default: Identity
  if (color_type & PNG_COLOR_MASK_COLOR) {
    if (G != 1.0) {
      dev_gamma = QPDFObjectHandle::newArray();
      dev_gamma.appendItem(QPDFObjectHandle::newReal(ROUND(G, 0.00001)));
      dev_gamma.appendItem(QPDFObjectHandle::newReal(ROUND(G, 0.00001)));
      dev_gamma.appendItem(QPDFObjectHandle::newReal(ROUND(G, 0.00001)));
      cal_param.replaceKey("/Gamma", dev_gamma);
    }

    matrix = QPDFObjectHandle::newArray();
    matrix.appendItem(QPDFObjectHandle::newReal(ROUND(Xr, 0.00001)));
    matrix.appendItem(QPDFObjectHandle::newReal(ROUND(Yr, 0.00001)));
    matrix.appendItem(QPDFObjectHandle::newReal(ROUND(Zr, 0.00001)));
    matrix.appendItem(QPDFObjectHandle::newReal(ROUND(Xg, 0.00001)));
    matrix.appendItem(QPDFObjectHandle::newReal(ROUND(Yg, 0.00001)));
    matrix.appendItem(QPDFObjectHandle::newReal(ROUND(Zg, 0.00001)));
    matrix.appendItem(QPDFObjectHandle::newReal(ROUND(Xb, 0.00001)));
    matrix.appendItem(QPDFObjectHandle::newReal(ROUND(Yb, 0.00001)));
    matrix.appendItem(QPDFObjectHandle::newReal(ROUND(Zb, 0.00001)));
    cal_param.replaceKey("/Matrix", matrix);
  } else { // Gray
    if (G != 1.0)
      cal_param.replaceKey("/Gamma",
		                       QPDFObjectHandle::newReal(ROUND(G, 0.00001)));
  }

  return  cal_param;
}

// Set up Indexed ColorSpace for color-type PALETTE:
//
// PNG allows only RGB color for base color space. If gAMA and/or cHRM
// chunk is available, we can use CalRGB color space instead of DeviceRGB
//  for base color space.
//
static QPDFObjectHandle
create_cspace_Indexed (png_structp png_ptr, png_infop info_ptr, QPDF& qpdf)
{
  png_byte  *data_ptr;
  png_colorp plte;
  int        num_plte;

  if (!png_get_valid(png_ptr, info_ptr, PNG_INFO_PLTE) ||
      !png_get_PLTE(png_ptr, info_ptr, &plte, &num_plte)) {
    std::cerr << "PNG does not have valid PLTE chunk." << std::endl;
    return  QPDFObjectHandle::newNull();
  }

  /* Order is important. */
  QPDFObjectHandle colorspace = QPDFObjectHandle::newArray();
  colorspace.appendItem(QPDFObjectHandle::newName("/Indexed"));
  // Base ColorSpace
  QPDFObjectHandle base;
  if (png_get_valid(png_ptr, info_ptr, PNG_INFO_iCCP))
    base = create_cspace_ICCBased(png_ptr, info_ptr, qpdf);
  else if (png_get_valid(png_ptr, info_ptr, PNG_INFO_sRGB))
    base = create_cspace_sRGB(png_ptr, info_ptr);
  else if (png_get_valid(png_ptr, info_ptr, PNG_INFO_cHRM))
    base = create_cspace_CalRGB(png_ptr, info_ptr);
  else {
    base = QPDFObjectHandle::newName("/DeviceRGB");
  }
  colorspace.appendItem(base);
  colorspace.appendItem(QPDFObjectHandle::newInteger(num_plte - 1));
  data_ptr = new png_byte[num_plte * 3];
  for (int i = 0; i < num_plte; i++) {
    data_ptr[3*i]   = plte[i].red;
    data_ptr[3*i+1] = plte[i].green;
    data_ptr[3*i+2] = plte[i].blue;
  }
  QPDFObjectHandle lookup =
      QPDFObjectHandle::newString(
        std::string(reinterpret_cast<const char *>(data_ptr), num_plte * 3));
  colorspace.appendItem(lookup);
  delete data_ptr;

  return colorspace;
}

//
// Colorkey Mask: array
//
//  [component_0_min component_0_max ... component_n_min component_n_max]
//

static QPDFObjectHandle
create_ckey_mask (png_structp png_ptr, png_infop info_ptr)
{
  png_byte  color_type;
  png_bytep trans;
  int       num_trans;
  png_color_16p colors;

  if (!png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS) ||
      !png_get_tRNS(png_ptr, info_ptr, &trans, &num_trans, &colors)) {
    std::cerr << "PNG does not have valid tRNS chunk!" << std::endl;
    return  QPDFObjectHandle::newNull();
  }

  QPDFObjectHandle colorkeys = QPDFObjectHandle::newArray();
  color_type = png_get_color_type(png_ptr, info_ptr);
  switch (color_type) {
  case PNG_COLOR_TYPE_PALETTE:
    for (int i = 0; i < num_trans; i++) {
      if (trans[i] == 0x00) {
        colorkeys.appendItem(QPDFObjectHandle::newReal(i));
        colorkeys.appendItem(QPDFObjectHandle::newReal(i));
      } else {
        assert(trans[i] == 0xff);
      }
    }
    break;
  case PNG_COLOR_TYPE_RGB:
    colorkeys.appendItem(QPDFObjectHandle::newReal(colors->red));
    colorkeys.appendItem(QPDFObjectHandle::newReal(colors->red));
    colorkeys.appendItem(QPDFObjectHandle::newReal(colors->green));
    colorkeys.appendItem(QPDFObjectHandle::newReal(colors->green));
    colorkeys.appendItem(QPDFObjectHandle::newReal(colors->blue));
    colorkeys.appendItem(QPDFObjectHandle::newReal(colors->blue));
    break;
  case PNG_COLOR_TYPE_GRAY:
    colorkeys.appendItem(QPDFObjectHandle::newReal(colors->gray));
    colorkeys.appendItem(QPDFObjectHandle::newReal(colors->gray));
    break;
  default:
    assert(1);
    break;
  }

  return  colorkeys;
}

static QPDFObjectHandle
create_colorspace (png_structp png_ptr, png_infop info_ptr, QPDF& qpdf)
{
  QPDFObjectHandle colorspace;
  png_byte         color_type;

  color_type = png_get_color_type(png_ptr, info_ptr);

  switch (color_type) {
  case PNG_COLOR_TYPE_PALETTE:
    colorspace = create_cspace_Indexed(png_ptr, info_ptr, qpdf);
    break;

  case PNG_COLOR_TYPE_RGB:
  case PNG_COLOR_TYPE_RGB_ALPHA:
    if (png_get_valid(png_ptr, info_ptr, PNG_INFO_iCCP)) {
      colorspace = create_cspace_ICCBased(png_ptr, info_ptr, qpdf);
    } else if (png_get_valid(png_ptr, info_ptr, PNG_INFO_sRGB)) {
      colorspace = create_cspace_sRGB(png_ptr, info_ptr);
    } else if (png_get_valid(png_ptr, info_ptr, PNG_INFO_cHRM)) {
      colorspace = create_cspace_CalRGB(png_ptr, info_ptr);
    } else {
      colorspace = QPDFObjectHandle::newName("/DeviceRGB");
    }
    break;

  case PNG_COLOR_TYPE_GRAY:
  case PNG_COLOR_TYPE_GRAY_ALPHA:
    if (png_get_valid(png_ptr, info_ptr, PNG_INFO_iCCP)) {
      colorspace = create_cspace_ICCBased(png_ptr, info_ptr, qpdf);
    } else if (png_get_valid(png_ptr, info_ptr, PNG_INFO_sRGB)) {
      colorspace = create_cspace_sRGB(png_ptr, info_ptr);
    } else if (png_get_valid(png_ptr, info_ptr, PNG_INFO_cHRM)) {
      colorspace = create_cspace_CalGray(png_ptr, info_ptr);
    } else {
      colorspace = QPDFObjectHandle::newName("/DeviceGray");
    }
    break;
  default:
    assert(1);
    break;
  }

  return  colorspace;
}

// Soft-Mask: stream
// Stream object want QPDF object...
static QPDFObjectHandle
create_soft_mask (png_structp png_ptr, png_infop info_ptr,
                  png_bytep image_data_ptr,
                  png_uint_32 width, png_uint_32 height, QPDF& qpdf)
{
  png_bytep   smask_data_ptr;
  png_bytep   trans;
  int         num_trans;

  if (!png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS) ||
      !png_get_tRNS(png_ptr, info_ptr, &trans, &num_trans, NULL)) {
    std::cerr << "PNG does not have valid tRNS chunk!" << std::endl;
    return  QPDFObjectHandle::newNull();
  }

  QPDFObjectHandle smask = QPDFObjectHandle::newStream(&qpdf);
  QPDFObjectHandle dict  = smask.getDict();
  smask_data_ptr = new png_byte[width * height];
  dict.replaceKey("/Type",    QPDFObjectHandle::newName("/XObjcect"));
  dict.replaceKey("/Subtype", QPDFObjectHandle::newName("/Image"));
  dict.replaceKey("/Width"  , QPDFObjectHandle::newInteger(width));
  dict.replaceKey("/Height" , QPDFObjectHandle::newInteger(height));
  dict.replaceKey("/ColorSpace", QPDFObjectHandle::newName("/DeviceGray"));
  dict.replaceKey("/BitsPerComponent", QPDFObjectHandle::newInteger(8));

  for (png_uint_32 i = 0; i < width * height; i++) {
    png_byte idx = image_data_ptr[i];
    smask_data_ptr[i] = (idx < num_trans) ? trans[idx] : 0xff;
  }
  std::string raster =
      std::string(reinterpret_cast<const char *>(smask_data_ptr),
                  width * height);
  smask.replaceStreamData(raster,
                          QPDFObjectHandle::newNull(),
                          QPDFObjectHandle::newNull());
  delete smask_data_ptr;

  return smask;
}

static QPDFObjectHandle
strip_soft_mask (png_structp png_ptr, png_infop info_ptr,
                 // next two values will be modified.
                 png_bytep image_data_ptr, png_uint_32p rowbytes_ptr,
                 png_uint_32 width, png_uint_32 height, QPDF& qpdf)
{
  png_byte color_type, bpc;

  color_type = png_get_color_type(png_ptr, info_ptr);
  bpc        = png_get_bit_depth (png_ptr, info_ptr);
  std::cerr << bpc << std::endl;
  if (color_type & PNG_COLOR_MASK_COLOR) {
    int bps = (bpc == 8) ? 4 : 8;
    if (*rowbytes_ptr != bps*width*sizeof(png_byte)) { // Something wrong
      std::cerr << "Inconsistent rowbytes value.";
      return  QPDFObjectHandle::newNull();
    }
  } else {
    int bps = (bpc == 8) ? 2 : 4;
    if (*rowbytes_ptr != bps*width*sizeof(png_byte)) { // Something wrong
      std::cerr << "Inconsistent rowbytes value.";
      return  QPDFObjectHandle::newNull();
    }
  }

  QPDFObjectHandle smask = QPDFObjectHandle::newStream(&qpdf);
  QPDFObjectHandle dict  = smask.getDict();
  dict.replaceKey("/Type",    QPDFObjectHandle::newName("/XObjcect"));
  dict.replaceKey("/Subtype", QPDFObjectHandle::newName("/Image"));
  dict.replaceKey("/Width"  , QPDFObjectHandle::newInteger(width));
  dict.replaceKey("/Height" , QPDFObjectHandle::newInteger(height));
  dict.replaceKey("/ColorSpace", QPDFObjectHandle::newName("/DeviceGray"));
  dict.replaceKey("/BitsPerComponent", QPDFObjectHandle::newInteger(bpc));

  png_bytep smask_data_ptr = new png_byte[(bpc / 8) * width * height];

  switch (color_type) {
  case PNG_COLOR_TYPE_RGB_ALPHA:
    if (bpc == 8) {
      for (png_uint_32 i = 0; i < width * height; i++) {
        memmove(image_data_ptr+(3*i), image_data_ptr+(4*i), 3);
        smask_data_ptr[i] = image_data_ptr[4*i+3];
      }
      *rowbytes_ptr = 3 * width * sizeof(png_byte);
    } else {
      for (png_uint_32 i = 0; i < width * height; i++) {
        memmove(image_data_ptr+(6*i), image_data_ptr+(8*i), 6);
        smask_data_ptr[2*i]   = image_data_ptr[8*i+6];
        smask_data_ptr[2*i+1] = image_data_ptr[8*i+7];
      }
      *rowbytes_ptr = 6 * width * sizeof(png_byte);
    }
    break;
  case PNG_COLOR_TYPE_GRAY_ALPHA:
    if (bpc == 8) {
      for (png_uint_32 i = 0; i < width*height; i++) {
        image_data_ptr[i] = image_data_ptr[2*i];
        smask_data_ptr[i] = image_data_ptr[2*i+1];
      }
      *rowbytes_ptr = width * sizeof(png_byte);
    } else {
      for (png_uint_32 i = 0; i < width * height; i++) {
        image_data_ptr[2*i]   = image_data_ptr[4*i];
        image_data_ptr[2*i+1] = image_data_ptr[4*i+1];
        smask_data_ptr[2*i]   = image_data_ptr[4*i+2];
        smask_data_ptr[2*i+1] = image_data_ptr[4*i+3];
      }
      *rowbytes_ptr = 2 * width * sizeof(png_byte);
    }
    break;
  default:
    assert(1);
    break;
  }

  std::string raster =
      std::string(reinterpret_cast<const char *>(smask_data_ptr),
                  (bpc / 8) * width * height);
  smask.replaceStreamData(raster,
                          QPDFObjectHandle::newNull(),
                          QPDFObjectHandle::newNull());
  delete smask_data_ptr;

  return  smask;
}

static void
read_image_data (png_structp png_ptr, png_bytep dest_ptr,
                 png_uint_32 height, png_uint_32 rowbytes)
{
  png_bytepp rows_p = new png_bytep[height];
  for (png_uint_32 i = 0; i < height; i++)
    rows_p[i] = dest_ptr + (rowbytes * i);
  png_read_image(png_ptr, rows_p);
  delete rows_p;
}


// Main
static const int EXIT_ERROR = 2;
static const int EXIT_OK    = 0;

static const char* const usage = "bals!";

static void error_exit (std::string const& msg)
{
  if (!msg.empty())
    std::cerr << msg << std::endl;
  exit(EXIT_ERROR);
}

static void warn (std::string const& msg)
{
  std::cerr << msg << std::endl;
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
        else if (kv[1] == "assembly" || kv[1] == "assem")
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
  Margins           margin;
  struct permission perm = {qpdf_r3p_full, qpdf_r3m_all, true, true};
  std::string       outfile, upasswd, opasswd;
  bool              nocompress = false, linearize = false, encrypt = false,
                    use_RC4 = false, is_min_version = false;
  int               keysize = 40;

  while ((opt = getopt(argc, argv, "o:V:m:ZlEK:U:O:P:R")) != -1) {
    switch (opt) {
    case 'o': /* output file */
      outfile = std::string(optarg);
      break;
    case 'V':
      if (strlen(optarg) > 1 &&
          optarg[strlen(optarg) - 1] == '+') {
        version = std::string(optarg, optarg + strlen(optarg) - 1);
        is_min_version = true;
      } else {
        version = std::string(optarg);
        is_min_version = false;
      }
      break;
    case 'm':
      error   = optarg_parse_margins(optarg, margin);
      break;
    case 'Z':
      nocompress = true;
      break;
    case 'l':
      linearize  = true;
      break;
    case 'E':
      encrypt = true;
      break;
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
    case 'R':
      use_RC4 = true;
      break;
    case ':': case '?':
      error_exit(usage);
      break;
    }
  }

  QPDF qpdf;

  qpdf.emptyPDF();
  if (argc - optind > 2 && outfile.empty()) {
    // require -o option
    error_exit(usage);
  } else if (outfile.empty() && argc > 1) { // last filename is for output
    argc--;
    outfile = std::string(argv[argc]);
  } else if (argc - optind == 0) {
    error_exit(usage);
  }

  for ( ;!error && optind < argc; optind++) {
    std::cerr << "Proccessing file: " << argv[optind] << std::endl;
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

  QPDFWriter w(qpdf, outfile.c_str());
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
        if (version < "1.4" && !is_min_version) {
          warn("Current encryption setting requires PDF ver. >= 1.4.\n" \
               "Encryption will be disabled.");
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
        if (version < "1.5" && !is_min_version) {
          warn("Current encryption setting requires PDF ver. >= 1.5.\n" \
               "Encryption will be disabled.");
          encrypt = false;
        } else {
          w.setR4EncryptionParameters(upasswd.c_str(), opasswd.c_str(),
                                      perm.accessibility, perm.extract,
                                      perm.print, perm.modify, true, true);
        }
      } else if (keysize == 256) {
        // Unencrypt Metadata unsupported yet
        if (version < "1.7" && !is_min_version) {
          warn("Current encryption setting requires PDF ver. > 1.7.\n" \
                "Encryption will be disabled.");
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
  if (is_min_version)
    w.setMinimumPDFVersion(version);
  else {
    w.forcePDFVersion(version);
  }
  w.setStreamDataMode(nocompress ? qpdf_s_uncompress : qpdf_s_compress);
  w.setLinearization(linearize);
  std::vector<QPDFExc> warns = qpdf.getWarnings();
  std::vector<QPDFExc>::iterator it;
  w.write();
  for (it = warns.begin(); it != warns.end(); it++) {
    warn(it->getMessageDetail());
  }
  std::cerr << "Output written to \""<< outfile << "\"." << std::endl;

  return EXIT_OK;
}
