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
#include <fstream>
#include <cassert>

#include <lcms2.h>

#include "Image.hh"
#include "PNGImage.hh"

static struct
{
  struct {
    std::string  version;
  } PDF;

  struct {
    std::string  defaultRGBProfilePath;
    std::string  sRGBProfilePath;
    std::string  CMYKProfilePath;
  } colorManagement;
} config = {
  { "1.7" },
  {"sRGB Color Space Profile.icm", "JapanColor2001Coated.icc"}
};

#include <qpdf/QPDF.hh>
#include <qpdf/QPDFExc.hh>
#include <qpdf/QPDFObjectHandle.hh>
#include <qpdf/QPDFWriter.hh>
#include <qpdf/QUtil.hh>

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

static std::pair<std::string, std::string> strip_soft_mask (PNGImage& src);

static QPDFObjectHandle
newStreamFromFile (QPDF& qpdf, const std::string filename)
{
  std::ifstream ifs(filename.c_str(), std::ios::binary | std::ios::in);
  std::string   data((std::istreambuf_iterator<char>(ifs)),
                      std::istreambuf_iterator<char>());
  return QPDFObjectHandle::newStream(&qpdf, data);
}

int
png_include_image (QPDF& qpdf, std::string filename, Margins margin)
{
  int32_t     page_width, page_height;
  std::string raster_cmyk, raster_rgb, raster_mask;
  bool        has_smask = false, use_cmyk = false;

  PNGImage src(filename);

  page_width  = 72 * src.getWidth()  / src.getResolutionX()
                + margin.left + margin.right ;
  page_height = 72 * src.getHeight() / src.getResolutionY()
                + margin.top  + margin.bottom;

  if (src.getNComps() == 2 || src.getNComps() == 4) {
    std::pair<std::string, std::string> data = strip_soft_mask(src);
    raster_rgb  = data.first;
    raster_mask = data.second;
    has_smask   = true;
  } else {
    raster_rgb = src.getPixelBytes();
    has_smask  = false;
  }


  // Color transform with LittleCMS library
  cmsHPROFILE hInProfile = NULL;
  switch (src.getColorSpaceType()) {
  case png_colorspace_iccp:
    {
      std::vector<unsigned char> profile = src.getICCProfile();
      hInProfile = cmsOpenProfileFromMem(profile.data(), profile.size());
    }
    break;
  case png_colorspace_srgb:
    {
      // first try loading sRGB ICC Profile from file.
      if (!config.colorManagement.sRGBProfilePath.empty()) {
        hInProfile =
            cmsOpenProfileFromFile(
              config.colorManagement.sRGBProfilePath.c_str(), "r");
      }
      if (!hInProfile) {
        hInProfile = cmsCreate_sRGBProfile();
        // Use approximate sRGB
      }
    }
    break;
  case png_colorspace_calibrated:
    {
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
  case png_colorspace_gamma_only:
    {
      // Assume D65
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
            config.colorManagement.defaultRGBProfilePath.c_str(), "r");
    } else {
      hInProfile = cmsCreate_sRGBProfile();
    }
    break;
  }
  // Fallback to built-in sRGB
  if (hInProfile == NULL)
    hInProfile = cmsCreate_sRGBProfile();

  cmsHPROFILE hOutProfile =
    cmsOpenProfileFromFile(config.colorManagement.CMYKProfilePath.c_str(), "r");
  cmsHTRANSFORM hTransform =
    cmsCreateTransform(
        hInProfile,  src.getBPC() == 8 ? TYPE_RGB_8  : TYPE_RGB_16,
        hOutProfile, src.getBPC() == 8 ? TYPE_CMYK_8 : TYPE_CMYK_16,
        INTENT_PERCEPTUAL, 0);
  cmsCloseProfile(hInProfile);
  cmsCloseProfile(hOutProfile);

  if (hTransform) {
    size_t pixel_size  = src.getWidth() * src.getHeight();
    size_t data_size   = (src.getBPC() / 8) * 4 * pixel_size;
    char *outputBuffer = new char[data_size];
    cmsDoTransform(hTransform, raster_rgb.data(), outputBuffer, pixel_size);
    raster_cmyk = std::string(outputBuffer, data_size);
    delete outputBuffer;
    cmsDeleteTransform(hTransform);
    use_cmyk = true;
  }

  // Creating an image XObject.
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

  QPDFObjectHandle colorspace;
  if (use_cmyk) {
    colorspace = QPDFObjectHandle::newArray();
    QPDFObjectHandle profile =
        newStreamFromFile(qpdf, config.colorManagement.CMYKProfilePath);
    profile.getDict().replaceKey("/N", QPDFObjectHandle::newInteger(4));
    colorspace.appendItem(QPDFObjectHandle::newName("/ICCBased"));
    colorspace.appendItem(qpdf.makeIndirectObject(profile));
  } else {
    colorspace = QPDFObjectHandle::newName("/DeviceRGB"); // FIXME
  }
  image_dict.replaceKey("/ColorSpace", colorspace);
  // Soft-mask
  if (has_smask) {
    QPDFObjectHandle smask = QPDFObjectHandle::newStream(&qpdf);
    QPDFObjectHandle dict  = smask.getDict();
    dict.replaceKey("/Type",    QPDFObjectHandle::newName("/XObjcect"));
    dict.replaceKey("/Subtype", QPDFObjectHandle::newName("/Image"));
    dict.replaceKey("/Width"  , QPDFObjectHandle::newInteger(src.getWidth()));
    dict.replaceKey("/Height" , QPDFObjectHandle::newInteger(src.getHeight()));
    dict.replaceKey("/ColorSpace", QPDFObjectHandle::newName("/DeviceGray"));
    dict.replaceKey("/BitsPerComponent",
                    QPDFObjectHandle::newInteger(src.getBPC()));
    smask.replaceStreamData(raster_mask,
                            QPDFObjectHandle::newNull(),
                            QPDFObjectHandle::newNull());
    image_dict.replaceKey("/SMask", smask);
  }

  image.replaceStreamData(use_cmyk ? raster_cmyk : raster_rgb,
                          QPDFObjectHandle::newNull(),
                          QPDFObjectHandle::newNull());

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
      "q "   + QUtil::double_to_string(72 * src.getWidth()  / src.getResolutionX()) + " " +
      "0 0 " + QUtil::double_to_string(72 * src.getHeight() / src.getResolutionY()) + " " +
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
  qpdf.addPage(page, false);

  return 0;
}

#if 0
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
#endif

// returns <image data without alpha channel, soft mask data>
static std::pair<std::string, std::string>
strip_soft_mask (PNGImage& src)
{
  // We must be sure that image has alpha channel.
  // Bit depth must be either 8 or 16 here.
  assert( src.getNComps() == 2 || src.getNComps() == 4 );
  assert( src.getBPC() == 8 || src.getBPC() == 16 );

  size_t size = (src.getBPC() / 8) * src.getWidth() * src.getHeight();
  char  *smask_data_ptr = new char[size];
  char  *image_data_ptr = new char[size*(src.getNComps()-1)];

  switch (src.getNComps()) {
  case 4:
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
  }

  std::string raster =
      std::string(image_data_ptr, size * (src.getNComps() - 1));
  std::string smask = std::string(smask_data_ptr, size);

  delete smask_data_ptr;
  delete image_data_ptr;

  return  std::make_pair(raster, smask);
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
        config.PDF.version = std::string(optarg, optarg + strlen(optarg) - 1);
        is_min_version = true;
      } else {
        config.PDF.version = std::string(optarg);
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
        if (config.PDF.version < "1.4" && !is_min_version) {
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
        if (config.PDF.version < "1.5" && !is_min_version) {
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
        if (config.PDF.version < "1.7" && !is_min_version) {
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
    w.setMinimumPDFVersion(config.PDF.version);
  else {
    w.forcePDFVersion(config.PDF.version);
  }
  w.setStreamDataMode(nocompress ? qpdf_s_uncompress : qpdf_s_compress);
  w.setLinearization(linearize);
  std::vector<QPDFExc> warns = qpdf.getWarnings();
  std::vector<QPDFExc>::iterator it;
  w.write();
  for (it = warns.begin(); it != warns.end(); it++) {
    warn(it->getMessageDetail());
  }

  return EXIT_OK;
}
