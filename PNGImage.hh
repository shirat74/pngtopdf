#ifndef __PNGIMAGE_HH__
#define __PNGIMAGE_HH__

#include <vector>
#include <stdio.h>

#include "Image.hh"

// This is not PNG specific so we put them here for a moment.
enum calibration_type_e
{
  calibration_none = 0,   // none -- use device dependent
  calibration_profile,
  calibration_srgb,
  calibration_matrix,     // cHRM specified possibly with gAMA
  calibration_gamma_only, // only gAMA specified
};


class PNGImage : public Image
{
public:
  enum eLoadOption {
    eLoadOptionNone      = 0x00,
    eLoadPrecomposeAlpha = 0x01,
    eLoadGammaCorrect    = 0x02, // Do gamma correction when gAMA chunk exist but
                                 // not cHRM exist.
    eLoadStrip16         = 0x04,
  };

public:
  PNGImage(int32_t width, int32_t height, int8_t nComps, int8_t bpc) :
    Image(width, height, nComps, bpc),
    isValid(true), isIndexed(false), hasMaskColor(false),
    calibrationType(calibration_none),
    dpi_x(72), dpi_y(72)
    {  };
  PNGImage(int32_t width, int32_t height, int8_t nComps, int8_t bpc,
           const std::string raster) :
      Image(width, height, nComps, bpc, raster),
      isValid(true), isIndexed(false), hasMaskColor(false),
      calibrationType(calibration_none),
      dpi_x(72), dpi_y(72)
      {  };
  // Read data from file and construct PNGImage object.
  PNGImage(const std::string filename,
           enum eLoadOption options = eLoadGammaCorrect);

  // Save to file.
  int  save(const std::string filename) const;
  // Check if load image succeeded.
  bool valid() const { return isValid; };

  enum calibration_type_e getCalibrationType() const
      { return calibrationType; };
  bool hasGamma() const { return colorSpace.hasGamma; };

  std::string getICCProfileName() const { return colorSpace.ICCP.name; };
  std::vector<unsigned char> getICCProfile() const
      { return colorSpace.ICCP.profile; };
  std::string getsRGBIntent() const { return colorSpace.sRGB; };
  float getGamma() const { return colorSpace.gamma; };
  // 8 values for chromaticity of white point ant primary phosphors
  // sorry for this misleading name, I don't know how to name it.
  std::vector<float> getChromaticity() const {
    std::vector<float> cal(8);
    cal[0] = colorSpace.calRGB.xw;
    cal[1] = colorSpace.calRGB.yw;
    cal[2] = colorSpace.calRGB.xr;
    cal[3] = colorSpace.calRGB.yr;
    cal[4] = colorSpace.calRGB.xg;
    cal[5] = colorSpace.calRGB.yg;
    cal[6] = colorSpace.calRGB.xb;
    cal[7] = colorSpace.calRGB.yb;
    return cal;
  }
  bool  hasColorKeyMask() const { return hasMaskColor; };
  Color getMaskColor() const { return maskColor; };

  bool hasPalette() const { return isIndexed; };
  std::vector<Color> getPalette() const { return palette; };

  // The following easily introduces inconsistency.
  void setCalibrationType(enum calibration_type_e type)
      { calibrationType = type; };
  void setICCProfileName(std::string name) { colorSpace.ICCP.name = name; };
  void setICCProfile(std::vector<unsigned char> profile)
      { colorSpace.ICCP.profile = profile; };
  void setsRGBIntent(std::string intent) { colorSpace.sRGB = intent; };
  void setGamma(float gamma)
      { colorSpace.gamma = gamma; colorSpace.hasGamma = true; };
  // 8 values for whitepoint and matrix
  void setChromaticity(float xw, float yw, float xr, float yr,
                       float xg, float yg, float xb, float yb) {
    colorSpace.calRGB.xw = xw;
    colorSpace.calRGB.yw = yw;
    colorSpace.calRGB.xr = xr;
    colorSpace.calRGB.yr = yr;
    colorSpace.calRGB.xg = xg;
    colorSpace.calRGB.yg = yg;
    colorSpace.calRGB.xb = xb;
    colorSpace.calRGB.yb = yb;
  }

  float getResolutionX() const { return dpi_x; };
  float getResolutionY() const { return dpi_y; };
  void  setResolutionX(float dpi) { dpi_x = dpi; };
  void  setResolutionY(float dpi) { dpi_y = dpi; };
  void  setResolution (float x, float y) { dpi_x = x; dpi_y = y; };

private:

  bool isValid;

  bool isIndexed;
  bool hasMaskColor;

  struct ICCP {
    std::string name;
    std::vector<unsigned char> profile;
  };
  struct CalRGB {
    float xw, yw, xr, yr, xg, yg, xb, yb;
  };

  // union { not working
  struct {
    std::string   sRGB;   // rendering intent stored
    struct ICCP   ICCP;
    float         gamma;  // Gama Only
    struct CalRGB calRGB; // Calibrated: with cHRM and gAMA
    bool          hasGamma;
  } colorSpace;

  // colorPalette may have alpha values from tRNS chunk
  std::vector<Color> palette;
  Color maskColor;

  enum calibration_type_e calibrationType;

  float dpi_x, dpi_y;

  static int check_for_png(FILE *fp);

};

#endif // __PNGIMAGE_HH__
