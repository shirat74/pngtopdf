#ifndef __TIFFFILTER_HH__
#define __TIFFFILTER_HH__

#include <string>
#include <qpdf/QPDFObjectHandle.hh>

// Utility functions for TIFF 2 predictor filtering
class TIFFFilter
{
public:
  static QPDFObjectHandle filter (std::string& raster,
                                  int32_t width, int32_t height,
                                  int8_t bpc, int8_t num_comp);
  static void apply_filter_1_2_4 (std::string& raster,
                                  int32_t width, int32_t height,
                                  int8_t bpc, int8_t num_comp);
};

#endif // __TIFFFILTER_HH__
