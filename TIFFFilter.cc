
#include <string>
#include <cassert>

#include <qpdf/QPDFObjectHandle.hh>

#include "TIFFFilter.hh"

// Many PDF viewers seems to have broken TIFF 2 predictor support?
// Ony GhostScript and MuPDF render 4bpc grayscale image with TIFF 2 predictor
// filter applied correctly.
//
//  Acrobat Reader DC  2015.007.20033  NG
//  Adobe Acrobat X    10.1.13         NG
//  Foxit Reader       4.1.5.425       NG
//  GhostScript        9.16            OK
//  SumatraPDF(MuPDF)  v3.0            OK
//  Evince(poppler)    2.32.0.145      NG (1bit and 4bit broken)
//
void
TIFFFilter::apply_filter_1_2_4 (std::string& raster,
                                int32_t width, int32_t height,
                                int8_t bpc, int8_t num_comp)
{
  assert( bpc > 0 && bpc <= 8 );

  int32_t  rowbytes = (bpc * num_comp * width + 7) / 8;
  uint8_t  mask     = (1 << bpc) - 1;
  std::vector<uint8_t> prev(num_comp);

  // Generic routine for 1 to 16 bit.
  // It supports, e.g., 7 bpc images too.
  // Actually, it is not necessary to have 16 bit inbuf and outbuf
  // since we only need 1, 2, and 4 bit support here. 8 bit is enough.
  for (int j = 0; j < height; j++) {
    int32_t  k, l, inbits, outbits;
    uint16_t inbuf, outbuf;

    std::fill(prev.begin(), prev.end(), 0);
    inbuf = outbuf = 0; inbits = outbits = 0;
    l = k = j * rowbytes;
    for (int i = 0; i < width; i++) {
      for (int c = 0; c < num_comp; c++) {
        if (inbits < bpc) { // need more byte
          inbuf   = (inbuf << 8) | raster[l]; l++;
          inbits += 8;
        }
        uint8_t cur = (inbuf >> (inbits - bpc)) & mask;
        inbits -= bpc; // consumed bpc bits
        int8_t  sub = cur - prev[c];
        prev[c] = cur;
        if (sub < 0)
          sub += (1 << bpc);
        // Append newly filtered component value
        outbuf   = (outbuf << bpc) | sub;
        outbits += bpc;
        // flush
        if (outbits >= 8) {
          raster[k] = (outbuf >> (outbits - 8)); k++;
          outbits  -= 8;
        }
      }
    }
    if (outbits > 0)
      raster[k] = (outbuf << (8 - outbits)); k++;
  }
}

// Returns DecodeParms dictionary
// NOTICE: "raster" modified.
QPDFObjectHandle
TIFFFilter::filter (std::string& raster,
                    int32_t width, int32_t height, int8_t bpc, int8_t num_comp)
{
  std::vector<uint16_t> prev(num_comp);

  if (num_comp > 4)
    return QPDFObjectHandle::newNull(); // Not supported yet

  switch (bpc) {
  case 1: case 2: case 4:
    apply_filter_1_2_4(raster, width, height, bpc, num_comp);
    break;

  case 8:
    for (int32_t j = 0; j < height; j++) {
      std::fill(prev.begin(), prev.end(), 0);
      for (int32_t i = 0; i < width; i++) {
        int32_t pos = num_comp * (width * j + i);
        for (int c = 0; c < num_comp; c++) {
          uint8_t cur   = raster[pos+c];
          int32_t sub   = cur - prev[c];
          prev[c]       = cur;
          raster[pos+c] = sub;
        }
      }
    }
    break;

  case 16:
    for (int32_t j = 0; j < height; j++) {
      std::fill(prev.begin(), prev.end(), 0);
      for (int32_t i = 0; i < width; i++) {
        int32_t pos = 2 * num_comp * (width * j + i);
        for (int c = 0; c < num_comp; c++) {
          uint16_t cur = ((uint8_t)raster[pos+2*c])*256 +
                           (uint8_t)raster[pos+2*c+1];
          uint16_t sub  = cur - prev[c];
          prev[c]       = cur;
          raster[pos+2*c  ] = (sub >> 8) & 0xff;
          raster[pos+2*c+1] = sub & 0xff;
        }
      }
    }
    break;

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
