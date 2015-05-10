#include <string>
#include <math.h>

#include <iostream>

#include "Image.hh"

// Only 8 and 16 bpc supported
Image::Image (int32_t width, int32_t height, int8_t num_comp, int8_t bpc)
{
  this->width  = width;
  this->height = height;
  this->nComps = num_comp;
  this->bpc    = bpc;

  data = new Color[width*height];
}

Image::Image (int32_t width, int32_t height, int8_t nComps, int8_t bpc,
              const std::string& raster)
{
  this->width  = width;
  this->height = height;
  this->nComps = nComps;
  this->bpc    = bpc;

  data = data_from_string(raster, width, height, nComps, bpc);
}

Image::~Image ()
{
  delete data;
}

Color *
Image::data_from_string (const std::string raster,
                         int32_t width, int32_t height,
                         int8_t nComps, int8_t bpc)
{
  Color *data     = new Color[width*height];
  size_t rowbytes = (width * nComps * bpc + 7) / 8;
  size_t expected = height * rowbytes;
  std::string str = (raster.size() < expected) ?
      raster + std::string(expected - raster.size(), 0) : raster;

  // nComps = 1 for bpc < 8 for the sake of simplicity
  if (bpc < 8 && nComps > 2)
    return  data; // FIXME; Actually fail

  switch (bpc) {
  case 1: case 2: case 4: // nComps = 1 here
    {
      static const uint8_t mask[5] = {0, 0x01, 0x03, 0, 0x0f};
      int maxval = (1 << bpc) - 1;
      int numpix = 8 / bpc; // number of dots contained in a byte
      for (int32_t j = 0; j < height; j++) {
        for (int32_t i = 0; i < width; i += numpix) {
          size_t  pos = width*j + i;
          uint8_t val = (uint8_t) str[j*rowbytes + i/numpix];
          for (int k = numpix-1; k >= 0; k--) {
            data[pos+k].v[0] =
              round(65535*((float)(val & mask[bpc]) / maxval));
            data[pos+k].n    = 1;
            val >>= bpc;
          }
        }
      }
    }
    break;
  case 8:
    for (int32_t i = 0; i < width * height; i++) {
      Color   pixel;
      int32_t pos = nComps * i;
      for (int c = 0; c < nComps; c++)
        pixel.v[c] = round(65535 * ((uint8_t) str.at(pos+c) / 255.0));
      data[i] = pixel;
    }
    break;
  case 16:
    for (int32_t i = 0; i < width * height; i++) {
      Color   pixel;
      int32_t pos = 2 * nComps * i;
      for (int c = 0; c < nComps; c++)
        pixel.v[c] = ((uint8_t) str.at(pos+2*c)) * 256  +
                                      (uint8_t) str.at(pos+2*c+1);
      data[i] = pixel;
    }
    break;
  }

  return  data;
}

// Periodic boundary...
Color
Image::getPixel (int32_t x, int32_t y) const
{
  if (x < 0) {
    while (x < 0)
      x += width;
  } else if (x >= width) {
    while (x >= width)
      x -= width;
  }
  if (y < 0) {
    while (y < 0)
      y += height;
  } else if (y >= height) {
    while (y >= height)
      y -= height;
  }
  return data[y * width + x];
}

void
Image::putPixel (int32_t x, int32_t y, Color value)
{
  if (x < 0) {
    while (x < 0)
      x += width;
  } else if (x >= width) {
    while (x >= width)
      x -= width;
  }
  if (y < 0) {
    while (y < 0)
      y += height;
  } else if (y >= height) {
    while (y >= height)
      y -= height;
  }
  data[y * width + x] = value;
}

std::string
Image::getPixelBytes () const
{
  std::string bytes;
  switch (bpc) {
  case 1: case 2: case 4:
    bytes = getPixelBytesN(bpc);
    break;
  case 8:
    bytes = getPixelBytes8();
    break;
  case 16:
    bytes = getPixelBytes16();
    break;
  }
  return bytes;
}

// For NComps = 1 and bpc = 1, 2, and 4
std::string
Image::getPixelBytesN (int N) const
{
  size_t  rowbytes = (width * nComps * N + 7) / 8;
  int     numpix   =  8 / N; // num dots per byte
  uint8_t maxval   = (1 << N) - 1;

  std::string bytes(rowbytes * height, 0); // returned value
  for (int32_t j = 0; j < height; j++) {
    for (int32_t i = 0; i < width; i += numpix) {
      size_t  pos = rowbytes*j + i/numpix;
      uint8_t val = 0;
      for (int k = 0; k < numpix; k++) {
        Color pixel = getPixel(i+k, j);
        val |=
          ((uint8_t)round(maxval*(pixel.v[0]/65535.))) << (8-N*k-N);
         // std::cerr << (pixel.v[0] ? "*" :  " ");
      }
      bytes[pos] = val;
    }
    // std::cerr << std::endl;
  }
  return bytes;
}

// Converted to bit depth 8
std::string
Image::getPixelBytes8 () const
{
  std::string bytes(width * height * nComps, 0);
  for (int32_t j = 0; j < height; j++) {
    for (int32_t i = 0; i < width; i++) {
      Color  pixel = getPixel(i, j);
      size_t pos   = nComps * (j * width + i);
      for (int c = 0; c < nComps; c++) {
        bytes[pos + c] = round(255 * (pixel.v[c] / 65535.0));
      }
    }
  }
  return bytes;
}

std::string
Image::getPixelBytes16 () const
{
  std::string bytes(2 * width * height * nComps, 0);
  for (int32_t j = 0; j < height; j++) {
    for (int32_t i = 0; i < width; i++) {
      Color  pixel = getPixel(i, j);
      size_t pos   = (bpc / 8) * nComps * (j * width + i);
      for (int c = 0; c < nComps; c++) {
        bytes[pos+2*c]   = (pixel.v[c] >> 8) & 0xff;
        bytes[pos+2*c+1] =  pixel.v[c] & 0xff;
      }
    }
  }
  return bytes;
}
