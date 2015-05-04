# pngtopdf

This is pngtopdf - a PNG to PDF converter.
To build pngtopdf, you need libz, libpng, and libqpdf libraries to be installed.

Usage

    pngtopdf [options] input output
  
or

    pngtopdf [options] -o output input1 input2 ... inputn

options:

- -o string     set output filename
- -V number     set PDF version
- -m dimensions set margins
 
 The margins are specified as 1-4 numbers with unit (pt if ommited).
 When four values are specified, then they represents top, right,
 bottom, left margin in that order. When three values are given,
 then they are taken as top, left and right side, and bottom margin.
 If only two values are supplied, they are interpreted as top-bottom
 and left-right margines. And single value means that all four side
 has same margin. "pt", "cm", "mm", and "in" can be used as unit
 specification.
