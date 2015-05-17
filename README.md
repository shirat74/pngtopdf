# pngtopdf

This is pngtopdf - a PNG to PDF converter.
To build pngtopdf, you need zlib, libpng, lcms2, and libqpdf libraries to be installed.

Usage
=====

For single PNG to PDF file conversion:

    pngtopdf [options] input output

or multiple PNGs to a multipage PDF:

    pngtopdf [options] -o output input1 input2 ... inputn

Options:

    -o string     set output filename
    -v number     set PDF version (1.7)
    -m dimensions set margins (0 0 0 0)

    -e, -E        enable/disable encryption (disabled)
    -c, -C        enable/disable RGB to CMYK conversion (disabled)
    -f, -F        enable/disable use of Flate predictor (disabled)
    -g, -G        enable/disable gamma pre-compensation (enabled)
    -l, -L        enable/disable Linearization (disabled)
    -s, -S        enable/disable use of Object Stream (enabled)
    -z, -Z        enable/disable compression (enabled)

Options for CMYK conversion:

    -b, -B        enable/disable balck point compensation (disabled)
    -i string     Set rendering intent
                  ["perceptual", "saturation", "absolute", "relative"]
                  ("perceptual")


Options for encryption:

    -K number     encryption key size [40-256] (40)
    -R            use RC4 encryption algorithm (false)
    -O string     owner password (empty)
    -U string     user password (empty)
    -P string     set permission flags

The margins are specified as 1-4 numbers with unit (pt if ommited). When four
values are specified, they represent top, right, bottom, left margin in that
order. When three values are given, then they are taken as top, left and right
side, and bottom margin. If only two values are supplied, they are interpreted
as top-bottom and left-right margines. And single value means that all four side
has same margin. "pt", "cm", "mm", and "in" can be used as unit specification.

Access permission flags are specified by space separated key-value pairs, e.g.,
"p:none m:all e:allow", meaning disallow any kind of printing, allow all type of
modification, and allow extraction. Valid keys for permission flags are "a"
(accessibility), "e" (extract), "m" (modify), "p" (print), and allowed values
are ("allow"|"yes"|"true") or "no" for "a" and "e", "all" (allow any), "annot"
(add or modify annotations), "form" (fill in forms), "assem" (assemble), or
"none" (nothing allowed) for key "m", and "full", "low", or "none" for "p".


Planned Features
================

Resampling, XMP metadata, JSON configuration file, PDF/X and PDF/A output,
more flexible way to place images onto page, add captions?, ... and many more!


Sending Feature Request and Bug Report
======================================

Concact via GitHub account or send e-mail to shunsaku.hirata74 at gmail.com
