# pngtopdf

This is pngtopdf - a PNG to PDF converter.
To build pngtopdf, you need zlib, libpng, and libqpdf libraries to be installed.

Usage

    pngtopdf [options] input output

or

    pngtopdf [options] -o output input1 input2 ... inputn

options:

    -o string     set output filename
    -V number     set PDF version (1.7)
    -m dimensions set margins (0 0 0 0)
    -l            linearize (false)
    -E            enable encryption (false)
    -K number     encryption key size [40-256] (40)
    -R            use RC4 encryption algorithm (false)
    -Z            disable compression (false)
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
