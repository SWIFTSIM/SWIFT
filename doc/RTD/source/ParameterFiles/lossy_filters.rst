.. Lossy compression filters

.. _Compression_filters:

Compression Filters
~~~~~~~~~~~~~~~~~~~

Filters to compress the data in snapshots can be applied to reduce the
disk footprint of the datasets. The filters provided by SWIFT are
filters natively provided by HDF5, implying that the library will
automatically and transparently apply the reverse filter when reading
the data stored on disk. They can be applied in combination with, or
instead of, the lossless gzip compression filter.

**These compression filters are lossy, meaning that they modify the
data written to disk**

*The filters will reduce the accuracy of the data stored. No check is
made inside SWIFT to verify that the applied filters make sense. Poor
choices can lead to all the values of a given array reduced to 0, Inf,
or to have lost too much accuracy to be useful. The onus is entirely
on the user to choose wisely how they want to compress their data.*

The filters are not applied when using parallel-hdf5.

The available filters are listed below.

N-bit filters for long long integers
------------------------------------

The N-bit filter takes a `long long` and saves only the most
significant N bits.

This can be used in cases similar to the particle IDs. For instance,
if they cover the range :math:`[1, 10^{10}]` then 64-bits is too many
and a lot of disk space is wasted storing the 0s. In this case
:math:`\left\lceil{\log_2(10^{10})}\right\rceil + 1 = 35` bits are
sufficient (The extra "+1" is for the sign bit).

SWIFT implements 5 variants of this filter:

 * ``Nbit36`` stores the 36 most significant bits (Numbers up to
   :math:`3.4\times10^{10}`, comp. ratio: 1.78)
 * ``Nbit40`` stores the 40 most significant bits (Numbers up to
   :math:`5.4\times10^{11}`, comp. ratio: 1.6)
 * ``Nbit44`` stores the 44 most significant bits (Numbers up to
   :math:`8.7\times10^{12}`, comp. ratio: 1.45)
 * ``Nbit48`` stores the 48 most significant bits (Numbers up to
   :math:`1.4\times10^{14}`, comp. ratio: 1.33)
 * ``Nbit56`` stores the 56 most significant bits (Numbers up to
   :math:`3.6\times10^{16}`, comp. ratio: 1.14)

Note that if the data written to disk is requiring more than the N
bits then part of the information written to the snapshot will
lost. SWIFT does not apply any verification before applying the
filter.

Scaling filters for floating-point numbers
------------------------------------------

The D-scale filters can be used to round floating-point values to a fixed
*absolute* accuracy.

They start by computing the minimum of an array that is then deducted from
all the values. The array is then multiplied by :math:`10^n` and truncated
to the nearest integer. These integers are stored with the minimal number
of bits required to store the values. That process is reversed when reading
the data.

For an array of values

+--------+--------+-------+
|  1.2345| -0.1267| 0.0897|
+--------+--------+-------+

and :math:`n=2`, we get stored on disk (but hidden to the user):

+--------+--------+-------+
|    136 |      0 |     22|
+--------+--------+-------+

This can be stored with 8 bits instead of the 32 bits needed to store the
original values in floating-point precision, realising a gain of 4x.

When reading the values (for example via ``h5py`` or ``swiftsimio``), that
process is transparently reversed and we get:

+--------+--------+-------+
|  1.2333| -0.1267| 0.0933|
+--------+--------+-------+

Using a scaling of :math:`n=2` hence rounds the values to two digits after
the decimal point.

SWIFT implements 4 variants of this filter:

 * ``DScale1`` scales by :math:`10^1`
 * ``DScale2`` scales by :math:`10^2`
 * ``DScale3`` scales by :math:`10^3`
 * ``DScale6`` scales by :math:`10^6`

An example application is to store the positions with ``pc`` accuracy in
simulations that use ``Mpc`` as their base unit by using the ``DScale6``
filter.

The compression rate of these filters depends on the data. On an
EAGLE-like simulation, compressing the positions from ``Mpc`` to ``pc`` (via
``Dscale6``) leads to rate of around 2.2x.

Modified floating-point representation filters
----------------------------------------------

These filters modify the bit-representation of floating point numbers
to get a different *relative* accuracy.

In brief, floating point (FP) numbers are represented in memory as
:math:`(\pm 1)\times a \times 2^b` with a certain number of bits used to store each
of :math:`a` (the mantissa) and :math:`b` (the exponent) as well as
one bit for the overall sign [#f1]_.  For example, a standard 4-bytes
`float` uses 23 bits for :math:`a` and 8 bits for :math:`b`. The
number of bits in the exponent mainly drives the range of values that
can be represented whilst the number of bits in the mantissa drives
the relative accuracy of the numbers.

Converting to the more familiar decimal notation, we get that the
number of decimal digits that are correctly represented is
:math:`\log_{10}(2^{n(a)+1})`, with :math:`n(x)` the number of bits in
:math:`x`. The range of positive numbers that can be represented is
given by :math:`[2^{-2^{n(b)-1}+2}, 2^{2^{n(b)-1}}]`. For a standard
`float`, this gives a relative accuracy of :math:`7.2` decimal digits
and a representable range of :math:`[1.17\times 10^{-38}, 3.40\times
10^{38}]`. Numbers above the upper limit are labeled as `Inf` and
below this range they default to zero.

The filters in this category change the number of bits in the mantissa and
exponent. When reading the values (for example via ``h5py`` or
``swiftsimio``) the numbers are transparently restored to regular ``float``
but with 0s in the bits of the mantissa that were not stored on disk, hence
changing the result from what was stored originally before compression.

These filters offer a fixed compression ratio and a fixed relative
accuracy. The available options in SWIFT are:


+-----------------+--------------+--------------+-------------+---------------------------------------------------+-------------------+
| Filter name     | :math:`n(a)` | :math:`n(b)` | Accuracy    | Range                                             | Compression ratio |
+=================+==============+==============+=============+===================================================+===================+
| No filter       | 23           | 8            | 7.22 digits | :math:`[1.17\times 10^{-38}, 3.40\times 10^{38}]` | ---               |
+-----------------+--------------+--------------+-------------+---------------------------------------------------+-------------------+
| ``FMantissa13`` | 13           | 8            | 4.21 digits | :math:`[1.17\times 10^{-38}, 3.40\times 10^{38}]` | 1.45x             |
+-----------------+--------------+--------------+-------------+---------------------------------------------------+-------------------+
| ``FMantissa9``  | 9            | 8            | 3.01 digits | :math:`[1.17\times 10^{-38}, 3.40\times 10^{38}]` | 1.78x             |
+-----------------+--------------+--------------+-------------+---------------------------------------------------+-------------------+
| ``BFloat16``    | 7            | 8            | 2.41 digits | :math:`[1.17\times 10^{-38}, 3.40\times 10^{38}]` | 2x                |
+-----------------+--------------+--------------+-------------+---------------------------------------------------+-------------------+
| ``HalfFloat``   | 10           | 5            | 3.31 digits | :math:`[6.1\times 10^{-5}, 6.5\times 10^{4}]`     | 2x                |
+-----------------+--------------+--------------+-------------+---------------------------------------------------+-------------------+

The accuracy given in the table corresponds to the number of decimal digits
that can be correctly stored. The "no filter" row is displayed for
comparison purposes.

The first two filters are useful to keep the same range as a standard
`float` but with a reduced accuracy of 3 or 4 decimal digits. The last two
are the two standard reduced precision options fitting within 16 bits: one
with a much reduced relative accuracy and one with a much reduced
representable range.

An example application is to store the densities with the ``FMantissa9``
filter as we rarely need more than 3 decimal digits of accuracy for this
quantity.

------------------------

.. [#f1] Note that the representation in memory of FP numbers is more
	 complicated than this simple picture. See for instance this
	 `Wikipedia
	 <https://en.wikipedia.org/wiki/Single-precision_floating-point_format>`_
	 article.

	    
