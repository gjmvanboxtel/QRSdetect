# QRSdetect-internal.R - internal or barely commented functions
# Copyright (C) 2019  Geert van Boxtel,
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# Version history
# 20190116  GvB       Initial setup for package QRSdetect
# 20190129  GvB       Added RMS function to compute Root Mean Square
#
#---------------------------------------------------------------------------------------------------------------------

#' QRSdetect-internal
#'
#' Internal or barely commented functions not exported from the namespace.
#' @name QRSdetect-internal
NULL
#'
#' \preformatted{
#'   is.scalar(x)                                      # test if x is a scalar
#'   is.posscal(x)                                     # test if x is a positive scalar
#'   is.wholenumber(x, tol = .Machine$double.eps^0.5)  # test if x is a whole number
#'   RMS(x)                                            # Root Mean Square
#' }
#' @param x object
#' @param tol tolerance
#' @author Geert van Boxtel

# test if x is a scalar
#' @rdname QRSdetect-internal
is.scalar <- function(x) is.atomic(x) && length(x) == 1L && !is.character(x) && Im(x)==0

# test if x is a positive scalar
#' @rdname QRSdetect-internal
is.posscal <- function (x) is.scalar(x) && x >= 0

# test if x is a whole number
#' @rdname QRSdetect-internal
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

# Root Mean Square
#' @rdname QRSdetect-internal
RMS <- function (x) sqrt(mean(x^2, na.rm = TRUE))
