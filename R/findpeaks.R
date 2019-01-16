# findpeaks.R - implementation of Matlab/Octave Signal toolbox 'findpeaks' function
# Copyright (C) 2019  Geert van Boxtel
# The Octave function is Copyright (C) 2012 Juan Pablo Carbajal <carbajal@ifi.uzh.ch>
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
#
#---------------------------------------------------------------------------------------------------------------------

#' Find local maxima
#'
#' Find local maxima (peaks) in an input signal vector.
#'
#' This function searches for peaks in a signal vector. Peaks of a positive array of data are defined as local maxima.
#' For double-sided data, they are maxima of the positive part and minima of the negative part. The function provides
#' various options to search for peaks in noisy data, such specifying a minimum peak height (\code{minh}), a minimum distance
#' between peaks (\code{mind}), and a minimum or maximum width of the peaks (\code{minw} and \code{maxw}).\cr\cr
#' The width of the peaks is estimated using a parabola fitted to the neighborhood of each peak. The width is caulculated
#' with the formula \eqn{a * (width - x0)^2 = 1}, where \eqn{a} is the the concavity of the parabola and \eqn{x0} its vertex.
#' The neighborhood size is equal to the value of \code{mind}.
#'
#' @param data    The input signal vector.
#' @param minh    Minimum peak height (non-negative scalar). Only peaks that exceed this value will be returned.
#'                For data taking positive and negative values use the option \code{ds}.
#'                Default: machine precision (.Machine$double.eps).
#' @param mind    Minimum separation between peaks (positive integer). Peaks separated by less than this distance are considered
#'                a single peak. This distance is also used to fit a second order polynomial to the peaks to estimate their width,
#'                (see \code{details}, therefore it acts as a smoothing parameter. Default value 1.
#' @param minw    Minimum width of peaks (positive integer). Default value 1.
#' @param maxw    Maximum width of peaks (positive integer). Default value \code{Inf}.
#' @param ds      Double-sided. Tells the function that data takes positive and negative values. The base-line for the peaks
#'                is taken as the mean value of the function. This is equivalent as passing the absolute value of the data
#'                after removing the mean.
#' @return When called with \code{minw = 0} and \code{maxw = Inf}, this function returns a \code{\link{list}} containing
#'         two values:
#'   \describe{
#'     \item{\code{$pks}: }{array containing the value of \code{data} at the peaks}
#'     \item{\code{$idx}: }{array containing the peak indices}
#'   }
#'   When called with either \code{minw > 0} or \code{maxw < Inf}, then the returned \code{\link{list}} contains these additional
#'   variables:
#'   \describe{
#'     \item{\code{$parabol}: }{a \code{\link{list}} containing additional information about the parabol fitted to the peak.
#'     The \code{\link{list}} \code{$pp} contains the coefficients of the 2nd degree polynomial (\code{a}, \code{b}, and \code{b2}),
#'     and \code{$x} the extrema of the interval where it was fitted (\code{$from}, \code{to}).}
#'     \item{\code{$height}:}{The estimated height of the returned peaks (in units of \code{data}).}
#'     \item{\code{$baseline}: }{The height at which the roots of the returned peaks were calculated (in units of \code{data}).}
#'     \item{\code{$roots}: }{The abscissa values (in index units) at which the parabola fitted to each of the returned peaks
#'       realizes its width.}
#'   }
#' @references \href{https://octave.sourceforge.io/signal/}{Octave Signal package}
#' @author Geert van Boxtel
#' @examples
#' # Example 1: Finding the peaks of smooth data is not a big deal
#'
#' t <- 2*pi*seq(0,1,length=1024)
#' y <- sin(3.14*t) + 0.5*cos(6.09*t) + 0.1*sin(10.11*t+1/6) + 0.1*sin(15.3*t+1/3)
#'
#' data1 <- abs(y) # Positive values
#' peaks1 <- findpeaks(data1)
#'
#' data2 <- y # Double-sided
#' peaks2 <- findpeaks(data2, ds=TRUE)
#' peaks3 <- findpeaks (data2, ds=TRUE, minh=0.5)
#'
#' \dontrun{
#'   op <- par(mfrow=c(1,2))
#'   plot(t,data1,type="l", xlab="", ylab="")
#'   points (t[peaks1$idx],peaks1$pks,col="red", pch=1)
#'   plot(t,data2,type="l", xlab="", ylab="")
#'   points (t[peaks2$idx],peaks2$pks,col="red", pch=1)
#'   points (t[peaks3$idx],peaks3$pks,col="red", pch=4)
#'   legend ("topleft", '0: >2*sd, x: >0.5', bty="n", text.col="red")
#'   par (op)}
#'
#' # Example 2: Noisy data may need tuning of the parameters. In this example,
#' # "mind" is used as a smoother of the peaks.
#'
#' t <- 2*pi*seq(0,1,length=1024)
#' y <- sin(3.14*t) + 0.5*cos(6.09*t) + 0.1*sin(10.11*t+1/6) + 0.1*sin(15.3*t+1/3)
#' data <- abs(y + 0.1*rnorm(length(y),1)); # Positive values + noise
#' peaks1 <- findpeaks(data, minh=1)
#' dt <- t[2]-t[1]
#' peaks2 <- findpeaks(data, minh=1, mind=round(0.5/dt))
#'
#' \dontrun{
#'   op <- par(mfrow=c(1,2))
#'   plot(t, data, type="l", xlab="", ylab="")
#'   points (t[peaks1$idx],peaks1$pks,col="red", pch=1)
#'   plot(t, data, type="l", xlab="", ylab="")
#'   points (t[peaks2$idx],peaks2$pks,col="red", pch=1)
#'   par (op)}
#'
#' @export

findpeaks <- function (data, minh = .Machine$double.eps, mind = 1, minw = 1, maxw = Inf, ds = FALSE) {

  # check function arguments
  ld <- length(data)
  if (!(is.vector (data) || is.array(data) || class(data) == 'ts') || ld < 3)
    stop ("findpeaks: data must be a vector of at least 3 elements")
  if (!is.posscal(minh)) stop ("findpeaks: minh must be a positive scalar")
  if (!is.posscal(mind)) stop ("findpeaks: mind must be a positive scalar")
  if (!is.posscal(minw)) stop ("findpeaks: minw must be a positive scalar")
  if (!is.logical(ds)) stop ("findpeaks: ds should a a logical value T/F")

  # always work with recified data
  .data <- abs(data-mean(data))
  if (ds) {
    tmp <- data
    data <- .data
    .data <- tmp
  } else {
    if (min(data,na.rm=T) < 0) warning("findpeaks: ds is FALSE but data contain negative values")
  }

  # Rough estimates of first and second derivative
  df1 <- diff (data, differences=1)[c(1,1:(length(data)-1))]
  df2 <- diff (data, differences=2)[c(1,1,1:(length(data)-2))]

  # check for changes of sign of 1st derivative and negativity of 2nd derivative.
  # <= in 1st derivative includes the case of oversampled signals.
  idx <- which(df1*c(df1[2:length(df1)],0) <= 0 & c(df2[2:length(df2)],0) < 0)

  # Get peaks that are beyond given height
  tf  <- which(data[idx] > minh)
  idx <- idx[tf]
  if (length(idx) <= 0) return (NULL)

  # sort according to magnitude
  tmp <- sort(data[idx],decreasing=T, index=T)
  idx.s <- idx[tmp$ix]

  ## Treat peaks separated less than mind as one
  D <- with(expand.grid(A=idx.s,B=t(idx.s)), abs(A-B))
  dim(D) <- c(length(idx.s),length(idx.s))
  if (any(D) < mind) {
    i <- 1
    peak <- NULL
    node2visit <- 1:length(idx.s)
    visited <- NULL
    idx.pruned <- idx.s
    while (length(node2visit) > 0) {
      d <- D[node2visit[1],]
      visited <- c(visited, node2visit[1])
      node2visit <- node2visit[-1]
      neighs <- setdiff(which(d < mind), visited)
      if (length(neighs) > 0) {
        idx.pruned <- setdiff (idx.pruned, idx.s[neighs])
        visited <- c(visited, neighs)
        node2visit <- setdiff (node2visit, visited)
      }
    }
    idx <- idx.pruned
  }


  # If minw or maxw are specified, then do the following
  # Estimate widths of peaks and filter for:
  # width smaller than given.
  # wrong concavity.
  # not high enough
  # data at peak is lower than parabola by 1%
  if (minw > 0 || maxw < Inf) {

    extra.x <- extra.pp <- extra.roots <- extra.height <- extra.baseline <- data.frame()

    idx.pruned <- idx
    n  <- length (idx)
    for (i in 1:n) {
      ind <- round(max(idx[i]-mind/2,1)) : round (min(idx[i]+mind/2,ld))
      pp <- stats::coef(stats::lm(data[ind] ~ ind + I(ind^2)))
      H  <- pp[1] - pp[2]^2/(4*pp[3])
      rz <- try(rev(Re(polyroot(c(pp[1]-mean(c(H,minh)), pp[2], pp[3])))), silent =TRUE)
      width <- try(abs (diff (rz)), silent=TRUE)
      if (!is.na(H) && !is.na(pp[3]) && (width < minw || width > maxw || pp[3] > 0 || H < minh || data[idx[i]] < 0.99*H)) {
        idx.pruned = setdiff (idx.pruned, idx[i])
      } else {
        extra.x <- rbind(extra.x, ind[c(1, length(ind))])
        extra.pp = rbind(extra.pp, rev(pp))
        if (class(rz)=="try-error") {
          extra.roots <- rbind(extra.roots, c(NA,NA))
        } else {
          extra.roots <- rbind(extra.roots, rz)
        }
        extra.height <- rbind(extra.height, H)
        extra.baseline <- rbind(extra.baseline, mean (c(H,minh)))
      }
    }
    idx <- idx.pruned
  }

  # check for double sided
  if (ds) pks <- .data[idx]
  else pks <- data[idx]

  # return values
  if (minw > 0 || maxw < Inf) {
    colnames(extra.x) <- c('from','to')
    colnames(extra.pp) <- c('b2','b', 'a')
    colnames(extra.roots) <- c('a0','a1')
    return (list(pks=pks, idx=idx, parabol=list(x=as.list(extra.x),pp=as.list(extra.pp)),
           height=as.numeric(extra.height[,1]), baseline=as.numeric(extra.baseline[,1]), roots=as.list(extra.roots)))
  } else {
    return (list(pks=pks, idx=idx))
  }
}


# #demo
# t <- 2*pi*seq(0,1,length=1024)
# y <- sin(3.14*t) + 0.5*cos(6.09*t) + 0.1*sin(10.11*t+1/6) + 0.1*sin(15.3*t+1/3)
#
# data1 <- abs(y) # Positive values
# peaks1 <- findpeaks(data1)
#
# data2 <- y # Double-sided
# peaks2 <- findpeaks(data2,ds=T)
# peaks3 <- findpeaks (data2, ds=T, minh=0.5)
#
# op <- par(mfrow=c(1,2))
# plot(t,data1,type="l", xlab="", ylab="")
# points (t[peaks1$idx],peaks1$pks,col="red", pch=1)
# plot(t,data2,type="l", xlab="", ylab="")
# points (t[peaks2$idx],peaks2$pks,col="red", pch=1)
# points (t[peaks3$idx],peaks3$pks,col="red", pch=4)
# legend ("topleft", '0: >2*sd, x: >0.5', bty="n", text.col="red")
# par (op)
#----------------------------------------------------------------------------
# Finding the peaks of smooth data is not a big deal!

# #demo
# t <- 2*pi*seq(0,1,length=1024)
# y <- sin(3.14*t) + 0.5*cos(6.09*t) + 0.1*sin(10.11*t+1/6) + 0.1*sin(15.3*t+1/3)
# data <- abs(y + 0.1*rnorm(length(y),1)); # Positive values + noise
# peaks1 <- findpeaks(data, minh=1)
# dt <- t[2]-t[1]
# peaks2 <- findpeaks(data, minh=1, mind=round(0.5/dt))
# op <- par(mfrow=c(1,2))
# plot(t, data, type="l", xlab="", ylab="")
# points (t[peaks1$idx],peaks1$pks,col="red", pch=1)
# plot(t, data, type="l", xlab="", ylab="")
# points (t[peaks2$idx],peaks2$pks,col="red", pch=1)
# par (op)
#----------------------------------------------------------------------------
# Noisy data may need tuning of the parameters. In the 2nd example,
# mind is used as a smoother of the peaks.

# findpeaks (c(1, 1, 1))
# findpeaks (t(c(1, 1, 1)))

## Test for bug #45056
## Test input vector is an oversampled sinusoid with clipped peaks
# x <- pmin (3, cos (2*pi*c(0:8000) / 600) + 2.01)
# findpeaks(x)

## Test input validation
# findpeaks ()
# findpeaks (1)
# findpeaks (c(1, 2))

