## morphology.R
##
## Morphological transformations with R
## (c) 2004 Vincent Zoonekynd
##
## Distributed under the GPL v2 or later
##
## (This work was mainly written before my arrival at BioRet,
## http://www.bioret.com/)
##

############################################################

# Translation
# Translates the image "x" by "i" units horizontally and "j" units
# vertically. The pixels that leave the image limits are discarded.
# The newly exposed pixels are set to "zero".
translate <- function (x,i,j,zero=0) {
  n <- dim(x)[1]
  m <- dim(x)[2]
  while (i>0) {
    x <- rbind( rep(zero,m), x[-n,] )
    i <- i - 1
  }
  while (i<0) {
    x <- rbind( x[-1,], rep(zero,m) )
    i <- i + 1
  }
  while (j>0) {
    x <- cbind( rep(zero,n), x[,-m] )
    j <- j - 1
  }
  while (j<0) {
    x <- cbind( x[,-1], rep(zero,n) )
    j <- j + 1
  }
  x
}

############################################################

# Structuring elements
new.structuring.element <- function (m, xmin, xmax, ymin, ymax) {
  m <- t(m)
  stopifnot( dim(m)[1] == xmax - xmin + 1 )
  stopifnot( dim(m)[2] == ymax - ymin + 1 )
  r <- list()
  r$xrange <- seq(xmin,xmax)
  r$yrange <- seq(ymin,ymax)
  r$value <- function (i,j) { m[1+i-xmin,1+j-ymin] }
  r$m <- m    # useless...
  class(r) <- "structuring.element"
  r
}
structuring.element.plus <-
  new.structuring.element( matrix(c(F,T,F,T,T,T,F,T,F),nr=3,nc=3), -1, 1, -1, 1)
structuring.element.times <-
  new.structuring.element( matrix(c(T,F,T,F,T,F,T,F,T),nr=3,nc=3), -1, 1, -1, 1)
structuring.element.square <-
  new.structuring.element( matrix(T,nr=3,nc=3), -1, 1, -1, 1)
structuring.element.disc <- function (radius) {
  r <- ceiling(radius)
  n <- 2*r+1
  m <- matrix(F, nr=n, nc=n)
  x <- (n+1)/2
  m <- (col(m)-x)^2 + (row(m)-x)^2 <= radius^2
  new.structuring.element(m, -r,r,-r,r)
}
structuring.element.square <- function (radius=1) {
  r <- ceiling(radius)
  n <- 2*r+1
  m <- matrix(F, nr=n, nc=n)
  x <- (n+1)/2
  m <- sup(abs(col(m)-x), abs(row(m)-x)) <= radius
  new.structuring.element(m, -r,r,-r,r)
}

############################################################

# Hit-and-Miss Transform (French: transformation tout-ou-rien)
TTR <- function (x, s) {
  res <- matrix(T, nr=dim(x)[1], nc=dim(x)[2])
  for (i in s$xrange) {
    for (j in s$yrange) {
      if(!is.na(s$value(i,j))) {
        if(s$value(i,j)) {
          res <- res & translate(x,i,j)
        } else {
          res <- res & !translate(x,i,j)
        }
      }
    }
  }
  res
}
se.rotation <- function (s) {
  M <- s$m
  n <- dim(M)[1]
  m <- dim(M)[2]
  y <- range(s$xrange)
  x <- range(-s$yrange)
  new.structuring.element(t(M)[m:1,], x[1],x[2],y[1],y[2])
}
rotations <- function (s) {
  res <- list()
  res <- append(res, list(s))
  s <- se.rotation(s)
  res <- append(res, list(s))
  s <- se.rotation(s)
  res <- append(res, list(s))
  s <- se.rotation(s)
  res <- append(res, list(s))
  res
}
TTRr <- function (x,s) {
  res <- matrix(F, nr=dim(x)[1], nc=dim(x)[2])
  for (ss in rotations(s)) {
    res <- res | TTR(x,ss)
  }
  res
}
iTTRr <- function (x,s) { x & ! TTRr(x,s) }

############################################################

# The "sup" and "inf" operations are the equivalent, for grayscale
# images, of the "union" and "intersection" operations for binary
# images.
sup <- function (x,y) {
  ifelse(x>y,x,y)
}
inf <- function (x,y) {
  ifelse(x>y,y,x)

}
complementaire <- function (image) {
  max(image) - image
}
prive.de <- function (a,b) {
  inf(a, complementaire(b))
}

############################################################

# Dilatation and erosion (here, for grayscale images) are the basic
# operations, from which you can define all the others.

dilatation <- function (image, kernel) {
  n <- dim(image)[1]
  m <- dim(image)[2]
  result <- image
  for (k in kernel$xrange)
    for (l in kernel$yrange)
      if (!is.na(kernel$value(k,l)) && kernel$value(k,l))
        result <- sup(result,translate(image,k,l))
  result
}
erosion <- function (image, kernel) {
  complementaire(dilatation(complementaire(image), kernel))
}

############################################################

# Opening and closing
ouverture <- function (image, kernel) {
  dilatation( erosion(image, kernel), kernel )
}
fermeture <- function (image, kernel) {
  erosion( dilatation(image, kernel), kernel )
}
gradient <- function (image, kernel) {
  dilatation(image,kernel) & ! erosion(image,kernel)
}
gradient.interne <- function (image, kernel) {
  prive.de(image, erosion(image,kernel))
}
gradient.externe <- function (image, kernel) {
  prive.de(dilatation(image,kernel), image)
}
laplacien <- function (image) {
  n <- dim(image)[1]
  m <- dim(image)[2]
  result <- matrix(FALSE, nr=n, nc=m)
  for (i in 2:(n-1)) {
    for (j in 2:(m-1)) {
      result[i,j] <- (
        4*image[i,j]
        - image[i-1,j] - image[i,j-1] - image[i+1,j] - image[i,j+1]
        ) != 0
    }
  }
  result
}
# Top hat
chapeau.haut.de.forme <- function (image, kernel) {
  prive.de( fermeture(image, kernel), image)
}
chapeau.haut.de.forme.conjugue <- function (image, kernel) {
  prive.de( image, ouverture(image, kernel) )
}
reconstruction <- function (x,y) { # x is binary, y need not be
  xx <- x*y
  repeat {
    x <- xx
    xx <- inf(  dilatation(x,structuring.element.plus),  y  )
    if(all(x==xx)) break
  }
  x
}
# Hysteresis threshold
seuillage.par.hysteresis <- function (image,a,b) {
  reconstruction( image >= b, image >= a )
}
# Hole filling
bouchage.de.trous <- function (x) {
  n <- dim(x)[1]
  m <- dim(x)[2]
  M <- matrix(F, nr=n, nc=m)
  M[1,] <- T
  M[,1] <- T
  M[n,] <- T
  M[,m] <- T
  !reconstruction(M,!x)
}

# Local maxima
maxima.locaux.4 <- function (d) {
  n <- dim(d)[1]
  m <- dim(d)[2]
  d.N <- rbind( d[-1,], rep(0,m) )
  d.S <- rbind( rep(0,m), d[-n,] )
  d.W <- cbind( d[,-1], rep(0,n) )
  d.E <- cbind( rep(0,n), d[,-m] )
  (d>d.N) & (d>d.E) & (d>d.W) & (d>d.S)
}
maxima.locaux.8 <- function (d) {
  n <- dim(d)[1]
  m <- dim(d)[2]
  t.N <- function (d) rbind( d[-1,], rep(0,dim(d)[2]) )
  t.S <- function (d) rbind( rep(0,dim(d)[2]), d[-dim(d)[1],] )
  t.W <- function (d) cbind( d[,-1], rep(0,dim(d)[1]) )
  t.E <- function (d) cbind( rep(0,dim(d)[1]), d[,-dim(d)[2]] )
  (d>t.N(d)) & (d>t.S(d)) & (d>t.E(d)) & (d>t.W(d)) &
  (d>t.N(t.E(d))) & (d>t.N(t.W(d))) &
  (d>t.S(t.E(d))) & (d>t.S(t.W(d)))
}
max.loc.hys <- function (x, r, h) {
  res <- matrix(T, nr=dim(x)[1], nc=dim(x)[2])
  for (k in r$xrange) {
    for (l in r$yrange) {
      res <- res & (translate(x,k,l)<x+h)
    }
  }
  res
}

# Regional maxima
dilatation.geodesique <- function (x, y, r) {
  for (i in 1:r) {
    # On prend la 4-distance et la 8-distance une fois sur deux :
    # Áa correspond a une boule octogonale.
    if ( i %% 2 == 0 ) {
      x <- dilatation(x, structuring.element.plus) & y
    } else {
      x <- dilatation(x, structuring.element.square()) & y
    }
  }
  x
}
r.h.maxima <- function (x,r,h) {
  res <- matrix(F, nr=dim(x)[1], nc=dim(x)[2])
  for (n in 1:255) {
    res <- res | (
      (x>n) & ! dilatation.geodesique(x>n+h, x>n, r)
    )
  }
  res
}

# A skeleton
MB.iteration <- function (x) {
  a <- translate(x,-1,0) & translate(x,0,1)
  a <- a & translate(a,1,-1) & x
  b <- (translate(a,1,0)  & ! translate(x,-1,0) ) |
       (translate(a,-1,0) & ! translate(x,1,0)  ) |
       (translate(a,0,1)  & ! translate(x,0,-1) ) |
       (translate(a,0,-1) & ! translate(x,0,1)  )
  c <- xor(x, translate(x,0,1))
  c <- c & translate(c,-1,0) & xor(x, translate(x,-1,0))
  c <- c | translate(c,0,-1)
  c <- c | translate(c,1,0)
  x <- x & ! ( b & ! c )
  x
}
MB <- function (x) {
  xx <- x
  repeat {
    x <- xx
    xx <- MB.iteration(x)
    if (all(x==xx)) break
  }
  x
}
# After computing a skeleton, we are often left with a lot of
# unwanted "hairs".
coupe <- function (x) {
  a <- ( x & translate(x,1,0) & !translate(x,-1,0) &
           !translate(x,0,1) & !translate(x,0,-1) &
           !translate(x,1,1) & !translate(x,1,-1) &
           !translate(x,-1,1) & !translate(x,-1,-1) ) |
       ( x & !translate(x,1,0) & !translate(x,-1,0) &
           translate(x,0,1) & !translate(x,0,-1) &
           !translate(x,1,1) & !translate(x,1,-1) &
           !translate(x,-1,1) & !translate(x,-1,-1) ) |
       ( x & !translate(x,1,0) & !translate(x,-1,0) &
           !translate(x,0,1) & !translate(x,0,-1) &
           translate(x,1,1) & !translate(x,1,-1) &
           !translate(x,-1,1) & !translate(x,-1,-1) ) |
       ( x & !translate(x,1,0) & !translate(x,-1,0) &
           !translate(x,0,1) & !translate(x,0,-1) &
           !translate(x,1,1) & !translate(x,1,-1) &
           translate(x,-1,1) & !translate(x,-1,-1) ) |
       ( x & !translate(x,1,0) & translate(x,-1,0) &
           !translate(x,0,1) & !translate(x,0,-1) &
           !translate(x,1,1) & !translate(x,1,-1) &
           !translate(x,-1,1) & !translate(x,-1,-1) ) |
       ( x & !translate(x,1,0) & !translate(x,-1,0) &
           !translate(x,0,1) & translate(x,0,-1) &
           !translate(x,1,1) & !translate(x,1,-1) &
           !translate(x,-1,1) & !translate(x,-1,-1) ) |
       ( x & !translate(x,1,0) & !translate(x,-1,0) &
           !translate(x,0,1) & !translate(x,0,-1) &
           !translate(x,1,1) & translate(x,1,-1) &
           !translate(x,-1,1) & !translate(x,-1,-1) ) |
       ( x & !translate(x,1,0) & !translate(x,-1,0) &
           !translate(x,0,1) & !translate(x,0,-1) &
           !translate(x,1,1) & !translate(x,1,-1) &
           !translate(x,-1,1) & translate(x,-1,-1) )
  x & ! a
}

# Homotopic kernel
axe.median <- function (x, kernel=structuring.element.plus) {
  n <- dim(x)[1]
  m <- dim(x)[2]
  r <- matrix(F, nr=n, nc=m)
  while(any(x)) {
    y <- erosion(x, kernel)
    r <- r | ( x & ! dilatation(y,kernel) )
    x <- y
  }
  r
}
# (Figure 3.11.1 page 96 from Manzanera's thesis)
noyau.homotopique.iteration <- function (x, zero) {
  y <- x & ! translate(x,0,1,zero)
  x <- x & ( y | translate(y,1,0,zero) | translate(y,-1,0,zero) | translate(x,0,-1,zero) )
   y <- x & ! translate(x,1,0,zero)
  x <- x & ( y | translate(y,0,1,zero) | translate(y,0,-1,zero) | translate(x,-1,0,zero) )
   y <- x & ! translate(x,0,-1,zero)
  x <- x & ( y | translate(y,1,0,zero) | translate(y,-1,0,zero) | translate(x,0,1,zero) )
   y <- x & ! translate(x,-1,0,zero)
  x <- x & ( y | translate(y,0,-1,zero) | translate(y,0,1,zero) | translate(x,1,0,zero) )
  x
}
noyau.homotopique <- function (x, zero=0) {
  xx <- x
  repeat {
    x <- xx
    xx <- noyau.homotopique.iteration(x, zero)
    if (all(xx==x)) break
  }
  x
}

# Printing functions
plot.structuring.element <- function (e, col='grey', ...) {
  plot.new()
  plot.window(xlim=range(e$xrange)+c(-.5,.5),
              ylim=range(e$yrange)+c(-.5,.5))
  for (x in e$xrange) {
    for (y in e$yrange) {
      v <- e$value(x,y)
      if (!is.na(v)) {
        if (v) {
          polygon(x+c(-.5,.5,.5,-.5), y+c(-.5,-.5,.5,.5), col=col, ...)
        } else {
          polygon(x+c(-.5,.5,.5,-.5), y+c(-.5,-.5,.5,.5), ...)
        }
      }
    }
  }
}

plot.image <- function (im, col=rainbow(255), border=T, ...) {
  stopifnot(length(dim(im))==2)
  n <- dim(im)[1]
  m <- dim(im)[2]
  N <- length(col)
  im <- floor(  N*(1+im-min(im))/( 1 + max(im) - min(im) )  )
  plot.new()
  plot.window(xlim=c(.5,n+.5), ylim=c(.5,m+.5))
  for (x in 1:n) {
    for (y in 1:m) {
      if (border) {
        polygon(x+c(-.5,.5,.5,-.5), y+c(-.5,-.5,.5,.5), col=col[im[x,y]], ...)
      } else {
        polygon(x+c(-.5,.5,.5,-.5), y+c(-.5,-.5,.5,.5), col=col[im[x,y]], border=NA,...)
      }
    }
  }
}

plot3d.image <- function (x,
        theta=30, phi=35,
        expand=.1,
        shade = .75, ltheta=120,
        border=NA, box=F
) {
  nx <- dim(x)[1]
  ny <- dim(x)[2]
  y <- 0.25 * (x[-nx,-ny] + x[-1,-ny] + x[-nx,-1] + x[-1,-1])
  y <- (y-min(x))/(max(y)-min(y))
  persp(x,
        col=heat.colors(255)[floor(254*y+1)],
        theta=theta, phi=phi,
        expand=expand,
        shade = shade, ltheta=ltheta,
        border=border, box=box)
}
