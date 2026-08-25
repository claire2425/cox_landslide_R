# 00_functions.R: shared functions used across the whole pipeline

# source this once in run_all.R, and each script also sources it on its own, so can still run any single script by itself without breaking.
# Keeping every function in one file means two scripts can't accidentally end up with slightly different versions of ahp_weights() or ahp_cr().

# AHP principal-eigenvector weights.
# A is a positive reciprocal pairwise matrix (a_ij = how much more important i is than j),
# and this returns a named weight vector that sums to 1.

ahp_weights<- function(A) {
  stopifnot(is.matrix(A), nrow(A) == ncol(A))   # A must be a square matrix
  ep<- eigen(A)
  idx<- which.max(Re(ep$values))     # take the principal eigenvalue -- don't just assume it's the first one
  v<- abs(Re(ep$vectors[, idx]))     # the eigenvector's sign is arbitrary, so abs() to keep weights positive
  w<- v / sum(v)                     # normalise so the weights sum to 1
  if (!is.null(rownames(A))) names(w)<- rownames(A)
  w
  }

# AHP consistency diagnostics: lambda_max, CI and CR.
# use the same which.max convention as ahp_weights() here, so the two functions
# can never disagree about which eigenvalue is the principal one.
ahp_cr<- function(A) {
  n<- nrow(A)
  RI<- c(0, 0, 0.58, 0.90, 1.12, 1.24, 1.32, 1.41, 1.45, 1.49)  # Saaty's random index, for n up to 10
  if (n > length(RI)) stop("ahp_cr: matrix too large; RI only defined up to n=10.")
  ep<- eigen(A)
  idx<- which.max(Re(ep$values))
  lambda_max<- Re(ep$values[idx])
  CI<- (lambda_max - n) / (n - 1)
  RI_n<- RI[n]
  CR<- if (RI_n == 0) 0 else CI / RI_n   # for n<=2 the RI is 0, so CR is defined as 0 to avoid dividing by zero
  list(lambda_max = lambda_max, CI = CI, RI = RI_n, CR = CR, n = n)
  }

# A quick sanity check that A really is a valid AHP matrix:
# diagonal should be all 1, it should be reciprocal (a_ij * a_ji = 1), and all entries positive.
check_ahp_matrix<- function(A, tol = 1e-8) {
  ok_diag<- all(abs(diag(A) - 1) < tol)
  ok_recip<- all(abs(A * t(A) - 1) < tol)
  ok_pos<- all(A > 0)
  if (!ok_diag)  warning("AHP matrix: diagonal is not all 1.")
  if (!ok_recip) warning("AHP matrix: reciprocal property a_ij * a_ji = 1 violated.")
  if (!ok_pos)   warning("AHP matrix: contains non-positive entries.")
  invisible(ok_diag && ok_recip && ok_pos)
  }

# Reclassify a raster by discrete integer code -> suitability score.
# Any code not listed in the map becomes NA by default (I want unmatched cells dropped, not guessed).
reclass_codes<- function(r, map, others_to_NA = TRUE) {
  from<- as.numeric(names(map)); to<- as.numeric(map)
  rcl<- cbind(from, to)
  other_val<- if (others_to_NA) NA else NULL
  out<- terra::classify(r, rcl = rcl, others = other_val)
  out
  }

# Resample raster r onto ref's grid, but only if they don't already line up
# (skipping the resample when they match avoids introducing needless interpolation error).
align_to<- function(r, ref, categorical = FALSE) {
  same<- terra::compareGeom(r, ref, stopOnError = FALSE, messages = FALSE)
  if (isTRUE(same)) return(r)
  method<- if (categorical) "near" else "bilinear"   # nearest for categorical layers, bilinear for continuous ones
  terra::resample(r, ref, method = method)
  }

# Test whether an sf object and a terra CRS are the same.
# I compare them as st_crs() objects rather than as WKT strings, because the raw strings
# can differ in whitespace/format even when the CRS is genuinely identical.
same_crs<- function(sf_obj, terra_crs_wkt) {
  isTRUE(sf::st_crs(sf_obj) == sf::st_crs(terra_crs_wkt))
  }

cat("00_functions.R loaded: ahp_weights, ahp_cr, check_ahp_matrix, reclass_codes, align_to, same_crs\n")

