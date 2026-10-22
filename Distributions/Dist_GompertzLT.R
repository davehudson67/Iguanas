## ============================================================
## CORRECT LEFT-TRUNCATED GOMPERTZ DISTRIBUTION FOR NIMBLE
##
## T is age at death from birth:
##
##   h(t) = a * exp(b * t)
##
## Conditioning on survival to lowerBound L gives:
##
##   f(t | T > L) = f(t) / S(L),  t >= L
##
## so
##
##   log f(t | T > L)
##     = log(a) + b*t
##       - (a/b) * [exp(b*t) - exp(b*L)]
##
## This is NOT the same as simply replacing t by t-L while
## leaving a unchanged.
## ============================================================


dGompertzLT <- nimbleFunction(
  run = function(
      x = double(0),
      a = double(0),
      b = double(0),
      lowerBound = double(0, default = 0),
      log = integer(0, default = 0)
  ) {

    returnType(double(0))

    if(
      a <= 0 |
      b <= 0 |
      lowerBound < 0
    ) {
      return(NaN)
    }

    if(
      x <
      lowerBound
    ) {
      if(log) {
        return(-Inf)
      } else {
        return(0.0)
      }
    }

    bx <-
      b *
      x

    bL <-
      b *
      lowerBound

    ## Numerical guard. Values this extreme have effectively
    ## zero density for the parameter ranges relevant here.
    if(
      bx >
      700
    ) {
      if(log) {
        return(-Inf)
      } else {
        return(0.0)
      }
    }

    exp_bx <-
      exp(
        bx
      )

    exp_bL <-
      exp(
        bL
      )

    logProb <-
      log(
        a
      ) +
      bx -
      (
        a /
        b
      ) *
      (
        exp_bx -
        exp_bL
      )

    if(log) {
      return(logProb)
    } else {
      return(
        exp(
          logProb
        )
      )
    }
  }
)


rGompertzLT <- nimbleFunction(
  run = function(
      n = integer(0),
      a = double(0),
      b = double(0),
      lowerBound = double(0, default = 0)
  ) {

    returnType(double(0))

    if(
      a <= 0 |
      b <= 0 |
      lowerBound < 0
    ) {
      return(NaN)
    }

    if(
      n !=
      1
    ) {
      nimPrint(
        "rGompertzLT only allows n = 1; using n = 1."
      )
    }

    u <-
      runif(
        1,
        0,
        1
      )

    exp_bL <-
      exp(
        b *
        lowerBound
      )

    inside <-
      exp_bL -
      (
        b /
        a
      ) *
      log(
        1 -
        u
      )

    out <-
      log(
        inside
      ) /
      b

    return(out)
  }
)


pGompertzLT <- nimbleFunction(
  run = function(
      q = double(0),
      a = double(0),
      b = double(0),
      lowerBound = double(0, default = 0),
      lower.tail = integer(0, default = 1),
      log.p = integer(0, default = 0)
  ) {

    returnType(double(0))

    if(
      a <= 0 |
      b <= 0 |
      lowerBound < 0
    ) {
      return(NaN)
    }

    if(
      q <
      lowerBound
    ) {

      if(lower.tail) {

        if(log.p) {
          return(-Inf)
        } else {
          return(0.0)
        }

      } else {

        if(log.p) {
          return(0.0)
        } else {
          return(1.0)
        }
      }
    }

    bq <-
      b *
      q

    bL <-
      b *
      lowerBound

    if(
      bq >
      700
    ) {

      if(lower.tail) {

        if(log.p) {
          return(0.0)
        } else {
          return(1.0)
        }

      } else {

        if(log.p) {
          return(-Inf)
        } else {
          return(0.0)
        }
      }
    }

    logS <-
      -(
        a /
        b
      ) *
      (
        exp(
          bq
        ) -
        exp(
          bL
        )
      )

    if(
      !lower.tail
    ) {

      if(log.p) {
        return(logS)
      } else {
        return(
          exp(
            logS
          )
        )
      }

    } else {

      cdf <-
        1 -
        exp(
          logS
        )

      if(log.p) {
        return(
          log(
            cdf
          )
        )
      } else {
        return(cdf)
      }
    }
  }
)


qGompertzLT <- nimbleFunction(
  run = function(
      p = double(0),
      a = double(0),
      b = double(0),
      lowerBound = double(0, default = 0),
      lower.tail = integer(0, default = 1),
      log.p = integer(0, default = 0)
  ) {

    returnType(double(0))

    if(
      a <= 0 |
      b <= 0 |
      lowerBound < 0
    ) {
      return(NaN)
    }

    if(log.p) {
      p <-
        exp(
          p
        )
    }

    if(
      !lower.tail
    ) {
      p <-
        1 -
        p
    }

    if(
      p <=
      0
    ) {
      return(lowerBound)
    }

    if(
      p >=
      1
    ) {
      return(Inf)
    }

    exp_bL <-
      exp(
        b *
        lowerBound
      )

    inside <-
      exp_bL -
      (
        b /
        a
      ) *
      log(
        1 -
        p
      )

    out <-
      log(
        inside
      ) /
      b

    return(out)
  }
)


registerDistributions(
  list(
    dGompertzLT = list(
      BUGSdist =
        "dGompertzLT(a, b, lowerBound)",
      Rdist =
        "dGompertzLT(a, b, lowerBound)",
      pqAvail =
        TRUE,
      range =
        c(
          0,
          Inf
        )
    )
  )
)
