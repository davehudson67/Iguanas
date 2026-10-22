## ============================================================
## Dist_GompertzCalendarLT_V2.R
##
## LEFT-TRUNCATED GOMPERTZ WITH PIECEWISE CALENDAR EXPOSURES
##
## Base age-specific hazard:
##
##   h0(age) = a10 * exp[b * (age - 10)]
##
## Calendar multiplier in each calendar year y:
##
##   exp(
##     betaFeed * feed[y] +
##     betaMatthew * matthew[y] +
##     sitePeriodSlope * period_decade[y]
##   )
##
## The exposure multiplier is piecewise constant within each
## calendar year. The integrated hazard is calculated exactly
## over each calendar-year segment.
##
## The density is conditional on survival to lowerBound:
##
##   f(x | T > L)
##     = h(x) * exp[- integral_L^x h(u) du]
##
## so the post-entry density integrates to 1.
##
## Exposure after the final observed year is held fixed at the
## final-year value. This only provides a proper tail for latent
## right-censored death times; likelihood contributions up to the
## censoring date do not depend on post-study exposure assumptions.
## ============================================================


dGompertzCalendarLT <- nimbleFunction(
  run = function(
      x = double(0),
      a10 = double(0),
      b = double(0),
      tB = double(0),
      lowerBound = double(0, default = 0),
      feed = double(1),
      matthew = double(1),
      betaFeed = double(0),
      betaMatthew = double(0),
      sitePeriodSlope = double(0),
      Origin = double(0),
      periodRefYear = double(0),
      tMax = integer(0),
      log = integer(0, default = 0)
  ) {

    returnType(double(0))

    if(
      a10 <= 0 |
      b <= 0 |
      lowerBound < 0 |
      tMax < 1
    ) {
      return(NaN)
    }

    if(
      x < lowerBound
    ) {
      if(log) {
        return(-Inf)
      } else {
        return(0.0)
      }
    }

    calStart <- tB + lowerBound
    calEnd <- tB + x
    currentCal <- calStart
    H <- 0.0
    studyEndBoundary <- tMax + 1.0

    while(
      currentCal < calEnd &
      currentCal < studyEndBoundary
    ) {

      yearIdx <- floor(currentCal)

      if(yearIdx < 1) {
        yearIdx <- 1
      }

      if(yearIdx > tMax) {
        yearIdx <- tMax
      }

      nextBoundary <- floor(currentCal) + 1.0

      ## NIMBLE's compiled min() is a unary reduction rather than
      ## base R's variadic min(), so choose the minimum explicitly.
      segEnd <- calEnd

      if(
        nextBoundary <
          segEnd
      ) {
        segEnd <-
          nextBoundary
      }

      if(
        studyEndBoundary <
          segEnd
      ) {
        segEnd <-
          studyEndBoundary
      }

      ageStart <- currentCal - tB
      ageEnd <- segEnd - tB

      calendarYear <- Origin + yearIdx - 1.0

      periodDecade <-
        (
          calendarYear -
          periodRefYear
        ) /
        10.0

      eta <-
        betaFeed * feed[yearIdx] +
        betaMatthew * matthew[yearIdx] +
        sitePeriodSlope * periodDecade

      Hseg <-
        exp(eta) *
        (
          a10 /
          b
        ) *
        (
          exp(
            b *
            (
              ageEnd -
              10.0
            )
          ) -
          exp(
            b *
            (
              ageStart -
              10.0
            )
          )
        )

      H <- H + Hseg
      currentCal <- segEnd
    }

    if(
      currentCal < calEnd
    ) {

      yearIdx <- tMax
      ageStart <- currentCal - tB
      ageEnd <- calEnd - tB
      calendarYear <- Origin + tMax - 1.0

      periodDecade <-
        (
          calendarYear -
          periodRefYear
        ) /
        10.0

      eta <-
        betaFeed * feed[yearIdx] +
        betaMatthew * matthew[yearIdx] +
        sitePeriodSlope * periodDecade

      Hseg <-
        exp(eta) *
        (
          a10 /
          b
        ) *
        (
          exp(
            b *
            (
              ageEnd -
              10.0
            )
          ) -
          exp(
            b *
            (
              ageStart -
              10.0
            )
          )
        )

      H <- H + Hseg
    }

    deathCal <- tB + x
    yearIdxDeath <- floor(deathCal)

    if(yearIdxDeath < 1) {
      yearIdxDeath <- 1
    }

    if(yearIdxDeath > tMax) {
      yearIdxDeath <- tMax
    }

    calendarYearDeath <- Origin + yearIdxDeath - 1.0

    periodDecadeDeath <-
      (
        calendarYearDeath -
        periodRefYear
      ) /
      10.0

    etaDeath <-
      betaFeed * feed[yearIdxDeath] +
      betaMatthew * matthew[yearIdxDeath] +
      sitePeriodSlope * periodDecadeDeath

    logHazard <-
      log(a10) +
      b *
      (
        x -
        10.0
      ) +
      etaDeath

    logProb <- logHazard - H

    if(log) {
      return(logProb)
    } else {
      return(
        exp(logProb)
      )
    }
  }
)


rGompertzCalendarLT <- nimbleFunction(
  run = function(
      n = integer(0),
      a10 = double(0),
      b = double(0),
      tB = double(0),
      lowerBound = double(0, default = 0),
      feed = double(1),
      matthew = double(1),
      betaFeed = double(0),
      betaMatthew = double(0),
      sitePeriodSlope = double(0),
      Origin = double(0),
      periodRefYear = double(0),
      tMax = integer(0)
  ) {

    returnType(double(0))

    if(
      a10 <= 0 |
      b <= 0 |
      lowerBound < 0 |
      tMax < 1
    ) {
      return(NaN)
    }

    if(
      n != 1
    ) {
      nimPrint(
        "rGompertzCalendarLT only allows n = 1; using n = 1."
      )
    }

    remainingH <- rexp(
      1,
      rate =
        1.0
    )

    currentCal <- tB + lowerBound
    studyEndBoundary <- tMax + 1.0

    while(
      currentCal < studyEndBoundary
    ) {

      yearIdx <- floor(currentCal)

      if(yearIdx < 1) {
        yearIdx <- 1
      }

      if(yearIdx > tMax) {
        yearIdx <- tMax
      }

      ## Avoid variadic min() for NIMBLE compilation.
      nextBoundary <-
        floor(
          currentCal
        ) +
        1.0

      if(
        studyEndBoundary <
          nextBoundary
      ) {
        nextBoundary <-
          studyEndBoundary
      }

      ageStart <- currentCal - tB
      ageEnd <- nextBoundary - tB
      calendarYear <- Origin + yearIdx - 1.0

      periodDecade <-
        (
          calendarYear -
          periodRefYear
        ) /
        10.0

      eta <-
        betaFeed * feed[yearIdx] +
        betaMatthew * matthew[yearIdx] +
        sitePeriodSlope * periodDecade

      multiplier <- exp(eta)

      Hseg <-
        multiplier *
        (
          a10 /
          b
        ) *
        (
          exp(
            b *
            (
              ageEnd -
              10.0
            )
          ) -
          exp(
            b *
            (
              ageStart -
              10.0
            )
          )
        )

      if(
        remainingH <= Hseg
      ) {

        baseStart <-
          exp(
            b *
            (
              ageStart -
              10.0
            )
          )

        inside <-
          baseStart +
          (
            remainingH *
            b
          ) /
          (
            multiplier *
            a10
          )

        deathAge <-
          10.0 +
          log(inside) /
          b

        return(deathAge)
      }

      remainingH <- remainingH - Hseg
      currentCal <- nextBoundary
    }

    yearIdx <- tMax
    ageStart <- currentCal - tB
    calendarYear <- Origin + tMax - 1.0

    periodDecade <-
      (
        calendarYear -
        periodRefYear
      ) /
      10.0

    eta <-
      betaFeed * feed[yearIdx] +
      betaMatthew * matthew[yearIdx] +
      sitePeriodSlope * periodDecade

    multiplier <- exp(eta)

    baseStart <-
      exp(
        b *
        (
          ageStart -
          10.0
        )
      )

    inside <-
      baseStart +
      (
        remainingH *
        b
      ) /
      (
        multiplier *
        a10
      )

    deathAge <-
      10.0 +
      log(inside) /
      b

    return(deathAge)
  }
)


registerDistributions(
  list(
    dGompertzCalendarLT =
      list(
        BUGSdist =
          "dGompertzCalendarLT(a10, b, tB, lowerBound, feed, matthew, betaFeed, betaMatthew, sitePeriodSlope, Origin, periodRefYear, tMax)",

        Rdist =
          "dGompertzCalendarLT(a10, b, tB, lowerBound, feed, matthew, betaFeed, betaMatthew, sitePeriodSlope, Origin, periodRefYear, tMax)",

        types =
          c(
            "value = double(0)",
            "a10 = double(0)",
            "b = double(0)",
            "tB = double(0)",
            "lowerBound = double(0)",
            "feed = double(1)",
            "matthew = double(1)",
            "betaFeed = double(0)",
            "betaMatthew = double(0)",
            "sitePeriodSlope = double(0)",
            "Origin = double(0)",
            "periodRefYear = double(0)",
            "tMax = integer(0)"
          ),

        pqAvail =
          FALSE,

        range =
          c(
            0,
            Inf
          )
      )
  )
)
