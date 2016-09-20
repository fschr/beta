package beta

import "math"

// BetaInc is the normalized incomplete beta function.
// This function has undefined behavior when
// a < 0.0 || b < 0.0.
func BetaInc(a, b, x float64) float64 {
	if x < 0.0 || x > 1.0 {
		panic("betaInc: x < 0.0 || x > 1.0")
	}
	if x == 1.0 {
		return 1.0
	}

	lnPreVal := -math.Log(Beta(a, b)) + a*math.Log(x) + b*math.Log(1.0-x)

	if x < (a+1.0)/(a+b+2.0) {
		/* Apply continued fraction directly. */
		return math.Exp(lnPreVal) * betaContFrac(a, b, x) / a
	} else {
		/* Apply continued fraction after hypogeometric transformation. */
		return 1.0 - math.Exp(lnPreVal)*betaContFrac(b, a, 1.0-x)/b
	}
}

// Beta is the beta function
func Beta(a, b float64) float64 {
	return math.Gamma(a) * math.Gamma(b) / math.Gamma(a+b)
}


// consts carried over from the GSL port, used in betaContFrac
const (
	GSL_DBL_MIN     = 2.2250738585072014e-308
	GSL_DBL_EPSILON = 2.2204460492503131e-16
)

// betaContFrac is a function that helps calculate
// the normalized/regularized incomplete beta function.
// It is so named because it applies continued fractions.
func betaContFrac(a, b, x float64) float64 {
	const maxIter uint = 512                 /* control iterations */
	const cutoff float64 = 2.0 * GSL_DBL_MIN /* control the zero cutoff */
	iterCount := uint(0)
	var cf float64

	/* standard initialization for continued fraction */
	numTerm := 1.0
	denTerm := 1.0 - (a+b)*x/(a+1.0)
	if math.Abs(denTerm) < cutoff {
		denTerm = cutoff
	}
	denTerm = 1.0 / denTerm
	cf = denTerm

	for iterCount < maxIter {
		k := float64(iterCount) + 1.0
		coeff := k * (b - k) * x / (((a - 1.0) + 2.0*k) * (a + 2.0*k))
		var deltaFrac float64

		/* first step */
		denTerm = 1.0 + coeff*denTerm
		numTerm = 1.0 + coeff/numTerm
		if math.Abs(denTerm) < cutoff {
			denTerm = cutoff
		}
		if math.Abs(numTerm) < cutoff {
			numTerm = cutoff
		}
		denTerm = 1.0 / denTerm
		deltaFrac = denTerm * numTerm
		cf *= deltaFrac
		coeff = -(a + k) * (a + b + k) * x / ((a + 2*k) * (a + 2*k + 1.0))

		/* second step */
		denTerm = 1.0 + coeff*denTerm
		numTerm = 1.0 + coeff/numTerm
		if math.Abs(denTerm) < cutoff {
			denTerm = cutoff
		}
		if math.Abs(numTerm) < cutoff {
			numTerm = cutoff
		}
		denTerm = 1.0 / denTerm
		deltaFrac = denTerm * numTerm
		cf *= deltaFrac
		if math.Abs(deltaFrac-1.0) < 2.0*GSL_DBL_EPSILON {
			break
		}

		iterCount += 1
	}
	return cf
}
