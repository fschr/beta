package beta

import "testing"

func TestBeta(t *testing.T) {
	pairs := []struct {
		x        float64
		y        float64
		expected float64
	}{
		{x: 1.5, y: 0.2, expected: 4.477609374347168},
		{x: 3.9, y: 4.2, expected: 0.006662554701084894},
	}
	for i, pair := range pairs {
		got := Beta(pair.x, pair.y)
		if got != pair.expected {
			t.Errorf("#%d: expected %v, got %v", i, pair.expected, got)
		}
	}
}

func TestBetaInc(t *testing.T) {
	pairs := []struct {
		a float64
		b float64
		x float64
		expected float64
	}{
		{a: 1.0, b: 3.0, x: 0.4, expected: 0.784},
	}
	for i, pair := range pairs {
		got := BetaInc(pair.a, pair.b, pair.x)
		if got != pair.expected {
			t.Errorf("#%d: expected %v, got %v", i, pair.expected, got)
		}
	}
}
