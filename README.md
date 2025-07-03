# Root Finding Methods Comparison

## Description
Comprehensive comparison of numerical methods for finding equation roots. Implements and analyzes Bisection, Newton's method (both analytic and numeric derivatives), and Secant method across three different test functions.

## Features
- **Bisection Method** - Guaranteed convergence for continuous functions
- **Newton's Method (Analytic)** - Fast convergence using exact derivatives
- **Newton's Method (Numeric)** - Uses finite difference approximation for derivatives
- **Secant Method** - Newton's method without requiring derivatives
- **Automatic Root Detection** - Systematically searches for all roots in interval
- **Error Analysis** - Compares results with Wolfram reference solutions
- **Tabular Output** - Professional formatting with convergence statistics

## Mathematical Background
- **Bisection**: Iteratively halves interval using intermediate value theorem
- **Newton's**: x_{n+1} = x_n - f(x_n)/f'(x_n), quadratic convergence
- **Secant**: Approximates derivative using secant line, superlinear convergence
- **Numeric Derivative**: f'(x) ≈ [f(x+h) - f(x-h)]/(2h)

## Test Functions
1. **f₁(x) = ln(x+1) - 1/(x+3)** - Logarithmic/rational function
2. **f₂(x) = x³ + 30cos(x) - 1/x** - Polynomial/trigonometric/rational mix
3. **f₃(x) = sin(3πx)/(x+2) + 1/(x+4)** - Oscillatory function with singularities

## Usage
Program automatically searches interval [-3, 4] for roots, applies all methods, and compares accuracy/efficiency against Wolfram reference solutions.
