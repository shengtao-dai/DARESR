# Changelog

## DARESR 0.1.1

- Fixed the confidence interval returned by
  [`series_eq_estimate()`](https://shengtao-dai.github.io/DARESR/reference/series_eq_estimate.md):
  it is now the estimate -/+ 1.96 standard errors (it previously used
  1.96 divided by the standard error, so its width grew with the sample
  size).
- Added full function documentation with examples, a pkgdown site and R
  CMD check on GitHub Actions.
- Added citation information for Dai (2025).

## DARESR 0.1.0

- Initial release.
