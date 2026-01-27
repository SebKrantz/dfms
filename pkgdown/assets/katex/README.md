# KaTeX Local Files

This directory contains local copies of KaTeX CSS and JS files to avoid loading them from the internet.

## Files:

- `katex.min.css` - ✅ KaTeX stylesheet
- `katex.min.js` - ✅ KaTeX JavaScript library

These files are automatically copied to `docs/katex/` when building the pkgdown site and loaded via `_pkgdown.yml` configuration.

## Note on fonts:

The CSS file references fonts with relative paths (`fonts/KaTeX_*.woff2`). The fonts will fall back to system fonts if not available, which is acceptable for most use cases. If you want to include fonts locally, you would need to download the entire KaTeX distribution and place the `fonts/` directory here as well.
