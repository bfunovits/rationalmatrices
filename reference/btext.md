# Write "Boxed" Text into the Margins of a Plot.

This tool writes text into a "colored" box in one of the margins of the
plot region.

## Usage

``` r
btext(
  text,
  side = 3,
  cex = 1,
  col = "black",
  size = cex,
  bg = "lightgray",
  border = "lightgray",
  parse_text = FALSE
)
```

## Arguments

- text:

  character or [`expression`](https://rdrr.io/r/base/expression.html).

- side:

  on which side of the plot (1=bottom, 2=left, 3=top, 4=right).

- cex:

  character expansion factor.

- col:

  color to use (for the text).

- size:

  determines the size of the box (in "line" units)

- bg, border:

  background and border color of the box.

- parse_text:

  if yes the procedure tries to coerce the text into an
  [`expression`](https://rdrr.io/r/base/expression.html).

## Value

(invisible) vector with the "box" corner and center coordinates.

## See also

[`mtext`](https://rdrr.io/r/graphics/mtext.html) and
[`plotmath`](https://rdrr.io/r/grDevices/plotmath.html) for details on
mathematical annotation.

## Examples

``` r
if (FALSE) { # \dontrun{
plot(1:10, axes = FALSE)
graphics::box()
btext('hallo', side = 1)
btext('Sigma[11]', side = 2, bg = 'lightblue')
btext('Sigma[11]', side = 3, bg = NA, border = 'black', size = 1.2, parse_text = TRUE)
btext(expression(Sigma[11]), side = 4, bg = 'orange', border = 'black', cex = 1.5)
} # }
```
