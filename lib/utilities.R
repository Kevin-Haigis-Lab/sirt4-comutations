
# Replace missing values in `x` with `replacement`.
replace_na <- function(x, replacement = FALSE) {
    y <- x
    y[is.na(y)] <- replacement
    return(y)
}