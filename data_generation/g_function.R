
g_indicator <- function(s, t) as.numeric(s >= t)

g_relu <- function(s, t) pmax(s - t, 0)

# Active g function, change by calling set_g_type("relu") or set_g_type("indicator")
g <- g_indicator

set_g_type <- function(type = c("indicator", "relu")) {
  type <- match.arg(type)
  g <<- switch(type, indicator = g_indicator, relu = g_relu)
  invisible(g)
}

g_type_int <- function() if (identical(g, g_indicator)) 0L else 1L

