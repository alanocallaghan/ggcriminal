#' @rdname geom_rbf
#' @importFrom ggplot2 layer has_flipped_aes GeomLine GeomRibbon flipped_names
#' @importFrom grid gList
#' @export
geom_rbf <- function(
        mapping = NULL,
        data = NULL,
        position = "identity",
        ...,
        n_rbfs = 12,
        rbf_variance = 1,
        se = TRUE,
        na.rm = FALSE,
        orientation = NA,
        show.legend = NA,
        inherit.aes = TRUE
    ) {

    params <- list(
        na.rm = na.rm,
        orientation = orientation,
        se = se,
        n_rbfs = n_rbfs,
        rbf_variance = rbf_variance,
        ...
    )

    ggplot2::layer(
        data = data,
        mapping = mapping,
        stat = StatRBF,
        geom = GeomRBF,
        position = position,
        show.legend = show.legend,
        inherit.aes = inherit.aes,
        params = params
    )
}

# #' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomRBF <- ggproto(
    "GeomRBF",
    Geom,
    setup_params = function(data, params) {
        params$flipped_aes <- ggplot2::has_flipped_aes(
            data,
            params,
            range_is_orthogonal = TRUE,
            ambiguous = TRUE
        )
        params
    },

    extra_params = c("na.rm", "orientation"),
    setup_data = function(data, params) {
        GeomLine$setup_data(data, params)
    },

    # The `se` argument is set to false here to make sure drawing the
    # geom and drawing the legend is in synch. If the geom is used by a
    # stat that doesn't set the `se` argument then `se` will be missing
    # and the legend key won't be drawn. With `se = FALSE` here the
    # ribbon won't be drawn either in that case, keeping the overall
    # behavior predictable and sensible. The user will realize that they
    # need to set `se = TRUE` to obtain the ribbon and the legend key.
    draw_group = function(data, panel_params, coord, se = FALSE, flipped_aes = FALSE) {
        ribbon <- transform(data, colour = NA)
        path <- transform(data, alpha = NA)

        ymin = flipped_names(flipped_aes)$ymin
        ymax = flipped_names(flipped_aes)$ymax
        has_ribbon <- se && !is.null(data[[ymax]]) && !is.null(data[[ymin]])

        gList(
            if (has_ribbon) GeomRibbon$draw_group(ribbon, panel_params, coord, flipped_aes = flipped_aes),
            GeomLine$draw_panel(path, panel_params, coord)
        )
    },

    draw_key = ggplot2::draw_key_smooth,

    required_aes = c("x", "y"),
    optional_aes = c("ymin", "ymax"),

    default_aes = aes(
        colour = "#3366FF",
        fill = "grey60",
        size = 1,
        linetype = 1, weight = 1, alpha = 0.4
    )
)


rbf_lm <- function(x, y, n_rbfs = 12, rbf_variance = 1.2) {
    X <- .designMat(x, n_reg_terms = n_rbfs, variance = rbf_variance)
    fit <- lm(y ~ 0 + ., data = as.data.frame(X))
    attr(fit, "rbf_variance") <- rbf_variance
    fit
}

predict_rbf <- function(model, xseq, se = TRUE, level = 0.95) {
    Xseq <- .designMat(xseq,
        n_reg_terms = length(coef(model)),
        variance = attr(model, "rbf_variance")
    )
    pred <- stats::predict(
        model,
        newdata = as.data.frame(Xseq),
        se.fit = se,
        level = level,
        interval = if (se) "confidence" else "none"
    )

    if (se) {
        fit <- as.data.frame(pred$fit)
        names(fit) <- c("y", "ymin", "ymax")
        base::data.frame(x = xseq, fit, se = pred$se.fit)
    } else {
        base::data.frame(x = xseq, y = as.vector(pred))
    }
}
















#' Smoothed conditional means using Gaussian radial basis function regression.
#' 
#' @param n_rbfs Number of Gaussian radial basis functions.
#' @param rbf_variance Variance of Gaussian radial basis functions.
#' @param mapping,data,position,se,n,span,fullrange,na.rm,orientation,show.legend,inherit.aes,... Geom parameters. See \code{\link[ggplot2]{geom_smooth}} for details.
#' @export
#' @importFrom ggplot2 ggproto Stat Geom draw_key_smooth flip_data
#' @rdname geom_rbf
stat_rbf <- function(
        mapping = NULL,
        data = NULL,
        position = "identity",
        n_rbfs = 12,
        rbf_variance = 1,
        ...,
        se = TRUE,
        n = 80,
        span = 0.75,
        fullrange = FALSE,
        level = 0.95,
        na.rm = FALSE,
        orientation = NA,
        show.legend = NA,
        inherit.aes = TRUE) {

  layer(
    data = data,
    mapping = mapping,
    stat = StatRBF,
    geom = GeomRBF,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      se = se,
      n = n,
      n_rbfs = n_rbfs,
      rbf_variance = rbf_variance,
      fullrange = fullrange,
      level = level,
      na.rm = na.rm,
      orientation = orientation,
      span = span,
      ...
    )
  )
}

# #' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
StatRBF <- ggplot2::ggproto(
    "StatRBF", ggplot2::Stat,
    setup_params = function(data, params) {
        params$flipped_aes <- ggplot2::has_flipped_aes(data, params, ambiguous = TRUE)
        params
    },

    extra_params = c("na.rm", "orientation"),

    compute_group = function(
            data,
            scales,
            n_rbfs = 12,
            rbf_variance = 1,
            se = TRUE,
            n = 80,
            span = 0.75,
            fullrange = FALSE,
            xseq = NULL,
            level = 0.95,
            na.rm = FALSE,
            flipped_aes = NA) {

        data <- ggplot2::flip_data(data, flipped_aes)
        if (length(unique(data$x)) < 2) {
            # Not enough data to perform fit
            return(data.frame())
        }

        if (is.null(data$weight)) data$weight <- 1

        if (is.null(xseq)) {
        if (is.integer(data$x)) {
            if (fullrange) {
                xseq <- scales$x$dimension()
            } else {
                xseq <- sort(unique(data$x))
            }
        } else {
            if (fullrange) {
                range <- scales$x$dimension()
            } else {
                range <- range(data$x, na.rm = TRUE)
            }
            xseq <- seq(range[1], range[2], length.out = n)
        }
        }

        model <- rbf_lm(
            data$x,
            data$y,
            n_rbfs = n_rbfs,
            rbf_variance = rbf_variance
        )

        prediction <- predict_rbf(model, xseq, se, level)
        prediction$flipped_aes <- flipped_aes
        ggplot2::flip_data(prediction, flipped_aes)
    },

    required_aes = c("x", "y")
)
