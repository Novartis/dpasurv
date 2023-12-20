# The below 'commented' code generates the logo.png stored in man/figures,
# which can be used as package logo on pkgdown page

# coord_dag <- list(
#   x = c(dNt = 3, Mt = 1.5, X = 0),
#   y = c(dNt = 1, Mt = -2, X = 1)
# )
#
# dg <- dagify(dNt ~ Mt + X, Mt ~ X, coords = coord_dag)
#
# p <- dg %>%
#   ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
#   geom_dag_node(size=10) +
#   geom_dag_edges(lineheight=10) +
#   geom_dag_text(size=8) +
#   ylim(c(-2.5, 1.5)) +
#   xlim(c(-0.5, 3.5)) +
#   theme_dag_blank() +
#   theme_transparent()
#
# sticker(p, package="dpa",
#         p_size=30, p_x = 1, p_y = 1.55, p_color="black",
#         s_x=0.98, s_y=0.73, s_width=1.7, s_height=1.3,
#         url = "https://cran.r-project.org/package=UCSCXenaTools", u_color = "white", u_size = 1,
#         h_fill="white", h_color="#1F65CC", h_size=3,
#         filename="logo.png", dpi=300, white_around_sticker = TRUE)

# The below 'commented' code generates the logo.png stored in man/figures,
# which can be used as package logo on pkgdown page

# coord_dag <- list(
#   x = c(dNt = 3, Mt = 1.5, X = 0),
#   y = c(dNt = 1, Mt = -2, X = 1)
# )
#
# dg <- ggdag::dagify(dNt ~ Mt + X, Mt ~ X, coords = coord_dag)
#
# p <- dg %>%
#   ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
#   ggdag::geom_dag_node(size=10) +
#   ggdag::geom_dag_edges(lineheight=10) +
#   ggdag::geom_dag_text(size=8) +
#   ylim(c(-2.5, 1.5)) +
#   xlim(c(-0.5, 3.5)) +
#   ggdag::theme_dag_blank() +
#   ggpubr::theme_transparent()
#
# hexSticker::sticker(p, package="dpasurv",
#                     p_size=23, p_x = 1, p_y = 1.5, p_color="black",
#                     s_x=0.98, s_y=0.73, s_width=1.7, s_height=1.3,
#                     url = "https://cran.r-project.org/package=UCSCXenaTools", u_color = "white", u_size = 1,
#                     h_fill="white", h_color="#1F65CC", h_size=3,
#                     filename="logo.png", dpi=300, white_around_sticker = TRUE)
