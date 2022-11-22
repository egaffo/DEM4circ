# https://github.com/EmilHvitfeldt/r-color-palettes

# https://emilhvitfeldt.github.io/r-color-palettes/discrete.html
# install.packages("paletteer")
library(paletteer)
ggprism_colors20 <- paletteer_d("ggprism::colors")
pals_polychrome <- rev(paletteer_d("pals::polychrome"))[1:20]

library(Polychrome)
# show_col(Polychrome::dark.colors(n = 20))
# kelly_poly <- rev(Polychrome::kelly.colors())[1:20]
kelly_poly <- Polychrome::kelly.colors()[2:21]

# https://hughjonesd.github.io/tweaking-colours-with-the-shades-package.html
# install.packages("shades")
suppressPackageStartupMessages(library(shades))

# install.packages("jrnold/ggthemes")
# https://github.com/jrnold/ggthemes/commit/9f8772ba89c7c9f022657fee1d9c112295490cc1
old_gdocs <- c('#3366cc', '#dc3912', '#ff9900', '#109618', '#990099',
               '#0099c6', '#dd4477', '#66aa00', '#b82e2e', '#316395',
               '#994499', '#22aa99', '#aaaa11', '#6633cc', '#e67300',
               '#8b0707', '#651067', '#329262', '#5574a6', '#3b3eac',
               '#b77322', '#16d620', '#b91383', '#f4359e', '#9c5935',
               '#a9c413', '#2a778d', '#668d1c', '#bea413', '#0c5922',
               '#743411', '#3366cc')

custom_old_gdocs <- c("#000000", old_gdocs[c(     1,  5,  8,  3,
                                             17,  6,  7, 16,  9,
                                             24, 13, 15, 31, 20,
                                             14, 30,  2, 18,  4)])
custom_old_gdocs <- as.character(brightness(custom_old_gdocs, delta(.1)))
# show_col(custom_old_gdocs)

get_color_shades <-
  function(nshades, baseColor){
    if(nshades > 0){
      # colors_shades <- shades::brightness(shades = baseColor,
      #                                     values = seq(from = 0.7, to = 1, length.out = 6))
      # colors_shades <- shades::lightness(shades = baseColor,
      #                                    values = c(seq(from = 50, to = 80, length.out = 6)))
      colors_shades <- shades::saturation(shades = baseColor,
                                         values = c(seq(from = 0.5, to = 1, length.out = 6)))

      colors_shades <- c(baseColor, colors_shades[1:nshades])
    }else{
      colors_shades <- baseColor
    }
    colors_shades
  }


# https://medialab.github.io/iwanthue/
# install.packages('hues')
library(hues)
# iwanthue

iwanthue_intense20_1 <-
  c("#edb7cf",
    "#f32600",
    "#0064f7",
    "#e4c447",
    "#c700bb",
    "#b1d262",
    "#111b69",
    "#59dda7",
    "#ff7cf0",
    "#008345",
    "#ff4c4d",
    "#006aa7",
    "#bf7200",
    "#cca3ff",
    "#0f4000",
    "#a2003d",
    "#cdbf9a",
    "#580037",
    "#ff8e7f",
    "#5d362f")
