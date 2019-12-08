library(shiny)
library(NLMR)
library(landscapetools)
library(landscapemetrics)
library(dplyr)
library(raster)
library(rgdal)
library(mobsim)
library(vegan)

cols <- rev(terrain.colors(255))[c(1,255)]

plot_map <- function(raster_map){
  
  
  par(mar = c(4,4,1,1), cex.lab = 1.4)
  
  if (max(raster_map[]) == 0){
    image(raster_map, asp = 1, xlab = "x [m]", ylab = "y [m]", col = cols[1])
    box()
    # plot(raster_map, col = cols[1], legend = F)
  } else if (min(raster_map[]) == 1) {
    image(raster_map, asp = 1, xlab = "x [m]", ylab = "y [m]", col = cols[2])
    box()
    # plot(raster_map, col = cols[255], legend = F)
  } else {
    # plot(raster_map, legend = F)
    image(raster_map, asp = 1, xlab = "x [m]", ylab = "y [m]", col = cols)
    box()
  }
}

ui <- fluidPage(
  
  titlePanel("Simulation of geometric fragmentation effects"),
  
  fluidRow(
    column(4,
           inputPanel(
             h4("Landscape parameter"),
             numericInput("hab_amount", "Habitat amount [0 - 100%]", 50,
                      min = 0, max = 100, step = 5),
             numericInput("hab_rough", "Fragmentation per se [0 - 2]", 0.0,
                      min = 0, max = 2, step = 0.1),
             actionButton("updateMap","Update Landscape")
           )
    ),
    column(4,
           inputPanel(
             h4("Biodiversity parameter"),
             numericInput("nInd","Total number of individuals [100 - 10000]",
                          min = 100, max = 10000, value = 1000, step = 500),
             numericInput("nSpec","Number of species [1 - 100]",
                          min = 1, max = 100, value = 10, step = 5),
             numericInput("sigma","Cluster size [1 - 1000]",
                          min = 1, max = 1000, value = 20, step = 10),
             actionButton("simSpecies","Simulate species distributions")
           )
    )
  ), # end row1
  
  fluidRow(
    column(4,
           wellPanel(
             h4("Distribution of habitat and matrix"),
             plotOutput("landscape"))),
    column(4,
           wellPanel(
           h4("Species distributions in continuous habitat"),
           plotOutput("community"))),
    column(4,
           wellPanel(
           h4("Species distributions in fragmented landscape"),
           plotOutput("com_frag")))
  ),
  
  fluidRow(
    column(4,
           wellPanel(
             h4("Landscape indices"),
             tableOutput("metrics"))),
    column(4,
           wellPanel(
             h4("Diversity indices in continuous habitat"),
             tableOutput("tab_div_cont"))),
    column(4,
           wellPanel(
             h4("Diversity indices in fragmented landscape"),
             tableOutput("tab_div_frag")))
  ),
  
  tags$h5("More information on geometric fragmentation effects:",
          tags$a(href = "https://onlinelibrary.wiley.com/doi/10.1002/ece3.4951","May et al. (2019) Ecology & Evolution")),
  
  tags$h5("Source code of this App in", tags$strong("R"),":",
          tags$a(href = "https://github.com/FelixMay/GeomFragApp","GitHub"))
)

server <- function(input, output) {
  
  map <- eventReactive(input$updateMap, {
    map1 <- nlm_mpd(ncol = 100,
                    nrow = 100,
                    roughness = input$hab_rough,
                    verbose = F,
                    resolution = 10)
    if (input$hab_amount < 0.1){
      map1_bin <- map1
      map1_bin[] <- 0
    } else if (input$hab_amount > 99.9){
      map1_bin <- map1
      map1_bin[] <- 1
    } else {
      map1_bin <- map1
      map1_bin[] <- 0
      map1_bin[map1 > quantile(map1[], probs = 1.0 - input$hab_amount/100)] <- 1
      #map1_bin <- util_binarize(map1, breaks = input$hab_amount/100)
    }
    projection(map1_bin) <- "+init=epsg:32631"
    # check_landscape(map1_bin)
    map1_bin
  })

  metrics <- reactive({
    class_metrics <- dplyr::bind_rows(
      lsm_c_np(map()),
      lsm_c_area_mn(map()), # hectares
      lsm_c_area_sd(map()), # hectares
      lsm_c_enn_mn(map()),  # meters
      lsm_c_enn_sd(map())   # meters
    ) %>% dplyr::filter(class == 1)

    table1 <- data.frame(Index = c("No. of fragments",
                                               "Fragment area (mean)",
                                               "Fragment area (standard dev.)",
                                               "Distance to nearest neighbor (mean)",
                                               "Distance to nearest neighbor (standard dev.)"
                                               ),
                         Value = class_metrics$value,
                         Unit = c("","ha","ha","m","m")
        )
    table1
  })

  output$landscape <- renderPlot({
    plot_map(map())
  }) # end output landscape

  output$metrics <- renderTable(metrics())

  com1 <- eventReactive(input$simSpecies,{
    sim1 <- sim_thomas_community(input$nSpec, input$nInd, sigma = input$sigma, fix_s_sim = F,
                                 xrange = c(0,1000), yrange = c(0,1000), mother_points = 1)
    points1 <- sim1$census
    coordinates(points1) = ~x + y
    projection(points1) <- "+init=epsg:32631"
    grid1 <- as(map(), "SpatialGridDataFrame")
    #grid1 <- as(map1_bin, "SpatialGridDataFrame")
    points1$Class <- as.factor(over(points1, grid1)[[1]])
    points1
  })
  
  div_cont <- reactive({
    # abund_vec <- table(points1$species)
    abund_vec <- table(com1()$species)
    abund_vec <- abund_vec[abund_vec > 0]
    N <- sum(abund_vec)
    S <- length(abund_vec)
    Shannon <- diversity(abund_vec, index = "shannon")
    J <- Shannon/log(S)
    data.frame(Index = c("No. of individuals","Species richness","Shannon-Diversity","Evenness"),
               Wert = c(N, S, Shannon, J))
  })
  
  output$tab_div_cont <- renderTable(div_cont())
  
  div_frag <- reactive({
    #points2 <- points1$[points1$Class == 1,]
    points2 <- com1()[com1()$Class == 1,]
    abund_vec <- table(points2$species)
    abund_vec <- abund_vec[abund_vec > 0]
    N <- sum(abund_vec)
    S <- length(abund_vec)
    Shannon <- diversity(abund_vec, index = "shannon")
    J <- Shannon/log(S)
    data.frame(Index = c("No. of individuals","Species richness","Shannon-Diversity","Evenness"),
               Value = c(N, S, Shannon, J))
  })
  
  output$tab_div_frag <- renderTable(div_frag())

  output$community <- renderPlot({
    rb_cols <- rainbow(length(unique(com1()$species)))
    map2 <- map()
    map2[] <- 1
    plot_map(map2)
    points(y ~ x, data = com1(), pch = 19, col = rb_cols[com1()$species])
  })
  
  output$com_frag <- renderPlot({
    plot_map(map())
    rb_cols <- rainbow(length(unique(com1()$species)))
    points2 <- com1()[com1()$Class == 1,]
    points(y ~ x, data = points2, pch = 19, col = rb_cols[points2$species])
  })

}

shinyApp(ui = ui, server = server)