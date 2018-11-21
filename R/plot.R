



plot_link_osm<-function(osmid){
    require(utils)
    a  = paste0('https://www.openstreetmap.org/way/', osmid)
    browseURL(a)

}
