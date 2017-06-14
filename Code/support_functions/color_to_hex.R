color_to_hex <- function(colorname, transparency = 0) {
  # transparency should be between 0 and 1
  
  color_256 = col2rgb(colorname, alpha = T)
  color_256[[4]] = round((1 - transparency) * 255)
  
  color_hex = paste(as.hexmode(color_256), collapse = "")
  
  # add octothorpe and capitalize letters
  color_hex = paste("#", toupper(color_hex), sep = "")
  
  return(color_hex)
}