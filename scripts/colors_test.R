#print colored text of each colors() matching key
plot0()
par(mai = rep(0,4))
key = "blue"
colors_key = colors()[grepl(key, colors())]
usr = par('usr')
x_min = usr[1]; x_max = usr[2]; y_min = usr[3]; y_max = usr[4]
ncol = 5
nrow = 20
x_step = (x_max - x_min) / ncol
y_step = (y_max - y_min) / nrow
x = x_min + x_step
y = y_max - y_step
for(c in colors_key){
  text(x,y, col = c, labels = c)
  y = y - y_step
  if(y < y_min){
    y = y_max - y_step
    x = x + x_step
  }
  if(x > x_max){
    print(paste("stopped after printing", c))
   break 
  }
}