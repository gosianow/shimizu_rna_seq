require(graphics)
data(crabs, package = "MASS")

lcrabs <- log(crabs[, 4:8])
crabs.grp <- factor(c("B", "b", "O", "o")[rep(1:4, rep(50,4))])
gr <- somgrid(topo = "hexagonal")
crabs.som <- SOM(lcrabs, gr)
plot(crabs.som)

## 2-phase training
crabs.som2 <- SOM(lcrabs, gr,
                  alpha = list(seq(0.05, 0, len = 1e4), seq(0.02, 0, len = 1e5)),
                  radii = list(seq(8, 1, len = 1e4), seq(4, 1, len = 1e5)))
plot(crabs.som2)


####################################################


# Load the kohonen package 
require(kohonen)

# Create a training data set (rows are samples, columns are variables
# Here I am selecting a subset of my variables available in "data"
# data_train <- data[, c(2,4,5,8)]

data_train <- crabs
data_train$sp <- factor(data_train$sp)
levels(data_train$sp) <- c(1, 2)
data_train$sp <- as.numeric(data_train$sp)
data_train$sex <- factor(data_train$sex)
levels(data_train$sex) <- c(1, 2)
data_train$sex <- as.numeric(data_train$sex)


# Change the data frame with training data to a matrix
# Also center and scale all variables to give them equal importance during
# the SOM training process. 
data_train_matrix <- scale(as.matrix(data_train))

# Create the SOM Grid - you generally have to specify the size of the 
# training grid prior to training the SOM. Hexagonal and Circular 
# topologies are possible
som_grid <- somgrid(xdim = 4, ydim = 3, topo="hexagonal")

# Finally, train the SOM, options for the number of iterations,
# the learning rates, and the neighbourhood are available
som_model <- som(data_train_matrix, 
                 grid=som_grid)


som_model <- som(data_train_matrix, 
                 grid=som_grid, 
                 rlen=100, 
                 alpha=c(0.05,0.01), 
                 keep.data = TRUE,
                 n.hood="circular")



plot(som_model, type="changes")

coolBlueHotRed <- function(n, alpha = 1) {
  rainbow(n, end=4/6, alpha=alpha)[n:1]
}

library(RColorBrewer)

colgrey <- function(n) colorRampPalette(c("white", "black"))(n)


plot(som_model, type="counts")

plot(som_model, type="dist.neighbours", palette.name = colgrey)

plot(som_model, type="codes")

var <- 3 
plot(som_model, type = "property", property = som_model$codes[,var], main=colnames(som_model$data)[var], palette.name=coolBlueHotRed)


var <- 2 #define the variable to plot 
var_unscaled <- aggregate(as.numeric(data_train[,var]), by=list(som_model$unit.classif), FUN=mean, simplify=TRUE)[,2] 
plot(som_model, type = "property", property=var_unscaled, main=names(data_train)[var], palette.name=coolBlueHotRed)


mydata <- som_model$codes 
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var)) 
for (i in 2:6) {
  wss[i] <- sum(kmeans(mydata, centers=i)$withinss)
}
plot(wss)



## use hierarchical clustering to cluster the codebook vectors
som_cluster <- cutree(hclust(dist(som_model$codes)), 3)
# plot these results:
plot(som_model, type="mapping", bgcol = coolBlueHotRed(3)[som_cluster] , main = "Clusters") 
add.cluster.boundaries(som_model, som_cluster)




#################################################### paraller coordianates

# Load required packages
require(GGally)

# Load datasets
data(state)
df <- data.frame(state.x77,
                 State = state.name,
                 Abbrev = state.abb,
                 Region = state.region,
                 Division = state.division
) 

# Generate basic parallel coordinate plot
p <- ggparcoord(data = df,                 
                # Which columns to use in the plot
                columns = 1:4,                 
                # Which column to use for coloring data
                groupColumn = 11,                 
                # Allows order of vertical bars to be modified
                order = "anyClass",                
                # Do not show points
                showPoints = FALSE,                
                # Turn on alpha blending for dense plots
                alphaLines = 0.6,                
                # Turn off box shading range
                shadeBox = NULL,                
                # Will normalize each column's values to [0, 1]
                scale = "uniminmax" # try "std" also
)

# Start with a basic theme
p <- p + theme_minimal()

# Decrease amount of margin around x, y values
p <- p + scale_y_continuous(expand = c(0.02, 0.02))
p <- p + scale_x_discrete(expand = c(0.02, 0.02))

# Remove axis ticks and labels
p <- p + theme(axis.ticks = element_blank())
p <- p + theme(axis.title = element_blank())
p <- p + theme(axis.text.y = element_blank())

# Clear axis lines
p <- p + theme(panel.grid.minor = element_blank())
p <- p + theme(panel.grid.major.y = element_blank())

# Darken vertical lines
p <- p + theme(panel.grid.major.x = element_line(color = "#bbbbbb"))

# Move label to bottom
p <- p + theme(legend.position = "bottom")

# Figure out y-axis range after GGally scales the data
min_y <- min(p$data$value)
max_y <- max(p$data$value)
pad_y <- (max_y - min_y) * 0.1

# Calculate label positions for each veritcal bar
lab_x <- rep(1:4, times = 2) # 2 times, 1 for min 1 for max
lab_y <- rep(c(min_y - pad_y, max_y + pad_y), each = 4)

# Get min and max values from original dataset
lab_z <- c(sapply(df[, 1:4], min), sapply(df[, 1:4], max))

# Convert to character for use as labels
lab_z <- as.character(lab_z)

# Add labels to plot
p <- p + annotate("text", x = lab_x, y = lab_y, label = lab_z, size = 3)

# Display parallel coordinate plot
print(p)


##################################### spline

require(graphics)

attach(cars)
plot(speed, dist, main = "data(cars)  &  smoothing splines")
cars.spl <- smooth.spline(speed, dist)
(cars.spl)
## This example has duplicate points, so avoid cv = TRUE

lines(cars.spl, col = "blue")
lines(smooth.spline(speed, dist, df = 10), lty = 2, col = "red")
legend(5,120,c(paste("default [C.V.] => df =",round(cars.spl$df,1)),
               "s( * , df = 10)"), col = c("blue","red"), lty = 1:2,
       bg = 'bisque')
detach()


## Residual (Tukey Anscombe) plot:
plot(residuals(cars.spl) ~ fitted(cars.spl))
abline(h = 0, col = "gray")

## consistency check:
stopifnot(all.equal(cars$dist,
                    fitted(cars.spl) + residuals(cars.spl)))

## Visualize the behavior of  .nknots.smspl()
nKnots <- Vectorize(.nknots.smspl) ; c.. <- adjustcolor("gray20",.5)
curve(nKnots, 1, 250, n=250)
abline(0,1, lty=2, col=c..); text(90,90,"y = x", col=c.., adj=-.25)
abline(h=100,lty=2); abline(v=200, lty=2)

n <- c(1:799, seq(800, 3490, by=10), seq(3500, 10000, by = 50))
plot(n, nKnots(n), type="l", main = "Vectorize(.nknots.smspl) (n)")
abline(0,1, lty=2, col=c..); text(180,180,"y = x", col=c..)
n0 <- c(50, 200, 800, 3200); c0 <- adjustcolor("blue3", .5)
lines(n0, nKnots(n0), type="h", col=c0)
axis(1, at=n0, line=-2, col.ticks=c0, col=NA, col.axis=c0)
axis(4, at=.nknots.smspl(10000), line=-.5, col=c..,col.axis=c.., las=1)

##-- artificial example
y18 <- c(1:3, 5, 4, 7:3, 2*(2:5), rep(10, 4))
xx  <- seq(1, length(y18), len = 201)
(s2  <- smooth.spline(y18)) # GCV
(s02  <- smooth.spline(y18, spar = 0.2))
(s02. <- smooth.spline(y18, spar = 0.2, cv = NA))
plot(y18, main = deparse(s2$call), col.main = 2)
lines(s2, col = "gray"); lines(predict(s2, xx), col = 2)
lines(predict(s02, xx), col = 3); mtext(deparse(s02$call), col = 3)



## The following shows the problematic behavior of 'spar' searching:
(s2  <- smooth.spline(y18, control =
                      list(trace = TRUE, tol = 1e-6, low = -1.5)))
(s2m <- smooth.spline(y18, cv = TRUE, control =
                      list(trace = TRUE, tol = 1e-6, low = -1.5)))
## both above do quite similarly (Df = 8.5 +- 0.2)


##################################### spline 2

## Linear Model
fit = lm(y ~ x )

## Polynomial
fit.3 = lm(y ~ poly(x,3) )

## Natural Spline
fit.ns.3 = lm(y ~ ns(x, 3) )
## Smoothing Spline
fit.sp = smooth.spline(y ~ x, nknots=15)



























