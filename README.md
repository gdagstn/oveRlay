<img src="https://user-images.githubusercontent.com/21171362/169328724-cf55b344-bee2-44c7-b4a8-eddac054dff8.png" align="right" alt="" width="200" />

# oveRlay
Annotate 2D point clouds using overlay polygons

# Motivation and methodology

This R package (in development) offers an alternative to the concave hull geom in `ggforce` by taking an entirely different approach:
- creaties a regularly spaced grid of squares, whose side is defined by a step size (a fraction of the range of your data)
- calculates which squares contain points, with user-defined tolerance, and keeps the vertices
- uses `isoband` to calculate a contour joining vertices on the grid
- calculates holes, i.e. which polygons are inside other polygons, so that their space can be "subtracted" when plotting by `ggplot2` 
- optional smoothing (via `smoothr`), offset (via `polyclip`) and joining of overlapping polygons (again via `polyclip`)

The resulting polygon (or sets of polygons) follow the shapes of the point cloud quite closely and, if a point or set of points is far enough, another disjointed polygon is created. This allows to deal with outliers nicely; the current `ggforce` implementation of `geom_mark_hull()` is super fast, but will enclose *all* points no matter how far they are. Moreover, and to the best of my knowledge, the concave hull in `ggforce` does not allow holes. 

The step size controls how granular the resulting polygon set is: small step sizes will create shapes that follow the data more closely (at the expense of computing speed and, eventually, usability). 

Here is an example of how different step sizes behave on the same set of 2D `rnorm(1000)`:

<img width="976" alt="Screenshot 2022-05-19 at 8 05 26 PM" src="https://user-images.githubusercontent.com/21171362/169289332-3b4d2c5e-5160-4ab6-8847-85318c1062c4.png">

# Installation

You can install `oveRlay` through `devtools`:

```{r}
library(devtools)
install_github("gdagstn/oveRlay")
```

# Examples

Using oveRlay is simple:


```{r}
library(oveRlay)

dat <- matrix(rnorm(1000), ncol = 2)

overlay <- makeOverlay(dat, min_pts = 1, stepsize = 0.06, minsize = 4)
plot(dat, pch = 16, cex = 0.5, xlim = range(overlay[,1:2]), ylim = range(overlay[,1:2]))
for(i in unique(overlay$cluster)) polygon(overlay[overlay$cluster == i, 1:2])
```
<img width="635" alt="overlay polygon" src="https://user-images.githubusercontent.com/21171362/169289869-45c337cd-7b97-4adf-bd76-dfa6bc52ada8.png">

Decreasing step size (increasing granularity)

```
overlay <- makeOverlay(dat, min_pts = 1, stepsize = 0.02, minsize = 4)
plot(dat, pch = 16, cex = 0.5, xlim = range(overlay[,1:2]), ylim = range(overlay[,1:2]))
for(i in unique(overlay$cluster)) polygon(overlay[overlay$cluster == i, 1:2])
```
<img width="650" alt="overlay with decreased step size" src="https://user-images.githubusercontent.com/21171362/169290024-3c730088-9deb-43d8-8ef7-ba4522eccef5.png">

Increasing offset

```
overlay <- makeOverlay(dat,min_pts = 1, stepsize = 0.02, minsize = 4, offset_prop = 0.08)
plot(dat, pch = 16, cex = 0.5, xlim = range(overlay[,1:2]), ylim = range(overlay[,1:2]))
for(i in unique(overlay$cluster)) polygon(overlay[overlay$cluster == i, 1:2])
```
<img width="652" alt="overlay with increased offset" src="https://user-images.githubusercontent.com/21171362/169290134-d2d68930-4f55-47af-889d-baa0af641348.png">

Increasing offset without joining polygons:

```
overlay <- makeOverlay(dat,min_pts = 1, stepsize = 0.02, minsize = 4, offset_prop = 0.08, join_polys = FALSE)
plot(dat, pch = 16, cex = 0.5, xlim = range(overlay[,1:2]), ylim = range(overlay[,1:2]))
for(i in unique(overlay$cluster)) polygon(overlay[overlay$cluster == i, 1:2])
```

<img width="641" alt="overlay without joining" src="https://user-images.githubusercontent.com/21171362/169290198-e2213731-fe43-463e-9241-4e9cf7771025.png">

However, the output of `makeOverlay()` is mostly designed to be used with `ggplot2`. When using it as an input to `geom_polygon()` it is important to use `group = "cluster_hole"` and `subgroup = "id_hole"` in a call to `aes_string()` (or `aes` without quotes). The use of `group` and `subgroup` allows drawing the polygon geom with holes. 

```
library(ggplot2)

dat =  as.data.frame(matrix(rnorm(1000), ncol = 2))
colnames(dat) = c("x", "y")

# Make a hole in the distribution
dat = dat[abs(dat$x) > 0.5 | abs(dat$y) > 0.5,]

overlay <- makeOverlay(dat, min_pts = 1, stepsize = 0.02, minsize = 4, offset_prop = 0.01, join_polys = TRUE)

ggplot(data = dat, aes(x = x, y = y)) + 
  geom_point() + 
  geom_polygon(data = overlay,
               aes_string(x = "x", y = "y", group = "cluster_hole", subgroup = "id_hole"),
               color = "red", fill = "red", alpha = 0.3) + 
  theme_bw()
```

<img width="656" alt="Screenshot 2022-05-19 at 10 46 00 PM" src="https://user-images.githubusercontent.com/21171362/169324584-f7fefb09-d34a-4cc3-a992-b78086a90f0a.png">

# To do

- Handling more than one group of points at the same time 
- Handling overlapping regions with different categories
- Integration with UMAP/tSNE visualizations in single cell packages
- Recentering of single-point polygons

# Acknowledgments

- Claus Wilke and Thomas Lin Pedersen for the [isoband](https://wilkelab.org/isoband/) package
- Thomas Lin Pedersen for the [ggforce](https://ggforce.data-imaginist.com/) package
- Adrian Baddeley for the [polyclip](https://github.com/baddstats/polyclip) package, and Angus Johnson for the original [Clipper](http://angusj.com/delphi/clipper.php) library
- Matthew Strimas Mackey for the [smoothr](https://strimas.com/smoothr/) package
- Hadley Wickham, Winston Chang, Lionel Henry, Thomas Lin Pedesen, Kohske Takahashi, Claus Wilke, Kara Woo, and Hiroaki Yutani for the [ggplot2](https://ggplot2.tidyverse.org/) package
