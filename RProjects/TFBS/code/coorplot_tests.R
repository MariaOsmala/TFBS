library("corrplot")
M = cor(mtcars)
corrplot(M, method = 'number') # colorful number

corrplot(M, method = 'color', order="alphabet") # colorful number

corrplot(M)
corrplot(M, order="AOE")

corrplot(M, method = 'shade', order = 'AOE', diag = FALSE)
corrplot(M, method = 'square', order = 'FPC', type = 'lower', diag = FALSE)

corrplot(M, method = 'ellipse', order = 'AOE', type = 'upper')

corrplot.mixed(M, order = 'AOE')
corrplot.mixed(M, lower = 'shade', upper = 'pie', order = 'hclust')

corrplot(M, order = 'AOE', 
         addCoef.col = 'black', #adds text
         tl.pos = 'd', #adds variable names to the diagomal
         cl.pos = 'n', #show color bar or not
         col = COL2('PiYG'))


N1 = matrix(runif(80, 0, 1), 8)
corrplot(N1, is.corr = FALSE, col.lim = c(0, 1), method = 'color', tl.pos = 'n',
         col = rev(COL1('YlOrRd')), cl.pos = 'b', addgrid.col = 'white', addCoef.col = 'grey50')


corrplot(N1, is.corr = FALSE, col.lim = c(0, 1), method = 'color',
#         tl.pos = 'n',
 col = rev(COL1('YlOrRd')), 
#cl.pos = 'b', #color bar position
        #addgrid.col = 'white', 
        addCoef.col = 'grey50',
        )


testRes = cor.mtest(mtcars, conf.level = 0.95)

## specialized the insignificant value according to the significant level
corrplot(M, p.mat = testRes$p, sig.level = 0.10, order = 'hclust', addrect = 2)

## leave blank on non-significant coefficient
## add significant correlation coefficients
corrplot(M, p.mat = testRes$p, method = 'circle', type = 'lower', insig='blank',
         addCoef.col ='black', number.cex = 0.8, order = 'AOE', diag=FALSE)

## leave blank on non-significant coefficient
## add all correlation coefficients
corrplot(M, p.mat = testRes$p, method = 'circle', type = 'lower', insig='blank',
         order = 'AOE', diag = FALSE)$corrPos -> p1
text(p1$x, p1$y, round(p1$corr, 2))

corrplot(M, p.mat = testRes$p, insig = 'p-value')

## add all p-values
corrplot(M, p.mat = testRes$p, insig = 'p-value', sig.level = -1)

## add significant level stars
corrplot(M, p.mat = testRes$p, method = 'color', diag = FALSE, type = 'upper',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'grey20', order = 'AOE')


set.seed(123)
mat <- matrix(rnorm(36), ncol=6)
rownames(mat) <- colnames(mat) <- as.character(1:6)
# Convert matrix to long format dataframe
df <- as.data.frame(as.table(mat))
df$Var1 <- as.numeric(as.character(df$Var1))
df$Var2 <- as.numeric(as.character(df$Var2))


# Filter out the upper triangle and diagonal
df <- subset(df, Var1 < Var2)

# Plot the heatmap
ggplot(df, aes(x=Var2, y=Var1)) +
  geom_tile(aes(fill=Freq), color="white") +
  scale_fill_viridis_c() + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


library(ggplot2)
library(grid)

# Sample ggplot object
data(mtcars)
p <- ggplot(mtcars, aes(mpg, disp)) + geom_point()

# Convert ggplot to a gList
gList_obj <- ggplotGrob(p)

# Sample grob (rectangle)
rect_grob <- rectGrob(width = 0.5, height = 0.5, gp = gpar(fill = "blue"))

# Sample gTree (circle)
circle_gtree <- gTree(children = gList(circleGrob( gp = gpar(fill = "red"))))

# Create viewports
vp1 <- viewport(width = 0.3, height = 0.3, x = 0.25, y = 0.75, name="VP1")
vp2 <- viewport(width = 0.3, height = 0.3, x = 0.75, y = 0.75, name="VP2")
vp3 <- viewport(width = 0.3, height = 0.3, x = 0.5, y = 0.25, name="VP3")

# Plot using the grid system
grid.newpage()
pushViewport(vp1)
grid.draw(gList_obj)
popViewport()
pushViewport(vp2)
grid.draw(rect_grob)
popViewport()
pushViewport(vp3)
grid.draw(circle_gtree)
