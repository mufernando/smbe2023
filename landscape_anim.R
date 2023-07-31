library(ggplot2)
library(gganimate)
library(gifski)


norm = function(x) return((x-mean(x))/sd(x))

sim_ar = function(n=10000, a1=0.18828, a2=0.05861, sd=2, normalize=F) {
  
  # generate n+2 standard normal variates
  
  E = rnorm(n+2)
  
  # create an autoregressive process and plot the first 200 observations,
  # the autocorrelation function, and the partial autocorrelation function
  
  Y = numeric(n)
  Y[1] = E[3] + a1*E[2] + a2*E[1]# + rnorm(1,sd=sd)
  Y[2] = E[4] + a1*Y[1] + a2*E[2]# + rnorm(1, sd=sd)
  for (i in 3:n) Y[i] = E[i+2] + a1*Y[i-1] + a2*Y[i-2]# + rnorm(1, sd=sd)
  if(normalize) {
    Y = norm(Y)
  }
  return (norm(Y))
}

n=500
set.seed(93384)



produce_df = function(n,gene_effect=0., invert=1, normalize=T) {
y1 = sim_ar(n, a1=0.3, a2=0.2, normalize=normalize)
y2 = sim_ar(n, a1=0.15, a2=0.1, normalize=normalize)
df = data.frame(win=rep(1:n, 2), y = c(y1,y1+rnorm(n, sd=0.5)), spp=factor(c(rep("v",n), rep("w",n))), time=0)
df2 = df
new_y1 = sim_ar(n, a1=0.3, a2=0.2)
new_y2 = sim_ar(n, a1=0.15, a2=0.1)
new_y1[100:200] = new_y1[100:200] - gene_effect
new_y1[300:400] = new_y1[300:400] - gene_effect
if (gene_effect > 0) {
  new_y2[100:200] = new_y1[100:200] + rnorm(n, sd=0.5)
  new_y2[300:400] = new_y1[300:400] + rnorm(n, sd=0.5)
}
df2$y = c(norm(new_y1),norm(new_y2))
df2$time = 4
df = rbind(df, df2)

return(df)
}


df = produce_df(n)

cpal = c("#BB4B49", "#6069c0")
anim = ggplot(df, aes(y=y, x=win, col=spp)) +
  geom_smooth(method="loess", span=0.1, se=F) + 
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        text=element_text(size=38), #change font size of all text
        axis.text=element_text(size=38), #change font size of axis text
        axis.title=element_text(size=38), #change font size of axis titles
        plot.title=element_text(size=38), #change font size of plot title
        legend.text=element_text(size=38), #change font size of legend text
        legend.title=element_text(size=38)) + #change font size of legend title   ) +
  labs(x="Chromosome position", y="Diversity (\U03C0)", title = 'Time: {round(frame_time,2)}N generations', col="Species", ) +
  scale_colour_manual(values=cpal) +
  transition_time(time)

animate(
  anim + enter_fade() + exit_fly(y_loc = 1), renderer = gifski_renderer(), fps=7, res=150, width="300", height="200", units="mm"
)
anim_save("~/projects/smbe23/two_spp_neutral_div_landscapes.gif")

set.seed(19928)

df = produce_df(n, gene_effect=1, normalize=T)

#cpal = c("#80b573", "#a8447d")
anim = ggplot(df, aes(y=y, x=win, col=spp)) +
  geom_smooth(method="loess", span=0.1, se=F) + 
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        text=element_text(size=38), #change font size of all text
        axis.text=element_text(size=38), #change font size of axis text
        axis.title=element_text(size=38), #change font size of axis titles
        plot.title=element_text(size=38), #change font size of plot title
        legend.text=element_text(size=38), #change font size of legend text
        legend.title=element_text(size=38)) + #change font size of legend title   ) +
  labs(x="Chromosome position", y="Diversity (\U03C0)", title = 'Time: {round(frame_time,2)}N generations', col="Species", ) +
  scale_colour_manual(values=cpal) +
  transition_time(time)

animate(
  anim + enter_fade() + exit_fly(y_loc = 1), renderer = gifski_renderer(), fps=7, res=150, width="300", height="200", units="mm"
)
anim_save("~/projects/smbe23/two_spp_sel_div_landscapes.gif")

