"""
Generate a 3x3 heatmap where each square is the average (over N runs) # of live cancer cells.
The x-axis has the 3 attachment rate values (0.033, 0.2, 1) and 
the y-axis has the 3 attachment lifetime values (15, 60, 120). 

Each of the 3x3 squares is annotated with N=# of runs (random seeds) that completed. 


Input argument:
  fixed migration bias (use the middle value = 0.5 for the heatmap in the BMC paper).

Output
  a .png at 300 dpi.

Module dependencies: numpy, matplotlib

Author:
Randy Heiland 
Intelligent Systems Engineering, Indiana University

"""

import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Rectangle
import numpy as np

mpl.rcParams['savefig.dpi'] = 300

print("len(sys.argv)=", len(sys.argv))
if len(sys.argv) < 2:
  print("Usage: %s <migration_bias: 0.25, 0.5, or 0.75>" % sys.argv[0])
  sys.exit(0)

migration_bias = float(sys.argv[1])
print("migration_bias=", migration_bias)

# -------- read in metrics file:
v = np.loadtxt("final_analysis.txt", delimiter=',')
run_num = v[0]
num_live_cancer_cells = v[1]
num_live_cancer_cells_above_threshold = v[2]
mean_oncoprotein = v[3]


# -------- read in params file:
#  -- on each line, we have 4 values: run#, rate, lifetime, bias
#     where:
# Attachment_rate = [0.033, 0.2, 1] - mean wait time, hours
# Attachment_lifetime = [15, 60, 120] - mins
# Migration_bias = [0.25, 0.5, 0.75] - random migration bias towards chemotactic gradients

run, rate, lifetime, bias = np.genfromtxt('run_rate_lifetime_bias.txt', delimiter=',').transpose()


# compute average # of <metric values>, e.g. # live cancer cells, mean oncoprotein, 
m1 = np.zeros([3,3])
N_val = np.zeros([3,3])

for dummy in [1]:  # each bias value generates entire set of 3 heatmaps (for 3 metrics)
  iy = 0
  for lval in [15, 60, 120]:  # from bottom-to-top (Y-axis)
    ix = 0
    for rval in [0.033, 0.2, 1]:  # from left-to-rightop (X-axis)
      idx = np.where( (bias==migration_bias) & (lifetime==lval) & (rate==rval) )
      metric_vals = num_live_cancer_cells[idx]   # NB! change RHS for other metrics
      N = len(metric_vals)  # length of array = # of (random seed) runs successful
      N_val[iy][ix] = N
      avg_metric_vals = metric_vals.sum() / N
      print("N, avg_metric_vals = ",N, avg_metric_vals)
      m1[iy][ix] = avg_metric_vals
      ix += 1
    iy += 1


#-------------------------------------------------
# plot coarse-grained (3x3) mesh, color-mapped for "Avg(<metric>) over N runs" (N=max 10 = random seed range)

# make these smaller to increase the resolution
dx, dy = 1.0, 1.0

# generate 2 2d grids for the x & y bounds
y, x = np.mgrid[slice(0, 3 + dy, dy),
                slice(0, 3 + dx, dx)]

# x and y are bounds, so m1 is the value *inside* those bounds.
levels = MaxNLocator(nbins=128).tick_values(m1.min(), m1.max())


# pick the desired colormap, sensible levels, and define a normalization
# instance which takes data values and translates those into levels.
#cmap = plt.get_cmap('PiYG')
cmap = plt.get_cmap('plasma')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

fig, ax = plt.subplots(nrows=1)

im = ax.pcolormesh(x, y, m1, cmap=cmap, norm=norm)
fig.colorbar(im, ax=ax)
title_string = "Avg(# live cancer cells) per N runs; Migration Bias=%s" % (migration_bias)
ax.set_title(title_string)
ax.set_xlabel("attachment rate")
ax.set_ylabel("attachment lifetime")

#ax.set_xticks([])
#ax.set_yticks([])
plt.xticks([0.5, 1.5, 2.5], ["0.033", "0.2", "1"])   # attachment rate
plt.yticks([0.5, 1.5, 2.5], ["15", "60", "120"])   # attachment lifetime

#-----------
t = ax.transData
canvas = ax.figure.canvas
for iy in range(3):
  for ix in range(3):
    sval = "N=%d" % N_val[iy][ix]
    text = ax.text(ix+0.4, iy+0.1, sval, backgroundcolor="white", color="black", transform=t)
    text.set_bbox(dict(facecolor='white'))
#    text.set_bbox(dict(boxstyle='round,pad=1'))
    text.draw(canvas.get_renderer())

#text = ax.text(1.4, 0.1, "N=9", color="white", transform=t)
#text.draw(canvas.get_renderer())
#text = ax.text(2.4, 0.1, "N=8", color="white", transform=t)
#text.draw(canvas.get_renderer())


# adjust spacing between subplots so `ax1` title and `ax0` tick labels
# don't overlap
fig.tight_layout()

png_filename = "heatmap_avg_num_live_cancer_cells_bias%s.png" % (migration_bias)
print("---> ", png_filename)
fig.savefig(png_filename, bbox_inches='tight')

#plt.show()
