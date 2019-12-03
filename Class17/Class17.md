Biological Network Analysis
================

Examine data from the Tara Ocean project on bacteria in ocean water samples from around the world.

Set up Cytoscape and R connection
=================================

We will use the **igraph** package from CRAN and the **RCy3** pacakage from bioconductor to build and visualize networks of these data.

``` r
#Load packages.
library(RCy3)
```

    ## Warning: package 'RCy3' was built under R version 3.6.1

``` r
library(igraph)
```

    ## 
    ## Attaching package: 'igraph'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     decompose, spectrum

    ## The following object is masked from 'package:base':
    ## 
    ##     union

``` r
library(RColorBrewer)
library(BBmisc)
```

    ## 
    ## Attaching package: 'BBmisc'

    ## The following object is masked from 'package:igraph':
    ## 
    ##     normalize

    ## The following object is masked from 'package:base':
    ## 
    ##     isFALSE

``` r
#Test the connection to Cytoscape.
cytoscapePing()
```

    ## [1] "You are connected to Cytoscape!"

We can test things further by making a small network (here in igraph format as used by the R igraph package) and sending it to Cytoscape.

``` r
g <- makeSimpleIgraph()
createNetworkFromIgraph(g,"myGraph")
```

    ## Loading data...
    ## Applying default style...
    ## Applying preferred layout...

    ## networkSUID 
    ##          52

We can include this Cytoscape rendered network image in our report. Cytoscape provides a number of canned visual styles. The code below explores some of these styles. For example check out the **marquee** style!

``` r
setVisualStyle("Marquee")
```

    ##                 message 
    ## "Visual Style applied."

``` r
fig <- exportImage(filename="demo_marquee", type="png", height=350)
```

    ## Warning: This file already exists. A Cytoscape popup 
    ##                 will be generated to confirm overwrite.

``` r
knitr::include_graphics("./demo_marquee.png")
```

<img src="./demo_marquee.png" width="640" />

Read our metagenomics data
==========================

We will read in a species co-occurrence matrix that was calculated using Spearman Rank coefficient.

``` r
#Scripts for processing located in "inst/data-raw/".
prok_vir_cor <- read.delim("virus_prok_cor_abundant.tsv", stringsAsFactors = FALSE)

#Have a peak at the first 6 rows.
head(prok_vir_cor)
```

    ##       Var1          Var2    weight
    ## 1  ph_1061 AACY020068177 0.8555342
    ## 2  ph_1258 AACY020207233 0.8055750
    ## 3  ph_3164 AACY020207233 0.8122517
    ## 4  ph_1033 AACY020255495 0.8487498
    ## 5 ph_10996 AACY020255495 0.8734617
    ## 6 ph_11038 AACY020255495 0.8740782

Here we will use the igraph package to convert the co-occurrence dataframe into a network that we can send to Cytoscape. In this case our graph is undirected (so we will set `directed = FALSE`) since we do not have any information about the direction of the interactions from this type of data.

``` r
g1 <- graph.data.frame(prok_vir_cor, directed = FALSE)
class(g1)
```

    ## [1] "igraph"

``` r
g1
```

    ## IGRAPH bb42d27 UNW- 845 1544 -- 
    ## + attr: name (v/c), weight (e/n)
    ## + edges from bb42d27 (vertex names):
    ##  [1] ph_1061 --AACY020068177 ph_1258 --AACY020207233
    ##  [3] ph_3164 --AACY020207233 ph_1033 --AACY020255495
    ##  [5] ph_10996--AACY020255495 ph_11038--AACY020255495
    ##  [7] ph_11040--AACY020255495 ph_11048--AACY020255495
    ##  [9] ph_11096--AACY020255495 ph_1113 --AACY020255495
    ## [11] ph_1208 --AACY020255495 ph_13207--AACY020255495
    ## [13] ph_1346 --AACY020255495 ph_14679--AACY020255495
    ## [15] ph_1572 --AACY020255495 ph_16045--AACY020255495
    ## + ... omitted several edges

In this case the first line of output (“UNW- 854 1544 –”) tells that our network graph has 845 vertices (i.e. nodes, which represent our bacteria and viruses) and 1544 edges (i.e. linking lines, which indicate their co-occurrence). Note that the first four characters (i.e. the “UNW-” part) tell us about the network setup. In this case our network is Undirected, Named (i.e. has the ‘name’ node/vertex attribute set) and Weighted (i.e. the ‘weight’ edge attribute is set).

``` r
#Plotting g1 alone is a hot mess. Vertex labels are automatically shown and the default vertex size is 15.
plot(g1, vertex.size=3, vertex.label=NA)
```

![](Class17_files/figure-markdown_github/unnamed-chunk-8-1.png)

Network querys
--------------

Note that we can query (and set) vertex and edge attributes with the V() and E() functions respectively.

``` r
V(g1)
```

    ## + 845/845 vertices, named, from bb42d27:
    ##   [1] ph_1061       ph_1258       ph_3164       ph_1033       ph_10996     
    ##   [6] ph_11038      ph_11040      ph_11048      ph_11096      ph_1113      
    ##  [11] ph_1208       ph_13207      ph_1346       ph_14679      ph_1572      
    ##  [16] ph_16045      ph_1909       ph_1918       ph_19894      ph_2117      
    ##  [21] ph_2231       ph_2363       ph_276        ph_2775       ph_2798      
    ##  [26] ph_3217       ph_3336       ph_3493       ph_3541       ph_3892      
    ##  [31] ph_4194       ph_4602       ph_4678       ph_484        ph_4993      
    ##  [36] ph_4999       ph_5001       ph_5010       ph_5286       ph_5287      
    ##  [41] ph_5302       ph_5321       ph_5643       ph_6441       ph_654       
    ##  [46] ph_6954       ph_7389       ph_7920       ph_8039       ph_8695      
    ## + ... omitted several vertices

``` r
E(g1)
```

    ## + 1544/1544 edges from bb42d27 (vertex names):
    ##  [1] ph_1061 --AACY020068177 ph_1258 --AACY020207233
    ##  [3] ph_3164 --AACY020207233 ph_1033 --AACY020255495
    ##  [5] ph_10996--AACY020255495 ph_11038--AACY020255495
    ##  [7] ph_11040--AACY020255495 ph_11048--AACY020255495
    ##  [9] ph_11096--AACY020255495 ph_1113 --AACY020255495
    ## [11] ph_1208 --AACY020255495 ph_13207--AACY020255495
    ## [13] ph_1346 --AACY020255495 ph_14679--AACY020255495
    ## [15] ph_1572 --AACY020255495 ph_16045--AACY020255495
    ## [17] ph_1909 --AACY020255495 ph_1918 --AACY020255495
    ## [19] ph_19894--AACY020255495 ph_2117 --AACY020255495
    ## + ... omitted several edges

Network community detection
---------------------------

Community structure detection algorithms try to find dense sub-graphs within larger network graphs (i.e. clusters of well connected nodes that are densely connected themselves but sparsely connected to other nodes outside the cluster) . Here we use the classic Girvan & Newman betweenness clustering method. The igraph package has lots of different community detection algorithms (i.e. different methods for finding communities).

``` r
cb <- cluster_edge_betweenness(g1)
```

    ## Warning in cluster_edge_betweenness(g1): At community.c:460 :Membership
    ## vector will be selected based on the lowest modularity score.

    ## Warning in cluster_edge_betweenness(g1): At community.c:467 :Modularity
    ## calculation with weighted edge betweenness community detection might not
    ## make sense -- modularity treats edge weights as similarities while edge
    ## betwenness treats them as distances

``` r
cb
```

    ## IGRAPH clustering edge betweenness, groups: 18, mod: 0.82
    ## + groups:
    ##   $`1`
    ##   [1] "ph_1061"       "AACY020068177"
    ##   
    ##   $`2`
    ##    [1] "ph_1258"       "ph_5861"       "ph_7172"       "ph_11569"     
    ##    [5] "ph_1291"       "ph_1600"       "ph_2702"       "ph_5790"      
    ##    [9] "ph_5858"       "ph_7594"       "ph_7816"       "ph_784"       
    ##   [13] "ph_1359"       "ph_1534"       "ph_1874"       "ph_2465"      
    ##   [17] "ph_5453"       "ph_900"        "ph_908"        "ph_811"       
    ##   [21] "ph_1367"       "ph_1452"       "ph_1458"       "ph_1723"      
    ##   + ... omitted several groups/vertices

``` r
plot(cb, y=g1, vertex.label=NA,  vertex.size=3)
```

![](Class17_files/figure-markdown_github/unnamed-chunk-11-1.png)

You can extract a cluster/community membership vector for further inspection with the `membership()` function:

``` r
head(membership(cb))
```

    ##  ph_1061  ph_1258  ph_3164  ph_1033 ph_10996 ph_11038 
    ##        1        2        3        4        4        4

Node degree
-----------

The degree of a node or vertex is its most basic structural property, the number of its adjacent edges. Here we calculate and plot the node degree distribution.

``` r
#Calculate and plot node degree of our network.
d <- degree(g1)
hist(d, breaks=30, col="lightblue", main ="Node Degree Distribution")
```

![](Class17_files/figure-markdown_github/unnamed-chunk-13-1.png)

For the `degree_distribution()` function a numeric vector of the same length as the maximum degree plus one is returned. The first element is the relative frequency zero degree vertices, the second vertices with degree one, et

``` r
plot(degree_distribution(g1), type="h")
```

![](Class17_files/figure-markdown_github/unnamed-chunk-14-1.png)

Centrality analysis
-------------------

Centrality gives an estimation on how important a node or edge is for the connectivity (or the information flow) of a network. It is a particularly useful parameter in signaling networks and it is often used when trying to find drug targets for example.

Centrality analysis often aims to answer the following question: Which nodes are the most important and why?

``` r
pr <- page_rank(g1)
head(pr$vector)
```

    ##      ph_1061      ph_1258      ph_3164      ph_1033     ph_10996 
    ## 0.0011834320 0.0011599483 0.0019042088 0.0005788564 0.0005769663 
    ##     ph_11038 
    ## 0.0005745460

Lets plot our network with nodes size scaled via this page rank centrality scores.

``` r
# Make a size vector btwn 2 and 20 for node plotting size.
v.size <- BBmisc::normalize(pr$vector, range=c(2,20), method="range")
plot(g1, vertex.size=v.size, vertex.label=NA)
```

![](Class17_files/figure-markdown_github/unnamed-chunk-16-1.png)

One of the simplest centrality scores is of course degree that we calculated previously and stored as the object `d`. Lets plot this one out also.

``` r
v.size <- BBmisc::normalize(d, range=c(2,20), method="range")
plot(g1, vertex.size=v.size, vertex.label=NA)
```

![](Class17_files/figure-markdown_github/unnamed-chunk-17-1.png)

Another very common centrality score is betweenness. The vertex and edge betweenness are (roughly) defined by the number of geodesics (shortest paths) going through a vertex or an edge.

``` r
b <- betweenness(g1)
v.size <- BBmisc::normalize(b, range=c(2,20), method="range")
plot(g1, vertex.size=v.size, vertex.label=NA)
```

![](Class17_files/figure-markdown_github/unnamed-chunk-18-1.png)
