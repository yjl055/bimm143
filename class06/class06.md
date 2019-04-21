Class 06: Why, when and how of writing my own R functons
================

Section 1: Improving analysis code by writing functions
-------------------------------------------------------

``` r
# Improve this analysis code
df <- data.frame(a=1:10, b=seq(200,400,length=10),c=11:20,d=NA)
df$a <- (df$a - min(df$a)) / (max(df$a) - min(df$a))
df$b <- (df$b - min(df$a)) / (max(df$b) - min(df$b))
df$c <- (df$c - min(df$c)) / (max(df$c) - min(df$c))
df$d <- (df$d - min(df$d)) / (max(df$a) - min(df$d)) 
```

Simplified the code first
=========================

x &lt;- (x-min(x))/(max(x)-min(x)) \# Reduced calculation duplication xmin &lt;- min(x) x &lt;- (x-xmin)/(max(x)-xmin) rng &lt;- range(x) x &lt;- (x-rng\[1\])/(rng\[2\]-rng\[1\])

Turned the reduced code into a function
=======================================

``` r
rescale <- function(x) {
  rng <- range(x)
  (x-rng[1]) / (rng[2]-rng[1])
}

rescale(1:10)
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

Improve the code for the analysis of protein drug interactions
--------------------------------------------------------------

``` r
library(bio3d)
s1 <- read.pdb("4AKE") # kinase with drug
```

    ##   Note: Accessing on-line PDB file

``` r
s2 <- read.pdb("1AKE") # kinase no drug
```

    ##   Note: Accessing on-line PDB file
    ##    PDB has ALT records, taking A only, rm.alt=TRUE

``` r
s3 <- read.pdb("1E4Y") # kinase with drug
```

    ##   Note: Accessing on-line PDB file

``` r
s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
s2.chainA <- trim.pdb(s2, chain="A", elety="CA")
s3.chainA <- trim.pdb(s1, chain="A", elety="CA")
s1.b <- s1.chainA$atom$b
s2.b <- s2.chainA$atom$b
s3.b <- s3.chainA$atom$b
plotb3(s1.b, sse=s1.chainA, typ="l", ylab="Bfactor")
```

![](class06_files/figure-markdown_github/unnamed-chunk-3-1.png)

``` r
plotb3(s2.b, sse=s2.chainA, typ="l", ylab="Bfactor")
```

![](class06_files/figure-markdown_github/unnamed-chunk-3-2.png)

``` r
plotb3(s3.b, sse=s3.chainA, typ="l", ylab="Bfactor")
```

![](class06_files/figure-markdown_github/unnamed-chunk-3-3.png)

``` r
plot_protein <- function (s) {
  s <- read.pdb(s) #s=name of the protein from data base
  s.chainA <- trim.pdb(s, chain = "A",elety="CA") #input = s
  s.b <- s.chainA$atom$b #input=s.chainA=trimmed s, assigning beta chain of residue from trimmed as s.b
  
  plotb3(s.b, sse=s.chainA, type="l", ylab="Bfactor") # inputs = s.b. & s.chianA, output = plot of residues
}
```

``` r
plot_protein
```

    ## function (s) {
    ##   s <- read.pdb(s) #s=name of the protein from data base
    ##   s.chainA <- trim.pdb(s, chain = "A",elety="CA") #input = s
    ##   s.b <- s.chainA$atom$b #input=s.chainA=trimmed s, assigning beta chain of residue from trimmed as s.b
    ##   
    ##   plotb3(s.b, sse=s.chainA, type="l", ylab="Bfactor") # inputs = s.b. & s.chianA, output = plot of residues
    ## }

``` r
plot_protein("4AKE")
```

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): C:
    ## \Users\yenaj\AppData\Local\Temp\RtmpAbUIaq/4AKE.pdb exists. Skipping
    ## download

![](class06_files/figure-markdown_github/unnamed-chunk-6-1.png)

``` r
plot_protein("1AKE")
```

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): C:
    ## \Users\yenaj\AppData\Local\Temp\RtmpAbUIaq/1AKE.pdb exists. Skipping
    ## download

    ##    PDB has ALT records, taking A only, rm.alt=TRUE

![](class06_files/figure-markdown_github/unnamed-chunk-7-1.png)

``` r
plot_protein("1E4Y")
```

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): C:
    ## \Users\yenaj\AppData\Local\Temp\RtmpAbUIaq/1E4Y.pdb exists. Skipping
    ## download

![](class06_files/figure-markdown_github/unnamed-chunk-8-1.png)
