---
title: "GWAS Tutorial"
subtitle: "For Anaerobic Germination in Rice"
author: "Yen-Hsiang (Teddy) Huang"
date: "2024-12-16"
knit: "bookdown::render_book"
site: bookdown::bookdown_site
output: bookdown::bs4_book
documentclass: book
bibliography: [book.bib]
biblio-style: apalike
link-citations: true
links-as-notes: true
colorlinks: true
github-repo: dgrtwo/tidy-text-mining
cover-image: images/cover.png
url: https://www.tidytextmining.com/
description: "A guide to GWAS analysis in rice, using the R and PLINK programs"
---



# Welcome to GWAS Tutorial {.unnumbered}

This is the [website](http://tidytextmining.com/) for ***GWAS analysis with R and PLINK programs***.

This work by **Yen-Hsiang (Teddy) Huang** is licensed under a <a rel="license" href="https://www.gnu.org/licenses/fdl-1.3.html.en#license-text">GNU Free Documentation License</a>.

<a href="https://www.gnu.org/licenses/fdl-1.3.html.en#license-text"> ![](images/clipboard-1489856658.png)</a>

```{=html}
<script>
  (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
  (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
  m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
  })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');

  ga('create', 'UA-68765210-2', 'auto');
  ga('send', 'pageview');

</script>
```
# Preface {.unnumbered}

If you work in analytics or data science, like we do, you are familiar with the fact that data is being generated all the time at ever faster rates. (You may even be a little weary of people pontificating about this fact.) Analysts are often trained to handle tabular or rectangular data that is mostly numeric, but much of the data proliferating today is unstructured and text-heavy. Many of us who work in analytical fields are not trained in even simple interpretation of natural language.

We developed the [tidytext](https://github.com/juliasilge/tidytext) [@R-tidytext] R package because we were familiar with many methods for data wrangling and visualization, but couldn't easily apply these same methods to text. We found that using tidy data principles can make many text mining tasks easier, more effective, and consistent with tools already in wide use. Treating text as data frames of individual words allows us to manipulate, summarize, and visualize the characteristics of text easily and integrate natural language processing into effective workflows we were already using.

This book serves as an introduction of text mining using the tidytext package and other tidy tools in R. The functions provided by the tidytext package are relatively simple; what is important are the possible applications. Thus, this book provides compelling examples of real text mining problems.

## Outline {.unnumbered}

We start by introducing the tidy text format, and some of the ways dplyr, tidyr, and tidytext allow informative analyses of this structure.

-   **Chapter** \@ref(tidytext) outlines the tidy text format and the `unnest_tokens()` function. It also introduces the gutenbergr and janeaustenr packages, which provide useful literary text datasets that we'll use throughout this book.
-   **Chapter** \@ref(sentiment) shows how to perform sentiment analysis on a tidy text dataset, using the `sentiments` dataset from tidytext and `inner_join()` from dplyr.
-   **Chapter** \@ref(tfidf) describes the tf-idf statistic (term frequency times inverse document frequency), a quantity used for identifying terms that are especially important to a particular document.
-   **Chapter** \@ref(ngrams) introduces n-grams and how to analyze word networks in text using the widyr and ggraph packages.

Text won't be tidy at all stages of an analysis, and it is important to be able to convert back and forth between tidy and non-tidy formats.

-   **Chapter** \@ref(dtm) introduces methods for tidying document-term matrices and corpus objects from the tm and quanteda packages, as well as for casting tidy text datasets into those formats.
-   **Chapter** \@ref(topicmodeling) explores the concept of topic modeling, and uses the `tidy()` method to interpret and visualize the output of the topicmodels package.

We conclude with several case studies that bring together multiple tidy text mining approaches we've learned.

-   **Chapter** \@ref(twitter) demonstrates an application of a tidy text analysis by analyzing the authors' own Twitter archives. How do Dave's and Julia's tweeting habits compare?
-   **Chapter** \@ref(nasa) explores metadata from over 32,000 NASA datasets (available in JSON) by looking at how keywords from the datasets are connected to title and description fields.
-   **Chapter** \@ref(usenet) analyzes a dataset of Usenet messages from a diverse set of newsgroups (focused on topics like politics, hockey, technology, atheism, and more) to understand patterns across the groups.

## About this book {.unnumbered}

This book is focused on practical software examples and data explorations. There are few equations, but a great deal of code. We especially focus on generating real insights from the literature, news, and social media that we analyze.

We don't assume any previous knowledge of text mining. Professional linguists and text analysts will likely find our examples elementary, though we are confident they can build on the framework for their own analyses.

We do assume that the reader is at least slightly familiar with dplyr, ggplot2, and the `%>%` "pipe" operator in R, and is interested in applying these tools to text data. For users who don't have this background, we recommend books such as [R for Data Science](http://r4ds.had.co.nz/). We believe that with a basic background and interest in tidy data, even a user early in their R career can understand and apply our examples.

## Using code examples {.unnumbered}

This book was written in [RStudio](http://www.rstudio.com/ide/) using [bookdown](http://bookdown.org/). The [website](https://www.tidytextmining.com/) is hosted via [Netlify](http://netlify.com/), and automatically built after every push by [GitHub Actions](https://help.github.com/actions). While we show the code behind the vast majority of the analyses, in the interest of space we sometimes choose not to show the code generating a particular visualization if we've already provided the code for several similar graphs. We trust the reader can learn from and build on our examples, and the code used to generate the book can be found in our [public GitHub repository](https://github.com/dgrtwo/tidy-text-mining). We generated all plots in this book using [ggplot2](https://ggplot2.tidyverse.org/) and its light theme (`theme_light()`).

This version of the book was built with R version 4.4.1 (2024-06-14 ucrt) and the following packages:



|package     |version   |source         |
|:-----------|:---------|:--------------|
|bookdown    |0.40      |CRAN (R 4.4.1) |
|dplyr       |1.1.4     |CRAN (R 4.4.1) |
|forcats     |1.0.0     |CRAN (R 4.4.1) |
|ggforce     |0.4.2     |CRAN (R 4.4.1) |
|ggplot2     |3.5.1     |CRAN (R 4.4.1) |
|ggraph      |2.2.1     |CRAN (R 4.4.1) |
|gutenbergr  |NA        |NA             |
|igraph      |2.1.1     |CRAN (R 4.4.1) |
|janeaustenr |1.0.0     |CRAN (R 4.4.1) |
|jsonlite    |1.8.9     |CRAN (R 4.4.2) |
|lubridate   |1.9.3     |CRAN (R 4.4.1) |
|mallet      |NA        |NA             |
|Matrix      |1.7-0     |CRAN (R 4.4.0) |
|quanteda    |NA        |NA             |
|readr       |2.1.5     |CRAN (R 4.4.1) |
|reshape2    |1.4.4     |CRAN (R 4.4.1) |
|sessioninfo |1.2.2     |CRAN (R 4.4.1) |
|stringr     |1.5.1     |CRAN (R 4.4.1) |
|styler      |1.10.3    |CRAN (R 4.4.1) |
|textdata    |NA        |NA             |
|tidyr       |1.3.1     |CRAN (R 4.4.1) |
|tidytext    |0.4.2     |CRAN (R 4.4.1) |
|tidyverse   |2.0.0     |CRAN (R 4.4.1) |
|tm          |NA        |NA             |
|topicmodels |NA        |NA             |
|widyr       |NA        |NA             |
|wordcloud   |2.6       |CRAN (R 4.4.2) |
|XML         |3.99-0.17 |CRAN (R 4.4.1) |



## Acknowledgements {.unnumbered}

We are so thankful for the contributions, help, and perspectives of people who have moved us forward in this project. There are several people and organizations we would like to thank in particular.

We would like to thank Os Keyes and Gabriela de Queiroz for their contributions to the tidytext package, Lincoln Mullen for his work on the [tokenizers](https://github.com/ropensci/tokenizers) package, Kenneth Benoit for his work on the [quanteda](https://github.com/kbenoit/quanteda) package, Thomas Pedersen for his work on the [ggraph](https://github.com/thomasp85/ggraph) package, and Hadley Wickham for his work in framing tidy data principles and building tidy tools. We would also like to thank [rOpenSci](https://ropensci.org/), which hosted us at the unconference where we began work, and the [NASA Datanauts](https://open.nasa.gov/explore/datanauts/) program, for the opportunities and support they have provided Julia during her time with them.

We received thoughtful, thorough technical reviews that improved the quality of this book significantly. We would like to thank Mara Averick, Carolyn Clayton, Simon Jackson, Sean Kross, and Lincoln Mullen for their investment of time and energy in these technical reviews.

Finally, we want to dedicate this book to our spouses, Robert and Dana. We both could produce a great deal of sentimental text on this subject but will restrict ourselves to heartfelt thanks.
