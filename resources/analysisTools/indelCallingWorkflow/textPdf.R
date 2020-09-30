#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
    text <- ""
} else {
    text <- args[1]
}

write(text, stdout())

pdf("/dev/stdout", height=3, width=6)
plot.new()
mtext(text)
ignore <- dev.off()
