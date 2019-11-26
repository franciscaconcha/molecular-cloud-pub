#!/bin/bash

# $1: framerate, 5 or less looks good
# $2: path for files
# $3: output name (including extension)

ffmpeg -framerate $1 -i $2/%04d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p $2/$3