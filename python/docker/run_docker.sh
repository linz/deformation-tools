#!/bin/sh
docker run -it --rm -v `realpath ..`:/data python3-gdal:latest
