#!/bin/bash
set -euo
TAG="annotate"
docker build --progress=plain -t fcunial/hapestry .
docker tag fcunial/hapestry fcunial/hapestry:${TAG}
docker push fcunial/hapestry:${TAG}
