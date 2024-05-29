#!/bin/bash
set -euo
TAG="annotate"
docker build --progress=plain -t fcunial/hapestry .
docker push fcunial/hapestry:${TAG}
