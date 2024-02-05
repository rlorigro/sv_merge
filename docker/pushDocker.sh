#!/bin/bash
#
TAG=""
docker build --progress=plain -t fcunial/hapestry .
docker push fcunial/hapestry${TAG}
