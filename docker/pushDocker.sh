#!/bin/bash
set -euxo

if [ $# -eq 0 ]; then
  echo "Need to provide tag for the docker image."
  exit 1
fi

TAG="$1"

docker build --platform linux/amd64 --progress=plain -t fcunial/hapestry .
docker tag fcunial/hapestry fcunial/hapestry:${TAG}
docker push fcunial/hapestry:${TAG}
