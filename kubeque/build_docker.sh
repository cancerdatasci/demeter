set +ex

docker build . -t us.gcr.io/broad-achilles/demeter
gcloud docker -- push us.gcr.io/broad-achilles/demeter
