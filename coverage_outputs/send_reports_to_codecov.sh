##Grab the uploader
curl -Os https://uploader.codecov.io/latest/linux/codecov

##make it executable
chmod +x codecov

##You'll need the WODEN_CODECOV_TOKEN environment variable. Read the
##contribution guidelines on how to obtain
./codecov --sha=$(git rev-parse HEAD) \
  -t ${WODEN_CODECOV_TOKEN} -f 'coverage_outputs/*' \
  -f '!*.sh' -f '!codecov'
