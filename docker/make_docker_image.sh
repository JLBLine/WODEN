
woden_version="2.5"

#Build some CUDA versions. Could we just do this with a multi-arch flag?
for arch in "60" "61" "70" "75" "80" "86"
# for arch in "80"
do
    docker build --no-cache --progress=plain --build-arg="CUDA_ARCH=${arch}" -t jlbline/woden-${woden_version}:cuda-$arch -f Dockerfile_cuda .
    docker push jlbline/woden-${woden_version}:cuda-$arch
done

##Build the setonix version. Someone could replace this with a HIP equivalent
##of the CUDA stuff above if they can work out how to do it.
docker build --no-cache --progress=plain -t jlbline/woden-${woden_version}:setonix -f Dockerfile_setonix .
docker push jlbline/woden-${woden_version}:setonix


